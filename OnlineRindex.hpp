/*!
 * Copyright (c) 2018 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file OnlineRindex.hpp
 * @brief Online r-index, index based on Run-length encoded Burrows–Wheeler transform (RLBWT).
 * @author Tomohiro I
 * @date 2018-05-04
 */
#ifndef INCLUDE_GUARD_OnlineRindex
#define INCLUDE_GUARD_OnlineRindex

#include <stdint.h>
#include <cassert>
#include <iostream>
#include <fstream>

namespace itmmti 
{
  /*!
   * @brief Online r-index, index based on Run-length encoded Burrows–Wheeler transform (RLBWT).
   * @note
   *   ::OnlineRindex wraps ::DynRLE to use it for representing dynamic RLE of BWT.
   *   In contrast to ::DynRle, OnlineRLWT has an implicit end marker (em_) at emPos_.
   */
  template<class DynRle, class DynSuccT>
  class OnlineRlbwtIndex
  {
  public:
    //// bwtintvl: Intervals are [left..right) : right bound is excluded.
    using bwtintvl = std::pair<uint64_t, uint64_t>;
    //// pattracker: The first two uints represent a bwt-interval [left..right),
    //// the third uint is tracking the idxS from which we can obtain sampled position for "BWT[left]",
    //// and the last uint represents how far we are from the last sampled position.
    using pattracker = std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>;
    using CharT = typename DynRle::CharT;
    using BTreeNodeT = typename DynRle::BTreeNodeT;
    static constexpr uintptr_t NOTFOUND{BTreeNodeT::NOTFOUND};
    static constexpr uint8_t kB{DynRle::kB};
    static constexpr uint8_t kBtmBM{DynRle::kBtmBM};
    static constexpr uint8_t kBtmBS{DynRle::kBtmBS};


  private:
    DynRle drle_; //!< Dynamic RLE data structure with one sampling position per run.
    DynSuccT succ_; //!< Dynamic data structure (supporting indel) for successor queries.
    uint64_t emPos_; //!< Current position (0base) of end marker.
    uint64_t prevSamplePos_; //!< Tracking txt-pos (0base) for character at bwt-position previous to current emPos_.
    uint64_t nextSamplePos_; //!< Tracking txt-pos (0base) for character at bwt-position next to current emPos_.
    CharT em_; //!< End marker. It is used only when bwt[emPos_] is accessed (and does not matter if em_ appears in the input text).


  public:
    OnlineRlbwtIndex
    (
     const size_t initNumBtms, //!< Initial size of DynRle to reserve.
     const uint64_t initSampleUb = 256, //!< Initial upper bound (excluding) of sample positions.
     CharT em = 0 //!< End marker (default 0).
     ) :
      drle_(),
      succ_(),
      emPos_(0),
      prevSamplePos_(0),
      nextSamplePos_(0),
      em_(em)
    {
      drle_.init_rindex(initNumBtms, initSampleUb);
    }


    ~OnlineRlbwtIndex()
    {
      //// drle_ and succ_ are deleted.
    }


    /*!
     * @brief Initialize data structure with empty elements.
     */
    void init
    (
     const size_t initNumBtms, //!< Initial size of DynRle to reserve.
     const uint64_t initSampleUb = 256 //!< Initial upper bound (excluding) of sample positions.
     ) {
      if (isReady()) {
        clearAll();
      }
      drle_.init_rindex(initNumBtms, initSampleUb);
      succ_.init();
    }


    /*!
     * @brief Free/delete all allocated objects.
     */
    void clearAll() {
      if (!isReady()) { // already cleared
        return;
      }
      drle_.clearAll();
      succ_.clearAll();
    }


    /*!
     * @brief Return if data structure is ready.
     */
    bool isReady() const noexcept {
      return drle_.isReady();
    }


    /*!
     * @brief Get end marker.
     */
    CharT getEm() const noexcept {
      return em_;
    }


    /*!
     * @brief Get current position of end marker.
     */
    uint64_t getEndmarkerPos() const noexcept {
      return emPos_;
    }


    /*!
     * @brief Return current length including end marker.
     */
    uint64_t getLenWithEndmarker() const noexcept {
      return drle_.getSumOfWeight() + 1; // +1 for end marker, which is not in drle_.
    }


    /*!
     * @brief Get current "prevSamplePos_".
     */
    uint64_t getPrevSamplePos() const noexcept {
      return prevSamplePos_;
    }


    /*!
     * @brief Get current "nextSamplePos_".
     */
    uint64_t getNextSamplePos() const noexcept {
      return nextSamplePos_;
    }

  
    /*!
     * @brief Get associated value at "idxM".
     */
    uint64_t getSample(uint64_t idxM) const noexcept {
      return drle_.getSampleFromIdxM(idxM);
    }


    /*!
     * @brief Set associated value at "idxM".
     */
    void setSample(uint64_t idxM, uint64_t val) noexcept {
      drle_.setSample(idxM, val);
    }


    /*!
     * @brief Extend RLBWT by appending one character while maintaining r-index.
     */
    void extend
    (
     const CharT ch //!< Char to append.
     ) {
      const uint64_t txtPos = drle_.getSumOfWeight(); //!< Txt-position of "ch" (0base).
      {
        const auto sampleUb = drle_.getSampleUb();
        if (sampleUb && txtPos >= sampleUb) {
          drle_.increaseSampleUb(txtPos + 1);
        }
      }

      uint64_t idxM;
      auto pos = emPos_;
      bool split = false;
      uint64_t originalSample = 0;
      if (pos < drle_.getSumOfWeight()) {
        idxM = drle_.searchPosM(pos);
        if ((pos > 0) && (drle_.getCharFromIdxM(idxM) != ch)) { // split a run
          split = true;
          originalSample = drle_.getSampleFromIdxM(idxM);
        }
        idxM = drle_.insertRun(idxM, pos, ch); // 'pos' is modified to be the relative pos in the run of 'idxM'.
      } else {
        idxM = drle_.pushbackRun_rindex(pos, ch);
      }
      emPos_ = drle_.rank(ch, idxM, pos, true);

      const uint64_t exponent = drle_.getWeightFromIdxM(idxM);
      // std::cout << __func__ << ": pos = " << pos << ", exponent = " << exponent << ", split = " << split
      //           << ", prev = " << prevSamplePos_ << ", next = " << nextSamplePos_ << std::endl;
      if (pos + 1 != exponent) {
        ++nextSamplePos_;
      } else {
        if (exponent > 1) {
          succ_.removeKey(prevSamplePos_);
        }
        setSample(drle_.getNextIdxM(idxM), txtPos);
        succ_.setKeyVal(txtPos, nextSamplePos_);
        //// Update nextSamplePos_
        nextSamplePos_ = 0;
        uint64_t idxS = drle_.idxM2S(idxM);
        const uint64_t nextIdxS = drle_.getNextIdxS(idxS);
        if (nextIdxS != BTreeNodeT::NOTFOUND) { // Succcessor with "ch" was found.
          const auto pos = drle_.getSampleFromIdxS(nextIdxS);
          nextSamplePos_ = succ_.nextPos(pos) + 1;
        } else { // Succcessor with "ch" was NOT found.
          const auto nextRootS = drle_.getNextRootS(drle_.getParentFromBtmS(idxS / kBtmBS));
          if (reinterpret_cast<uintptr_t>(nextRootS) != BTreeNodeT::NOTFOUND) {
            nextSamplePos_ = 1;
            const auto btmS = reinterpret_cast<uint64_t>(nextRootS->getLmBtm_DirectJump());
            if (btmS) {
              const auto pos = drle_.getSampleFromIdxS(btmS * kBtmBS + 1);
              nextSamplePos_ += succ_.nextPos(pos);
            }
          }
        }
      }

      if (pos == 0) {
        if (exponent == 1) {
          setSample(idxM, prevSamplePos_);
        }
        if (split) {
          setSample(drle_.getPrevIdxM(idxM), originalSample);
        }
        succ_.setKeyVal(prevSamplePos_, txtPos);
        //// Update prevSamplePos_
        prevSamplePos_ = 0;
        uint64_t idxS = drle_.getPrevIdxS(drle_.idxM2S(idxM));
        uint64_t prevIdxM = drle_.idxS2M(idxS);
        if (prevIdxM == 0) { // Predecessor with "ch" was NOT found.
          const auto prevRootS = drle_.getPrevRootS(drle_.getParentFromBtmS(idxS / kBtmBS));
          if (!prevRootS->isDummy()) {
            const uint64_t btmS = reinterpret_cast<uint64_t>(prevRootS->getRmBtm());
            idxS = btmS * kBtmBS + drle_.getNumChildrenFromBtmS(btmS) - 1;
            prevIdxM = drle_.idxS2M(idxS);
          }
        }
        if (prevIdxM) {
          prevIdxM = drle_.getNextIdxM(prevIdxM);
          prevSamplePos_ = drle_.getSampleFromIdxM(prevIdxM) + 1;
        }
      } else {
        ++prevSamplePos_;
      }
    }


    /*!
     * @brief Access to the current RLBWT by [] operator.
     */
    CharT operator[]
    (
     uint64_t pos //!< in [0..getLenWithEndmarker()].
     ) const noexcept {
      assert(pos < getLenWithEndmarker());

      if (pos == emPos_) {
        return em_;
      }
      pos -= (pos > emPos_);
      uint64_t idxM = drle_.searchPosM(pos);
      return drle_.getCharFromIdxM(idxM);
    }


    /*!
     * @brief Return 'rank of ch at pos' + 'num of total occ of characters smaller than ch'.
     */
    uint64_t totalRank
    (
     const CharT ch,
     uint64_t pos //!< in [0..getLenWithEndmarker()].
     ) const noexcept {
      assert(pos < getLenWithEndmarker());

      pos -= (pos > emPos_);
      return drle_.rank(ch, pos, true);
    }


    
    /*!
     * @brief Function to get text-position for "BWT[bwtpos + 1]", where "BWT[bwtpos]" corresponds to "T[txtpos]".
     */
    uint64_t nextPos
    (
     uint64_t txtpos //!< Text-position for currently focused character.
     ) const noexcept {
      return succ_.nextPos(txtpos);
    }


    /*!
     * @brief Compute the first occ from a valid pattracker.
     */
    uint64_t calcFstOcc
    (
     const pattracker & tracker
     ) const noexcept {
      const uint64_t idxS = std::get<2>(tracker);
      const uint64_t key = drle_.getSampleFromIdxS(drle_.getNextIdxS(idxS));
      return succ_.nextPos(key) + std::get<3>(tracker);
    }


    /*!
     * @brief Compute bwt-interval for cW from bwt-interval for W.
     * @note Intervals are [left..right), where right bound is excluded.
     */
    bool lfMap
    (
     pattracker & tracker,
     const CharT ch
     ) const noexcept {
      // {//debug
      //   std::cerr << __func__ << ": ch = " << ch << ", tracker = {" << std::get<0>(tracker) << ", " << std::get<1>(tracker) << ", " << std::get<2>(tracker) << "}" << std::endl;
      // }
      assert(std::get<0>(tracker) <= getLenWithEndmarker() && std::get<1>(tracker) <= getLenWithEndmarker());
      assert(std::get<0>(tracker) < std::get<1>(tracker));

      //// If "ch" is not in the alphabet or empty interval, return empty interval.
      const auto * retRootS = drle_.searchCharA(ch);
      if (retRootS->isDummy() || drle_.getCharFromNodeS(retRootS) != ch) {
        return false;
      }

      uint64_t r_in_drle = std::get<1>(tracker) - (std::get<1>(tracker) > emPos_); // Taking (implicit) end-marker into account.
      //// +1 because in F we are not taking into account the end-marker,
      //// which is in position 0 but not explicitly stored in F. */
      r_in_drle = drle_.rank(ch, r_in_drle - 1, true) + 1;
      if (r_in_drle <= 1) {
        return false;
      }

      uint64_t l_in_drle = std::get<0>(tracker) - (std::get<0>(tracker) > emPos_); // Taking (implicit) end-marker into account.
      const uint64_t idxM = drle_.searchPosM(l_in_drle); // l_in_drle is modified to relative pos.
      //// In order to get idxS, replicate variant of rank function,
      //// where pos is specified by 'idxM' and 'relativePos' with several modifications.
      const auto chNow = drle_.getCharFromIdxM(idxM);
      uint64_t idxS;
      {
        if (ch == chNow) {
          idxS = drle_.idxM2S(idxM);
        } else {
          l_in_drle = 0;
          idxS = drle_.calcPredIdxSFromIdxM(retRootS, ch, idxM); // We know that "ch" already exists by "r_in_dre > 1".
        }
        const auto btmS = idxS / kBtmBS;
        for (auto tmpIdxS = btmS * kBtmBS; tmpIdxS < idxS + (ch != chNow); ++tmpIdxS) {
          l_in_drle += drle_.getWeightFromIdxS(tmpIdxS);
        }
        BTreeNodeT * rootS;
        l_in_drle += drle_.getParentFromBtmS(btmS)->calcPSum(drle_.getIdxInSiblingFromBtmS(btmS), rootS);
        l_in_drle += rootS->getParent()->calcPSum(rootS->getIdxInSibling());
      }
      ++l_in_drle; // +1 for (implicit) end-marker in F.

      if (l_in_drle >= r_in_drle) {
        return false;
      }

      //// Update tracker.
      std::get<0>(tracker) = l_in_drle;
      std::get<1>(tracker) = r_in_drle;
      if (l_in_drle == emPos_) {
        std::get<2>(tracker) = drle_.getSumOfWeight();
      } else if (ch == chNow) {
        std::get<3>(tracker) += 1;
      } else {
        std::get<2>(tracker) = idxS;
        std::get<3>(tracker) = 1;
      }

      return true;
    }


    /*!
     * @brief Compute bwt-interval for cW from bwt-interval for W
     * @note Intervals are [left, right) : right bound is excluded
     */
    bwtintvl lfMap
    (
     const bwtintvl intvl,
     const CharT ch
     ) const noexcept {
      assert(intvl.first <= getLenWithEndmarker() && intvl.second <= getLenWithEndmarker());

      //// If "ch" is not in the alphabet or empty interval, return empty interval.
      const auto * retRootS = drle_.searchCharA(ch);
      if (retRootS->isDummy() || drle_.getCharFromNodeS(retRootS) != ch || intvl.first >= intvl.second) {
        return {0, 0};
      }

      uint64_t l = intvl.first - (intvl.first > emPos_);
      uint64_t r = intvl.second - (intvl.second > emPos_);
      const uint64_t idxM = drle_.searchPosM(l); // l is modified to relative pos.
      //// +1 because in F.select(0, ch) we are not taking into account the end-marker,
      //// which is in position 0 but not explicitly stored in F.
      return {
        drle_.rank(ch, idxM, l, true) - (drle_.getCharFromIdxM(idxM) == ch) + 1,
          drle_.rank(ch, r-1, true) + 1
          };
    }


    /*!
     * @brief LF map.
     */
    uint64_t lfMap(uint64_t i) const noexcept {
      assert(i < getLenWithEndmarker());

      i -= (i > emPos_);
      const uint64_t idxM = drle_.searchPosM(i);
      const auto ch = drle_.getCharFromIdxM(idxM);
      return drle_.rank(ch, idxM, i, true);
    }


    /*!
     * @brief Output original text to std::ofstream.
     */
    void invert
    (
     std::ofstream & ofs
     ) const noexcept {
      uint64_t pos = 0;
      for (uint64_t i = 0; i < this->getLenWithEndmarker() - 1; ++i) {
        pos -= (pos > emPos_);
        const uint64_t idxM = drle_.searchPosM(pos);
        const auto ch = drle_.getCharFromIdxM(idxM);
        ofs.put(static_cast<unsigned char>(ch));
        pos = drle_.rank(ch, idxM, pos, true);
      }
    }


    //////////////////////////////// statistics
    /*!
     * @brief Calculate total memory usage in bytes.
     */
    size_t calcMemBytes
    (
     bool includeThis = true
     ) const noexcept {
      size_t size = sizeof(*this) * includeThis;
      size += drle_.calcMemBytes(false);
      size += succ_.calcMemBytes(false);
      return size;
    }


    /*!
     * @brief Print statistics of ::DynRLE (not of ::OnlineRLBWT).
     */
    void printStatistics
    (
     std::ostream & os, //!< Output stream, e.g., std::cout.
     const bool verbose
     ) const noexcept {
      os << "OnlineRindex object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
      os << "emPos_ = " << emPos_ << ", em_ = " << em_ << std::endl;
      if (isReady()) {
        const size_t totalBytes = drle_.calcMemBytes(true) + succ_.calcMemBytes(true);
        os << "Total memory usage: " << totalBytes << " bytes = "
           << (double)(totalBytes) / 1024 << " KiB = "
           << ((double)(totalBytes) / 1024) / 1024 << " MiB" << std::endl;
        drle_.printStatistics(os, verbose);
        succ_.printStatistics(os, verbose);
      } else {
        os << "Data structure is empty (not ready)." << std::endl;
      }
      os << "OnlineRindex object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
    }


    void printDebugInfo
    (
     std::ostream & os //!< Output stream, e.g., std::cout.
     ) const noexcept {
      os << "emPos_ = " << emPos_ << ", em_ = " << em_ << std::endl;
      // drle_.printDebugInfo(os);
      succ_.printDebugInfo(os);

      //// Check correctness of sampled position links.
      uint64_t bwtpos = lfMap(0);
      for (uint64_t len = 1; len < getLenWithEndmarker() - 1; ++len) {
        uint64_t temp = bwtpos - (bwtpos > emPos_);
        uint64_t idxM = drle_.searchPosM(temp); // temp is modified
        // os << "check: len = " << len << ", bwtpos = " << bwtpos << ", fbwtpos = " << bwtpos - (bwtpos > emPos_)
        //    << ", temp = " << temp << ", weight = " << drle_.getWeightFromIdxM(idxM) << ", sample = " << drle_.getSampleFromIdxM(idxM) << std::endl;
        if (temp == 0) {
          const uint64_t prevPos = drle_.getSampleFromIdxM(idxM);
          const uint64_t nextPos = succ_.nextPos(prevPos);
          if (len != nextPos) {
            os << "error: bwtpos = " << bwtpos - (bwtpos > emPos_) << ", len = " << len << ", nextPos = " << nextPos << std::endl;
            drle_.printDebugInfoOfBtmM(idxM / kBtmBM, os);
            succ_.printStatistics(os, true);
          }
        }
        if (temp + 1 == drle_.getWeightFromIdxM(idxM)) {
          idxM = drle_.getNextIdxM(idxM);
          const uint64_t sample = drle_.getSampleFromIdxM(idxM);
          if (len != sample) {
            os << "error: bwtpos = " << bwtpos - (bwtpos > emPos_) << ", len = " << len << ", sample = " << sample << std::endl;
          }
        }
        bwtpos = lfMap(bwtpos);
      }
    }
  };
};

#endif

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
    //// BwtIntvl: Intervals are [left..right) : right bound is excluded.
    using BwtIntvl = std::pair<uint64_t, uint64_t>;
    //// PatTracker: The first two uints represent a bwt-interval [left..right),
    //// the third uint is tracking the idxS from which we can obtain sampled position for "BWT[left]",
    //// and the last uint represents how far we are from the last sampled position.
    using PatTracker = std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>;
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
    uint64_t lastSamplePos_; //!< Tracking txt-pos (0base) for last bwt-position.
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
      em_(em)
    {
      if (initNumBtms) {
        drle_.init(initNumBtms, initSampleUb);
      }
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
      if (initNumBtms) {
        drle_.init(initNumBtms, initSampleUb);
        succ_.init();
        emPos_ = prevSamplePos_ = nextSamplePos_ = lastSamplePos_ = 0;
      }
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
     * @brief Return current length excluding end marker.
     */
    uint64_t getLenWithoutEndmarker() const noexcept {
      return drle_.getSumOfWeight();
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
      // std::cerr << __func__ << ": idxM = " << idxM << ", val = " << val << std::endl;

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
      auto bwtPos = emPos_;
      bool split = false;
      uint64_t originalSample = 0;
      if (bwtPos < drle_.getSumOfWeight()) {
        idxM = drle_.searchPosM(bwtPos);
        if ((bwtPos > 0) && (drle_.getCharFromIdxM(idxM) != ch)) { // split a run
          split = true;
          originalSample = drle_.getSampleFromIdxM(idxM);
        }
        idxM = drle_.insertRun(idxM, bwtPos, ch); // 'bwtPos' is modified to be the relative pos in the run of 'idxM'.
      } else {
        idxM = drle_.pushbackRun(bwtPos, ch);
      }
      emPos_ = drle_.rank(ch, idxM, bwtPos, true);

      const uint64_t exponent = drle_.getWeightFromIdxM(idxM);

      const auto prevSamplePos_old = prevSamplePos_;
      const auto nextSamplePos_old = nextSamplePos_;
      if (bwtPos == 0) {
        if (bwtPos == 0 && idxM > 1) {
          if (exponent == 1) {
            setSample(idxM, prevSamplePos_old);
          }
          if (split) {
            setSample(drle_.getPrevIdxM(idxM), originalSample);
          }
          succ_.setKeyVal(prevSamplePos_old, txtPos);
        }
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
          if (prevIdxM != BTreeNodeT::NOTFOUND) {
            prevSamplePos_ = drle_.getSampleFromIdxM(prevIdxM) + 1;
          } else {
            prevSamplePos_ = lastSamplePos_ + 1;
          }
        }
      } else {
        ++prevSamplePos_;
      }

      //// Update nextSamplePos_
      if (bwtPos + 1 < exponent) {
        ++nextSamplePos_;
      } else if (emPos_ == txtPos + 1) {
        nextSamplePos_ = txtPos + 2;
      } else {
        nextSamplePos_ = succ_.calcNextPos(prevSamplePos_, txtPos, prevSamplePos_old, nextSamplePos_old);
      }

      //// Update successor data structure.
      if (bwtPos + 1 == exponent) {
        const auto nextIdxM = drle_.getNextIdxM(idxM);
        if (nextIdxM != BTreeNodeT::NOTFOUND) {
          if (exponent > 1) {
            succ_.removeKey(prevSamplePos_old);
          }
          setSample(nextIdxM, txtPos);
          succ_.setKeyVal(txtPos, nextSamplePos_old);
        } else {
          lastSamplePos_ = txtPos;
        }
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
     * @brief Function to get text-position for "BWT[bwtPos + 1]", where "BWT[bwtPos]" corresponds to "T[txtPos]".
     */
    uint64_t calcNextPos
    (
     const uint64_t txtPos //!< Text-position for currently focused character.
     ) const noexcept {
      return succ_.calcNextPos(txtPos, getLenWithoutEndmarker(), prevSamplePos_, nextSamplePos_);
      // if (txtPos < tailBeg_) {
      //   return succ_.calcNextPos(txtPos, getLenWithoutEndmarker(), prevSamplePos_, nextSamplePos_);
      // } else { // tail part
      //   return nextSamplePos_ - (drle_.getSumOfWeight() - txtPos);
      // }
    }


    /*!
     * @brief Compute the first occ (ending position) from valid PatTracker.
     */
    uint64_t calcFstOcc
    (
     const PatTracker & tracker
     ) const noexcept {
      if (std::get<0>(tracker) == emPos_) {
        return getLenWithoutEndmarker();
      }
      const uint64_t idxS = std::get<2>(tracker);
      if (idxS) {
        const uint64_t key = drle_.getSampleFromIdxS(drle_.getNextIdxS(idxS));
        return calcNextPos(key) + std::get<3>(tracker);
      }
      return std::get<3>(tracker);
    }


    /*!
     * @brief Compute bwt-interval for cW from bwt-interval for W.
     * @note Intervals are [left..right), where right bound is exclusive.
     */
    bool lfMap
    (
     PatTracker & tracker, //!< Valid PatTracker.
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
      //// which is in position 0 but not explicitly stored in F.
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
      if (ch == chNow) {
        std::get<3>(tracker) += 1;
      } else {
        std::get<2>(tracker) = idxS;
        std::get<3>(tracker) = 1;
      }

      return true;
    }


    /*!
     * @brief Compute bwt-interval for cW from bwt-interval for W.
     * @note Intervals are [left..right), where right bound is exclusive.
     */
    BwtIntvl lfMap
    (
     const BwtIntvl intvl,
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


    /*!
     * @brief Get initial patter tracker based on current rindex.
     */
    PatTracker getInitialPatTracker() const noexcept {
      assert(isReady());

      return {0, getLenWithEndmarker(), 0, 0};
    }


    /*!
     * @brief Get number of occurrences from valid PatTracker.
     */
    bool includeEmPos(const PatTracker & tracker) const noexcept {
      assert(isReady());

      return (std::get<0>(tracker) <= emPos_ && emPos_ < std::get<1>(tracker));
    }


    /*!
     * @brief Get number of occurrences from valid PatTracker.
     */
    uint64_t getNumOcc(const PatTracker & tracker) const noexcept {
      assert(isReady());

      return std::get<1>(tracker) - std::get<0>(tracker);
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
      os << "emPos_ = " << emPos_ << ", em_ = " << em_
         << ", prevSamplePos_ = " << prevSamplePos_ << ", nextSamplePos_ = " << nextSamplePos_ << std::endl;
      if (isReady()) {
        const size_t totalBytes = drle_.calcMemBytes(true) + succ_.calcMemBytes(true);
        os << "Total: " << totalBytes << " bytes = "
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
      os << "emPos_ = " << emPos_ << ", em_ = " << em_ << ", prevSamplePos_ = " << prevSamplePos_
         << ", nextSamplePos_ = " << nextSamplePos_ << ", lastSamplePos_ = " << lastSamplePos_ << std::endl;
      // drle_.printDebugInfo(os);
      succ_.printDebugInfo(os);

      //// Check correctness of sampled position links.
      const auto len = getLenWithoutEndmarker();
      for (uint64_t txtPos = 0, bwtPos = 0; bwtPos < getLenWithEndmarker(); ++bwtPos) {
        if (bwtPos == emPos_) {
          if (txtPos != len) {
            os << "!error!: bwtPos = " << bwtPos << ", txtPos = " << txtPos << ", len = " << len << std::endl;
          }
        } else if (bwtPos == emPos_ + 1) {
          if (txtPos != nextSamplePos_) {
            os << "!error!: bwtPos = " << bwtPos << ", txtPos = " << txtPos << ", nextSamplePos_ = " << nextSamplePos_ << std::endl;
          }
        } else if (bwtPos == emPos_ - 1) {
          if (txtPos != prevSamplePos_) {
            os << "!error!: bwtPos = " << bwtPos << ", txtPos = " << txtPos << ", prevSamplePos_ = " << prevSamplePos_ << std::endl;
          }
        } else {
          uint64_t temp = bwtPos - (bwtPos > emPos_);
          auto idxM = drle_.searchPosM(temp);
          if (temp + 1 == drle_.getWeightFromIdxM(idxM)) {
            idxM = drle_.getNextIdxM(idxM);
            if (idxM != BTreeNodeT::NOTFOUND) {
              const uint64_t sample = drle_.getSampleFromIdxM(idxM);
              if (txtPos != sample) {
                os << "!error!: bwtPos = " << bwtPos << ", txtPos = " << txtPos << ", sample = " << sample << std::endl;
              }
            } else {
              if (txtPos != lastSamplePos_) {
                os << "!error!: bwtPos = " << bwtPos << ", txtPos = " << txtPos << ", lastSamplePos_ = " << lastSamplePos_ << std::endl;
              }
            }
          }
        }
        txtPos = calcNextPos(txtPos);
      }
    }


    bool checkDecompress
    (
     std::ifstream & ifs
     ) const noexcept {
      // {//debug
      //   std::cerr << __func__ << std::endl;
      // }

      uint64_t pos = 0;
      for (uint64_t i = 0; i < this->getLenWithEndmarker() - 1; ++i) {
        pos -= (pos > emPos_);
        // std::cout << "T[" << i << "] searchPos(" << pos << ") ";
        const uint64_t idxM = drle_.searchPosM(pos);
        const auto ch = static_cast<unsigned char>(drle_.getCharFromIdxM(idxM));
        // std::cout << "idxM = " << idxM << ", pos = " << pos << ", ch = " << (int)ch << "(" << ch << ")" << std::endl;
        char c; // Assume that the input character fits in char.
        ifs.get(c);
        unsigned char uc = static_cast<unsigned char>(c);
        if (uc != ch) {
          std::cerr << "!error!: bad expansion at i = " << i << ", (" << ch << ") should be (" << uc << ")" << ", idxM = " << idxM << ", pos = " << pos << std::endl;
          return false;
        }
        pos = drle_.rank(ch, idxM, pos, true);
      }
      return true;
    }
  };
};

#endif

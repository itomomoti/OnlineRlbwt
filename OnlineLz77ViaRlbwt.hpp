/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file OnlineLz77ViaRlbwt.hpp
 * @brief Online LZ77 computation via online RLBWT.
 * @author Tomohiro I
 * @date 2017-10-12
 */
#ifndef INCLUDE_GUARD_OnlineLz77ViaRlbwt
#define INCLUDE_GUARD_OnlineLz77ViaRlbwt

#include <stdint.h>
#include <cassert>
#include <iostream>
#include <fstream>

namespace itmmti 
{
  /*!
   * @brief Online LZ77 computation via online RLBWT.
   * @note
   *   ::OnlineLz77ViaRlbwt wraps ::DynRle to use it for representing dynamic RLE of BWT.
   *   In contrast to ::DynRle, OnlineLz77ViaRlbwt has an implicit end marker (em_) at emPos_.
   */
  template<class DynRle>
  class OnlineLz77ViaRlbwt
  {
  public:
    //// bwtintvl: [left..right), where right bound is excluded.
    using bwtintvl = std::pair<uint64_t, uint64_t>;
    //// bwttracker: The first two uints represent a bwt-interval [left..right),
    //// and the last uint is tracking the text-position of BWT[left].
    using bwttracker = std::tuple<uint64_t, uint64_t, uint64_t>;
    using CharT = typename DynRle::CharT;
    using BTreeNodeT = typename DynRle::BTreeNodeT;
    static constexpr uintptr_t NOTFOUND{BTreeNodeT::NOTFOUND};
    static constexpr uint8_t kB{DynRle::kB};
    static constexpr uint8_t kBtmBM{DynRle::kBtmBM};
    static constexpr uint8_t kBtmBS{DynRle::kBtmBS};


  private:
    DynRle drle_; //!< Dynamic RLE data structure with one sampling position per run.
    uint64_t emPos_; //!< Current position (0base) of end marker.
    uint64_t succSamplePos_; //!< Tracking txt-pos (0base) for bwt-position next to current emPos_.
    CharT em_; //!< End marker. It is used only when bwt[emPos_] is accessed (and does not matter if em_ appears in the input text).


  public:
    OnlineLz77ViaRlbwt
    (
     const size_t initNumBtms, //!< Initial size of DynRle to reserve.
     const uint64_t initSampleUb = 256, //!< Initial upper bound (exclusive) of sample positions.
     CharT em = 0 //!< End marker (default 0).
     ) :
      drle_(initNumBtms, initSampleUb),
      emPos_(0),
      succSamplePos_(0),
      em_(em)
    {}


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
     * @brief Get current "succSamplePos_".
     */
    uint64_t getSuccSamplePos() const noexcept {
      return succSamplePos_;
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
     * @brief Extend RLBWT by appending one character while maintaining sample pos for leftend of each run.
     */
    void extend
    (
     const CharT ch //!< Character to append.
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
        idxM = drle_.pushbackRun(pos, ch);
      }
      emPos_ = drle_.rank(ch, idxM, pos, true);

      const uint64_t exponent = drle_.getWeightFromIdxM(idxM);
      if (pos == 0) {
        drle_.setSample(idxM, txtPos);
        if (split) {
          drle_.setSample(drle_.getPrevIdxM(idxM), originalSample);
          drle_.setSample(drle_.getNextIdxM(idxM), succSamplePos_);
        }
      }
      if (pos + 1 != exponent) {
        ++succSamplePos_;
      } else {
        uint64_t idxS = drle_.idxM2S(idxM);
        const uint64_t nextIdxS = drle_.getNextIdxS(idxS);
        if (nextIdxS != BTreeNodeT::NOTFOUND) { // Succcessor with "ch" was found.
          succSamplePos_ = drle_.getSampleFromIdxS(nextIdxS) + 1;
        } else { // Succcessor with "ch" was NOT found.
          //// Take the smallest character larger "next_ch" than "ch".
          //// If such character exists, set "succSamplePos_" to the
          //// sampled position for the first BWT-run of "next_ch" minus one.
          const auto nextRootS = drle_.getNextRootS(drle_.getParentFromBtmS(idxS / kBtmBS));
          if (reinterpret_cast<uintptr_t>(nextRootS) != BTreeNodeT::NOTFOUND) {
            succSamplePos_ = drle_.getSampleFromIdxS(reinterpret_cast<uint64_t>(nextRootS->getLmBtm_DirectJump()) * kBtmBS + 1) + 1;
          }
        }
      }
    }


    /*!
     * @brief Access to the current RLBWT by [] operator.
     */
    CharT operator[]
    (
     uint64_t pos //!< in [0, OnlineRLBWT::getLenWithEndmarker()].
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
     * @brief Compute bwt-interval for cW from bwt-interval for W.
     * @note Intervals are [left..right), where right bound is excluded.
     */
    bool lfMap
    (
     bwttracker & tracker,
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
      if (l_in_drle == emPos_) {
        std::get<2>(tracker) = drle_.getSumOfWeight();
      } else if (ch == chNow) {
        std::get<2>(tracker) += 1;
      } else {
        std::get<2>(tracker) = drle_.getSampleFromIdxS(drle_.getNextIdxS(idxS)) + 1;
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
    uint64_t lfMap(uint64_t i) {
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
      size += drle_.calcMemBytes();
      return size;
    }


    /*!
     * @brief Print statistics
     */
    void printStatistics
    (
     std::ostream & os, //!< std::ostream (e.g., std::cout).
     const bool verbose
     ) const noexcept {
      os << "OnlineLz77ViaRlbwt object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
      os << "Len with endmarker = " << getLenWithEndmarker() << std::endl;
      os << "emPos_ = " << emPos_ << ", succSamplePos_ = " << succSamplePos_ << ", em_ = " << em_ << std::endl;
      drle_.printStatistics(os, verbose);
      os << "OnlineLz77ViaRlbwt object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
    }


    void printDebugInfo
    (
     std::ostream & os //!< std::ostream (e.g., std::cout).
     ) const noexcept {
      os << "Len with endmarker = " << getLenWithEndmarker() << std::endl;
      os << "emPos_ = " << emPos_ << ", succSamplePos_ = " << succSamplePos_ << ", em_ = " << em_ << std::endl;
      drle_.printDebugInfo(os);
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
          std::cerr << "error: bad expansion at i = " << i << ", (" << ch << ") should be (" << uc << ")" << ", idxM = " << idxM << ", pos = " << pos << std::endl;
          return false;
        }
        pos = drle_.rank(ch, idxM, pos, true);
      }
      return true;
    }
  };
};

#endif

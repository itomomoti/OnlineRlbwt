/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file OnlineRLBWT.hpp
 * @brief Online Run-length encoded Burrows–Wheeler transform (RLBWT).
 * @author Tomohiro I
 * @date 2017-01-29
 */
#ifndef INCLUDE_GUARD_OnlineRLBWT
#define INCLUDE_GUARD_OnlineRLBWT

#include <stdint.h>
#include <cassert>
#include <ostream>
#include <fstream>

#include "BTree.hpp"
#include "DynRleWithValue.hpp"

namespace itmmti 
{
  using bwtintvl = std::pair<uint64_t, uint64_t>;
  using bwttracker = std::tuple<uint64_t, uint64_t, uint64_t>;

  /*!
   * @brief Online Run-length encoded Burrows–Wheeler transform (RLBWT).
   * @note
   *   ::OnlineRlbwt wraps ::DynRLE to use it for representing dynamic RLE of BWT.
   *   In contrast to ::DynRle, OnlineRLWT has an implicit end marker (em_) at emPos_.
   */
  template<class DynRle>
  class OnlineRlbwt
  {
  public:
    using BTreeNodeT = typename DynRle::BTreeNodeT;
    static constexpr uint8_t kB{DynRle::kB};
    static constexpr uint8_t kBtmB{DynRle::kBtmB};


  private:
    DynRle drle_;
    uint64_t emPos_; //!< Current position (0base) of end marker.
    uint64_t em_; //!< End marker that should not appear in the input text.
    uint64_t succSamplePos_; // Tracking txt-pos (0base) for bwt-position next to current emPos_.


  public:
    OnlineRlbwt
    (
     const size_t initNumBtms, //!< Give initial size of DynRle to reserve.
     uint64_t em = UINT64_MAX //!< Give end marker (default UINT64_MAX).
     ) :
      drle_(initNumBtms),
      emPos_(0),
      em_(em),
      succSamplePos_(0)
    {}


    /*!
     * @brief Get end marker.
     */
    uint64_t getEndmarker() const noexcept {
      return em_;
    }


    /*!
     * @brief Get current position of end marker.
     */
    uint64_t getEndmarkerPos() const noexcept {
      return emPos_;
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
    uint64_t getLeafVal(uint64_t idxM) const noexcept {
      return drle_.getLeafValFromIdxM(idxM);
    }


    /*!
     * @brief Set associated value at "idxM".
     */
    void setLeafVal(uint64_t val, uint64_t idxM) noexcept {
      drle_.setLeafVal(val, idxM);
    }


    /*!
     * @brief Extend RLBWT by appending one.
     */
    void extend
    (
     const uint64_t ch //!< 64bit-char to append.
     ) {
      const uint64_t txtPos = drle_.getSumOfWeight(); //!< Txt-position of "ch" (0base).
      uint64_t idxM;
      auto pos = emPos_;
      idxM = drle_.insertRun(ch, 1, txtPos, succSamplePos_, pos); // 'pos' is modified to be the relative pos in the run of 'idxM'.
      emPos_ = drle_.rank(ch, idxM, pos, true);

      const uint64_t exponent = drle_.getWeightFromIdxM(idxM);
      if (pos == 0 && exponent > 1) {
        drle_.setLeafVal(idxM, txtPos);
      }
      if (pos + 1 != exponent) {
        ++succSamplePos_;
      } else {
        uint64_t idxS = drle_.idxM2S(idxM);
        const uint64_t nextIdxS = drle_.getNextIdxS(idxS);
        if (nextIdxS != DynRle::BTreeNodeT::NOTFOUND) { // Succcessor with "ch" was found.
          succSamplePos_ = drle_.getLeafValFromIdxS(nextIdxS) + 1;
        } else { // Succcessor with "ch" was NOT found.
          /* Take the smallest character larger "next_ch" than "ch".
             If such character exists, set "succSamplePos_" to the
             sampled position for the first BWT-run of "next_ch" minus one. */
          const auto nextRootS = drle_.getNextRootS(drle_.getParentFromBtmS(idxS / kBtmB));
          if (reinterpret_cast<uintptr_t>(nextRootS) != BTreeNodeT::NOTFOUND) {
            succSamplePos_ = drle_.getLeafValFromBtmNodeS(nextRootS->getLmBtm_DirectJump(), 1) + 1;
          }
        }
      }
    }


    /*!
     * @brief Access to the current RLBWT by [] operator.
     */
    uint64_t operator[]
    (
     uint64_t pos //!< in [0, OnlineRLBWT::getLenWithEndmarker()].
     ) const noexcept {
      assert(pos < getLenWithEndmarker());

      if (pos == emPos_) {
        return em_;
      } else if (pos > emPos_) {
        --pos;
      }
      uint64_t idxM = drle_.searchPosM(pos);
      return drle_.getCharFromIdxM(idxM);
    }


    /*!
     * @brief Return current length including end marker.
     */
    uint64_t getLenWithEndmarker() const noexcept {
      return drle_.getSumOfWeight() + 1; // +1 for end marker, which is not in drle_.
    }


    /*!
     * @brief Return 'rank of ch at pos' + 'num of total occ of characters smaller than ch'.
     */
    uint64_t totalRank
    (
     uint64_t ch,
     uint64_t pos //!< in [0, OnlineRLBWT::getLenWithEndmarker()].
     ) const noexcept {
      assert(pos < getLenWithEndmarker());

      if (pos > emPos_) {
        --pos;
      }
      return drle_.rank(ch, pos, true);
    }


    /*!
     * @brief Compute bwt-interval for cW from bwt-interval for W
     * @note Intervals are [left, right) : right bound is excluded
     */
    bool lfMap
    (
     bwttracker & tracker,
     const uint64_t ch
     ) const noexcept {
      assert(ch != getEndmarker());
      assert(std::get<0>(tracker) <= getLenWithEndmarker() && std::get<1>(tracker) <= getLenWithEndmarker());
      assert(std::get<0>(tracker) < std::get<1>(tracker));

      // If "ch" is not in the alphabet or empty interval, return empty interval.
      const auto * retRootS = drle_.searchCharA(ch);
      if (retRootS->isDummy() || drle_.getCharFromNodeS(retRootS) != ch) {
        return false;
      }

      uint64_t r_in_drle = std::get<1>(tracker) - (std::get<1>(tracker) > emPos_); // Taking (implicit) end-marker into account.
      /* +1 because in F we are not taking into account the end-marker,
         which is in position 0 but not explicitly stored in F. */
      r_in_drle = drle_.rank(ch, r_in_drle - 1, true) + 1;
      if (r_in_drle <= 1) {
        return false;
      }

      uint64_t l_in_drle = std::get<0>(tracker) - (std::get<0>(tracker) > emPos_); // Taking (implicit) end-marker into account.
      const uint64_t idxM = drle_.searchPosM(l_in_drle); // l_in_drle is modified to relative pos.
      /*
       * In order to get idxS, replicate variant of rank function,
       * where pos is specified by 'idxM' and 'relativePos' with several modifications.
       */
      const auto chNow = drle_.getCharFromIdxM(idxM);
      uint64_t idxS;
      {
        if (ch == chNow) {
          idxS = drle_.idxM2S(idxM);
        } else {
          l_in_drle = 0;
          idxS = drle_.getPredIdxSFromIdxM(retRootS, ch, idxM); // We know that "ch" already exists by "r_in_dre > 1".
        }
        const auto btmS = idxS / kBtmB;
        for (auto tmpIdxS = btmS * kBtmB; tmpIdxS < idxS + (ch != chNow); ++tmpIdxS) {
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

      // Update tracker.
      std::get<0>(tracker) = l_in_drle;
      std::get<1>(tracker) = r_in_drle;
      if (l_in_drle == emPos_) {
        std::get<2>(tracker) = drle_.getSumOfWeight();
      } else if (ch == chNow) {
        std::get<2>(tracker) += 1;
      } else {
        std::get<2>(tracker) = drle_.getLeafValFromIdxS(drle_.getNextIdxS(idxS)) + 1;
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
     const uint64_t ch
     ) const noexcept {
      assert(ch != getEndmarker());
      assert(intvl.first <= getLenWithEndmarker() && intvl.second <= getLenWithEndmarker());

      // If "ch" is not in the alphabet or empty interval, return empty interval.
      const auto * retRootS = drle_.searchCharA(ch);
      if (retRootS->isDummy() || drle_.getCharFromNodeS(retRootS) != ch || intvl.first >= intvl.second) {
        return {0, 0};
      }

      uint64_t l = intvl.first - (intvl.first > emPos_);
      uint64_t r = intvl.second - (intvl.second > emPos_);
      const uint64_t idxM = drle_.searchPosM(l); // l is modified to relative pos.
      // +1 because in F.select(0, ch) we are not taking into account the end-marker,
      // which is in position 0 but not explicitly stored in F.
      return {
        drle_.rank(ch, idxM, l, true) - (drle_.getCharFromIdxM(idxM) == ch) + 1,
          drle_.rank(ch, r-1, true) + 1
          };
    }


    /*!
     * @brief LF map.
     */
    uint64_t lfMap(uint64_t i){
      assert(i < getLenWithEndmarker());

      if (i > emPos_) {
        --i;
      }
      const uint64_t idxM = drle_.searchPosM(i);
      const unsigned char ch = drle_.getCharFromIdxM(idxM);
      return drle_.rank(ch, idxM, i, true);
    }


    /*!
     * @brief Print statistics of ::DynRLE (not of ::OnlineRLBWT).
     */
    void printStatictics
    (
     std::ostream & os //!< std::ostream (e.g., std::cout).
     ) const noexcept {
      drle_.printStatictics(os);
    }


    void printDebugInfo(std::ostream & os) const noexcept {
      os <<  "emPos_ = " << emPos_ << ", em_ = " << em_ << ", succSamplePos_ = " << succSamplePos_ << std::endl;
      drle_.printDebugInfo(os);
    }


    /*!
     * @brief Calculate total memory usage in bytes.
     */
    size_t calcMemBytes() const noexcept {
      return sizeof(*this) + drle_.calcMemBytes();
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
        if (pos > emPos_) {
          --pos;
        }
        const uint64_t idxM = drle_.searchPosM(pos);
        const unsigned char ch = drle_.getCharFromIdxM(idxM);
        ofs.put(ch);
        pos = drle_.rank(ch, idxM, pos, true);
      }
    }


    void checkDecompress
    (
     std::ifstream & ifs
     ) const noexcept {
      uint64_t pos = 0;
      for (uint64_t i = 0; i < this->getLenWithEndmarker() - 1; ++i) {
        if (pos > emPos_) {
          --pos;
        }
        // std::cout << "T[" << i << "] searchPos(" << pos << ") ";
        const uint64_t idxM = drle_.searchPosM(pos);
        const unsigned char ch = drle_.getCharFromIdxM(idxM);
        // std::cout << "idxM = " << idxM << ", pos = " << pos << ", ch = " << (int)ch << "(" << ch << ")" << std::endl;
        char c; // Assume that the input character fits in char.
        ifs.get(c);
        unsigned char uc = static_cast<unsigned char>(c);
        if (uc != ch) {
          std::cout << "error: bad expansion at i = " << i << ", (" << ch << ") should be (" << uc << ")" << ", idxM = " << idxM << ", pos = " << pos << std::endl;
          return;
        }
        pos = drle_.rank(ch, idxM, pos, true);
      }
      std::cout << "correctly decompressed" << std::endl;
    }
  };
};

#endif

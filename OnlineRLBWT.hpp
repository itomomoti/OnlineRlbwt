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
#include <ostream>
#include <fstream>

using bwtintvl = std::pair<uint64_t, uint64_t>;

/*!
 * @brief Online Run-length encoded Burrows–Wheeler transform (RLBWT).
 * @note
 *   ::OnlineRLBWT wraps ::DynRLE to use it for representing dynamic RLE of BWT.
 *   In contrast to ::DynRLE, OnlineRLWT has a vertial end marker (em_) at emPos_.
 */
template <class DynRLE>
class OnlineRLBWT
{
  DynRLE drle_;
  uint64_t emPos_; //!< Current position (0base) of end marker.
  uint64_t em_; //!< End marker that should not appear in the input text.


public:
  OnlineRLBWT
  (
   const size_t initNumBtms, //!< Give initial size of DynRLE to reserve.
   uint64_t em = UINT64_MAX //!< Give end marker (default UINT64_MAX).
   ) :
    drle_(initNumBtms),
    emPos_(0),
    em_(em)
  {}


  /*!
   * @brief Get end marker.
   */
  uint64_t getEm() const noexcept {
    return em_;
  }


  /*!
   * @brief Get current position of end marker.
   */
  uint64_t getEmPos() const noexcept {
    return emPos_;
  }


  /*!
   * @brief Extend RLBWT by appending one.
   */
  void extend
  (
   uint64_t ch //!< 64bit-char to append.
   ) {
    uint64_t idxM = drle_.insertRun(ch, 1, emPos_);
    emPos_ = drle_.rank(ch, idxM, emPos_, true);
  }


  /*!
   * @brief Access to the current RLBWT by [] operator.
   */
  uint64_t operator[]
  (
   uint64_t pos //!< in [0, OnlineRLBWT::getLenWithEm()].
   ) const noexcept {
    assert(pos < getLenWithEm());

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
  uint64_t getLenWithEm() const noexcept {
    return drle_.getSumOfWeight() + 1; // +1 for end marker, which is not in drle_.
  }


  /*!
   * @brief Return 'rank of ch at pos' + 'num of total occ of characters smaller than ch'.
   */
  uint64_t totalRank
  (
   uint64_t ch,
   uint64_t pos //!< in [0, OnlineRLBWT::getLenWithEm()].
   ) const noexcept {
    assert(pos < getLenWithEm());

    if (pos > emPos_) {
      --pos;
    }
    return drle_.rank(ch, pos, true);
  }


  /*!
   * @brief Compute bwt-interval for cW from bwt-interval for W
   * @note Intervals are [left, right) : right bound is excluded
   */
  bwtintvl lfMap
  (
   bwtintvl intvl,
   uint64_t ch
   ) const noexcept {
    assert(ch != getEm());
    assert(intvl.first <= getLenWithEm() && intvl.second <= getLenWithEm());

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
    assert(i < getLenWithEm());

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
    for (uint64_t i = 0; i < this->getLenWithEm() - 1; ++i) {
      if (pos > emPos_) {
        --pos;
      }
      const uint64_t idxM = drle_.searchPosM(pos);
      const unsigned char ch = drle_.getCharFromIdxM(idxM);
      ofs.put(ch);
      pos = drle_.rank(ch, idxM, pos, true);
    }
  }
};

#endif

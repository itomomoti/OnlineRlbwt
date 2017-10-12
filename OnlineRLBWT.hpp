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
#ifndef INCLUDE_GUARD_DynRLBWT
#define INCLUDE_GUARD_DynRLBWT

#include <stdint.h>
#include <ostream>
#include <fstream>


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

    if (pos == emPos_ && ch == em_) {
      pos = 0;
      return 1;
    }
    if (pos > emPos_) {
      --pos;
    }
    return drle_.rank(ch, pos, true);
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


  // /*!
  //  * @brief Output string represented by current RLE to std::ofstream.
  //  */
  // void printString
  // (
  //  std::ofstream & ofs
  //  ) const noexcept {
  //   drle_.printString(ofs);
  // }
};

#endif

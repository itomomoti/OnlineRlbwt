/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file TagRelabelAlgo.hpp
 * @brief Utilities to implement linked-lists that can answer order queries in O(1) time using Tag-range Relabeling Algorithm (TRA).
 * @author Tomohiro I
 * @date 2017-04-08
 */
#ifndef INCLUDE_GUARD_TagRelabelAlgo
#define INCLUDE_GUARD_TagRelabelAlgo

#include <stdint.h>
#include <cassert>

#include "BitsUtil.hpp"

/*!
 * @brief Utilities to implement linked-lists that can answer order queries in O(1) time using Tag-range Relabeling Algorithm (TRA).
 * @par Reference
 *   [1] Bender, M. A.; Cole, R.; Demaine, E. D.; Farach-Colton, M. & Zito, J. Two Simplified Algorithms for Maintaining Order in a List, ESA, 2002, 152-164.
 * @par Notation
 *   A traCode is uint ranging from 9 to 15, which determines the threshold of overflow (it corresponds to parameter T of [1] but not exactly the same).
 *   More precisely, a tag-range sharing '64 - x' most significant bits is said to be overflow iff the # of elements in the tag-range is larger than '(traCode/8)^x'.
 */
struct TagRelabelAlgo
{
  TagRelabelAlgo() = delete;

  static constexpr uint8_t NUM_TRACODE{7};
  static constexpr uint64_t MAX_LABEL{ctcbits::UINTW_MAX(63)};
  /*!
   * @brief TBL_Capacities[i] represents the max # of the elements handled with traCode of 'i + 9'.
   * @note
   *   TBL_Capacities[i] is roughly '((i+9)/8)^63' (not exactly so due to rounding up for each multiplication of '((i+9)/8)').
   */
  static const uint64_t TBL_Capacities[NUM_TRACODE];


  /*!
   * @brief Get the smallest traCode to handle 'num' elements.
   */
  static uint8_t getSmallestTraCode(uint64_t num) noexcept {
    assert(TBL_Capacities[NUM_TRACODE-1] >= num);

    uint8_t i = 0;
    while (TBL_Capacities[i] < num) {
      ++i;
    }
    return i + 9;
  }

  
  /*!
   * @brief Multiply '(traCode/8)' to 'overflowNum'.
   */
  static uint64_t getNextOverflowNum(uint64_t overflowNum, uint8_t traCode) noexcept {
    return ((overflowNum * traCode) >> 3) + 1;
  }
};


#endif

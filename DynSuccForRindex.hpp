/*!
 * Copyright (c) 2018 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file DynSuccWithVal.hpp
 * @brief Dynamic succsessor data structure (supporting indel) used in r-index.
 * @note Implementation by wrapping PSumWithValue.
 * @author Tomohiro I
 * @date 2018-05-02
 */
#ifndef INCLUDE_GUARD_DynSuccForRindex
#define INCLUDE_GUARD_DynSuccForRindex

#include <stdint.h>

#include <iostream>
#include <algorithm>
#include <cassert>
#include <string>
#include <sstream>

// include from Basics
#include "BitsUtil.hpp"
#include "StepCode.hpp"
#include "MemUtil.hpp"
//

#include "PSumWithValue.hpp"


namespace itmmti
{
  /*!
   * @brief
   * @tparam tparam_BTreeNodeT Use ::BTreeNode<kB> for now.
   * @tparam tparam_BtmNodeT Use ::BtmNodeForPSumWithVal<kBtmB> for now.
   */
  template<typename tparam_BTreeNodeT, typename tparam_BtmNodeT>
  class DynSuccForRindex
  {
  public:
    //// Public constant, alias etc.
    using BTreeNodeT = tparam_BTreeNodeT;
    using BtmNodeT = tparam_BtmNodeT;
    using PSumT = PSumWithValue<BTreeNodeT, BtmNodeT>;
    static constexpr uint8_t kB{BTreeNodeT::kB};
    static constexpr uint8_t kBtmB{BtmNodeT::kBtmB};
    // //// Key (uint stored in successor data structure) and associated value.
    // using KeyVal = std::pair<uint64_t, uint64_t>;


  private:
    //// Private member variables.
    PSumT psum_; //!< Dynamic PSum data structure.


  public:
    DynSuccForRindex()
      : psum_()
    {
      init();
    }


    ~DynSuccForRindex()
    {
      // psum_ is deleted.
    }


    /*!
     * @brief Initialize data structure with empty elements.
     */
    void init() {
      if (isReady()) {
        clearAll();
      }
      psum_.init();
      // Insert a sentinel having UINT64_MAX weight at the end.
      const uint64_t dummyArray[] = {UINT64_MAX};
      BTreeNodeT * parent;
      uint8_t idxInSib;
      auto btm = psum_.getRmBtm(parent, idxInSib);
      btm->insert(dummyArray, dummyArray, 1, 1, 0, parent, idxInSib);
    }


    /*!
     * @brief Free/delete all allocated objects.
     */
    void clearAll() {
      if (!isReady()) { // already cleared
        return;
      }
      psum_.clearAll();
    }


    /*!
     * @brief Return if data structure is ready.
     */
    bool isReady() const noexcept {
      return psum_.isReady();
    }


    /*!
     * @brief Function to get text-position for "BWT[bwtpos + 1]", where "BWT[bwtpos]" corresponds to "T[txtpos]".
     */
    // uint64_t calcNextPos
    // (
    //  uint64_t txtPos //!< Text-position for currently focused character.
    //  ) const noexcept {
    //   std::cerr << __func__ << ": txtPos = " << txtPos << std::endl;

    //   BTreeNodeT * parent;
    //   uint8_t idxInSib;
    //   auto btm = psum_.searchBtm(txtPos, parent, idxInSib);

    //   uint64_t retWeight, retVal;
    //   btm->searchPos(txtPos, retWeight, retVal);
    //   return retVal - (retWeight - txtPos - 1);
    // }


    /*!
     * @brief Function to get text-position for "BWT[bwtpos + 1]", where "BWT[bwtpos]" corresponds to "T[txtpos]".
     */
    uint64_t calcNextPos
    (
     const uint64_t txtPos, //!< Text-position for currently focused character.
     const uint64_t txtLen, //!< Text length without end marker.
     const uint64_t emPrev, //!< Text-position for previous character of implicit end marker.
     const uint64_t emNext //!< Text-position for next character of implicit end marker.
     ) const noexcept {
      // std::cerr << __func__ << ": txtPos = " << txtPos << ", txtLen = " << txtLen
      //           << ", emPrev = " << emPrev << ", emNext = " << emNext << std::endl;
      assert(txtPos <= txtLen);

      uint64_t q = txtPos;
      BTreeNodeT * parent;
      uint8_t idxInSib;
      auto btm = psum_.searchBtm(q, parent, idxInSib);

      uint64_t retWeight, retVal;
      btm->searchPos(q, retWeight, retVal);
      const uint64_t dist = retWeight - q - 1; // distance to sampled position
      if (txtPos <= emPrev && emPrev - txtPos <= dist) {
        return txtLen - (emPrev - txtPos);
      } else if (txtLen - txtPos <= dist) {
        return emNext - (txtLen - txtPos);
      } else {
        return retVal - dist;
      }
    }


    /*!
     * @brief Register (key, val) in successor data structure.
     */
    void setKeyVal
    (
     const uint64_t key, //!< Key
     const uint64_t val //!< Value
     ) noexcept {
      // std::cerr << __func__ << "(" << this << "): key = " << key << ", val = " << val << std::endl;

      uint64_t q = key;
      BTreeNodeT * parent;
      uint8_t idxInSib;
      auto btm = psum_.searchBtm(q, parent, idxInSib); // q is modified.

      uint64_t retWeight, retVal;
      uint16_t childIdx = btm->searchPos(q, retWeight, retVal); // q is modified.
      const auto fstHalf = q + 1;
      if (fstHalf == retWeight) {
        const uint64_t array[] = {val};
        btm->replace(array, 1, 2 * childIdx + 1, parent, idxInSib);
      } else {
        const uint64_t weights[] = {fstHalf, retWeight - fstHalf};
        const uint64_t vals[] = {val, retVal};
        btm->insert(weights, vals, 2, childIdx, 1, parent, idxInSib);
      }
    }


    /*!
     * @brief Remove key from successor data structure.
     * @note key must exist.
     */
    void removeKey
    (
     const uint64_t key //!< Key
     ) noexcept {
      // std::cerr << __func__ << "(" << this << "): key = " << key << std::endl;

      uint64_t q = key;
      BTreeNodeT * parent;
      uint8_t idxInSib;
      auto btm = psum_.searchBtm(q, parent, idxInSib); // q is modified

      uint64_t retWeight, retVal, bitPos;
      uint16_t childIdx = btm->searchPos(q, retWeight, retVal, bitPos); // q is modified
      assert(q + 1 == retWeight); // Assume that key exists.
      if (childIdx + 1 < btm->getNumChildren()) {
        const uint16_t i = 2 * (childIdx + 1);
        const uint64_t weights[] = {retWeight + btm->readNext(i, bitPos)};
        const uint64_t vals[] = {btm->readNext(i + 1, bitPos)};
        btm->insert_shrink(weights, vals, 1, childIdx, 2, parent, idxInSib);
      } else {
        uint64_t weights[] = {retWeight};
        btm->insert_shrink(weights, weights, 0, childIdx, 1, parent, idxInSib); // Nothing is inserted, weights is used as dummy array.
        //// search again to get element that is next to deleted element.
        q = key - retWeight + 1;
        btm = psum_.searchBtm(q, parent, idxInSib); // q is modified
        childIdx = btm->searchPos(q, retWeight); // q is modified
        weights[0] += retWeight;
        btm->replace(weights, 1, 2 * childIdx, parent, idxInSib);
      }
    }


  public:
    //// statistics
    size_t calcMemBytes
    (
     bool includeThis = true
     ) const noexcept {
      return psum_.calcMemBytes(includeThis);
    }


    void printStatistics
    (
     std::ostream & os,
     const bool verbose
     ) const noexcept {
      os << "DynSuccForRindex object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
      psum_.printStatistics(os, verbose);
      os << "DynSuccForRindex object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
    }


    void printDebugInfo
    (
     std::ostream & os
     ) const noexcept {
      os << __func__ << " DynSuccForRindex: kB = " << (int)kB << ", kBtmB = " << (int)kBtmB << std::endl;
      psum_.printDebugInfo(os);
    }
  };
} // namespace itmmti

#endif

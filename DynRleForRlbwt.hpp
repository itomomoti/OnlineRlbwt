/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file DynRleForRlbwt.hpp
 * @brief Dynamic Run-length encoding with value.
 * @author Tomohiro I
 * @date 2017-12-16
 */
#ifndef INCLUDE_GUARD_DynRleForRlbwt
#define INCLUDE_GUARD_DynRleForRlbwt

#include <stdint.h>

#include <iostream>
#include <algorithm>
#include <cassert>
#include <string>
#include <fstream>
#include <sstream>

// include from Basics
#include "BitsUtil.hpp"
#include "StepCode.hpp"
#include "MemUtil.hpp"

// include from BTree
#include "BTree.hpp"
#include "TagRelabelAlgo.hpp"


namespace itmmti
{
  template<uint8_t kB> class BTreeNode;


  /*!
   * @brief Dynamic partial sum data structure implemented by B+tree, where each leaf is associated with value.
   * @tparam kB Arity for internal node of B+tree, which should be in {32, 64, 128}.
   * @tparam kBtmB Arity for bottom node of B+tree, which should be in {32, 64, 128}.
   * @par Notation
   *   - T: Current string represented by RLE.
   *   - Mixed tree: B+tree representing RLE of T.
   *     - btmM: Index of bottom node of B+tree on the array-based implementation.
   *             Each bottom node "btmM" can have "kBtmB" children, which correspond to indexes [btmM * kBtmB, (btmM+1) * kBtmB).
   *     - idxM: Indexes that are corresponding to children of "btmM".
   *   - Separated tree: B+tree separately representing runs for each character.
   *     - btmS: Index of bottom node of B+tree on the array-based implementation (all separated trees share arrays).
   *             Each bottom node "btmS" can have "kBtmB" children, which correspond to indexes [btmS * kBtmB, (btmS+1) * kBtmB).
   *     - idxS: Indexes that are corresponding to children of "btmS".
   */
  template<uint8_t param_kB, uint8_t param_kBtmB>
  class DynRleForRlbwt
  {
  public:
    // Public constant, alias etc.
    static constexpr uint8_t kB{param_kB};
    static constexpr uint8_t kBtmB{param_kBtmB};
    using BTreeNodeT = BTreeNode<kB>;

    static constexpr size_t kUnitBits{kBtmB * 8};
    static constexpr bool kM{0};
    static constexpr bool kS{1};


  private:
    // Inner class
    class BtmNodeM
    {
    public:
      //// Private member variables.
      StepCodeCore<kBtmB> stcc_; //!< Storing weights
      uint64_t btmVal_; //!< TRA label of btmM.
      BTreeNodeT * parent_; //!< Pointer to parent of btmM.
      uint16_t stccCapacity_; //!< Current bit capacity of stcc_
      uint16_t stccSize_; //!< Current bit size of stcc_
      uint8_t numChildren_; //!< Current size (number of elements).
      uint8_t idxInSibling_; //!< idxInSibling of btmM.
      uint8_t wCodesAuxM_[kBtmB / StepCodeUtil::kWCNum - 1];


      BtmNodeM
      (
       uint16_t initStccCapacity = 0
       ) :
        stcc_(),
        btmVal_(0),
        parent_(nullptr),
        stccCapacity_(0),
        stccSize_(0),
        numChildren_(0),
        idxInSibling_(0)
      {
        assert(initStccCapacity <= UINT16_MAX - 64);

        reserveBitCapacity(initStccCapacity);
      }


      ~BtmNodeM() {
        // stcc_ is freed.
      }


      uint16_t getStccSize() const noexcept {
        return stccSize_;
      }


      uint16_t getStccCapacity() const noexcept {
        return stccCapacity_;
      }


      uint16_t getNumChildren() const noexcept {
        return numChildren_;
      }


      uint16_t getIdxInSibling() const noexcept {
        return idxInSibling_;
      }


      uint64_t getBtmVal() const noexcept {
        return btmVal_;
      }


      void setBtmVal
      (
       uint64_t val
       ) noexcept {
        btmVal_ = val;
      }


      BTreeNodeT * getParent() const noexcept {
        return parent_;
      }


      uint64_t readStccVal
      (
       const uint8_t idx //!< in [0..numChildren_)
       ) const noexcept {
        assert(idx < numChildren_);

        return stcc_.read(idx);
      }


      uint64_t calcSumOfWeight() const noexcept {
        const auto & stcc = this->getConstRef_stcc();
        const auto num = numChildren_;
        uint64_t sum = 0;
        uint64_t bitPos = 0;
        for (uint64_t i = 0; i < num; ++i) {
          const auto w = stcc.readW(i);
          sum += stcc.readWBits(bitPos, w);
          bitPos += w;
        }

        return sum;
      }


      /*!
       * @brief Calculate the beginning bit-pos of "idx"-th value in stcc_.
       */
      uint16_t calcBitPos
      (
       const uint8_t idx //!< in [0..numChildren_]
       ) const noexcept {
        assert(idx <= numChildren_);

        if (idx < numChildren_) {
          return static_cast<uint16_t>(stcc_.calcBitPos(idx, wCodesAuxM_));
        } else {
          return stccSize_;
        }
      }


      /*!
       * @brief Get read-only array pointer.
       */
      const StepCodeCore<kBtmB> & getConstRef_stcc() const noexcept
      {
        return stcc_;
      }


      /*!
       * @brief Get read-only array pointer.
       */
      const uint64_t * getConstPtr_vals() const noexcept
      {
        return stcc_.getConstPtr_vals();
      }


      /*!
       * @brief Get read-only wCodes_ array pointer.
       */
      const uint64_t * getConstPtr_wCodes() const noexcept
      {
        return stcc_.getConstPtr_wCodes();
      }


      void reserveBitCapacity
      (
       uint16_t givenBitCapacity
       ) {
        if (givenBitCapacity > this->stccCapacity_) {
          size_t newSize = (static_cast<size_t>(givenBitCapacity) / kUnitBits + 2) * kUnitBits;
          this->stccCapacity_ = static_cast<uint16_t>(this->stcc_.setBitCapacity(static_cast<size_t>(givenBitCapacity)));
        }
      }


      void shrinkBitCapacity() {
        if (this->stccCapacity_ - this->stccSize_ > kUnitBits) {
          this->stccCapacity_ = static_cast<uint16_t>(stcc_.setBitCapacity(static_cast<size_t>(this->stccSize_)));
        }
      }


      void setParentRef
      (
       BTreeNodeT * newParent,
       uint8_t newIdxInSibling
       ) noexcept {
        this->parent_ = newParent;
        this->idxInSibling_ = newIdxInSibling;
      }


      /*!
       * @brief Resize "numChildren_" to "newSize".
       * @note
       *   It does not change bitCapacity.
       */
      void resize
      (
       const uint8_t newSize
       ) noexcept {
        assert(newSize <= kBtmB);

        numChildren_ = newSize;
      }


      /*!
       * @brief update wCodesAuxM.
       */
      void updateWCodesAuxM
      (
       const uint16_t idxBeg,
       const uint16_t idxEnd
       ) noexcept {
        assert(idxBeg < idxEnd);

        const uint64_t beg = idxBeg / StepCodeUtil::kWCNum;
        const uint64_t end = (idxEnd - 1) / StepCodeUtil::kWCNum + 1;
        // std::cerr << __FUNCTION__ << ": " << idxBeg << "->" << beg << ", " << idxEnd << "->" << end << std::endl;
        stcc_.updateWCodesAuxM(wCodesAuxM_, beg, end);
      }


      /*!
       * @brief Replace values.
       */
      void replace
      (
       const uint64_t * newVals, //!< Storing stcc values that replace existing stcc values
       const uint8_t num, //!< Number of elements to replace.
       const uint8_t idx //!< in [0..numChildren_). Beginning idx of tgt.
       ) {
        assert(idx + num <= numChildren_);

        const uint16_t bitPos0 = this->calcBitPos(idx);
        uint16_t bitPos = bitPos0;
        uint16_t sumW_ins = 0;
        uint16_t sumW_del = 0;
        for (uint16_t i = idx; i < idx + num; ++i) {
          const uint8_t w_old = stcc_.readW(i);
          const uint8_t w_new = StepCodeUtil::calcSteppedW(newVals[i - idx]);
          sumW_del += w_old;
          sumW_ins += w_new;
          stcc_.writeWCode(StepCodeUtil::calcWCodeFromSteppedW(w_new), i);
          bitPos += w_old;
        }
        this->updateWCodesAuxM(idx, idx + num);

        if (sumW_ins != sumW_del) {
          const uint16_t newStccSize = this->stccSize_ + sumW_ins - sumW_del;
          if (newStccSize != this->stccSize_) {
            this->reserveBitCapacity(newStccSize);
            this->stcc_.mvVals(this->stcc_.getConstPtr_vals(), bitPos, bitPos0 + sumW_ins, this->stccSize_ - bitPos);
          }
          this->stccSize_ = newStccSize;
        }
        bitPos = bitPos0;
        for (uint16_t i = idx; i < idx + num; ++i) {
          uint8_t w = this->stcc_.readW(i);
          this->stcc_.writeWBits(newVals[i - idx], bitPos, w);
          bitPos += w;
        }
      }


      //////////////////////////////// statistics
      size_t calcMemBytes
      (
       bool includeThis = true
       ) const noexcept {
        size_t size = sizeof(*this) * includeThis;
        return size + calcMemBytesStccDynArray();
      }


      size_t calcMemBytesStccDynArray() const noexcept {
        return stccCapacity_ / 8;
      }


      void printStatistics
      (
       std::ostream & os,
       const bool verbose
       ) const noexcept {
        os << "DynRleWithValue::BtmNode object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
        os << "BTree arity for bottom node = " << static_cast<int>(kBtmB) << ", btmVal = " << btmVal_ << std::endl;
        os << "parent = " << parent_ << ", idxInSibling = " << (int)idxInSibling_ << ", numChildren = " << static_cast<uint64_t>(numChildren_) << std::endl;
        os << "bit size = " << stccSize_ << ", bit capacity = " << stccCapacity_ << std::endl;
        os << "Total: " << calcMemBytes() << " bytes" << std::endl;
        os << "dynamic array of step code: " << calcMemBytesStccDynArray() << " bytes" << std::endl;
        os << "DynRleWithValue::BtmNode object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
      }


      void printDebugInfo
      (
       std::ostream & os,
       const bool verbose
       ) const noexcept {
        os << "DynRleWithValue::BtmNodeM object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
        os << "BTree arity for bottom node = " << static_cast<int>(kBtmB);
        os << ", SumOfWeight = " << this->calcSumOfWeight() << ", label = " << btmVal_ << std::endl;
        os << "parent = " << parent_ << ", idxInSibling = " << (int)idxInSibling_ << ", numChildren = " << static_cast<uint64_t>(numChildren_) << std::endl;
        os << "bit size = " << stccSize_ << ", bit capacity = " << stccCapacity_ << std::endl;
        os << "Total: " << calcMemBytes() << " bytes, dynamic array of step code: " << calcMemBytesStccDynArray() << " bytes" << std::endl;
        {
          os << "dump bit witdth stored in wCodes (" << stcc_.getConstPtr_wCodes() << ")" << std::endl;
          for (uint8_t i = 0; i < numChildren_; ++i) {
            os << (uint64_t)(stcc_.readW(i)) << " ";
          }
          os << std::endl;
        }
        {
          os << "dump values" << std::endl;
          for (uint8_t i = 0; i < numChildren_; ++i) {
            os << stcc_.read(i) << " ";
          }
          os << std::endl;
        }
        {
          os << "dump bits in vals_ (" << stcc_.getConstPtr_vals() << ")" << std::endl;
          for (uint64_t i = 0; i < (stccSize_ + 63) / 64; ++i) {
            os << "(" << i << ")";
            for (uint64_t j = 0; j < 64; ++j) {
              os << bits::readWBits_S(stcc_.getConstPtr_vals(), 64 * i + 63 - j, ctcbits::UINTW_MAX(1));
            }
            os << " ";
          }
          os << std::endl;
        }
        os << "DynRleWithValue::BtmNodeM object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
      }
    }; // end of BtmNode


  private:
    WBitsVec sample_; //!< Sampled positions.
    typename BTreeNodeT::SuperRootT srootM_; //!< Super root of mixed tree
    // Information for leaves and elements for mixed tree.
    WBitsVec idxM2S_; //!< Packed array mapping idxM to corresponding idxS.
    BtmNodeM ** btmPtrsM_; //!< Storing weights
    // Alphabet tree: the bottoms of alphabet tree are roots of separated trees.
    typename BTreeNodeT::SuperRootT srootA_; //!< Super root of alphabet tree
    // Information for leaves and elements for separated tree.
    WBitsVec idxS2M_; //!< Packed array mapping idxS to corresponding idxM.
    BTreeNodeT ** parentS_; //!< Pointer to parent of btmS.
    uint64_t * charS_; //!< 64bit-char of btmS.
    uint8_t * idxInSiblingS_; //!< idxInSibling of btmS.
    uint8_t * numChildrenS_; //!< Num of children of btmS.

    uint8_t traCode_; //!< traCode in [9..16).
    uint64_t sampleUb_; //!< current sample upper bound. Sample pos (0base) must be less than "sampleUb_". If 0, it means samples are not needed.


  public:
    DynRleForRlbwt
    (
     const size_t initNumBtms = 0,
     const uint64_t initSampleUb = 0
     ) :
      sample_(),
      srootM_(),
      idxM2S_(),
      btmPtrsM_(nullptr),
      srootA_(),
      idxS2M_(),
      parentS_(nullptr),
      charS_(nullptr),
      idxInSiblingS_(nullptr),
      numChildrenS_(nullptr),
      traCode_(9),
      sampleUb_(0)
    {
      if (initNumBtms) {
        init(initNumBtms, initSampleUb);
      }
    }


    ~DynRleForRlbwt() {
      clearAll();
    }


    /*!
     * @brief Reserve space to accomodate 'initNumBtms' bottoms, and init.
     */
    void init
    (
     const size_t initNumBtms,
     const uint64_t initSampleUb
     ) {
      // {//debug
      //   std::cerr << __func__ << "initNumBtms = " << initNumBtms << ", initSampleUb = " << initSampleUb << std::endl;
      // }
      assert(initNumBtms > 0);

      if (isReady()) {
        clearAll();
      }
      if (initSampleUb) {
        increaseWOfSample(initSampleUb - 1);
      }
      reserveBtm(initNumBtms);

      srootM_.setRoot(new BTreeNodeT(reinterpret_cast<void *>(0), true, true, true, true));
      srootM_.root_->putFirstBtm(reinterpret_cast<void *>(0), 0);
      // sentinel
      auto firstBtmNodeM = new BtmNodeM(kUnitBits);
      firstBtmNodeM->setParentRef(srootM_.root_, 0);
      setNewBtmNodeM(firstBtmNodeM);
      idxM2S_.resize(kBtmB);
      sample_.resize(kBtmB);
      {
        const uint64_t newVals[] = {0};
        const uint64_t newLinks[] = {0};
        insertNewElemM_simple(0, 0, newVals, newLinks, 1, 0);
      }

      auto * dummyRootS = new BTreeNodeT(nullptr, true, true, true, false, true);
      dummyRootS->putFirstBtm(nullptr, 0);
      srootA_.setRoot(new BTreeNodeT(dummyRootS, true, true, true, true));
      srootA_.root_->pushbackBTreeNode(dummyRootS);
    }


    /*!
     * @brief Free/delete all allocated objects.
     */
    void clearAll() {
      if (!isReady()) { // already cleared
        return;
      }

      for (uint64_t i = 0; i < idxM2S_.size() / kBtmB; ++i) {
        delete btmPtrsM_[i];
      }
      if (sampleUb_) {
        sample_.clear(); 
        sample_.convert(1, 0, true);
        sampleUb_ = 256;
      }
      idxM2S_.clear(); 
      idxM2S_.changeCapacity(0);
      idxS2M_.clear(); 
      idxS2M_.changeCapacity(0);
      memutil::safefree(btmPtrsM_);
      memutil::safefree(parentS_);
      memutil::safefree(charS_);
      memutil::safefree(idxInSiblingS_);
      memutil::safefree(numChildrenS_);

      { // delete separated tree
        auto * rootS = srootA_.root_->getLmBtm_DirectJump();
        while (reinterpret_cast<uintptr_t>(rootS) != BTreeNodeT::NOTFOUND) {
          auto * next = getNextRootS(rootS);
          delete rootS;
          rootS = next;
        }
      }
      memutil::safedelete(srootM_.root_);
      memutil::safedelete(srootA_.root_);

      traCode_ = 9;
    }


    void increaseWOfSample
    (
     const uint64_t minSupportSampleVal
     ) {
      // {// debug
      //   std::cerr << __func__ << ": minSupportSampleVal = " << minSupportSampleVal << std::endl;
      // }
      if (minSupportSampleVal >= sampleUb_) {
        const uint8_t newW = bits::bitSize(minSupportSampleVal);
        sample_.convert(newW, idxM2S_.capacity());
        sampleUb_ = UINT64_C(1) << newW;
      }
    }


    uint64_t getSampleUb() const noexcept {
      return sampleUb_;
    }


    /*!
     * @brief Return if data structure is ready.
     */
    bool isReady() const noexcept {
      return (srootM_.root_ != nullptr);
    }


    /*!
     * @brief Return if given "idxM" corresponds to valid run.
     */
    bool isValidIdxM(const uint64_t idxM) const noexcept {
      return (isReady() &&
              (idxM / kBtmB) < idxM2S_.size() / kBtmB &&
              (idxM % kBtmB) < getNumChildrenFromBtmM(idxM / kBtmB));
    }


    /*!
     * @brief Return if given "idxS" corresponds to valid run.
     */
    bool isValidIdxS(const uint64_t idxS) const noexcept {
      return (isReady() &&
              (idxS / kBtmB) < idxS2M_.size() / kBtmB &&
              (idxS % kBtmB) < getNumChildrenFromBtmS(idxS / kBtmB));
    }


    /*!
     * @brief Return if given "btmM" is valid
     */
    bool isValidBtmM(const uint64_t btmM) const noexcept {
      if (!(btmM < idxM2S_.size() / kBtmB)) {
        std::cerr << __func__ << ": btmM = " << btmM << std::endl;
      }
      return (isReady() && btmM < idxM2S_.size() / kBtmB);
    }


    /*!
     * @brief Return if given "idxS" corresponds to valid run.
     */
    bool isValidBtmS(const uint64_t btmS) const noexcept {
      return (isReady() && btmS < idxS2M_.size() / kBtmB);
    }


    /*!
     * @brief Get rootM_.
     */
    const auto rootM() const noexcept {
      return srootM_.root_;
    }


    auto rootM() noexcept {
      return srootM_.root_;
    }


    //////////////////////////////// Get values stored for each bottom node via "btm" and "btmNode"
    //////////////// Get parent
    /*!
     * @brief Get parent of "btmM"
     */
    const auto getParentFromBtmM(const uint64_t btmM) const noexcept {
      assert(isValidBtmM(btmM));

      return btmPtrsM_[btmM]->getParent();
    }


    auto getParentFromBtmM(const uint64_t btmM) noexcept {
      assert(isValidBtmM(btmM));

      return btmPtrsM_[btmM]->getParent();
    }


    /*!
     * @brief Get parent of "btmS"
     */
    const auto getParentFromBtmS(const uint64_t btmS) const noexcept {
      assert(isValidBtmS(btmS));

      return parentS_[btmS];
    }


    auto getParentFromBtmS(const uint64_t btmS) noexcept {
      assert(isValidBtmS(btmS));

      return parentS_[btmS];
    }


    //////////////// Get idxInSibling
    /*!
     * @brief Get idxInSibling of "btmM"
     */
    auto getIdxInSiblingFromBtmM(const uint64_t btmM) const noexcept {
      assert(isValidBtmM(btmM));

      return btmPtrsM_[btmM]->getIdxInSibling();
    }


    /*!
     * @brief Get idxInSibling of "btmS"
     */
    auto getIdxInSiblingFromBtmS(const uint64_t btmS) const noexcept {
      assert(isValidBtmS(btmS));

      return idxInSiblingS_[btmS];
    }


    //////////////// Get numChildren
    /*!
     * @brief Get numChildren of "btmM"
     */
    uint8_t getNumChildrenFromBtmM(const uint64_t btmM) const noexcept {
      assert(isValidBtmM(btmM));

      return btmPtrsM_[btmM]->getNumChildren();
    }


    /*!
     * @brief Get numChildren of "btmS"
     */
    uint8_t getNumChildrenFromBtmS(const uint64_t btmS) const noexcept {
      assert(isValidBtmS(btmS));

      return numChildrenS_[btmS];
    }


    //////////////////////////////// Get value stored for each bottom node via "idx" or "btm"
    /*!
     * @brief Get char of run corresponding to "idxS"
     */
    uint64_t getCharFromBtmS(uint64_t btmS) const noexcept {
      assert(isValidBtmS(btmS));

      return charS_[btmS];
    }


    /*!
     * @brief Get char of run corresponding to "idxM"
     */
    uint64_t getCharFromIdxM(uint64_t idxM) const noexcept {
      // {//debug
      //   std::cerr << __func__ << ": idxM = " << idxM << std::endl;
      // }
      assert(isValidIdxM(idxM));

      return getCharFromIdxS(idxM2S(idxM));
    }


    /*!
     * @brief Get char of run corresponding to "idxS"
     */
    uint64_t getCharFromIdxS(uint64_t idxS) const noexcept {
      // {//debug
      //   std::cerr << __func__ << ": idxS = " << idxS << std::endl;
      // }
      assert(isValidIdxS(idxS));

      return charS_[idxS / kBtmB];
    }


    /*!
     * @brief Get TRA-label of run corresponding to "btmM"
     */
    uint64_t getLabelFromBtmM(uint64_t btmM) const noexcept {
      assert(isValidBtmM(btmM));

      return btmPtrsM_[btmM]->getBtmVal();
    }


    //////////////////////////////// Get value stored for each leaf
    //////////////// Get link
    /*!
     * @brief Get "idxS" from "idxM".
     */
    uint64_t idxM2S(const uint64_t idxM) const noexcept {
      assert(isValidIdxM(idxM));

      return idxM2S_.read(idxM);
    }


    /*!
     * @brief Get "idxM" from "idxS".
     */
    uint64_t idxS2M(const uint64_t idxS) const noexcept {
      assert(isValidIdxS(idxS));

      return idxS2M_.read(idxS);
    }


    //////////////// Get leafVal
    /*!
     * @brief Get leafVal from "idxM"
     */
    uint64_t getSampleFromIdxM(uint64_t idxM) const noexcept {
      assert(isValidIdxM(idxM));
      assert(sampleUb_);

      return sample_.read(idxM);
    }


    /*!
     * @brief Get leafVal from "idxS"
     */
    uint64_t getSampleFromIdxS(uint64_t idxS) const noexcept {
      assert(isValidIdxS(idxS));
      assert(sampleUb_);

      return getSampleFromIdxM(idxS2M(idxS));
    }


    //////////////// Get weight
    /*!
     * @brief Get length of run corresponding to "idxM"
     */
    uint64_t getWeightFromIdxM(uint64_t idxM) const noexcept {
      assert(isValidIdxM(idxM));

      return btmPtrsM_[idxM / kBtmB]->readStccVal(idxM % kBtmB);
    }


    /*!
     * @brief Get length of run corresponding to "idxS"
     */
    uint64_t getWeightFromIdxS(uint64_t idxS) const noexcept {
      // {//debug
      //   std::cerr << __func__ << ": idxS = " << idxS << std::endl;
      // }
      assert(isValidIdxS(idxS));

      return getWeightFromIdxM(idxS2M(idxS));
    }


    //////////////////////////////// Get something from BTreeNode
    /*!
     * @brief Get character corresponding to a node of separated tree.
     */
    uint64_t getCharFromNodeS(const BTreeNodeT * nodeS) const noexcept {
      assert(isReady());
      assert(nodeS); // nodeS should be valid node

      return getCharFromBtmS(reinterpret_cast<uint64_t>(nodeS->getLmBtm_DirectJump()));
    }


    //////////////////////////////// Get sum of weights
    /*!
     * @brief Return |T|.
     */
    size_t getSumOfWeight() const noexcept {
      assert(isReady());

      return srootM_.root_->getSumOfWeight();
    }


    /*!
     * @brief Compute num of occ of "ch" in T
     */
    size_t getSumOfWeight(const uint64_t ch) const noexcept {
      assert(isReady());

      const auto * retRootS = searchCharA(ch);
      if (retRootS->isDummy() || getCharFromNodeS(retRootS) != ch) {
        return 0;
      }
      return retRootS->getSumOfWeight();
    }


    uint64_t calcSumOfWeightOfBtmM
    (
     const uint64_t btmM //!< valid btmM
     ) const noexcept {
      assert(isValidBtmM(btmM));

      const auto & stcc = btmPtrsM_[btmM]->getConstRef_stcc();
      const auto num = getNumChildrenFromBtmM(btmM);
      uint64_t sum = 0;
      uint64_t bitPos = 0;
      for (uint64_t i = 0; i < num; ++i) {
        const auto w = stcc.readW(i);
        sum += stcc.readWBits(bitPos, w);
        bitPos += w;
      }

      return sum;
    }


    uint64_t calcSumOfWeightOfBtmM
    (
     const uint64_t btmM,
     const uint8_t childIdx_beg,
     const uint8_t childIdx_end
     ) const noexcept {
      assert(isValidBtmM(btmM));
      assert(childIdx_beg <= childIdx_end);
      assert(childIdx_end <= getNumChildrenFromBtmM(btmM));

      const auto & stcc = btmPtrsM_[btmM]->getConstRef_stcc();
      uint64_t bitPos = stcc.calcBitPos(childIdx_beg);
      uint64_t sum = 0;
      for (uint16_t i = childIdx_beg; i < childIdx_end; ++i) {
        const auto w = stcc.readW(i);
        sum += stcc.readWBits(bitPos, w);
        bitPos += w;
      }

      return sum;
    }


    uint64_t calcSumOfWeightOfBtmS
    (
     const uint64_t btmS
     ) const noexcept {
      assert(isValidBtmS(btmS));

      uint64_t sum = 0;
      for (uint64_t i = 0; i < getNumChildrenFromBtmS(btmS); ++i) { // use OpenMP?
        const auto idxM = idxS2M(btmS * kBtmB + i);
        sum += getWeightFromIdxM(idxM);
      }
      return sum;
    }


    uint64_t calcSumOfWeightOfBtmS
    (
     const uint64_t btmS,
     uint8_t childIdx_beg,
     uint8_t childIdx_end
     ) const noexcept {
      assert(isValidBtmS(btmS));
      assert(childIdx_beg <= childIdx_end);
      assert(childIdx_end <= getNumChildrenFromBtmS(btmS));

      uint64_t sum = 0;
      for (uint16_t i = childIdx_beg; i < childIdx_end; ++i) { // use OpenMP?
        const auto idxM = idxS2M(btmS * kBtmB + i);
        sum += getWeightFromIdxM(idxM);
      }

      return sum;
    }


    //////////////////////////////// rank/select
    /*!
     * @brief Compute rank_{ch}[0..pos], i.e., num of ch in T[0..pos].
     */
    uint64_t rank
    (
     const uint64_t ch, //!< 64bit-char.
     uint64_t pos, //!< Pos (0base) < |T|.
     const bool calcTotalRank //!< If true, compute 'rank_{ch}[0..pos] + num of occ of characters in T smaller than ch'.
     ) const noexcept {
      assert(isReady());
      assert(pos < srootM_.root_->getSumOfWeight());

      const auto idxM = searchPosM(pos); // pos is modified to relative pos
      return rank(ch, idxM, pos, calcTotalRank);
    }


    /*!
     * @brief Variant of rank function, where pos is specified by 'idxM' and 'relativePos'.
     */
    uint64_t rank
    (
     const uint64_t ch, //!< 64bit-char.
     const uint64_t idxM, //!< Valid idxM.
     const uint64_t relativePos, //!< Relative pos (0base) < |T|.
     const bool calcTotalRank //!< If true, compute 'rank_{ch}[0..pos] + num of occ of characters in T smaller than ch'.
     ) const noexcept {
      // {//debug
      //   std::cerr << __func__ << ": ch = " << ch << ", idxM = " << idxM << ", relativePos = " << relativePos << std::endl;
      // }
      assert(isValidIdxM(idxM));
      assert(relativePos < calcSumOfWeightOfBtmM(idxM / kBtmB));

      const auto chNow = getCharFromIdxM(idxM);
      uint64_t ret = 0;
      uint64_t idxS;
      if (ch == chNow) {
        ret = relativePos + 1;
        idxS = idxM2S(idxM);
      } else {
        const auto * retRootS = searchCharA(ch);
        if (retRootS->isDummy() || getCharFromNodeS(retRootS) != ch) {
          return 0;
        }
        idxS = getPredIdxSFromIdxM(retRootS, ch, idxM);
      }
      ret += calcSumOfWeightOfBtmS(idxS / kBtmB, 0, (idxS % kBtmB) + (ch != chNow)); // TODO: calc sum of weights for smaller part
      if (calcTotalRank) {
        BTreeNodeT * root;
        ret += getParentFromBtmS(idxS / kBtmB)->calcPSum(getIdxInSiblingFromBtmS(idxS / kBtmB), root);
        return ret + root->getParent()->calcPSum(root->getIdxInSibling());
      } else {
        return ret + getParentFromBtmS(idxS / kBtmB)->calcPSum(getIdxInSiblingFromBtmS(idxS / kBtmB));
      }
    }


    /*!
     * @brief Compute smallest pos (0base) s.t. 'rank == rank_{ch}[0..pos]'.
     * @attention Rank is 1base.
     */
    uint64_t select
    (
     const BTreeNodeT * rootS, //!< Root of separated tree for 'ch'.
     const uint64_t rank //!< Rank > 0.
     ) const noexcept {
      assert(rank > 0);
      assert(rootS); // rootS should be valid node

      if (rank > rootS->getSumOfWeight()) {
        return BTreeNodeT::NOTFOUND;
      }
      auto pos = rank - 1; // -1 for translating rank into 0base pos.
      const auto idxS = searchPosS(pos, rootS); // pos is modified to the relative pos
      const auto idxM = idxS2M(idxS);
      pos += calcSumOfWeightOfBtmM(idxM / kBtmB, 0, idxM % kBtmB); // TODO: do this and next in parallel
      return pos + getParentFromBtmM(idxM / kBtmB)->calcPSum(getIdxInSiblingFromBtmM(idxM / kBtmB));
    }


    /*!
     * @brief Compute smallest pos s.t. 'rank == rank_{ch}[0..pos]'.
     * @attention Rank is 1base.
     */
    uint64_t select
    (
     const uint64_t ch, //!< character for select query.
     const uint64_t rank //!< Rank > 0.
     ) const noexcept {
      assert(rank > 0);

      const auto * retRootS = searchCharA(ch);
      if (retRootS->isDummy() || getCharFromNodeS(retRootS) != ch) {
        return BTreeNodeT::NOTFOUND;
      }
      return select(retRootS, rank);
    }


    /*!
     * @brief Compute smallest pos s.t. 'totalRank == totalRank_{ch}[0..pos]'.
     * @attention TotalRank is 1base.
     */
    uint64_t select
    (
     const uint64_t totalRank //!< TotalRank > 0.
     ) const noexcept {
      assert(totalRank > 0);

      if (totalRank > srootA_.root_->getSumOfWeight()) {
        return BTreeNodeT::NOTFOUND;
      }
      auto pos = totalRank - 1;
      const auto * retRootS = searchPosA(pos);
      return select(retRootS, pos + 1); // +1 for 1base rank
    }


    /*!
     * @brief Output string represented by current RLE to std::ofstream.
     */
    void printString(std::ofstream & ofs) const noexcept {
      assert(isReady());

      uint64_t pos = 0;
      for (auto idxM = searchPosM(pos); idxM != BTreeNodeT::NOTFOUND; idxM = getNextIdxM(idxM)) {
        const size_t exponent = getWeightFromIdxM(idxM);
        char ch = getCharFromIdxM(idxM);
        for (size_t i = 0; i < exponent; ++i) {
          ofs.put(ch);
        }
      }
    }


  public:
    //////////////////////////////// Public search functions
    /*!
     * @brief Return 'idxM' corresponding to the run containing 'pos'-th character (0base).
     * @attention 'pos' is modified to be the relative position (0base) from the beginning of the run.
     */
    uint64_t searchPosM
    (
     uint64_t & pos //!< [in,out] Give position to search (< |T|). It is modified to relative position.
     ) const noexcept {
      // {//debug
      //   std::cerr << __func__ << ": pos = " << pos << std::endl;
      // }
      assert(isReady());
      assert(pos < srootM_.root_->getSumOfWeight());

      const auto btmM = reinterpret_cast<uint64_t>(srootM_.root_->searchPos(pos));
      const auto & stcc = btmPtrsM_[btmM]->getConstRef_stcc();
      uint8_t i = 0;
      uint64_t bitPos = 0;
      while (true) {
        const auto w = stcc.readW(i);
        const auto weight = stcc.readWBits(bitPos, w);
        if (pos < weight) {
          break;
        }
        ++i;
        bitPos += w;
        pos -= weight;
      }

      return btmM * kBtmB + i;
    }


    /*!
     * @brief Search root of separated tree of the largest character that is smaller or equal to 'ch'.
     */
    const BTreeNodeT * searchCharA
    (
     const uint64_t ch
     ) const noexcept {
      assert(isReady());

      auto * nodeA = srootA_.root_;
      while (true) {
        const bool nowOnBorder = nodeA->isBorder();
        uint8_t lb = 0;
        uint8_t ub = nodeA->getNumChildren();
        while (lb+1 != ub) { // invariant: the answer is in [lb..ub)
          uint8_t mid = (lb + ub) / 2;
          if (ch < getCharFromNodeA(nodeA->getChildPtr(mid), nowOnBorder)) {
            ub = mid;
          } else {
            lb = mid;
          }
        }
        nodeA = nodeA->getChildPtr(lb);
        if (nowOnBorder) {
          return nodeA;
        }
      }
    }


    /*!
     * @brief Search root of separated tree of the largest character that is smaller or equal to 'ch'.
     */
    BTreeNodeT * searchCharA
    (
     const uint64_t ch
     ) noexcept {
      // {//debug
      //   std::cerr << __func__ << ": ch = " << ch << std::endl;
      // }
      assert(isReady());

      return const_cast<BTreeNodeT *>(static_cast<const DynRleForRlbwt &>(*this).searchCharA(ch));
    }


    uint64_t searchPosS
    (
     uint64_t & pos,
     const BTreeNodeT * rootS
     ) const noexcept {
      assert(isReady());
      assert(rootS); // rootS should be valid node
      assert(pos < rootS->getSumOfWeight());

      const auto btmS = reinterpret_cast<uint64_t>(rootS->searchPos(pos));

      uint8_t idx = 0;
      while (true) {
        const uint64_t idxM = idxS2M(btmS * kBtmB + idx);
        auto weight = getWeightFromIdxM(idxM);
        if (pos >= weight) {
          pos -= weight;
          ++idx;
        } else {
          return btmS * kBtmB + idx;
        }
      }
    }


    /*!
     * @brief Search idxS having the largest label that is smaller or equal to "label".
     */
    uint64_t searchLabelS
    (
     const uint64_t label,
     const BTreeNodeT * rootS
     ) const noexcept {
      assert(isReady());
      assert(rootS); // rootS should be valid node

      const auto * nodeS = rootS;
      while (true) {
        const bool nowOnBorder = nodeS->isBorder();
        uint8_t lb = 0;
        uint8_t ub = nodeS->getNumChildren();
        while (lb+1 != ub) {
          uint8_t mid = (lb + ub) / 2;
          if (label < getLabelFromNodeS(nodeS->getChildPtr(mid), nowOnBorder)) {
            ub = mid;
          } else {
            lb = mid;
          }
        }
        nodeS = nodeS->getChildPtr(lb);
        if (nowOnBorder) {
          break;
        }
      }
      const auto btmS = reinterpret_cast<uint64_t>(nodeS);
      uint8_t lb = 0;
      uint8_t ub = getNumChildrenFromBtmS(btmS);
      while (lb+1 != ub) {
        uint8_t mid = (lb + ub) / 2;
        if (label < getLabelFromBtmM(idxS2M(btmS * kBtmB + mid) / kBtmB)) {
          ub = mid;
        } else {
          lb = mid;
        }
      }
      return btmS * kBtmB + lb;
    }


    uint64_t getPredIdxSFromIdxM
    (
     const BTreeNodeT * rootS,
     const uint64_t ch,
     const uint64_t idxM
     ) const noexcept {
      // {//debug
      //   std::cerr << __func__ << ": ch = " << ch << "(" << (char)(ch) << "), idxM = " << idxM << std::endl;
      // }
      assert(isValidIdxM(idxM));

      const uint64_t btmM = idxM / kBtmB;
      const auto btmNodeM = btmPtrsM_[btmM];
      uint8_t i = (idxM % kBtmB) - 1;
      if (btmM) { // If btmM is not 0 (0 means btmM is the first btm in the mixed tree).
        while (i < kBtmB && getCharFromIdxM(btmM * kBtmB + i) != ch) { // "i < kBtmB" breaks when "i becomes below 0"
          --i;
        }
        if (i < kBtmB) {
          return idxM2S(btmM * kBtmB + i);
        } else {
          return searchLabelS(btmNodeM->getBtmVal() - 1, rootS); // -1 is needed.
        }
      } else { // btmM == 0: dummy idx (== 0) should be ignored.
        while (i > 0 && getCharFromIdxM(btmM * kBtmB + i) != ch) {
          --i;
        }
        if (i > 0) {
          return idxM2S(i);
        } else {
          return reinterpret_cast<uintptr_t>(rootS->getLmBtm_DirectJump()) * kBtmB;
        }
      }
    }


  public:
    //////////////////////////////// Iterator like functions
    /*!
     * @brief Get previous btmM
     */
    uint64_t getPrevBtmM
    (
     const uint64_t btmM //!< Valid btmM
     ) const noexcept {
      assert(isValidBtmM(btmM));

      const auto prevBtmM = reinterpret_cast<uint64_t>
        (getParentFromBtmM(btmM)->getPrevBtm(getIdxInSiblingFromBtmM(btmM)));
      if (prevBtmM != BTreeNodeT::NOTFOUND) {
        return prevBtmM;
      }
      return BTreeNodeT::NOTFOUND;
    }


    /*!
     * @brief Get next btmM
     */
    uint64_t getNextBtmM
    (
     const uint64_t btmM //!< Valid btmM
     ) const noexcept {
      assert(isValidBtmM(btmM));

      const auto nextBtmM = reinterpret_cast<uint64_t>
        (getParentFromBtmM(btmM)->getNextBtm_DirectJump(getIdxInSiblingFromBtmM(btmM)));
      if (nextBtmM != BTreeNodeT::NOTFOUND) {
        return nextBtmM;
      }
      return BTreeNodeT::NOTFOUND;
    }


    /*!
     * @brief Get previous btmS
     */
    uint64_t getPrevBtmS
    (
     const uint64_t btmS //!< Valid btmS
     ) const noexcept {
      assert(isValidBtmS(btmS));

      const auto prevBtmS = reinterpret_cast<uint64_t>
        (getParentFromBtmS(btmS)->getPrevBtm(getIdxInSiblingFromBtmS(btmS)));
      if (prevBtmS != BTreeNodeT::NOTFOUND) {
        return prevBtmS;
      }
      return BTreeNodeT::NOTFOUND;
    }


    /*!
     * @brief Get next btmS
     */
    uint64_t getNextBtmS
    (
     const uint64_t btmS //!< Valid btmS
     ) const noexcept {
      assert(isValidBtmS(btmS));

      const auto nextBtmS = reinterpret_cast<uint64_t>
        (getParentFromBtmS(btmS)->getNextBtm_DirectJump(getIdxInSiblingFromBtmS(btmS)));
      if (nextBtmS != BTreeNodeT::NOTFOUND) {
        return nextBtmS;
      }
      return BTreeNodeT::NOTFOUND;
    }


    /*!
     * @brief Get previous idxM.
     */
    uint64_t getPrevIdxM
    (
     const uint64_t idxM //!< Valid idxM.
     ) const noexcept {
      assert(isValidIdxM(idxM));

      if (idxM % kBtmB) {
        return idxM - 1;
      }
      const auto prevBtmM = reinterpret_cast<uint64_t>
        (getParentFromBtmM(idxM / kBtmB)->getPrevBtm(getIdxInSiblingFromBtmM(idxM / kBtmB)));
      if (prevBtmM != BTreeNodeT::NOTFOUND) {
        return prevBtmM * kBtmB + getNumChildrenFromBtmM(prevBtmM) - 1;
      }
      return BTreeNodeT::NOTFOUND;
    }


    /*!
     * @brief Get next idxM.
     */
    uint64_t getNextIdxM
    (
     const uint64_t idxM //!< Valid idxM.
     ) const noexcept {
      assert(isValidIdxM(idxM));

      if ((idxM % kBtmB) + 1 < getNumChildrenFromBtmM(idxM / kBtmB)) {
        return idxM + 1;
      }
      const auto nextBtmM = reinterpret_cast<uint64_t>
        (getParentFromBtmM(idxM / kBtmB)->getNextBtm_DirectJump(getIdxInSiblingFromBtmM(idxM / kBtmB)));
      if (nextBtmM != BTreeNodeT::NOTFOUND) {
        return nextBtmM * kBtmB;
      }
      return BTreeNodeT::NOTFOUND;
    }


    /*!
     * @brief Get first root of separated tree that is not dummy.
     */
    const BTreeNodeT * getFstRootS() const noexcept {
      assert(isReady());

      return getNextRootS(srootA_.root_->getLmBtm_DirectJump());
    }


    /*!
     * @brief Get first root of separated tree that is not dummy.
     */
    BTreeNodeT * getFstRootS() noexcept {
      assert(isReady());

      return const_cast<BTreeNodeT *>(static_cast<const DynRleForRlbwt &>(*this).getFstRootS());
    }


    /*!
     * @brief Get root of separated tree for previous character.
     */
    const BTreeNodeT * getPrevRootS(const BTreeNodeT * node) const noexcept {
      assert(isReady());
      assert(node); // rootS should be valid node

      while (!node->isRoot()) {
        node = node->getParent();
      }
      return node->getParent()->getPrevBtm(node->getIdxInSibling());
    }


    BTreeNodeT * getPrevRootS(BTreeNodeT * nodeS) noexcept {
      assert(isReady());
      assert(nodeS); // nodeS should be valid node

      return const_cast<BTreeNodeT *>(static_cast<const DynRleForRlbwt &>(*this).getPrevRootS(nodeS));
    }


    /*!
     * @brief Get root of separated tree for next character.
     */
    const BTreeNodeT * getNextRootS(const BTreeNodeT * nodeS) const noexcept {
      assert(isReady());
      assert(nodeS); // nodeS should be valid node

      while (!nodeS->isRoot()) {
        nodeS = nodeS->getParent();
      }
      return nodeS->getParent()->getNextBtm_DirectJump(nodeS->getIdxInSibling());
    }


    /*!
     * @brief Get root of separated tree for next character.
     */
    BTreeNodeT * getNextRootS(BTreeNodeT * nodeS) noexcept {
      assert(isReady());
      assert(nodeS); // nodeS should be valid node

      return const_cast<BTreeNodeT *>(getNextRootS(static_cast<const BTreeNodeT *>(nodeS)));
    }


    uint64_t getPrevIdxS(const uint64_t idxS) const noexcept {
      assert(isValidIdxS(idxS));

      if (idxS % kBtmB) {
        return idxS - 1;
      }
      const auto prevBtmS = reinterpret_cast<uint64_t>
        (getParentFromBtmS(idxS / kBtmB)->getPrevBtm(getIdxInSiblingFromBtmS(idxS / kBtmB)));
      if (prevBtmS != BTreeNodeT::NOTFOUND) {
        return prevBtmS * kBtmB + getNumChildrenFromBtmS(prevBtmS) - 1;
      }
      return BTreeNodeT::NOTFOUND;
    }


    uint64_t getNextIdxS(const uint64_t idxS) const noexcept {
      assert(isValidIdxS(idxS));

      if ((idxS % kBtmB) + 1 < getNumChildrenFromBtmS(idxS / kBtmB)) {
        return idxS + 1;
      }
      const auto nextBtmS = reinterpret_cast<uint64_t>
        (getParentFromBtmS(idxS / kBtmB)->getNextBtm_DirectJump(getIdxInSiblingFromBtmS(idxS / kBtmB)));
      if (nextBtmS != BTreeNodeT::NOTFOUND) {
        return nextBtmS * kBtmB;
      }
      return BTreeNodeT::NOTFOUND;
    }


  private:
    //////////////////////////////// private functions (utilities)
    uint64_t getLabelFromNodeS(const BTreeNodeT * nodeS, const bool isChildOfBorder) const noexcept {
      uint64_t btmS;
      if (!isChildOfBorder) {
        btmS = reinterpret_cast<uintptr_t>(nodeS->getLmBtm_DirectJump());
      } else {
        btmS = reinterpret_cast<uintptr_t>(nodeS);
      }
      return getLabelFromBtmM(idxS2M(btmS * kBtmB) / kBtmB);
    }


    uint64_t getCharFromNodeA(const BTreeNodeT * nodeA, const bool isChildOfBorder) const noexcept {
      uint64_t btmS;
      if (!isChildOfBorder) {
        btmS = reinterpret_cast<uintptr_t>(nodeA->getLmBtm_DirectJump()->getLmBtm_DirectJump());
      } else {
        btmS = reinterpret_cast<uintptr_t>(nodeA->getLmBtm_DirectJump());
      }
      return getCharFromBtmS(btmS);
    }


    /*!
     * @brief Return root of separated tree that contains the position 'pos' (0based) in alphabetically sorted array
     */
    BTreeNodeT * searchPosA(uint64_t & pos) const noexcept {
      return srootA_.root_->searchPos(pos);
    }


    void reserveBtm(const size_t numBtms) {
      const uint64_t numIdx = numBtms * kBtmB;
      const uint8_t w = bits::bitSize(numIdx - 1);
      if (sampleUb_) {
        sample_.convert(sample_.getW(), numIdx);
      }
      idxM2S_.convert(w, numIdx);
      idxS2M_.convert(w, numIdx);
      memutil::realloc_AbortOnFail(btmPtrsM_, numBtms);
      memutil::realloc_AbortOnFail(parentS_, numBtms);
      memutil::realloc_AbortOnFail(charS_, numBtms);
      memutil::realloc_AbortOnFail(idxInSiblingS_, numBtms);
      memutil::realloc_AbortOnFail(numChildrenS_, numBtms);
      traCode_ = TagRelabelAlgo::getSmallestTraCode(numBtms);
    }


    void expandBtm() {
      const uint64_t newNum = 2 * idxM2S_.capacity() / kBtmB; // number of capacity of bottoms is doubled
      reserveBtm(newNum);
    }


    uint64_t setNewBtmNodeM
    (
     BtmNodeM * btmNodeM
     ) noexcept {
      assert(btmNodeM != nullptr);

      const uint64_t retBtmIdx = idxM2S_.size() / kBtmB;
      if (retBtmIdx == idxM2S_.capacity() / kBtmB) {
        expandBtm();
      }
      btmPtrsM_[retBtmIdx] = btmNodeM;
      idxM2S_.resize((retBtmIdx + 1) * kBtmB);
      sample_.resize((retBtmIdx + 1) * kBtmB);

      return retBtmIdx;
    }


    uint64_t setNewBtmNodeS
    (
     const uint64_t ch
     ) noexcept {
      const uint64_t retBtmIdx = idxS2M_.size() / kBtmB;
      if (retBtmIdx == idxS2M_.capacity() / kBtmB) {
        expandBtm();
      }
      idxS2M_.resize((retBtmIdx + 1) * kBtmB);
      charS_[retBtmIdx] = ch;
      numChildrenS_[retBtmIdx] = 0;

      return retBtmIdx;
    }


    void asgnLabel
    (
     const uint64_t btmM
     ) noexcept {
      uint64_t next = getNextBtmM(btmM);
      uint64_t prev = getPrevBtmM(btmM); // assume that prev alwarys exists
      uint64_t base = (next == BTreeNodeT::NOTFOUND) ? TagRelabelAlgo::MAX_LABEL : getLabelFromBtmM(next);
      if (getLabelFromBtmM(prev) < base - 1) {
        btmPtrsM_[btmM]->setBtmVal((getLabelFromBtmM(prev) + base) / 2);
        return;
      }

      base >>= 1;
      uint64_t tmpBtmM = btmM;
      uint8_t l = 1;
      uint64_t num = 1;
      uint64_t overflowNum = 2;
      while (true) {
        while (prev != BTreeNodeT::NOTFOUND && (getLabelFromBtmM(prev) >> l) == base) { // expand backward
          ++num;
          tmpBtmM = prev;
          prev = getPrevBtmM(prev);
        }
        while (next != BTreeNodeT::NOTFOUND && (getLabelFromBtmM(next) >> l) == base){ // expand forward
          ++num;
          next = getNextBtmM(next);
        }
        if (overflowNum >= num) {
          break;
        }
        ++l;
        base >>= 1;
        overflowNum = TagRelabelAlgo::getNextOverflowNum(overflowNum, traCode_);
      }

      // relabel num labels
      uint64_t tmpLabel = base << l;
      const uint64_t interval = (UINT64_C(1) << l) / num;
      while (true) {
        btmPtrsM_[tmpBtmM]->setBtmVal(tmpLabel);
        if (--num == 0) {
          return;
        }
        tmpLabel += interval;
        tmpBtmM = getNextBtmM(tmpBtmM);
      }
    }


    void changePSumFromParentM
    (
     const uint64_t btmM,
     const int64_t change
     ) const noexcept {
      // {//debug
      //   std::cerr << __func__ << ": btmM = " << btmM << ", change = " << change << std::endl;
      // }
      assert(isValidBtmM(btmM));

      getParentFromBtmM(btmM)->changePSumFrom(getIdxInSiblingFromBtmM(btmM), change);
    }


    void changePSumFromParentS
    (
     const uint64_t btmS,
     const int64_t change
     ) const noexcept {
      getParentFromBtmS(btmS)->changePSumFrom(getIdxInSiblingFromBtmS(btmS), change);
    }


    uint64_t setupNewSTree
    (
     BTreeNodeT * predNode,
     const uint64_t ch
     ) {
      // {//debug
      //   std::cerr << __func__ << ": predNode = " << predNode << ", ch = " << ch << std::endl;
      // }

      const uint64_t newIdxS = idxS2M_.size();
      const uint64_t btmS = newIdxS / kBtmB;
      if (newIdxS + kBtmB >= idxS2M_.capacity()) {
        expandBtm();
      }
      idxS2M_.resize(newIdxS + kBtmB);
    
      auto * newRootS = new BTreeNodeT(reinterpret_cast<void *>(btmS), true, true, true, false);
      newRootS->putFirstBtm(reinterpret_cast<void *>(btmS), 0);
      parentS_[btmS] = newRootS;
      idxInSiblingS_[btmS] = 0;
      charS_[btmS] = ch;
      numChildrenS_[btmS] = 1; // only dummy idxS exists
      idxS2M_.write(0, newIdxS); // link to dummy idxM of weight 0

      predNode->getParent()->handleSplitOfChild(newRootS, predNode->getIdxInSibling());
      return newIdxS;
    }


    /*!
     * @brief Handle split btm node
     * @post
     *   This function will do the following:
     *   - update
     *     - links from upper nodes (through handleSplitOfBtm())
     *     - links to upper nodes (in this function)
     *     - labels (by asgnLabel()) for kM
     *     - character for kS
     */
    void handleSplitOfBtmInBtmM
    (
     const uint64_t btmM1, //!< First half of splitted node
     const uint64_t btmM2 //!< Second half of splitted node
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << std::endl;
      // }

      auto * uNode = getParentFromBtmM(btmM1);
      const auto idxInSib = getIdxInSiblingFromBtmM(btmM1);
      const auto oriNum = uNode->getNumChildren();
      const uint8_t numToL = uNode->handleSplitOfBtm(reinterpret_cast<BTreeNodeT *>(btmM2), calcSumOfWeightOfBtmM(btmM2), idxInSib);
      if (numToL == 0) {
        for (uint8_t i = idxInSib + 1; i < uNode->getNumChildren(); ++i) {
          auto tmpBtmM = reinterpret_cast<uint64_t>(uNode->getChildPtr(i));
          btmPtrsM_[tmpBtmM]->setParentRef(uNode, i);
        }
        if (oriNum == kB) {
          auto * nextNode = uNode->getNextSib();
          for (uint8_t i = 0; i < nextNode->getNumChildren(); ++i) {
            auto tmpBtmM = reinterpret_cast<uint64_t>(nextNode->getChildPtr(i));
            btmPtrsM_[tmpBtmM]->setParentRef(nextNode, i);
          }
        }
      } else {
        for (uint8_t i = 0; i < uNode->getNumChildren(); ++i) {
          auto tmpBtmM = reinterpret_cast<uint64_t>(uNode->getChildPtr(i));
          btmPtrsM_[tmpBtmM]->setParentRef(uNode, i);
        }
        auto * prevNode = uNode->getPrevSib();
        const uint8_t numL = prevNode->getNumChildren();
        for (uint8_t i = numL - (numToL + (idxInSib < numToL)); i < numL; ++i) {
          auto tmpBtmM = reinterpret_cast<uint64_t>(prevNode->getChildPtr(i));
          btmPtrsM_[tmpBtmM]->setParentRef(prevNode, i);
        }
      }

      asgnLabel(btmM2);
    }


    /*!
     * @brief Handle split btm node
     * @post
     *   This function will do the following:
     *   - update
     *     - links from upper nodes (through handleSplitOfBtm())
     *     - links to upper nodes (in this function)
     *     - labels (by asgnLabel()) for kM
     *     - character for kS
     */
    void handleSplitOfBtmInBtmS
    (
     const uint64_t btmS1, //!< First half of splitted node
     const uint64_t btmS2 //!< Second half of splitted node
     ) noexcept {
      // {//debug
      //   std::cerr << __func__ << ": btmS1 = " << btmS1 << ", btmS2 = " << btmS2 << std::endl;
      // }

      auto * uNode = getParentFromBtmS(btmS1);
      const auto idxInSib = getIdxInSiblingFromBtmS(btmS1);
      const auto oriNum = uNode->getNumChildren();
      const uint8_t numToL = uNode->handleSplitOfBtm(reinterpret_cast<BTreeNodeT *>(btmS2), calcSumOfWeightOfBtmS(btmS2), idxInSib);
      if (numToL == 0) {
        for (uint8_t i = idxInSib + 1; i < uNode->getNumChildren(); ++i) {
          uint64_t tmpBtmS = reinterpret_cast<uintptr_t>(uNode->getChildPtr(i));
          parentS_[tmpBtmS] = uNode;
          idxInSiblingS_[tmpBtmS] = i;
        }
        if (oriNum == kB) {
          auto * nextNode = uNode->getNextSib();
          for (uint8_t i = 0; i < nextNode->getNumChildren(); ++i) {
            uint64_t tmpBtmS = reinterpret_cast<uintptr_t>(nextNode->getChildPtr(i));
            parentS_[tmpBtmS] = nextNode;
            idxInSiblingS_[tmpBtmS] = i;
          }
        }
      } else {
        for (uint8_t i = 0; i < uNode->getNumChildren(); ++i) {
          uint64_t tmpBtmS = reinterpret_cast<uintptr_t>(uNode->getChildPtr(i));
          parentS_[tmpBtmS] = uNode;
          idxInSiblingS_[tmpBtmS] = i;
        }
        auto * prevNode = uNode->getPrevSib();
        const uint8_t numL = prevNode->getNumChildren();
        for (uint8_t i = numL - (numToL + (idxInSib < numToL)); i < numL; ++i) {
          uint64_t tmpBtmS = reinterpret_cast<uintptr_t>(prevNode->getChildPtr(i));
          parentS_[tmpBtmS] = prevNode;
          idxInSiblingS_[tmpBtmS] = i;
        }
      }

      charS_[btmS2] = charS_[btmS1]; // set character of the btm node
    }


    void mvIdxRL
    (
     WBitsVec & wba,
     const uint64_t srcIdx,
     const uint64_t tgtIdx,
     const uint8_t num,
     WBitsVec & wbaOther
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": tgtIdxBase = " << tgtIdxBase << ", srcIdx = " << (int)srcIdx << ", tgtIdx = " << (int)tgtIdx << ", num = " << (int)num
      //             << ", minSupportW = " << (int)minSupportW << ", kM? = " << (btmPtrs_other == btmPtrs_[kS]) << std::endl;
      // }
      assert(srcIdx % kBtmB + num <= static_cast<uint64_t>(kBtmB));
      assert(tgtIdx % kBtmB + num <= static_cast<uint64_t>(kBtmB));

      for (uint8_t i = num; i > 0; --i) {
        const uint64_t idxOther = wba.read(srcIdx + i - 1);
        wba.write(idxOther, tgtIdx + i - 1);
        wbaOther.write(tgtIdx + i - 1, idxOther);
      }
    }


    void mvIdxLR
    (
     WBitsVec & wba,
     const uint64_t srcIdx,
     const uint64_t tgtIdx,
     const uint8_t num,
     WBitsVec & wbaOther
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": tgtIdxBase = " << tgtIdxBase << ", srcIdx = " << (int)srcIdx << ", tgtIdx = " << (int)tgtIdx << ", num = " << (int)num
      //             << ", minSupportW = " << (int)minSupportW << ", kM? = " << (btmPtrs_other == btmPtrs_[kS]) << std::endl;
      // }
      assert(srcIdx % kBtmB + num <= static_cast<uint64_t>(kBtmB));
      assert(tgtIdx % kBtmB + num <= static_cast<uint64_t>(kBtmB));

      for (uint8_t i = 0; i < num; ++i) {
        const uint64_t idxOther = wba.read(srcIdx + i);
        wba.write(idxOther, tgtIdx + i);
        wbaOther.write(tgtIdx + i, idxOther);
      }
    }


    void mvSample
    (
     const uint64_t srcIdx,
     const uint64_t tgtIdx,
     const uint8_t num
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": tgtIdxBase = " << tgtIdxBase << ", srcIdx = " << (int)srcIdx << ", tgtIdx = " << (int)tgtIdx << ", num = " << (int)num
      //             << ", minSupportW = " << (int)minSupportW << ", kM? = " << (btmPtrs_other == btmPtrs_[kS]) << std::endl;
      // }
      if (sampleUb_) {
        mvWBA_SameW(sample_.getItrAt(srcIdx), sample_.getItrAt(tgtIdx), num);
      }
    }


    void makeSpaceInOneBtmNodeM
    (
     const uint64_t idxM,
     const uint64_t * srcWCodes,
     const uint16_t sumW_ins,
     const uint8_t numChild_ins,
     const uint8_t numChild_del //!< Length of wCodes of tgt to delete.
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": idxBase = " << idxBase << ", childIdx = " << (int)childIdx << ", sumW_ins = " << (int)sumW_ins
      //             << ", numChild_ins = " << (int)numChild_ins << ", numChild_del = " << (int)numChild_del << ", insertMorS = " << insertMorS << std::endl;
      //   // btmPtrs_[insertMorS][idxBase / kBtmB]->printDebugInfo(std::cerr);
      // }

      auto btmNodeM = btmPtrsM_[idxM / kBtmB];
      const uint8_t childIdx = idxM % kBtmB;
      uint16_t sumW_del = 0;
      for (uint8_t i = 0; i < numChild_del; ++i) {
        sumW_del = btmNodeM->stcc_.readW(childIdx + i);
      }

      const uint8_t tailNum = btmNodeM->numChildren_ - (childIdx + numChild_del); // at least 0 by assumption.
      uint16_t tailW = 0;
      if (tailNum) {
        tailW = btmNodeM->stccSize_ - btmNodeM->calcBitPos(childIdx + numChild_del);
        btmNodeM->stcc_.mvWCodes(btmNodeM->stcc_.getConstPtr_wCodes(), childIdx + numChild_del, childIdx + numChild_ins, tailNum);
        if (numChild_ins > numChild_del) {
          mvIdxRL(idxM2S_, idxM + numChild_del, idxM + numChild_ins, tailNum, idxS2M_);
        } else {
          mvIdxLR(idxM2S_, idxM + numChild_del, idxM + numChild_ins, tailNum, idxS2M_);
        }
        mvSample(idxM + numChild_del, idxM + numChild_ins, tailNum);
      }
      btmNodeM->stcc_.mvWCodes(srcWCodes, 0, childIdx, numChild_ins);
      btmNodeM->numChildren_ += numChild_ins - numChild_del;
      btmNodeM->updateWCodesAuxM(childIdx, btmNodeM->numChildren_);
      if (sumW_ins != sumW_del) {
        const uint16_t newBitSize = btmNodeM->stccSize_ + sumW_ins - sumW_del;
        btmNodeM->reserveBitCapacity(newBitSize);
        if (tailNum) {
          btmNodeM->stcc_.mvVals(btmNodeM->stcc_.getConstPtr_vals(), btmNodeM->stccSize_ - tailW, newBitSize - tailW, tailW);
        }
        btmNodeM->stccSize_ = newBitSize;
      }
    }


    uint64_t overflowToLeftM
    (
     const uint64_t lIdxBase,
     const uint64_t rIdxBase,
     const uint8_t childIdx,
     const uint64_t * srcWCodes,
     const uint8_t numChild_ins,
     const uint8_t numChild_del //!< Length of wCodes of tgt to delete.
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": lIdxBase = " << lIdxBase << ", rIdxBase = " << rIdxBase << ", childIdx = "
      //             << (int)childIdx << ", numChild_ins = " << (int)numChild_ins << ", numChild_del = " << (int)numChild_del << ", insertMorS = " << insertMorS << std::endl;
      // }
      assert(childIdx + numChild_del <= getNumChildrenFromBtmM(rIdxBase / kBtmB));

      auto lnode = btmPtrsM_[lIdxBase / kBtmB];
      auto rnode = btmPtrsM_[rIdxBase / kBtmB];

      std::tuple<uint64_t, uint64_t, uint64_t> changeList[2];
      uint8_t clSize = 0;
      const uint8_t numL_old = lnode->getNumChildren();
      const uint8_t numR_old = rnode->getNumChildren();
      const uint8_t numTotal = numL_old + numR_old + numChild_ins - numChild_del;
      const uint8_t numL_new = numTotal / 2;
      const uint8_t numR_new = numTotal - numL_new;
      const uint8_t numToLeft = numL_new - numL_old;
      const bool isNewElemInL = childIdx < numToLeft;
      const bool isNewElemInR = childIdx + numChild_ins > numToLeft;
      uint16_t sumWL = lnode->stccSize_;
      uint8_t numL = numL_old;
      uint8_t curNumAfterDel = 0;
      uint8_t curNumSrcWCodes = 0;
      {
        const auto num = (isNewElemInL)? childIdx : numToLeft;
        if (num) {
          const auto w = rnode->calcBitPos(num);
          changeList[clSize++] = {0, sumWL, w};
          sumWL += w;
          lnode->stcc_.mvWCodes(rnode->getConstPtr_wCodes(), 0, numL, num);
          mvIdxLR(idxM2S_, rIdxBase, lIdxBase + numL, num, idxS2M_);
          mvSample(rIdxBase, lIdxBase + numL, num);
          numL += num;
        }
      }
      if (isNewElemInL) {
        curNumSrcWCodes = std::min(static_cast<uint8_t>(numToLeft - childIdx), numChild_ins);
        sumWL += StepCodeUtil::sumW(srcWCodes, 0, curNumSrcWCodes);
        lnode->stcc_.mvWCodes(srcWCodes, 0, numL, curNumSrcWCodes);
        numL += curNumSrcWCodes;
        if (numL < numL_new) { // Still need to move elements to left after inserting srcWCodes
          curNumAfterDel = numL_new - numL;
          const auto w = rnode->stcc_.sumW(childIdx + numChild_del, childIdx + numChild_del + curNumAfterDel);
          changeList[clSize++] = {rnode->calcBitPos(childIdx + numChild_del), sumWL, w};
          sumWL += w;
          lnode->stcc_.mvWCodes(rnode->getConstPtr_wCodes(), childIdx + numChild_del, numL, curNumAfterDel);
          mvIdxLR(idxM2S_, rIdxBase + childIdx + numChild_del, lIdxBase + numL, curNumAfterDel, idxS2M_);
          mvSample(rIdxBase + childIdx + numChild_del, lIdxBase + numL, curNumAfterDel);
        }
      }
      lnode->numChildren_ = numL_new;
      lnode->updateWCodesAuxM(numL_old, numL_new);
      { // Update vals of lnode.
        lnode->reserveBitCapacity(sumWL);
        for (uint8_t i = 0; i < clSize; ++i) {
          lnode->stcc_.mvVals(rnode->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
        }
        lnode->stccSize_ = sumWL;
      }

      // Update rnode.
      clSize = 0;
      const uint16_t bitPosOfLastChunk = rnode->calcBitPos(childIdx + numChild_del + curNumAfterDel);
      const uint16_t bitSizeOfLastChunk = rnode->stccSize_ - bitPosOfLastChunk;
      uint16_t sumWR = bitSizeOfLastChunk;
      if (numToLeft < childIdx) {
        const uint8_t num = childIdx - numToLeft;
        const uint16_t bitPos = rnode->calcBitPos(numToLeft);
        const uint16_t w = rnode->calcBitPos(childIdx) - bitPos;
        sumWR += w;
        changeList[clSize++] = {bitPos, 0, w};
        rnode->stcc_.mvWCodes(rnode->getConstPtr_wCodes(), numToLeft, 0, num);
        mvIdxLR(idxM2S_, rIdxBase + numToLeft, rIdxBase, num, idxS2M_);
        mvSample(rIdxBase + numToLeft, rIdxBase, num);
      }
      if (isNewElemInR) {
        sumWR += StepCodeUtil::sumW(srcWCodes, curNumSrcWCodes, numChild_ins);
      }
      if (numR_old != childIdx + numChild_del) { // There are remaining children in tail.
        if (numR_old != numR_new) { // Need shift wCodes of "this" node.
          const uint8_t srcBeg = childIdx + numChild_del + curNumAfterDel;
          const uint8_t num = numR_old - srcBeg;
          const uint8_t tgtBeg = numR_new - num;
          rnode->stcc_.mvWCodes(rnode->getConstPtr_wCodes(), srcBeg, tgtBeg, num);
          if (tgtBeg < srcBeg) {
            mvIdxLR(idxM2S_, rIdxBase + srcBeg, rIdxBase + tgtBeg, num, idxS2M_);
          } else {
            mvIdxRL(idxM2S_, rIdxBase + srcBeg, rIdxBase + tgtBeg, num, idxS2M_);
          }
          mvSample(rIdxBase + srcBeg, rIdxBase + tgtBeg, num);
        }
        changeList[clSize++] = {bitPosOfLastChunk, sumWR - bitSizeOfLastChunk, bitSizeOfLastChunk};
      }
      if (isNewElemInR) {
        const uint8_t num = numChild_ins - curNumSrcWCodes;
        rnode->stcc_.mvWCodes(srcWCodes, curNumSrcWCodes, childIdx + curNumSrcWCodes - numToLeft, num);
      }
      rnode->numChildren_ = numR_new;
      rnode->updateWCodesAuxM(0, numR_new);
      { // Update vals of rnode
        rnode->reserveBitCapacity(sumWR);
        for (uint8_t i = 0; i < clSize; ++i) {
          rnode->stcc_.mvVals(rnode->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
        }
        rnode->stccSize_ = sumWR;
      }

      return (isNewElemInL) ? lIdxBase + numL_old + childIdx : rIdxBase + childIdx - numToLeft;
    }


    uint64_t overflowToRightM
    (
     const uint64_t lIdxBase,
     const uint64_t rIdxBase,
     const uint8_t childIdx,
     const uint64_t * srcWCodes,
     const uint8_t numChild_ins,
     const uint8_t numChild_del //!< Length of wCodes of tgt to delete.
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": lIdxBase = " << lIdxBase << ", rIdxBase = " << rIdxBase << ", childIdx = "
      //             << (int)childIdx << ", numChild_ins = " << (int)numChild_ins << ", numChild_del = " << (int)numChild_del << ", insertMorS = " << insertMorS << std::endl;
      // }
      assert(childIdx + numChild_del <= getNumChildrenFromBtmM(lIdxBase / kBtmB));

      auto lnode = btmPtrsM_[lIdxBase / kBtmB];
      auto rnode = btmPtrsM_[rIdxBase / kBtmB];

      std::tuple<uint64_t, uint64_t, uint64_t> changeList[2];
      uint8_t clSize = 0;
      const uint8_t numL_old = lnode->getNumChildren();
      const uint8_t numR_old = rnode->getNumChildren();
      const uint8_t numTotal = numL_old + numR_old + numChild_ins - numChild_del;
      const uint8_t numL_new = numTotal / 2;
      const uint8_t numR_new = numTotal - numL_new;
      const uint8_t numToRight = numR_new - numR_old;

      uint8_t numSrcWCodesInL = 0;
      uint8_t numToRight1 = 0;
      uint8_t numToRight2 = 0;
      if (childIdx < numL_new) { // new elements are in L
        if (childIdx + numChild_ins <= numL_new) { // new elements are only in L
          numSrcWCodesInL = numChild_ins;
          numToRight2 = numToRight;
        } else { // new elements are also in R
          numSrcWCodesInL = numL_new - childIdx;
          numToRight2 = numL_old - (childIdx + numChild_del);
          // {
          //   std::cerr << "koko: numToRight2 = " << numToRight2 << std::endl;
          // }
        }
      } else { // new elements are in R
        numToRight1 = childIdx - numL_new;
        numToRight2 = numL_old - (childIdx + numChild_del);
      }

      if (numR_old) { // shift wCodes of R to make space
        rnode->stcc_.mvWCodes(rnode->getConstPtr_wCodes(), 0, numToRight, numR_old);
        mvIdxRL(idxM2S_, rIdxBase, rIdxBase + numToRight, numR_old, idxS2M_);
        mvSample(rIdxBase, rIdxBase + numToRight, numR_old);
      }

      uint8_t numR_increment = 0;
      uint16_t sumWR_increment = 0;
      if (numToRight1) {
        const uint16_t bitPos = lnode->calcBitPos(childIdx - numToRight1);
        const uint16_t w = lnode->calcBitPos(childIdx) - bitPos;
        changeList[clSize++] = {bitPos, 0, w};
        sumWR_increment += w;
        rnode->stcc_.mvWCodes(lnode->getConstPtr_wCodes(), childIdx - numToRight1, 0, numToRight1);
        mvIdxLR(idxM2S_, lIdxBase + childIdx - numToRight1, rIdxBase, numToRight1, idxS2M_);
        mvSample(lIdxBase + childIdx - numToRight1, rIdxBase, numToRight1);
        numR_increment += numToRight1;
      }
      if (numSrcWCodesInL != numChild_ins) {
        sumWR_increment += StepCodeUtil::sumW(srcWCodes, numSrcWCodesInL, numChild_ins);
        rnode->stcc_.mvWCodes(srcWCodes, numSrcWCodesInL, numR_increment, numChild_ins - numSrcWCodesInL);
        numR_increment += (numChild_ins - numSrcWCodesInL);
      }
      if (numToRight2) {
        const uint16_t bitPos = lnode->calcBitPos(numL_old - numToRight2);
        const uint16_t w = lnode->stccSize_ - bitPos;
        changeList[clSize++] = {bitPos, sumWR_increment, w};
        sumWR_increment += w;
        rnode->stcc_.mvWCodes(lnode->getConstPtr_wCodes(), numL_old - numToRight2, numR_increment, numToRight2);
        mvIdxLR(idxM2S_, lIdxBase + numL_old - numToRight2, rIdxBase + numR_increment, numToRight2, idxS2M_);
        mvSample(lIdxBase + numL_old - numToRight2, rIdxBase + numR_increment, numToRight2);
      }
      rnode->numChildren_ = numR_new;
      rnode->updateWCodesAuxM(0, numR_new);
      { // Update vals of "rnode".
        rnode->reserveBitCapacity(rnode->stccSize_ + sumWR_increment);
        if (numR_old) {
          rnode->stcc_.mvVals(rnode->getConstPtr_vals(), 0, sumWR_increment, rnode->stccSize_);
        }
        for (uint8_t i = 0; i < clSize; ++i) {
          rnode->stcc_.mvVals(lnode->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
        }
        rnode->stccSize_ += sumWR_increment;
      }

      if (numSrcWCodesInL) {
        // {
        //   std::cerr << "numSrcWCodesInL = " << numSrcWCodesInL << std::endl;
        // }
        const uint16_t sumWL_ins = static_cast<uint16_t>(StepCodeUtil::sumW(srcWCodes, 0, numSrcWCodesInL));
        const uint16_t tailBitPos_new = lnode->calcBitPos(childIdx) + sumWL_ins;
        lnode->stccSize_ = tailBitPos_new;
        const uint8_t numTail = numL_new - (childIdx + numSrcWCodesInL);
        if (numTail) {
          const uint16_t tailBitPos_old = lnode->calcBitPos(childIdx + numChild_del);
          const uint16_t w = lnode->calcBitPos(childIdx + numChild_del + numTail) - tailBitPos_old;
          lnode->stccSize_ += w;
          if (tailBitPos_new != tailBitPos_old) {
            lnode->reserveBitCapacity(lnode->stccSize_);
            lnode->stcc_.mvVals(lnode->getConstPtr_vals(), tailBitPos_old, tailBitPos_new, w);
          }
          if (numChild_ins != numChild_del) {
            lnode->stcc_.mvWCodes(lnode->getConstPtr_wCodes(), childIdx + numChild_del, childIdx + numChild_ins, numTail);
            if (numChild_ins > numChild_del) {
              mvIdxRL(idxM2S_, lIdxBase + childIdx + numChild_del, lIdxBase + childIdx + numChild_ins, numTail, idxS2M_);
            } else {
              mvIdxLR(idxM2S_, lIdxBase + childIdx + numChild_del, lIdxBase + childIdx + numChild_ins, numTail, idxS2M_);
            }
            mvSample(lIdxBase + childIdx + numChild_del, lIdxBase + childIdx + numChild_ins, numTail);
          }
        } else {
          lnode->reserveBitCapacity(lnode->stccSize_);
        }
        lnode->stcc_.mvWCodes(srcWCodes, 0, childIdx, numSrcWCodesInL);
        lnode->updateWCodesAuxM(childIdx, numL_new);
      } else { // shrink
        lnode->stccSize_ = lnode->calcBitPos(numL_new); // shrink (just change bitSize)
        lnode->shrinkBitCapacity();
        lnode->updateWCodesAuxM(numL_new - 1, numL_new);
      }
      lnode->numChildren_ = numL_new;

      return (numSrcWCodesInL) ? lIdxBase + childIdx : rIdxBase + childIdx - numL_new;
    }


    uint64_t overflowToLeftM_simple
    (
     const uint64_t lIdxBase,
     const uint64_t rIdxBase,
     const uint8_t childIdx,
     const uint64_t * srcWCodes,
     const uint8_t numChild_ins,
     const uint8_t numChild_del //!< Length of wCodes of tgt to delete.
     ) noexcept {
      // {//debug
      //   std::cerr << __func__ << ": lIdxBase = " << lIdxBase << ", rIdxBase = " << rIdxBase << ", childIdx = " << (int)childIdx
      //             << ", numChild_ins = " << (int)numChild_ins << ", numChild_del = " << (int)numChild_del << std::endl;
      // }
      assert(childIdx + numChild_del <= getNumChildrenFromBtmM(rIdxBase / kBtmB));

      auto lnode = btmPtrsM_[lIdxBase / kBtmB];
      auto rnode = btmPtrsM_[rIdxBase / kBtmB];

      std::tuple<uint64_t, uint64_t, uint64_t> changeList[2];
      uint8_t clSize = 0;
      const uint8_t numL_old = lnode->getNumChildren();
      const uint8_t numR_old = rnode->getNumChildren();
      const uint8_t numTotal = numL_old + numR_old + numChild_ins - numChild_del;
      const uint8_t numOldMid = (numL_old + numR_old) / 2;
      const bool isNewElemInL = numL_old + childIdx < numOldMid;
      const uint8_t numL_new = (isNewElemInL) ? numOldMid + numChild_ins - numChild_del : numOldMid;
      const uint8_t numR_new = numTotal - numL_new;
      const uint8_t numToLeft = numL_new - numL_old;

      uint16_t sumWL = lnode->stccSize_;
      uint8_t numL = numL_old;
      uint8_t curNumAfterDel = 0;
      {
        const auto num = (isNewElemInL)? childIdx : numToLeft;
        if (num) {
          const auto w = rnode->calcBitPos(num);
          changeList[clSize++] = {0, sumWL, w};
          sumWL += w;
          lnode->stcc_.mvWCodes(rnode->getConstPtr_wCodes(), 0, numL, num);
          mvIdxLR(idxM2S_, rIdxBase, lIdxBase + numL, num, idxS2M_);
          mvSample(rIdxBase, lIdxBase + numL, num);
          numL += num;
        }
      }
      if (isNewElemInL) {
        sumWL += StepCodeUtil::sumW(srcWCodes, 0, numChild_ins);
        lnode->stcc_.mvWCodes(srcWCodes, 0, numL, numChild_ins);
        numL += numChild_ins;
        if (numL < numL_new) { // Still need to move elements to left after inserting srcWCodes
          curNumAfterDel = numL_new - numL;
          const auto w = rnode->stcc_.sumW(childIdx + numChild_del, childIdx + numChild_del + curNumAfterDel);
          changeList[clSize++] = {rnode->calcBitPos(childIdx + numChild_del), sumWL, w};
          sumWL += w;
          lnode->stcc_.mvWCodes(rnode->getConstPtr_wCodes(), childIdx + numChild_del, numL, curNumAfterDel);
          mvIdxLR(idxM2S_, rIdxBase + childIdx + numChild_del, lIdxBase + numL, curNumAfterDel, idxS2M_);
          mvSample(rIdxBase + childIdx + numChild_del, lIdxBase + numL, curNumAfterDel);
        }
      }
      lnode->numChildren_ = numL_new;
      lnode->updateWCodesAuxM(numL_old, numL_new);
      { // Update vals of lnode.
        lnode->reserveBitCapacity(sumWL);
        for (uint8_t i = 0; i < clSize; ++i) {
          lnode->stcc_.mvVals(rnode->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
        }
        lnode->stccSize_ = sumWL;
      }

      // Update rnode.
      clSize = 0;
      const uint16_t bitPosOfLastChunk = rnode->calcBitPos(childIdx + numChild_del + curNumAfterDel);
      const uint16_t bitSizeOfLastChunk = rnode->stccSize_ - bitPosOfLastChunk;
      uint16_t sumWR = bitSizeOfLastChunk;
      if (numToLeft < childIdx) {
        const uint8_t num = childIdx - numToLeft;
        const uint16_t bitPos = rnode->calcBitPos(numToLeft);
        const uint16_t w = rnode->calcBitPos(childIdx) - bitPos;
        sumWR += w;
        changeList[clSize++] = {bitPos, 0, w};
        rnode->stcc_.mvWCodes(rnode->getConstPtr_wCodes(), numToLeft, 0, num);
        mvIdxLR(idxM2S_, rIdxBase + numToLeft, rIdxBase, num, idxS2M_);
        mvSample(rIdxBase + numToLeft, rIdxBase, num);
      }
      if (!isNewElemInL) {
        sumWR += StepCodeUtil::sumW(srcWCodes, 0, numChild_ins);
      }
      if (numR_old != childIdx + numChild_del) { // There are remaining children in tail.
        if (numR_old != numR_new) { // Need shift wCodes of "this" node.
          const uint8_t srcBeg = childIdx + numChild_del + curNumAfterDel;
          const uint8_t num = numR_old - srcBeg;
          const uint8_t tgtBeg = numR_new - num;
          rnode->stcc_.mvWCodes(rnode->getConstPtr_wCodes(), srcBeg, tgtBeg, num);
          if (tgtBeg < srcBeg) {
            mvIdxLR(idxM2S_, rIdxBase + srcBeg, rIdxBase + tgtBeg, num, idxS2M_);
          } else {
            mvIdxRL(idxM2S_, rIdxBase + srcBeg, rIdxBase + tgtBeg, num, idxS2M_);
          }
          mvSample(rIdxBase + srcBeg, rIdxBase + tgtBeg, num);
        }
        changeList[clSize++] = {bitPosOfLastChunk, sumWR - bitSizeOfLastChunk, bitSizeOfLastChunk};
      }
      if (!isNewElemInL) {
        rnode->stcc_.mvWCodes(srcWCodes, 0, childIdx - numToLeft, numChild_ins);
      }
      rnode->numChildren_ = numR_new;
      rnode->updateWCodesAuxM(0, numR_new);
      { // Update vals of rnode
        rnode->reserveBitCapacity(sumWR);
        for (uint8_t i = 0; i < clSize; ++i) {
          rnode->stcc_.mvVals(rnode->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
        }
        rnode->stccSize_ = sumWR;
        rnode->shrinkBitCapacity();
      }

      return (isNewElemInL) ? lIdxBase + numL_old + childIdx : rIdxBase + childIdx - numToLeft;
    }


    uint64_t overflowToRightM_simple
    (
     const uint64_t lIdxBase,
     const uint64_t rIdxBase,
     const uint8_t childIdx,
     const uint64_t * srcWCodes,
     const uint8_t numChild_ins,
     const uint8_t numChild_del //!< Length of wCodes of tgt to delete.
     ) noexcept {
      // {//debug
      //   std::cerr << __func__ << ": lIdxBase = " << lIdxBase << ", rIdxBase = " << rIdxBase << ", childIdx = " << (int)childIdx
      //             << ", numChild_ins = " << (int)numChild_ins << ", numChild_del = " << (int)numChild_del << std::endl;
      // }
      assert(childIdx + numChild_del <= getNumChildrenFromBtmM(lIdxBase / kBtmB));

      auto lnode = btmPtrsM_[lIdxBase / kBtmB];
      auto rnode = btmPtrsM_[rIdxBase / kBtmB];

      std::tuple<uint64_t, uint64_t, uint64_t> changeList[2];
      uint8_t clSize = 0;
      const uint8_t numL_old = lnode->getNumChildren();
      const uint8_t numR_old = rnode->getNumChildren();
      const uint8_t numOldTotal = numL_old + numR_old;
      const uint8_t numOldMid = numOldTotal / 2;
      const bool isNewElemInL = (childIdx < numOldMid);

      uint8_t numToRight1 = 0;
      uint8_t numToRight2 = 0;
      uint8_t numL_new = numOldMid;
      uint8_t numR_new = numOldTotal - numOldMid;
      if (isNewElemInL) { // new elements are in L
        numToRight2 = numL_old - numOldMid;
        numL_new += numChild_ins - numChild_del;
      } else { // new elements are in R
        numToRight1 = childIdx - numOldMid;
        numToRight2 = numL_old - (childIdx + numChild_del);
        numR_new += numChild_ins - numChild_del;
      }
      const uint8_t numToRight = numR_new - numR_old;

      if (numR_old) { // shift wCodes of R to make space
        rnode->stcc_.mvWCodes(rnode->getConstPtr_wCodes(), 0, numToRight, numR_old);
        mvIdxRL(idxM2S_, rIdxBase, rIdxBase + numToRight, numR_old, idxS2M_);
        mvSample(rIdxBase, rIdxBase + numToRight, numR_old);
      }

      uint8_t numR_increment = 0;
      uint16_t sumWR_increment = 0;
      if (numToRight1) {
        const uint16_t bitPos = lnode->calcBitPos(childIdx - numToRight1);
        const uint16_t w = lnode->calcBitPos(childIdx) - bitPos;
        changeList[clSize++] = {bitPos, 0, w};
        sumWR_increment += w;
        rnode->stcc_.mvWCodes(lnode->getConstPtr_wCodes(), childIdx - numToRight1, 0, numToRight1);
        mvIdxLR(idxM2S_, lIdxBase + childIdx - numToRight1, rIdxBase, numToRight1, idxS2M_);
        mvSample(lIdxBase + childIdx - numToRight1, rIdxBase, numToRight1);
        numR_increment += numToRight1;
      }
      if (!isNewElemInL) {
        sumWR_increment += StepCodeUtil::sumW(srcWCodes, 0, numChild_ins);
        rnode->stcc_.mvWCodes(srcWCodes, 0, numToRight1, numChild_ins);
        numR_increment += numChild_ins;
      }
      if (numToRight2) {
        const uint16_t bitPos = lnode->calcBitPos(numL_old - numToRight2);
        const uint16_t w = lnode->stccSize_ - bitPos;
        changeList[clSize++] = {bitPos, sumWR_increment, w};
        sumWR_increment += w;
        rnode->stcc_.mvWCodes(lnode->getConstPtr_wCodes(), numL_old - numToRight2, numR_increment, numToRight2);
        mvIdxLR(idxM2S_, lIdxBase + numL_old - numToRight2, rIdxBase + numR_increment, numToRight2, idxS2M_);
        mvSample(lIdxBase + numL_old - numToRight2, rIdxBase + numR_increment, numToRight2);
      }
      rnode->numChildren_ = numR_new;
      rnode->updateWCodesAuxM(0, numR_new);
      { // Update vals of "rnode".
        rnode->reserveBitCapacity(rnode->stccSize_ + sumWR_increment);
        if (numR_old) {
          rnode->stcc_.mvVals(rnode->getConstPtr_vals(), 0, sumWR_increment, rnode->stccSize_);
        }
        for (uint8_t i = 0; i < clSize; ++i) {
          rnode->stcc_.mvVals(lnode->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
        }
        rnode->stccSize_ += sumWR_increment;
      }

      if (isNewElemInL) {
        // {
        //   std::cerr << "numSrcWCodesInL = " << numSrcWCodesInL << std::endl;
        // }
        const uint16_t sumWL_ins = static_cast<uint16_t>(StepCodeUtil::sumW(srcWCodes, 0, numChild_ins));
        const uint16_t tailBitPos_new = lnode->calcBitPos(childIdx) + sumWL_ins;
        lnode->stccSize_ = tailBitPos_new;
        const uint8_t numTail = numL_new - (childIdx + numChild_ins);
        if (numTail) {
          const uint16_t tailBitPos_old = lnode->calcBitPos(childIdx + numChild_del);
          const uint16_t w = lnode->calcBitPos(childIdx + numChild_del + numTail) - tailBitPos_old;
          lnode->stccSize_ += w;
          if (tailBitPos_new != tailBitPos_old) {
            lnode->reserveBitCapacity(lnode->stccSize_);
            lnode->stcc_.mvVals(lnode->getConstPtr_vals(), tailBitPos_old, tailBitPos_new, w);
          }
          lnode->stcc_.mvWCodes(lnode->getConstPtr_wCodes(), childIdx + numChild_del, childIdx + numChild_ins, numTail);
          mvIdxRL(idxM2S_, lIdxBase + childIdx + numChild_del, lIdxBase + childIdx + numChild_ins, numTail, idxS2M_);
          mvSample(lIdxBase + childIdx + numChild_del, lIdxBase + childIdx + numChild_ins, numTail);
        } else {
          lnode->reserveBitCapacity(lnode->stccSize_);
        }
        lnode->stcc_.mvWCodes(srcWCodes, 0, childIdx, numChild_ins);
        lnode->updateWCodesAuxM(childIdx, numL_new);
      } else { // shrink
        lnode->stccSize_ = lnode->calcBitPos(numL_new); // shrink (just change bitSize)
        lnode->shrinkBitCapacity();
        lnode->updateWCodesAuxM(numL_new - 1, numL_new);
      }
      lnode->numChildren_ = numL_new;

      return (isNewElemInL) ? lIdxBase + childIdx : rIdxBase + childIdx - numL_new;
    }


    void writeNewElemInOneBtmM
    (
     const uint64_t idxM,
     const uint64_t * newWeights, //!< Storing weights to insert to stcc
     const uint64_t * newLinks, //!< Storing new links to insert
     const uint8_t numChild_ins
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": lnode = " << lnode << ", rnode = " << rnode
      //             << ", childIdx = " << (int)childIdx << ", numChild_ins = " << (int)numChild_ins << std::endl;
      // }
      const uint64_t btmM = idxM / kBtmB;
      const uint8_t childIdx = idxM % kBtmB;
      auto btmNodeM = btmPtrsM_[btmM];
      uint64_t bitPos = btmNodeM->calcBitPos(childIdx);
      for (uint8_t i = childIdx; i < childIdx + numChild_ins; ++i) {
        idxM2S_.write(newLinks[i - childIdx], btmM * kBtmB + i);
        uint8_t w = btmNodeM->stcc_.readW(i);
        btmNodeM->stcc_.writeWBits(newWeights[i - childIdx], bitPos, w);
        bitPos += w;
      }
    }


    void writeNewElemInTwoBtmM
    (
     const uint64_t lIdxBase,
     const uint64_t rIdxBase,
     uint8_t childIdx, //!< Relative idx to write counting from left-end of lnode
     const uint64_t * newWeights, //!< Storing weights to insert to stcc
     const uint64_t * newLinks, //!< Storing new links to insert
     const uint8_t numChild_ins
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": lnode = " << lnode << ", rnode = " << rnode
      //             << ", childIdx = " << (int)childIdx << ", numChild_ins = " << (int)numChild_ins << std::endl;
      // }

      auto lnode = btmPtrsM_[lIdxBase / kBtmB];
      auto rnode = btmPtrsM_[rIdxBase / kBtmB];
      uint8_t numL = lnode->getNumChildren();
      uint8_t numNewElemInL = 0;
      if (childIdx < numL) { // Insert new elements to lnode
        uint64_t bitPos = lnode->calcBitPos(childIdx);
        numNewElemInL = std::min(numChild_ins, static_cast<uint8_t>(numL - childIdx));
        // {//debug
        //   std::cerr << "insert to l: (" << lnode << ")" << numNewElemInL << std::endl;
        //   // lnode->printStatistics(std::cerr, true);
        // }
        for (uint8_t i = childIdx; i < childIdx + numNewElemInL; ++i) {
          idxM2S_.write(newLinks[i - childIdx], lIdxBase + i);
          uint8_t w = lnode->stcc_.readW(i);
          lnode->stcc_.writeWBits(newWeights[i - childIdx], bitPos, w);
          bitPos += w;
        }
      }
      if (numNewElemInL < numChild_ins) { // Insert new elements to rnode
        // {//debug
        //   std::cerr << "insert to r:" << std::endl;
        //   rnode->printStatistics(std::cerr, true);
        // }
        childIdx += numNewElemInL - numL;
        uint64_t bitPos = rnode->calcBitPos(childIdx);
        for (uint8_t i = childIdx; i < childIdx + numChild_ins - numNewElemInL; ++i) {
          idxM2S_.write(newLinks[i - childIdx + numNewElemInL], rIdxBase + i);
          uint8_t w = rnode->stcc_.readW(i);
          // std::cerr << "ci = " << ci
          //           << ", w = " << (int)w
          //           << ", weight = " << weights[ci - childIdx_ins] << std::endl;
          rnode->stcc_.writeWBits(newWeights[i - childIdx + numNewElemInL], bitPos, w);
          bitPos += w;
        }
      }
    }


    /*!
     * @brief Insert stcc values
     * @note Weights of BTreeNodes should be changed in advance
     * @node For simplicity, assume the following two cases (which are necessary for online construction):
     *   - numChild_ins == 1 and numChild_del == 0, or
     *   - numChild_ins == 3 and numChild_del == 1.
     */
    uint64_t insertNewElemM_simple
    (
     const uint64_t idxBase,
     const uint8_t childIdx,
     const uint64_t * newVals, //!< Storing stcc vals to insert
     const uint64_t * newLinks, //!< Storing new links to insert
     const uint8_t numChild_ins,
     const uint8_t numChild_del //!< Length of wCodes of tgt to delete
     ) noexcept {
      assert(numChild_ins <= kBtmB);
      assert(childIdx + numChild_del <= getNumChildrenFromBtmM(idxBase / kBtmB)); // could be equal. Especialy "childIdx" could be "numChildren"
      // {//debug
      //   std::cerr << __func__ << " idxBase = " << idxBase << ", childIdx = " << (int)childIdx
      //             << ", numChild_ins = " << (int)numChild_ins << ", numChild_del = " << (int)numChild_del << std::endl;
      //   // btmPtrs_[insertMorS][idxBase / kBtmB]->printDebugInfo(std::cerr);
      // }

      auto btmNodeM = btmPtrsM_[idxBase / kBtmB];
      uint64_t wCodesTemp[1];
      uint16_t sumW_ins = 0;
      for (uint8_t i = 0; i < numChild_ins; ++i) {
        uint8_t w = StepCodeUtil::calcSteppedW(newVals[i]);
        sumW_ins += w;
        StepCodeUtil::writeWCode(StepCodeUtil::calcWCodeFromSteppedW(w), wCodesTemp, i);
      }

      const uint8_t num_old = btmNodeM->getNumChildren();
      const uint16_t num = static_cast<uint16_t>(num_old) + numChild_ins - numChild_del;
      if (num <= static_cast<uint16_t>(kBtmB)) { // Easy case: This node can accommodate inserting elements.
        makeSpaceInOneBtmNodeM(idxBase + childIdx, wCodesTemp, sumW_ins, numChild_ins, numChild_del);
        writeNewElemInOneBtmM(idxBase + childIdx, newVals, newLinks, numChild_ins);
        return idxBase + childIdx;
      }

      const uint8_t excess = static_cast<uint8_t>(num - kBtmB);
      auto parent = btmNodeM->getParent();
      const auto idxInSib = btmNodeM->getIdxInSibling();
      if (idxInSib) { // Check previous sibling.
        uint64_t lBtmM = reinterpret_cast<uintptr_t>(parent->getChildPtr(idxInSib - 1));
        auto lBtmNodeM = btmPtrsM_[lBtmM];
        const auto numL = lBtmNodeM->getNumChildren();
        if (kBtmB - numL >= excess + 2) { // +2 for simplisity
          const auto retIdx = overflowToLeftM_simple(lBtmM * kBtmB, idxBase, childIdx, wCodesTemp, numChild_ins, numChild_del);
          writeNewElemInOneBtmM(retIdx, newVals, newLinks, numChild_ins);
          parent->changePSumAt(idxInSib - 1, parent->getPSum(idxInSib) + calcSumOfWeightOfBtmM(lBtmM, numL, lBtmNodeM->getNumChildren()));
          return retIdx;
        }
      }

      if (idxInSib + 1 < parent->getNumChildren()) { // Check next sibling.
        uint64_t rBtmM = reinterpret_cast<uintptr_t>(parent->getChildPtr(idxInSib + 1));
        auto rnode = btmPtrsM_[rBtmM];
        const auto numR = rnode->getNumChildren();
        if (kBtmB - numR >= excess + 2) { // +2 for simplisity
          const auto retIdx = overflowToRightM_simple(idxBase, rBtmM * kBtmB, childIdx, wCodesTemp, numChild_ins, numChild_del);
          writeNewElemInOneBtmM(retIdx, newVals, newLinks, numChild_ins);
          parent->changePSumAt(idxInSib, parent->getPSum(idxInSib + 1) - calcSumOfWeightOfBtmM(rBtmM, 0, rnode->getNumChildren() - numR));
          return retIdx;
        }
      }

      { // This bottom node has to be split
        auto rnode = new BtmNodeM();
        const auto rBtmM = setNewBtmNodeM(rnode);
        const auto retIdx = overflowToRightM_simple(idxBase, rBtmM * kBtmB, childIdx, wCodesTemp, numChild_ins, numChild_del);
        writeNewElemInOneBtmM(retIdx, newVals, newLinks, numChild_ins);
        handleSplitOfBtmInBtmM(idxBase / kBtmB, rBtmM);
        return retIdx;
      }
    }


    uint64_t insertRunAfterM
    (
     const uint64_t idxM,
     const uint64_t link
     ) noexcept {
      // {//debug
      //   std::cerr << __func__ << " idxM = " << idxM << ", link = " << link << std::endl;
      // }
      const uint8_t childIdx = (idxM % kBtmB) + 1; // +1 is needed to put new run AFTER "idx". "childIdx" could be "kBtmB"
      const uint64_t newVals[] = {1};
      const uint64_t newLinks[] = {link};
      return insertNewElemM_simple(idxM / kBtmB * kBtmB, childIdx, newVals, newLinks, 1, 0);
    }


    uint64_t insertRunWithSplitM
    (
     const uint64_t idxM,
     const uint64_t splitPos
     ) noexcept {
      // {//debug
      //   std::cerr << __func__ << " idxM = " << idxM << ", splitPos = " << splitPos << std::endl;
      // }

      const uint64_t btmM = idxM / kBtmB;
      auto btmNodeM = btmPtrsM_[btmM];
      const uint8_t childIdx = idxM % kBtmB;
      changePSumFromParentM(btmM, 1);
      const uint64_t weight2 = btmNodeM->readStccVal(childIdx) - splitPos;
      const uint64_t newVals[] = {splitPos, 1, weight2};
      const uint64_t newLinks[] = {idxM2S_.read(btmM * kBtmB + childIdx), 0, 0}; // 0, 0 are dummy
      return insertNewElemM_simple(idxM / kBtmB * kBtmB, childIdx, newVals, newLinks, 3, 1);
    }


    uint64_t insertRunAfterS
    (
     const uint64_t idxS,
     const uint64_t link,
     const uint64_t ch
     ) noexcept {
      // {//debug
      //   std::cerr << __func__ << " idxS = " << idxS << ", link = " << link << ", ch = " << ch << std::endl;
      // }
      const uint8_t childIdx = (idxS % kBtmB) + 1; // +1 is needed to put new run AFTER "idx". "childIdx" could be "kBtmB"
      return insertNewElemS_simple(idxS / kBtmB * kBtmB, childIdx, link, ch);
    }


    /*!
     * @brief Insert stcc values
     * @note Weights of BTreeNodes should be changed in advance
     * @node For simplicity, assume the following two cases (which are necessary for online construction):
     *   - numChild_ins == 1 and numChild_del == 0, or
     *   - numChild_ins == 3 and numChild_del == 1.
     */
    uint64_t insertNewElemS_simple
    (
     const uint64_t idxBase,
     const uint8_t childIdx,
     const uint64_t newLink, //!< new link to insert
     const uint64_t ch
     ) noexcept {
      // {//debug
      //   std::cerr << __func__ << " idxBase = " << idxBase << ", childIdx = " << (int)childIdx
      //             << ", newLink = " << newLink << ", ch = " << ch << std::endl;
      //   // btmPtrs_[insertMorS][idxBase / kBtmB]->printDebugInfo(std::cerr);
      // }
      assert(childIdx <= getNumChildrenFromBtmS(idxBase / kBtmB)); // could be equal

      const uint64_t idxS = idxBase + childIdx;
      const uint64_t btmS = idxBase / kBtmB;
      const uint8_t num_old = getNumChildrenFromBtmS(btmS);
      if (num_old + 1 <= kBtmB) { // Easy case: This node can accommodate inserting elements.
        numChildrenS_[btmS] = num_old + 1;
        mvIdxRL(idxS2M_, idxS, idxS + 1, num_old - childIdx, idxM2S_);
        idxS2M_.write(newLink, idxS);
        return idxS;
      }

      auto parent = getParentFromBtmS(btmS);
      const auto idxInSib = getIdxInSiblingFromBtmS(btmS);
      if (idxInSib) { // Check previous sibling.
        auto lBtmS = reinterpret_cast<uint64_t>(parent->getChildPtr(idxInSib - 1));
        const auto numL = getNumChildrenFromBtmS(lBtmS);
        if (numL < kBtmB - 1) { // -1 for simplisity
          const auto retIdx = overflowToLeftS_simple(lBtmS * kBtmB, idxBase, childIdx);
          idxS2M_.write(newLink, retIdx);
          parent->changePSumAt(idxInSib - 1, parent->getPSum(idxInSib) + calcSumOfWeightOfBtmS(lBtmS, numL, getNumChildrenFromBtmS(lBtmS)));
          return retIdx;
        }
      }

      if (idxInSib + 1 < parent->getNumChildren()) { // Check next sibling.
        auto rBtmS = reinterpret_cast<uint64_t>(parent->getChildPtr(idxInSib + 1));
        const auto numR = getNumChildrenFromBtmS(rBtmS);
        if (numR < kBtmB - 1) { // -1 for simplisity
          const auto retIdx = overflowToRightS_simple(idxBase, rBtmS * kBtmB, childIdx);
          idxS2M_.write(newLink, retIdx);
          parent->changePSumAt(idxInSib, parent->getPSum(idxInSib + 1) - calcSumOfWeightOfBtmS(rBtmS, 0, getNumChildrenFromBtmS(rBtmS) - numR));
          return retIdx;
        }
      }

      { // This bottom node has to be split
        const auto rBtmS = setNewBtmNodeS(ch);
        const auto retIdx = overflowToRightS_simple(idxBase, rBtmS * kBtmB, childIdx);
        idxS2M_.write(newLink, retIdx);
        handleSplitOfBtmInBtmS(idxBase / kBtmB, rBtmS);
        return retIdx;
      }
    }


    uint64_t overflowToLeftS_simple
    (
     const uint64_t lIdxBase,
     const uint64_t rIdxBase,
     const uint8_t childIdx
     ) noexcept {
      // {//debug
      //   std::cerr << __func__ << ": lIdxBase = " << lIdxBase << ", rIdxBase = " << rIdxBase << ", childIdx = " << (int)childIdx << std::endl;
      //   // btmPtrs_[insertMorS][idxBase / kBtmB]->printDebugInfo(std::cerr);
      // }
      const uint8_t numL = getNumChildrenFromBtmS(lIdxBase / kBtmB);
      const uint8_t numToLeft = numL / 2 + kBtmB / 2 - numL;
      const bool isNewElemInL = (childIdx < numToLeft);
      const uint8_t numL_new = numL + numToLeft + isNewElemInL;
      const uint8_t numR_new = kBtmB - numToLeft + !isNewElemInL;
      { // update left node (possibly first part)
        const uint8_t num = (isNewElemInL) ? childIdx : numToLeft;
        if (num) {
          mvIdxLR(idxS2M_, rIdxBase, lIdxBase + numL, num, idxM2S_);
        }
      }
      if (isNewElemInL) { // update left node (second part)
        const uint8_t num = numToLeft - childIdx;
        if (num) {
          mvIdxLR(idxS2M_, rIdxBase + childIdx, lIdxBase + numL_new - num, num, idxM2S_);
        }
      }

      { // update this node (possibly first part)
        const uint8_t num = (isNewElemInL) ? kBtmB - numToLeft : childIdx - numToLeft;
        if (num) {
          mvIdxLR(idxS2M_, rIdxBase + numToLeft, rIdxBase, num, idxM2S_);
        }
      }
      if (!isNewElemInL) { // update this node (second part)
        const uint8_t num = kBtmB - childIdx;
        if (num && kBtmB != numR_new) {
          mvIdxLR(idxS2M_, rIdxBase + kBtmB - num, rIdxBase + numR_new - num, num, idxM2S_);
        }
      }

      numChildrenS_[lIdxBase / kBtmB] = numL_new;
      numChildrenS_[rIdxBase / kBtmB] = numR_new;
      return (isNewElemInL) ? lIdxBase + numL + childIdx : rIdxBase + childIdx - numToLeft;
    }


    uint64_t overflowToRightS_simple
    (
     const uint64_t lIdxBase,
     const uint64_t rIdxBase,
     const uint8_t childIdx
     ) noexcept {
      // {//debug
      //   std::cerr << __func__ << ": lIdxBase = " << lIdxBase << ", rIdxBase = " << rIdxBase << ", childIdx = " << (int)childIdx << std::endl;
      //   // btmPtrs_[insertMorS][idxBase / kBtmB]->printDebugInfo(std::cerr);
      // }
      assert(childIdx <= kBtmB);

      const uint8_t numR = getNumChildrenFromBtmS(rIdxBase / kBtmB);
      const uint8_t leftEnd = numR / 2 + kBtmB /2;
      const bool isNewElemInL = (childIdx < leftEnd);
      const uint8_t shift = kBtmB - leftEnd + !isNewElemInL;
      mvIdxRL(idxS2M_, rIdxBase, rIdxBase + shift, numR, idxM2S_);
      const uint8_t numR_new = numR + shift;

      if (isNewElemInL) {
        mvIdxRL(idxS2M_, lIdxBase + leftEnd, rIdxBase, kBtmB - leftEnd, idxM2S_);
        const uint8_t num = leftEnd - childIdx;
        if (num) {
          mvIdxRL(idxS2M_, lIdxBase + childIdx, lIdxBase + childIdx + 1, num, idxM2S_);
        }
      } else {
        {
          const uint8_t num = childIdx - leftEnd;
          if (num) {
            mvIdxLR(idxS2M_, lIdxBase + leftEnd, rIdxBase, num, idxM2S_);
          }
        }
        {
          const uint8_t num = kBtmB - childIdx;
          if (num) {
            mvIdxLR(idxS2M_, lIdxBase + childIdx, rIdxBase + shift - num, num, idxM2S_);
          }
        }
      }

      numChildrenS_[lIdxBase / kBtmB] = leftEnd + isNewElemInL;
      numChildrenS_[rIdxBase / kBtmB] = numR_new;
      return (isNewElemInL) ? lIdxBase + childIdx : rIdxBase + childIdx - leftEnd;
    }        


    /*!
     * @brief Change (increase/decrease) length of run at "idxM".
     */
    void changeWeight
    (
     const uint64_t idxM,
     const int64_t change
     ) noexcept {
      // {//debug
      //   std::cerr << __func__ << ": idxM = " << idxM << ", change = " << change << std::endl;
      // }
      // update btm node
      auto btmNodeM = btmPtrsM_[idxM / kBtmB];
      const uint64_t curWeight = btmNodeM->readStccVal(idxM % kBtmB);
      assert(curWeight + change > 0);
      const uint64_t newVals[] = {static_cast<uint64_t>(curWeight + change)};
      btmNodeM->replace(newVals, 1, idxM % kBtmB);
      // update mixed tree
      changePSumFromParentM(idxM / kBtmB, change);
      // update separated tree AND alphabet tree (they are connected seamlessly)
      changePSumFromParentS(idxM2S(idxM) / kBtmB, change);
    }


    /*!
     * @brief Insert new run of character 'ch' and length '1' splitting run at "idxM".
     * @return IdxM of the inserted run.
     * @note Assume that "ch" is different from the one for the splitted run.
     */
    uint64_t insertRunWithSplit
    (
     const uint64_t idxM,
     const uint64_t splitPos,
     const uint64_t ch
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << " idxM = " << idxM << ", splitPos = " << splitPos
      //             << ", weight = " << weight << ", leafVal1 = " << leafVal1 << ", leafVal2 = " << leafVal2 << std::endl;
      // }

      const uint64_t idxS0 = idxM2S(idxM);
      uint64_t tempIdxM = insertRunWithSplitM(idxM, splitPos);
      if (idxM != tempIdxM) {
        idxS2M_.write(tempIdxM, idxS0);
      }
      auto * retRootS = searchCharA(ch);
      uint64_t idxS;
      if (retRootS->isDummy() || getCharFromNodeS(retRootS) != ch) {
        idxS = setupNewSTree(retRootS, ch);
      } else {
        idxS = getPredIdxSFromIdxM(retRootS, ch, tempIdxM);
      }
      const auto newIdxM = getNextIdxM(tempIdxM);
      { // insert new run with character "ch"
        tempIdxM = newIdxM;
        changePSumFromParentS(idxS / kBtmB, 1);
        idxS = insertRunAfterS(idxS, tempIdxM, ch);
        idxM2S_.write(idxS, tempIdxM);
      }
      { // insert second half of splitted run
        tempIdxM = getNextIdxM(tempIdxM);
        idxS = insertRunAfterS(idxS0, tempIdxM, ch);
        idxM2S_.write(idxS, tempIdxM);
      }
      return newIdxM;
    }


    /*!
     * @brief Insert new run of character 'ch' and length '1' after 'idxM'.
     * @return IdxM of the inserted run.
     */
    uint64_t insertRunAfter
    (
     const uint64_t idxM,
     const uint64_t ch
     ) noexcept {
      // {//debug
      //   std::cerr << __func__ << ": idxM = " << idxM << ", ch = " << ch << std::endl;
      //   // std::cerr << __func__ << ": BEFORE idxM = " << idxM << std::endl;
      //   // btmPtrs_[kM][idxM / kBtmB]->printDebugInfo(std::cerr);
      // }

      changePSumFromParentM(idxM / kBtmB, 1);
      const auto newIdxM = insertRunAfterM(idxM, 0);
      BTreeNodeT * retRootS = searchCharA(ch);
      uint64_t idxS;
      if (retRootS->isDummy() || getCharFromNodeS(retRootS) != ch) {
        idxS = setupNewSTree(retRootS, ch);
      } else {
        idxS = getPredIdxSFromIdxM(retRootS, ch, newIdxM);
      }
      // {//debug
      //   std::cerr << __FUNCTION__ << ": BEFORE idxS = " << idxS << std::endl;
      //   printDebugInfoOfBtmS(idxS / kBtmB, std::cerr);
      // }
      changePSumFromParentS(idxS / kBtmB, 1);
      const auto newIdxS = insertRunAfterS(idxS, newIdxM, ch);
      idxM2S_.write(newIdxS, newIdxM);
      // parent->changePSumAt(idxInSib - 1, parent->getPSum(idxInSib) + calcSumOfWeightOfBtmNode(lnode, numL, lnode->getNumChildren(), insertMorS));
      // parent->changePSumAt(idxInSib, parent->getPSum(idxInSib + 1) - calcSumOfWeightOfBtmNode(rnode, 0, rnode->getNumChildren() - numR, insertMorS));

      // {//debug
      //   std::cerr << __func__ << ": AFTER idxS = " << idxS << ", newIdxS = " << newIdxS << std::endl;
      //   std::cerr << __func__ << ": AFTER show idxS = " << idxS << std::endl;
      //   printDebugInfoOfBtmS(idxS / kBtmB, std::cerr);
      //   // if (idxM / kBtmB != newIdxM / kBtmB) {
      //   //   std::cerr << __FUNCTION__ << ": AFTER show newIdxS = " << newIdxS << std::endl;
      //   //   btmPtrs_[kS][newIdxS / kBtmB]->printDebugInfo(std::cerr);
      //   // }
      // }
      // {//debug
      //   std::cerr << __FUNCTION__ << ": AFTER idxM = " << idxM << ", newIdxM = " << newIdxM << std::endl;
      //   std::cerr << __FUNCTION__ << ": AFTER show idxM = " << idxM << std::endl;
      //   btmPtrs_[kM][idxM / kBtmB]->printDebugInfo(std::cerr);
      //   if (idxM / kBtmB != newIdxM / kBtmB) {
      //     std::cerr << __FUNCTION__ << ": AFTER show newIdxM = " << newIdxM << std::endl;
      //     btmPtrs_[kM][newIdxM / kBtmB]->printDebugInfo(std::cerr);
      //   }
      // }
      return newIdxM;
    }


  public:
    //////////////////////////////// Public functions (interface)
    /*!
     * @brief Change (increase/decrease) length of run at "idxM".
     */
    void setSample
    (
     const uint64_t idxM,
     const uint64_t newSample
     ) noexcept {
      assert(isValidIdxM(idxM));
      assert(sampleUb_);

      sample_.write(newSample, idxM);
    }


    /*!
     * @brief Pushback a run, merging into the last run if possible.
     */
    uint64_t pushbackRun
    (
     uint64_t & pos, //!< [out] It is set to relative position of a run.
     const uint64_t ch //!< 64bit-char.
     ) {
      // {//debug
      //   std::cerr << __func__ << ": pos = " << pos << ", ch = " << ch << std::endl;
      // }

      const auto btmM = reinterpret_cast<uint64_t>(srootM_.root_->getRmBtm());
      const auto idxM = btmM * kBtmB + getNumChildrenFromBtmM(btmM) - 1;
      const auto idxS = idxM2S(idxM);
      if (idxS == 0) { // dummy
        pos = 0;
        return insertRunAfter(0, ch);
      }
      if (getCharFromBtmS(idxS / kBtmB) != ch) {
        pos = 0;
        return insertRunAfter(idxM, ch);
      } else { // merge into the last run
        pos = getWeightFromIdxM(idxM);
        changeWeight(idxM, 1);
        return idxM;
      }
    }


    /*!
     * @brief Pushback a run without merge.
     */
    // uint64_t pushbackRunWithoutMerge
    // (
    //  const uint64_t ch, //!< 64bit-char.
    //  const uint64_t weight, //!< Weight (exponent) of new run.
    //  const uint64_t leafVal
    //  ) {
    //   const auto btmNodeM = reinterpret_cast<BtmNode *>(srootM_.root_->getRmBtm());
    //   return insertRunAfter(calcIdxBase(btmNodeM, btmPtrs_[kS]) + btmNodeM->getNumChildrenFromBtmM() - 1, weight, leafVal, ch);
    // }


    /*!
     * @brief Insert run of "ch^{1}" at relative position in "idxM", merging into adjacent runs if possible.
     */
    uint64_t insertRun
    (
     uint64_t idxM,
     uint64_t & pos, //!< [in,out] 0base position where inserted run will start. It is modified to relative position in a run.
     const uint64_t ch //!< 64bit-char.
     ) {
      // {//debug
      //   std::cerr << __func__ << ": idxM = " << idxM << ", pos = " << pos << ", ch = " << ch << std::endl;
      // }
      auto chNow = getCharFromIdxM(idxM);
      if (ch == chNow) {
        changeWeight(idxM, 1);
      } else if (pos == 0) {
        idxM = getPrevIdxM(idxM); // Move to previous idxM.
        if (idxM > 0 && ch == getCharFromIdxM(idxM)) { // Check if 'ch' can be merged with the previous run.
          pos = getWeightFromIdxM(idxM);
          changeWeight(idxM, 1);
        } else {
          idxM = insertRunAfter(idxM, ch);
        }
      } else { // Current run is split with fstHalf of weight 'pos'.
        idxM = insertRunWithSplit(idxM, pos, ch);
        pos = 0;
      }
      return idxM;
    }


    /*!
     * @brief Insert run of "ch^{1}" at "pos", merging into adjacent runs if possible.
     */
    uint64_t insertRun
    (
     uint64_t & pos, //!< [in,out] 0base position where inserted run will start. It is modified to relative position in a run.
     const uint64_t ch //!< 64bit-char.
     ) {
      if (pos > srootM_.root_->getSumOfWeight()) {
        return BTreeNodeT::NOTFOUND;
      } else if (pos == srootM_.root_->getSumOfWeight()) {
        return pushbackRun(pos, ch);
      }
      auto idxM = searchPosM(pos); // 'pos' is modified to be the relative pos in the run of 'idxM'.
      return insertRun(idxM, pos, ch);
    }


    /*!
     * @brief Variant of DynRLE::insertRun for rvalue pos.
     */
    uint64_t insertRun
    (
     uint64_t && pos, //!< 0base position where inserted run will start.
     const uint64_t ch //!< 64bit-char.
     ) {
      auto tmp = pos;
      return insertRun(tmp, ch);
    }


  public:
    //////////////////////////////// statistics
    size_t calcMemBytesMTree() const noexcept {
      if (isReady()) {
        return srootM_.root_->calcMemBytes();
      } else {
        return 0;
      }
    }


    size_t calcMemBytesATree() const noexcept {
      if (isReady()) {
        return srootA_.root_->calcMemBytes();
      } else {
        return 0;
      }
    }


    size_t calcMemBytesSTree() const noexcept {
      size_t size = 0;
      if (isReady()) {
        for (const auto * rootS = getFstRootS(); reinterpret_cast<uintptr_t>(rootS) != BTreeNodeT::NOTFOUND; rootS = getNextRootS(rootS)) {
          size += rootS->calcMemBytes();
        }
      }
      return size;
    }


    /*!
     * @brief Return memory for bottom part of MTree
     * @note idxM2S_ and sample_ are excluded
     */
    size_t calcMemBytesBtmM() const noexcept {
      size_t size = (idxM2S_.capacity() / kBtmB) * sizeof(btmPtrsM_[0]);
      const uint64_t numUsed = idxM2S_.size() / kBtmB;
      for (uint64_t i = 0; i < numUsed; ++i) {
        size += btmPtrsM_[i]->calcMemBytes();
      }
      return size;
    }


    /*!
     * @brief Return actual memory used for bottom part of MTree
     * @note idxM2S_ and sample_ are excluded
     * @node Unused (over reserved) part is not counted
     */
    size_t calcMemBytesBtmM_used() const noexcept {
      const uint64_t numUsed = idxM2S_.size() / kBtmB;
      size_t size = numUsed * sizeof(btmPtrsM_[0]);
      for (uint64_t i = 0; i < numUsed; ++i) {
        size += btmPtrsM_[i]->calcMemBytes();
      }
      return size;
    }


    /*!
     * @brief Return memory for bottom part of STree
     * @note idxS2M_ and sample_ are excluded
     */
    size_t calcMemBytesBtmS() const noexcept {
      return (idxS2M_.capacity() / kBtmB) * (sizeof(parentS_[0]) +
                                             sizeof(charS_[0]) +
                                             sizeof(idxInSiblingS_[0]) +
                                             sizeof(numChildrenS_[0]));
    }


    /*!
     * @brief Return actual memory used for bottom part of STree
     * @note idxS2M_ and sample_ are excluded
     * @node Unused (over reserved) part is not counted
     */
    size_t calcMemBytesBtmS_used() const noexcept {
      return (idxS2M_.size() / kBtmB) * (sizeof(parentS_[0]) +
                                         sizeof(charS_[0]) +
                                         sizeof(idxInSiblingS_[0]) +
                                         sizeof(numChildrenS_[0]));
    }


    size_t calcMemBytesLinks() const noexcept {
      return idxM2S_.calcMemBytes() + idxS2M_.calcMemBytes();
    }


    size_t calcMemBytesLinks_used() const noexcept {
      return (idxM2S_.getW() * idxM2S_.size()) / 8 + (idxS2M_.getW() * idxS2M_.size()) / 8;
    }


    size_t calcMemBytesSample() const noexcept {
      return sample_.calcMemBytes();
    }


    size_t calcMemBytesSample_used() const noexcept {
      return (sample_.getW() * sample_.size()) / 8;
    }


    size_t calcMemBytesBtmWeights() const noexcept {
      size_t size = 0;
      const uint64_t numUsed = idxM2S_.size() / kBtmB;
      for (uint64_t i = 0; i < numUsed; ++i) {
        size += btmPtrsM_[i]->calcMemBytesStccDynArray();
      }
      return size;
    }


    size_t calcMemBytesOverReserved() const noexcept {
      size_t size = calcMemBytesBtmM() - calcMemBytesBtmM_used();
      size += calcMemBytesBtmS() - calcMemBytesBtmS_used();
      size += calcMemBytesLinks() - calcMemBytesLinks_used();
      size += calcMemBytesSample() - calcMemBytesSample_used();
      return size;
    }


    size_t calcMemBytes
    (
     bool includeThis = true
     ) const noexcept {
      size_t size = sizeof(*this) * includeThis;
      size += calcMemBytesBtmM();
      size += calcMemBytesBtmS();
      size += calcMemBytesMTree();
      size += calcMemBytesATree();
      size += calcMemBytesSTree();
      size += calcMemBytesLinks();
      size += calcMemBytesSample();
      return size;
    }


    size_t calcNumUsedSTree() const noexcept {
      size_t numUsed = 0;
      if (isReady()) {
        for (const auto * rootS = getFstRootS(); reinterpret_cast<uintptr_t>(rootS) != BTreeNodeT::NOTFOUND; rootS = getNextRootS(rootS)) {
          numUsed += rootS->calcNumUsed();
        }
      }
      return numUsed;
    }


    size_t calcNumSlotsSTree() const noexcept {
      size_t numSlots = 0;
      if (isReady()) {
        for (const auto * rootS = getFstRootS(); reinterpret_cast<uintptr_t>(rootS) != BTreeNodeT::NOTFOUND; rootS = getNextRootS(rootS)) {
          numSlots += rootS->calcNumSlots();
        }
      }
      return numSlots;
    }


    size_t calcNumUsedBtmM() const noexcept {
      size_t numUsed = 0;
      for (uint64_t i = 0; i < idxM2S_.size() / kBtmB; ++i) {
        numUsed += getNumChildrenFromBtmM(i);
      }
      return numUsed;
    }


    size_t calcNumSlotsBtmM() const noexcept {
      return idxM2S_.size();
    }


    size_t calcNumUsedBtmS() const noexcept {
      size_t numUsed = 0;
      for (uint64_t i = 0; i < idxS2M_.size() / kBtmB; ++i) {
        numUsed += getNumChildrenFromBtmS(i);
      }
      return numUsed;
    }


    size_t calcNumSlotsBtmS() const noexcept {
      return idxS2M_.size();
    }


    size_t calcNumRuns() const noexcept {
      size_t numRuns = 0;
      for (size_t i = 0; i < idxM2S_.size() / kBtmB; ++i) {
        numRuns += getNumChildrenFromBtmM(i);
      }
      return numRuns - 1; // -1 due to the first dummy
    }


    size_t calcNumAlph() const noexcept {
      size_t numAlph = 0;
      if (isReady()) {
        for (const auto * rootS = getFstRootS(); reinterpret_cast<uintptr_t>(rootS) != BTreeNodeT::NOTFOUND; rootS = getNextRootS(rootS)) {
          ++numAlph;
        }
      }
      return numAlph;
    }


    void printStatictics(std::ostream & os) const noexcept {
      if (isReady()) {
        const size_t totalLen = getSumOfWeight();
        const size_t numRuns = calcNumRuns();
        const size_t numSlotsM = srootM_.root_->calcNumSlots();
        const size_t numUsedM = srootM_.root_->calcNumUsed();
        const size_t numSlotsA = srootA_.root_->calcNumSlots();
        const size_t numUsedA = srootA_.root_->calcNumUsed();
        const size_t numSlotsS = calcNumSlotsSTree();
        const size_t numUsedS = calcNumUsedSTree();
        const size_t numSlotsBtm = calcNumSlotsBtmM() + calcNumSlotsBtmS();
        const size_t numUsedBtm = calcNumUsedBtmM() + calcNumUsedBtmS();
        const size_t totalBytes = calcMemBytes();
        os << "TotalLen = " << totalLen << ", #Runs = " << numRuns << ", Alphabet Size = " << calcNumAlph()
           << ", BTreeNode arity kB = " << static_cast<uint64_t>(kB) << " BtmNode arity kBtmB = " << static_cast<uint64_t>(kBtmB) << std::endl;
        os << "MTree bottom array size = " << idxM2S_.size() / kBtmB << ", capacity = " << idxM2S_.capacity() / kBtmB << std::endl;
        os << "STree bottom array size = " << idxS2M_.size() / kBtmB << ", capacity = " << idxS2M_.capacity() / kBtmB << std::endl;
        os << "Total: " << totalBytes << " bytes = " << totalBytes / 1024 << " KiB = " << (totalBytes / 1024) / 1024 << " MiB" << std::endl;
        os << "MTree: " << calcMemBytesMTree() << " bytes, OccuRate = " << ((numSlotsM) ? 100.0 * numUsedM / numSlotsM : 0)
           << " (= 100*" << numUsedM << "/" << numSlotsM << ")" << std::endl;
        os << "ATree: " << calcMemBytesATree() << " bytes, OccuRate = " << ((numSlotsA) ? 100.0 * numUsedA / numSlotsA : 0)
           << " (= 100*" << numUsedA << "/" << numSlotsA << ")" << std::endl;
        os << "STree: " << calcMemBytesSTree() << " bytes, OccuRate = " << ((numSlotsS) ? 100.0 * numUsedS / numSlotsS : 0)
           << " (= 100*" << numUsedS << "/" << numSlotsS << ")" << std::endl;
        os << "BtmNodes: " << calcMemBytesBtmM() + calcMemBytesBtmS() << " bytes, OccuRate = " << ((numSlotsBtm) ? 100.0 * numUsedBtm / numSlotsBtm : 0)
           << " (= 100*" << numUsedBtm << "/" << numSlotsBtm << ")" << std::endl;
        os << "Links: " << calcMemBytesLinks() << " bytes" << std::endl;
        os << "Samples: " << calcMemBytesSample() << " bytes" << std::endl;
        os << "---------------- additional details ----------------" << std::endl;
        os << "Dynamic array for weights: " << calcMemBytesBtmWeights() << " bytes" << std::endl;
        os << "Over reserved: " << calcMemBytesOverReserved() << " bytes" << std::endl;
      }
    }


    void printDebugInfoOfBtmS
    (
     const uint64_t btmS,
     std::ostream & os
     ) const noexcept {
      os << __func__
         << ": btmS = " << btmS
         << ", ch = " << getCharFromBtmS(btmS)
         << ", parent_ = " << getParentFromBtmS(btmS)
         << ", idxInSibling = " << (int)getIdxInSiblingFromBtmS(btmS)
         << ", numChildren = " << (int)getNumChildrenFromBtmS(btmS)
         << std::endl;
      os << "idxS2M:";
      for (uint8_t i = 0; i < getNumChildrenFromBtmS(btmS); ++i) {
        const uint64_t idxS = btmS * kBtmB + i;
        os << " [" << idxS << "]" << idxS2M_.read(idxS);
      }
      os << std::endl;
    }


    void printDebugInfoOfBtmM
    (
     const uint64_t btmM,
     std::ostream & os
     ) const noexcept {
      os << __func__ << ": btmM = " << btmM << std::endl;
      btmPtrsM_[btmM]->printDebugInfo(os, true);
      os << "idxM2S:";
      for (uint8_t i = 0; i < getNumChildrenFromBtmM(btmM); ++i) {
        const uint64_t idxM = btmM * kBtmB + i;
        os << " [" << idxM << "]" << idxM2S_.read(idxM);
      }
      os << std::endl;
    }


    void printDebugInfo
    (
     std::ostream & os
     ) const noexcept {
      const uint64_t btmSizeM = idxM2S_.size() / kBtmB;
      const uint64_t btmSizeS = idxS2M_.size() / kBtmB;
      os << "btmSizeM = " << btmSizeM << ", btmCapacityM = " << idxM2S_.capacity() / kBtmB
         << ", btmSizeS = " << btmSizeS << ", btmCapacityS = " << idxS2M_.capacity() / kBtmB
         << ", traCode_ = " << (int)traCode_ << std::endl;
      if (isReady() && getSumOfWeight() > 0) {
        {
          os << "dump btmPtrsM_" << std::endl;
          for (uint64_t i = 0; i < btmSizeM; ++i) {
            os << "[" << i << "]" << btmPtrsM_[i] << " ";
          }
          os << std::endl;
        }

        // {
        //   os << "dump btmPtrs_[kM] debugInfo" << std::endl;
        //   for (uint64_t i = 0; i < size_[kM]; ++i) {
        //     btmPtrs_[kM][i]->printDebugInfo(os);
        //   }
        //   os << "dump btmPtrs_[kS] debugInfo" << std::endl;
        //   for (uint64_t i = 0; i < size_[kS]; ++i) {
        //     btmPtrs_[kS][i]->printDebugInfo(os);
        //   }
        // }

        { // check links of idxM2S and idxS2M
          for (uint64_t i = 0; i < btmSizeM; ++i) {
            for (uint64_t j = 0; j < getNumChildrenFromBtmM(i); ++j) {
              uint64_t idxM = kBtmB * i + j;
              if (idxM != idxS2M(idxM2S(idxM))) {
                os << "error!! links of idxM2S and idxS2M: idxM = " << idxM
                   << ", idxS = " << idxM2S(idxM) << std::endl; // WARNING, links are not maintained correctly
                printDebugInfoOfBtmM(idxM / kBtmB, std::cerr);
                printDebugInfoOfBtmS(idxM2S(idxM) / kBtmB, std::cerr);
              }
            }
          }
        }

        { // check links of parent-child for M
          for (uint64_t i = 0; i < btmSizeM; ++i) {
            uint8_t idx = getIdxInSiblingFromBtmM(i);
            auto node = getParentFromBtmM(i);
            bool islmbtm = (idx == 0);
            if (reinterpret_cast<uint64_t>(node->getChildPtr(idx)) != i) {
              os << "error!! " << "parent-child for btmM = " << i << std::endl;
            }
            if (islmbtm && reinterpret_cast<uint64_t>(node->getLmJumpNode()) != i) {
              os << "error!! lmJumNode for btmM = " << i << std::endl;
            }
            while (!(node->isRoot())) {
              idx = node->getIdxInSibling();
              islmbtm &= (idx == 0);
              if (node->getParent()->getChildPtr(idx) != node) {
                os << "error!! " << "parent-child for child node = " << node << std::endl;
              }
              if (islmbtm && reinterpret_cast<uint64_t>(node->getLmJumpNode()) != i) {
                os << "error!! lmJumNode for btmM = " << i << std::endl;
              }
              node = node->getParent();
            }
          }
        }

        { // check links of parent-child for S
          for (uint64_t i = 0; i < btmSizeS; ++i) {
            uint8_t idx = getIdxInSiblingFromBtmS(i);
            auto node = getParentFromBtmS(i);
            bool islmbtm = (idx == 0);
            if (reinterpret_cast<uint64_t>(node->getChildPtr(idx)) != i) {
              os << "error!! " << "parent-child for btmS = " << i << std::endl;
            }
            if (islmbtm && reinterpret_cast<uint64_t>(node->getLmJumpNode()) != i) {
              os << "error!! lmJumpNode for btmS = " << i << std::endl;
            }
            while (!(node->isRoot())) {
              idx = node->getIdxInSibling();
              islmbtm &= (idx == 0);
              if (node->getParent()->getChildPtr(idx) != node) {
                os << "error!! " << "parent-child for child node = " << node << std::endl;
              }
              if (islmbtm && reinterpret_cast<uint64_t>(node->getLmJumpNode()) != i) {
                os << "error!! lmJumNode for btmM = " << i << std::endl;
              }
              node = node->getParent();
            }
          }
        }

        { // check correctness of runs
          uint64_t c = UINT64_MAX;
          os << "check runs:" << std::endl;
          // std::cerr << srootM_.root_ << " " << srootA_.root_ << std::endl;
          uint64_t pos = 0;
          uint64_t len = 0;
          for (auto idxM = searchPosM(pos); idxM != BTreeNodeT::NOTFOUND; idxM = getNextIdxM(idxM)) {
            ++pos;
            len += getWeightFromIdxM(idxM);
            if (getWeightFromIdxM(idxM) == 0) {
              os << "error!! detected 0 length run: " << idxM << ", " << pos << std::endl;
            }
            if (c == getCharFromIdxM(idxM)) {
              auto idxM0 = getPrevIdxM(idxM);
              os << "error!! detected consecutive runs having the same char: " 
                 << idxM << ", " << pos << ", (" << c << ", " << getWeightFromIdxM(idxM0) << ")" << ", (" << c << ", " << getWeightFromIdxM(idxM) << ")" << std::endl;
            }
            c = getCharFromIdxM(idxM);
          }
          std::cerr << "run: " << pos << ", len: " << len << std::endl;
        }

        // {
        //   uint64_t pos = 0;
        //   for (auto idxM = searchPosM(pos); idxM != BTreeNodeT::NOTFOUND; idxM = getNextIdxM(idxM)) {
        //     os << "(" << idxM << ":" << getCharFromIdxM(idxM) << "^" << getWeightFromIdxM(idxM) << ", " << getLeafValFromIdxM(idxM) << ") ";
        //   }
        //   os << std::endl;
        // }

        // {//MTree
        //   srootM_.root_->printStatistics(std::cerr, true);
        // }

        {
          os << "Information on M" << std::endl;
          uint64_t pos = 0;
          for (uint64_t btmM = reinterpret_cast<uintptr_t>(srootM_.root_->getLmBtm_DirectJump());
               btmM != BTreeNodeT::NOTFOUND;
               btmM = getNextBtmM(btmM)) {
            os << "[" << btmM * kBtmB << "~" << (btmM+1) * kBtmB - 1<< "]";
            printDebugInfoOfBtmM(btmM, std::cerr);
          }
          os << std::endl;
        }

        {
          os << "Alphabet: " << std::endl;
          for (const auto * rootS = getFstRootS();
               reinterpret_cast<uintptr_t>(rootS) != BTreeNodeT::NOTFOUND;
               rootS = getNextRootS(rootS)) {
            os << "(" << getCharFromNodeS(rootS) << ", " << rootS->getSumOfWeight() << ") ";
          }
          os << std::endl;
          os << std::endl;
        }

        {
          os << "Information on S" << std::endl;
          for (const auto * rootS = getFstRootS();
               reinterpret_cast<uintptr_t>(rootS) != BTreeNodeT::NOTFOUND;
               rootS = getNextRootS(rootS)) {
            for (uint64_t btmS = reinterpret_cast<uintptr_t>(rootS->getLmBtm_DirectJump());
                 btmS != BTreeNodeT::NOTFOUND;
                 btmS = getNextBtmS(btmS)) {
              os << "[" << btmS * kBtmB << "~" << (btmS+1) * kBtmB - 1 << "]";
              printDebugInfoOfBtmS(btmS, std::cerr);
            }
          }
        }
      }
    }
  };
} // namespace itmmti

#endif

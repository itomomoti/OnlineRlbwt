/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file DynRleWithValue.hpp
 * @brief Dynamic Run-length encoding with value.
 * @author Tomohiro I
 * @date 2017-12-16
 */
#ifndef INCLUDE_GUARD_DynRleWithValue
#define INCLUDE_GUARD_DynRleWithValue

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
#include "PSumWithValue.hpp"
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
  class DynRleWithValue
  {
  public:
    // Public constant, alias etc.
    static constexpr uint8_t kB{param_kB};
    static constexpr uint8_t kBtmB{param_kBtmB};
    using BTreeNodeT = BTreeNode<kB>;


  private:
    // Inner class
    class BtmNode
    {
    public:
      //// Private member variables.
      uint64_t btmVal_; //! TRA-label for btmM and character for btmS.
      BTreeNodeT * parent_;
      uint64_t * links_; //!< w-bit packed links storing idxM or idxS
      StepCodeCore<kBtmB> stcc_; //!< Storing weights (for MTree) or leaf values (for STree)
      uint16_t stccCapacity_; //!< Current bit capacity of stcc_
      uint16_t stccSize_; //!< Current bit size of stcc_
      uint8_t numChildren_; //!< Current size (number of elements).
      uint8_t w_; //!< Bit-width for "links_"
      uint8_t idxInSibling_; //!< idxInSibling
      uint8_t wCodesAuxM_[kBtmB / StepCodeUtil::kWCNum - 1];


      BtmNode
      (
       uint16_t initStccCapacity = 0
       ) :
        btmVal_(0),
        parent_(nullptr),
        links_(nullptr),
        stcc_(),
        stccCapacity_(0),
        stccSize_(0),
        numChildren_(0),
        w_(0),
        idxInSibling_(0)
      {
        assert(initStccCapacity <= UINT16_MAX - 64);

        reserveBitCapacity(initStccCapacity);
      }


      ~BtmNode() {
        memutil::safefree(links_);
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


      uint8_t getW() const noexcept {
        return w_;
      }


      BTreeNodeT * getParent() const noexcept {
        return parent_;
      }


      BtmNode * getPrevBtmNode() const noexcept {
        auto ret = this->getParent()->getPrevBtm(this->getIdxInSibling());
        if (reinterpret_cast<uintptr_t>(ret) != BTreeNodeT::NOTFOUND) {
          return reinterpret_cast<BtmNode *>(ret);
        } else {
          return nullptr;
        }
      }


      BtmNode * getNextBtmNode() const noexcept {
        auto ret = this->getParent()->getNextBtm_DirectJump(this->getIdxInSibling());
        if (reinterpret_cast<uintptr_t>(ret) != BTreeNodeT::NOTFOUND) {
          return reinterpret_cast<BtmNode *>(ret);
        } else {
          return nullptr;
        }
      }


      uint64_t readStccVal
      (
       const uint8_t idx //!< in [0..numChildren_)
       ) const noexcept {
        assert(idx < numChildren_);

        return stcc_.read(idx);
      }


      uint64_t readLink
      (
       const uint8_t idx //!< in [0..numChildren_)
       ) const noexcept {
        assert(idx < numChildren_);

        return bits::readWBits(links_, static_cast<uint64_t>(idx) * w_, w_, bits::UINTW_MAX(w_));
      }


      void writeLink
      (
       const uint64_t val, //!< in [0, 2^{w_}).
       const uint8_t idx //!< in [0, capacity_).
       ) noexcept {
        assert(idx < kBtmB);
        assert(val <= bits::UINTW_MAX(w_));

        bits::writeWBits(val, links_, static_cast<uint64_t>(idx) * w_, w_, bits::UINTW_MAX(w_));
      }


      void increaseW
      (
       uint8_t minSupportW //!< New bit-width is at least minSupportW
       ) noexcept {
        assert(0 < minSupportW && minSupportW <= 64);

        if (minSupportW <= w_) {
          return;
        }
        const size_t oldLen = (kBtmB * w_ + 63) / 64; // +63 for roundup
        const size_t minLen = (kBtmB * minSupportW + 63) / 64; // +63 for roundup
        if (minLen > oldLen) {
          memutil::realloc_AbortOnFail(links_, minLen);
          minSupportW = minLen * 64 / kBtmB; // Set minSupportW here
        }

        for (uint64_t i = numChildren_ - 1; i != UINT64_MAX; --i) {
          bits::writeWBits(bits::readWBits(links_, i * w_, w_, bits::UINTW_MAX(w_)), links_, i * minSupportW, minSupportW, bits::UINTW_MAX(minSupportW));
        }
        w_ = minSupportW;
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
        const uint64_t end = (idxEnd - 1) / StepCodeUtil::kWCNum + (idxEnd <= (kBtmB - StepCodeUtil::kWCNum));
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
        return size + calcMemBytesStccDynArray() + calcMemBytesLinksArray();
      }


      size_t calcMemBytesStccDynArray() const noexcept {
        return stccCapacity_ / 8;
      }


      size_t calcMemBytesLinksArray() const noexcept {
        return (w_ * kBtmB) / 8;
      }


      void printStatistics
      (
       std::ostream & os,
       const bool verbose
       ) const noexcept {
        os << "DynRleWithValue::BtmNode object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
        os << "BTree arity for bottom node = " << static_cast<int>(kBtmB) << ", btmVal = " << btmVal_ << std::endl;
        os << "parent = " << parent_ << ", idxInSibling = " << (int)idxInSibling_ << ", numChildren = " << static_cast<uint64_t>(numChildren_) << std::endl;
        os << "bit size = " << stccSize_ << ", bit capacity = " << stccCapacity_ << ", bit-width for links = " << (int)w_ << std::endl;
        os << "Total: " << calcMemBytes() << " bytes" << std::endl;
        os << "dynamic array of step code: " << calcMemBytesStccDynArray() << " bytes" << std::endl;
        os << "dynamic array of links: " << calcMemBytesLinksArray() << " bytes" << std::endl;
        os << "DynRleWithValue::BtmNode object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
      }


      void printDebugInfo
      (
       std::ostream & os,
       const bool verbose
       ) const noexcept {
        os << "DynRleWithValue::BtmNode object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
        os << "BTree arity for bottom node = " << static_cast<int>(kBtmB) << ", btmVal = " << btmVal_ << std::endl;
        os << "parent = " << parent_ << ", idxInSibling = " << (int)idxInSibling_ << ", numChildren = " << static_cast<uint64_t>(numChildren_) << std::endl;
        os << "bit size = " << stccSize_ << ", bit capacity = " << stccCapacity_ << ", bit-width for links = " << (int)w_ << std::endl;
        os << "Total: " << calcMemBytes() << " bytes" << std::endl;
        os << "dynamic array of step code: " << calcMemBytesStccDynArray() << " bytes" << std::endl;
        os << "dynamic array of links: " << calcMemBytesLinksArray() << " bytes" << std::endl;
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
        {
          os << "dump links_ (" << links_ << ")" << std::endl;
          for (uint64_t i = 0; i < this->numChildren_; ++i) {
            os << this->readLink(i) << ", ";
          }
          os << std::endl;
        }
        os << "DynRleWithValue::BtmNode object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
      }


      void printDebugInfo
      (
       std::ostream & os,
       const bool verbose,
       const DynRleWithValue<kB, kBtmB> & obj,
       const bool nodeMorS
       ) const noexcept {
        uint64_t btmIdx = obj.calcIdxBase(this, obj.btmPtrs_[!nodeMorS]);
        os << "[" << btmIdx << "~" << btmIdx + kBtmB << ")DynRleWithValue::BtmNode object "
           << ((nodeMorS == obj.kM) ? "kM" : "kS") << " (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
        os << "BTree arity for bottom node = " << static_cast<int>(kBtmB);
        if (nodeMorS == obj.kM) {
          os << ", SumOfWeight = " << obj.calcSumOfWeightOfBtmNodeM(this) << ", label = " << btmVal_ << std::endl;
        } else {
          os << ", SumOfWeight = " << obj.calcSumOfWeightOfBtmNodeS(this) << ", ch = " << btmVal_ << "(" << (char)btmVal_ << ")" << std::endl;
        }
        os << "parent = " << parent_ << ", idxInSibling = " << (int)idxInSibling_ << ", numChildren = " << static_cast<uint64_t>(numChildren_) << std::endl;
        os << "bit size = " << stccSize_ << ", bit capacity = " << stccCapacity_ << ", bit-width for links = " << (int)w_ << std::endl;
        os << "Total: " << calcMemBytes() << " bytes, dynamic array of step code: " << calcMemBytesStccDynArray()
           << " bytes, dynamic array of links: " << calcMemBytesLinksArray() << " bytes" << std::endl;
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
        {
          os << "dump links_ (" << links_ << ")" << std::endl;
          for (uint64_t i = 0; i < this->numChildren_; ++i) {
            os << this->readLink(i) << ", ";
          }
          os << std::endl;
        }
        os << "[" << btmIdx << "~" << btmIdx + kBtmB << ")DynRleWithValue::BtmNode object "
           << ((nodeMorS == obj.kM) ? "kM" : "kS") << " (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
      }
    }; // end of BtmNode


  public:
    // Public constant, alias etc.
    static constexpr size_t kUnitBits{kBtmB * 8};
    static constexpr bool kM{0};
    static constexpr bool kS{1};


  private:
    typename BTreeNodeT::SuperRootT srootM_; //!< Super root of mixed tree
    typename BTreeNodeT::SuperRootT srootA_; //!< Super root of alphabet tree
    BtmNode ** btmPtrs_[2]; //!< Pointers to BtmNodes "btmPtrs_[0]" for MTree and "btmPtrs_[1]" for STree
    uint64_t capacity_[2]; //!< capacities of btmPtrs_ arrays
    uint64_t size_[2]; //!< sizes of btmPtrs_ arrays
    uint8_t traCode_; //!< traCode in [9..16)


  public:
    DynRleWithValue
    (
     const size_t initNumBtms = 0
     ) :
      srootM_(),
      srootA_(),
      traCode_(9)
    {
      btmPtrs_[0] = nullptr;
      btmPtrs_[1] = nullptr;
      capacity_[0] = 0;
      capacity_[1] = 0;
      size_[0] = 0;
      size_[1] = 0;

      if (initNumBtms) {
        init(initNumBtms);
      }
    }


    ~DynRleWithValue() {
      clearAll();
    }


    /*!
     * @brief Reserve space to accomodate 'initNumBtms' bottoms, and init.
     */
    void init
    (
     const size_t initNumBtms
     ) {
      assert(initNumBtms > 0);

      if (isReady()) {
        clearAll();
      }
      reserveBtmM(initNumBtms);
      reserveBtmS(initNumBtms);

      auto firstBtmNodeM = new BtmNode();
      setNewBtmNode(firstBtmNodeM, kM);
      srootM_.setRoot(new BTreeNodeT(firstBtmNodeM, true, true, true, true));
      srootM_.root_->putFirstBtm(firstBtmNodeM, 0);
      firstBtmNodeM->setParentRef(srootM_.root_, 0);
      firstBtmNodeM->increaseW(8);
      const uint64_t newVals[] = {0};
      const uint64_t newLinks[] = {0};
      insertNewElem(0, 0, newVals, newLinks, 1, 0, kM); // sentinel
      // isRoot = true, isBorder = true, isJumpToBtm = true, isUnderSuperRoot = false, isDummy = true
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

      for (uint64_t i = 0; i < size_[kM]; ++i) {
        delete btmPtrs_[kM][i];
      }
      for (uint64_t i = 0; i < size_[kS]; ++i) {
        delete btmPtrs_[kS][i];
      }
      memutil::safefree(btmPtrs_[kM]);
      memutil::safefree(btmPtrs_[kS]);
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
              (idxM / kBtmB) < size_[kM] &&
              (idxM % kBtmB) < btmPtrs_[kM][idxM / kBtmB]->getNumChildren());
    }


    /*!
     * @brief Return if given "idxS" corresponds to valid run.
     */
    bool isValidIdxS(const uint64_t idxS) const noexcept {
      return (isReady() &&
              (idxS / kBtmB) < size_[kS] &&
              (idxS % kBtmB) < btmPtrs_[kS][idxS / kBtmB]->getNumChildren());
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
      return btmPtrs_[kM][btmM]->getParent();
    }


    auto getParentFromBtmM(const uint64_t btmM) noexcept {
      return btmPtrs_[kM][btmM]->getParent();
    }


    /*!
     * @brief Get parent of "btmS"
     */
    const auto getParentFromBtmS(const uint64_t btmS) const noexcept {
      return btmPtrs_[kS][btmS]->getParent();
    }


    auto getParentFromBtmS(const uint64_t btmS) noexcept {
      return btmPtrs_[kS][btmS]->getParent();
    }


    /*!
     * @brief Get parent of "btmNode"
     */
    const auto getParentFromBtmNode
    (
     const void * btmNode //!< Pointer to BtmNode.
     ) const noexcept {
      return static_cast<BtmNode *>(btmNode)->getParent();
    }


    auto getParentFromBtmNode
    (
     void * btmNode //!< Pointer to BtmNode.
     ) noexcept {
      return static_cast<BtmNode *>(btmNode)->getParent();
    }


    //////////////// Get idxInSibling
    /*!
     * @brief Get idxInSibling of "btmM"
     */
    auto getIdxInSiblingFromBtmM(const uint64_t btmM) const noexcept {
      return btmPtrs_[kM][btmM]->getIdxInSibling();
    }


    /*!
     * @brief Get idxInSibling of "btmS"
     */
    auto getIdxInSiblingFromBtmS(const uint64_t btmS) const noexcept {
      return btmPtrs_[kS][btmS]->getIdxInSibling();
    }


    /*!
     * @brief Get idxInSibling of "btmNode"
     */
    auto getIdxInSiblingFromBtmNode
    (
     const void * btmNode //!< Pointer to BtmNode.
     ) const noexcept {
      return static_cast<BtmNode *>(btmNode)->getIdxInSibling();
    }


    //////////////// Get numChildren
    /*!
     * @brief Get numChildren of "btmM"
     */
    auto getNumChildrenFromBtmM(const uint64_t btmM) const noexcept {
      return btmPtrs_[kM][btmM]->getNumChildren();
    }


    /*!
     * @brief Get numChildren of "btmS"
     */
    auto getNumChildrenFromBtmS(const uint64_t btmS) const noexcept {
      return btmPtrs_[kS][btmS]->getNumChildren();
    }


    /*!
     * @brief Get num of "btmNode"
     */
    auto getNumChildrenFromBtmNode
    (
     const void * btmNode //!< Pointer to BtmNode.
     ) const noexcept {
      return static_cast<const BtmNode *>(btmNode)->getNumChildren();
    }


    //////////////////////////////// Get value stored for each bottom node via "idx"
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
      assert(isValidIdxS(idxS));

      return btmPtrs_[kS][idxS / kBtmB]->getBtmVal();
    }


    /*!
     * @brief Get TRA-label of run corresponding to "idxM"
     */
    uint64_t getLabelFromIdxM(uint64_t idxM) const noexcept {
      assert(isValidIdxM(idxM));

      return btmPtrs_[kM][idxM / kBtmB]->getBtmVal();
    }


    //////////////////////////////// Get value stored for each leaf
    //////////////// Get link
    /*!
     * @brief Get "idxS" from "idxM".
     */
    uint64_t idxM2S(const uint64_t idxM) const noexcept {
      assert(isValidIdxM(idxM));

      return btmPtrs_[kM][idxM / kBtmB]->readLink(idxM % kBtmB);
    }


    /*!
     * @brief Get "idxM" from "idxS".
     */
    uint64_t idxS2M(const uint64_t idxS) const noexcept {
      assert(isValidIdxS(idxS));

      return btmPtrs_[kS][idxS / kBtmB]->readLink(idxS % kBtmB);
    }


    /*!
     * @brief Get link from "btmNodeM"
     */
    uint64_t getLinkFromBtmNodeM
    (
     const void * btmNodeM, //!< Pointer to BtmNode for MTree
     const uint8_t childIdx //!< in [0..numChildren of btmNode)
     ) const noexcept {
      return static_cast<BtmNode *>(btmNodeM)->readLink(childIdx);
    }


    /*!
     * @brief Get link from "btmNodeS"
     */
    uint64_t getLinkFromBtmNodeS
    (
     const void * btmNodeS, //!< Pointer to BtmNode for STree
     const uint8_t childIdx //!< in [0..numChildren of btmNode)
     ) const noexcept {
      return static_cast<BtmNode *>(btmNodeS)->readLink(childIdx);
    }


    //////////////// Get leafVal
    /*!
     * @brief Get leafVal from "idxM"
     */
    uint64_t getLeafValFromIdxM(uint64_t idxM) const noexcept {
      assert(isValidIdxM(idxM));

      return getLeafValFromIdxS(idxM2S(idxM));
    }


    /*!
     * @brief Get leafVal from "idxS"
     */
    uint64_t getLeafValFromIdxS(uint64_t idxS) const noexcept {
      assert(isValidIdxS(idxS));

      return btmPtrs_[kS][idxS / kBtmB]->readStccVal(idxS % kBtmB);
    }


    /*!
     * @brief Get leafVal from "btmNodeM"
     */
    uint64_t getLeafValFromBtmNodeM
    (
     const void * btmNodeM, //!< Pointer to BtmNode for MTree
     const uint8_t childIdx //!< in [0..numChildren of btmNode)
     ) const noexcept {
      return getLeafValFromIdxS(static_cast<BtmNode *>(btmNodeM)->readLink(childIdx));
    }


    /*!
     * @brief Get leafVal from "btmNodeS"
     */
    uint64_t getLeafValFromBtmNodeS
    (
     const void * btmNodeS, //!< Pointer to BtmNode for STree
     const uint8_t childIdx //!< in [0..numChildren of btmNode)
     ) const noexcept {
      return static_cast<const BtmNode *>(btmNodeS)->readStccVal(childIdx);
    }


    //////////////// Get weight
    /*!
     * @brief Get length of run corresponding to "idxM"
     */
    uint64_t getWeightFromIdxM(uint64_t idxM) const noexcept {
      assert(isValidIdxM(idxM));

      return btmPtrs_[kM][idxM / kBtmB]->readStccVal(idxM % kBtmB);
    }


    /*!
     * @brief Get length of run corresponding to "idxS"
     */
    uint64_t getWeightFromIdxS(uint64_t idxS) const noexcept {
      assert(isValidIdxS(idxS));

      return getWeightFromIdxM(idxS2M(idxS));
    }


    /*!
     * @brief Get weight from "btmNodeM"
     */
    uint64_t getWeightFromBtmNodeM
    (
     const void * btmNodeM, //!< Pointer to BtmNode for MTree
     const uint8_t childIdx //!< in [0..numChildren of btmNode)
     ) const noexcept {
      return static_cast<BtmNode *>(btmNodeM)->readStccVal(childIdx);
    }


    /*!
     * @brief Get weight from "btmNodeS"
     */
    uint64_t getWeightFromBtmNodeS
    (
     const void * btmNodeS, //!< Pointer to BtmNode for STree
     const uint8_t childIdx //!< in [0..numChildren of btmNode)
     ) const noexcept {
      return getWeightFromIdxM(static_cast<BtmNode *>(btmNodeS)->readLink(childIdx));
    }


    //////////////////////////////// Get something from BTreeNode
    /*!
     * @brief Get character corresponding to a node of separated tree.
     */
    uint64_t getCharFromNodeS(const BTreeNodeT * nodeS) const noexcept {
      assert(isReady());
      assert(nodeS); // nodeS should be valid node

      return reinterpret_cast<const BtmNode *>(nodeS->getLmBtm_DirectJump())->getBtmVal();
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


    uint64_t calcSumOfWeightOfBtmNodeM
    (
     const void * btmNodeM //!< Pointer to BtmNode for MTree 
     ) const noexcept {
      const auto & stcc = static_cast<const BtmNode *>(btmNodeM)->getConstRef_stcc();
      const auto num = static_cast<const BtmNode *>(btmNodeM)->getNumChildren();

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
     const uint64_t btmM
     ) const noexcept {
      assert(isValidIdxM(kBtmB * btmM));

      return calcSumOfWeightOfBtmNodeM(btmPtrs_[kM][btmM]);
    }


    uint64_t calcSumOfWeightOfBtmNodeM
    (
     const void * btmNodeM, //!< Pointer to BtmNode for MTree 
     const uint8_t childIdx_beg,
     const uint8_t childIdx_end
     ) const noexcept {
      const auto & stcc = static_cast<const BtmNode *>(btmNodeM)->getConstRef_stcc();
      uint64_t bitPos = static_cast<const BtmNode *>(btmNodeM)->calcBitPos(childIdx_beg);

      uint64_t sum = 0;
      for (uint16_t i = childIdx_beg; i < childIdx_end; ++i) {
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
      assert(isValidIdxM(kBtmB * btmM));
      assert(childIdx_beg < childIdx_end);
      assert(childIdx_end <= getNumChildrenFromBtmM(btmM));

      return calcSumOfWeightOfBtmNodeM(btmPtrs_[kM][btmM], childIdx_beg, childIdx_end);
    }


    uint64_t calcSumOfWeightOfBtmNodeS
    (
     const void * btmNodeS //!< Pointer to BtmNode for STree 
     ) const noexcept {
      uint64_t sum = 0;
      for (uint64_t i = 0; i < getNumChildrenFromBtmNode(btmNodeS); ++i) {
        const auto idxM = static_cast<const BtmNode *>(btmNodeS)->readLink(i);
        sum += getWeightFromIdxM(idxM);
      }

      return sum;
    }


    uint64_t calcSumOfWeightOfBtmS
    (
     const uint64_t btmS
     ) const noexcept {
      assert(isValidIdxS(kBtmB * btmS));

      return calcSumOfWeightOfBtmNodeS(btmPtrs_[kS][btmS]);
    }


    uint64_t calcSumOfWeightOfBtmNodeS
    (
     const void * btmNodeS, //!< Pointer to BtmNode for STree 
     uint8_t childIdx_beg,
     uint8_t childIdx_end
     ) const noexcept {
      uint64_t sum = 0;
      for (uint16_t i = childIdx_beg; i < childIdx_end; ++i) {
        const auto idxM = static_cast<const BtmNode *>(btmNodeS)->readLink(i);
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
      assert(isValidIdxM(kBtmB * btmS));
      assert(childIdx_beg < childIdx_end);
      assert(childIdx_end <= getNumChildrenFromBtmS(btmS));

      return calcSumOfWeightBtmNodeS(btmPtrs_[kS][btmS], childIdx_beg, childIdx_end);
    }


    uint64_t calcSumOfWeightOfBtmNode
    (
     const void * btmNode, //!< Pointer to BtmNode
     const bool nodeMorS
     ) const noexcept {
      return (nodeMorS == kM) ?
        calcSumOfWeightOfBtmNodeM(btmNode) :
        calcSumOfWeightOfBtmNodeS(btmNode);
    }


    uint64_t calcSumOfWeightOfBtmNode
    (
     const void * btmNode, //!< Pointer to BtmNode
     uint8_t childIdx_beg,
     uint8_t childIdx_end,
     const bool nodeMorS
     ) const noexcept {
      return (nodeMorS == kM) ?
        calcSumOfWeightOfBtmNodeM(btmNode, childIdx_beg, childIdx_end) :
        calcSumOfWeightOfBtmNodeS(btmNode, childIdx_beg, childIdx_end);
    }


    uint64_t calcSumOfWeightOfBtm
    (
     const uint64_t btm,
     const bool nodeMorS
     ) const noexcept {
      return (nodeMorS == kM) ?
        calcSumOfWeightOfBtmM(btm) :
        calcSumOfWeightOfBtmS(btm);
    }


    uint64_t calcSumOfWeightOfBtm
    (
     const uint64_t btm,
     uint8_t childIdx_beg,
     uint8_t childIdx_end,
     const bool nodeMorS
     ) const noexcept {
      return (nodeMorS == kM) ?
        calcSumOfWeightOfBtmM(btm, childIdx_beg, childIdx_end) :
        calcSumOfWeightOfBtmS(btm, childIdx_beg, childIdx_end);
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
      const auto btmNodeS = btmPtrs_[kS][idxS / kBtmB];
      for (uint8_t i = 0; i < (idxS % kBtmB) + (ch != chNow); ++i) {
        ret += getWeightFromIdxM(btmNodeS->readLink(i));
      }
      if (calcTotalRank) {
        BTreeNodeT * root;
        ret += btmNodeS->getParent()->calcPSum(btmNodeS->getIdxInSibling(), root);
        return ret + root->getParent()->calcPSum(root->getIdxInSibling());
      } else {
        return ret + btmNodeS->getParent()->calcPSum(btmNodeS->getIdxInSibling());
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
      const auto btmNodeM = btmPtrs_[kM][idxM / kBtmB];
      pos += calcSumOfWeightOfBtmM(btmNodeM, 0, idxM % kBtmB);
      return pos + btmNodeM->getParent()->calcPSum(btmNodeM->getIdxInSibling());
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
      //   std::cerr << __func__ << ": pos = " << posidxM << std::endl;
      // }
      assert(isReady());
      assert(pos < srootM_.root_->getSumOfWeight());

      const auto btmNodeM = reinterpret_cast<BtmNode *>(srootM_.root_->searchPos(pos));

      const auto & stcc = btmNodeM->getConstRef_stcc();
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

      return calcIdxBase(btmNodeM, btmPtrs_[kS]) + i;
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
      assert(isReady());

      return const_cast<BTreeNodeT *>(static_cast<const DynRleWithValue &>(*this).searchCharA(ch));
    }


    uint64_t searchPosS
    (
     uint64_t & pos,
     const BTreeNodeT * rootS
     ) const noexcept {
      assert(isReady());
      assert(rootS); // rootS should be valid node
      assert(pos < rootS->getSumOfWeight());

      const auto btmNodeS = reinterpret_cast<BtmNode *>(rootS->searchPos(pos));

      uint8_t idx = 0;
      while (true) {
        const uint64_t idxM = btmNodeS->readLink(idx);
        auto weight = getWeightFromIdxM(idxM);
        if (pos >= weight) {
          pos -= weight;
          ++idx;
        } else {
          return calcIdxBase(btmNodeS, btmPtrs_[kM]) + idx;
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
      const auto btmNodeS = reinterpret_cast<const BtmNode *>(nodeS);
      uint8_t lb = 0;
      uint8_t ub = btmNodeS->getNumChildren();
      while (lb+1 != ub) {
        uint8_t mid = (lb + ub) / 2;
        if (label < getLabelFromIdxM(btmNodeS->readLink(mid))) {
          ub = mid;
        } else {
          lb = mid;
        }
      }
      return calcIdxBase(btmNodeS, btmPtrs_[kM]) + lb;
    }


    uint64_t getPredIdxSFromIdxM
    (
     const BTreeNodeT * rootS,
     const uint64_t ch,
     const uint64_t idxM
     ) const noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << " ch = " << ch << "(" << (char)(ch) << "), idxM = " << idxM << std::endl;
      // }
      const uint64_t btmM = idxM / kBtmB;
      const auto btmNodeM = btmPtrs_[kM][btmM];
      uint8_t i = (idxM % kBtmB) - 1;
      if (btmM) { // If btmM is not 0 (0 means btmM is the first btm in the mixed tree).
        while (i < kBtmB && getCharFromIdxS(btmNodeM->readLink(i)) != ch) { // "i < kBtmB" holds when "i becomes below 0"
          --i;
        }
        if (i < kBtmB) {
          return btmNodeM->readLink(i);
        } else {
          return searchLabelS(btmNodeM->getBtmVal() - 1, rootS); // -1 is needed.
        }
      } else { // btmM == 0: dummy idx (== 0) should be ignored.
        while (i > 0 && getCharFromIdxS(btmNodeM->readLink(i)) != ch) {
          --i;
        }
        if (i > 0) {
          return btmNodeM->readLink(i);
        } else {
          return calcIdxBase(reinterpret_cast<const BtmNode *>(rootS->getLmBtm_DirectJump()), btmPtrs_[kM]);
        }
      }
    }


  public:
    //////////////////////////////// Iterator like functions
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
      const auto btmNodeM = btmPtrs_[kM][idxM / kBtmB];
      const auto prev = btmNodeM->getPrevBtmNode();
      if (prev != nullptr) {
        return calcIdxBase(prev, btmPtrs_[kS]) + prev->getNumChildren() - 1;
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

      const auto btmNodeM = btmPtrs_[kM][idxM / kBtmB];
      if ((idxM % kBtmB) + 1 < btmNodeM->getNumChildren()) {
        return idxM + 1;
      }
      const auto next = btmNodeM->getNextBtmNode();
      if (next != nullptr) {
        return calcIdxBase(next, btmPtrs_[kS]);
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

      return const_cast<BTreeNodeT *>(static_cast<const DynRleWithValue &>(*this).getFstRootS());
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

      return const_cast<BTreeNodeT *>(static_cast<const DynRleWithValue &>(*this).getPrevRootS(nodeS));
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
      const auto prev = btmPtrs_[kS][idxS / kBtmB]->getNextBtmNode();
      if (prev != nullptr) {
        return calcIdxBase(prev, btmPtrs_[kM]) + prev->getNumChildren() - 1;
      }
      return BTreeNodeT::NOTFOUND;
    }


    uint64_t getNextIdxS(const uint64_t idxS) const noexcept {
      assert(isValidIdxS(idxS));

      const auto btmNodeS = btmPtrs_[kS][idxS / kBtmB];
      if ((idxS % kBtmB) + 1 < btmNodeS->getNumChildren()) {
        return idxS + 1;
      }
      const auto next = btmNodeS->getNextBtmNode();
      if (next != nullptr) {
        return calcIdxBase(next, btmPtrs_[kM]);
      }
      return BTreeNodeT::NOTFOUND;
    }


  private:
    //////////////////////////////// private functions (utilities)
    /*!
     * @brief Get idxBase of "btmNode"
     */
    uint64_t calcIdxBase
    (
     const BtmNode * btmNode,
     BtmNode ** btmPtrs_other
     ) const noexcept {
      uint64_t otherIdx = btmNode->readLink(btmNode->getNumChildren() - 1);
      return btmPtrs_other[otherIdx / kBtmB]->readLink(otherIdx % kBtmB) & ~(kBtmB - 1);
    }


    uint64_t getLabelFromNodeS(const BTreeNodeT * nodeS, const bool isChildOfBorder) const noexcept {
      uint64_t idxM;
      if (!isChildOfBorder) {
        idxM = reinterpret_cast<const BtmNode *>(nodeS->getLmBtm_DirectJump())->readLink(0);
      } else {
        idxM = reinterpret_cast<const BtmNode *>(nodeS)->readLink(0);
      }
      return getLabelFromIdxM(idxM);
    }


    uint64_t getCharFromNodeA(const BTreeNodeT * nodeA, const bool isChildOfBorder) const noexcept {
      if (!isChildOfBorder) {
        return reinterpret_cast<const BtmNode *>(nodeA->getLmBtm_DirectJump()->getLmBtm_DirectJump())->getBtmVal();
      } else {
        return reinterpret_cast<const BtmNode *>(nodeA->getLmBtm_DirectJump())->getBtmVal();
      }
    }


    /*!
     * @brief Return root of separated tree that contains the position 'pos' (0based) in alphabetically sorted array
     */
    BTreeNodeT * searchPosA(uint64_t & pos) const noexcept {
      return srootA_.root_->searchPos(pos);
    }


    void reserveBtmM(const size_t numBtms) {
      memutil::realloc_AbortOnFail(btmPtrs_[kM], numBtms);
      capacity_[kM] = numBtms;
      traCode_ = TagRelabelAlgo::getSmallestTraCode(numBtms);
    }


    void reserveBtmS(const size_t numBtms) {
      memutil::realloc_AbortOnFail(btmPtrs_[kS], numBtms);
      capacity_[kS] = numBtms;
    }


    void expandBtmM() {
      const uint64_t newNum = 2 * size_[kM]; // number of capacity of bottoms is doubled
      reserveBtmM(newNum);
    }


    void expandBtmS() {
      const uint64_t newNum = 2 * size_[kS]; // number of capacity of bottoms is doubled
      reserveBtmS(newNum);
    }


    uint64_t setNewBtmNode
    (
     BtmNode * btmNode,
     const bool nodeMorS //!< kM or kS
     ) noexcept {
      assert(btmNode != nullptr);

      uint64_t retBtmIdx;
      if (nodeMorS == kM) {
        retBtmIdx = size_[kM];
        if (retBtmIdx == capacity_[kM]) {
          expandBtmM();
        }
        btmPtrs_[kM][retBtmIdx] = btmNode;
        size_[kM] = retBtmIdx + 1;
      } else {
        retBtmIdx = size_[kS];
        if (retBtmIdx == capacity_[kS]) {
          expandBtmS();
        }
        btmPtrs_[kS][retBtmIdx] = btmNode;
        size_[kS] = retBtmIdx + 1;
      }

      return retBtmIdx;
    }


    void asgnLabel
    (
     BtmNode * btmNodeM
     ) noexcept {
      auto nextNodeM = btmNodeM->getNextBtmNode();
      auto prevNodeM = btmNodeM->getPrevBtmNode(); // assume that prev alwarys exists
      uint64_t base = (nextNodeM == nullptr) ? TagRelabelAlgo::MAX_LABEL : nextNodeM->getBtmVal();
      if (btmNodeM->getBtmVal() < base - 1) {
        btmNodeM->setBtmVal((prevNodeM->getBtmVal() + base) / 2);
        return;
      }

      base >>= 1;
      auto tmpNodeM = btmNodeM;
      uint8_t l = 1;
      uint64_t num = 1;
      uint64_t overflowNum = 2;
      while (true) {
        while (prevNodeM != nullptr && (prevNodeM->getBtmVal() >> l) == base) { // expand backward
          ++num;
          tmpNodeM = prevNodeM;
          prevNodeM = prevNodeM->getPrevBtmNode();
        }
        while (nextNodeM != nullptr && (nextNodeM->getBtmVal() >> l) == base){ // expand forward
          ++num;
          nextNodeM = nextNodeM->getNextBtmNode();
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
        tmpNodeM->setBtmVal(tmpLabel);
        if (--num == 0) {
          return;
        }
        tmpLabel += interval;
        tmpNodeM = tmpNodeM->getNextBtmNode();
      }
    }


    void changePSumFromParent
    (
     BtmNode * btmNode,
     const int64_t change
     ) const noexcept {
      btmNode->getParent()->changePSumFrom(btmNode->getIdxInSibling(), change);
    }


    uint64_t setupNewSTree
    (
     BTreeNodeT * predNode,
     const uint64_t ch
     ) {
      // {//debug
      //   std::cerr << __FUNCTION__ << " ch = " << ch << "(" << (char)(ch) << ")" << std::endl;
      // }

      auto newBtmNodeS = new BtmNode();
      auto * newRootS = new BTreeNodeT(newBtmNodeS, true, true, true, false);
      uint64_t newIdxS = setNewBtmNode(newBtmNodeS, kS) * kBtmB;
      newRootS->putFirstBtm(newBtmNodeS, 0);
      newBtmNodeS->setParentRef(newRootS, 0);
      newBtmNodeS->setBtmVal(ch);
      newBtmNodeS->increaseW(8);
      const uint64_t newVals[] = {0};
      const uint64_t newLinks[] = {0};
      insertNewElem(newIdxS, 0, newVals, newLinks, 1, 0, kS); // dummy idxS
      newBtmNodeS->writeLink(0, 0);

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
    void handleSplitOfBtmInBtm
    (
     BtmNode * btmNode1, // First half of splitted node
     BtmNode * btmNode2, // Second half of splitted node
     const bool nodeMorS // kM or kS
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << std::endl;
      // }

      auto * uNode = btmNode1->getParent();
      const auto idxInSib = btmNode1->getIdxInSibling();
      const auto oriNum = uNode->getNumChildren();
      const uint8_t numToL = uNode->handleSplitOfBtm(reinterpret_cast<BTreeNodeT *>(btmNode2), calcSumOfWeightOfBtmNode(btmNode2, nodeMorS), idxInSib);
      if (numToL == 0) {
        for (uint8_t i = idxInSib + 1; i < uNode->getNumChildren(); ++i) {
          auto tmpBtmNode = reinterpret_cast<BtmNode *>(uNode->getChildPtr(i));
          tmpBtmNode->setParentRef(uNode, i);
        }
        if (oriNum == kB) {
          auto * nextNode = uNode->getNextSib();
          for (uint8_t i = 0; i < nextNode->getNumChildren(); ++i) {
            auto tmpBtmNode = reinterpret_cast<BtmNode *>(nextNode->getChildPtr(i));
            tmpBtmNode->setParentRef(nextNode, i);
          }
        }
      } else {
        for (uint8_t i = 0; i < uNode->getNumChildren(); ++i) {
          auto tmpBtmNode = reinterpret_cast<BtmNode *>(uNode->getChildPtr(i));
          tmpBtmNode->setParentRef(uNode, i);
        }
        auto * prevNode = uNode->getPrevSib();
        const uint8_t numL = prevNode->getNumChildren();
        for (uint8_t i = numL - (numToL + (idxInSib < numToL)); i < numL; ++i) {
          auto tmpBtmNode = reinterpret_cast<BtmNode *>(prevNode->getChildPtr(i));
          tmpBtmNode->setParentRef(prevNode, i);
        }
      }

      if (nodeMorS == kM) {
        asgnLabel(btmNode2);
      } else {
        btmNode2->setBtmVal(btmNode1->getBtmVal()); // set character of the btm node
      }
    }


    void mvIdxRL
    (
     BtmNode * srcBtmNode,
     BtmNode * tgtBtmNode,
     const uint64_t tgtIdxBase,
     const uint8_t srcIdx,
     const uint8_t tgtIdx,
     const uint8_t num,
     const uint8_t minSupportW,
     BtmNode ** btmPtrs_other
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": tgtIdxBase = " << tgtIdxBase << ", srcIdx = " << (int)srcIdx << ", tgtIdx = " << (int)tgtIdx << ", num = " << (int)num
      //             << ", minSupportW = " << (int)minSupportW << ", kM? = " << (btmPtrs_other == btmPtrs_[kS]) << std::endl;
      // }
      assert(srcIdx + num <= kBtmB);
      assert(tgtIdx + num <= kBtmB);

      for (uint8_t i = num; i > 0; --i) {
        const uint8_t src = srcIdx + i - 1;
        const uint8_t tgt = tgtIdx + i - 1;
        const uint64_t idx_other = srcBtmNode->readLink(src);
        tgtBtmNode->writeLink(idx_other, tgt);
        auto btmNode_other = btmPtrs_other[idx_other / kBtmB];
        // {
        //   std::cerr << __FUNCTION__ << ": idx_other = " << idx_other << std::endl;
        //   btmNode_other->printDebugInfo(std::cerr);
        // }
        btmNode_other->increaseW(minSupportW);
        btmNode_other->writeLink(tgtIdxBase + tgt, idx_other % kBtmB);
        // {
        //   std::cerr << __FUNCTION__ << ": idx_other after = " << idx_other << std::endl;
        //   btmNode_other->printDebugInfo(std::cerr);
        // }
      }
    }


    void mvIdxLR
    (
     BtmNode * srcBtmNode,
     BtmNode * tgtBtmNode,
     const uint64_t tgtIdxBase,
     const uint8_t srcIdx,
     const uint8_t tgtIdx,
     const uint8_t num,
     const uint8_t minSupportW,
     BtmNode ** btmPtrs_other
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": tgtIdxBase = " << tgtIdxBase << ", srcIdx = " << (int)srcIdx << ", tgtIdx = " << (int)tgtIdx << ", num = " << (int)num
      //             << ", minSupportW = " << (int)minSupportW << ", kM? = " << (btmPtrs_other == btmPtrs_[kS]) << std::endl;
      // }
      assert(srcIdx + num <= kBtmB);
      assert(tgtIdx + num <= kBtmB);

      for (uint64_t i = 0; i < num; ++i) {
        const uint8_t src = srcIdx + i;
        const uint8_t tgt = tgtIdx + i;
        const uint64_t idx_other = srcBtmNode->readLink(src);
        tgtBtmNode->writeLink(idx_other, tgt);
        auto btmNode_other = btmPtrs_other[idx_other / kBtmB];
        btmNode_other->increaseW(minSupportW);
        btmNode_other->writeLink(tgtIdxBase + tgt, idx_other % kBtmB);
      }
    }


    void makeSpaceInOneBtmNode
    (
     const uint64_t idxBase,
     const uint8_t childIdx,
     const uint64_t * srcWCodes,
     const uint16_t sumW_ins,
     const uint8_t numChild_ins,
     const uint8_t numChild_del, //!< Length of wCodes of tgt to delete.
     const bool insertMorS //!< If "insertMorS == kM", insert to M. If "insertMorS == kS", insert to S.
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": idxBase = " << idxBase << ", childIdx = " << (int)childIdx << ", sumW_ins = " << (int)sumW_ins
      //             << ", numChild_ins = " << (int)numChild_ins << ", numChild_del = " << (int)numChild_del << ", insertMorS = " << insertMorS << std::endl;
      //   // btmPtrs_[insertMorS][idxBase / kBtmB]->printDebugInfo(std::cerr);
      // }

      auto btmNode = btmPtrs_[insertMorS][idxBase / kBtmB];
      btmNode->increaseW(bits::bitSize(size_[!insertMorS] * kBtmB));

      uint16_t sumW_del = 0;
      for (uint8_t i = 0; i < numChild_del; ++i) {
        sumW_del = btmNode->stcc_.readW(childIdx + i);
      }

      const uint8_t tailNum = btmNode->numChildren_ - (childIdx + numChild_del); // at least 0 by assumption.
      uint16_t tailW = 0;
      if (tailNum) {
        const uint8_t minSupportW = bits::bitSize(size_[insertMorS] * kBtmB);
        auto btmPtrs_other = btmPtrs_[!insertMorS];
        tailW = btmNode->stccSize_ - btmNode->calcBitPos(childIdx + numChild_del);
        btmNode->stcc_.mvWCodes(btmNode->stcc_.getConstPtr_wCodes(), childIdx + numChild_del, childIdx + numChild_ins, tailNum);
        if (numChild_ins > numChild_del) {
          mvIdxRL(btmNode, btmNode, idxBase, childIdx + numChild_del, childIdx + numChild_ins, tailNum, minSupportW, btmPtrs_other);
        } else {
          mvIdxLR(btmNode, btmNode, idxBase, childIdx + numChild_del, childIdx + numChild_ins, tailNum, minSupportW, btmPtrs_other);
        }
      }
      btmNode->stcc_.mvWCodes(srcWCodes, 0, childIdx, numChild_ins);
      btmNode->numChildren_ += numChild_ins - numChild_del;
      btmNode->updateWCodesAuxM(childIdx, btmNode->numChildren_);
      if (sumW_ins != sumW_del) {
        const uint16_t newBitSize = btmNode->stccSize_ + sumW_ins - sumW_del;
        btmNode->reserveBitCapacity(newBitSize);
        if (tailNum) {
          btmNode->stcc_.mvVals(btmNode->stcc_.getConstPtr_vals(), btmNode->stccSize_ - tailW, newBitSize - tailW, tailW);
        }
        btmNode->stccSize_ = newBitSize;
      }
    }


    uint64_t overflowToL
    (
     const uint64_t lIdxBase,
     const uint64_t rIdxBase,
     const uint8_t childIdx,
     const uint64_t * srcWCodes,
     const uint8_t numChild_ins,
     const uint8_t numChild_del, //!< Length of wCodes of tgt to delete.
     const bool insertMorS //!< If "insertMorS == kM", insert to M. If "insertMorS == kS", insert to S.
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": lIdxBase = " << lIdxBase << ", rIdxBase = " << rIdxBase << ", childIdx = "
      //             << (int)childIdx << ", numChild_ins = " << (int)numChild_ins << ", numChild_del = " << (int)numChild_del << ", insertMorS = " << insertMorS << std::endl;
      // }
      assert(childIdx + numChild_del <= btmPtrs_[insertMorS][rIdxBase / kBtmB]->getNumChildren());

      auto lnode = btmPtrs_[insertMorS][lIdxBase / kBtmB];
      auto rnode = btmPtrs_[insertMorS][rIdxBase / kBtmB];
      lnode->increaseW(bits::bitSize(size_[!insertMorS] * kBtmB));
      rnode->increaseW(bits::bitSize(size_[!insertMorS] * kBtmB));
      const uint8_t minSupportW = bits::bitSize(size_[insertMorS] * kBtmB);
      auto btmPtrs_other = btmPtrs_[!insertMorS];

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
          mvIdxLR(rnode, lnode, lIdxBase, 0, numL, num, minSupportW, btmPtrs_other);
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
          mvIdxLR(rnode, lnode, lIdxBase, childIdx + numChild_del, numL, curNumAfterDel, minSupportW, btmPtrs_other);
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
        mvIdxLR(rnode, rnode, rIdxBase, numToLeft, 0, num, minSupportW, btmPtrs_other);
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
            mvIdxLR(rnode, rnode, rIdxBase, srcBeg, tgtBeg, num, minSupportW, btmPtrs_other);
          } else {
            mvIdxRL(rnode, rnode, rIdxBase, srcBeg, tgtBeg, num, minSupportW, btmPtrs_other);
          }
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


    uint64_t overflowToR
    (
     const uint64_t lIdxBase,
     const uint64_t rIdxBase,
     const uint8_t childIdx,
     const uint64_t * srcWCodes,
     const uint8_t numChild_ins,
     const uint8_t numChild_del, //!< Length of wCodes of tgt to delete.
     const bool insertMorS //!< If "insertMorS == kM", insert to M. If "insertMorS == kS", insert to S.
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": lIdxBase = " << lIdxBase << ", rIdxBase = " << rIdxBase << ", childIdx = "
      //             << (int)childIdx << ", numChild_ins = " << (int)numChild_ins << ", numChild_del = " << (int)numChild_del << ", insertMorS = " << insertMorS << std::endl;
      // }
      assert(childIdx + numChild_del <= btmPtrs_[insertMorS][lIdxBase / kBtmB]->getNumChildren());

      auto lnode = btmPtrs_[insertMorS][lIdxBase / kBtmB];
      auto rnode = btmPtrs_[insertMorS][rIdxBase / kBtmB];
      lnode->increaseW(bits::bitSize(size_[!insertMorS] * kBtmB));
      rnode->increaseW(bits::bitSize(size_[!insertMorS] * kBtmB));
      const uint8_t minSupportW = bits::bitSize(size_[insertMorS] * kBtmB);
      auto btmPtrs_other = btmPtrs_[!insertMorS];

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
        mvIdxRL(rnode, rnode, rIdxBase, 0, numToRight, numR_old, minSupportW, btmPtrs_other);
      }

      uint8_t numR_increment = 0;
      uint16_t sumWR_increment = 0;
      if (numToRight1) {
        const uint16_t bitPos = lnode->calcBitPos(childIdx - numToRight1);
        const uint16_t w = lnode->calcBitPos(childIdx) - bitPos;
        changeList[clSize++] = {bitPos, 0, w};
        sumWR_increment += w;
        rnode->stcc_.mvWCodes(lnode->getConstPtr_wCodes(), childIdx - numToRight1, 0, numToRight1);
        mvIdxLR(lnode, rnode, rIdxBase, childIdx - numToRight1, 0, numToRight1, minSupportW, btmPtrs_other);
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
        mvIdxLR(lnode, rnode, rIdxBase, numL_old - numToRight2, numR_increment, numToRight2, minSupportW, btmPtrs_other);
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
              mvIdxRL(lnode, lnode, lIdxBase, childIdx + numChild_del, childIdx + numChild_ins, numTail, minSupportW, btmPtrs_other);
            } else {
              mvIdxLR(lnode, lnode, lIdxBase, childIdx + numChild_del, childIdx + numChild_ins, numTail, minSupportW, btmPtrs_other);
            }
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


    void writeNewElemInTwo
    (
     BtmNode * lnode,
     BtmNode * rnode,
     uint8_t childIdx, //!< Relative idx to write counting from left-end of lnode
     const uint64_t * newVals, //!< Storing stcc vals to insert
     const uint64_t * newLinks, //!< Storing new links to insert
     const uint8_t numChild_ins
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": lnode = " << lnode << ", rnode = " << rnode
      //             << ", childIdx = " << (int)childIdx << ", numChild_ins = " << (int)numChild_ins << std::endl;
      // }

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
          lnode->writeLink(newLinks[i - childIdx], i);
          uint8_t w = lnode->stcc_.readW(i);
          lnode->stcc_.writeWBits(newVals[i - childIdx], bitPos, w);
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
          rnode->writeLink(newLinks[i - childIdx + numNewElemInL], i);
          uint8_t w = rnode->stcc_.readW(i);
          // std::cerr << "ci = " << ci
          //           << ", w = " << (int)w
          //           << ", weight = " << weights[ci - childIdx_ins] << std::endl;
          rnode->stcc_.writeWBits(newVals[i - childIdx + numNewElemInL], bitPos, w);
          bitPos += w;
        }
      }
    }


    /*!
     * @brief Insert stcc values
     * @note Weights of BTreeNodes should be changed in advance
     */
    uint64_t insertNewElem
    (
     const uint64_t idxBase,
     const uint8_t childIdx,
     const uint64_t * newVals, //!< Storing stcc vals to insert
     const uint64_t * newLinks, //!< Storing new links to insert
     const uint8_t numChild_ins,
     const uint8_t numChild_del, //!< Length of wCodes of tgt to delete
     const bool insertMorS //!< If "insertMorS == kM", insert to M. If "insertMorS == kS", insert to S.
     ) noexcept {
      assert(numChild_ins <= kBtmB);
      assert(childIdx + numChild_del <= btmPtrs_[insertMorS][idxBase / kBtmB]->getNumChildren()); // could be equal. Especialy "childIdx" could be "numChildren"
      // {//debug
      //   std::cerr << __FUNCTION__ << " idxBase = " << idxBase << ", childIdx = " << (int)childIdx
      //             << ", numChild_ins = " << (int)numChild_ins << ", numChild_del = " << (int)numChild_del
      //             << ", insertMorS = " << (int)insertMorS << std::endl;
      //   // btmPtrs_[insertMorS][idxBase / kBtmB]->printDebugInfo(std::cerr);
      // }

      auto btmNode = btmPtrs_[insertMorS][idxBase / kBtmB];
      auto btmPtrs_other = btmPtrs_[!insertMorS];
      uint64_t wCodesTemp[kBtmB / StepCodeUtil::kWCNum];
      uint16_t sumW_ins = 0;
      for (uint8_t i = 0; i < numChild_ins; ++i) {
        uint8_t w = StepCodeUtil::calcSteppedW(newVals[i]);
        sumW_ins += w;
        StepCodeUtil::writeWCode(StepCodeUtil::calcWCodeFromSteppedW(w), wCodesTemp, i);
      }

      const uint16_t num = static_cast<uint16_t>(btmNode->getNumChildren()) + numChild_ins - numChild_del;
      if (num <= static_cast<uint16_t>(kBtmB)) { // Easy case: This node can accommodate inserting elements.
        makeSpaceInOneBtmNode(idxBase, childIdx, wCodesTemp, sumW_ins, numChild_ins, numChild_del, insertMorS);
        writeNewElemInTwo(btmNode, btmNode, childIdx, newVals, newLinks, numChild_ins);
        return idxBase + childIdx;
      }

      const uint8_t excess = static_cast<uint8_t>(num - kBtmB);
      auto parent = btmNode->getParent();
      const auto idxInSib = btmNode->getIdxInSibling();
      if (idxInSib) { // Check previous sibling.
        auto lnode = reinterpret_cast<BtmNode *>(parent->getChildPtr(idxInSib - 1));
        const auto numL = lnode->getNumChildren();
        if (kBtmB - numL >= excess) { // Previous sibling can accommodate overflowed elements.
          const auto retIdx = overflowToL(calcIdxBase(lnode, btmPtrs_other), idxBase, childIdx, wCodesTemp, numChild_ins, numChild_del, insertMorS);
          writeNewElemInTwo(lnode, btmNode, numL + childIdx, newVals, newLinks, numChild_ins);
          parent->changePSumAt(idxInSib - 1, parent->getPSum(idxInSib) + calcSumOfWeightOfBtmNode(lnode, numL, lnode->getNumChildren(), insertMorS));
          return retIdx;
        }
      }
      if (idxInSib + 1 < parent->getNumChildren()) { // Check next sibling.
        auto rnode = reinterpret_cast<BtmNode *>(parent->getChildPtr(idxInSib + 1));
        const auto numR = rnode->getNumChildren();
        if (kBtmB - numR >= excess) { // Next sibling can accommodate overflowed elements.
          const auto retIdx = overflowToR(idxBase, calcIdxBase(rnode, btmPtrs_other), childIdx, wCodesTemp, numChild_ins, numChild_del, insertMorS);
          writeNewElemInTwo(btmNode, rnode, childIdx, newVals, newLinks, numChild_ins);
          parent->changePSumAt(idxInSib, parent->getPSum(idxInSib + 1) - calcSumOfWeightOfBtmNode(rnode, 0, rnode->getNumChildren() - numR, insertMorS));
          return retIdx;
        }
      }

      { // This bottom node has to be split
        auto rnode = new BtmNode();
        const auto rBtmIdx = setNewBtmNode(rnode, insertMorS);
        const auto retIdx = overflowToR(idxBase, rBtmIdx * kBtmB, childIdx, wCodesTemp, numChild_ins, numChild_del, insertMorS);
        writeNewElemInTwo(btmNode, rnode, childIdx, newVals, newLinks, numChild_ins);
        handleSplitOfBtmInBtm(btmNode, rnode, insertMorS);
        return retIdx;
      }
    }


    uint64_t insertRunAfter_each
    (
     const uint64_t idx,
     const uint64_t val,
     const uint64_t link,
     const bool insertMorS
     ) noexcept {
      // {//debug
      //   std::cerr << __func__ << " idx = " << idx << ", val = " << val << ", insertMorS = " << (int)insertMorS << std::endl;
      // }

      const uint8_t childIdx = (idx % kBtmB) + 1; // +1 is needed to put new run AFTER "idx". "childIdx" could be "kBtmB"
      const uint64_t newVals[] = {val};
      const uint64_t newLinks[] = {link};
      return insertNewElem(idx / kBtmB * kBtmB, childIdx, newVals, newLinks, 1, 0, insertMorS);
    }


    uint64_t insertRunWithSplitM
    (
     const uint64_t idxM,
     const uint64_t splitPos,
     const uint64_t weight
     ) noexcept {
      // {//debug
      //   std::cerr << __func__ << " idxM = " << idxM << ", splitPos = " << splitPos << ", weight = " << weight << std::endl;
      // }

      auto btmNodeM = btmPtrs_[kM][idxM / kBtmB];
      const uint8_t childIdx = idxM % kBtmB;
      changePSumFromParent(btmNodeM, weight);
      const uint64_t weight2 = btmNodeM->readStccVal(childIdx) - splitPos;
      const uint64_t newVals[] = {splitPos, weight, weight2};
      const uint64_t newLinks[] = {btmNodeM->readLink(childIdx), 0, 0}; // 0, 0 are dummy
      return insertNewElem(idxM / kBtmB * kBtmB, childIdx, newVals, newLinks, 3, 1, kM);
    }


  public:
    //////////////////////////////// Public functions (interface)
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
      auto btmNodeM = btmPtrs_[kM][idxM / kBtmB];
      const uint64_t curWeight = btmNodeM->readStccVal(idxM % kBtmB);
      assert(curWeight + change > 0);
      const uint64_t newVals[] = {static_cast<uint64_t>(curWeight + change)};
      btmNodeM->replace(newVals, 1, idxM % kBtmB);
      // update mixed tree
      changePSumFromParent(btmNodeM, change);
      // update separated tree AND alphabet tree (they are connected seamlessly)
      auto btmNodeS = btmPtrs_[kS][idxM2S(idxM) / kBtmB];
      changePSumFromParent(btmNodeS, change);
    }


    /*!
     * @brief Change (increase/decrease) length of run at "idxM".
     */
    void setLeafVal
    (
     const uint64_t idxM,
     const uint64_t newLeafVal
     ) noexcept {
      // update btm node
      const auto idxS = idxM2S(idxM);
      const uint64_t newVals[] = {newLeafVal};
      btmPtrs_[kS][idxS / kBtmB]->replace(newVals, 1, idxS % kBtmB);
    }


    /*!
     * @brief Pushback a run, merging into the last run if possible.
     */
    uint64_t pushbackRun
    (
     const uint64_t ch, //!< 64bit-char.
     const uint64_t weight, //!< Weight (exponent) of new run.
     const uint64_t leafVal,
     uint64_t & pos //!< [out] It is set to relative position of a run.
     ) {
      const auto btmNodeM = reinterpret_cast<BtmNode *>(srootM_.root_->getRmBtm());
      const auto idxS = btmNodeM->readLink(btmNodeM->getNumChildren() - 1);
      if (idxS == 0) { // dummy
        pos = 0;
        return insertRunAfter(0, weight, leafVal, ch);
      }
      const auto btmNodeS = btmPtrs_[kS][idxS / kBtmB];
      const auto idxM = btmNodeS->readLink(idxS % kBtmB);
      if (btmNodeS->getBtmVal() != ch) {
        pos = 0;
        return insertRunAfter(idxM, weight, leafVal, ch);
      } else { // merge into the last run
        pos = getWeightFromIdxM(idxM);
        changeWeight(idxM, weight);
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
     * @brief Insert run of "ch^{weight}" at "pos", merging into adjacent runs if possible.
     */
    uint64_t insertRun
    (
     const uint64_t ch, //!< 64bit-char.
     const uint64_t weight, //!< Weight (exponent) of new run.
     const uint64_t leafVal1,
     const uint64_t leafVal2,
     uint64_t & pos //!< [in,out] 0base position where inserted run will start. It is modified to relative position in a run.
     ) {
      // {//debug
      //   std::cerr << __func__ << ": ch = " << ch << ", weight = " << weight
      //             << ", leafVal1 = " << leafVal1 << ", leafVal2 = " << leafVal2 << ", pos = " << pos << std::endl;
      // }
      if (pos > srootM_.root_->getSumOfWeight()) {
        return BTreeNodeT::NOTFOUND;
      } else if (pos == srootM_.root_->getSumOfWeight()) {
        return pushbackRun(ch, weight, leafVal1, pos);
      }
      auto idxM = searchPosM(pos); // 'pos' is modified to be the relative pos in the run of 'idxM'.
      auto chNow = getCharFromIdxM(idxM);
      if (ch == chNow) {
        changeWeight(idxM, weight);
      } else if (pos == 0) {
        idxM = getPrevIdxM(idxM); // Move to previous idxM.
        if (idxM > 0 && ch == getCharFromIdxM(idxM)) { // Check if 'ch' can be merged with the previous run.
          pos = getWeightFromIdxM(idxM);
          changeWeight(idxM, weight);
        } else {
          idxM = insertRunAfter(idxM, weight, leafVal1, ch);
        }
      } else { // Current run is split with fstHalf of weight 'pos'.
        idxM = insertRunWithSplit(idxM, pos, weight, leafVal1, leafVal2, ch);
        pos = 0;
      }
      return idxM;
    }


    /*!
     * @brief Variant of DynRLE::insertRun for rvalue pos.
     */
    uint64_t insertRun
    (
     const uint64_t ch, //!< 64bit-char.
     const uint64_t weight, //!< Weight (exponent) of new run.
     const uint64_t leafVal1,
     const uint64_t leafVal2,
     uint64_t && pos //!< 0base position where inserted run will start.
     ) {
      auto tmp = pos;
      return insertRun(ch, weight, leafVal1, leafVal2, tmp);
    }


    /*!
     * @brief Insert new run of character 'ch' and length 'weight' after 'idxM'.
     * @return IdxM of the inserted run.
     */
    uint64_t insertRunAfter
    (
     const uint64_t idxM,
     const uint64_t weight,
     const uint64_t leafVal,
     const uint64_t ch
     ) noexcept {
      // {//debug
      //   std::cerr << __func__ << " idxM = " << idxM << ", weight = " << weight << ", leafVal = " << leafVal << std::endl;
      //   std::cerr << __func__ << ": BEFORE idxM = " << idxM << std::endl;
      //   // btmPtrs_[kM][idxM / kBtmB]->printDebugInfo(std::cerr);
      // }

      changePSumFromParent(btmPtrs_[kM][idxM / kBtmB], weight);
      const auto newIdxM = insertRunAfter_each(idxM, weight, 0, kM);
      BTreeNodeT * retRootS = searchCharA(ch);
      uint64_t idxS;
      if (retRootS->isDummy() || getCharFromNodeS(retRootS) != ch) {
        idxS = setupNewSTree(retRootS, ch);
      } else {
        idxS = getPredIdxSFromIdxM(retRootS, ch, newIdxM);
      }
      // {//debug
      //   std::cerr << __FUNCTION__ << ": BEFORE idxS = " << idxS << std::endl;
      //   btmPtrs_[kS][idxS / kBtmB]->printDebugInfo(std::cerr);
      // }
      changePSumFromParent(btmPtrs_[kS][idxS / kBtmB], weight);
      const auto newIdxS = insertRunAfter_each(idxS, leafVal, newIdxM, kS);
      btmPtrs_[kM][newIdxM / kBtmB]->writeLink(newIdxS, newIdxM % kBtmB);
      btmPtrs_[kS][newIdxS / kBtmB]->writeLink(newIdxM, newIdxS % kBtmB);
      // parent->changePSumAt(idxInSib - 1, parent->getPSum(idxInSib) + calcSumOfWeightOfBtmNode(lnode, numL, lnode->getNumChildren(), insertMorS));
      // parent->changePSumAt(idxInSib, parent->getPSum(idxInSib + 1) - calcSumOfWeightOfBtmNode(rnode, 0, rnode->getNumChildren() - numR, insertMorS));

      // {//debug
      //   std::cerr << __FUNCTION__ << ": AFTER idxS = " << idxS << ", newIdxS = " << newIdxS << std::endl;
      //   std::cerr << __FUNCTION__ << ": AFTER show idxS = " << idxS << std::endl;
      //   btmPtrs_[kS][idxS / kBtmB]->printDebugInfo(std::cerr);
      //   if (idxM / kBtmB != newIdxM / kBtmB) {
      //     std::cerr << __FUNCTION__ << ": AFTER show newIdxS = " << newIdxS << std::endl;
      //     btmPtrs_[kS][newIdxS / kBtmB]->printDebugInfo(std::cerr);
      //   }
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


    /*!
     * @brief Insert new run of character 'ch' and length 'weight' splitting run at "idxM".
     * @return IdxM of the inserted run.
     * @note Assume that "ch" is different from the one for the splitted run.
     */
    uint64_t insertRunWithSplit
    (
     const uint64_t idxM,
     const uint64_t splitPos,
     const uint64_t weight,
     const uint64_t leafVal1,
     const uint64_t leafVal2,
     const uint64_t ch
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << " idxM = " << idxM << ", splitPos = " << splitPos
      //             << ", weight = " << weight << ", leafVal1 = " << leafVal1 << ", leafVal2 = " << leafVal2 << std::endl;
      // }

      const uint64_t idxS0 = idxM2S(idxM);
      uint64_t tempIdxM = insertRunWithSplitM(idxM, splitPos, weight);
      if (idxM != tempIdxM) {
        btmPtrs_[kS][idxS0 / kBtmB]->increaseW(bits::bitSize(size_[kM] * kBtmB));
        btmPtrs_[kS][idxS0 / kBtmB]->writeLink(tempIdxM, idxS0 % kBtmB);
      }
      auto * retRootS = searchCharA(ch);
      uint64_t idxS;
      if (retRootS->isDummy() || reinterpret_cast<BtmNode *>(retRootS->getLmBtm_DirectJump())->getBtmVal() != ch) {
        idxS = setupNewSTree(retRootS, ch);
      } else {
        idxS = getPredIdxSFromIdxM(retRootS, ch, tempIdxM);
      }
      const auto newIdxM = getNextIdxM(tempIdxM);
      { // insert new run with character "ch"
        tempIdxM = newIdxM;
        changePSumFromParent(btmPtrs_[kS][idxS / kBtmB], weight);
        idxS = insertRunAfter_each(idxS, leafVal1, tempIdxM, kS);
        btmPtrs_[kM][tempIdxM / kBtmB]->writeLink(idxS, tempIdxM % kBtmB);
      }
      { // insert second half of splitted run
        tempIdxM = getNextIdxM(tempIdxM);
        idxS = insertRunAfter_each(idxS0, leafVal2, tempIdxM, kS);
        btmPtrs_[kM][tempIdxM / kBtmB]->writeLink(idxS, tempIdxM % kBtmB);
      }
      return newIdxM;
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


    size_t calcMemBytesBtmM() const noexcept {
      size_t size = 0;
      for (uint64_t i = 0; i < size_[kM]; ++i) {
        size += btmPtrs_[kM][i]->calcMemBytes();
      }
      return size;
    }


    size_t calcMemBytesBtmS() const noexcept {
      size_t size = 0;
      for (uint64_t i = 0; i < size_[kS]; ++i) {
        size += btmPtrs_[kS][i]->calcMemBytes();
      }
      return size;
    }


    size_t calcMemBytesLinks() const noexcept {
      size_t size = 0;
      for (uint64_t i = 0; i < size_[kM]; ++i) {
        size += btmPtrs_[kM][i]->calcMemBytesLinksArray();
      }
      for (uint64_t i = 0; i < size_[kS]; ++i) {
        size += btmPtrs_[kS][i]->calcMemBytesLinksArray();
      }
      return size;
    }


    size_t calcMemBytesBtmWeights() const noexcept {
      size_t size = 0;
      for (uint64_t i = 0; i < size_[kM]; ++i) {
        size += btmPtrs_[kM][i]->calcMemBytesStccDynArray();
      }
      return size;
    }


    size_t calcMemBytesLeafVals() const noexcept {
      size_t size = 0;
      for (uint64_t i = 0; i < size_[kS]; ++i) {
        size += btmPtrs_[kS][i]->calcMemBytesStccDynArray();
      }
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
      for (uint64_t i = 0; i < size_[kM]; ++i) {
        numUsed += btmPtrs_[kM][i]->getNumChildren();
      }
      return numUsed;
    }


    size_t calcNumSlotsBtmM() const noexcept {
      return size_[kM] * kBtmB;
    }


    size_t calcNumUsedBtmS() const noexcept {
      size_t numUsed = 0;
      for (uint64_t i = 0; i < size_[kS]; ++i) {
        numUsed += btmPtrs_[kS][i]->getNumChildren();
      }
      return numUsed;
    }


    size_t calcNumSlotsBtmS() const noexcept {
      return size_[kS] * kBtmB;
    }


    size_t calcNumRuns() const noexcept {
      size_t numRuns = 0;
      for (size_t i = 0; i < size_[kM]; ++i) {
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
        os << "TotalLen = " << totalLen << ", #Runs = " << numRuns << ", Alphabet Size = " << calcNumAlph()
           << ", BTreeNode arity kB = " << static_cast<uint64_t>(kB) << " BtmNode arity kBtmB = " << static_cast<uint64_t>(kBtmB) << std::endl;
        os << "MTree bottom array size = " << size_[kM] << ", capacity = " << capacity_[kM] << std::endl;
        os << "STree bottom array size = " << size_[kS] << ", capacity = " << capacity_[kS] << std::endl;
        os << "Total: " << calcMemBytes() << " bytes" << std::endl;
        os << "MTree: " << calcMemBytesMTree() << " bytes, OccuRate = " << ((numSlotsM) ? 100.0 * numUsedM / numSlotsM : 0)
           << " (= 100*" << numUsedM << "/" << numSlotsM << ")" << std::endl;
        os << "ATree: " << calcMemBytesATree() << " bytes, OccuRate = " << ((numSlotsA) ? 100.0 * numUsedA / numSlotsA : 0)
           << " (= 100*" << numUsedA << "/" << numSlotsA << ")" << std::endl;
        os << "STree: " << calcMemBytesSTree() << " bytes, OccuRate = " << ((numSlotsS) ? 100.0 * numUsedS / numSlotsS : 0)
           << " (= 100*" << numUsedS << "/" << numSlotsS << ")" << std::endl;
        os << "BtmNodes: " << calcMemBytesBtmM() + calcMemBytesBtmS() << " bytes, OccuRate = " << ((numSlotsBtm) ? 100.0 * numUsedBtm / numSlotsBtm : 0)
           << " (= 100*" << numUsedBtm << "/" << numSlotsBtm << ")" << std::endl;
        os << "Links: " << calcMemBytesLinks() << " bytes" << std::endl;
        os << "Weights: " << calcMemBytesBtmWeights() << " bytes" << std::endl;
        os << "LeafVals: " << calcMemBytesLeafVals() << " bytes" << std::endl;
      }
    }


    void printDebugInfo
    (
     std::ostream & os
     ) const noexcept {
      os << "size_[kM] = " << size_[kM] << ", capacity_[kM] = " << capacity_[kM]
         << ", size_[kS] = " << size_[kS] << ", capacity_[kS] = " << capacity_[kS]
         << ", traCode_ = " << (int)traCode_ << std::endl;
      if (isReady() && getSumOfWeight() > 0) {
        {
          os << "dump btmPtrs_[kM]" << std::endl;
          for (uint64_t i = 0; i < size_[kM]; ++i) {
            os << "[" << i << "]" << btmPtrs_[kM][i] << " ";
          }
          os << std::endl;

          os << "dump btmPtrs_[kS]" << std::endl;
          for (uint64_t i = 0; i < size_[kS]; ++i) {
            os << "[" << i << "]" << btmPtrs_[kS][i] << " ";
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
          for (uint64_t i = 0; i < size_[kM]; ++i) {
            for (uint64_t j = 0; j < getNumChildrenFromBtmM(i); ++j) {
              uint64_t idxM = kBtmB * i + j;
              if (idxM != idxS2M(idxM2S(idxM))) {
                os << "error!! links of idxM2S and idxS2M: idxM = " << idxM
                   << ", idxS = " << idxM2S(idxM) << std::endl; // WARNING, links are not maintained correctly
                btmPtrs_[kM][idxM /kBtmB]->printDebugInfo(std::cerr, true, *this, kM);
                btmPtrs_[kS][idxM2S(idxM) /kBtmB]->printDebugInfo(std::cerr, true, *this, kS);
              }
            }
          }
        }

        { // check links of parent-child for M
          for (uint64_t i = 0; i < size_[kM]; ++i) {
            uint8_t idx = getIdxInSiblingFromBtmM(i);
            auto node = getParentFromBtmM(i);
            bool islmbtm = (idx == 0);
            if (static_cast<void *>(node->getChildPtr(idx)) != btmPtrs_[kM][i]) {
              os << "error!! " << "parent-child for btmM = " << i << std::endl;
            }
            if (islmbtm && static_cast<void *>(node->getLmJumpNode()) != btmPtrs_[kM][i]) {
              os << "error!! lmJumNode for btmM = " << i << std::endl;
            }
            while (!(node->isRoot())) {
              idx = node->getIdxInSibling();
              islmbtm &= (idx == 0);
              if (node->getParent()->getChildPtr(idx) != node) {
                os << "error!! " << "parent-child for child node = " << node << std::endl;
              }
              if (islmbtm && static_cast<void *>(node->getLmJumpNode()) != btmPtrs_[kM][i]) {
                os << "error!! lmJumNode for btmM = " << i << std::endl;
              }
              node = node->getParent();
            }
          }
        }

        { // check links of parent-child for S
          for (uint64_t i = 0; i < size_[kS]; ++i) {
            uint8_t idx = getIdxInSiblingFromBtmS(i);
            auto node = getParentFromBtmS(i);
            bool islmbtm = (idx == 0);
            if (static_cast<void *>(node->getChildPtr(idx)) != btmPtrs_[kS][i]) {
              os << "error!! " << "parent-child for btmS = " << i << std::endl;
            }
            if (islmbtm && static_cast<void *>(node->getLmJumpNode()) != btmPtrs_[kS][i]) {
              os << "error!! lmJumpNode for btmS = " << i << std::endl;
            }
            while (!(node->isRoot())) {
              idx = node->getIdxInSibling();
              islmbtm &= (idx == 0);
              if (static_cast<void *>(node->getParent()->getChildPtr(idx)) != node) {
                os << "error!! " << "parent-child for child node = " << node << std::endl;
              }
              if (islmbtm && static_cast<void *>(node->getLmJumpNode()) != btmPtrs_[kS][i]) {
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

        // {
        //   os << "Information on M" << std::endl;
        //   uint64_t pos = 0;
        //   for (auto btmNodeM = reinterpret_cast<const BtmNode *>(srootM_.root_->getLmBtm_DirectJump());
        //        btmNodeM != nullptr;
        //        btmNodeM = btmNodeM->getNextBtmNode()) {
        //     btmNodeM->printDebugInfo(std::cerr, true, *this, kM);
        //   }
        //   os << std::endl;
        // }

        // {
        //   os << "Alphabet: " << std::endl;
        //   for (const auto * rootS = getFstRootS();
        //        reinterpret_cast<uintptr_t>(rootS) != BTreeNodeT::NOTFOUND;
        //        rootS = getNextRootS(rootS)) {
        //     os << "(" << getCharFromNodeS(rootS) << ", " << rootS->getSumOfWeight() << ") ";
        //   }
        //   os << std::endl;
        //   os << std::endl;
        // }

        // {
        //   os << "Information on S" << std::endl;
        //   for (const auto * rootS = getFstRootS();
        //        reinterpret_cast<uintptr_t>(rootS) != BTreeNodeT::NOTFOUND;
        //        rootS = getNextRootS(rootS)) {
        //     for (auto btmNodeS = reinterpret_cast<const BtmNode *>(rootS->getLmBtm_DirectJump());
        //          btmNodeS != nullptr;
        //          btmNodeS = btmNodeS->getNextBtmNode()) {
        //       btmNodeS->printDebugInfo(std::cerr, true, *this, kS);
        //     }
        //   }
        // }

      }
    }
  };
} // namespace itmmti

#endif

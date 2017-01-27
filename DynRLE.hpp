/**
 * @file DynRLE.hpp
 * @brief Dynamic run-length encoding supporting access, rank, select, and insert (TODO delete).
 * @author Tomohiro I
 * @date 2017-01-15
 */
#ifndef INCLUDE_GUARD_DynRLE
#define INCLUDE_GUARD_DynRLE

#include <stdint.h>

#include <iostream>
#include <algorithm>
#include <cassert>
#include <string>
#include <fstream>
#include <sstream>

#include "../../Basics/BitsUtil.hpp"
#include "../../Basics/WBitsArray.hpp"
#include "../../Basics/MemUtil.hpp"
#include "BTree.hpp"

#define XT 1.5  //ラベル振り直しの基準値
#define X 20
uint64_t numRelabel = 0;
uint64_t numInRelabel = 0;


template <uint8_t B> class BTreeUpperNode;

template <uint8_t B = 64> // B should be in {4, 8, 16, 32, 64, 128}. B/2 <= 'numChildren_' <= B
class DynRLE
{
  // mixed tree
  BTreeUpperNode<B> * rootM_;
  // information for leaves and elements for mixed tree
  WBitsArray idxM2S_;
  BTreeUpperNode<B> ** parentM_;
  uint64_t * labelM_;
  uint8_t * idxInSiblingM_;
  WBitsArray ** weightArrayVec_;
  // alphabet tree: the bottoms of alphabet tree are roots of separated trees
  BTreeUpperNode<B> * rootA_;
  // information for leaves and elements for separated tree
  WBitsArray idxS2M_;
  BTreeUpperNode<B> ** parentS_;
  uint64_t * charS_;
  uint8_t * idxInSiblingS_;
  uint8_t * numChildrenS_;

  // uint8_t labelingParam_;

public:
  DynRLE() :
    rootM_(NULL),
    idxM2S_(),
    parentM_(NULL),
    labelM_(NULL),
    idxInSiblingM_(NULL),
    weightArrayVec_(NULL),
    rootA_(NULL),
    idxS2M_(),
    parentS_(NULL),
    charS_(NULL),
    idxInSiblingS_(NULL),
    numChildrenS_(NULL)
    // labelingParam_(0)
  {}


  DynRLE(const size_t initNumBtms) {
    init(initNumBtms);
  }


  ~DynRLE() {
    clearAll();
  }


  void init(const size_t initNumBtms) {
    assert(initNumBtms > 0);
    if (isReady()) {
      clearAll();
    }
    reserveBtms(initNumBtms);

    rootM_ = new BTreeUpperNode<B>(true, true, reinterpret_cast<BTreeUpperNode<B> *>(0));
    // sentinel
    parentM_[0] = rootM_;
    idxInSiblingM_[0] = 0;
    labelM_[0] = 0;
    idxM2S_.resize(B);
    idxM2S_.write(0, 0); // sentinel: should not be used
    weightArrayVec_[0] = new WBitsArray(8, B);
    weightArrayVec_[0]->resize(1);
    weightArrayVec_[0]->write(0, 0);
    rootM_->pushbackBtm(reinterpret_cast<BTreeUpperNode<B> *>(0), 0);

    auto * dummyRootS = new BTreeUpperNode<B>(true, true, NULL, true);
    dummyRootS->pushbackBtm(NULL, 0);
    rootA_ = new BTreeUpperNode<B>(true, true, dummyRootS);
    rootA_->pushbackUNode(dummyRootS);
  }


  void clearAll() {
    if (!isReady()) { // already cleared
      return;
    }
    for (uint64_t i = 0; i < idxM2S_.size() / B; ++i) {
      delete weightArrayVec_[i];
    }
    memutil::myfree(weightArrayVec_);
    { // delete separated tree
      auto * rootS = rootA_->getLmBtm();
      while (reinterpret_cast<uintptr_t>(rootS) != BTreeUpperNode<B>::NOTFOUND) {
        auto * next = getNextRootS(rootS);
        delete rootS;
        rootS = next;
      }
    }
    idxM2S_.clear();
    idxM2S_.shrink_to_fit();
    idxS2M_.clear();
    idxS2M_.shrink_to_fit();
    memutil::myfree(parentM_);
    memutil::myfree(parentS_);
    memutil::myfree(labelM_);
    memutil::myfree(charS_);
    memutil::myfree(idxInSiblingM_);
    memutil::myfree(idxInSiblingS_);
    memutil::myfree(numChildrenS_);

    memutil::mydelete(rootM_);
    memutil::mydelete(rootA_);
  }


  bool isReady() const noexcept {
    return (rootM_ != NULL);
  }


  size_t getSumOfWeight() const noexcept {
    return rootM_->getSumOfWeight();
  }


  size_t getSumOfWeight(const uint64_t ch) const noexcept {
    const auto * retRootS = searchCharA(ch);
    if (retRootS->isDummy() || charS_[reinterpret_cast<uintptr_t>(retRootS->getLmBtm())] != ch) {
      return 0;
    }
    return retRootS->getSumOfWeight();
  }


  uint64_t getWeightFromIdxM(uint64_t idxM) const noexcept {
    return weightArrayVec_[idxM / B]->read(idxM % B);
  }


  uint64_t getCharFromIdxM(const uint64_t idxM) const noexcept {
    return charS_[idxM2S_.read(idxM) / B];
  }


  uint64_t getCharFromNodeS(const BTreeUpperNode<B> * nodeS) const noexcept {
    return charS_[reinterpret_cast<uintptr_t>(nodeS->getLmBtm())];
  }


  uint64_t rank(const uint64_t ch, uint64_t pos, const bool calcTotalRank) const noexcept {
    assert(pos < rootM_->getSumOfWeight());
    auto idxM = searchPosM(pos); // pos is modified to relative pos
    return rank(ch, idxM, pos, calcTotalRank);
  }


  uint64_t rank(const uint64_t ch, const uint64_t idxM, const uint64_t relativePos, const bool calcTotalRank) const noexcept {
    auto chNow = getCharFromIdxM(idxM);
    uint64_t ret = 0;
    uint64_t idxS;
    if (ch == chNow) {
      ret = relativePos + 1;
      idxS = idxM2S_.read(idxM);
    } else {
      const auto * retRootS = searchCharA(ch);
      if (retRootS->isDummy() || getCharFromNodeS(retRootS) != ch) {
        return 0;
      }
      idxS = getPredIdxSFromIdxM(retRootS, ch, idxM);
    }
    const auto btmS = idxS / B;
    for (auto tmpIdxS = btmS * B; tmpIdxS < idxS + (ch != chNow); ++tmpIdxS) {
      ret += getWeightFromIdxS(tmpIdxS);
    }
    return ret + parentS_[btmS]->calcPSum(idxInSiblingS_[btmS], calcTotalRank);
  }


  /**
   * NOTE: rank is 1base
   * @return the smallest pos (0base) s.t. rank_{ch}[0..pos] (pos inclusive).
   */
  uint64_t select(const BTreeUpperNode<B> * rootS, const uint64_t rank) const noexcept {
    assert(rank > 0);
    auto pos = rank - 1;
    const auto idxS = searchPosS(pos, rootS); // pos is modified to the relative pos
    const auto idxM = idxS2M_.read(idxS);
    const auto btmM = idxM / B;
    for (auto tmpIdxM = btmM * B; tmpIdxM < idxM; ++tmpIdxM) {
      pos += getWeightFromIdxM(tmpIdxM);
    }
    return pos + parentM_[btmM]->calcPSum(idxInSiblingM_[btmM], false);
  }


  uint64_t select(const uint64_t ch, const uint64_t rank) const noexcept {
    assert(rank > 0);
    const auto * retRootS = searchCharA(ch);
    if (retRootS->isDummy() || getCharFromNodeS(retRootS) != ch) {
      return BTreeUpperNode<B>::NOTFOUND;
    }
    if (rank > retRootS->getSumOfWeight()) {
      return BTreeUpperNode<B>::NOTFOUND;
    }
    return select(retRootS, rank);
  }


  uint64_t select(const uint64_t totalRank) const noexcept {
    assert(totalRank > 0);
    if (totalRank > rootA_->getSumOfWeight()) {
      return BTreeUpperNode<B>::NOTFOUND;
    }
    auto pos = totalRank - 1;
    const auto * retRootS = searchPosA(pos);
    return select(retRootS, pos + 1); // +1 for 1base rank
  }


  void printString(std::ofstream & ofs) const noexcept {
    uint64_t pos = 0;
    for (auto idxM = searchPosM(pos); idxM != BTreeUpperNode<B>::NOTFOUND; idxM = getNextIdxM(idxM)) {
      const size_t exponent = getWeightFromIdxM(idxM);
      char ch = getCharFromIdxM(idxM);
      for (size_t i = 0; i < exponent; ++i) {
        ofs.put(ch);
      }
    }
  }


  //// public search functions
public:
  /**
     Return pointer ('idxM') to the run containing 'sum'-th character (0base) in RLE.
     sum is modified to be the relative position (0base) from the beginning of the run.
   */
  uint64_t searchPosM(uint64_t & pos) const noexcept {
    assert(pos < rootM_->getSumOfWeight());
    // if (sum >= rootM_->getSumOfWeight()) {
    //   return BTreeUpperNode<B>::NOTFOUND;
    // }
    uint64_t btmM = reinterpret_cast<uintptr_t>(rootM_->searchPos(pos));

    const auto * wArray = weightArrayVec_[btmM];
    uint8_t i = 0;
    while (pos >= wArray->read(i)) {
      pos -= wArray->read(i);
      ++i;
    }
    return btmM * B + i;
  }


  /**
     search root of separated tree of the largest character that is smaller or equal to 'ch'.
  */
  BTreeUpperNode<B> * searchCharA(const uint64_t ch) const noexcept {
    auto * nodeA = rootA_;
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


  //// private search functions
private:
  uint64_t searchPosS(uint64_t & pos, const BTreeUpperNode<B> * rootS) const noexcept {
    // check this outside
    // if (pos >= rootS->getSumOfWeight()) {
    //   return BTreeUpperNode<B>::NOTFOUND;
    // }
    uint64_t idxS = B * reinterpret_cast<uintptr_t>(rootS->searchPos(pos));

    while (true) {
      auto weight = getWeightFromIdxS(idxS);
      if (pos >= weight) {
        pos -= weight;
        ++idxS;
      } else {
        return idxS;
      }
    }
  }


  /**
     search idxS with the largest label that is smaller or equal to 'label'
  */
  uint64_t searchLabelS(const uint64_t label, const BTreeUpperNode<B> * rootS) const noexcept {
    const auto * nodeS = rootS;
    while (true) {
      const bool nowOnBorder = nodeS->isBorder();
      uint8_t lb = 0;
      uint8_t ub = nodeS->getNumChildren();
      while (lb+1 != ub) {
        uint8_t mid = (lb + ub) / 2;
        if (label < getLabelFromNodeU(nodeS->getChildPtr(mid), nowOnBorder)) {
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
    const uint64_t idxS = B * reinterpret_cast<uintptr_t>(nodeS);
    uint8_t lb = 0;
    uint8_t ub = numChildrenS_[idxS / B];
    while (lb+1 != ub) {
      uint8_t mid = (lb + ub) / 2;
      if (label < labelM_[idxS2M_.read(idxS + mid) / B]) {
        ub = mid;
      } else {
        lb = mid;
      }
    }
    return idxS + lb;
  }


  //// iterator like functions
public:
  uint64_t getPrevIdxM(const uint64_t idxM) const noexcept {
    if (idxM % B) {
      return idxM - 1;
    }
    const uint64_t prevBtmM
      = reinterpret_cast<uintptr_t>(parentM_[idxM / B]->getPrevBtm(idxInSiblingM_[idxM / B]));
    if (prevBtmM != BTreeUpperNode<B>::NOTFOUND) {
      return prevBtmM * B + getNumChildrenM(prevBtmM) - 1;
    }
    return BTreeUpperNode<B>::NOTFOUND;
  }


  uint64_t getNextIdxM(const uint64_t idxM) const noexcept {
    if ((idxM % B) + 1 < getNumChildrenM(idxM / B)) {
      return idxM + 1;
    }
    const uint64_t nextBtmM
      = reinterpret_cast<uintptr_t>(parentM_[idxM / B]->getNextBtm(idxInSiblingM_[idxM / B]));
    if (nextBtmM != BTreeUpperNode<B>::NOTFOUND) {
      return nextBtmM * B;
    }
    return BTreeUpperNode<B>::NOTFOUND;
  }


  BTreeUpperNode<B> * getFstRootS() const noexcept {
    return getNextRootS(rootA_->getLmBtm());
  }


  BTreeUpperNode<B> * getPrevRootS(const BTreeUpperNode<B> * node) const noexcept {
    uint8_t idxInSib;
    do {
      idxInSib = node->getIdxInSibling();
      node = node->getParent();
    } while (!node->isBorder());
    return node->getPrevBtm(idxInSib);
  }


  BTreeUpperNode<B> * getNextRootS(const BTreeUpperNode<B> * node) const noexcept {
    uint8_t idxInSib;
    do {
      idxInSib = node->getIdxInSibling();
      node = node->getParent();
    } while (!(node->isBorder()));
    return node->getNextBtm(idxInSib);
  }


private:
  uint64_t getPrevIdxS(const size_t idxS) const noexcept {
    if (idxS % B) {
      return idxS - 1;
    }
    const uint64_t prevBtmS
      = reinterpret_cast<uintptr_t>(parentS_[idxS / B]->getPrevBtm(idxInSiblingS_[idxS / B]));
    if (prevBtmS != BTreeUpperNode<B>::NOTFOUND) {
      return prevBtmS * B + numChildrenS_[prevBtmS] - 1;
    }
    return BTreeUpperNode<B>::NOTFOUND;
  }


  uint64_t getNextIdxS(const size_t idxS) const noexcept {
    if ((idxS % B) + 1 < numChildrenS_[idxS / B]) {
      return idxS + 1;
    }
    const uint64_t nextBtmS
      = reinterpret_cast<uintptr_t>(parentS_[idxS / B]->getNextBtm(idxInSiblingS_[idxS / B]));
    if (nextBtmS != BTreeUpperNode<B>::NOTFOUND) {
      return nextBtmS * B;
    }
    return BTreeUpperNode<B>::NOTFOUND;
  }


  //// private getter functions (utilities)
private:
  uint8_t getNumChildrenM(const uint64_t btmM) const noexcept {
    return weightArrayVec_[btmM]->size();
  }


  uint64_t getWeightFromIdxS(uint64_t idxS) const noexcept {
    const uint64_t idxM = idxS2M_.read(idxS);
    return weightArrayVec_[idxM / B]->read(idxM % B);
  }


  uint64_t getLabelFromNodeU(const BTreeUpperNode<B> * nodeU, const bool isChildOfBorder)  const noexcept {
    uint64_t idxS;
    if (!isChildOfBorder) {
      idxS = B * reinterpret_cast<uintptr_t>(nodeU->getLmBtm());
    } else {
      idxS = B * reinterpret_cast<uintptr_t>(nodeU);
    }
    uint64_t idxM = idxS2M_.read(idxS); // idxM corresponding to the left most idxS in btmS
    return labelM_[idxM / B];
  }


  uint64_t getCharFromNodeA(const BTreeUpperNode<B> * nodeA, const bool isChildOfBorder) const noexcept {
    uint64_t btmS;
    if (!isChildOfBorder) {
      btmS = reinterpret_cast<uintptr_t>(nodeA->getLmBtm()->getLmBtm());
    } else {
      btmS = reinterpret_cast<uintptr_t>(nodeA->getLmBtm());
    }
    return charS_[btmS];
  }


  uint64_t getPrevBtmM(const uint64_t btmM) const noexcept {
    return reinterpret_cast<uintptr_t>(parentM_[btmM]->getPrevBtm(idxInSiblingM_[btmM]));
  }


  uint64_t getNextBtmM(const uint64_t btmM) const noexcept {
    //    debugstream2 << __LINE__ << ": getNextBtmM: btmM " << btmM << ", parentM_[btmM] " << parentM_[btmM] << std::endl;
    return reinterpret_cast<uintptr_t>(parentM_[btmM]->getNextBtm(idxInSiblingM_[btmM]));
  }


  /**
     returns root of separated tree that contains the position 'pos' (0based) in alphabetically sorted array
  */
  BTreeUpperNode<B> * searchPosA(uint64_t & pos) const noexcept {
    return rootA_->searchPos(pos);
  }


  void reserveBtms(const size_t numBtms) {
    const uint8_t w = bits::bitSize(numBtms * B - 1);
    idxM2S_.convert(w, numBtms * B);
    idxS2M_.convert(w, numBtms * B);
    memutil::myrealloc(weightArrayVec_, numBtms);
    memutil::myrealloc(parentM_, numBtms);
    memutil::myrealloc(parentS_, numBtms);
    memutil::myrealloc(labelM_, numBtms);
    memutil::myrealloc(charS_, numBtms);
    memutil::myrealloc(idxInSiblingM_, numBtms);
    memutil::myrealloc(idxInSiblingS_, numBtms);
    memutil::myrealloc(numChildrenS_, numBtms);
  }


  void expandBtms() {
    const uint64_t newNumBtms = 2 * (idxM2S_.capacity() / B); // number of capacity of bottoms is doubled
    reserveBtms(newNumBtms);
  }


  ////
  void changeWeight(const uint64_t idxM, const int64_t change) {
    // update Leaf
    auto * wArray = weightArrayVec_[idxM / B];
    const uint64_t val = wArray->read(idxM % B) + change;
    const uint8_t w = wArray->getW();
    const uint8_t needW = bits::bitSize(val);
    if (needW > w) {
      wArray->convert(needW, B);
    }
    wArray->write(val, idxM % B);
    // update mixed tree
    parentM_[idxM / B]->changePSumFrom(idxInSiblingM_[idxM / B], change);
    // update separated tree AND alphabet tree (they are connected seamlessly)
    const uint64_t btmS = idxM2S_.read(idxM) / B;
    parentS_[btmS]->changePSumFrom(idxInSiblingS_[btmS], change);
  }


  void asgnLabel(const uint64_t btmM) {
    uint64_t next = getNextBtmM(btmM);
    uint64_t prev = getPrevBtmM(btmM);
    if (next == BTreeUpperNode<B>::NOTFOUND) {
      labelM_[btmM] = labelM_[prev] + X;
      return;
    } else if (labelM_[prev] < labelM_[next] - 1){
      labelM_[btmM] = (labelM_[prev] + labelM_[next]) / 2;
      return;
    }
    ++numInRelabel; // 4test

    double criterion = 1;
    uint64_t criterionLabel = labelM_[next];
    uint64_t num = 2;
    uint64_t btmM0 = btmM;
    uint8_t l = 0;
    next = getNextBtmM(next);
    do {
      assert(l < 63);
      ++l;
      criterion *= (2/XT); // *= T giwaku
      criterionLabel >>= 1;
      while (prev != BTreeUpperNode<B>::NOTFOUND && (labelM_[prev] >> l) == criterionLabel) { // expand backward
        ++num;
        btmM0 = prev;
        prev = getPrevBtmM(prev);
      }
      while (next != BTreeUpperNode<B>::NOTFOUND && (labelM_[next] >> l) == criterionLabel){ // expand forward
        ++num;
        next = getNextBtmM(next);
      }
    } while (criterion <= num);

    numRelabel += num; // 4test
    // relable num labels
    uint64_t tLabel = criterionLabel << l;
    uint64_t interval = (1ULL << l) / num;
    while (true) {
      labelM_[btmM0] = tLabel;
      if (--num == 0) {
        return;
      }
      tLabel += interval;
      btmM0 = getNextBtmM(btmM0);
    }
  }


  /**
     setup
       parentM_[retBtmM] (by handleSplitBtmM())
       idxInSiblingM_[retBtmM] (by handleSplitBtmM())
       labelM_[retBtmM] (by asgnLabel())
     resize
       idxM2S_ to use range [endIdxM, endIdxM + B)
     reserve
       weightArrayVec_[retBtmM]
     update
       upper nodes (through handleSplitBtmM())
       labels
  */
  uint64_t splitBtmM(const uint8_t width, const uint64_t btmM, const uint64_t weight) {
    const uint64_t endIdxM = idxM2S_.size();
    const uint64_t retBtmM = endIdxM / B;
    if (!(idxM2S_.resizeWithoutReserve(endIdxM + B))) {
      expandBtms();
    }
    idxM2S_.resize(endIdxM + B);
    // reserve
    weightArrayVec_[retBtmM] = new WBitsArray(width, B);
    // setup and update
    handleSplitOfBtmInBtm(btmM, retBtmM, weight, parentM_, idxInSiblingM_);
    if (!(rootM_->isRoot())) { // root needs update
      rootM_ = rootM_->getParent();
    }
    asgnLabel(retBtmM);
    return retBtmM;
  }


  /**
     setup
       parentS_[retBtmS] (by handleSplitBtmS())
       idxInSiblingS_[retBtmS] (by handleSplitBtmS())
       charS_[retBtmS]
     resize
       idxS2M_ to use range [endIdxS, endIdxS + B)
     update
       upper nodes (through handleSplitBtmS())
  */
  uint64_t splitBtmS(const uint64_t btmS, const uint64_t weight) {
    const uint64_t endIdxS = idxS2M_.size();
    const uint64_t retBtmS = endIdxS / B;
    if (!(idxS2M_.resizeWithoutReserve(endIdxS + B))) {
      expandBtms();
    }
    idxS2M_.resize(endIdxS + B);
    // setup and update
    handleSplitOfBtmInBtm(btmS, retBtmS, weight, parentS_, idxInSiblingS_);
    charS_[retBtmS] = charS_[btmS];
    return retBtmS;
  }


  uint64_t setupNewSTree(BTreeUpperNode<B> * predNode, const uint64_t ch) {
    const uint64_t endIdxS = idxS2M_.size();
    const uint64_t btmS = endIdxS / B;
    if (!(idxS2M_.resizeWithoutReserve(endIdxS + B))) {
      expandBtms();
    }
    idxS2M_.resize(endIdxS + B);
    
    auto * newRootS = new BTreeUpperNode<B>(true, true, reinterpret_cast<BTreeUpperNode<B> *>(btmS));
    parentS_[btmS] = newRootS;
    idxInSiblingS_[btmS] = 0;
    charS_[btmS] = ch;
    numChildrenS_[btmS] = 1; // only dummy idxS exists
    idxS2M_.write(0, btmS * B); // link to dummy idxM of weight 0

    newRootS->pushbackBtm(reinterpret_cast<BTreeUpperNode<B> *>(btmS), 0);
    const auto idxInSib = predNode->getIdxInSibling();
    auto * parent = predNode->getParent();
    parent->handleSplitOfChild(newRootS, idxInSib);
    if (!(rootA_->isRoot())) {
      rootA_ = rootA_->getParent();
    }
    return endIdxS;
  }


  void mvIdxFwd(WBitsArray & wba, uint64_t srcIdx, uint64_t tgtIdx, uint64_t num, WBitsArray & wbaOther) {
    for (uint64_t i = num; i > 0; --i) {
      const uint64_t idxOther = wba.read(srcIdx + i - 1);
      wba.write(idxOther, tgtIdx + i - 1);
      wbaOther.write(tgtIdx + i - 1, idxOther);
    }
  }


  uint64_t makeSpaceAfterIdxM(const uint64_t idxM) {
    const uint8_t remIdxM = idxM % B;
    const uint64_t btmM = idxM / B;
    WBitsArray * wArray0 = weightArrayVec_[btmM];
    const uint8_t oriNum = wArray0->size();
    if (oriNum < B) {
      wArray0->resize(oriNum + 1);
      const uint8_t mvNum = oriNum - remIdxM - 1;
      if (mvNum) {
        mvWBA_SameW(wArray0->getItrAt(remIdxM + 1), wArray0->getItrAt(remIdxM + 2), mvNum);
        mvIdxFwd(idxM2S_, idxM + 1, idxM + 2, mvNum, idxS2M_);
      }
      wArray0->write(0, remIdxM + 1);
      return idxM + 1;
    }
    // split
    uint64_t sum = 0;
    for (uint8_t i = B/2; i < B; ++i) {
      sum += wArray0->read(i);
    }
    const auto newBtmM = splitBtmM(wArray0->getW(), btmM, sum);
    WBitsArray * wArray1 = weightArrayVec_[newBtmM];
    mvWBA_SameW(wArray0->getItrAt(B/2), wArray1->getItrAt(0), B/2);
    wArray0->resize(B/2);
    wArray1->resize(B/2);
    mvIdxFwd(idxM2S_, btmM*B + B/2, newBtmM*B, B/2, idxS2M_);
    if (remIdxM < B/2) {
      return makeSpaceAfterIdxM(idxM);
    } else {
      return makeSpaceAfterIdxM(newBtmM * B + remIdxM - B/2);
    }
  }


  uint64_t makeSpaceAfterIdxS(const uint64_t idxS) {
    const uint8_t remIdxS = idxS % B;
    const uint64_t btmS = idxS / B;
    const uint8_t oriNum = numChildrenS_[btmS];
    if (oriNum < B) {
      numChildrenS_[btmS] = oriNum + 1;
      const uint8_t mvNum = oriNum - remIdxS - 1;
      if (mvNum) {
        mvIdxFwd(idxS2M_, idxS + 1, idxS + 2, mvNum, idxM2S_);
      }
      return idxS + 1;
    }
    uint64_t sum = 0;
    for (uint8_t i = B/2; i < B; ++i) {
      sum += getWeightFromIdxS(btmS * B + i);
    }
    const auto newBtmS = splitBtmS(btmS, sum);
    numChildrenS_[btmS] = B/2;
    numChildrenS_[newBtmS] = B/2;
    mvIdxFwd(idxS2M_, btmS*B + B/2, newBtmS*B, B/2, idxM2S_);
    if (remIdxS < B/2) {
      return makeSpaceAfterIdxS(idxS);
    } else {
      return makeSpaceAfterIdxS(newBtmS * B + remIdxS - B/2);
    }
  }


  void handleSplitOfBtmInBtm(const uint64_t btm, const uint64_t newBtm, const uint64_t weight, BTreeUpperNode<B> ** parentArray, uint8_t * idxInSibArray) {
    auto * uNode = parentArray[btm];
    const auto idxInSib = idxInSibArray[btm];
    const auto oriNum = uNode->getNumChildren();
    BTreeUpperNode<B> * newNode = uNode->handleSplitOfBtm(reinterpret_cast<BTreeUpperNode<B> *>(newBtm), weight, idxInSib);
    const uint8_t newNum = (oriNum < B) ? oriNum + 1 : B/2 + (idxInSib < B/2);
    // update links to upper nodes
    for (uint8_t i = idxInSib + 1; i < newNum; ++i) {
      const uint64_t tmp = reinterpret_cast<uintptr_t>(uNode->getChildPtr(i));
      parentArray[tmp] = uNode;
      idxInSibArray[tmp] = i;
    }
    if (oriNum == B) { // split
      for (uint8_t i = 0; i < B/2 + (idxInSib >= B/2); ++i) { // for children of newNode
        const uint64_t tmp = reinterpret_cast<uintptr_t>(newNode->getChildPtr(i));
        parentArray[tmp] = newNode;
        idxInSibArray[tmp] = i;
      }
    }
  }


  uint64_t getPredIdxSFromIdxM(const BTreeUpperNode<B> * rootS, const uint64_t ch, const uint64_t idxM) const noexcept {
    const uint64_t btmM = idxM / B;
    if (btmM) { // if btmM is not 0 (0 means btmM is the first btm in the mixed tree)
      uint64_t i = idxM - 1;
      for ( ; i >= btmM * B && getCharFromIdxM(i) != ch; --i) {}
      if (i >= btmM * B) {
        return idxM2S_.read(i);
      } else {
        return searchLabelS(labelM_[btmM] - 1, rootS); // -1 is needed
      }
    } else { // btmM == 0: dummy idx (== 0) should be ignored.
      uint64_t i = idxM - 1;
      for ( ; i > 0 && getCharFromIdxM(i) != ch; --i) {}
      if (i > 0) {
        return idxM2S_.read(i);
      } else {
        return B * reinterpret_cast<uintptr_t>(rootS->getLmBtm());
      }
    }
  }


  /**
   * insert new run of character 'ch' and length 'weight' after 'idxM'.
   * return idx to the inserted run.
   */
  uint64_t insertNewRunAfter(const uint64_t ch, const uint64_t weight, const uint64_t idxM) {
    const auto newIdxM = makeSpaceAfterIdxM(idxM);
    auto * retRootS = searchCharA(ch);
    uint64_t idxS;
    if (retRootS->isDummy() || charS_[reinterpret_cast<uintptr_t>(retRootS->getLmBtm())] != ch) {
      idxS = setupNewSTree(retRootS, ch);
    } else {
      idxS = getPredIdxSFromIdxM(retRootS, ch, newIdxM);
    }
    const auto newIdxS = makeSpaceAfterIdxS(idxS);
    idxM2S_.write(newIdxS, newIdxM);
    idxS2M_.write(newIdxM, newIdxS);
    changeWeight(newIdxM, weight);
    return newIdxM;
  }


public:
  uint64_t pushbackRun(const uint64_t ch, const uint64_t weight) {
    const uint64_t btmM = reinterpret_cast<uintptr_t>(rootM_->getRmBtm());
    return insertNewRunAfter(ch, weight, btmM * B + getNumChildrenM(btmM) - 1);
  }


  ////
  uint64_t insertRun(const uint64_t ch, const uint64_t weight, uint64_t & pos) {
    if (pos > rootM_->getSumOfWeight()) {
      return BTreeUpperNode<B>::NOTFOUND;
    } else if (pos == rootM_->getSumOfWeight()) {
      pos = 0;
      return pushbackRun(ch, weight);
    }
    auto idxM = searchPosM(pos); // 'pos' is modified to be the relative pos in the run of 'idxM'
    auto chNow = getCharFromIdxM(idxM);
    if (ch == chNow) {
      changeWeight(idxM, weight);
    } else if (pos == 0) {
      idxM = getPrevIdxM(idxM); // move to previous idxM
      if (idxM > 0 && ch == getCharFromIdxM(idxM)) { // check if 'ch' can be merged with the previous run
        pos = getWeightFromIdxM(idxM);
        changeWeight(idxM, weight);
      } else {
        idxM = insertNewRunAfter(ch, weight, idxM);
      }
    } else { // current run is split with fstHalf of weight 'pos'
      const auto weightSndHalf = getWeightFromIdxM(idxM) - pos;
      pos = 0;
      changeWeight(idxM, -1 * weightSndHalf);
      idxM = insertNewRunAfter(ch, weight, idxM);
      insertNewRunAfter(chNow, weightSndHalf, idxM);
    }
    return idxM;
  }


  void insertRunWithoutReturn(const uint64_t ch, const uint64_t weight, const uint64_t pos) {
    auto tmp = pos;
    insertRun(ch, weight, tmp);
  }


  //// statistics
public:
  size_t calcMemBytesMTree() const noexcept {
    return rootM_->calcMemBytes();
  }


  size_t calcMemBytesATree() const noexcept {
    return rootA_->calcMemBytes();
  }


  size_t calcMemBytesSTree() const noexcept {
    size_t size = 0;
    for (const auto * rootS = getFstRootS();
         reinterpret_cast<uintptr_t>(rootS) != BTreeUpperNode<B>::NOTFOUND;
         rootS = getNextRootS(rootS)) {
      size += rootS->calcMemBytes();
    }
    return size;
  }


  size_t calcMemBytesWeightArrays() const noexcept {
    size_t size = 0;
    for (uint64_t i = 0; i < idxM2S_.size() / B; ++i) {
      size += weightArrayVec_[i]->calcMemBytes();
    }
    return size;
  }


  size_t calcMemBytesIdxConvertArrays() const noexcept {
    size_t size = 0;
    size += idxM2S_.calcMemBytes();
    size += idxS2M_.calcMemBytes();
    return size;
  }


  size_t calcMemBytesBtmArrays() const noexcept {
    return (idxM2S_.capacity() / B) * (sizeof(parentM_[0]) + sizeof(parentS_[0]) +
                                       sizeof(labelM_[0]) + sizeof(charS_[0]) +
                                       sizeof(idxInSiblingM_[0]) + sizeof(idxInSiblingS_[0]) +
                                       sizeof(numChildrenS_[0]) + sizeof(weightArrayVec_[0]));
  }


  size_t calcMemBytes() const noexcept {
    size_t size = sizeof(*this);
    size += calcMemBytesMTree();
    size += calcMemBytesATree();
    size += calcMemBytesSTree();
    size += calcMemBytesWeightArrays();
    size += calcMemBytesIdxConvertArrays();
    size -= sizeof(idxM2S_); // minus double counted part
    size -= sizeof(idxS2M_); // minus double counted part
    size += calcMemBytesBtmArrays();
    return size;
  }


  size_t calcNumUsedSTree() const noexcept {
    size_t numUsed = 0;
    for (const auto * rootS = getFstRootS();
         reinterpret_cast<uintptr_t>(rootS) != BTreeUpperNode<B>::NOTFOUND;
         rootS = getNextRootS(rootS)) {
      numUsed += rootS->calcNumUsed();
    }
    return numUsed;
  }


  size_t calcNumSlotsSTree() const noexcept {
    size_t numSlots = 0;
    for (const auto * rootS = getFstRootS();
         reinterpret_cast<uintptr_t>(rootS) != BTreeUpperNode<B>::NOTFOUND;
         rootS = getNextRootS(rootS)) {
      numSlots += rootS->calcNumSlots();
    }
    return numSlots;
  }


  size_t calcNumRuns() const noexcept {
    size_t numRuns = 0;
    for (size_t i = 0; i < idxM2S_.size() / B; ++i) {
      numRuns += getNumChildrenM(i);
    }
    return numRuns - 1; // -1 due to the first dummy
  }


  size_t calcNumAlph() const noexcept {
    size_t numAlph = 0;
    for (const auto * rootS = getFstRootS();
         reinterpret_cast<uintptr_t>(rootS) != BTreeUpperNode<B>::NOTFOUND;
         rootS = getNextRootS(rootS)) {
      ++numAlph;
    }
    return numAlph;
  }


  void printStatictics(std::ostream & os) const noexcept {
    const size_t totalLen = getSumOfWeight();
    const size_t numRuns = calcNumRuns();
    os << "TotalLen = " << totalLen << ", #Runs = " << numRuns << ", Alphabet Size = " << calcNumAlph() << ", BTree arity param B = " << static_cast<int>(B) << std::endl;
    os << "Total: " << calcMemBytes() << " bytes" << std::endl;
    os << "MTree: " << calcMemBytesMTree() << " bytes, OccuRate = " << ((rootM_->calcNumSlots()) ? 100.0 * rootM_->calcNumUsed() / rootM_->calcNumSlots() : 0)
       << " (= 100*" << rootM_->calcNumUsed() << "/" << rootM_->calcNumSlots() << ")" << std::endl;
    os << "ATree: " << calcMemBytesATree() << " bytes, OccuRate = " << ((rootA_->calcNumSlots()) ? 100.0 * rootA_->calcNumUsed() / rootA_->calcNumSlots() : 0)
       << " (= 100*" << rootA_->calcNumUsed() << "/" << rootA_->calcNumSlots() << ")" << std::endl;
    os << "STree: " << calcMemBytesSTree() << " bytes, OccuRate = " << ((calcNumSlotsSTree()) ? 100.0 * calcNumUsedSTree() / calcNumSlotsSTree() : 0)
       << " (= 100*" << calcNumUsedSTree() << "/" << calcNumSlotsSTree() << ")" << std::endl;
    os << "IdxConverArrays: " << calcMemBytesIdxConvertArrays() << " bytes ~ "
       << "(2*" << static_cast<int>(idxM2S_.getW()) << "(bitwidth)*" << idxM2S_.capacity() << "(capacity each))/8, "
       << "OccuRate = " << ((idxM2S_.capacity() + idxS2M_.capacity()) ? 100.0 * 2 * numRuns / (idxM2S_.capacity() + idxS2M_.capacity()) : 0)
       << " (= 100*2*" << numRuns << "/" << (idxM2S_.capacity() + idxS2M_.capacity()) << ")" << std::endl;
    os << "WeightArrays: " << calcMemBytesWeightArrays() << " bytes" << std::endl;
    os << "BtmArrays: " << calcMemBytesBtmArrays() << " bytes, "
       << "OccuRate = " << ((idxM2S_.capacity() + idxS2M_.capacity()) ? 100.0 * (idxM2S_.size() + idxS2M_.size()) / (idxM2S_.capacity() + idxS2M_.capacity()) : 0)
       << " (= 100*" << (idxM2S_.size() + idxS2M_.size())/B << "/" << (idxM2S_.capacity() + idxS2M_.capacity())/B << "), "
       << "OccuRate (btmM) = " << ((idxM2S_.capacity()) ? 100.0 * idxM2S_.size() / idxM2S_.capacity() : 0)
       << " (= 100*" << idxM2S_.size()/B << "/" << idxM2S_.capacity()/B << "), "
       << "OccuRate (btmS) = " << ((idxS2M_.capacity()) ? 100.0 * idxS2M_.size() / idxS2M_.capacity() : 0)
       << " (= 100*" << idxS2M_.size()/B << "/" << idxS2M_.capacity()/B << ")" << std::endl;
    // 4test
    os << "numInRelabel = " << numInRelabel << ", numRelabel = " << numRelabel << ", ratio = " << 1.0 * numRelabel / numInRelabel << std::endl;
  }


  void printDebugInfo(std::ostream & os) const noexcept {
    {
      uint64_t pos = 0;
      for (auto idxM = searchPosM(pos); idxM != BTreeUpperNode<B>::NOTFOUND; idxM = getNextIdxM(idxM)) {
        os << "(" << idxM << ":" << getCharFromIdxM(idxM) << "^" << getWeightFromIdxM(idxM) << ") ";
      }
      os << std::endl;
    }

    {
      const uint64_t numBtmM = idxM2S_.size() / B;
      os << "information on M" << std::endl;
      for (uint64_t i = 0; i < numBtmM; ++i) {
        const auto nextBtmM = getNextBtmM(i);
        os << "[" << i*B << "-" << (i+1)*B-1 << "] (num=" << (int)getNumChildrenM(i) << " lbl=" 
                  << labelM_[i] << " par=" << parentM_[i] << " sib=" << (int)idxInSiblingM_[i] << ") "
                  << "=> " << nextBtmM * B << std::endl;
        for (uint64_t j = 0; j < B; ++j) {
          if (j < getNumChildrenM(i) && B*i+j != idxS2M_.read(idxM2S_.read(B*i+j))) {
            os << "!!"; // WARNING, links are not maintained correctly
          }
          os << idxM2S_.read(B*i+j) << "(" << getWeightFromIdxM(B*i+j) << ")  ";
        }
        os << std::endl;
      }
    }

    {
      const uint64_t numBtmS = idxS2M_.size() / B;
      os << "information on S" << std::endl;
      for (uint64_t i = 0; i < numBtmS; ++i) {
        const auto nextIdxS = getNextIdxS(i*B + numChildrenS_[i] - 1);
        os << "[" << i*B << "-" << (i+1)*B-1 << "] (num=" << (int)numChildrenS_[i] << " ch=" << charS_[i] << " par=" 
                  << parentS_[i] << " sib=" << (int)idxInSiblingS_[i] << ") "
                  << "=> " << nextIdxS << std::endl;
        for (uint64_t j = 0; j < B; ++j) {
          os << idxS2M_.read(B*i+j) << "  ";
        }
        os << std::endl;
      }
    }

    os << "Alphabet: " << std::endl;
    for (const auto * rootS = getFstRootS();
         reinterpret_cast<uintptr_t>(rootS) != BTreeUpperNode<B>::NOTFOUND;
         rootS = getNextRootS(rootS)) {
      const uint64_t btmS = reinterpret_cast<uintptr_t>(rootS->getLmBtm());
      os << "(" << charS_[btmS] << ", " << rootS->getSumOfWeight() << ") ";
    }
    os << std::endl;
  }
};

#endif

/**
 * @file BTree.hpp
 * @author Tomohiro I
 * @date 2017-01-26
 */
#ifndef INCLUDE_GUARD_BTree
#define INCLUDE_GUARD_BTree

#include <algorithm>
#include <cassert>



template <uint8_t B = 64> // B should be in {4, 8, 16, 32, 64, 128}. B/2 <= 'numChildren_' <= B
class BTreeNode
{
  uint64_t psum_[B]; // partial sum: psum_[i] = sum_{i = 0}^{i} [weight of i-th child (0base)]
  BTreeNode<B> * parent_;
  uint8_t idxInSibling_; // this node is the 'idxInSibling_'-th child (0base) of its parent
  uint8_t numChildren_;
  uint8_t flags_;
  BTreeNode<B> * children_[B];
  BTreeNode<B> * lmBtm_;

  enum { // for 'flags_'
    isBorderBit = 1,
    isRootBit = 2,
    isDummyBit = 4,
  };

public:
  enum {
    NOTFOUND = UINTPTR_MAX
  };

public:
  BTreeNode<B>(bool isBorder, bool isRoot, BTreeNode<B> * lmBtm, bool isDummy = false)
  : parent_(NULL),
    numChildren_(0),
    flags_(isBorder * isBorderBit | isRoot * isRootBit | isDummy * isDummyBit),
    lmBtm_(lmBtm)
  {}
  ~BTreeNode<B>() = default;
  BTreeNode<B>(const BTreeNode<B> &) = delete;
  BTreeNode<B> & operator=(const BTreeNode<B> &) = delete;


  //// simple getter
public:
  uint64_t getPSum(uint8_t i) const noexcept {
    assert(i < numChildren_);
    return psum_[i];
  }


  uint64_t getWeightOfChild(uint8_t i) const noexcept {
    assert(i < numChildren_);
    return (i > 0) ? psum_[i] - psum_[i-1] : psum_[0];
  }


  uint64_t getSumOfWeight() const noexcept {
    assert(numChildren_ > 0);
    return psum_[numChildren_ - 1];
  }


  BTreeNode<B> * getChildPtr(uint8_t i) const noexcept {
    return children_[i];
  }


  BTreeNode<B> * getParent() const noexcept {
    return parent_;
  }


  uint8_t getIdxInSibling() const noexcept {
    return idxInSibling_;
  }


  uint8_t getNumChildren() const noexcept {
    return numChildren_;
  }


  bool isBorder() const noexcept {
    return flags_ & isBorderBit;
  }


  bool isRoot() const noexcept {
    return flags_ & isRootBit;
  }


  bool isDummy() const noexcept {
    return flags_ & isDummyBit;
  }


  //// function to traverse tree
  BTreeNode<B> * getLmBtm() const noexcept {
    return lmBtm_;
  }


  BTreeNode<B> * getRmBtm() const noexcept {
    const auto * node = this;
    while (!(node->isBorder())) {
      node = node->children_[node->getNumChildren() - 1];
    }
    return node->children_[node->getNumChildren() - 1];
  }


  BTreeNode<B> * getNextBtm(uint8_t idxInSib) const noexcept {
    const auto * node = this;
    while (idxInSib + 1 == node->getNumChildren() && !(node->isRoot())) {
      idxInSib = node->getIdxInSibling();
      node = node->getParent();
    }
    if (idxInSib + 1 < node->getNumChildren()) {
      if (node->isBorder()) {
        return node->getChildPtr(idxInSib + 1);
      } else {
        return node->getChildPtr(idxInSib + 1)->getLmBtm();
      }
    }
    return reinterpret_cast<BTreeNode<B> *>(NOTFOUND);
  }


  BTreeNode<B> * getPrevBtm(uint8_t idxInSib) const noexcept {
    const auto * node = this;
    while (idxInSib == 0 && !(node->isRoot())) {
      idxInSib = node->getIdxInSibling();
      node = node->getParent();
    }
    if (idxInSib) {
      if (node->isBorder()) {
        return node->getChildPtr(idxInSib - 1);
      } else {
        return node->getChildPtr(idxInSib - 1)->getRmBtm();
      }
    }
    return reinterpret_cast<BTreeNode<B> *>(NOTFOUND);
  }


  /**
   * @param calcTotalPSum stopping criteria: If false, goes up until node 'isRoot() == true'. If true, goes up until node with 'parent_ == NULL' (convenient when BTrees are stacked).
   * @return partial sum up to the node (exclusive) indicated by idx-th child (0base) of this node
   */
  uint64_t calcPSum(uint8_t idx, const bool calcTotalPSum) const noexcept {
    assert(isBorder());
    const auto * node = this;
    uint64_t ret = 0;
    while (true) {
      ret += (idx > 0) ? node->getPSum(idx - 1) : 0;
      if (node->getParent() == NULL || (!calcTotalPSum && node->isRoot())) {
        return ret;
      }
      idx = node->getIdxInSibling();
      node = node->getParent();
    }
  }


  ////
  BTreeNode<B> * searchPos(uint64_t & pos) const noexcept {
    const auto * ptr = std::upper_bound(psum_, psum_ + numChildren_, pos);
    const uint8_t i = ptr - psum_;
    if (i) {
      pos -= *(--ptr);
    }
    if (isBorder()) {
      return children_[i];
    }
    return children_[i]->searchPos(pos);
  }


  //// private modifier (intend to call from member function of BTreeNode)
private:
  void setLmBtm(BTreeNode<B> * btmRoot) noexcept {
    lmBtm_ = btmRoot;
  }


  void updateLmBtm(BTreeNode<B> * btmRoot) noexcept {
    auto * node = this;
    while (true) {
      node->setLmBtm(btmRoot);
      if (node->isRoot() || node->getIdxInSibling() > 0) {
        break;
      }
      node = node->getParent();
    }
  }


  void unroot() noexcept {
    flags_ &= ~isRootBit;
  }


  void setParentRef(BTreeNode<B> * newParent, uint8_t newIdxInSibling) noexcept {
    this->parent_ = newParent;
    this->idxInSibling_ = newIdxInSibling;
  }


  void setChildPtr(BTreeNode * child, uint8_t idx) noexcept {
    assert(idx < numChildren_);
    children_[idx] = child;
  }


  void makeNewRoot(BTreeNode<B> * fstHalf, BTreeNode<B> * sndHalf) {
    auto newRoot = new BTreeNode<B>(false, true, fstHalf->getLmBtm());
    auto * parent = fstHalf->getParent();
    if (parent != NULL) { // BTrees are stacked
      const auto idxInSib = fstHalf->getIdxInSibling();
      parent->setChildPtr(newRoot, idxInSib); // parent points to newRoot
      newRoot->setParentRef(parent, idxInSib); // newRoot points to parent
      if (idxInSib == 0) {
        parent->updateLmBtm(newRoot);
      }
    }
    newRoot->pushbackUNode(fstHalf);
    newRoot->pushbackUNode(sndHalf);
  }


  //// public modifier
public:
  void pushbackUNode(BTreeNode<B> * child) noexcept {
    assert(numChildren_ < B);
    children_[numChildren_] = child;
    psum_[numChildren_] = (numChildren_ > 0) ?
      psum_[numChildren_ - 1] + child->getSumOfWeight() : child->getSumOfWeight();
    child->setParentRef(this, numChildren_);
    ++numChildren_;
  }


  void pushbackBtm(BTreeNode<B> * child, const uint64_t psumVal) noexcept {
    assert(isBorder());
    assert(numChildren_ < B);
    children_[numChildren_] = child;
    psum_[numChildren_] = psumVal;
    ++numChildren_;
  }


  //// 'children_[idx]' is split to 'children_[idx]' and 'sndHalf'
  void handleSplitOfChild(BTreeNode<B> * sndHalf, const uint8_t idx) {
    // assert(!isBorder()); // If node is on border and a child is not treated as BTreeNode<B> type, it should be processed differently
    const auto end = numChildren_;
    assert(idx <= end);
    if (end < B) {
      numChildren_ = idx;
      this->pushbackUNode(children_[idx]);
      auto * pushC = sndHalf;
      for (uint8_t i = idx+1; i <= end; ++i) {
        auto * tmp = children_[i];
        this->pushbackUNode(pushC);
        pushC = tmp;
      }
      return;
    }
    // this node has to be split
    auto * lmBtm = (this->isBorder()) ? children_[B/2] : children_[B/2]->getLmBtm();
    auto newNode = new BTreeNode<B>(this->isBorder(), false, lmBtm);
    for (uint8_t i = B/2; i < B; ++i) {
      newNode->pushbackUNode(children_[i]);
    }
    numChildren_ = B/2;
    if (idx < B/2) {
      this->handleSplitOfChild(sndHalf, idx);
    } else {
      newNode->handleSplitOfChild(sndHalf, idx - B/2);
    }
    if (!isRoot()) {
      parent_->handleSplitOfChild(newNode, idxInSibling_);
    } else {
      this->unroot();
      makeNewRoot(this, newNode);
    }
  }


  BTreeNode<B> * handleSplitOfBtm(BTreeNode<B> * sndHalf, const uint64_t weight, const uint8_t idx) {
    assert(isBorder());
    assert(idx <= numChildren_);
    if (numChildren_ < B) {
      for (uint8_t i = numChildren_; idx + 1 < i; --i) {
        children_[i] = children_[i-1];
        psum_[i] = psum_[i-1];
      }
      ++numChildren_;
      children_[idx+1] = sndHalf;
      psum_[idx+1] = psum_[idx];
      psum_[idx] -= weight;
      return NULL;
    }
    // this node has to be split
    auto newNode = new BTreeNode(true, false, children_[B/2]);
    const auto minus = psum_[B/2 - 1];
    for (uint8_t i = B/2; i < B; ++i) {
      newNode->pushbackBtm(children_[i], psum_[i] - minus);
    }
    numChildren_ = B/2;
    if (idx < B/2) {
      this->handleSplitOfBtm(sndHalf, weight, idx);
    } else {
      newNode->handleSplitOfBtm(sndHalf, weight, idx - B/2);
    }
    if (!isRoot()) {
      parent_->handleSplitOfChild(newNode, idxInSibling_);
    } else {
      this->unroot();
      makeNewRoot(this, newNode);
    }
    return newNode;
  }


  void changePSumFrom(const uint8_t idx, const int64_t change) noexcept {
    for (uint8_t i = idx; i < numChildren_; ++i) {
      psum_[i] += change;
    }
    if (parent_ != NULL) { // we do not use isRoot() here for convenience. That is, when we stack two or more BTrees, the change will be propagated.
      parent_->changePSumFrom(idxInSibling_, change);
    }
  }


  //// calculate statistics
public:
  size_t calcMemBytes() const noexcept {
    size_t sumOfSize = sizeof(*this);
    if (!isBorder()) {
      for (uint8_t i = 0; i < numChildren_; ++i) {
        sumOfSize += children_[i]->calcMemBytes();
      }
    }
    return sumOfSize;
  }


  size_t calcNumUsed() const noexcept {
    size_t numOfUsed = numChildren_;
    if (!isBorder()) {
      for (uint8_t i = 0; i < numChildren_; ++i) {
        numOfUsed += children_[i]->calcNumUsed();
      }
    }
    return numOfUsed;
  }


  size_t calcNumSlots() const noexcept {
    size_t numOfSlots = B;
    if (!isBorder()) {
      for (uint8_t i = 0; i < numChildren_; ++i) {
        numOfSlots += children_[i]->calcNumSlots();
      }
    }
    return numOfSlots;
  }
};

#endif

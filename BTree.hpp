/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file BTree.hpp
 * @brief Pointer-based implementation of upper part of B+tree.
 * @attention Bottom part of B+tree can be implemented in more space efficient way (e.g., array-based implementation).
 * @author Tomohiro I
 * @date 2017-01-26
 * @todo Support delete.
 */
#ifndef INCLUDE_GUARD_BTree
#define INCLUDE_GUARD_BTree

#include <algorithm>
#include <cassert>


/*!
 * @brief Pointer-based implementation of upper part of B+tree.
 * @attention Bottom part of B+tree can be implemented in more space efficient way (e.g., array-based implementation).
 * @tparam B should be in {4, 8, 16, 32, 64, 128}. B/2 <= 'numChildren_' <= B.
 */
template <uint8_t B = 64>
class BTreeNode
{
  uint64_t psum_[B]; //!< Partial sum: psum_[i] = sum_{i = 0}^{i} [weight of i-th child (0base)]
  BTreeNode<B> * parent_; //!< Pointer to parent node.
  uint8_t idxInSibling_; //!< This node is 'idxInSibling_'-th child (0base) of its parent.
  uint8_t numChildren_; //!< Current num of children.
  uint8_t flags_;
  BTreeNode<B> * children_[B]; //!< Pointers to children. Note that the children might not BTreeNode type when this node is in border.
  BTreeNode<B> * lmBtm_; //!< To remember leftmost bottom of this node.

  /*!
   * @brief For representing current status of node by 'flags_'.
   */
  enum {
    isBorderBit = 1, //!< This node lies in border, i.e., its children are bottom nodes.
    isRootBit = 2, //!< This node is a root.
    isDummyBit = 4, //!< This node is a dummy node.
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
  /*!
   * @brief Get the sum of weights of child subtrees numbered from 0 to i.
   */
  uint64_t getPSum
  (
   uint8_t i //!< in [0, 'numChildren_').
   ) const noexcept {
    assert(i < numChildren_);

    return psum_[i];
  }


  /*!
   * @brief Get the weight of a child subtree.
   */
  uint64_t getWeightOfChild
  (
   uint8_t i //!< in [0, 'numChildren_').
   ) const noexcept {
    assert(i < numChildren_);

    return (i > 0) ? psum_[i] - psum_[i-1] : psum_[0];
  }


  /*!
   * @brief Get the weight of this node.
   * @pre There is at least one child.
   */
  uint64_t getSumOfWeight() const noexcept {
    assert(numChildren_ > 0);

    return psum_[numChildren_ - 1];
  }


  /*!
   * @brief Get pointer to the i-th child (0base).
   */
  BTreeNode<B> * getChildPtr
  (
   uint8_t i //!< in [0, 'numChildren_').
   ) const noexcept {
    assert(i < numChildren_);

    return children_[i];
  }


  /*!
   * @brief Get pointer to parent.
   */
  BTreeNode<B> * getParent() const noexcept {
    return parent_;
  }


  /*!
   * @brief Get BTreeNode<B>::idxInSibling_.
   */
  uint8_t getIdxInSibling() const noexcept {
    return idxInSibling_;
  }


  /*!
   * @brief Get num of children.
   */
  uint8_t getNumChildren() const noexcept {
    return numChildren_;
  }


  /*!
   * @brief Return if this node is in border.
   */
  bool isBorder() const noexcept {
    return flags_ & isBorderBit;
  }


  /*!
   * @brief Return if this node is root.
   */
  bool isRoot() const noexcept {
    return flags_ & isRootBit;
  }


  /*!
   * @brief Return if this node is a dummy node.
   */
  bool isDummy() const noexcept {
    return flags_ & isDummyBit;
  }


  //// function to traverse tree
  /*!
   * @brief Get leftmost bottom. It takes O(1) time.
   */
  BTreeNode<B> * getLmBtm() const noexcept {
    return lmBtm_;
  }


  /*!
   * @brief Get leftmost bottom. It takes O(h) time, where h is the height of this node.
   */
  BTreeNode<B> * getRmBtm() const noexcept {
    const auto * node = this;
    while (!(node->isBorder())) {
      node = node->children_[node->getNumChildren() - 1];
    }
    return node->children_[node->getNumChildren() - 1];
  }


  /*!
   * @brief Return next bottm starting from 'idxInSib'-th child (0base) of this node.
   */
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


  /*!
   * @brief Return previous bottm starting from 'idxInSib'-th child (0base) of this node.
   */
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


  /*!
   * @brief Return partial sum up to the node (exclusive) indicated by 'idx'-th child (0base) of this node.
   */
  uint64_t calcPSum
  (
   uint8_t idx,
   const bool calcTotalPSum //!< Stopping criterion: If false, goes up until node 'isRoot() == true'. If true, goes up until node with 'parent_ == NULL' (convenient when BTrees are stacked).
   ) const noexcept {
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
  /*!
   * @brief Traverse tree looking for 'pos'.
   * @return Pointer to bottom node where partial sum of 'pos' is achieved, where weight-0 nodes (e.g. dummy nodes) are skipped.
   */
  BTreeNode<B> * searchPos
  (
   uint64_t & pos //!< [in,out] Give global position to search. It is modified to relative position in bottom node.
   ) const noexcept {
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
    newRoot->pushbackBTreeNode(fstHalf);
    newRoot->pushbackBTreeNode(sndHalf);
  }


  //// public modifier
public:
  /*!
   * @brief Pushback node of ::BTreeNode type.
   * @note
   *   Set psum_, fix links between this node and child node, and increment 'numChildren_'.
   */
  void pushbackBTreeNode(BTreeNode<B> * child) noexcept {
    assert(numChildren_ < B);

    children_[numChildren_] = child;
    psum_[numChildren_] = (numChildren_ > 0) ?
      psum_[numChildren_ - 1] + child->getSumOfWeight() : child->getSumOfWeight();
    child->setParentRef(this, numChildren_);
    ++numChildren_;
  }


  /*!
   * @brief Pushback bottom node other than ::BTreeNode type.
   * @note
   *   Set psum_ to given value, set link from this node to child node, and increment 'numChildren_'.
   *   Link from child node to this node should be set outside this function.
   */
  void pushbackBtm(BTreeNode<B> * child, const uint64_t psumVal) noexcept {
    assert(isBorder());
    assert(numChildren_ < B);

    children_[numChildren_] = child;
    psum_[numChildren_] = psumVal;
    ++numChildren_;
  }


  /*!
   * @brief Handle the situation where 'children_[idx]' is split to 'children_[idx]' and 'sndHalf' when child node is ::BTreeNode type.
   */
  void handleSplitOfChild(BTreeNode<B> * sndHalf, const uint8_t idx) {
    assert(idx <= numChildren_);

    const auto end = numChildren_;
    if (end < B) {
      numChildren_ = idx;
      this->pushbackBTreeNode(children_[idx]);
      auto * pushC = sndHalf;
      for (uint8_t i = idx+1; i <= end; ++i) {
        auto * tmp = children_[i];
        this->pushbackBTreeNode(pushC);
        pushC = tmp;
      }
      return;
    }
    // this node has to be split
    auto * lmBtm = (this->isBorder()) ? children_[B/2] : children_[B/2]->getLmBtm();
    auto newNode = new BTreeNode<B>(this->isBorder(), false, lmBtm);
    for (uint8_t i = B/2; i < B; ++i) {
      newNode->pushbackBTreeNode(children_[i]);
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


  /*!
   * @brief Handle the situation where 'children_[idx]' is split to 'children_[idx]' and 'sndHalf' when child node is not ::BTreeNode type.
   */
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


  /*!
   * @brief Change weight of 'idx'-th child of this node, and accordingly all 'psum_' values needed to be fixed.
   */
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

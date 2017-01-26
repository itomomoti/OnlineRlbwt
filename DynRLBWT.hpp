#ifndef INCLUDE_GUARD_DynRLBWT
#define INCLUDE_GUARD_DynRLBWT

#include <stdint.h>
#include <ostream>


/*
 * In contrast to DynRLE, DynRLWT has a vertial end mark at insertPos_
*/
template <class DynRLE>
class DynRLBWT
{
  // mixed tree
  DynRLE drle_;
  uint64_t emPos_;
  uint64_t em_; // end marker


public:
  DynRLBWT(const size_t initNumBtms, uint64_t em = UINT64_MAX)
    : drle_(initNumBtms),
      emPos_(0),
      em_(em)
  {}


  uint64_t getEm() const noexcept {
    return em_;
  }


  uint64_t getEmPos() const noexcept {
    return emPos_;
  }


  void extend(uint64_t ch) {
    uint64_t idxM = drle_.insertRun(ch, 1, emPos_);
    emPos_ = drle_.rank(ch, idxM, emPos_, true);
  }


  uint64_t operator[](uint64_t pos) const noexcept {
    if (pos == emPos_) {
      return em_;
    } else if (pos > emPos_) {
      --pos;
    }
    uint64_t idxM = drle_.searchPosM(pos);
    return drle_.getCharFromIdxM(idxM);
  }


  uint64_t getLenWithEm() const noexcept {
    return drle_.getSumOfWeight() + 1;
  }


  uint64_t totalRank(uint64_t ch, uint64_t pos) const noexcept {
    if (pos == emPos_) {
      return em_;
    } else if (pos > emPos_) {
      --pos;
    }
    return drle_.rank(ch, pos, true);
  }


  void printStatictics(std::ostream & os) const noexcept {
    drle_.printStatictics(os);
  }


  size_t calcMemBytes() const noexcept {
    return sizeof(*this) + drle_.calcMemBytes();
  }
};

#endif

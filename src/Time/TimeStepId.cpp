// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Time/TimeStepId.hpp"

#include <boost/functional/hash.hpp>
#include <ostream>
#include <pup.h>

#include "Time/EvolutionOrdering.hpp"
#include "Time/Slab.hpp"
#include "Time/Time.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"

UnsizedTimeStepId::UnsizedTimeStepId(const bool time_runs_forward,
                                     const int64_t slab_number,
                                     const Time& step_time,
                                     const uint64_t substep)
    : time_runs_forward_(time_runs_forward),
      slab_number_(slab_number),
      step_time_(step_time),
      substep_(substep) {
  if (time_runs_forward_ ? step_time_.is_at_slab_end()
                         : step_time_.is_at_slab_start()) {
    ASSERT(substep_ == 0,
           "Time needs to be advanced, but step already started");
    const Slab new_slab = time_runs_forward_ ? step_time_.slab().advance()
                                             : step_time_.slab().retreat();
    ++slab_number_;
    step_time_ = step_time_.with_slab(new_slab);
  }
}

bool UnsizedTimeStepId::is_at_slab_boundary() const {
  return substep_ == 0 and step_time_.is_at_slab_boundary();
}

void UnsizedTimeStepId::pup(PUP::er& p) {
  p | time_runs_forward_;
  p | slab_number_;
  p | step_time_;
  p | substep_;
}

TimeStepId::TimeStepId(const bool time_runs_forward, const int64_t slab_number,
                       const Time& time)
    : UnsizedTimeStepId(time_runs_forward, slab_number, time),
      substep_time_(time.value()) {}

TimeStepId::TimeStepId(const bool time_runs_forward, const int64_t slab_number,
                       const Time& step_time, const uint64_t substep,
                       const TimeDelta& step_size, const double substep_time)
    : UnsizedTimeStepId(time_runs_forward, slab_number, step_time, substep),
      step_size_(step_size),
      substep_time_(substep_time) {
  ASSERT(substep_ != 0 or step_time_.value() == substep_time_,
         "Initial substep must align with the step.");
  if (time_runs_forward_) {
    ASSERT(substep_time_ >= step_time_.value(),
           "Substep must be within the step.");
  } else {
    ASSERT(substep_time_ <= step_time_.value(),
           "Substep must be within the step.");
  }
}

const TimeDelta& TimeStepId::step_size() const {
  ASSERT(substep_ != 0, "Step size not available at substep 0.");
  return step_size_;
}

TimeStepId TimeStepId::next_step(const TimeDelta& step_size) const {
  ASSERT(substep_ == 0 or step_size == step_size_,
         "Step size inconsistent: " << step_size_ << " " << step_size);
  return TimeStepId(time_runs_forward_, slab_number_, step_time_ + step_size);
}

TimeStepId TimeStepId::next_substep(const TimeDelta& step_size,
                                    const double step_fraction) const {
  ASSERT(substep_ == 0 or step_size == step_size_,
         "Step size inconsistent: " << step_size_ << " " << step_size);
  ASSERT(step_fraction >= 0.0 and step_fraction <= 1.0,
         "Substep must be within the step.");
  const double new_time = (1.0 - step_fraction) * step_time_.value() +
                          step_fraction * (step_time_ + step_size).value();
  return TimeStepId(time_runs_forward_, slab_number_, step_time_, substep_ + 1,
                    step_size, new_time);
}

void TimeStepId::pup(PUP::er& p) {
  UnsizedTimeStepId::pup(p);
  p | step_size_;
  p | substep_time_;
}

bool operator==(const UnsizedTimeStepId& a, const UnsizedTimeStepId& b) {
  ASSERT(a.time_runs_forward() == b.time_runs_forward(),
         "Time is not running in a consistent direction");
  return a.slab_number() == b.slab_number() and
         a.step_time() == b.step_time() and a.substep() == b.substep();
}

bool operator!=(const UnsizedTimeStepId& a, const UnsizedTimeStepId& b) {
  return not(a == b);
}

bool operator<(const UnsizedTimeStepId& a, const UnsizedTimeStepId& b) {
  ASSERT(a.time_runs_forward() == b.time_runs_forward(),
         "Time is not running in a consistent direction");
  if (a.slab_number() != b.slab_number()) {
    return a.slab_number() < b.slab_number();
  }
  if (a.step_time() != b.step_time()) {
    return evolution_less<Time>{a.time_runs_forward()}(a.step_time(),
                                                       b.step_time());
  }
  return a.substep() < b.substep();
}
bool operator<=(const UnsizedTimeStepId& a, const UnsizedTimeStepId& b) {
  return not(b < a);
}
bool operator>(const UnsizedTimeStepId& a, const UnsizedTimeStepId& b) {
  return b < a;
}
bool operator>=(const UnsizedTimeStepId& a, const UnsizedTimeStepId& b) {
  return not(a < b);
}

bool operator==(const TimeStepId& a, const TimeStepId& b) {
  bool equal = static_cast<const UnsizedTimeStepId&>(a) ==
               static_cast<const UnsizedTimeStepId&>(b);
  if (equal and a.substep() != 0) {
    equal = a.step_size() == b.step_size();
  }
  ASSERT(not equal or a.substep_time() == b.substep_time(),
         "IDs at same step and substep but different times");
  return equal;
}

bool operator!=(const TimeStepId& a, const TimeStepId& b) {
  return not(a == b);
}

std::ostream& operator<<(std::ostream& s, const UnsizedTimeStepId& id) {
  return s << id.slab_number() << ':' << id.step_time() << ':' << id.substep();
}

std::ostream& operator<<(std::ostream& s, const TimeStepId& id) {
  return s << id.slab_number() << ':' << id.step_time() << ':' << id.substep()
           << ':' << id.substep_time();
}

size_t hash_value(const UnsizedTimeStepId& id) {
  size_t h = 0;
  boost::hash_combine(h, id.slab_number());
  boost::hash_combine(h, id.step_time());
  boost::hash_combine(h, id.substep());
  return h;
}

size_t hash_value(const TimeStepId& id) {
  size_t h = 0;
  boost::hash_combine(h, static_cast<const UnsizedTimeStepId&>(id));
  if (id.substep() != 0) {
    boost::hash_combine(h, id.step_size());
    boost::hash_combine(h, id.substep_time());
  }
  return h;
}

// clang-tidy: do not modify std namespace (okay for hash)
namespace std {  // NOLINT
size_t hash<UnsizedTimeStepId>::operator()(const UnsizedTimeStepId& id) const {
  return boost::hash<UnsizedTimeStepId>{}(id);
}
size_t hash<TimeStepId>::operator()(const TimeStepId& id) const {
  return boost::hash<TimeStepId>{}(id);
}
}  // namespace std

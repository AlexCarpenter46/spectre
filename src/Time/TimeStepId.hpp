// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines class TimeStepId.

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <limits>

#include "Time/Time.hpp"

namespace PUP {
class er;
}  // namespace PUP

// FIXME doc, test, etc.
class UnsizedTimeStepId {
 public:
  UnsizedTimeStepId() = default;
  /// Create an UnsizedTimeStepId at a substep at time `substep_time`
  /// in a step starting at time `step_time`.
  UnsizedTimeStepId(bool time_runs_forward, int64_t slab_number,
                    const Time& step_time, uint64_t substep = 0);

  bool time_runs_forward() const { return time_runs_forward_; }
  int64_t slab_number() const { return slab_number_; }
  /// Time at the start of the current step
  const Time& step_time() const { return step_time_; }
  uint64_t substep() const { return substep_; }

  bool is_at_slab_boundary() const;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

 protected:
  bool time_runs_forward_{false};
  int64_t slab_number_{std::numeric_limits<int64_t>::lowest()};
  Time step_time_{};
  uint64_t substep_{0};
};

/// \ingroup TimeGroup
///
/// A unique identifier for the temporal state of an integrated
/// system.
class TimeStepId : public UnsizedTimeStepId {
 public:
  TimeStepId() = default;
  /// Create a TimeStepId at the start of a step.  If that step is at the
  /// (evolution-defined) end of the slab the TimeStepId will be advanced
  /// to the next slab.
  TimeStepId(bool time_runs_forward, int64_t slab_number, const Time& time);
  /// Create a TimeStepId at a substep at time `substep_time` in a step
  /// starting at time `step_time`.
  TimeStepId(bool time_runs_forward, int64_t slab_number, const Time& step_time,
             uint64_t substep, const TimeDelta& step_size, double substep_time);

  /// Current step size.  Only available when the substep is nonzero,
  /// because on the full step the state is valid for any step size.
  const TimeDelta& step_size() const;
  /// Time of the current substep
  double substep_time() const { return substep_time_; }

  /// Returns a new TimeStepId representing the start of the next step.
  TimeStepId next_step(const TimeDelta& step_size) const;

  /// Returns a new TimeStepId representing the next substep, given
  /// the position of the substep within the step.
  TimeStepId next_substep(const TimeDelta& step_size,
                          double step_fraction) const;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

 private:
  TimeDelta step_size_{};
  double substep_time_{};
};

bool operator==(const UnsizedTimeStepId& a, const UnsizedTimeStepId& b);
bool operator!=(const UnsizedTimeStepId& a, const UnsizedTimeStepId& b);
bool operator<(const UnsizedTimeStepId& a, const UnsizedTimeStepId& b);
bool operator<=(const UnsizedTimeStepId& a, const UnsizedTimeStepId& b);
bool operator>(const UnsizedTimeStepId& a, const UnsizedTimeStepId& b);
bool operator>=(const UnsizedTimeStepId& a, const UnsizedTimeStepId& b);

bool operator==(const TimeStepId& a, const TimeStepId& b);
bool operator!=(const TimeStepId& a, const TimeStepId& b);

std::ostream& operator<<(std::ostream& s, const UnsizedTimeStepId& id);
std::ostream& operator<<(std::ostream& s, const TimeStepId& id);

size_t hash_value(const UnsizedTimeStepId& id);
size_t hash_value(const TimeStepId& id);

namespace std {
template <>
struct hash<UnsizedTimeStepId> {
  size_t operator()(const UnsizedTimeStepId& id) const;
};
template <>
struct hash<TimeStepId> {
  size_t operator()(const TimeStepId& id) const;
};
}  // namespace std

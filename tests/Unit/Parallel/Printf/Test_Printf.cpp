// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <ostream>
#include <vector>

#include "Parallel/Printf/Printf.hpp"

namespace {
struct TestStream {
  double a{1.0};
  std::vector<int> b{0, 4, 8, -7};
};

std::ostream& operator<<(std::ostream& os, const TestStream& t) {
  os << t.a << " (";
  for (size_t i = 0; i < t.b.size() - 1; ++i) {
    os << t.b[i] << ",";
  }
  os << t.b[t.b.size() - 1] << ")";
  return os;
}

enum class TestEnum { Value1, Value2 };

std::ostream& operator<<(std::ostream& os, const TestEnum& t) {
  switch (t) {
    case TestEnum::Value1:
      return os << "Value 1";
    case TestEnum::Value2:
      return os << "Value 2";
    default:
      return os;
  }
}

}  // namespace

// [output_test_example]
// [[OutputRegex, -100 3000000000 1.0000000000000000000e\+00 \(0,4,8,-7\) test 1
// 2 3 abf a o e u Value 2]]
SPECTRE_TEST_CASE("Unit.Parallel.printf", "[Unit][Parallel]") {
  OUTPUT_TEST();
  // clang-tidy doesn't want c-style arrays, but here we are trying
  // to test them explicitly.
  const char c_string0[40] = {"test 1 2 3"}; // NOLINT
  // clang-tidy doesn't want raw pointers, wants gsl::owner<>.
  auto* c_string1 = new char[80]; // NOLINT
  // clang-tidy: do not use pointer arithmetic
  c_string1[0] = 'a';   // NOLINT
  c_string1[1] = 'b';   // NOLINT
  c_string1[2] = 'f';   // NOLINT
  c_string1[3] = '\0';  // NOLINT
  constexpr const char* const c_string2 = {"a o e u"};
  Parallel::printf("%d %lld %s %s %s %s %s\n", -100, 3000000000, TestStream{},
                   c_string0, c_string1, c_string2, TestEnum::Value2);
  // clang-tidy doesn't want delete on anything without gsl::owner<>.
  delete[] c_string1; // NOLINT
}
// [output_test_example]

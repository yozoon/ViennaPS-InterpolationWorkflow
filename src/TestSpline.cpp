#include <set>
#include <vector>

#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include "NaturalCubicSplineInterpolation.hpp"

template <typename KeyContainerType, typename ValueContainerType>
void printCSV(const KeyContainerType &x, const ValueContainerType &y,
              const std::string &filename) {
  auto out = fmt::output_file(filename);
  for (size_t i = 0; i < x.size() && i < y.size(); ++i) {
    auto xi = x.at(i);
    auto vi = y.at(i);
    out.print("{},{}\n", xi, fmt::join(vi, ","));
    fmt::print("{},{}\n", xi, fmt::join(vi, ","));
  }
}

int main() {
  using NumericType = double;

  static constexpr NumericType PI = std::acos(NumericType{-1});

  static constexpr int outputDimension = 2;
  std::vector<NumericType> knots(5);
  std::generate_n(
      knots.begin(), knots.size(),
      [i = 0, n = knots.size()]() mutable { return PI * i++ / (n - 1); });

  std::vector<std::array<NumericType, outputDimension>> f;
  std::transform(knots.begin(), knots.end(), std::back_inserter(f), [](auto v) {
    return std::array<NumericType, outputDimension>{std::sin(v), std::cos(v)};
  });

  printCSV(knots, f, "points.csv");

  std::vector<NumericType> x(50);
  NumericType min = -PI / 4;
  NumericType max = 5. * PI / 4;

  std::generate_n(x.begin(), x.size(), [&, i = 0, n = x.size()]() mutable {
    return min + (max - min) * i++ / (n - 1);
  });

  NaturalCubicSplineInterpolation spline(knots, f);

  std::vector<std::array<NumericType, outputDimension>> interpolated;
  interpolated.reserve(x.size());
  std::transform(x.begin(), x.end(), std::back_inserter(interpolated),
                 [&](auto x) { return spline(x); });

  printCSV(x, interpolated, "spline.csv");
}
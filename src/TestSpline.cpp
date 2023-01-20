#include <set>
#include <vector>

#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include "NaturalCubicSplineInterpolation.hpp"
#include "NaturalCubicSplineInterpolationEIGEN.hpp"

template <typename KeyContainerType, typename ValueContainerType>
void printCSV(const KeyContainerType &x, const ValueContainerType &y,
              const std::string &filename) {
  auto out = fmt::output_file(filename);
  for (size_t i = 0; i < x.size() && i < y.size(); ++i) {
    auto xi = x.at(i);
    auto vi = y.at(i);
    out.print("{},{}\n", xi, fmt::join(vi, ","));
    // fmt::print("{},{}\n", xi, fmt::join(vi, ","));
  }
}

template <typename KeyContainerType, typename ValueContainerType>
void printCSV(const KeyContainerType &x, const ValueContainerType &r,
              const ValueContainerType &y, const std::string &filename) {
  auto out = fmt::output_file(filename);
  for (size_t i = 0; i < x.size() && i < r.size() && i < y.size(); ++i) {
    auto xi = x.at(i);
    auto ri = r.at(i);
    auto vi = y.at(i);
    out.print("{},{},{}\n", xi, fmt::join(ri, ","), fmt::join(vi, ","));
    // fmt::print("{},{},{}\n", xi, fmt::join(ri, ","), fmt::join(vi, ","));
  }
}

int main() {
  using NumericType = double;

  static constexpr NumericType PI = std::acos(NumericType{-1});

  std::vector<NumericType> knots(5);
  std::generate_n(
      knots.begin(), knots.size(),
      [i = 0, n = knots.size()]() mutable { return PI * i++ / (n - 1); });

  std::vector<NumericType (*)(NumericType)> functions;
  functions.push_back([](NumericType x) { return std::sin(x); });
  functions.push_back([](NumericType x) { return std::cos(x); });
  functions.push_back([](NumericType x) {
    return std::exp(-2. * (x - PI / 2) * (x - PI / 2));
  });

  unsigned outputDimension = functions.size();

  std::vector<std::vector<NumericType>> f(
      knots.size(), std::vector<NumericType>(outputDimension));
  for (unsigned i = 0; i < knots.size(); ++i) {
    for (unsigned j = 0; j < outputDimension; ++j)
      f[i][j] = functions[j](knots[i]);
  }

  printCSV(knots, f, "points.csv");

  unsigned numberOfSamples = 50;
  std::vector<NumericType> x(numberOfSamples);
  NumericType min = -PI / 4;
  NumericType max = 5. * PI / 4;

  std::generate_n(x.begin(), x.size(), [&, i = 0, n = x.size()]() mutable {
    return min + (max - min) * i++ / (n - 1);
  });

  std::vector<std::vector<NumericType>> reference;
  reference.reserve(numberOfSamples);
  std::transform(x.begin(), x.end(), std::back_inserter(reference),
                 [&functions](NumericType x) {
                   std::vector<NumericType> r;
                   for (unsigned i = 0; i < functions.size(); ++i)
                     r.push_back(functions[i](x));
                   return r;
                 });

  {
    NaturalCubicSplineInterpolation<NumericType> spline(knots, f);

    std::vector<std::vector<NumericType>> interpolated;
    interpolated.reserve(x.size());
    std::transform(x.begin(), x.end(), std::back_inserter(interpolated),
                   [&](auto x) { return spline(x); });

    printCSV(x, reference, interpolated, "spline_lapack.csv");
  }

  {
    NaturalCubicSplineInterpolationEIGEN<NumericType> spline(knots, f);

    std::vector<std::vector<NumericType>> interpolated;
    interpolated.reserve(x.size());
    std::transform(x.begin(), x.end(), std::back_inserter(interpolated),
                   [&](auto x) { return spline(x); });

    printCSV(x, reference, interpolated, "spline_eigen.csv");
  }
}
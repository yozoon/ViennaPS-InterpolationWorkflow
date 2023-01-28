#include <functional>
#include <set>
#include <vector>

#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <psSmartPointer.hpp>

#include "CubicSplineInterpolation.hpp"

#include "SplineGridInterpolation.hpp"

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

  std::vector<NumericType> knots(7);

  NumericType min = -1.;
  NumericType max = 1.;
  std::generate_n(knots.begin(), knots.size(),
                  [&, i = 0, n = knots.size()]() mutable {
                    return min + (max - min) * i++ / (n - 1);
                  });

  std::vector<NumericType (*)(NumericType)> functions;
  functions.push_back([](NumericType x) { return 1 + x * x; });
  functions.push_back([](NumericType x) { return std::sin(x); });
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
  auto [minIt, maxIt] = std::minmax_element(knots.begin(), knots.end());
  auto range = *maxIt - *minIt;
  NumericType overshoot = 0.5;

  std::generate_n(
      x.begin(), x.size(),
      [min = min - overshoot * range, max = max + overshoot * range, i = 0,
       n = x.size()]() mutable { return min + (max - min) * i++ / (n - 1); });

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
    CubicSplineInterpolation<NumericType> spline(
        knots, f, SplineBoundaryConditionType::NOT_A_KNOT);

    std::vector<std::vector<NumericType>> interpolated;
    interpolated.reserve(x.size());
    std::transform(x.begin(), x.end(), std::back_inserter(interpolated),
                   [&](auto x) { return spline(x); });

    printCSV(x, reference, interpolated, "spline_lapack.csv");
  }

  {
    // The function(s) to interpolate
    std::vector<std::function<NumericType(NumericType, NumericType)>> functions;
    functions.push_back([](NumericType x, NumericType) { return 1. + x * x; });
    functions.push_back([&](NumericType x, NumericType y) {
      NumericType xp = x - PI / 2;
      NumericType yp = y - PI / 2;
      return std::exp(-(xp * xp + yp * yp) / 2);
    });

    functions.push_back([](NumericType x, NumericType y) {
      return std::sin(x) * std::cos(3. * y);
    });

    // Create sample points and save them to a file
    int nx = 4, ny = 7;
    NumericType minX = 0.;
    NumericType maxX = PI;
    NumericType minY = 0.;
    NumericType maxY = PI;

    auto out = fmt::output_file("points_nd.csv");
    std::vector<std::vector<NumericType>> data;
    for (int i = 0; i < nx; ++i) {
      NumericType x = minX + (maxX - minX) * i / (nx - 1);
      for (int j = 0; j < ny; ++j) {
        NumericType y = minY + (maxY - minY) * j / (ny - 1);
        std::vector<NumericType> tmp;
        tmp.push_back(x);
        tmp.push_back(y);
        out.print("{},{}\n", x, y);
        for (unsigned i = 0; i < functions.size(); ++i)
          tmp.push_back(functions[i](x, y));
        data.push_back(tmp);
      }
    }

    // Instantiate the interpolation class ans set the data
    SplineGridInterpolation<NumericType> sgi;
    sgi.setDataDimensions(2, functions.size());
    sgi.setBCType(SplineBoundaryConditionType::NOT_A_KNOT);
    sgi.setData(
        psSmartPointer<const std::vector<std::vector<NumericType>>>::New(data));

    // Interpolate along a grid
    int resolution = 30;

    NumericType overshoot = 0.2;
    NumericType minXi = minX - overshoot * (maxX - minX);
    NumericType maxXi = maxX + overshoot * (maxX - minX);
    NumericType minYi = minY - overshoot * (maxY - minY);
    NumericType maxYi = maxY + overshoot * (maxY - minY);

    auto out2 = fmt::output_file("spline_nd.csv");
    for (int i = 0; i < resolution; ++i) {
      NumericType x = minXi + (maxXi - minXi) * i / (resolution - 1);
      for (int j = 0; j < resolution; ++j) {
        NumericType y = minYi + (maxYi - minYi) * j / (resolution - 1);
        out2.print("{}, {}", x, y);
        for (auto &f : functions) {
          out2.print(",{}", f(x, y));
        }
        auto [pred, _] = sgi.estimate(std::vector{x, y}).value();
        out2.print(",{}\n", fmt::join(pred, ","));
      }
    }
  }
}
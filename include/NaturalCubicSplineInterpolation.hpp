#ifndef NATURAL_CUBIC_SPLINE_INTERPOLATION_HPP
#define NATURAL_CUBIC_SPLINE_INTERPOLATION_HPP

#include <algorithm>
#include <cassert>
#include <exception>
#include <iterator>
#include <memory>
#include <set>
#include <vector>

extern "C" {
void dgtsv_(int *, int *, double *, double *, double *, double *, int *, int *);
void sgtsv_(int *, int *, float *, float *, float *, float *, int *, int *);
}

template <typename NumericType> class NaturalCubicSplineInterpolation {
  static_assert(std::is_same_v<NumericType, float> ||
                    std::is_same_v<NumericType, double>,
                "NumericType is neither double nor float. No "
                "matching LAPACK function found.");

public:
  NaturalCubicSplineInterpolation(
      const std::vector<NumericType> &passedX,
      const std::vector<std::vector<NumericType>> &passedY)
      : N(passedX.size()) {

    if (passedX.size() != passedY.size())
      throw std::invalid_argument("NaturalCubicSplineInterpolation: The "
                                  "number of elements in the X "
                                  "and Y vector do not match.");

    // Make sure that the input values are sorted and not duplicate.
    // This is achieved by creating a permutation index set, whose comparison
    // function is based on the comparison of the passed X values. Thus, upon
    // inserting incrementing keys, the key will be placed in the set at the
    // position where the X value with that index would be located. If two or
    // more of the provided X values are the same, only the first value will be
    // inserted.
    const auto comparator = [&passedX](int i, int j) {
      return passedX[i] < passedX[j];
    };
    std::set<int, decltype(comparator)> indices(comparator);
    for (int i = 0; i < N; ++i)
      indices.insert(i);

    if (indices.size() < 3)
      throw std::invalid_argument(
          "NaturalCubicSplineInterpolation: Not enough unique X "
          "values were provided to apply cubic spline interpolation.");

    // The first element determins the output dimension
    outputDimension = passedY[0].size();

    knots.reserve(indices.size());
    y.resize(indices.size(), std::vector<NumericType>(outputDimension));

    // Use the permutation index array to actually populate our local vectors
    // holding knots and values.
    for (auto index : indices) {
      knots.push_back(passedX[index]);
      y[index] = passedY[index];
    }
  }

  std::vector<NumericType> operator()(NumericType x) {
    if (!initialized)
      initialize();

    auto lowerBound = std::lower_bound(knots.begin(), knots.end(), x);

    int i = 1;
    int d = std::distance(knots.begin(), lowerBound);
    i = std::clamp(d, 1, N - 1);

    NumericType t = (x - knots[i - 1]) / (knots[i] - knots[i - 1]);

    std::vector<NumericType> result(outputDimension, 0.);
    for (int j = 0; j < outputDimension; ++j) {
      result[j] = (1. - t) * y[i - 1][j] + t * y[i][j] +
                  t * (1. - t) *
                      ((1. - t) * a[i * outputDimension + j] +
                       t * b[i * outputDimension + j]);
    }
    return result;
  }

private:
  void initialize() {
    // Lower diagonal elements
    auto lowerDiag = std::make_unique<NumericType[]>(N - 1);
    for (int i = 1; i < N - 1; ++i)
      lowerDiag[i - 1] = 1.0 / (knots[i] - knots[i - 1]);
    lowerDiag[N - 2] = 1.0 / (knots[N - 1] - knots[N - 2]);

    // Diagonal elements
    auto diag = std::make_unique<NumericType[]>(N);
    diag[0] = 2.0 / (knots[1] - knots[0]);
    for (int i = 1; i < N - 1; ++i)
      diag[i] = 2.0 * (1.0 / (knots[i] - knots[i - 1]) +
                       1.0 / (knots[i + 1] - knots[i]));
    diag[N - 1] = 2.0 / (knots[N - 1] - knots[N - 2]);

    // Upper diagonal elements
    auto upperDiag = std::make_unique<NumericType[]>(N - 1);
    upperDiag[0] = 1.0 / (knots[1] - knots[0]);
    for (int i = 1; i < N - 1; ++i)
      upperDiag[i] = 1.0 / (knots[i + 1] - knots[i]);

    // RHS matrix (Column major layout!)
    auto B = std::make_unique<NumericType[]>(outputDimension * N);

    for (int j = 0; j < outputDimension; ++j) {
      // First row
      B[j * N] = 3.0 * (y[1][j] - y[0][j]) / (knots[1] - knots[0]);
      // Last row
      B[j * N + N - 1] =
          3.0 * (y[N - 1][j] - y[N - 2][j]) / (knots[N - 1] - knots[N - 2]);
    }

    for (int i = 1; i < N - 1; ++i) {
      for (int j = 0; j < outputDimension; ++j) {
        B[j * N + i] =
            3.0 * ((y[i][j] - y[i - 1][j]) / (knots[i] - knots[i - 1]) +
                   (y[i + 1][j] - y[i][j]) / (knots[i + 1] - knots[i]));
      }
    }

    int info;
    // Solve the system(s)
    if constexpr (std::is_same_v<NumericType, float>) {
      sgtsv_(&N /* N */, &outputDimension /* NRHS */, lowerDiag.get() /* DL */,
             diag.get() /* D */, upperDiag.get() /* DU */, B.get(),
             &N /* LDB */, &info /* INFO */);
    } else if constexpr (std::is_same_v<NumericType, double>) {
      dgtsv_(&N /* N */, &outputDimension /* outputDimension */,
             lowerDiag.get() /* DL */, diag.get() /* D */,
             upperDiag.get() /* DU */, B.get(), &N /* LDB */, &info /* INFO */);
    }

    a = std::make_unique<NumericType[]>(N * outputDimension);
    b = std::make_unique<NumericType[]>(N * outputDimension);

    for (int i = 1; i < N; ++i) {
      for (int j = 0; j < outputDimension; ++j) {
        a[i * outputDimension + j] =
            B[j * N + i - 1] * (knots[i] - knots[i - 1]) -
            (y[i][j] - y[i - 1][j]);
        b[i * outputDimension + j] =
            -B[j * N + i] * (knots[i] - knots[i - 1]) + (y[i][j] - y[i - 1][j]);
      }
    }

    initialized = true;
  }

  int N;
  int outputDimension = 0;
  bool initialized = false;

  std::vector<NumericType> knots;
  std::vector<std::vector<NumericType>> y;

  std::unique_ptr<NumericType[]> a = nullptr;
  std::unique_ptr<NumericType[]> b = nullptr;
};

#endif
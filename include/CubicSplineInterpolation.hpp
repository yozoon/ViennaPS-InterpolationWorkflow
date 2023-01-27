#ifndef CUBIC_SPLINE_INTERPOLATION_HPP
#define CUBIC_SPLINE_INTERPOLATION_HPP

#include <algorithm>
#include <cassert>
#include <exception>
#include <iterator>
#include <memory>
#include <set>
#include <vector>

extern "C" {
// Tridiagonal matrix solver routines
void dgtsv_(int *, int *, double *, double *, double *, double *, int *, int *);
void sgtsv_(int *, int *, float *, float *, float *, float *, int *, int *);
// General banded matrix solver routines
void dgbsv_(int *, int *, int *, int *, double *, int *, int *, double *, int *,
            int *);
void sgbsv_(int *, int *, int *, int *, float *, int *, int *, float *, int *,
            int *);
}

// #include <fmt/core.h>
// #include <fmt/ranges.h>

// template <typename NumericType>
// void printMatrix(NumericType *mat, int numRows, int numCols) {
//   for (int i = 0; i < numRows; ++i) {
//     for (int j = 0; j < numCols; ++j) {
//       fmt::print("{:.3f}, ", mat[j * numRows + i]);
//     }
//     fmt::print("\n");
//   }
// }

// template <typename NumericType>
// void printTridiagMatrix(NumericType *upperDiag, NumericType *diag,
//                         NumericType *lowerDiag, int size) {
//   for (int i = 0; i < size; ++i) {
//     for (int j = 0; j < size; ++j) {
//       if (i == j) {
//         // diag
//         fmt::print("{:.3f}, ", diag[i]);
//       } else if (j == i - 1) {
//         // Lower diag
//         fmt::print("{:.3f}, ", lowerDiag[j]);
//       } else if (j == i + 1) {
//         // upper diag
//         fmt::print("{:.3f}, ", upperDiag[i]);
//       } else {
//         fmt::print("0.000, ");
//       }
//     }
//     fmt::print("\n");
//   }
// }

enum SplineBoundaryConditionType : unsigned {
  NOT_A_KNOT,
  NATURAL,
};

template <typename NumericType> class CubicSplineInterpolation {
  static_assert(std::is_same_v<NumericType, float> ||
                    std::is_same_v<NumericType, double>,
                "NumericType is neither double nor float. No "
                "matching LAPACK function found.");

  const SplineBoundaryConditionType bcType;

public:
  CubicSplineInterpolation(const std::vector<NumericType> &passedX,
                           const std::vector<std::vector<NumericType>> &passedY,
                           SplineBoundaryConditionType passedBCType =
                               SplineBoundaryConditionType::NOT_A_KNOT)
      : bcType(passedBCType), N(passedX.size()) {

    if (passedX.size() != passedY.size())
      throw std::invalid_argument("CubicSplineInterpolation: The "
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

    if (indices.size() < 4)
      throw std::invalid_argument(
          "CubicSplineInterpolation: Not enough unique X "
          "values were provided to apply cubic spline interpolation. Cubic "
          "spline interpolation requires at least four points.");

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
    // fmt::print("knots={}\n", knots);
    // fmt::print("{}\n", y);
  }

  std::vector<NumericType> operator()(NumericType x) {
    if (!initialized)
      initialize();

    auto lowerBound = std::lower_bound(knots.begin(), knots.end(), x);

    int d = std::distance(knots.begin(), lowerBound);
    int i = std::clamp(d, 1, N - 1);

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

    // RHS matrix (Column major layout!)
    auto B = std::make_unique<NumericType[]>(outputDimension * N);

    auto toColMajIndex = [&](int nrows, int ncols, int row, int col) {
      assert(row < nrows && "Matrix row index too large!");
      assert(col < ncols && "Matrix column index too large!");
      return col * nrows + row;
    };

    for (int i = 1; i < N - 1; ++i) {
      for (int j = 0; j < outputDimension; ++j) {
        B[toColMajIndex(N, outputDimension, i, j)] =
            3.0 * ((y[i][j] - y[i - 1][j]) /
                       ((knots[i] - knots[i - 1]) * (knots[i] - knots[i - 1])) +
                   (y[i + 1][j] - y[i][j]) /
                       ((knots[i + 1] - knots[i]) * (knots[i + 1] - knots[i])));
      }
    }

    if (bcType == SplineBoundaryConditionType::NATURAL) {
      // RHS entries for natural BC
      for (int j = 0; j < outputDimension; ++j) {
        // First row
        NumericType dx1 = knots[1] - knots[0];
        B[toColMajIndex(N, outputDimension, 0, j)] =
            3.0 * (y[1][j] - y[0][j]) / (dx1 * dx1);
        // Last row
        NumericType dxn1 = knots[N - 1] - knots[N - 2];
        B[toColMajIndex(N, outputDimension, N - 1, j)] =
            3.0 * (y[N - 1][j] - y[N - 2][j]) / (dxn1 * dxn1);
      }

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

      // fmt::print("A=\n");
      // printTridiagMatrix(upperDiag.get(), diag.get(), lowerDiag.get(), N);
      // fmt::print("B=\n");
      // printMatrix(B.get(), N, outputDimension);

      int info;
      // Solve the system(s)
      if constexpr (std::is_same_v<NumericType, float>) {
        sgtsv_(&N /* N */, &outputDimension /* NRHS */,
               lowerDiag.get() /* DL */, diag.get() /* D */,
               upperDiag.get() /* DU */, B.get(), &N /* LDB */,
               &info /* INFO */);
      } else if constexpr (std::is_same_v<NumericType, double>) {
        dgtsv_(&N /* N */, &outputDimension /* outputDimension */,
               lowerDiag.get() /* DL */, diag.get() /* D */,
               upperDiag.get() /* DU */, B.get(), &N /* LDB */,
               &info /* INFO */);
      }
    } else if (bcType == SplineBoundaryConditionType::NOT_A_KNOT) {
      // RHS entries for natural BC
      for (int j = 0; j < outputDimension; ++j) {
        // First row
        NumericType dx1 = knots[1] - knots[0];
        NumericType dx2 = knots[2] - knots[1];
        B[j * N] = 2.0 * ((y[1][j] - y[0][j]) / (dx1 * dx1 * dx1) -
                          (y[2][j] - y[1][j]) / (dx2 * dx2 * dx2));
        // Last row
        NumericType dxn2 = knots[N - 2] - knots[N - 3];
        NumericType dxn1 = knots[N - 1] - knots[N - 2];
        B[j * N + N - 1] =
            2.0 * ((y[N - 2][j] - y[N - 3][j]) / (dxn2 * dxn2 * dxn2) -
                   (y[N - 1][j] - y[N - 2][j]) / (dxn1 * dxn1 * dxn1));
      }

      // The system matrix
      int KL = 2;
      int KU = 2;
      int LDAB = std::max(2 * KL + KU + 1, N);
      auto A = std::make_unique<NumericType[]>(LDAB * N);

      auto toBandedIndex = [&](int nrows, int ncols, int row, int col) {
        // return toColMajIndex(nrows, ncols, row, col);
        return toColMajIndex(nrows, ncols, KL + KU + row - col, col);
      };

      for (int i = 1; i < N - 1; ++i) {
        // Diagonal
        A[toBandedIndex(LDAB, N, i, i)] =
            2.0 *
            (1.0 / (knots[i] - knots[i - 1]) + 1.0 / (knots[i + 1] - knots[i]));

        // Diagonal with offet -1 (lower diagonal)
        A[toBandedIndex(LDAB, N, i, i - 1)] = 1.0 / (knots[i] - knots[i - 1]);

        // Diagonal with offet -2
        if (i > 1)
          A[toBandedIndex(LDAB, N, i, i - 2)] = 0.;

        // Diagonal with offet +1 (upper diagonal)
        A[toBandedIndex(LDAB, N, i, i + 1)] = 1.0 / (knots[i + 1] - knots[i]);

        // Diagonal with offet +2
        if (i < N - 2)
          A[toBandedIndex(LDAB, N, i, i + 2)] = 0.;
      }

      NumericType dx1 = knots[1] - knots[0];
      NumericType dx2 = knots[2] - knots[1];
      A[toBandedIndex(LDAB, N, 0, 0)] = 1.0 / (dx1 * dx1);
      A[toBandedIndex(LDAB, N, 0, 1)] = 1.0 / (dx1 * dx1) - 1.0 / (dx2 * dx2);
      A[toBandedIndex(LDAB, N, 0, 2)] = -1.0 / (dx2 * dx2);

      NumericType dxn1 = knots[N - 1] - knots[N - 2];
      NumericType dxn2 = knots[N - 2] - knots[N - 3];
      A[toBandedIndex(LDAB, N, N - 1, N - 3)] = 1.0 / (dxn2 * dxn2);
      A[toBandedIndex(LDAB, N, N - 1, N - 2)] =
          1.0 / (dxn2 * dxn2) - 1.0 / (dxn1 * dxn1);
      A[toBandedIndex(LDAB, N, N - 1, N - 1)] = -1.0 / (dxn1 * dxn1);

      // fmt::print("A=\n");
      // printMatrix(A.get(), LDAB, N);

      auto ipiv = std::make_unique<int[]>(N);
      int info;
      if constexpr (std::is_same_v<NumericType, float>) {
        sgbsv_(&N /* N */, &KL /* KL */, &KU /* KU*/,
               &outputDimension /* output dimension*/, A.get() /* AB */,
               &LDAB /* LDAB */, ipiv.get() /* IPIV */, B.get() /* B */,
               &N /* LDB */, &info /* INFO */);
      } else if constexpr (std::is_same_v<NumericType, double>) {
        dgbsv_(&N /* N */, &KL /* KL */, &KU /* KU*/,
               &outputDimension /* output dimension*/, A.get() /* AB */,
               &LDAB /* LDAB */, ipiv.get() /* IPIV */, B.get() /* B */,
               &N /* LDB */, &info /* INFO */);
      }
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
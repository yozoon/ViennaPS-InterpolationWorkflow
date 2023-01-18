#ifndef NATURAL_CUBIC_SPLINE_INTERPOLATION_HPP
#define NATURAL_CUBIC_SPLINE_INTERPOLATION_HPP

#include <algorithm>
#include <cassert>
#include <exception>
#include <iterator>
#include <set>
#include <vector>

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

template <typename KnotContainerType, typename ValueContainerType>
class NaturalCubicSplineInterpolation {
  using NumericType = typename KnotContainerType::value_type;
  using OutputType = typename ValueContainerType::value_type;
  using SizeType = size_t;

  static constexpr SizeType outputDimension = std::tuple_size_v<OutputType>;
  bool initialized = false;

public:
  NaturalCubicSplineInterpolation(const KnotContainerType &passedX,
                                  const ValueContainerType &passedY)
      : N(passedX.size()) {

    if (passedX.size() != passedY.size())
      throw std::invalid_argument(
          "NaturalCubicSplineInterpolation: The number of elements in the X "
          "and Y vector do not match.");

    // Make sure that the input values are sorted and not duplicate.
    // This is achieved by creating a permutation index set, whose comparison
    // function is based on the comparison of the passed X values. Thus, upon
    // inserting incrementing keys, the key will be placed in the set at the
    // position where the X value with that index would be located. If two or
    // more of the provided X values are the same, only the first value will be
    // inserted.
    const auto comparator = [&passedX](SizeType i, SizeType j) {
      return passedX[i] < passedX[j];
    };
    std::set<SizeType, decltype(comparator)> indices(comparator);
    for (SizeType i = 0; i < N; ++i)
      indices.insert(i);

    if (indices.size() < 3)
      throw std::invalid_argument(
          "NaturalCubicSplineInterpolation: Not enough unique X "
          "values were provided to apply cubic spline interpolation.");

    knots.reserve(indices.size());
    f.reserve(indices.size());

    // Use the permutation index array to actually populate our local vectors
    // holding knots and values.
    for (auto index : indices) {
      knots.push_back(passedX[index]);
      f.push_back(passedY[index]);
    }
  }

  OutputType operator()(NumericType x) {
    if (!initialized)
      initialize();

    auto lowerBound = std::lower_bound(knots.begin(), knots.end(), x);

    SizeType i = 1;
    SizeType d = std::distance(knots.begin(), lowerBound);
    i = std::clamp(d, SizeType{1}, SizeType{N - 1});

    NumericType t = (x - knots.at(i - 1)) / (knots.at(i) - knots.at(i - 1));

    OutputType result;
    for (SizeType j = 0; j < outputDimension; ++j) {
      result[j] = (1. - t) * f.at(i - 1).at(j) + t * f.at(i).at(j) +
                  t * (1. - t) * ((1. - t) * a(i, j) + t * b(i, j));
    }
    return result;
  }

private:
  void initialize() {
    using TripletType = Eigen::Triplet<NumericType>;

    std::vector<TripletType> triplets;
    triplets.reserve(3 * N - 2);

    // Initialize the system matrix
    for (SizeType i = 1; i < N - 1; ++i) {
      for (SizeType j = std::min(i - 1, SizeType{0}); j < std::max(i + 1, N);
           ++j) {
        if (j == i) {
          // Diagonal element
          triplets.push_back(
              TripletType(i, j,
                          2.0 * (1.0 / (knots.at(i) - knots.at(i - 1)) +
                                 1.0 / (knots.at(i + 1) - knots.at(i)))));
        } else if (j - 1 == i) {
          // Left to diagonal
          triplets.push_back(
              TripletType(i, j, 1.0 / (knots.at(i) - knots.at(i - 1))));
        } else if (j + 1 == i) {
          // Right to diagonal
          triplets.push_back(
              TripletType(i, j, 1.0 / (knots.at(i + 1) - knots.at(i))));
        }
      }
    }

    // Now also set the values required for the "natural spline" conditions
    triplets.push_back(TripletType(0, 0, 2.0 / (knots.at(1) - knots.at(0))));
    triplets.push_back(TripletType(0, 1, 1.0 / (knots.at(1) - knots.at(0))));
    triplets.push_back(
        TripletType(N - 1, N - 1, 2.0 / (knots.at(N - 1) - knots.at(N - 2))));
    triplets.push_back(
        TripletType(N - 1, N - 2, 1.0 / (knots.at(N - 1) - knots.at(N - 2))));

    Eigen::SparseMatrix<NumericType> A(N, N);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();

    // Create the RHS matrix
    Eigen::Matrix<NumericType, Eigen::Dynamic, Eigen::Dynamic> rhs(
        N, outputDimension);

    for (SizeType i = 1; i < N - 1; ++i)
      for (SizeType j = 0; j < outputDimension; ++j) {
        rhs(i, j) = 3.0 * ((f.at(i).at(j) - f.at(i - 1).at(j)) /
                               (knots.at(i) - knots.at(i - 1)) +
                           (f.at(i + 1).at(j) - f.at(i).at(j)) /
                               (knots.at(i + 1) - knots.at(i)));
      }

    // Now also set the values required for the "natural spline" conditions
    for (SizeType j = 0; j < outputDimension; ++j) {
      rhs(0, j) =
          3.0 * (f.at(1).at(j) - f.at(0).at(j)) / (knots.at(1) - knots.at(0));
      rhs(N - 1, j) = 3.0 * (f.at(N - 1).at(j) - f.at(N - 2).at(j)) /
                      (knots.at(N - 1) - knots.at(N - 2));
    }

    // Solve the system using SparseLU method
    Eigen::SparseLU<Eigen::SparseMatrix<NumericType>,
                    Eigen::COLAMDOrdering<int>>
        solver;
    solver.compute(A);
    Eigen::Matrix<NumericType, Eigen::Dynamic, Eigen::Dynamic> k =
        solver.solve(rhs);

    a.resize(N, outputDimension);
    b.resize(N, outputDimension);

    for (SizeType i = 1; i < N; ++i) {
      for (SizeType j = 0; j < outputDimension; ++j) {
        a(i, j) = k(i - 1, j) * (knots.at(i) - knots.at(i - 1)) -
                  (f.at(i).at(j) - f.at(i - 1).at(j));
      }
      for (SizeType j = 0; j < outputDimension; ++j) {
        b(i, j) = -k(i, j) * (knots.at(i) - knots.at(i - 1)) +
                  (f.at(i).at(j) - f.at(i - 1).at(j));
      }
    }

    initialized = true;
  }
  const SizeType N;
  KnotContainerType knots;
  ValueContainerType f;

  Eigen::Matrix<NumericType, Eigen::Dynamic, Eigen::Dynamic> a;
  Eigen::Matrix<NumericType, Eigen::Dynamic, Eigen::Dynamic> b;
};

#endif
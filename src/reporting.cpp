#define R_NO_REMAP
#include <cpp11.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace cpp11;

[[cpp11::register]]
writable::list mfrm_cpp_row_predictions(
    doubles eta, doubles_matrix<> steps, integers group, integers lower,
    integers upper, double rating_min, doubles alpha, doubles theta,
    doubles offset, bool ability_only) {
  const int n = eta.size();
  const int categories = steps.ncol();
  writable::doubles_matrix<> probs(n, categories);
  writable::doubles expected(n), variance(n);
  std::vector<double> mass(categories);

  for (int i = 0; i < n; ++i) {
    expected[i] = variance[i] = NA_REAL;
    const double discrimination = alpha.size() ? alpha[i] : 1.0;
    const int g = group[i];
    if (std::isnan(eta[i]) || !std::isfinite(discrimination) || discrimination <= 0 ||
        g < 0 || g >= steps.nrow()) {
      for (int k = 0; k < categories; ++k) probs(i, k) = NA_REAL;
      continue;
    }
    const int first = lower[g], last = upper[g];
    if (first < 0 || last >= categories || first > last) stop("Invalid category support.");
    for (int k = 0; k < categories; ++k) probs(i, k) = 0;
    if (std::isinf(eta[i])) {
      const int extreme = eta[i] > 0 ? last : first;
      probs(i, extreme) = 1;
      expected[i] = rating_min + extreme;
      variance[i] = 0;
      continue;
    }

    const bool shift_ability = ability_only && discrimination != 1 && theta.size() && offset.size();
    const double location = shift_ability ? discrimination * theta[i] + offset[i] : eta[i];
    const double multiplier = shift_ability ? 1.0 : discrimination;
    const double minus_infinity = -std::numeric_limits<double>::infinity();
    double maximum = minus_infinity;
    for (int k = first; k <= last; ++k) {
      const double threshold = steps(g, k);
      mass[k] = std::isfinite(threshold) ? multiplier * ((k - first) * location - threshold) : minus_infinity;
      maximum = std::max(maximum, mass[k]);
    }
    if (maximum == minus_infinity) continue;
    double total = 0;
    for (int k = first; k <= last; ++k) {
      mass[k] = std::exp(mass[k] - maximum);
      total += mass[k];
    }
    double mean = 0;
    for (int k = first; k <= last; ++k) {
      const double p = mass[k] / total;
      probs(i, k) = p;
      mean += k * p;
    }
    expected[i] = rating_min + mean;
    double squared_deviations = 0;
    for (int k = first; k <= last; ++k) {
      const double difference = k - mean;
      squared_deviations += difference * difference * probs(i, k);
    }
    variance[i] = squared_deviations;
  }

  return writable::list({"probs"_nm = probs, "expected"_nm = expected, "variance"_nm = variance});
}

[[cpp11::register]]
writable::doubles mfrm_cpp_expected_agreement(
    list groups, doubles scores, doubles_matrix<> probs) {
  double units = 0, pairs = 0, exact = 0, adjacent = 0;
  double expected_exact = 0, expected_adjacent = 0;
  for (SEXP group : groups) {
    integers rows(group);
    bool counted = false;
    for (R_xlen_t i = 0; i < rows.size(); ++i) {
      const int a = rows[i] - 1;
      for (R_xlen_t j = i + 1; j < rows.size(); ++j) {
        const int b = rows[j] - 1;
        double exact_pair = 0, adjacent_pair = 0;
        double total_a = 0, total_b = 0;
        for (int k = 0; k < probs.ncol(); ++k) {
          const double pa = probs(a, k), pb = probs(b, k);
          total_a += pa;
          total_b += pb;
          exact_pair += pa * pb;
          adjacent_pair += pa * (pb + (k > 0 ? probs(b, k - 1) : 0) +
                                      (k + 1 < probs.ncol() ? probs(b, k + 1) : 0));
        }
        if (!std::isfinite(total_a) || !std::isfinite(total_b) ||
            std::abs(total_a - 1) > 1e-8 || std::abs(total_b - 1) > 1e-8) continue;
        counted = true;
        ++pairs;
        exact += scores[a] == scores[b];
        adjacent += std::abs(scores[a] - scores[b]) <= 1.0;
        expected_exact += exact_pair;
        expected_adjacent += adjacent_pair;
      }
    }
    if (counted) ++units;
  }
  return writable::doubles({units, pairs, exact, adjacent, expected_exact, expected_adjacent});
}

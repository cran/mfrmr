#define R_NO_REMAP
#include <cpp11.hpp>
#include <Rmath.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

using namespace cpp11;

namespace {

// Share only bit-identical finite inputs on a common quadrature grid. The
// mapping is local to this call, so no parameter-dependent cache survives it.
struct ProbabilityRows {
  std::vector<int> group;
  std::vector<int> first;
};

ProbabilityRows probability_representatives(doubles base_eta,
                                             doubles_matrix<> person_nodes,
                                             const double* steps,
                                             int k_cat,
                                             bool observation_steps,
                                             const double* slopes = nullptr) {
  const int n = base_eta.size();
  if (n < 32 || person_nodes.ncol() < 2 || k_cat < 2 ||
      person_nodes.nrow() == 0) {
    return {};
  }
  for (int q = 0; q < person_nodes.ncol(); ++q) {
    const double first = person_nodes(0, q);
    if (!std::isfinite(first)) return {};
    for (int p = 1; p < person_nodes.nrow(); ++p) {
      const double node = person_nodes(p, q);
      if (std::memcmp(&first, &node, sizeof(double)) != 0) return {};
    }
  }
  if (!observation_steps) {
    for (int k = 0; k < k_cat; ++k) {
      if (!std::isfinite(steps[k])) return {};
    }
  }

  std::unordered_map<std::string, int> first_rows;
  ProbabilityRows rows;
  rows.group.resize(static_cast<size_t>(n));
  std::string key;
  const size_t key_doubles = 1 + (slopes != nullptr ? 1 : 0) +
                             (observation_steps ? static_cast<size_t>(k_cat) : 0);
  key.reserve(sizeof(double) * key_doubles);
  for (int i = 0; i < n; ++i) {
    const double eta = base_eta[i];
    if (!std::isfinite(eta)) return {};
    key.assign(reinterpret_cast<const char*>(&eta), sizeof(double));
    if (slopes != nullptr) {
      const double sl = slopes[i];
      if (!std::isfinite(sl)) return {};
      key.append(reinterpret_cast<const char*>(&sl), sizeof(double));
    }
    if (observation_steps) {
      for (int k = 0; k < k_cat; ++k) {
        const double step = steps[static_cast<size_t>(i) + static_cast<size_t>(n) * k];
        if (!std::isfinite(step)) return {};
        key.append(reinterpret_cast<const char*>(&step), sizeof(double));
      }
    }
    const auto inserted = first_rows.emplace(key, static_cast<int>(rows.first.size()));
    if (inserted.second) rows.first.push_back(i);
    rows.group[static_cast<size_t>(i)] = inserted.first->second;
    // With few repetitions, retain the original loop and its smaller workspace.
    if (first_rows.size() > static_cast<size_t>(n / 8)) return {};
  }
  // Bound the extra category cache by the already allocated log-probability matrix.
  if (rows.first.size() * static_cast<size_t>(k_cat) >
      static_cast<size_t>(n) * person_nodes.ncol()) return {};
  return rows;
}

ProbabilityRows probability_representatives_indexed(int model_code,
                                                    doubles base_eta,
                                                    doubles_matrix<> person_nodes,
                                                    const double* step_rsm,
                                                    const std::vector<double>& step_rm,
                                                    int n_criteria,
                                                    const int* step_idx_ptr,
                                                    const double* slopes_ptr,
                                                    int n_slopes,
                                                    const int* slope_idx_ptr,
                                                    int k_cat) {
  const int n = base_eta.size();
  if (n < 32 || person_nodes.ncol() < 2 || k_cat < 2 ||
      person_nodes.nrow() == 0) {
    return {};
  }
  for (int q = 0; q < person_nodes.ncol(); ++q) {
    const double first = person_nodes(0, q);
    if (!std::isfinite(first)) return {};
    for (int p = 1; p < person_nodes.nrow(); ++p) {
      const double node = person_nodes(p, q);
      if (std::memcmp(&first, &node, sizeof(double)) != 0) return {};
    }
  }
  if (model_code == 0) {
    for (int k = 0; k < k_cat; ++k) {
      if (!std::isfinite(step_rsm[k])) return {};
    }
  } else {
    for (size_t idx = 0; idx < step_rm.size(); ++idx) {
      if (!std::isfinite(step_rm[idx])) return {};
    }
    if (model_code == 2) {
      for (int s = 0; s < n_slopes; ++s) {
        if (!std::isfinite(slopes_ptr[s])) return {};
      }
    }
  }

  std::unordered_map<std::string, int> first_rows;
  ProbabilityRows rows;
  rows.group.resize(static_cast<size_t>(n));
  std::string key;
  const size_t key_doubles = 1 + (model_code == 2 ? 1 : 0) +
                             (model_code != 0 ? static_cast<size_t>(k_cat) : 0);
  key.reserve(sizeof(double) * key_doubles);
  for (int i = 0; i < n; ++i) {
    const double eta = base_eta[i];
    if (!std::isfinite(eta)) return {};
    key.assign(reinterpret_cast<const char*>(&eta), sizeof(double));
    if (model_code == 2) {
      const int sl_i = slope_idx_ptr[i] - 1;
      if (sl_i < 0 || sl_i >= n_slopes) return {};
      const double sl = slopes_ptr[sl_i];
      key.append(reinterpret_cast<const char*>(&sl), sizeof(double));
    }
    if (model_code != 0) {
      const int st_i = step_idx_ptr[i] - 1;
      if (st_i < 0 || st_i >= n_criteria) return {};
      const double* st_row = step_rm.data() + static_cast<size_t>(st_i) * k_cat;
      key.append(reinterpret_cast<const char*>(st_row), sizeof(double) * static_cast<size_t>(k_cat));
    }
    const auto inserted = first_rows.emplace(key, static_cast<int>(rows.first.size()));
    if (inserted.second) rows.first.push_back(i);
    rows.group[static_cast<size_t>(i)] = inserted.first->second;
    if (first_rows.size() > static_cast<size_t>(n / 8)) return {};
  }
  if (rows.first.size() * static_cast<size_t>(k_cat) >
      static_cast<size_t>(n) * person_nodes.ncol()) return {};
  return rows;
}

inline writable::doubles_matrix<> build_p_geq(const writable::doubles_matrix<>& probs) {
  const int n = probs.nrow();
  const int k_cat = probs.ncol();
  const int n_steps = k_cat - 1;
  writable::doubles_matrix<> out(n, std::max(n_steps, 0));
  if (n == 0 || n_steps <= 0) {
    return out;
  }

  for (int i = 0; i < n; ++i) {
    out(i, n_steps - 1) = probs(i, k_cat - 1);
  }
  for (int s = n_steps - 2; s >= 0; --s) {
    for (int i = 0; i < n; ++i) {
      out(i, s) = out(i, s + 1) + probs(i, s + 1);
    }
  }

  return out;
}

// Helper for 1-person standardized adaptive posterior evaluation at scalar z.
struct AdaptivePointMoments {
  double log_likelihood = 0.0;
  double score = 0.0;
  double information = 1.0;
  bool finite_ok = true;
};

inline AdaptivePointMoments eval_adaptive_point_moments(
    double z,
    double mu_p,
    double sigma,
    int m_obs,
    int k_cat,
    int model_code,
    const double* base_eta_sub,
    const int* score_sub,
    const double* weight_sub,
    const double* slope_sub,
    const double* step_rsm,
    const double* step_rm,
    const int* step_idx_sub,
    std::vector<double>& log_num_buf,
    std::vector<double>& prob_buf) {
  AdaptivePointMoments res;
  const double theta = mu_p + sigma * z;
  double score_sum = 0.0;
  double info_sum = 0.0;
  double ll_sum = 0.0;

  for (int r = 0; r < m_obs; ++r) {
    const double eta = base_eta_sub[r] + theta;
    const double sl = (model_code == 2) ? slope_sub[r] : 1.0;
    const double* st_ptr = (model_code == 0)
      ? step_rsm
      : (step_rm + static_cast<size_t>(step_idx_sub[r] - 1) * k_cat);

    double row_max = -std::numeric_limits<double>::infinity();
    for (int k = 0; k < k_cat; ++k) {
      const double lin = eta * static_cast<double>(k) - st_ptr[k];
      const double val = (model_code == 2) ? (lin * sl) : lin;
      log_num_buf[static_cast<size_t>(k)] = val;
      if (val > row_max) row_max = val;
    }

    double denom = 0.0;
    for (int k = 0; k < k_cat; ++k) {
      const double e = std::exp(log_num_buf[static_cast<size_t>(k)] - row_max);
      prob_buf[static_cast<size_t>(k)] = e;
      denom += e;
    }
    const double log_denom = std::log(denom);
    const int y = score_sub[r];
    const double w = weight_sub[r];
    ll_sum += w * ((log_num_buf[static_cast<size_t>(y)] - row_max) - log_denom);

    double expected = 0.0;
    for (int k = 0; k < k_cat; ++k) {
      prob_buf[static_cast<size_t>(k)] /= denom;
      expected += prob_buf[static_cast<size_t>(k)] * static_cast<double>(k);
    }
    double variance = 0.0;
    for (int k = 0; k < k_cat; ++k) {
      const double centered = static_cast<double>(k) - expected;
      variance += prob_buf[static_cast<size_t>(k)] * centered * centered;
    }
    score_sum += (w * sl) * (static_cast<double>(y) - expected);
    info_sum += (w * sl * sl) * variance;
  }

  const double m1 = sigma * score_sum;
  const double m2 = sigma * sigma * info_sum;
  if (!std::isfinite(m1) || !std::isfinite(m2)) {
    res.finite_ok = false;
    return res;
  }
  res.log_likelihood = ll_sum;
  res.score = m1 - z;
  res.information = m2 + 1.0;
  return res;
}

struct AdaptiveModeResult {
  double mode = 0.0;
  double scale = 1.0;
  double information = 1.0;
  int error_code = 0;
  std::string error_message;
};

// Safeguarded Newton-Raphson with bracketed bisection fallback on the strictly
// decreasing standardized posterior score S(z), matching uniroot(tol = 1e-12)
// while reaching machine-precision stationarity in ~4-6 evaluations.
inline AdaptiveModeResult find_adaptive_posterior_mode(
    double mu_p,
    double sigma,
    int m_obs,
    int k_cat,
    int model_code,
    const double* base_eta_sub,
    const int* score_sub,
    const double* weight_sub,
    const double* slope_sub,
    const double* step_rsm,
    const double* step_rm,
    const int* step_idx_sub,
    std::vector<double>& log_num_buf,
    std::vector<double>& prob_buf) {
  AdaptiveModeResult out;
  auto eval_z = [&](double z) {
    return eval_adaptive_point_moments(
      z, mu_p, sigma, m_obs, k_cat, model_code,
      base_eta_sub, score_sub, weight_sub, slope_sub,
      step_rsm, step_rm, step_idx_sub, log_num_buf, prob_buf);
  };

  AdaptivePointMoments m0 = eval_z(0.0);
  if (!m0.finite_ok) {
    out.error_code = 1;
    out.error_message = "Non-finite adaptive posterior score or information.";
    return out;
  }
  if (std::abs(m0.score) <= 1e-15) {
    const double sc = 1.0 / std::sqrt(m0.information);
    if (!std::isfinite(sc) || sc <= 0.0) {
      out.error_code = 2;
      out.error_message = "Invalid local posterior scale.";
      return out;
    }
    out.mode = 0.0;
    out.scale = sc;
    out.information = m0.information;
    return out;
  }

  double low = -1.0;
  double high = 1.0;
  double z = 0.0;
  AdaptivePointMoments cur = m0;
  bool bracketed = false;

  if (m0.score > 0.0) {
    low = 0.0;
    double step = 1.0;
    high = 1.0;
    for (int expand = 0; expand < 64; ++expand) {
      AdaptivePointMoments mh = eval_z(high);
      if (!mh.finite_ok) {
        out.error_code = 1;
        out.error_message = "Non-finite adaptive posterior score or information.";
        return out;
      }
      if (mh.score <= 0.0) { bracketed = true; break; }
      low = high;
      step *= 2.0;
      high += step;
    }
  } else {
    high = 0.0;
    double step = 1.0;
    low = -1.0;
    for (int expand = 0; expand < 64; ++expand) {
      AdaptivePointMoments ml = eval_z(low);
      if (!ml.finite_ok) {
        out.error_code = 1;
        out.error_message = "Non-finite adaptive posterior score or information.";
        return out;
      }
      if (ml.score >= 0.0) { bracketed = true; break; }
      high = low;
      step *= 2.0;
      low -= step;
    }
  }

  // Preserve the R solver as a fallback when the native search cannot certify
  // a bracket or termination; a finite endpoint is not a posterior mode.
  if (!bracketed) {
    out.error_code = 5;
    out.error_message = "Unable to bracket adaptive posterior mode.";
    return out;
  }

  // Start Newton from the initial step at z = 0, safeguarded inside (low, high)
  double z_next = z + cur.score / cur.information;
  if (!std::isfinite(z_next) || z_next <= low || z_next >= high) {
    z_next = 0.5 * (low + high);
  }
  z = z_next;

  bool converged = false;
  for (int iter = 0; iter < 80; ++iter) {
    cur = eval_z(z);
    if (!cur.finite_ok) {
      out.error_code = 1;
      out.error_message = "Non-finite adaptive posterior score or information.";
      return out;
    }
    if (std::abs(cur.score) <= 1e-14) {
      converged = true;
      break;
    }
    if (cur.score > 0.0) {
      low = z;
    } else {
      high = z;
    }
    const double scale_z = std::max(1.0, std::abs(z));
    if ((high - low) <= 4.0 * std::numeric_limits<double>::epsilon() * scale_z) {
      converged = true;
      break;
    }
    double cand = z + cur.score / cur.information;
    if (!std::isfinite(cand) || cand <= low || cand >= high) {
      cand = 0.5 * (low + high);
    }
    const double dz = std::abs(cand - z);
    z = cand;
    if (dz <= 1e-14 * scale_z && std::abs(cur.score) <= 1e-12) {
      cur = eval_z(z);
      if (!cur.finite_ok) {
        out.error_code = 1;
        out.error_message = "Non-finite adaptive posterior score or information.";
        return out;
      }
      converged = true;
      break;
    }
  }

  if (!converged) {
    out.error_code = 5;
    out.error_message = "Adaptive posterior mode did not converge.";
    return out;
  }

  const double sc = 1.0 / std::sqrt(cur.information);
  if (!std::isfinite(sc) || sc <= 0.0) {
    out.error_code = 2;
    out.error_message = "Invalid local posterior scale.";
    return out;
  }
  out.mode = z;
  out.scale = sc;
  out.information = cur.information;
  return out;
}

inline void sparse_crossprod_add(int pat,
                                 const int* pat_col_starts,
                                 const int* active_col,
                                 const int* col_val_starts,
                                 const int* row_i,
                                 const double* val_x,
                                 const double* vec,
                                 double alpha,
                                 std::vector<double>& target) {
  const int c_begin = pat_col_starts[pat];
  const int c_end = pat_col_starts[pat + 1];
  for (int c_pos = c_begin; c_pos < c_end; ++c_pos) {
    const int col = active_col[c_pos];
    const int v_begin = col_val_starts[c_pos];
    const int v_end = col_val_starts[c_pos + 1];
    double acc = 0.0;
    for (int idx = v_begin; idx < v_end; ++idx) {
      acc += val_x[idx] * vec[row_i[idx]];
    }
    target[static_cast<size_t>(col)] += alpha * acc;
  }
}

}  // namespace

[[cpp11::register]]
writable::doubles_matrix<> mfrm_cpp_compute_p_geq(doubles_matrix<> probs) {
  const int n = probs.nrow();
  const int k_cat = probs.ncol();
  const int n_steps = k_cat - 1;
  writable::doubles_matrix<> out(n, std::max(n_steps, 0));
  if (n == 0 || n_steps <= 0) {
    return out;
  }

  for (int i = 0; i < n; ++i) {
    out(i, n_steps - 1) = probs(i, k_cat - 1);
  }
  for (int s = n_steps - 2; s >= 0; --s) {
    for (int i = 0; i < n; ++i) {
      out(i, s) = out(i, s + 1) + probs(i, s + 1);
    }
  }

  return out;
}

[[cpp11::register]]
writable::list mfrm_cpp_rsm_logprob_bundle(integers score_k,
                                           integers person_int,
                                           doubles base_eta,
                                           doubles_matrix<> person_nodes,
                                           doubles step_cum,
                                           sexp weight,
                                           bool include_probs) {
  const int n = base_eta.size();
  const int n_nodes = person_nodes.ncol();
  const int k_cat = step_cum.size();
  const bool has_weight = !Rf_isNull(weight);
  doubles weight_vec;
  if (has_weight) {
    weight_vec = as_doubles(weight);
  }

  writable::doubles_matrix<> log_prob_mat(n, n_nodes);
  writable::list prob_list;
  if (include_probs) {
    prob_list = writable::list(n_nodes);
  }

  std::vector<double> log_num(static_cast<size_t>(k_cat));
  std::vector<double> prob_row(static_cast<size_t>(k_cat));
  const ProbabilityRows probability_rows = probability_representatives(
    base_eta, person_nodes, REAL(step_cum), k_cat, false);
  const bool reuse_probabilities = !probability_rows.group.empty();
  std::vector<double> cached_log_num(probability_rows.first.size() * static_cast<size_t>(k_cat));
  std::vector<double> cached_row_max(probability_rows.first.size());
  std::vector<double> cached_log_denom(probability_rows.first.size());

  for (int q = 0; q < n_nodes; ++q) {
    if (include_probs) {
      writable::doubles_matrix<> probs_q(n, k_cat);
      for (int i = 0; i < n; ++i) {
        const int person_row = person_int[i] - 1;
        const int obs_col = score_k[i];
        const int group = reuse_probabilities ? probability_rows.group[static_cast<size_t>(i)] : 0;
        const int source_row = reuse_probabilities ? probability_rows.first[static_cast<size_t>(group)] : i;
        if (source_row != i) {
          double lp = (cached_log_num[static_cast<size_t>(group) * k_cat + obs_col] -
                       cached_row_max[static_cast<size_t>(group)]) -
                      cached_log_denom[static_cast<size_t>(group)];
          if (has_weight) lp *= weight_vec[i];
          log_prob_mat(i, q) = lp;
          for (int k = 0; k < k_cat; ++k) {
            probs_q(i, k) = probs_q(source_row, k);
          }
          continue;
        }
        const double eta_q = base_eta[i] + person_nodes(person_row, q);
        double row_max = -std::numeric_limits<double>::infinity();

        for (int k = 0; k < k_cat; ++k) {
          const double value = eta_q * static_cast<double>(k) - step_cum[k];
          log_num[static_cast<size_t>(k)] = value;
          if (value > row_max) {
            row_max = value;
          }
        }

        double denom_sum = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          const double prob = std::exp(log_num[static_cast<size_t>(k)] - row_max);
          prob_row[static_cast<size_t>(k)] = prob;
          denom_sum += prob;
        }
        const double log_denom = std::log(denom_sum);
        if (reuse_probabilities) {
          for (int k = 0; k < k_cat; ++k) {
            cached_log_num[static_cast<size_t>(group) * k_cat + k] = log_num[static_cast<size_t>(k)];
          }
          cached_row_max[static_cast<size_t>(group)] = row_max;
          cached_log_denom[static_cast<size_t>(group)] = log_denom;
        }

        double lp = (log_num[static_cast<size_t>(obs_col)] - row_max) - log_denom;
        if (has_weight) {
          lp *= weight_vec[i];
        }
        log_prob_mat(i, q) = lp;
        for (int k = 0; k < k_cat; ++k) {
          probs_q(i, k) = prob_row[static_cast<size_t>(k)] / denom_sum;
        }
      }
      prob_list[q] = probs_q;
    } else {
      for (int i = 0; i < n; ++i) {
        const int person_row = person_int[i] - 1;
        const int obs_col = score_k[i];
        const int group = reuse_probabilities ? probability_rows.group[static_cast<size_t>(i)] : 0;
        const int source_row = reuse_probabilities ? probability_rows.first[static_cast<size_t>(group)] : i;
        if (source_row != i) {
          double lp = (cached_log_num[static_cast<size_t>(group) * k_cat + obs_col] -
                       cached_row_max[static_cast<size_t>(group)]) -
                      cached_log_denom[static_cast<size_t>(group)];
          if (has_weight) lp *= weight_vec[i];
          log_prob_mat(i, q) = lp;
          continue;
        }
        const double eta_q = base_eta[i] + person_nodes(person_row, q);
        double row_max = -std::numeric_limits<double>::infinity();

        for (int k = 0; k < k_cat; ++k) {
          const double value = eta_q * static_cast<double>(k) - step_cum[k];
          log_num[static_cast<size_t>(k)] = value;
          if (value > row_max) {
            row_max = value;
          }
        }

        double denom_sum = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          denom_sum += std::exp(log_num[static_cast<size_t>(k)] - row_max);
        }
        const double log_denom = std::log(denom_sum);
        if (reuse_probabilities) {
          for (int k = 0; k < k_cat; ++k) {
            cached_log_num[static_cast<size_t>(group) * k_cat + k] = log_num[static_cast<size_t>(k)];
          }
          cached_row_max[static_cast<size_t>(group)] = row_max;
          cached_log_denom[static_cast<size_t>(group)] = log_denom;
        }

        double lp = (log_num[static_cast<size_t>(obs_col)] - row_max) - log_denom;
        if (has_weight) {
          lp *= weight_vec[i];
        }
        log_prob_mat(i, q) = lp;
      }
    }
  }

  if (include_probs) {
    return writable::list({
      "log_prob_mat"_nm = log_prob_mat,
      "prob_list"_nm = prob_list
    });
  }

  return writable::list({
    "log_prob_mat"_nm = log_prob_mat
  });
}

[[cpp11::register]]
writable::list mfrm_cpp_pcm_logprob_bundle(integers score_k,
                                           integers person_int,
                                           doubles base_eta,
                                           doubles_matrix<> person_nodes,
                                           doubles_matrix<> step_cum_obs,
                                           sexp weight,
                                           bool include_probs) {
  const int n = base_eta.size();
  const int n_nodes = person_nodes.ncol();
  const int k_cat = step_cum_obs.ncol();
  const bool has_weight = !Rf_isNull(weight);
  doubles weight_vec;
  if (has_weight) {
    weight_vec = as_doubles(weight);
  }

  writable::doubles_matrix<> log_prob_mat(n, n_nodes);
  writable::list prob_list;
  if (include_probs) {
    prob_list = writable::list(n_nodes);
  }

  std::vector<double> log_num(static_cast<size_t>(k_cat));
  std::vector<double> prob_row(static_cast<size_t>(k_cat));
  const ProbabilityRows probability_rows = probability_representatives(
    base_eta, person_nodes, REAL(step_cum_obs), k_cat, true);
  const bool reuse_probabilities = !probability_rows.group.empty();
  std::vector<double> cached_log_num(probability_rows.first.size() * static_cast<size_t>(k_cat));
  std::vector<double> cached_row_max(probability_rows.first.size());
  std::vector<double> cached_log_denom(probability_rows.first.size());

  for (int q = 0; q < n_nodes; ++q) {
    if (include_probs) {
      writable::doubles_matrix<> probs_q(n, k_cat);
      for (int i = 0; i < n; ++i) {
        const int person_row = person_int[i] - 1;
        const int obs_col = score_k[i];
        const int group = reuse_probabilities ? probability_rows.group[static_cast<size_t>(i)] : 0;
        const int source_row = reuse_probabilities ? probability_rows.first[static_cast<size_t>(group)] : i;
        if (source_row != i) {
          double lp = (cached_log_num[static_cast<size_t>(group) * k_cat + obs_col] -
                       cached_row_max[static_cast<size_t>(group)]) -
                      cached_log_denom[static_cast<size_t>(group)];
          if (has_weight) lp *= weight_vec[i];
          log_prob_mat(i, q) = lp;
          for (int k = 0; k < k_cat; ++k) {
            probs_q(i, k) = probs_q(source_row, k);
          }
          continue;
        }
        const double eta_q = base_eta[i] + person_nodes(person_row, q);
        double row_max = -std::numeric_limits<double>::infinity();

        for (int k = 0; k < k_cat; ++k) {
          const double value = eta_q * static_cast<double>(k) - step_cum_obs(i, k);
          log_num[static_cast<size_t>(k)] = value;
          if (value > row_max) {
            row_max = value;
          }
        }

        double denom_sum = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          const double prob = std::exp(log_num[static_cast<size_t>(k)] - row_max);
          prob_row[static_cast<size_t>(k)] = prob;
          denom_sum += prob;
        }
        const double log_denom = std::log(denom_sum);
        if (reuse_probabilities) {
          for (int k = 0; k < k_cat; ++k) {
            cached_log_num[static_cast<size_t>(group) * k_cat + k] = log_num[static_cast<size_t>(k)];
          }
          cached_row_max[static_cast<size_t>(group)] = row_max;
          cached_log_denom[static_cast<size_t>(group)] = log_denom;
        }

        double lp = (log_num[static_cast<size_t>(obs_col)] - row_max) - log_denom;
        if (has_weight) {
          lp *= weight_vec[i];
        }
        log_prob_mat(i, q) = lp;
        for (int k = 0; k < k_cat; ++k) {
          probs_q(i, k) = prob_row[static_cast<size_t>(k)] / denom_sum;
        }
      }
      prob_list[q] = probs_q;
    } else {
      for (int i = 0; i < n; ++i) {
        const int person_row = person_int[i] - 1;
        const int obs_col = score_k[i];
        const int group = reuse_probabilities ? probability_rows.group[static_cast<size_t>(i)] : 0;
        const int source_row = reuse_probabilities ? probability_rows.first[static_cast<size_t>(group)] : i;
        if (source_row != i) {
          double lp = (cached_log_num[static_cast<size_t>(group) * k_cat + obs_col] -
                       cached_row_max[static_cast<size_t>(group)]) -
                      cached_log_denom[static_cast<size_t>(group)];
          if (has_weight) lp *= weight_vec[i];
          log_prob_mat(i, q) = lp;
          continue;
        }
        const double eta_q = base_eta[i] + person_nodes(person_row, q);
        double row_max = -std::numeric_limits<double>::infinity();

        for (int k = 0; k < k_cat; ++k) {
          const double value = eta_q * static_cast<double>(k) - step_cum_obs(i, k);
          log_num[static_cast<size_t>(k)] = value;
          if (value > row_max) {
            row_max = value;
          }
        }

        double denom_sum = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          denom_sum += std::exp(log_num[static_cast<size_t>(k)] - row_max);
        }
        const double log_denom = std::log(denom_sum);
        if (reuse_probabilities) {
          for (int k = 0; k < k_cat; ++k) {
            cached_log_num[static_cast<size_t>(group) * k_cat + k] = log_num[static_cast<size_t>(k)];
          }
          cached_row_max[static_cast<size_t>(group)] = row_max;
          cached_log_denom[static_cast<size_t>(group)] = log_denom;
        }

        double lp = (log_num[static_cast<size_t>(obs_col)] - row_max) - log_denom;
        if (has_weight) {
          lp *= weight_vec[i];
        }
        log_prob_mat(i, q) = lp;
      }
    }
  }

  if (include_probs) {
    return writable::list({
      "log_prob_mat"_nm = log_prob_mat,
      "prob_list"_nm = prob_list
    });
  }

  return writable::list({
    "log_prob_mat"_nm = log_prob_mat
  });
}

[[cpp11::register]]
writable::list mfrm_cpp_gpcm_logprob_bundle(integers score_k,
                                            integers person_int,
                                            doubles base_eta,
                                            doubles_matrix<> person_nodes,
                                            doubles_matrix<> step_cum_obs,
                                            doubles slope_obs,
                                            sexp weight,
                                            bool include_probs,
                                            bool include_linear_part) {
  const int n = base_eta.size();
  const int n_nodes = person_nodes.ncol();
  const int k_cat = step_cum_obs.ncol();
  const bool has_weight = !Rf_isNull(weight);
  doubles weight_vec;
  if (has_weight) {
    weight_vec = as_doubles(weight);
  }

  writable::doubles_matrix<> log_prob_mat(n, n_nodes);
  writable::list prob_list;
  writable::list linear_part_list;
  if (include_probs) {
    prob_list = writable::list(n_nodes);
  }
  if (include_linear_part) {
    linear_part_list = writable::list(n_nodes);
  }

  std::vector<double> linear_row(static_cast<size_t>(k_cat));
  std::vector<double> log_num(static_cast<size_t>(k_cat));
  std::vector<double> prob_row(static_cast<size_t>(k_cat));
  const ProbabilityRows probability_rows = probability_representatives(
    base_eta, person_nodes, REAL(step_cum_obs), k_cat, true, REAL(slope_obs));
  const bool reuse_probabilities = !probability_rows.group.empty();
  std::vector<double> cached_log_num(probability_rows.first.size() * static_cast<size_t>(k_cat));
  std::vector<double> cached_row_max(probability_rows.first.size());
  std::vector<double> cached_log_denom(probability_rows.first.size());

  for (int q = 0; q < n_nodes; ++q) {
    if (include_probs || include_linear_part) {
      writable::doubles_matrix<> probs_q(include_probs ? n : 0, include_probs ? k_cat : 0);
      writable::doubles_matrix<> linear_q(include_linear_part ? n : 0, include_linear_part ? k_cat : 0);
      for (int i = 0; i < n; ++i) {
        const int person_row = person_int[i] - 1;
        const int obs_col = score_k[i];
        const int group = reuse_probabilities ? probability_rows.group[static_cast<size_t>(i)] : 0;
        const int source_row = reuse_probabilities ? probability_rows.first[static_cast<size_t>(group)] : i;
        if (source_row != i) {
          double lp = (cached_log_num[static_cast<size_t>(group) * k_cat + obs_col] -
                       cached_row_max[static_cast<size_t>(group)]) -
                      cached_log_denom[static_cast<size_t>(group)];
          if (has_weight) lp *= weight_vec[i];
          log_prob_mat(i, q) = lp;
          if (include_probs) {
            for (int k = 0; k < k_cat; ++k) {
              probs_q(i, k) = probs_q(source_row, k);
            }
          }
          if (include_linear_part) {
            for (int k = 0; k < k_cat; ++k) {
              linear_q(i, k) = linear_q(source_row, k);
            }
          }
          continue;
        }
        const double eta_q = base_eta[i] + person_nodes(person_row, q);
        const double sl = slope_obs[i];
        double row_max = -std::numeric_limits<double>::infinity();

        for (int k = 0; k < k_cat; ++k) {
          const double lin = eta_q * static_cast<double>(k) - step_cum_obs(i, k);
          const double value = lin * sl;
          linear_row[static_cast<size_t>(k)] = lin;
          log_num[static_cast<size_t>(k)] = value;
          if (value > row_max) {
            row_max = value;
          }
        }

        double denom_sum = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          const double prob = std::exp(log_num[static_cast<size_t>(k)] - row_max);
          prob_row[static_cast<size_t>(k)] = prob;
          denom_sum += prob;
        }
        const double log_denom = std::log(denom_sum);
        if (reuse_probabilities) {
          for (int k = 0; k < k_cat; ++k) {
            cached_log_num[static_cast<size_t>(group) * k_cat + k] = log_num[static_cast<size_t>(k)];
          }
          cached_row_max[static_cast<size_t>(group)] = row_max;
          cached_log_denom[static_cast<size_t>(group)] = log_denom;
        }

        double lp = (log_num[static_cast<size_t>(obs_col)] - row_max) - log_denom;
        if (has_weight) {
          lp *= weight_vec[i];
        }
        log_prob_mat(i, q) = lp;
        if (include_probs) {
          for (int k = 0; k < k_cat; ++k) {
            probs_q(i, k) = prob_row[static_cast<size_t>(k)] / denom_sum;
          }
        }
        if (include_linear_part) {
          for (int k = 0; k < k_cat; ++k) {
            linear_q(i, k) = linear_row[static_cast<size_t>(k)];
          }
        }
      }
      if (include_probs) prob_list[q] = probs_q;
      if (include_linear_part) linear_part_list[q] = linear_q;
    } else {
      for (int i = 0; i < n; ++i) {
        const int person_row = person_int[i] - 1;
        const int obs_col = score_k[i];
        const int group = reuse_probabilities ? probability_rows.group[static_cast<size_t>(i)] : 0;
        const int source_row = reuse_probabilities ? probability_rows.first[static_cast<size_t>(group)] : i;
        if (source_row != i) {
          double lp = (cached_log_num[static_cast<size_t>(group) * k_cat + obs_col] -
                       cached_row_max[static_cast<size_t>(group)]) -
                      cached_log_denom[static_cast<size_t>(group)];
          if (has_weight) lp *= weight_vec[i];
          log_prob_mat(i, q) = lp;
          continue;
        }
        const double eta_q = base_eta[i] + person_nodes(person_row, q);
        const double sl = slope_obs[i];
        double row_max = -std::numeric_limits<double>::infinity();

        for (int k = 0; k < k_cat; ++k) {
          const double lin = eta_q * static_cast<double>(k) - step_cum_obs(i, k);
          const double value = lin * sl;
          log_num[static_cast<size_t>(k)] = value;
          if (value > row_max) {
            row_max = value;
          }
        }

        double denom_sum = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          denom_sum += std::exp(log_num[static_cast<size_t>(k)] - row_max);
        }
        const double log_denom = std::log(denom_sum);
        if (reuse_probabilities) {
          for (int k = 0; k < k_cat; ++k) {
            cached_log_num[static_cast<size_t>(group) * k_cat + k] = log_num[static_cast<size_t>(k)];
          }
          cached_row_max[static_cast<size_t>(group)] = row_max;
          cached_log_denom[static_cast<size_t>(group)] = log_denom;
        }

        double lp = (log_num[static_cast<size_t>(obs_col)] - row_max) - log_denom;
        if (has_weight) {
          lp *= weight_vec[i];
        }
        log_prob_mat(i, q) = lp;
      }
    }
  }

  if (include_probs && include_linear_part) {
    return writable::list({
      "log_prob_mat"_nm = log_prob_mat,
      "prob_list"_nm = prob_list,
      "linear_part_list"_nm = linear_part_list
    });
  }
  if (include_probs) {
    return writable::list({
      "log_prob_mat"_nm = log_prob_mat,
      "prob_list"_nm = prob_list
    });
  }
  if (include_linear_part) {
    return writable::list({
      "log_prob_mat"_nm = log_prob_mat,
      "linear_part_list"_nm = linear_part_list
    });
  }
  return writable::list({
    "log_prob_mat"_nm = log_prob_mat
  });
}

// Fused fixed-grid MML (and EM M-step) likelihood, posterior, and gradient
// sufficient-statistics kernel for RSM (0), PCM (1), and GPCM (2).
// Eliminates R allocation of per-node n x K probability matrices, linear-part
// matrices, and n x (K-1) P_geq matrices during optimization and optimHess.
[[cpp11::register]]
writable::list mfrm_cpp_mml_eval_bundle(int model_code,
                                        integers score_k,
                                        integers person_int,
                                        doubles base_eta,
                                        doubles_matrix<> person_nodes,
                                        doubles_matrix<> person_log_weights,
                                        doubles step_cum_rsm,
                                        doubles_matrix<> step_cum_mat,
                                        integers step_idx,
                                        doubles slopes,
                                        integers slope_idx,
                                        sexp weight,
                                        sexp external_obs_posterior,
                                        bool include_grad,
                                        bool pop_active,
                                        doubles_matrix<> pop_design_person,
                                        doubles quad_std_nodes,
                                        double sigma_pop) {
  const int n = base_eta.size();
  const int n_nodes = person_nodes.ncol();
  const int k_cat = (model_code == 0) ? step_cum_rsm.size() : step_cum_mat.ncol();
  const int n_steps = std::max(k_cat - 1, 0);
  const int n_criteria = (model_code == 0) ? 1 : step_cum_mat.nrow();
  const int n_slopes = (model_code == 2) ? slopes.size() : 0;
  const bool has_weight = !Rf_isNull(weight);
  const bool has_ext_post = !Rf_isNull(external_obs_posterior);
  doubles weight_vec;
  if (has_weight) weight_vec = as_doubles(weight);
  doubles ext_post_vec;
  const double* ext_post_ptr = nullptr;
  if (has_ext_post) {
    ext_post_vec = as_doubles(external_obs_posterior);
    ext_post_ptr = REAL(ext_post_vec);
  }

  // Build row-major step table for cache-local category traversal.
  std::vector<double> step_rm;
  if (model_code != 0) {
    step_rm.resize(static_cast<size_t>(n_criteria) * k_cat);
    for (int c = 0; c < n_criteria; ++c) {
      for (int k = 0; k < k_cat; ++k) {
        step_rm[static_cast<size_t>(c) * k_cat + k] = step_cum_mat(c, k);
      }
    }
  }
  const double* step_rsm_ptr = (model_code == 0 && k_cat > 0) ? REAL(step_cum_rsm) : nullptr;
  const int* step_idx_ptr = (model_code != 0 && n > 0) ? INTEGER(step_idx) : nullptr;
  const double* slopes_ptr = (model_code == 2 && n_slopes > 0) ? REAL(slopes) : nullptr;
  const int* slope_idx_ptr = (model_code == 2 && n > 0) ? INTEGER(slope_idx) : nullptr;

  // Map person_int (1-based) to first-observed person rows (matching rowsum(..., reorder = FALSE)).
  const int n_person_rows = person_nodes.nrow();
  std::vector<int> person_to_row(static_cast<size_t>(n_person_rows), -1);
  std::vector<int> person_ids_vec;
  person_ids_vec.reserve(static_cast<size_t>(n_person_rows));
  std::vector<int> obs_person_row(static_cast<size_t>(n), 0);
  if (!has_ext_post) {
    for (int i = 0; i < n; ++i) {
      const int p_id = person_int[i];
      int r = person_to_row[static_cast<size_t>(p_id - 1)];
      if (r < 0) {
        r = static_cast<int>(person_ids_vec.size());
        person_to_row[static_cast<size_t>(p_id - 1)] = r;
        person_ids_vec.push_back(p_id);
      }
      obs_person_row[static_cast<size_t>(i)] = r;
    }
  }
  const int n_active_persons = static_cast<int>(person_ids_vec.size());

  const ProbabilityRows rep_rows = probability_representatives_indexed(
    model_code, base_eta, person_nodes, step_rsm_ptr, step_rm, n_criteria,
    step_idx_ptr, slopes_ptr, n_slopes, slope_idx_ptr, k_cat);
  const bool reuse_rep = !rep_rows.group.empty();
  const size_t n_rep = rep_rows.first.size();

  // When representative deduplication applies, precompute category log-probs,
  // expected category, cumulative P_geq, and (for GPCM) category slope contrasts
  // across the G * Q unique pairs once.
  std::vector<double> rep_lp;
  std::vector<double> rep_expected;
  std::vector<double> rep_p_geq;
  std::vector<double> rep_slope_contrast;
  std::vector<double> linear_buf(static_cast<size_t>(k_cat));
  std::vector<double> log_num_buf(static_cast<size_t>(k_cat));
  std::vector<double> prob_buf(static_cast<size_t>(k_cat));
  std::vector<double> pgeq_buf(static_cast<size_t>(std::max(n_steps, 1)));

  if (reuse_rep) {
    rep_lp.resize(n_rep * static_cast<size_t>(n_nodes) * k_cat);
    if (include_grad) {
      rep_expected.resize(n_rep * static_cast<size_t>(n_nodes));
      if (n_steps > 0) {
        rep_p_geq.resize(n_rep * static_cast<size_t>(n_nodes) * n_steps);
      }
      if (model_code == 2) {
        rep_slope_contrast.resize(n_rep * static_cast<size_t>(n_nodes) * k_cat);
      }
    }
    for (int q = 0; q < n_nodes; ++q) {
      const double node_q = person_nodes(0, q);
      for (size_t g = 0; g < n_rep; ++g) {
        const int i0 = rep_rows.first[g];
        const double eta_q = base_eta[i0] + node_q;
        const double sl = (model_code == 2) ? slopes_ptr[slope_idx_ptr[i0] - 1] : 1.0;
        const double* st_ptr = (model_code == 0)
          ? step_rsm_ptr
          : (step_rm.data() + static_cast<size_t>(step_idx_ptr[i0] - 1) * k_cat);

        double row_max = -std::numeric_limits<double>::infinity();
        for (int k = 0; k < k_cat; ++k) {
          const double lin = eta_q * static_cast<double>(k) - st_ptr[k];
          const double val = (model_code == 2) ? (lin * sl) : lin;
          linear_buf[static_cast<size_t>(k)] = lin;
          log_num_buf[static_cast<size_t>(k)] = val;
          if (val > row_max) row_max = val;
        }
        double denom_sum = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          const double e = std::exp(log_num_buf[static_cast<size_t>(k)] - row_max);
          prob_buf[static_cast<size_t>(k)] = e;
          denom_sum += e;
        }
        const double log_denom = std::log(denom_sum);
        const size_t gq = g * static_cast<size_t>(n_nodes) + static_cast<size_t>(q);
        double* lp_dst = rep_lp.data() + gq * k_cat;
        for (int k = 0; k < k_cat; ++k) {
          lp_dst[k] = (log_num_buf[static_cast<size_t>(k)] - row_max) - log_denom;
        }
        if (include_grad) {
          double expected = 0.0;
          for (int k = 0; k < k_cat; ++k) {
            prob_buf[static_cast<size_t>(k)] /= denom_sum;
            expected += prob_buf[static_cast<size_t>(k)] * static_cast<double>(k);
          }
          rep_expected[gq] = expected;
          if (n_steps > 0) {
            double* pg_dst = rep_p_geq.data() + gq * n_steps;
            pg_dst[n_steps - 1] = prob_buf[static_cast<size_t>(k_cat - 1)];
            for (int s = n_steps - 2; s >= 0; --s) {
              pg_dst[s] = pg_dst[s + 1] + prob_buf[static_cast<size_t>(s + 1)];
            }
          }
          if (model_code == 2) {
            double* sc_dst = rep_slope_contrast.data() + gq * k_cat;
            for (int c = 0; c < k_cat; ++c) {
              const double obs_lin = linear_buf[static_cast<size_t>(c)];
              double acc = 0.0;
              for (int k = 0; k < k_cat; ++k) {
                acc += prob_buf[static_cast<size_t>(k)] * (obs_lin - linear_buf[static_cast<size_t>(k)]);
              }
              sc_dst[c] = sl * acc;
            }
          }
        }
      }
    }
  }

  // Pass 1: accumulate person log-likelihoods (or complete-data objective for EM M-step).
  double neg_loglik = 0.0;
  std::vector<double> ll_by_person;
  if (!has_ext_post) {
    ll_by_person.assign(static_cast<size_t>(n_active_persons) * n_nodes, 0.0);
  }

  for (int q = 0; q < n_nodes; ++q) {
    if (reuse_rep) {
      for (int i = 0; i < n; ++i) {
        const size_t g = static_cast<size_t>(rep_rows.group[static_cast<size_t>(i)]);
        const size_t gq = g * static_cast<size_t>(n_nodes) + static_cast<size_t>(q);
        double lp = rep_lp[gq * k_cat + static_cast<size_t>(score_k[i])];
        if (has_weight) lp *= weight_vec[i];
        if (has_ext_post) {
          neg_loglik -= lp * ext_post_ptr[static_cast<size_t>(q) * n + i];
        } else {
          const int r = obs_person_row[static_cast<size_t>(i)];
          ll_by_person[static_cast<size_t>(r) * n_nodes + q] += lp;
        }
      }
    } else {
      for (int i = 0; i < n; ++i) {
        const int p_row = person_int[i] - 1;
        const double eta_q = base_eta[i] + person_nodes(p_row, q);
        const double sl = (model_code == 2) ? slopes_ptr[slope_idx_ptr[i] - 1] : 1.0;
        const double* st_ptr = (model_code == 0)
          ? step_rsm_ptr
          : (step_rm.data() + static_cast<size_t>(step_idx_ptr[i] - 1) * k_cat);

        double row_max = -std::numeric_limits<double>::infinity();
        for (int k = 0; k < k_cat; ++k) {
          const double lin = eta_q * static_cast<double>(k) - st_ptr[k];
          const double val = (model_code == 2) ? (lin * sl) : lin;
          log_num_buf[static_cast<size_t>(k)] = val;
          if (val > row_max) row_max = val;
        }
        double denom_sum = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          denom_sum += std::exp(log_num_buf[static_cast<size_t>(k)] - row_max);
        }
        double lp = (log_num_buf[static_cast<size_t>(score_k[i])] - row_max) - std::log(denom_sum);
        if (has_weight) lp *= weight_vec[i];
        if (has_ext_post) {
          neg_loglik -= lp * ext_post_ptr[static_cast<size_t>(q) * n + i];
        } else {
          const int r = obs_person_row[static_cast<size_t>(i)];
          ll_by_person[static_cast<size_t>(r) * n_nodes + q] += lp;
        }
      }
    }
  }

  writable::doubles log_marginal_out(has_ext_post ? 0 : n_active_persons);
  writable::integers person_ids_out(has_ext_post ? 0 : n_active_persons);
  std::vector<double> posterior_flat;
  if (!has_ext_post) {
    if (include_grad) {
      posterior_flat.resize(static_cast<size_t>(n_active_persons) * n_nodes);
    }
    std::vector<double> shifted_row(static_cast<size_t>(n_nodes));
    for (int r = 0; r < n_active_persons; ++r) {
      person_ids_out[r] = person_ids_vec[static_cast<size_t>(r)];
      const int p_id0 = person_ids_vec[static_cast<size_t>(r)] - 1;
      const double* ll_row = ll_by_person.data() + static_cast<size_t>(r) * n_nodes;

      double offset = -std::numeric_limits<double>::infinity();
      for (int q = 0; q < n_nodes; ++q) {
        if (ll_row[q] > offset) offset = ll_row[q];
      }
      double r_max = -std::numeric_limits<double>::infinity();
      for (int q = 0; q < n_nodes; ++q) {
        const double sh = (ll_row[q] - offset) + person_log_weights(p_id0, q);
        shifted_row[static_cast<size_t>(q)] = sh;
        if (sh > r_max) r_max = sh;
      }
      double sum_exp = 0.0;
      for (int q = 0; q < n_nodes; ++q) {
        shifted_row[static_cast<size_t>(q)] -= r_max;
        sum_exp += std::exp(shifted_row[static_cast<size_t>(q)]);
      }
      const double log_norm = std::log(sum_exp);
      const double lm = offset + (r_max + log_norm);
      log_marginal_out[r] = lm;
      neg_loglik -= lm;
      if (include_grad) {
        double* post_row = posterior_flat.data() + static_cast<size_t>(r) * n_nodes;
        for (int q = 0; q < n_nodes; ++q) {
          post_row[q] = std::exp(shifted_row[static_cast<size_t>(q)] - log_norm);
        }
      }
    }
  }

  if (!include_grad) {
    return writable::list({
      "value"_nm = neg_loglik,
      "log_marginal"_nm = log_marginal_out,
      "person_ids"_nm = person_ids_out
    });
  }

  // Pass 2: accumulate gradient sufficient statistics in identical q-first order.
  writable::doubles marginal_residual(n);
  for (int i = 0; i < n; ++i) marginal_residual[i] = 0.0;

  writable::doubles_matrix<> grad_step_mat(n_criteria, n_steps);
  for (int c = 0; c < n_criteria; ++c) {
    for (int s = 0; s < n_steps; ++s) {
      grad_step_mat(c, s) = 0.0;
    }
  }
  std::vector<double> marginal_step_residual;
  if (model_code != 0 && n_steps > 0) {
    marginal_step_residual.assign(static_cast<size_t>(n) * n_steps, 0.0);
  }

  writable::doubles grad_log_slope_exp(n_slopes);
  for (int s = 0; s < n_slopes; ++s) grad_log_slope_exp[s] = 0.0;
  std::vector<double> marginal_slope_score;
  if (model_code == 2) {
    marginal_slope_score.assign(static_cast<size_t>(n), 0.0);
  }

  const bool do_pop = pop_active && !has_ext_post;
  const int n_beta = do_pop ? pop_design_person.ncol() : 0;
  writable::doubles grad_beta(n_beta);
  for (int b = 0; b < n_beta; ++b) grad_beta[b] = 0.0;
  double grad_log_sigma2 = 0.0;
  std::vector<double> person_score_q(do_pop ? static_cast<size_t>(n_active_persons) : 0);

  for (int q = 0; q < n_nodes; ++q) {
    if (do_pop) {
      std::fill(person_score_q.begin(), person_score_q.end(), 0.0);
    }
    std::vector<double> rsm_step_q(model_code == 0 ? static_cast<size_t>(n_steps) : 0, 0.0);

    if (reuse_rep) {
      for (int i = 0; i < n; ++i) {
        const size_t g = static_cast<size_t>(rep_rows.group[static_cast<size_t>(i)]);
        const size_t gq = g * static_cast<size_t>(n_nodes) + static_cast<size_t>(q);
        const int y = score_k[i];
        const double sl = (model_code == 2) ? slopes_ptr[slope_idx_ptr[i] - 1] : 1.0;
        const double w = has_weight ? weight_vec[i] : 1.0;

        double res_q = (model_code == 2)
          ? (sl * (static_cast<double>(y) - rep_expected[gq]))
          : (static_cast<double>(y) - rep_expected[gq]);
        if (has_weight) res_q *= w;

        const double post_iq = has_ext_post
          ? ext_post_ptr[static_cast<size_t>(q) * n + i]
          : posterior_flat[static_cast<size_t>(obs_person_row[static_cast<size_t>(i)]) * n_nodes + q];

        if (do_pop) {
          person_score_q[static_cast<size_t>(obs_person_row[static_cast<size_t>(i)])] += res_q;
        }
        marginal_residual[i] += res_q * post_iq;

        if (n_steps > 0) {
          const double* pg = rep_p_geq.data() + gq * n_steps;
          if (model_code == 0) {
            const double scale_step = has_weight ? (post_iq * w) : post_iq;
            for (int s = 0; s < n_steps; ++s) {
              const double ind = (y >= (s + 1)) ? 1.0 : 0.0;
              rsm_step_q[static_cast<size_t>(s)] += (pg[s] - ind) * scale_step;
            }
          } else {
            double scale_step = (model_code == 2) ? (sl * post_iq) : post_iq;
            if (has_weight) scale_step *= w;
            double* msr = marginal_step_residual.data() + static_cast<size_t>(i) * n_steps;
            for (int s = 0; s < n_steps; ++s) {
              const double ind = (y >= (s + 1)) ? 1.0 : 0.0;
              msr[s] += (pg[s] - ind) * scale_step;
            }
          }
        }

        if (model_code == 2) {
          double sc = rep_slope_contrast[gq * k_cat + static_cast<size_t>(y)] * post_iq;
          if (has_weight) sc *= w;
          marginal_slope_score[static_cast<size_t>(i)] += sc;
        }
      }
    } else {
      for (int i = 0; i < n; ++i) {
        const int p_row = person_int[i] - 1;
        const int y = score_k[i];
        const double eta_q = base_eta[i] + person_nodes(p_row, q);
        const double sl = (model_code == 2) ? slopes_ptr[slope_idx_ptr[i] - 1] : 1.0;
        const double w = has_weight ? weight_vec[i] : 1.0;
        const double* st_ptr = (model_code == 0)
          ? step_rsm_ptr
          : (step_rm.data() + static_cast<size_t>(step_idx_ptr[i] - 1) * k_cat);

        double row_max = -std::numeric_limits<double>::infinity();
        for (int k = 0; k < k_cat; ++k) {
          const double lin = eta_q * static_cast<double>(k) - st_ptr[k];
          const double val = (model_code == 2) ? (lin * sl) : lin;
          linear_buf[static_cast<size_t>(k)] = lin;
          log_num_buf[static_cast<size_t>(k)] = val;
          if (val > row_max) row_max = val;
        }
        double denom_sum = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          const double e = std::exp(log_num_buf[static_cast<size_t>(k)] - row_max);
          prob_buf[static_cast<size_t>(k)] = e;
          denom_sum += e;
        }
        double expected = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          prob_buf[static_cast<size_t>(k)] /= denom_sum;
          expected += prob_buf[static_cast<size_t>(k)] * static_cast<double>(k);
        }

        double res_q = (model_code == 2)
          ? (sl * (static_cast<double>(y) - expected))
          : (static_cast<double>(y) - expected);
        if (has_weight) res_q *= w;

        const double post_iq = has_ext_post
          ? ext_post_ptr[static_cast<size_t>(q) * n + i]
          : posterior_flat[static_cast<size_t>(obs_person_row[static_cast<size_t>(i)]) * n_nodes + q];

        if (do_pop) {
          person_score_q[static_cast<size_t>(obs_person_row[static_cast<size_t>(i)])] += res_q;
        }
        marginal_residual[i] += res_q * post_iq;

        if (n_steps > 0) {
          pgeq_buf[static_cast<size_t>(n_steps - 1)] = prob_buf[static_cast<size_t>(k_cat - 1)];
          for (int s = n_steps - 2; s >= 0; --s) {
            pgeq_buf[static_cast<size_t>(s)] = pgeq_buf[static_cast<size_t>(s + 1)] +
                                               prob_buf[static_cast<size_t>(s + 1)];
          }
          if (model_code == 0) {
            const double scale_step = has_weight ? (post_iq * w) : post_iq;
            for (int s = 0; s < n_steps; ++s) {
              const double ind = (y >= (s + 1)) ? 1.0 : 0.0;
              rsm_step_q[static_cast<size_t>(s)] += (pgeq_buf[static_cast<size_t>(s)] - ind) * scale_step;
            }
          } else {
            double scale_step = (model_code == 2) ? (sl * post_iq) : post_iq;
            if (has_weight) scale_step *= w;
            double* msr = marginal_step_residual.data() + static_cast<size_t>(i) * n_steps;
            for (int s = 0; s < n_steps; ++s) {
              const double ind = (y >= (s + 1)) ? 1.0 : 0.0;
              msr[s] += (pgeq_buf[static_cast<size_t>(s)] - ind) * scale_step;
            }
          }
        }

        if (model_code == 2) {
          const double obs_lin = linear_buf[static_cast<size_t>(y)];
          double acc = 0.0;
          for (int k = 0; k < k_cat; ++k) {
            acc += prob_buf[static_cast<size_t>(k)] * (obs_lin - linear_buf[static_cast<size_t>(k)]);
          }
          double sc = sl * acc * post_iq;
          if (has_weight) sc *= w;
          marginal_slope_score[static_cast<size_t>(i)] += sc;
        }
      }
    }

    if (model_code == 0) {
      for (int s = 0; s < n_steps; ++s) {
        grad_step_mat(0, s) += rsm_step_q[static_cast<size_t>(s)];
      }
    }

    if (do_pop) {
      double sum_sigma_q = 0.0;
      const double z_q = quad_std_nodes[q];
      for (int r = 0; r < n_active_persons; ++r) {
        const double w_pq = posterior_flat[static_cast<size_t>(r) * n_nodes + q] *
                            person_score_q[static_cast<size_t>(r)];
        const int p_id0 = person_ids_vec[static_cast<size_t>(r)] - 1;
        for (int b = 0; b < n_beta; ++b) {
          grad_beta[b] += pop_design_person(p_id0, b) * w_pq;
        }
        sum_sigma_q += w_pq * z_q;
      }
      grad_log_sigma2 += sum_sigma_q * (0.5 * sigma_pop);
    }
  }

  if (model_code != 0 && n_steps > 0) {
    for (int i = 0; i < n; ++i) {
      const int c = step_idx_ptr[i] - 1;
      const double* msr = marginal_step_residual.data() + static_cast<size_t>(i) * n_steps;
      for (int s = 0; s < n_steps; ++s) {
        grad_step_mat(c, s) += msr[s];
      }
    }
  }
  if (model_code == 2) {
    for (int i = 0; i < n; ++i) {
      const int sl_i = slope_idx_ptr[i] - 1;
      grad_log_slope_exp[sl_i] += marginal_slope_score[static_cast<size_t>(i)];
    }
  }

  return writable::list({
    "value"_nm = neg_loglik,
    "log_marginal"_nm = log_marginal_out,
    "person_ids"_nm = person_ids_out,
    "marginal_residual"_nm = marginal_residual,
    "grad_step_mat"_nm = grad_step_mat,
    "grad_log_slope_exp"_nm = grad_log_slope_exp,
    "grad_beta"_nm = grad_beta,
    "grad_log_sigma2"_nm = grad_log_sigma2
  });
}

[[cpp11::register]]
writable::doubles mfrm_cpp_indexed_sum(doubles x,
                                       integers idx,
                                       int n_levels,
                                       double scale) {
  const int n = x.size();
  writable::doubles out(n_levels);
  for (int l = 0; l < n_levels; ++l) out[l] = 0.0;
  if (scale == 1.0) {
    for (int i = 0; i < n; ++i) {
      out[idx[i] - 1] += x[i];
    }
  } else if (scale == -1.0) {
    for (int i = 0; i < n; ++i) {
      out[idx[i] - 1] -= x[i];
    }
  } else {
    for (int i = 0; i < n; ++i) {
      out[idx[i] - 1] += scale * x[i];
    }
  }
  return out;
}

[[cpp11::register]]
writable::list mfrm_cpp_jml_eval_bundle(int model_code,
                                        integers score_k,
                                        integers person_int,
                                        int n_theta,
                                        doubles eta,
                                        doubles step_cum_rsm,
                                        doubles_matrix<> step_cum_mat,
                                        integers step_idx,
                                        doubles slopes,
                                        integers slope_idx,
                                        sexp weight,
                                        bool include_grad) {
  const int n = eta.size();
  const int k_cat = (model_code == 0) ? step_cum_rsm.size() : step_cum_mat.ncol();
  const int n_steps = std::max(k_cat - 1, 0);
  const int n_criteria = (model_code == 0) ? 1 : step_cum_mat.nrow();
  const int n_slopes = (model_code == 2) ? slopes.size() : 0;
  const bool has_weight = !Rf_isNull(weight);
  doubles weight_vec;
  if (has_weight) weight_vec = as_doubles(weight);

  std::vector<double> step_rm;
  if (model_code != 0) {
    step_rm.resize(static_cast<size_t>(n_criteria) * k_cat);
    for (int c = 0; c < n_criteria; ++c) {
      for (int k = 0; k < k_cat; ++k) {
        step_rm[static_cast<size_t>(c) * k_cat + k] = step_cum_mat(c, k);
      }
    }
  }
  const double* step_rsm_ptr = (model_code == 0 && k_cat > 0) ? REAL(step_cum_rsm) : nullptr;
  const int* step_idx_ptr = (model_code != 0 && n > 0) ? INTEGER(step_idx) : nullptr;
  const double* slopes_ptr = (model_code == 2 && n_slopes > 0) ? REAL(slopes) : nullptr;
  const int* slope_idx_ptr = (model_code == 2 && n > 0) ? INTEGER(slope_idx) : nullptr;

  std::vector<double> linear_buf(static_cast<size_t>(k_cat));
  std::vector<double> log_num_buf(static_cast<size_t>(k_cat));
  std::vector<double> prob_buf(static_cast<size_t>(k_cat));
  std::vector<double> pgeq_buf(static_cast<size_t>(std::max(n_steps, 1)));

  double neg_loglik = 0.0;
  writable::doubles residual(include_grad ? n : 0);
  writable::doubles grad_theta_exp(include_grad ? n_theta : 0);
  for (int t = 0; t < (include_grad ? n_theta : 0); ++t) grad_theta_exp[t] = 0.0;
  writable::doubles_matrix<> grad_step_mat(include_grad ? n_criteria : 0, include_grad ? n_steps : 0);
  if (include_grad) {
    for (int c = 0; c < n_criteria; ++c) {
      for (int s = 0; s < n_steps; ++s) {
        grad_step_mat(c, s) = 0.0;
      }
    }
  }
  writable::doubles grad_log_slope_exp(include_grad ? n_slopes : 0);
  for (int s = 0; s < (include_grad ? n_slopes : 0); ++s) grad_log_slope_exp[s] = 0.0;

  for (int i = 0; i < n; ++i) {
    const int y = score_k[i];
    const double eta_i = eta[i];
    const double sl = (model_code == 2) ? slopes_ptr[slope_idx_ptr[i] - 1] : 1.0;
    const double w = has_weight ? weight_vec[i] : 1.0;
    const double* st_ptr = (model_code == 0)
      ? step_rsm_ptr
      : (step_rm.data() + static_cast<size_t>(step_idx_ptr[i] - 1) * k_cat);

    double row_max = -std::numeric_limits<double>::infinity();
    for (int k = 0; k < k_cat; ++k) {
      const double lin = eta_i * static_cast<double>(k) - st_ptr[k];
      const double val = (model_code == 2) ? (lin * sl) : lin;
      linear_buf[static_cast<size_t>(k)] = lin;
      log_num_buf[static_cast<size_t>(k)] = val;
      if (val > row_max) row_max = val;
    }
    double denom_sum = 0.0;
    for (int k = 0; k < k_cat; ++k) {
      const double e = std::exp(log_num_buf[static_cast<size_t>(k)] - row_max);
      prob_buf[static_cast<size_t>(k)] = e;
      denom_sum += e;
    }
    double lp = (log_num_buf[static_cast<size_t>(y)] - row_max) - std::log(denom_sum);
    if (has_weight) lp *= w;
    neg_loglik -= lp;

    if (include_grad) {
      double expected = 0.0;
      for (int k = 0; k < k_cat; ++k) {
        prob_buf[static_cast<size_t>(k)] /= denom_sum;
        expected += prob_buf[static_cast<size_t>(k)] * static_cast<double>(k);
      }
      double res_i = (model_code == 2)
        ? (sl * (static_cast<double>(y) - expected))
        : (static_cast<double>(y) - expected);
      if (has_weight) res_i *= w;
      residual[i] = res_i;
      if (n_theta > 0) {
        grad_theta_exp[person_int[i] - 1] += res_i;
      }

      if (n_steps > 0) {
        pgeq_buf[static_cast<size_t>(n_steps - 1)] = prob_buf[static_cast<size_t>(k_cat - 1)];
        for (int s = n_steps - 2; s >= 0; --s) {
          pgeq_buf[static_cast<size_t>(s)] = pgeq_buf[static_cast<size_t>(s + 1)] +
                                             prob_buf[static_cast<size_t>(s + 1)];
        }
        const int c_row = (model_code == 0) ? 0 : (step_idx_ptr[i] - 1);
        double scale_step = (model_code == 2) ? sl : 1.0;
        if (has_weight) scale_step *= w;
        for (int s = 0; s < n_steps; ++s) {
          const double ind = (y >= (s + 1)) ? 1.0 : 0.0;
          grad_step_mat(c_row, s) += (pgeq_buf[static_cast<size_t>(s)] - ind) * scale_step;
        }
      }

      if (model_code == 2) {
        const double obs_lin = linear_buf[static_cast<size_t>(y)];
        double acc = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          acc += prob_buf[static_cast<size_t>(k)] * (obs_lin - linear_buf[static_cast<size_t>(k)]);
        }
        double sc = sl * acc;
        if (has_weight) sc *= w;
        grad_log_slope_exp[slope_idx_ptr[i] - 1] += sc;
      }
    }
  }

  if (!include_grad) {
    return writable::list({"value"_nm = neg_loglik});
  }
  return writable::list({
    "value"_nm = neg_loglik,
    "residual"_nm = residual,
    "grad_theta_exp"_nm = grad_theta_exp,
    "grad_step_mat"_nm = grad_step_mat,
    "grad_log_slope_exp"_nm = grad_log_slope_exp
  });
}

[[cpp11::register]]
writable::list mfrm_cpp_adaptive_quad_basis(int model_code,
                                            int k_cat,
                                            integers persons,
                                            integers group_starts,
                                            integers group_rows,
                                            integers score_k,
                                            doubles weights,
                                            doubles base_eta,
                                            doubles step_cum_rsm,
                                            doubles_matrix<> step_cum_mat,
                                            integers step_idx,
                                            doubles slopes,
                                            integers slope_idx,
                                            doubles mu,
                                            double sigma,
                                            doubles quad_nodes,
                                            doubles log_quad_weights,
                                            doubles log_dnorm_quad_nodes,
                                            doubles_matrix<> init_nodes,
                                            doubles_matrix<> init_log_weights) {
  const int n_groups = persons.size();
  const int n_nodes = quad_nodes.size();
  const int n_person_total = init_nodes.nrow();
  const int n_criteria = (model_code == 0) ? 1 : step_cum_mat.nrow();

  writable::doubles_matrix<> out_nodes(n_person_total, n_nodes);
  writable::doubles_matrix<> out_log_weights(n_person_total, n_nodes);
  for (int p = 0; p < n_person_total; ++p) {
    for (int q = 0; q < n_nodes; ++q) {
      out_nodes(p, q) = init_nodes(p, q);
      out_log_weights(p, q) = init_log_weights(p, q);
    }
  }

  std::vector<double> step_rm;
  if (model_code != 0) {
    step_rm.resize(static_cast<size_t>(n_criteria) * k_cat);
    for (int c = 0; c < n_criteria; ++c) {
      for (int k = 0; k < k_cat; ++k) {
        step_rm[static_cast<size_t>(c) * k_cat + k] = step_cum_mat(c, k);
      }
    }
  }
  const double* step_rsm_ptr = (model_code == 0 && k_cat > 0) ? REAL(step_cum_rsm) : nullptr;

  std::vector<double> base_sub;
  std::vector<int> score_sub;
  std::vector<double> weight_sub;
  std::vector<double> slope_sub;
  std::vector<int> step_sub;
  std::vector<double> log_num_buf(static_cast<size_t>(k_cat));
  std::vector<double> prob_buf(static_cast<size_t>(k_cat));

  for (int g = 0; g < n_groups; ++g) {
    const int p_id0 = persons[g] - 1;
    const int r_begin = group_starts[g];
    const int r_end = group_starts[g + 1];
    const int m_obs = r_end - r_begin;

    base_sub.resize(static_cast<size_t>(m_obs));
    score_sub.resize(static_cast<size_t>(m_obs));
    weight_sub.resize(static_cast<size_t>(m_obs));
    if (model_code != 0) step_sub.resize(static_cast<size_t>(m_obs));
    if (model_code == 2) slope_sub.resize(static_cast<size_t>(m_obs));

    for (int j = 0; j < m_obs; ++j) {
      const int obs = group_rows[r_begin + j] - 1;
      base_sub[static_cast<size_t>(j)] = base_eta[obs];
      score_sub[static_cast<size_t>(j)] = score_k[obs];
      weight_sub[static_cast<size_t>(j)] = weights[obs];
      if (model_code != 0) step_sub[static_cast<size_t>(j)] = step_idx[obs];
      if (model_code == 2) {
        const double sl = slopes[slope_idx[obs] - 1];
        if (!std::isfinite(sl) || sl <= 0.0) {
          return writable::list({
            "error_code"_nm = 3,
            "error_message"_nm = std::string("The retained slopes are not positive finite values.")
          });
        }
        slope_sub[static_cast<size_t>(j)] = sl;
      }
    }

    const AdaptiveModeResult loc = find_adaptive_posterior_mode(
      mu[p_id0], sigma, m_obs, k_cat, model_code,
      base_sub.data(), score_sub.data(), weight_sub.data(),
      model_code == 2 ? slope_sub.data() : nullptr,
      step_rsm_ptr, step_rm.data(),
      model_code != 0 ? step_sub.data() : nullptr,
      log_num_buf, prob_buf);

    if (loc.error_code != 0) {
      return writable::list({
        "error_code"_nm = loc.error_code,
        "error_message"_nm = loc.error_message
      });
    }

    const double log_scale = std::log(loc.scale);
    for (int q = 0; q < n_nodes; ++q) {
      const double z_q = loc.mode + loc.scale * quad_nodes[q];
      out_nodes(p_id0, q) = mu[p_id0] + sigma * z_q;
      out_log_weights(p_id0, q) = log_quad_weights[q] +
                                  Rf_dnorm4(z_q, 0.0, 1.0, 1) -
                                  log_dnorm_quad_nodes[q] +
                                  log_scale;
    }
  }

  return writable::list({
    "error_code"_nm = 0,
    "error_message"_nm = std::string(""),
    "nodes"_nm = out_nodes,
    "log_weights"_nm = out_log_weights
  });
}

[[cpp11::register]]
writable::list mfrm_cpp_adaptive_mml_eval(int model_code,
                                          int k_cat,
                                          int p_dim,
                                          integers pat_person_ids,
                                          doubles pat_multiplicity,
                                          integers pat_row_starts,
                                          integers pat_obs_idx,
                                          integers score_k,
                                          doubles weights,
                                          doubles base_eta,
                                          doubles step_cum_rsm,
                                          doubles_matrix<> step_cum_mat,
                                          integers step_idx,
                                          doubles slopes,
                                          integers slope_idx,
                                          doubles mu,
                                          double sigma,
                                          doubles rule_nodes,
                                          doubles log_rule_weights,
                                          doubles log_dnorm_rule_nodes,
                                          integers cum_pat_col_starts,
                                          integers cum_active_col,
                                          integers cum_col_val_starts,
                                          integers cum_row_i,
                                          doubles cum_val_x,
                                          bool has_log_slopes,
                                          integers ls_pat_col_starts,
                                          integers ls_active_col,
                                          integers ls_col_val_starts,
                                          integers ls_row_i,
                                          doubles ls_val_x,
                                          integers beta_cols0,
                                          doubles_matrix<> pat_pop_design,
                                          int log_sigma2_col0,
                                          bool include_person_summaries) {
  const int n_pat = pat_person_ids.size();
  const int n_nodes = rule_nodes.size();
  const int n_steps = k_cat - 1;
  const int n_criteria = (model_code == 0) ? 1 : step_cum_mat.nrow();
  const int n_beta = beta_cols0.size();

  std::vector<double> step_rm;
  if (model_code != 0) {
    step_rm.resize(static_cast<size_t>(n_criteria) * k_cat);
    for (int c = 0; c < n_criteria; ++c) {
      for (int k = 0; k < k_cat; ++k) {
        step_rm[static_cast<size_t>(c) * k_cat + k] = step_cum_mat(c, k);
      }
    }
  }
  const double* step_rsm_ptr = (model_code == 0 && k_cat > 0) ? REAL(step_cum_rsm) : nullptr;
  const int* cum_pcs = INTEGER(cum_pat_col_starts);
  const int* cum_ac = cum_active_col.size() > 0 ? INTEGER(cum_active_col) : nullptr;
  const int* cum_cvs = INTEGER(cum_col_val_starts);
  const int* cum_ri = cum_row_i.size() > 0 ? INTEGER(cum_row_i) : nullptr;
  const double* cum_vx = cum_val_x.size() > 0 ? REAL(cum_val_x) : nullptr;

  const int* ls_pcs = INTEGER(ls_pat_col_starts);
  const int* ls_ac = ls_active_col.size() > 0 ? INTEGER(ls_active_col) : nullptr;
  const int* ls_cvs = INTEGER(ls_col_val_starts);
  const int* ls_ri = ls_row_i.size() > 0 ? INTEGER(ls_row_i) : nullptr;
  const double* ls_vx = ls_val_x.size() > 0 ? REAL(ls_val_x) : nullptr;

  double value = 0.0;
  std::vector<double> gradient(static_cast<size_t>(p_dim), 0.0);
  writable::doubles_matrix<> pat_summaries(include_person_summaries ? n_pat : 0,
                                           include_person_summaries ? 3 : 0);

  std::vector<double> base_sub;
  std::vector<int> score_sub;
  std::vector<double> weight_sub;
  std::vector<double> slope_sub;
  std::vector<int> step_sub;

  std::vector<double> log_num_buf(static_cast<size_t>(k_cat));
  std::vector<double> prob_buf(static_cast<size_t>(k_cat));
  std::vector<double> mode_prob;
  std::vector<double> mode_linear;
  std::vector<double> mode_expected;
  std::vector<double> mode_variance;

  std::vector<double> vec_cum;
  std::vector<double> vec_ls;
  std::vector<double> mixed(static_cast<size_t>(p_dim));
  std::vector<double> curv_deriv(static_cast<size_t>(p_dim));
  std::vector<double> mode_deriv(static_cast<size_t>(p_dim));
  std::vector<double> log_scale_deriv(static_cast<size_t>(p_dim));
  std::vector<double> fixed_deriv(static_cast<size_t>(p_dim));

  std::vector<double> z_nodes(static_cast<size_t>(n_nodes));
  std::vector<double> node_ll(static_cast<size_t>(n_nodes));
  std::vector<double> node_score(static_cast<size_t>(n_nodes));
  std::vector<double> joint(static_cast<size_t>(n_nodes));
  std::vector<double> posterior(static_cast<size_t>(n_nodes));
  std::vector<double> node_probs_all;
  std::vector<double> node_linear_all;
  std::vector<double> residual_mean;
  std::vector<double> slope_mean;

  const double inv_sigma2 = 1.0 / (sigma * sigma);

  for (int pat = 0; pat < n_pat; ++pat) {
    const int p_id0 = pat_person_ids[pat] - 1;
    const double mult = pat_multiplicity[pat];
    const double mu_p = mu[p_id0];
    const int r_begin = pat_row_starts[pat];
    const int r_end = pat_row_starts[pat + 1];
    const int m_obs = r_end - r_begin;

    base_sub.resize(static_cast<size_t>(m_obs));
    score_sub.resize(static_cast<size_t>(m_obs));
    weight_sub.resize(static_cast<size_t>(m_obs));
    slope_sub.resize(static_cast<size_t>(m_obs));
    if (model_code != 0) step_sub.resize(static_cast<size_t>(m_obs));

    for (int r = 0; r < m_obs; ++r) {
      const int obs = pat_obs_idx[r_begin + r] - 1;
      base_sub[static_cast<size_t>(r)] = base_eta[obs];
      score_sub[static_cast<size_t>(r)] = score_k[obs];
      weight_sub[static_cast<size_t>(r)] = weights[obs];
      if (model_code != 0) step_sub[static_cast<size_t>(r)] = step_idx[obs];
      if (model_code == 2) {
        const double sl = slopes[slope_idx[obs] - 1];
        if (!std::isfinite(sl) || sl <= 0.0) {
          return writable::list({
            "error_code"_nm = 3,
            "error_message"_nm = std::string("The retained slopes are not positive finite values.")
          });
        }
        slope_sub[static_cast<size_t>(r)] = sl;
      } else {
        slope_sub[static_cast<size_t>(r)] = 1.0;
      }
    }

    const AdaptiveModeResult loc = find_adaptive_posterior_mode(
      mu_p, sigma, m_obs, k_cat, model_code,
      base_sub.data(), score_sub.data(), weight_sub.data(),
      model_code == 2 ? slope_sub.data() : nullptr,
      step_rsm_ptr, step_rm.data(),
      model_code != 0 ? step_sub.data() : nullptr,
      log_num_buf, prob_buf);

    if (loc.error_code != 0) {
      return writable::list({
        "error_code"_nm = loc.error_code,
        "error_message"_nm = loc.error_message
      });
    }

    const double mode = loc.mode;
    const double scale = loc.scale;

    // Evaluate mode probabilities, linear parts, moments, and third moment.
    mode_prob.resize(static_cast<size_t>(m_obs) * k_cat);
    mode_linear.resize(static_cast<size_t>(m_obs) * k_cat);
    mode_expected.resize(static_cast<size_t>(m_obs));
    mode_variance.resize(static_cast<size_t>(m_obs));
    const double theta_mode = mu_p + sigma * mode;
    double mode_info_sum = 0.0;
    double curvature_theta = 0.0;

    for (int r = 0; r < m_obs; ++r) {
      const double eta = base_sub[static_cast<size_t>(r)] + theta_mode;
      const double sl = slope_sub[static_cast<size_t>(r)];
      const double w = weight_sub[static_cast<size_t>(r)];
      const double* st_ptr = (model_code == 0)
        ? step_rsm_ptr
        : (step_rm.data() + static_cast<size_t>(step_sub[static_cast<size_t>(r)] - 1) * k_cat);

      double row_max = -std::numeric_limits<double>::infinity();
      for (int k = 0; k < k_cat; ++k) {
        const double lin = eta * static_cast<double>(k) - st_ptr[k];
        const double val = (model_code == 2) ? (lin * sl) : lin;
        mode_linear[static_cast<size_t>(r) * k_cat + k] = lin;
        log_num_buf[static_cast<size_t>(k)] = val;
        if (val > row_max) row_max = val;
      }
      double denom = 0.0;
      for (int k = 0; k < k_cat; ++k) {
        const double e = std::exp(log_num_buf[static_cast<size_t>(k)] - row_max);
        prob_buf[static_cast<size_t>(k)] = e;
        denom += e;
      }
      double exp_r = 0.0;
      for (int k = 0; k < k_cat; ++k) {
        const double p_rk = prob_buf[static_cast<size_t>(k)] / denom;
        mode_prob[static_cast<size_t>(r) * k_cat + k] = p_rk;
        exp_r += p_rk * static_cast<double>(k);
      }
      mode_expected[static_cast<size_t>(r)] = exp_r;
      double var_r = 0.0;
      double m3_r = 0.0;
      for (int k = 0; k < k_cat; ++k) {
        const double c_rk = static_cast<double>(k) - exp_r;
        const double c2 = c_rk * c_rk;
        const double p_rk = mode_prob[static_cast<size_t>(r) * k_cat + k];
        var_r += p_rk * c2;
        m3_r += p_rk * c2 * c_rk;
      }
      mode_variance[static_cast<size_t>(r)] = var_r;
      mode_info_sum += (w * sl * sl) * var_r;
      curvature_theta += (w * sl * sl * sl) * m3_r;
    }
    const double at_mode_information = sigma * sigma * mode_info_sum + 1.0;
    const double curvature = at_mode_information * inv_sigma2;

    vec_cum.resize(static_cast<size_t>(m_obs) * std::max(n_steps, 0));
    vec_ls.resize(static_cast<size_t>(m_obs));

    // mixed = -project(prob * centered * (w * slope), linear) +
    //         crossprod(log_slope, w * slope * (score_k - expected)) +
    //         mean_derivative / sigma^2 + (mode / sigma) * variance_derivative
    std::fill(mixed.begin(), mixed.end(), 0.0);
    for (int r = 0; r < m_obs; ++r) {
      const double ws = weight_sub[static_cast<size_t>(r)] * slope_sub[static_cast<size_t>(r)];
      const double sl = slope_sub[static_cast<size_t>(r)];
      const double exp_r = mode_expected[static_cast<size_t>(r)];
      for (int k = 1; k < k_cat; ++k) {
        const double coeff = mode_prob[static_cast<size_t>(r) * k_cat + k] *
                             (static_cast<double>(k) - exp_r) * ws;
        vec_cum[static_cast<size_t>(k - 1) * m_obs + r] = coeff * sl;
      }
    }
    if (n_steps > 0 && cum_ac != nullptr) {
      sparse_crossprod_add(pat, cum_pcs, cum_ac, cum_cvs, cum_ri, cum_vx,
                           vec_cum.data(), -1.0, mixed);
    }
    if (has_log_slopes && ls_ac != nullptr) {
      for (int r = 0; r < m_obs; ++r) {
        const double ws = weight_sub[static_cast<size_t>(r)] * slope_sub[static_cast<size_t>(r)];
        const double sl = slope_sub[static_cast<size_t>(r)];
        const double exp_r = mode_expected[static_cast<size_t>(r)];
        double dot_lin = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          const double coeff = mode_prob[static_cast<size_t>(r) * k_cat + k] *
                               (static_cast<double>(k) - exp_r) * ws;
          dot_lin += coeff * mode_linear[static_cast<size_t>(r) * k_cat + k];
        }
        // Combine -project's log_slope term with +crossprod(log_slope, w * slope * (y - expected))
        vec_ls[static_cast<size_t>(r)] =
          -dot_lin * sl + ws * (static_cast<double>(score_sub[static_cast<size_t>(r)]) - exp_r);
      }
      sparse_crossprod_add(pat, ls_pcs, ls_ac, ls_cvs, ls_ri, ls_vx,
                           vec_ls.data(), 1.0, mixed);
    }
    for (int b = 0; b < n_beta; ++b) {
      mixed[static_cast<size_t>(beta_cols0[b])] += pat_pop_design(pat, b) * inv_sigma2;
    }
    if (log_sigma2_col0 >= 0) {
      mixed[static_cast<size_t>(log_sigma2_col0)] += mode / sigma;
    }

    // curvature_derivative = project(prob * (centered^2 - variance) * (w * slope^2), linear) +
    //                        crossprod(log_slope, 2 * w * slope^2 * variance) -
    //                        variance_derivative / sigma^2
    std::fill(curv_deriv.begin(), curv_deriv.end(), 0.0);
    for (int r = 0; r < m_obs; ++r) {
      const double sl = slope_sub[static_cast<size_t>(r)];
      const double ws2 = weight_sub[static_cast<size_t>(r)] * sl * sl;
      const double exp_r = mode_expected[static_cast<size_t>(r)];
      const double var_r = mode_variance[static_cast<size_t>(r)];
      for (int k = 1; k < k_cat; ++k) {
        const double c_rk = static_cast<double>(k) - exp_r;
        const double coeff = mode_prob[static_cast<size_t>(r) * k_cat + k] *
                             (c_rk * c_rk - var_r) * ws2;
        vec_cum[static_cast<size_t>(k - 1) * m_obs + r] = coeff * sl;
      }
    }
    if (n_steps > 0 && cum_ac != nullptr) {
      sparse_crossprod_add(pat, cum_pcs, cum_ac, cum_cvs, cum_ri, cum_vx,
                           vec_cum.data(), 1.0, curv_deriv);
    }
    if (has_log_slopes && ls_ac != nullptr) {
      for (int r = 0; r < m_obs; ++r) {
        const double sl = slope_sub[static_cast<size_t>(r)];
        const double ws2 = weight_sub[static_cast<size_t>(r)] * sl * sl;
        const double exp_r = mode_expected[static_cast<size_t>(r)];
        const double var_r = mode_variance[static_cast<size_t>(r)];
        double dot_lin = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          const double c_rk = static_cast<double>(k) - exp_r;
          const double coeff = mode_prob[static_cast<size_t>(r) * k_cat + k] *
                               (c_rk * c_rk - var_r) * ws2;
          dot_lin += coeff * mode_linear[static_cast<size_t>(r) * k_cat + k];
        }
        vec_ls[static_cast<size_t>(r)] = dot_lin * sl + 2.0 * ws2 * var_r;
      }
      sparse_crossprod_add(pat, ls_pcs, ls_ac, ls_cvs, ls_ri, ls_vx,
                           vec_ls.data(), 1.0, curv_deriv);
    }
    if (log_sigma2_col0 >= 0) {
      curv_deriv[static_cast<size_t>(log_sigma2_col0)] -= inv_sigma2;
    }

    const double inv_curv = 1.0 / curvature;
    const double inv_2curv = 0.5 * inv_curv;
    for (int d = 0; d < p_dim; ++d) {
      const double md = mixed[static_cast<size_t>(d)] * inv_curv;
      mode_deriv[static_cast<size_t>(d)] = md;
      log_scale_deriv[static_cast<size_t>(d)] =
        -(curv_deriv[static_cast<size_t>(d)] + curvature_theta * md) * inv_2curv;
    }

    // Evaluate adaptive Gauss-Hermite nodes z = mode + scale * rule_nodes
    node_probs_all.resize(static_cast<size_t>(n_nodes) * m_obs * k_cat);
    if (has_log_slopes) {
      node_linear_all.resize(static_cast<size_t>(n_nodes) * m_obs * k_cat);
    }
    double offset = -std::numeric_limits<double>::infinity();

    for (int q = 0; q < n_nodes; ++q) {
      const double z_q = mode + scale * rule_nodes[q];
      z_nodes[static_cast<size_t>(q)] = z_q;
      const double theta_q = mu_p + sigma * z_q;
      double ll_q = 0.0;
      double score_sum_q = 0.0;
      double info_sum_q = 0.0;
      double* probs_q_ptr = node_probs_all.data() + static_cast<size_t>(q) * m_obs * k_cat;
      double* lin_q_ptr = has_log_slopes
        ? (node_linear_all.data() + static_cast<size_t>(q) * m_obs * k_cat)
        : nullptr;

      for (int r = 0; r < m_obs; ++r) {
        const double eta = base_sub[static_cast<size_t>(r)] + theta_q;
        const double sl = slope_sub[static_cast<size_t>(r)];
        const double w = weight_sub[static_cast<size_t>(r)];
        const int y = score_sub[static_cast<size_t>(r)];
        const double* st_ptr = (model_code == 0)
          ? step_rsm_ptr
          : (step_rm.data() + static_cast<size_t>(step_sub[static_cast<size_t>(r)] - 1) * k_cat);

        double row_max = -std::numeric_limits<double>::infinity();
        for (int k = 0; k < k_cat; ++k) {
          const double lin = eta * static_cast<double>(k) - st_ptr[k];
          const double val = (model_code == 2) ? (lin * sl) : lin;
          if (has_log_slopes) lin_q_ptr[static_cast<size_t>(r) * k_cat + k] = lin;
          log_num_buf[static_cast<size_t>(k)] = val;
          if (val > row_max) row_max = val;
        }
        double denom = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          const double e = std::exp(log_num_buf[static_cast<size_t>(k)] - row_max);
          prob_buf[static_cast<size_t>(k)] = e;
          denom += e;
        }
        ll_q += w * ((log_num_buf[static_cast<size_t>(y)] - row_max) - std::log(denom));
        double exp_r = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          const double p_rk = prob_buf[static_cast<size_t>(k)] / denom;
          probs_q_ptr[static_cast<size_t>(r) * k_cat + k] = p_rk;
          exp_r += p_rk * static_cast<double>(k);
        }
        double var_r = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          const double c_rk = static_cast<double>(k) - exp_r;
          var_r += probs_q_ptr[static_cast<size_t>(r) * k_cat + k] * c_rk * c_rk;
        }
        score_sum_q += (w * sl) * (static_cast<double>(y) - exp_r);
        info_sum_q += (w * sl * sl) * var_r;
      }

      const double m1_q = sigma * score_sum_q;
      const double m2_q = sigma * sigma * info_sum_q;
      if (!std::isfinite(m1_q) || !std::isfinite(m2_q)) {
        return writable::list({
          "error_code"_nm = 1,
          "error_message"_nm = std::string("Non-finite adaptive posterior score or information.")
        });
      }
      node_ll[static_cast<size_t>(q)] = ll_q;
      node_score[static_cast<size_t>(q)] = m1_q - z_q;
      if (ll_q > offset) offset = ll_q;
    }

    const double log_scale = std::log(scale);
    double max_joint = -std::numeric_limits<double>::infinity();
    for (int q = 0; q < n_nodes; ++q) {
      const double z_q = z_nodes[static_cast<size_t>(q)];
      const double j_q = log_rule_weights[q] +
                         (node_ll[static_cast<size_t>(q)] - offset) +
                         Rf_dnorm4(z_q, 0.0, 1.0, 1) -
                         log_dnorm_rule_nodes[q] +
                         log_scale;
      if (!std::isfinite(j_q)) {
        return writable::list({
          "error_code"_nm = 4,
          "error_message"_nm = std::string("Non-finite integration terms.")
        });
      }
      joint[static_cast<size_t>(q)] = j_q;
      if (j_q > max_joint) max_joint = j_q;
    }

    double sum_post = 0.0;
    for (int q = 0; q < n_nodes; ++q) {
      const double p_q = std::exp(joint[static_cast<size_t>(q)] - max_joint);
      posterior[static_cast<size_t>(q)] = p_q;
      sum_post += p_q;
    }
    const double normalizer = max_joint + std::log(sum_post);
    for (int q = 0; q < n_nodes; ++q) {
      posterior[static_cast<size_t>(q)] /= sum_post;
    }

    value -= mult * (offset + normalizer);

    if (include_person_summaries) {
      double center = 0.0;
      for (int q = 0; q < n_nodes; ++q) {
        center += posterior[static_cast<size_t>(q)] * z_nodes[static_cast<size_t>(q)];
      }
      double z_var = 0.0;
      for (int q = 0; q < n_nodes; ++q) {
        const double dz = z_nodes[static_cast<size_t>(q)] - center;
        z_var += posterior[static_cast<size_t>(q)] * dz * dz;
      }
      pat_summaries(pat, 0) = offset + normalizer;
      pat_summaries(pat, 1) = mu_p + sigma * center;
      pat_summaries(pat, 2) = sigma * std::sqrt(z_var);
    }

    residual_mean.assign(static_cast<size_t>(m_obs) * k_cat, 0.0);
    if (has_log_slopes) {
      slope_mean.assign(static_cast<size_t>(m_obs), 0.0);
    }
    for (int q = 0; q < n_nodes; ++q) {
      const double post_q = posterior[static_cast<size_t>(q)];
      const double* probs_q_ptr = node_probs_all.data() + static_cast<size_t>(q) * m_obs * k_cat;
      const double* lin_q_ptr = has_log_slopes
        ? (node_linear_all.data() + static_cast<size_t>(q) * m_obs * k_cat)
        : nullptr;
      for (int r = 0; r < m_obs; ++r) {
        const int y = score_sub[static_cast<size_t>(r)];
        double sl_acc = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          const double res_rk = ((k == y) ? 1.0 : 0.0) - probs_q_ptr[static_cast<size_t>(r) * k_cat + k];
          residual_mean[static_cast<size_t>(r) * k_cat + k] += post_q * res_rk;
          if (has_log_slopes) {
            sl_acc += res_rk * lin_q_ptr[static_cast<size_t>(r) * k_cat + k];
          }
        }
        if (has_log_slopes) {
          slope_mean[static_cast<size_t>(r)] += post_q * sl_acc;
        }
      }
    }

    std::fill(fixed_deriv.begin(), fixed_deriv.end(), 0.0);
    for (int r = 0; r < m_obs; ++r) {
      const double ws = weight_sub[static_cast<size_t>(r)] * slope_sub[static_cast<size_t>(r)];
      for (int k = 1; k < k_cat; ++k) {
        vec_cum[static_cast<size_t>(k - 1) * m_obs + r] =
          residual_mean[static_cast<size_t>(r) * k_cat + k] * ws;
      }
    }
    if (n_steps > 0 && cum_ac != nullptr) {
      sparse_crossprod_add(pat, cum_pcs, cum_ac, cum_cvs, cum_ri, cum_vx,
                           vec_cum.data(), 1.0, fixed_deriv);
    }
    if (has_log_slopes && ls_ac != nullptr) {
      for (int r = 0; r < m_obs; ++r) {
        vec_ls[static_cast<size_t>(r)] =
          slope_mean[static_cast<size_t>(r)] *
          (weight_sub[static_cast<size_t>(r)] * slope_sub[static_cast<size_t>(r)]);
      }
      sparse_crossprod_add(pat, ls_pcs, ls_ac, ls_cvs, ls_ri, ls_vx,
                           vec_ls.data(), 1.0, fixed_deriv);
    }

    double z_post_mean = 0.0;
    double z2m1_post_mean = 0.0;
    double score_post_mean = 0.0;
    double scale_post_mean = 0.0;
    for (int q = 0; q < n_nodes; ++q) {
      const double post_q = posterior[static_cast<size_t>(q)];
      const double z_q = z_nodes[static_cast<size_t>(q)];
      const double sc_q = node_score[static_cast<size_t>(q)];
      z_post_mean += post_q * z_q;
      z2m1_post_mean += post_q * (z_q * z_q - 1.0);
      score_post_mean += post_q * sc_q;
      scale_post_mean += post_q * (scale * sc_q * rule_nodes[q] + 1.0);
    }
    for (int b = 0; b < n_beta; ++b) {
      fixed_deriv[static_cast<size_t>(beta_cols0[b])] +=
        (z_post_mean / sigma) * pat_pop_design(pat, b);
    }
    if (log_sigma2_col0 >= 0) {
      fixed_deriv[static_cast<size_t>(log_sigma2_col0)] += 0.5 * z2m1_post_mean;
    }

    const double coef_mode = score_post_mean / sigma;
    const double coef_scale = scale_post_mean;
    for (int d = 0; d < p_dim; ++d) {
      gradient[static_cast<size_t>(d)] -= mult * (
        fixed_deriv[static_cast<size_t>(d)] +
        coef_mode * mode_deriv[static_cast<size_t>(d)] +
        coef_scale * log_scale_deriv[static_cast<size_t>(d)]
      );
    }
  }

  writable::doubles grad_out(p_dim);
  for (int d = 0; d < p_dim; ++d) {
    grad_out[d] = gradient[static_cast<size_t>(d)];
  }
  return writable::list({
    "error_code"_nm = 0,
    "error_message"_nm = std::string(""),
    "value"_nm = value,
    "gradient"_nm = grad_out,
    "pat_summaries"_nm = pat_summaries
  });
}

[[cpp11::register]]
writable::list mfrm_cpp_cml_esf(list features,
                                doubles par,
                                bool information) {
  const int p = par.size();
  const int n_items = features.size();

  int nz = 1;
  std::vector<double> logz(1, 0.0);
  // Row-major means: nz x p -> s * p + d
  std::vector<double> means(static_cast<size_t>(p), 0.0);
  // Slice-major symmetric covariances: s * p * p + d1 * p + d2
  std::vector<double> covs(information ? static_cast<size_t>(p) * p : 0, 0.0);

  std::vector<double> eta;
  std::vector<double> logs;
  std::vector<double> weights;
  std::vector<double> sqrt_weights;
  std::vector<double> cand;
  std::vector<double> centered;
  std::vector<double> mu(static_cast<size_t>(p));

  for (int item = 0; item < n_items; ++item) {
    doubles_matrix<> x(features[item]);
    const int k_max_item = x.nrow() - 1;
    eta.resize(static_cast<size_t>(k_max_item + 1));
    for (int k = 0; k <= k_max_item; ++k) {
      double acc = 0.0;
      for (int d = 0; d < p; ++d) {
        acc += x(k, d) * par[d];
      }
      if (!std::isfinite(acc)) {
        stop("Nonfinite CML category predictors.");
      }
      eta[static_cast<size_t>(k)] = acc;
    }

    const int next_nz = nz + k_max_item;
    std::vector<double> next_logz(static_cast<size_t>(next_nz), 0.0);
    std::vector<double> next_means(static_cast<size_t>(next_nz) * p, 0.0);
    std::vector<double> next_covs(information ? static_cast<size_t>(next_nz) * p * p : 0, 0.0);

    for (int s = 0; s < next_nz; ++s) {
      const int k_lo = std::max(0, s - nz + 1);
      const int k_hi = std::min(k_max_item, s);
      const int m_h = k_hi - k_lo + 1;

      logs.resize(static_cast<size_t>(m_h));
      weights.resize(static_cast<size_t>(m_h));
      double peak = -std::numeric_limits<double>::infinity();
      for (int h = 0; h < m_h; ++h) {
        const int k = k_lo + h;
        const int from = s - k;
        const double lv = logz[static_cast<size_t>(from)] + eta[static_cast<size_t>(k)];
        logs[static_cast<size_t>(h)] = lv;
        if (lv > peak) peak = lv;
      }
      double w_sum = 0.0;
      for (int h = 0; h < m_h; ++h) {
        const double w = std::exp(logs[static_cast<size_t>(h)] - peak);
        weights[static_cast<size_t>(h)] = w;
        w_sum += w;
      }
      next_logz[static_cast<size_t>(s)] = peak + std::log(w_sum);
      for (int h = 0; h < m_h; ++h) {
        weights[static_cast<size_t>(h)] /= w_sum;
      }

      cand.resize(static_cast<size_t>(m_h) * p);
      for (int h = 0; h < m_h; ++h) {
        const int k = k_lo + h;
        const int from = s - k;
        const double* m_from = means.data() + static_cast<size_t>(from) * p;
        double* c_row = cand.data() + static_cast<size_t>(h) * p;
        for (int d = 0; d < p; ++d) {
          c_row[d] = m_from[d] + x(k, d);
        }
      }

      double* dst_mu = next_means.data() + static_cast<size_t>(s) * p;
      for (int d = 0; d < p; ++d) {
        double acc = 0.0;
        for (int h = 0; h < m_h; ++h) {
          acc += weights[static_cast<size_t>(h)] * cand[static_cast<size_t>(h) * p + d];
        }
        mu[static_cast<size_t>(d)] = acc;
        dst_mu[d] = acc;
      }

      if (information) {
        sqrt_weights.resize(static_cast<size_t>(m_h));
        centered.resize(static_cast<size_t>(m_h) * p);
        for (int h = 0; h < m_h; ++h) {
          const double sw = std::sqrt(weights[static_cast<size_t>(h)]);
          sqrt_weights[static_cast<size_t>(h)] = sw;
          const double* c_row = cand.data() + static_cast<size_t>(h) * p;
          double* cen_row = centered.data() + static_cast<size_t>(h) * p;
          for (int d = 0; d < p; ++d) {
            cen_row[d] = (c_row[d] - mu[static_cast<size_t>(d)]) * sw;
          }
        }
        double* dst_cov = next_covs.data() + static_cast<size_t>(s) * p * p;
        for (int d1 = 0; d1 < p; ++d1) {
          for (int d2 = d1; d2 < p; ++d2) {
            double v_acc = 0.0;
            for (int h = 0; h < m_h; ++h) {
              v_acc += centered[static_cast<size_t>(h) * p + d1] *
                       centered[static_cast<size_t>(h) * p + d2];
            }
            for (int h = 0; h < m_h; ++h) {
              const int from = s - (k_lo + h);
              v_acc += weights[static_cast<size_t>(h)] *
                       covs[static_cast<size_t>(from) * p * p + static_cast<size_t>(d1) * p + d2];
            }
            dst_cov[static_cast<size_t>(d1) * p + d2] = v_acc;
            dst_cov[static_cast<size_t>(d2) * p + d1] = v_acc;
          }
        }
      }
    }

    nz = next_nz;
    logz.swap(next_logz);
    means.swap(next_means);
    if (information) covs.swap(next_covs);
  }

  writable::doubles logz_out(nz);
  for (int s = 0; s < nz; ++s) logz_out[s] = logz[static_cast<size_t>(s)];

  writable::doubles_matrix<> means_out(nz, p);
  for (int s = 0; s < nz; ++s) {
    for (int d = 0; d < p; ++d) {
      means_out(s, d) = means[static_cast<size_t>(s) * p + d];
    }
  }

  if (!information) {
    return writable::list({
      "logz"_nm = logz_out,
      "means"_nm = means_out,
      "covariance"_nm = R_NilValue
    });
  }

  writable::doubles cov_out(static_cast<R_xlen_t>(p) * p * nz);
  for (int s = 0; s < nz; ++s) {
    const double* src = covs.data() + static_cast<size_t>(s) * p * p;
    const R_xlen_t offset = static_cast<R_xlen_t>(s) * p * p;
    for (int d2 = 0; d2 < p; ++d2) {
      for (int d1 = 0; d1 < p; ++d1) {
        cov_out[offset + static_cast<R_xlen_t>(d2) * p + d1] =
          src[static_cast<size_t>(d1) * p + d2];
      }
    }
  }
  cov_out.attr("dim") = writable::integers({p, p, nz});

  return writable::list({
    "logz"_nm = logz_out,
    "means"_nm = means_out,
    "covariance"_nm = cov_out
  });
}

[[cpp11::register]]
writable::list mfrm_cpp_expected_category_bundle(list prob_list,
                                                 doubles_matrix<> obs_posterior,
                                                 bool include_p_geq) {
  const int n_nodes = prob_list.size();
  if (n_nodes == 0) {
    stop("`prob_list` must contain at least one probability matrix.");
  }

  doubles_matrix<> first_prob(prob_list[0]);
  const int n = first_prob.nrow();
  const int k_cat = first_prob.ncol();
  writable::doubles_matrix<> posterior_prob(n, k_cat);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < k_cat; ++k) {
      posterior_prob(i, k) = 0.0;
    }
  }

  for (int q = 0; q < n_nodes; ++q) {
    doubles_matrix<> probs_q(prob_list[q]);
    for (int i = 0; i < n; ++i) {
      const double post_q = obs_posterior(i, q);
      for (int k = 0; k < k_cat; ++k) {
        posterior_prob(i, k) += probs_q(i, k) * post_q;
      }
    }
  }

  writable::doubles expected_k(n);
  writable::doubles var_k(n);
  for (int i = 0; i < n; ++i) {
    double first_moment = 0.0;
    for (int k = 0; k < k_cat; ++k) {
      const double p = posterior_prob(i, k);
      first_moment += p * static_cast<double>(k);
    }
    expected_k[i] = first_moment;
    double variance = 0.0;
    for (int k = 0; k < k_cat; ++k) {
      const double centered = static_cast<double>(k) - first_moment;
      variance += posterior_prob(i, k) * centered * centered;
    }
    var_k[i] = variance;
  }

  if (include_p_geq && k_cat > 1) {
    return writable::list({
      "posterior_prob"_nm = posterior_prob,
      "expected_k"_nm = expected_k,
      "var_k"_nm = var_k,
      "p_geq"_nm = build_p_geq(posterior_prob)
    });
  }

  return writable::list({
    "posterior_prob"_nm = posterior_prob,
    "expected_k"_nm = expected_k,
    "var_k"_nm = var_k
  });
}

[[cpp11::register]]
writable::list mfrm_cpp_testlet_kernel(doubles_matrix<> design_X,
                                       int ns,
                                       doubles par,
                                       doubles theta,
                                       doubles rule_nodes,
                                       doubles rule_weights,
                                       bool gradient,
                                       integers persons,
                                       integers person_block_ptr,
                                       integers block_entry_ptr,
                                       integers entry_cell,
                                       integers entry_count) {
  const int nb = design_X.ncol();
  const int nd = design_X.nrow();
  const int npar = nb + ns + 1;
  const double v = par[nb + ns];
  const int nq = theta.size();
  const int np = persons.size();
  const int ng = (v == 0.0) ? 1 : rule_nodes.size();

  writable::doubles_matrix<> loglik_out(nq, np);
  for (int pi = 0; pi < np; ++pi) {
    for (int q = 0; q < nq; ++q) {
      loglik_out(q, pi) = 0.0;
    }
  }

  const R_xlen_t score_len = gradient
    ? (static_cast<R_xlen_t>(nq) * static_cast<R_xlen_t>(np) * static_cast<R_xlen_t>(npar))
    : 0;
  writable::doubles scores_out(score_len);
  writable::doubles_matrix<> theta_score_out(gradient ? nq : 0, gradient ? np : 0);
  writable::doubles_matrix<> theta_curv_out(gradient ? nq : 0, gradient ? np : 0);

  if (gradient) {
    double* sp = score_len > 0 ? REAL(scores_out) : nullptr;
    for (R_xlen_t idx = 0; idx < score_len; ++idx) {
      sp[idx] = 0.0;
    }
    scores_out.attr("dim") = writable::integers({nq, np, npar});
    for (int pi = 0; pi < np; ++pi) {
      for (int q = 0; q < nq; ++q) {
        theta_score_out(q, pi) = 0.0;
        theta_curv_out(q, pi) = 0.0;
      }
    }
  }

  if (nq == 0 || np == 0 || nd == 0) {
    return writable::list({
      "loglik"_nm = loglik_out,
      "score"_nm = gradient ? static_cast<SEXP>(scores_out) : R_NilValue,
      "theta_score"_nm = gradient ? static_cast<SEXP>(theta_score_out) : R_NilValue,
      "theta_curvature"_nm = gradient ? static_cast<SEXP>(theta_curv_out) : R_NilValue
    });
  }

  std::vector<double> gamma(static_cast<size_t>(ng), 0.0);
  std::vector<double> logweights(static_cast<size_t>(ng), 0.0);
  if (v > 0.0) {
    const double sqrt_v = std::sqrt(v);
    for (int m = 0; m < ng; ++m) {
      gamma[static_cast<size_t>(m)] = sqrt_v * rule_nodes[m];
      logweights[static_cast<size_t>(m)] = std::log(rule_weights[m]);
    }
  }

  const int k_cat = ns + 1;
  std::vector<double> steps_cum(static_cast<size_t>(k_cat), 0.0);
  for (int s = 1; s <= ns; ++s) {
    steps_cum[static_cast<size_t>(s)] =
      steps_cum[static_cast<size_t>(s - 1)] + par[nb + s - 1];
  }

  const int* p_blk_ptr = INTEGER(person_block_ptr);
  const int* blk_ent_ptr = INTEGER(block_entry_ptr);
  const int* ent_cell_ptr = entry_cell.size() > 0 ? INTEGER(entry_cell) : nullptr;
  const int* ent_cnt_ptr = entry_count.size() > 0 ? INTEGER(entry_count) : nullptr;

  std::vector<uint8_t> active_design(static_cast<size_t>(nd), 0);
  for (int pi = 0; pi < np; ++pi) {
    const int p0 = persons[pi] - 1;
    const int b_start = p_blk_ptr[p0];
    const int b_end = p_blk_ptr[p0 + 1];
    for (int b = b_start; b < b_end; ++b) {
      const int e_start = blk_ent_ptr[b];
      const int e_end = blk_ent_ptr[b + 1];
      for (int e = e_start; e < e_end; ++e) {
        const int d = ent_cell_ptr[e] % nd;
        active_design[static_cast<size_t>(d)] = 1;
      }
    }
  }

  const size_t n_qm = static_cast<size_t>(nq) * static_cast<size_t>(ng);
  std::vector<double> design_logp(static_cast<size_t>(nd) * k_cat * n_qm, 0.0);
  std::vector<double> design_mean;
  std::vector<double> design_var;
  std::vector<double> design_pgeq;
  if (gradient) {
    design_mean.assign(static_cast<size_t>(nd) * n_qm, 0.0);
    design_var.assign(static_cast<size_t>(nd) * n_qm, 0.0);
    design_pgeq.assign(static_cast<size_t>(nd) * std::max(ns, 1) * n_qm, 0.0);
  }

  std::vector<double> logw_buf(static_cast<size_t>(k_cat));
  std::vector<double> prob_buf(static_cast<size_t>(k_cat));

  for (int d = 0; d < nd; ++d) {
    if (!active_design[static_cast<size_t>(d)]) continue;
    double shift_d = 0.0;
    for (int b = 0; b < nb; ++b) {
      shift_d += design_X(d, b) * par[b];
    }
    for (int m = 0; m < ng; ++m) {
      const double gm = gamma[static_cast<size_t>(m)];
      for (int q = 0; q < nq; ++q) {
        const size_t iqm = static_cast<size_t>(m) * nq + static_cast<size_t>(q);
        const double eta = (theta[q] + gm) - shift_d;
        double rmax = -std::numeric_limits<double>::infinity();
        for (int k = 0; k < k_cat; ++k) {
          const double lw = eta * static_cast<double>(k) - steps_cum[static_cast<size_t>(k)];
          logw_buf[static_cast<size_t>(k)] = lw;
          if (lw > rmax) rmax = lw;
        }
        double sum_exp = 0.0;
        for (int k = 0; k < k_cat; ++k) {
          const double ev = std::exp(logw_buf[static_cast<size_t>(k)] - rmax);
          prob_buf[static_cast<size_t>(k)] = ev;
          sum_exp += ev;
        }
        const double den = std::log(sum_exp);
        for (int k = 0; k < k_cat; ++k) {
          design_logp[(static_cast<size_t>(d) * k_cat + k) * n_qm + iqm] =
            (logw_buf[static_cast<size_t>(k)] - rmax) - den;
        }
        if (gradient) {
          double mean_val = 0.0;
          for (int k = 0; k < k_cat; ++k) {
            const double pk = prob_buf[static_cast<size_t>(k)] / sum_exp;
            prob_buf[static_cast<size_t>(k)] = pk;
            mean_val += pk * static_cast<double>(k);
          }
          double var_val = 0.0;
          for (int k = 0; k < k_cat; ++k) {
            const double ck = static_cast<double>(k) - mean_val;
            var_val += prob_buf[static_cast<size_t>(k)] * ck * ck;
          }
          design_mean[static_cast<size_t>(d) * n_qm + iqm] = mean_val;
          design_var[static_cast<size_t>(d) * n_qm + iqm] = var_val;
          if (ns > 0) {
            double tail = 0.0;
            for (int s = ns - 1; s >= 0; --s) {
              tail += prob_buf[static_cast<size_t>(s + 1)];
              design_pgeq[(static_cast<size_t>(d) * ns + s) * n_qm + iqm] = tail;
            }
          }
        }
      }
    }
  }

  std::vector<double> lc(n_qm);
  std::vector<double> S_mat(gradient ? n_qm : 0);
  std::vector<double> V_mat(gradient ? n_qm : 0);
  std::vector<double> lw_mat(n_qm);
  std::vector<double> weights_mat(gradient ? n_qm : 0);
  std::vector<double> offset_vec(static_cast<size_t>(nq));
  std::vector<double> marginal_vec(static_cast<size_t>(nq));
  std::vector<double> r_bar(gradient ? static_cast<size_t>(nq) : 0);
  double* scores_ptr = (gradient && score_len > 0) ? REAL(scores_out) : nullptr;

  for (int pi = 0; pi < np; ++pi) {
    const int p0 = persons[pi] - 1;
    const int b_start = p_blk_ptr[p0];
    const int b_end = p_blk_ptr[p0 + 1];
    for (int b_idx = b_start; b_idx < b_end; ++b_idx) {
      std::fill(lc.begin(), lc.end(), 0.0);
      if (gradient) {
        std::fill(S_mat.begin(), S_mat.end(), 0.0);
        std::fill(V_mat.begin(), V_mat.end(), 0.0);
      }
      const int e_start = blk_ent_ptr[b_idx];
      const int e_end = blk_ent_ptr[b_idx + 1];
      for (int e_idx = e_start; e_idx < e_end; ++e_idx) {
        const int cell0 = ent_cell_ptr[e_idx];
        const int d = cell0 % nd;
        const int y = cell0 / nd;
        const double cnt = static_cast<double>(ent_cnt_ptr[e_idx]);
        const double* lp_src = design_logp.data() + (static_cast<size_t>(d) * k_cat + y) * n_qm;
        for (size_t iqm = 0; iqm < n_qm; ++iqm) {
          lc[iqm] += cnt * lp_src[iqm];
        }
        if (gradient) {
          const double yd = static_cast<double>(y);
          const double* mn_src = design_mean.data() + static_cast<size_t>(d) * n_qm;
          const double* vr_src = design_var.data() + static_cast<size_t>(d) * n_qm;
          for (size_t iqm = 0; iqm < n_qm; ++iqm) {
            S_mat[iqm] += cnt * (yd - mn_src[iqm]);
            V_mat[iqm] += cnt * vr_src[iqm];
          }
        }
      }

      for (int q = 0; q < nq; ++q) {
        double off = -std::numeric_limits<double>::infinity();
        for (int m = 0; m < ng; ++m) {
          const double val = lc[static_cast<size_t>(m) * nq + q];
          if (val > off) off = val;
        }
        offset_vec[static_cast<size_t>(q)] = off;
        double largest = -std::numeric_limits<double>::infinity();
        for (int m = 0; m < ng; ++m) {
          const size_t iqm = static_cast<size_t>(m) * nq + q;
          const double lw_val = (lc[iqm] - off) + logweights[static_cast<size_t>(m)];
          lw_mat[iqm] = lw_val;
          if (lw_val > largest) largest = lw_val;
        }
        double sum_lw = 0.0;
        for (int m = 0; m < ng; ++m) {
          sum_lw += std::exp(lw_mat[static_cast<size_t>(m) * nq + q] - largest);
        }
        const double marg = largest + std::log(sum_lw);
        marginal_vec[static_cast<size_t>(q)] = marg;
        loglik_out(q, pi) += off + marg;
      }

      if (gradient) {
        for (int m = 0; m < ng; ++m) {
          for (int q = 0; q < nq; ++q) {
            const size_t iqm = static_cast<size_t>(m) * nq + q;
            weights_mat[iqm] = std::exp(lw_mat[iqm] - marginal_vec[static_cast<size_t>(q)]);
          }
        }

        const R_xlen_t var_col_base =
          (static_cast<R_xlen_t>(npar - 1) * static_cast<R_xlen_t>(np) + pi) * nq;
        for (int q = 0; q < nq; ++q) {
          double b_score = 0.0;
          double sv_mean = 0.0;
          double var_score_q = 0.0;
          for (int m = 0; m < ng; ++m) {
            const size_t iqm = static_cast<size_t>(m) * nq + q;
            const double w_qm = weights_mat[iqm];
            const double s_qm = S_mat[iqm];
            const double sv_qm = s_qm * s_qm - V_mat[iqm];
            b_score += w_qm * s_qm;
            sv_mean += w_qm * sv_qm;
            if (v > 0.0) {
              var_score_q += w_qm * (s_qm * (gamma[static_cast<size_t>(m)] / (2.0 * v)));
            }
          }
          if (v == 0.0) {
            var_score_q = 0.5 * sv_mean;
          }
          theta_score_out(q, pi) += b_score;
          theta_curv_out(q, pi) += sv_mean - b_score * b_score;
          scores_ptr[var_col_base + q] += var_score_q;
        }

        for (int e_idx = e_start; e_idx < e_end; ++e_idx) {
          const int cell0 = ent_cell_ptr[e_idx];
          const int d = cell0 % nd;
          const int y = cell0 / nd;
          const double yd = static_cast<double>(y);
          const double cnt = static_cast<double>(ent_cnt_ptr[e_idx]);
          const double* mn_src = design_mean.data() + static_cast<size_t>(d) * n_qm;

          if (nb > 0) {
            for (int q = 0; q < nq; ++q) {
              double acc = 0.0;
              for (int m = 0; m < ng; ++m) {
                const size_t iqm = static_cast<size_t>(m) * nq + q;
                acc += weights_mat[iqm] * (yd - mn_src[iqm]);
              }
              r_bar[static_cast<size_t>(q)] = acc;
            }
            for (int b = 0; b < nb; ++b) {
              const double xb = -cnt * design_X(d, b);
              if (xb == 0.0) continue;
              const R_xlen_t col_base =
                (static_cast<R_xlen_t>(b) * static_cast<R_xlen_t>(np) + pi) * nq;
              for (int q = 0; q < nq; ++q) {
                scores_ptr[col_base + q] += xb * r_bar[static_cast<size_t>(q)];
              }
            }
          }

          for (int s = 0; s < ns; ++s) {
            const double ind = (y >= (s + 1)) ? 1.0 : 0.0;
            const double* pg_src = design_pgeq.data() + (static_cast<size_t>(d) * ns + s) * n_qm;
            const R_xlen_t col_base =
              (static_cast<R_xlen_t>(nb + s) * static_cast<R_xlen_t>(np) + pi) * nq;
            for (int q = 0; q < nq; ++q) {
              double pg_bar = 0.0;
              for (int m = 0; m < ng; ++m) {
                const size_t iqm = static_cast<size_t>(m) * nq + q;
                pg_bar += weights_mat[iqm] * pg_src[iqm];
              }
              scores_ptr[col_base + q] += cnt * (pg_bar - ind);
            }
          }
        }
      }
    }
  }

  return writable::list({
    "loglik"_nm = loglik_out,
    "score"_nm = gradient ? static_cast<SEXP>(scores_out) : R_NilValue,
    "theta_score"_nm = gradient ? static_cast<SEXP>(theta_score_out) : R_NilValue,
    "theta_curvature"_nm = gradient ? static_cast<SEXP>(theta_curv_out) : R_NilValue
  });
}

[[cpp11::register]]
writable::list mfrm_cpp_marginal_pairwise_expectations(list prob_list,
                                                       doubles_matrix<> obs_posterior,
                                                       integers row_i,
                                                       integers row_j,
                                                       integers adj_col1,
                                                       integers adj_col2) {
  const int n_pairs = row_i.size();
  writable::doubles exp_exact(n_pairs);
  writable::doubles exp_adjacent(n_pairs);
  if (n_pairs == 0) {
    return writable::list({
      "exp_exact"_nm = exp_exact,
      "exp_adjacent"_nm = exp_adjacent
    });
  }

  const int n_nodes = prob_list.size();
  doubles_matrix<> first_prob(prob_list[0]);
  const int n_obs = first_prob.nrow();
  const int k_cat = first_prob.ncol();
  const int n_adj = adj_col1.size();

  const int* ri_ptr = INTEGER(row_i);
  const int* rj_ptr = INTEGER(row_j);
  std::vector<int> slot_of_row(static_cast<size_t>(n_obs), -1);
  std::vector<int> used_rows;
  used_rows.reserve(static_cast<size_t>(n_obs));
  for (int p = 0; p < n_pairs; ++p) {
    const int i0 = ri_ptr[p] - 1;
    const int j0 = rj_ptr[p] - 1;
    if (slot_of_row[static_cast<size_t>(i0)] < 0) {
      slot_of_row[static_cast<size_t>(i0)] = static_cast<int>(used_rows.size());
      used_rows.push_back(i0);
    }
    if (slot_of_row[static_cast<size_t>(j0)] < 0) {
      slot_of_row[static_cast<size_t>(j0)] = static_cast<int>(used_rows.size());
      used_rows.push_back(j0);
    }
  }

  const size_t n_used = used_rows.size();
  std::vector<double> probs_flat(n_used * static_cast<size_t>(n_nodes) * static_cast<size_t>(k_cat));
  for (int q = 0; q < n_nodes; ++q) {
    doubles_matrix<> probs_q(prob_list[q]);
    for (size_t s = 0; s < n_used; ++s) {
      const int r = used_rows[s];
      double* dst = probs_flat.data() +
        (s * static_cast<size_t>(n_nodes) + static_cast<size_t>(q)) * static_cast<size_t>(k_cat);
      for (int k = 0; k < k_cat; ++k) {
        dst[k] = probs_q(r, k);
      }
    }
  }

  std::vector<int> adj1(static_cast<size_t>(n_adj));
  std::vector<int> adj2(static_cast<size_t>(n_adj));
  for (int a = 0; a < n_adj; ++a) {
    adj1[static_cast<size_t>(a)] = adj_col1[a] - 1;
    adj2[static_cast<size_t>(a)] = adj_col2[a] - 1;
  }

  for (int p = 0; p < n_pairs; ++p) {
    const int i0 = ri_ptr[p] - 1;
    const int j0 = rj_ptr[p] - 1;
    const size_t si = static_cast<size_t>(slot_of_row[static_cast<size_t>(i0)]);
    const size_t sj = static_cast<size_t>(slot_of_row[static_cast<size_t>(j0)]);
    double ee = 0.0;
    double ea = 0.0;
    for (int q = 0; q < n_nodes; ++q) {
      const double* pi_vec = probs_flat.data() +
        (si * static_cast<size_t>(n_nodes) + static_cast<size_t>(q)) * static_cast<size_t>(k_cat);
      const double* pj_vec = probs_flat.data() +
        (sj * static_cast<size_t>(n_nodes) + static_cast<size_t>(q)) * static_cast<size_t>(k_cat);
      double exact_q = 0.0;
      for (int k = 0; k < k_cat; ++k) {
        exact_q += pi_vec[k] * pj_vec[k];
      }
      double adjacent_q = 0.0;
      for (int a = 0; a < n_adj; ++a) {
        adjacent_q += pi_vec[adj1[static_cast<size_t>(a)]] * pj_vec[adj2[static_cast<size_t>(a)]];
      }
      const double post_q = obs_posterior(i0, q);
      ee += post_q * exact_q;
      ea += post_q * adjacent_q;
    }
    exp_exact[p] = ee;
    exp_adjacent[p] = ea;
  }

  return writable::list({
    "exp_exact"_nm = exp_exact,
    "exp_adjacent"_nm = exp_adjacent
  });
}

namespace {

struct PosteriorPersonObs {
  double delta;
  double alpha;
  const double* cum_steps;
  int y;
};

inline void eval_person_logpost_z(double z,
                                  double mu,
                                  double sigma,
                                  int n_cat,
                                  const std::vector<PosteriorPersonObs>& obs,
                                  std::vector<double>& log_num_buf,
                                  std::vector<double>& prob_buf,
                                  double& log_post,
                                  double* grad_z,
                                  double* neg_hess_z) {
  const double theta = mu + sigma * z;
  const int m_obs = static_cast<int>(obs.size());
  double ll = 0.0;
  double s_grad = 0.0;
  double s_info = 0.0;

  for (int j = 0; j < m_obs; ++j) {
    const PosteriorPersonObs& o = obs[static_cast<size_t>(j)];
    const double eta = theta + o.delta;
    const double alpha = o.alpha;
    const double* cs = o.cum_steps;
    double row_max = -std::numeric_limits<double>::infinity();
    for (int k = 0; k < n_cat; ++k) {
      const double u = alpha * (static_cast<double>(k) * eta - cs[k]);
      log_num_buf[static_cast<size_t>(k)] = u;
      if (u > row_max) row_max = u;
    }
    double denom = 0.0;
    for (int k = 0; k < n_cat; ++k) {
      const double e = std::exp(log_num_buf[static_cast<size_t>(k)] - row_max);
      prob_buf[static_cast<size_t>(k)] = e;
      denom += e;
    }
    ll += (log_num_buf[static_cast<size_t>(o.y)] - row_max) - std::log(denom);
    if (grad_z != nullptr || neg_hess_z != nullptr) {
      double ek = 0.0;
      double ek2 = 0.0;
      for (int k = 0; k < n_cat; ++k) {
        const double pk = prob_buf[static_cast<size_t>(k)] / denom;
        const double dk = static_cast<double>(k);
        ek += pk * dk;
        if (neg_hess_z != nullptr) ek2 += pk * dk * dk;
      }
      if (grad_z != nullptr) {
        s_grad += alpha * (static_cast<double>(o.y) - ek);
      }
      if (neg_hess_z != nullptr) {
        const double vk = std::max(0.0, ek2 - ek * ek);
        s_info += alpha * alpha * vk;
      }
    }
  }

  // -0.5 * log(2 * pi) = -0.91893853320467274178
  log_post = -0.5 * z * z - 0.91893853320467274178 + ll;
  if (grad_z != nullptr) {
    *grad_z = -z + sigma * s_grad;
  }
  if (neg_hess_z != nullptr) {
    *neg_hess_z = 1.0 + sigma * sigma * s_info;
  }
}

inline double integrate_gl_panel_post(double a,
                                      double b,
                                      double h_star,
                                      double mu,
                                      double sigma,
                                      int n_cat,
                                      const std::vector<PosteriorPersonObs>& obs,
                                      const double* gl_nodes,
                                      const double* gl_weights,
                                      int n_gl,
                                      std::vector<double>& log_num_buf,
                                      std::vector<double>& prob_buf) {
  if (!(b > a)) return 0.0;
  const double half = 0.5 * (b - a);
  const double mid = 0.5 * (a + b);
  double acc = 0.0;
  for (int q = 0; q < n_gl; ++q) {
    const double z = mid + half * gl_nodes[q];
    double lp = 0.0;
    eval_person_logpost_z(z, mu, sigma, n_cat, obs, log_num_buf, prob_buf, lp, nullptr, nullptr);
    acc += gl_weights[q] * std::exp(lp - h_star);
  }
  return half * acc;
}

inline double integrate_gl_panel_stdnorm(double a,
                                         double b,
                                         const double* gl_nodes,
                                         const double* gl_weights,
                                         int n_gl) {
  if (!(b > a)) return 0.0;
  const double half = 0.5 * (b - a);
  const double mid = 0.5 * (a + b);
  double acc = 0.0;
  for (int q = 0; q < n_gl; ++q) {
    const double z = mid + half * gl_nodes[q];
    acc += gl_weights[q] * Rf_dnorm4(z, 0.0, 1.0, 0);
  }
  return half * acc;
}

} // namespace

[[cpp11::register]]
writable::list mfrm_cpp_posterior_person_score(int n_cat,
                                               doubles_matrix<> loc_draws,
                                               doubles_matrix<> cum_step_draws,
                                               doubles_matrix<> slope_draws,
                                               doubles mu_draws,
                                               doubles sigma_draws,
                                               integers_matrix<> row_loc_idx,
                                               doubles facet_signs,
                                               integers row_owner,
                                               integers row_slope_idx,
                                               integers row_y,
                                               integers person_row_starts,
                                               integers person_rows,
                                               doubles probs,
                                               doubles cut_scores,
                                               doubles score_values,
                                               doubles gl_nodes,
                                               doubles gl_weights,
                                               bool compute_predictive,
                                               bool keep_person_draws) {
  const int n_draws = mu_draws.size();
  const int n_person = person_row_starts.size() - 1;
  const int n_facets = facet_signs.size();
  const int n_prob = probs.size();
  const int k_cut = cut_scores.size();
  const int n_dec_cols = (k_cut > 0) ? (2 * k_cut + 1) : 0;
  const int n_gl = gl_nodes.size();
  const int n_pred_rows = compute_predictive ? person_rows.size() : 0;
  const int k_max = n_cat - 1;

  writable::integers status_code(n_person);
  writable::strings error_message(n_person);
  writable::integers observed_count(n_person);
  writable::integers omitted_count(n_person);
  writable::doubles out_mean(n_person);
  writable::doubles out_within_sd(n_person);
  writable::doubles out_between_sd(n_person);
  writable::doubles out_total_sd(n_person);
  writable::doubles out_log_ml(n_person);
  writable::doubles_matrix<> out_quantiles(n_person, n_prob);
  writable::doubles_matrix<> out_dec_probs(n_person, n_dec_cols);
  writable::doubles_matrix<> out_dec_between_sd(n_person, n_dec_cols);
  writable::doubles_matrix<> out_pred_probs(n_pred_rows, compute_predictive ? n_cat : 0);
  writable::doubles out_pred_exp_sd(n_pred_rows);
  writable::doubles_matrix<> out_draw_mean(keep_person_draws ? n_person : 0, keep_person_draws ? n_draws : 0);
  writable::doubles_matrix<> out_draw_sd(keep_person_draws ? n_person : 0, keep_person_draws ? n_draws : 0);
  writable::doubles_matrix<> out_draw_log_z(keep_person_draws ? n_person : 0, keep_person_draws ? n_draws : 0);

  const double* gl_x = REAL(gl_nodes);
  const double* gl_w = REAL(gl_weights);
  const double* sc_vals = REAL(score_values);

  std::vector<double> log_num_buf(static_cast<size_t>(n_cat));
  std::vector<double> prob_buf(static_cast<size_t>(n_cat));
  std::vector<PosteriorPersonObs> obs_buf;

  std::vector<double> draw_mu_p(static_cast<size_t>(n_draws));
  std::vector<double> draw_var_p(static_cast<size_t>(n_draws));
  std::vector<double> draw_log_z_p(static_cast<size_t>(n_draws));
  std::vector<double> draw_z_star(static_cast<size_t>(n_draws));
  std::vector<double> draw_d1(static_cast<size_t>(n_draws));
  std::vector<double> draw_d2(static_cast<size_t>(n_draws));
  std::vector<double> draw_d3(static_cast<size_t>(n_draws));
  std::vector<double> draw_h_star(static_cast<size_t>(n_draws));
  std::vector<double> draw_i0(static_cast<size_t>(n_draws));
  std::vector<double> draw_panel0(static_cast<size_t>(n_draws));
  std::vector<double> draw_panel1(static_cast<size_t>(n_draws));
  std::vector<double> draw_panel2(static_cast<size_t>(n_draws));
  std::vector<double> draw_panel5(static_cast<size_t>(n_draws));
  std::vector<double> draw_panel6(static_cast<size_t>(n_draws));
  std::vector<double> draw_panel7(static_cast<size_t>(n_draws));
  std::vector<double> cum_steps_flat(static_cast<size_t>(cum_step_draws.nrow()) *
                                     static_cast<size_t>(cum_step_draws.ncol()));
  {
    const int nr_cs = cum_step_draws.nrow();
    const int nc_cs = cum_step_draws.ncol();
    for (int m = 0; m < nr_cs; ++m) {
      for (int c = 0; c < nc_cs; ++c) {
        cum_steps_flat[static_cast<size_t>(m) * static_cast<size_t>(nc_cs) + static_cast<size_t>(c)] =
          cum_step_draws(m, c);
      }
    }
  }
  const int nc_cs_stride = cum_step_draws.ncol();

  std::vector<double> dec_sum(static_cast<size_t>(n_dec_cols));
  std::vector<double> dec_draws(static_cast<size_t>(n_draws) * static_cast<size_t>(n_dec_cols));
  std::vector<double> band_tmp(static_cast<size_t>(k_cut + 1));

  std::vector<double> quad_z(static_cast<size_t>(8 * n_gl));
  std::vector<double> quad_w_norm(static_cast<size_t>(8 * n_gl));

  for (int p = 0; p < n_person; ++p) {
    const int r_start = person_row_starts[p];
    const int r_end = person_row_starts[p + 1];
    const int n_p_rows = r_end - r_start;

    std::vector<int> obs_rows_p;
    obs_rows_p.reserve(static_cast<size_t>(n_p_rows));
    int n_omit_p = 0;
    for (int idx = r_start; idx < r_end; ++idx) {
      const int r0 = person_rows[idx] - 1;
      if (row_y[r0] >= 0) {
        obs_rows_p.push_back(r0);
      } else {
        ++n_omit_p;
      }
    }
    const int m_obs = static_cast<int>(obs_rows_p.size());
    observed_count[p] = m_obs;
    omitted_count[p] = n_omit_p;

    auto mark_unavailable = [&](const std::string& msg) {
      status_code[p] = 2;
      error_message[p] = msg;
      out_mean[p] = NA_REAL;
      out_within_sd[p] = NA_REAL;
      out_between_sd[p] = NA_REAL;
      out_total_sd[p] = NA_REAL;
      out_log_ml[p] = NA_REAL;
      for (int qi = 0; qi < n_prob; ++qi) out_quantiles(p, qi) = NA_REAL;
      for (int dc = 0; dc < n_dec_cols; ++dc) {
        out_dec_probs(p, dc) = NA_REAL;
        out_dec_between_sd(p, dc) = NA_REAL;
      }
      if (compute_predictive) {
        for (int idx = r_start; idx < r_end; ++idx) {
          for (int k = 0; k < n_cat; ++k) out_pred_probs(idx, k) = NA_REAL;
          out_pred_exp_sd[idx] = NA_REAL;
        }
      }
      if (keep_person_draws) {
        for (int m = 0; m < n_draws; ++m) {
          out_draw_mean(p, m) = NA_REAL;
          out_draw_sd(p, m) = NA_REAL;
          out_draw_log_z(p, m) = NA_REAL;
        }
      }
    };

    if (n_draws < 1) {
      mark_unavailable("No calibration draws were supplied for scoring.");
      continue;
    }

    // Per-person predictive accumulators across rows r_start .. r_end - 1
    std::vector<double> pred_prob_acc;
    std::vector<double> pred_exp_draws;
    if (compute_predictive && n_p_rows > 0) {
      pred_prob_acc.assign(static_cast<size_t>(n_p_rows) * static_cast<size_t>(n_cat), 0.0);
      pred_exp_draws.assign(static_cast<size_t>(n_p_rows) * static_cast<size_t>(n_draws), 0.0);
    }
    if (n_dec_cols > 0) {
      std::fill(dec_sum.begin(), dec_sum.end(), 0.0);
    }

    bool failed = false;
    std::string fail_msg;

    if (m_obs == 0) {
      // Prior-only Person: exact normal mixture across calibration draws
      for (int m = 0; m < n_draws; ++m) {
        const double mu_m = mu_draws[m];
        const double sig_m = sigma_draws[m];
        if (!std::isfinite(mu_m) || !std::isfinite(sig_m) || sig_m <= 0.0) {
          failed = true;
          fail_msg = "Population mean or SD draw is not finite and positive.";
          break;
        }
        draw_mu_p[static_cast<size_t>(m)] = mu_m;
        draw_var_p[static_cast<size_t>(m)] = sig_m * sig_m;
        draw_log_z_p[static_cast<size_t>(m)] = 0.0;

        if (k_cut > 0) {
          for (int b = 0; b <= k_cut; ++b) {
            const double lo_x = (b == 0) ? -std::numeric_limits<double>::infinity() : cut_scores[b - 1];
            const double hi_x = (b == k_cut) ? std::numeric_limits<double>::infinity() : cut_scores[b];
            const double ua = std::isfinite(lo_x) ? ((lo_x - mu_m) / sig_m) : -std::numeric_limits<double>::infinity();
            const double ub = std::isfinite(hi_x) ? ((hi_x - mu_m) / sig_m) : std::numeric_limits<double>::infinity();
            double prob_b = 0.0;
            if (ub <= 0.0) {
              const double pb = Rf_pnorm5(ub, 0.0, 1.0, 1, 0);
              const double pa = std::isfinite(ua) ? Rf_pnorm5(ua, 0.0, 1.0, 1, 0) : 0.0;
              prob_b = pb - pa;
              if (std::isfinite(ua) && prob_b < 0.5 * pb) {
                prob_b = integrate_gl_panel_stdnorm(std::max(-12.0, ua), std::min(12.0, ub), gl_x, gl_w, n_gl);
              }
            } else if (ua >= 0.0) {
              const double qa = Rf_pnorm5(ua, 0.0, 1.0, 0, 0);
              const double qb = std::isfinite(ub) ? Rf_pnorm5(ub, 0.0, 1.0, 0, 0) : 0.0;
              prob_b = qa - qb;
              if (std::isfinite(ub) && prob_b < 0.5 * qa) {
                prob_b = integrate_gl_panel_stdnorm(std::max(-12.0, ua), std::min(12.0, ub), gl_x, gl_w, n_gl);
              }
            } else {
              const double left_p = std::isfinite(ua) ? (0.5 - Rf_pnorm5(ua, 0.0, 1.0, 1, 0)) : 0.5;
              const double right_p = std::isfinite(ub) ? (0.5 - Rf_pnorm5(ub, 0.0, 1.0, 0, 0)) : 0.5;
              prob_b = left_p + right_p;
            }
            band_tmp[static_cast<size_t>(b)] = std::max(0.0, std::min(1.0, prob_b));
          }
          double bsum = 0.0;
          for (int b = 0; b <= k_cut; ++b) bsum += band_tmp[static_cast<size_t>(b)];
          if (!std::isfinite(bsum) || std::abs(bsum - 1.0) > 1e-8) {
            failed = true;
            fail_msg = "Prior-only decision probability normalization failed.";
            break;
          }
          for (int b = 0; b <= k_cut; ++b) band_tmp[static_cast<size_t>(b)] /= bsum;
          double running_above = 0.0;
          for (int c = k_cut - 1; c >= 0; --c) {
            running_above += band_tmp[static_cast<size_t>(c + 1)];
            const double val = std::min(1.0, running_above);
            dec_draws[static_cast<size_t>(m) * static_cast<size_t>(n_dec_cols) + static_cast<size_t>(c)] = val;
            dec_sum[static_cast<size_t>(c)] += val;
          }
          for (int b = 0; b <= k_cut; ++b) {
            const int col = k_cut + b;
            const double val = band_tmp[static_cast<size_t>(b)];
            dec_draws[static_cast<size_t>(m) * static_cast<size_t>(n_dec_cols) + static_cast<size_t>(col)] = val;
            dec_sum[static_cast<size_t>(col)] += val;
          }
        }

        if (compute_predictive && n_p_rows > 0) {
          const double panels[9] = {-10.0, -6.0, -3.5, -1.5, 0.0, 1.5, 3.5, 6.0, 10.0};
          double i0_norm = 0.0;
          int q_slot = 0;
          for (int pan = 0; pan < 8; ++pan) {
            const double half = 0.5 * (panels[pan + 1] - panels[pan]);
            const double mid = 0.5 * (panels[pan + 1] + panels[pan]);
            for (int q = 0; q < n_gl; ++q) {
              const double z = mid + half * gl_x[q];
              const double u = half * gl_w[q] * Rf_dnorm4(z, 0.0, 1.0, 0);
              quad_z[static_cast<size_t>(q_slot)] = z;
              quad_w_norm[static_cast<size_t>(q_slot)] = u;
              i0_norm += u;
              ++q_slot;
            }
          }
          for (int q = 0; q < 8 * n_gl; ++q) quad_w_norm[static_cast<size_t>(q)] /= i0_norm;
          for (int r_local = 0; r_local < n_p_rows; ++r_local) {
            const int r0 = person_rows[r_start + r_local] - 1;
            double delta_r = 0.0;
            for (int f = 0; f < n_facets; ++f) {
              delta_r += facet_signs[f] * loc_draws(m, row_loc_idx(r0, f) - 1);
            }
            const double alpha_r = slope_draws(m, row_slope_idx[r0] - 1);
            const double* cs_r = cum_steps_flat.data() +
              static_cast<size_t>(m) * static_cast<size_t>(nc_cs_stride) +
              static_cast<size_t>(row_owner[r0] - 1) * static_cast<size_t>(n_cat);
            double exp_m = 0.0;
            for (int q = 0; q < 8 * n_gl; ++q) {
              const double wq = quad_w_norm[static_cast<size_t>(q)];
              const double eta = (mu_m + sig_m * quad_z[static_cast<size_t>(q)]) + delta_r;
              double rmax = -std::numeric_limits<double>::infinity();
              for (int k = 0; k < n_cat; ++k) {
                const double u = alpha_r * (static_cast<double>(k) * eta - cs_r[k]);
                log_num_buf[static_cast<size_t>(k)] = u;
                if (u > rmax) rmax = u;
              }
              double den = 0.0;
              for (int k = 0; k < n_cat; ++k) {
                const double e = std::exp(log_num_buf[static_cast<size_t>(k)] - rmax);
                prob_buf[static_cast<size_t>(k)] = e;
                den += e;
              }
              for (int k = 0; k < n_cat; ++k) {
                const double pk = prob_buf[static_cast<size_t>(k)] / den;
                pred_prob_acc[static_cast<size_t>(r_local) * static_cast<size_t>(n_cat) + static_cast<size_t>(k)] += wq * pk;
                exp_m += wq * pk * sc_vals[k];
              }
            }
            pred_exp_draws[static_cast<size_t>(r_local) * static_cast<size_t>(n_draws) + static_cast<size_t>(m)] = exp_m;
          }
        }
      }

      if (failed) {
        mark_unavailable(fail_msg);
        continue;
      }

      double mean_p = 0.0;
      double within_var = 0.0;
      for (int m = 0; m < n_draws; ++m) {
        mean_p += draw_mu_p[static_cast<size_t>(m)];
        within_var += draw_var_p[static_cast<size_t>(m)];
      }
      mean_p /= static_cast<double>(n_draws);
      within_var /= static_cast<double>(n_draws);
      double between_var = 0.0;
      for (int m = 0; m < n_draws; ++m) {
        const double diff = draw_mu_p[static_cast<size_t>(m)] - mean_p;
        between_var += diff * diff;
      }
      between_var /= static_cast<double>(n_draws);
      const double total_sd = std::sqrt(within_var + between_var);

      // Exact normal-mixture quantile inversion
      double x_min = std::numeric_limits<double>::infinity();
      double x_max = -std::numeric_limits<double>::infinity();
      for (int m = 0; m < n_draws; ++m) {
        const double mu_m = draw_mu_p[static_cast<size_t>(m)];
        const double sig_m = std::sqrt(draw_var_p[static_cast<size_t>(m)]);
        x_min = std::min(x_min, mu_m - 12.0 * sig_m);
        x_max = std::max(x_max, mu_m + 12.0 * sig_m);
      }

      for (int qi = 0; qi < n_prob; ++qi) {
        const double pr = probs[qi];
        if (n_draws == 1) {
          const double mu0 = draw_mu_p[0];
          const double sig0 = std::sqrt(draw_var_p[0]);
          out_quantiles(p, qi) = (pr <= 0.5)
            ? (mu0 + sig0 * Rf_qnorm5(pr, 0.0, 1.0, 1, 0))
            : (mu0 - sig0 * Rf_qnorm5(1.0 - pr, 0.0, 1.0, 1, 0));
          continue;
        }
        double lo = x_min;
        double hi = x_max;
        double x = mean_p + total_sd * Rf_qnorm5(pr, 0.0, 1.0, 1, 0);
        if (!(x > lo && x < hi)) x = 0.5 * (lo + hi);
        double prev_width = hi - lo;
        const double tol_err = std::min(1e-12, 1e-8 * std::min(pr, 1.0 - pr));
        for (int it = 0; it < 60; ++it) {
          double sum_tail = 0.0;
          double sum_pdf = 0.0;
          for (int m = 0; m < n_draws; ++m) {
            const double mu_m = draw_mu_p[static_cast<size_t>(m)];
            const double sig_m = std::sqrt(draw_var_p[static_cast<size_t>(m)]);
            const double z = (x - mu_m) / sig_m;
            if (pr <= 0.5) {
              sum_tail += Rf_pnorm5(z, 0.0, 1.0, 1, 0);
            } else {
              sum_tail += Rf_pnorm5(z, 0.0, 1.0, 0, 0);
            }
            sum_pdf += Rf_dnorm4(z, 0.0, 1.0, 0) / sig_m;
          }
          const double err = (pr <= 0.5)
            ? (sum_tail / static_cast<double>(n_draws) - pr)
            : ((1.0 - pr) - sum_tail / static_cast<double>(n_draws));
          const double pdf = sum_pdf / static_cast<double>(n_draws);
          if (err < 0.0 && x > lo) lo = x;
          if (err > 0.0 && x < hi) hi = x;
          if (std::abs(err) <= tol_err || (hi - lo) <= 1e-13 * std::max(1.0, std::abs(x))) {
            break;
          }
          double x_new = (pdf > 0.0) ? (x - err / pdf) : (0.5 * (lo + hi));
          if (!(x_new > lo && x_new < hi) || std::abs(x_new - x) > 0.5 * prev_width) {
            x_new = 0.5 * (lo + hi);
          }
          prev_width = hi - lo;
          x = x_new;
        }
        out_quantiles(p, qi) = x;
      }

      status_code[p] = 1;
      error_message[p] = "";
      out_mean[p] = mean_p;
      out_within_sd[p] = std::sqrt(within_var);
      out_between_sd[p] = std::sqrt(between_var);
      out_total_sd[p] = total_sd;
      out_log_ml[p] = 0.0;
    } else {
      // Observed Person: exact bracketed Newton mode + 6-panel GL quadrature
      obs_buf.resize(static_cast<size_t>(m_obs));
      for (int m = 0; m < n_draws; ++m) {
        const double mu_m = mu_draws[m];
        const double sig_m = sigma_draws[m];
        if (!std::isfinite(mu_m) || !std::isfinite(sig_m) || sig_m <= 0.0) {
          failed = true;
          fail_msg = "Population mean or SD draw is not finite and positive.";
          break;
        }

        double s1 = 0.0;
        double a_sum = 0.0;
        double a_max = 0.0;
        for (int j = 0; j < m_obs; ++j) {
          const int r0 = obs_rows_p[static_cast<size_t>(j)];
          double delta_j = 0.0;
          for (int f = 0; f < n_facets; ++f) {
            delta_j += facet_signs[f] * loc_draws(m, row_loc_idx(r0, f) - 1);
          }
          const double alpha_j = slope_draws(m, row_slope_idx[r0] - 1);
          if (!std::isfinite(delta_j) || !std::isfinite(alpha_j) || alpha_j <= 0.0) {
            failed = true;
            fail_msg = "Facet location or slope draw is not finite and positive.";
            break;
          }
          const double* cs_j = cum_steps_flat.data() +
            static_cast<size_t>(m) * static_cast<size_t>(nc_cs_stride) +
            static_cast<size_t>(row_owner[r0] - 1) * static_cast<size_t>(n_cat);
          const int y_j = row_y[r0];
          obs_buf[static_cast<size_t>(j)] = PosteriorPersonObs{delta_j, alpha_j, cs_j, y_j};
          s1 += alpha_j * static_cast<double>(y_j);
          a_sum += alpha_j;
          if (alpha_j > a_max) a_max = alpha_j;
        }
        if (failed) break;

        double z_lo = sig_m * (s1 - static_cast<double>(k_max) * a_sum);
        double z_hi = sig_m * s1;
        if (!(z_hi > z_lo)) {
          z_lo -= 1.0;
          z_hi += 1.0;
        }
        double z = (0.0 > z_lo && 0.0 < z_hi) ? 0.0 : (0.5 * (z_lo + z_hi));
        double h_val = 0.0;
        double grad = 0.0;
        double neg_hess = 1.0;

        for (int it = 0; it < 60; ++it) {
          eval_person_logpost_z(z, mu_m, sig_m, n_cat, obs_buf, log_num_buf, prob_buf,
                                h_val, &grad, &neg_hess);
          if (grad > 0.0 && z > z_lo) z_lo = z;
          if (grad < 0.0 && z < z_hi) z_hi = z;
          if (std::abs(grad) <= 1e-13 || (z_hi - z_lo) <= 1e-14 * std::max(1.0, std::abs(z))) {
            break;
          }
          double z_new = z + grad / neg_hess;
          if (!(z_new > z_lo && z_new < z_hi)) {
            z_new = 0.5 * (z_lo + z_hi);
          }
          z = z_new;
        }
        eval_person_logpost_z(z, mu_m, sig_m, n_cat, obs_buf, log_num_buf, prob_buf,
                              h_val, &grad, &neg_hess);
        if (!std::isfinite(h_val) || !std::isfinite(grad) || std::abs(grad) > 1e-8 ||
            !std::isfinite(neg_hess) || neg_hess < 1.0) {
          failed = true;
          fail_msg = "Conditional posterior mode optimization did not converge.";
          break;
        }

        const double s_eff = std::min(1.0 / std::sqrt(neg_hess), 1.5 / (sig_m * a_max));
        const double d1 = std::max(1e-4, std::min(1.0, 2.5 * s_eff));
        const double d2 = std::max(2.0 * d1, std::min(2.6, 7.0 * s_eff));
        const double d3 = std::max(1.6 * d2, std::min(5.2, 4.5 * std::sqrt(s_eff)));
        const double panels[9] = {
          z - 10.0, z - d3, z - d2, z - d1, z, z + d1, z + d2, z + d3, z + 10.0
        };
        double panel_mass[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double i1 = 0.0;
        double i2 = 0.0;
        int q_slot = 0;
        for (int pan = 0; pan < 8; ++pan) {
          const double half = 0.5 * (panels[pan + 1] - panels[pan]);
          const double mid = 0.5 * (panels[pan + 1] + panels[pan]);
          double pm = 0.0;
          for (int q = 0; q < n_gl; ++q) {
            const double zq = mid + half * gl_x[q];
            double lpq = 0.0;
            eval_person_logpost_z(zq, mu_m, sig_m, n_cat, obs_buf, log_num_buf, prob_buf,
                                  lpq, nullptr, nullptr);
            const double u = half * gl_w[q] * std::exp(lpq - h_val);
            pm += u;
            const double dz = zq - z;
            i1 += dz * u;
            i2 += dz * dz * u;
            quad_z[static_cast<size_t>(q_slot)] = zq;
            quad_w_norm[static_cast<size_t>(q_slot)] = u;
            ++q_slot;
          }
          panel_mass[pan] = pm;
        }
        const double i0 = panel_mass[0] + panel_mass[1] + panel_mass[2] + panel_mass[3] +
          panel_mass[4] + panel_mass[5] + panel_mass[6] + panel_mass[7];
        if (!std::isfinite(i0) || i0 <= 0.0) {
          failed = true;
          fail_msg = "Conditional posterior normalizing constant is non-positive or non-finite.";
          break;
        }
        const double dz_mean = i1 / i0;
        const double var_z = i2 / i0 - dz_mean * dz_mean;
        if (!std::isfinite(var_z) || var_z <= 0.0) {
          failed = true;
          fail_msg = "Conditional posterior variance is non-positive or non-finite.";
          break;
        }

        draw_mu_p[static_cast<size_t>(m)] = mu_m + sig_m * (z + dz_mean);
        draw_var_p[static_cast<size_t>(m)] = sig_m * sig_m * var_z;
        draw_log_z_p[static_cast<size_t>(m)] = h_val + std::log(i0);
        draw_z_star[static_cast<size_t>(m)] = z;
        draw_d1[static_cast<size_t>(m)] = d1;
        draw_d2[static_cast<size_t>(m)] = d2;
        draw_d3[static_cast<size_t>(m)] = d3;
        draw_h_star[static_cast<size_t>(m)] = h_val;
        draw_i0[static_cast<size_t>(m)] = i0;
        draw_panel0[static_cast<size_t>(m)] = panel_mass[0];
        draw_panel1[static_cast<size_t>(m)] = panel_mass[1];
        draw_panel2[static_cast<size_t>(m)] = panel_mass[2];
        draw_panel5[static_cast<size_t>(m)] = panel_mass[5];
        draw_panel6[static_cast<size_t>(m)] = panel_mass[6];
        draw_panel7[static_cast<size_t>(m)] = panel_mass[7];

        // Positive-weight 8-panel cut-score band integration
        if (k_cut > 0) {
          for (int b = 0; b <= k_cut; ++b) {
            const double za = (b == 0)
              ? panels[0]
              : std::max(panels[0], std::min(panels[8], (cut_scores[b - 1] - mu_m) / sig_m));
            const double zb = (b == k_cut)
              ? panels[8]
              : std::max(panels[0], std::min(panels[8], (cut_scores[b] - mu_m) / sig_m));
            double b_mass = 0.0;
            if (zb > za) {
              for (int pan = 0; pan < 8; ++pan) {
                const double sub_a = std::max(za, panels[pan]);
                const double sub_b = std::min(zb, panels[pan + 1]);
                if (sub_b > sub_a) {
                  if (sub_a == panels[pan] && sub_b == panels[pan + 1]) {
                    b_mass += panel_mass[pan];
                  } else {
                    b_mass += integrate_gl_panel_post(sub_a, sub_b, h_val, mu_m, sig_m, n_cat,
                                                      obs_buf, gl_x, gl_w, n_gl, log_num_buf, prob_buf);
                  }
                }
              }
            }
            band_tmp[static_cast<size_t>(b)] = std::max(0.0, b_mass / i0);
          }
          double bsum = 0.0;
          for (int b = 0; b <= k_cut; ++b) bsum += band_tmp[static_cast<size_t>(b)];
          if (!std::isfinite(bsum) || std::abs(bsum - 1.0) > 1e-8) {
            failed = true;
            fail_msg = "Decision probability band partition failed accuracy check.";
            break;
          }
          for (int b = 0; b <= k_cut; ++b) band_tmp[static_cast<size_t>(b)] /= bsum;
          double running_above = 0.0;
          for (int c = k_cut - 1; c >= 0; --c) {
            running_above += band_tmp[static_cast<size_t>(c + 1)];
            const double val = std::min(1.0, running_above);
            dec_draws[static_cast<size_t>(m) * static_cast<size_t>(n_dec_cols) + static_cast<size_t>(c)] = val;
            dec_sum[static_cast<size_t>(c)] += val;
          }
          for (int b = 0; b <= k_cut; ++b) {
            const int col = k_cut + b;
            const double val = band_tmp[static_cast<size_t>(b)];
            dec_draws[static_cast<size_t>(m) * static_cast<size_t>(n_dec_cols) + static_cast<size_t>(col)] = val;
            dec_sum[static_cast<size_t>(col)] += val;
          }
        }

        // Posterior predictive category probabilities across all rows of Person p
        if (compute_predictive && n_p_rows > 0) {
          for (int q = 0; q < 8 * n_gl; ++q) quad_w_norm[static_cast<size_t>(q)] /= i0;
          for (int r_local = 0; r_local < n_p_rows; ++r_local) {
            const int r0 = person_rows[r_start + r_local] - 1;
            double delta_r = 0.0;
            for (int f = 0; f < n_facets; ++f) {
              delta_r += facet_signs[f] * loc_draws(m, row_loc_idx(r0, f) - 1);
            }
            const double alpha_r = slope_draws(m, row_slope_idx[r0] - 1);
            const double* cs_r = cum_steps_flat.data() +
              static_cast<size_t>(m) * static_cast<size_t>(nc_cs_stride) +
              static_cast<size_t>(row_owner[r0] - 1) * static_cast<size_t>(n_cat);
            double exp_m = 0.0;
            for (int q = 0; q < 8 * n_gl; ++q) {
              const double wq = quad_w_norm[static_cast<size_t>(q)];
              const double eta = (mu_m + sig_m * quad_z[static_cast<size_t>(q)]) + delta_r;
              double rmax = -std::numeric_limits<double>::infinity();
              for (int k = 0; k < n_cat; ++k) {
                const double u = alpha_r * (static_cast<double>(k) * eta - cs_r[k]);
                log_num_buf[static_cast<size_t>(k)] = u;
                if (u > rmax) rmax = u;
              }
              double den = 0.0;
              for (int k = 0; k < n_cat; ++k) {
                const double e = std::exp(log_num_buf[static_cast<size_t>(k)] - rmax);
                prob_buf[static_cast<size_t>(k)] = e;
                den += e;
              }
              for (int k = 0; k < n_cat; ++k) {
                const double pk = prob_buf[static_cast<size_t>(k)] / den;
                pred_prob_acc[static_cast<size_t>(r_local) * static_cast<size_t>(n_cat) + static_cast<size_t>(k)] += wq * pk;
                exp_m += wq * pk * sc_vals[k];
              }
            }
            pred_exp_draws[static_cast<size_t>(r_local) * static_cast<size_t>(n_draws) + static_cast<size_t>(m)] = exp_m;
          }
        }
      }

      if (failed) {
        mark_unavailable(fail_msg);
        continue;
      }

      // Mixture moments and log marginal likelihood
      double mean_p = 0.0;
      double within_var = 0.0;
      double max_lz = -std::numeric_limits<double>::infinity();
      for (int m = 0; m < n_draws; ++m) {
        mean_p += draw_mu_p[static_cast<size_t>(m)];
        within_var += draw_var_p[static_cast<size_t>(m)];
        if (draw_log_z_p[static_cast<size_t>(m)] > max_lz) max_lz = draw_log_z_p[static_cast<size_t>(m)];
      }
      mean_p /= static_cast<double>(n_draws);
      within_var /= static_cast<double>(n_draws);

      double between_var = 0.0;
      double sum_exp_lz = 0.0;
      for (int m = 0; m < n_draws; ++m) {
        const double diff = draw_mu_p[static_cast<size_t>(m)] - mean_p;
        between_var += diff * diff;
        sum_exp_lz += std::exp(draw_log_z_p[static_cast<size_t>(m)] - max_lz);
      }
      between_var /= static_cast<double>(n_draws);
      const double total_sd = std::sqrt(within_var + between_var);
      const double log_ml = max_lz + std::log(sum_exp_lz) - std::log(static_cast<double>(n_draws));

      // Continuous mixture CDF/PDF inversion across calibration draws
      double x_min = std::numeric_limits<double>::infinity();
      double x_max = -std::numeric_limits<double>::infinity();
      for (int m = 0; m < n_draws; ++m) {
        const double mu_m = mu_draws[m];
        const double sig_m = sigma_draws[m];
        const double zs = draw_z_star[static_cast<size_t>(m)];
        x_min = std::min(x_min, mu_m + sig_m * (zs - 10.0));
        x_max = std::max(x_max, mu_m + sig_m * (zs + 10.0));
      }

      for (int qi = 0; qi < n_prob; ++qi) {
        const double pr = probs[qi];
        double lo = x_min;
        double hi = x_max;
        double x = mean_p + total_sd * Rf_qnorm5(pr, 0.0, 1.0, 1, 0);
        if (!(x > lo && x < hi)) x = 0.5 * (lo + hi);
        double prev_width = hi - lo;
        const double tol_err = std::min(1e-12, 1e-8 * std::min(pr, 1.0 - pr));
        double final_err = 1.0;

        for (int it = 0; it < 60; ++it) {
          double sum_tail = 0.0;
          double sum_pdf = 0.0;

          for (int m = 0; m < n_draws; ++m) {
            const double mu_m = mu_draws[m];
            const double sig_m = sigma_draws[m];
            const double zs = draw_z_star[static_cast<size_t>(m)];
            const double zx = (x - mu_m) / sig_m;
            if (zx <= zs - 10.0) {
              if (pr > 0.5) sum_tail += 1.0;
              continue;
            }
            if (zx >= zs + 10.0) {
              if (pr <= 0.5) sum_tail += 1.0;
              continue;
            }

            // Reconstruct obs_buf for draw m to evaluate single active sub-panel + exact PDF
            for (int j = 0; j < m_obs; ++j) {
              const int r0 = obs_rows_p[static_cast<size_t>(j)];
              double delta_j = 0.0;
              for (int f = 0; f < n_facets; ++f) {
                delta_j += facet_signs[f] * loc_draws(m, row_loc_idx(r0, f) - 1);
              }
              const double alpha_j = slope_draws(m, row_slope_idx[r0] - 1);
              const double* cs_j = cum_steps_flat.data() +
                static_cast<size_t>(m) * static_cast<size_t>(nc_cs_stride) +
                static_cast<size_t>(row_owner[r0] - 1) * static_cast<size_t>(n_cat);
              obs_buf[static_cast<size_t>(j)] = PosteriorPersonObs{delta_j, alpha_j, cs_j, row_y[r0]};
            }

            const double d1 = draw_d1[static_cast<size_t>(m)];
            const double d2 = draw_d2[static_cast<size_t>(m)];
            const double d3 = draw_d3[static_cast<size_t>(m)];
            const double hs = draw_h_star[static_cast<size_t>(m)];
            const double i0 = draw_i0[static_cast<size_t>(m)];

            double lp_x = 0.0;
            eval_person_logpost_z(zx, mu_m, sig_m, n_cat, obs_buf, log_num_buf, prob_buf,
                                  lp_x, nullptr, nullptr);
            sum_pdf += std::exp(lp_x - hs) / (sig_m * i0);

            double f_m = 0.0;
            double u_m = 0.0;
            if (zx <= zs - d3) {
              const double l_mass = integrate_gl_panel_post(zs - 10.0, zx, hs, mu_m, sig_m, n_cat,
                                                            obs_buf, gl_x, gl_w, n_gl, log_num_buf, prob_buf);
              f_m = std::min(1.0, std::max(0.0, l_mass / i0));
              u_m = 1.0 - f_m;
            } else if (zx <= zs - d2) {
              const double l_mass = draw_panel0[static_cast<size_t>(m)] +
                integrate_gl_panel_post(zs - d3, zx, hs, mu_m, sig_m, n_cat,
                                        obs_buf, gl_x, gl_w, n_gl, log_num_buf, prob_buf);
              f_m = std::min(1.0, std::max(0.0, l_mass / i0));
              u_m = 1.0 - f_m;
            } else if (zx <= zs - d1) {
              const double l_mass = draw_panel0[static_cast<size_t>(m)] +
                draw_panel1[static_cast<size_t>(m)] +
                integrate_gl_panel_post(zs - d2, zx, hs, mu_m, sig_m, n_cat,
                                        obs_buf, gl_x, gl_w, n_gl, log_num_buf, prob_buf);
              f_m = std::min(1.0, std::max(0.0, l_mass / i0));
              u_m = 1.0 - f_m;
            } else if (zx <= zs) {
              const double l_mass = draw_panel0[static_cast<size_t>(m)] +
                draw_panel1[static_cast<size_t>(m)] +
                draw_panel2[static_cast<size_t>(m)] +
                integrate_gl_panel_post(zs - d1, zx, hs, mu_m, sig_m, n_cat,
                                        obs_buf, gl_x, gl_w, n_gl, log_num_buf, prob_buf);
              f_m = std::min(1.0, std::max(0.0, l_mass / i0));
              u_m = 1.0 - f_m;
            } else if (zx < zs + d1) {
              const double r_mass = integrate_gl_panel_post(zx, zs + d1, hs, mu_m, sig_m, n_cat,
                                                            obs_buf, gl_x, gl_w, n_gl, log_num_buf, prob_buf) +
                draw_panel5[static_cast<size_t>(m)] +
                draw_panel6[static_cast<size_t>(m)] +
                draw_panel7[static_cast<size_t>(m)];
              u_m = std::min(1.0, std::max(0.0, r_mass / i0));
              f_m = 1.0 - u_m;
            } else if (zx < zs + d2) {
              const double r_mass = integrate_gl_panel_post(zx, zs + d2, hs, mu_m, sig_m, n_cat,
                                                            obs_buf, gl_x, gl_w, n_gl, log_num_buf, prob_buf) +
                draw_panel6[static_cast<size_t>(m)] +
                draw_panel7[static_cast<size_t>(m)];
              u_m = std::min(1.0, std::max(0.0, r_mass / i0));
              f_m = 1.0 - u_m;
            } else if (zx < zs + d3) {
              const double r_mass = integrate_gl_panel_post(zx, zs + d3, hs, mu_m, sig_m, n_cat,
                                                            obs_buf, gl_x, gl_w, n_gl, log_num_buf, prob_buf) +
                draw_panel7[static_cast<size_t>(m)];
              u_m = std::min(1.0, std::max(0.0, r_mass / i0));
              f_m = 1.0 - u_m;
            } else {
              const double r_mass = integrate_gl_panel_post(zx, zs + 10.0, hs, mu_m, sig_m, n_cat,
                                                            obs_buf, gl_x, gl_w, n_gl, log_num_buf, prob_buf);
              u_m = std::min(1.0, std::max(0.0, r_mass / i0));
              f_m = 1.0 - u_m;
            }

            if (pr <= 0.5) {
              sum_tail += f_m;
            } else {
              sum_tail += u_m;
            }
          }

          const double err = (pr <= 0.5)
            ? (sum_tail / static_cast<double>(n_draws) - pr)
            : ((1.0 - pr) - sum_tail / static_cast<double>(n_draws));
          const double pdf = sum_pdf / static_cast<double>(n_draws);
          final_err = std::abs(err);

          if (err < 0.0 && x > lo) lo = x;
          if (err > 0.0 && x < hi) hi = x;
          if (final_err <= tol_err || (hi - lo) <= 1e-13 * std::max(1.0, std::abs(x))) {
            break;
          }
          double x_new = (pdf > 0.0) ? (x - err / pdf) : (0.5 * (lo + hi));
          if (!(x_new > lo && x_new < hi) || std::abs(x_new - x) > 0.5 * prev_width) {
            x_new = 0.5 * (lo + hi);
          }
          prev_width = hi - lo;
          x = x_new;
        }

        if (final_err > 1e-7 && (hi - lo) > 1e-9 * std::max(1.0, std::abs(x))) {
          failed = true;
          fail_msg = "Continuous posterior quantile root-finding did not reach target accuracy.";
          break;
        }
        out_quantiles(p, qi) = x;
      }

      if (failed) {
        mark_unavailable(fail_msg);
        continue;
      }

      status_code[p] = 0;
      error_message[p] = "";
      out_mean[p] = mean_p;
      out_within_sd[p] = std::sqrt(within_var);
      out_between_sd[p] = std::sqrt(between_var);
      out_total_sd[p] = total_sd;
      out_log_ml[p] = log_ml;
    }

    // Store decision probabilities & between-draw SDs
    if (n_dec_cols > 0) {
      for (int dc = 0; dc < n_dec_cols; ++dc) {
        const double mean_dec = std::min(1.0, std::max(0.0, dec_sum[static_cast<size_t>(dc)] / static_cast<double>(n_draws)));
        out_dec_probs(p, dc) = mean_dec;
        double var_dec = 0.0;
        for (int m = 0; m < n_draws; ++m) {
          const double dv = dec_draws[static_cast<size_t>(m) * static_cast<size_t>(n_dec_cols) + static_cast<size_t>(dc)] - mean_dec;
          var_dec += dv * dv;
        }
        out_dec_between_sd(p, dc) = std::sqrt(var_dec / static_cast<double>(n_draws));
      }
    }

    // Store posterior predictive probabilities & expected-score between-draw SD
    if (compute_predictive && n_p_rows > 0) {
      for (int r_local = 0; r_local < n_p_rows; ++r_local) {
        const int out_r = r_start + r_local;
        double psum = 0.0;
        for (int k = 0; k < n_cat; ++k) {
          const double pk = pred_prob_acc[static_cast<size_t>(r_local) * static_cast<size_t>(n_cat) + static_cast<size_t>(k)] /
            static_cast<double>(n_draws);
          out_pred_probs(out_r, k) = pk;
          psum += pk;
        }
        if (psum > 0.0) {
          for (int k = 0; k < n_cat; ++k) out_pred_probs(out_r, k) /= psum;
        }
        double mean_exp = 0.0;
        for (int m = 0; m < n_draws; ++m) {
          mean_exp += pred_exp_draws[static_cast<size_t>(r_local) * static_cast<size_t>(n_draws) + static_cast<size_t>(m)];
        }
        mean_exp /= static_cast<double>(n_draws);
        double var_exp = 0.0;
        for (int m = 0; m < n_draws; ++m) {
          const double dv = pred_exp_draws[static_cast<size_t>(r_local) * static_cast<size_t>(n_draws) + static_cast<size_t>(m)] - mean_exp;
          var_exp += dv * dv;
        }
        out_pred_exp_sd[out_r] = std::sqrt(var_exp / static_cast<double>(n_draws));
      }
    }

    // Optional per-draw person moments
    if (keep_person_draws) {
      for (int m = 0; m < n_draws; ++m) {
        out_draw_mean(p, m) = draw_mu_p[static_cast<size_t>(m)];
        out_draw_sd(p, m) = std::sqrt(draw_var_p[static_cast<size_t>(m)]);
        out_draw_log_z(p, m) = draw_log_z_p[static_cast<size_t>(m)];
      }
    }
  }

  return writable::list({
    "status_code"_nm = status_code,
    "error_message"_nm = error_message,
    "observed_count"_nm = observed_count,
    "omitted_count"_nm = omitted_count,
    "mean"_nm = out_mean,
    "within_sd"_nm = out_within_sd,
    "between_sd"_nm = out_between_sd,
    "total_sd"_nm = out_total_sd,
    "log_marginal_lik"_nm = out_log_ml,
    "quantiles"_nm = out_quantiles,
    "decision_probs"_nm = out_dec_probs,
    "decision_between_sd"_nm = out_dec_between_sd,
    "pred_probs"_nm = out_pred_probs,
    "pred_expected_mean_sd"_nm = out_pred_exp_sd,
    "draw_mean"_nm = out_draw_mean,
    "draw_sd"_nm = out_draw_sd,
    "draw_log_z"_nm = out_draw_log_z
  });
}

[[cpp11::register]]
writable::list mfrm_cpp_drift_mml_eval(doubles tilde_psi,
                                       integers person_ptr,
                                       integers r_obs0,
                                       integers x_obs0,
                                       integers loc_idx0,
                                       int n_loc_blocks,
                                       int R,
                                       int K,
                                       int u_tilde_offset,
                                       bool step_rater,
                                       doubles quad_nodes,
                                       doubles quad_weights,
                                       double D,
                                       bool need_hessian) {
  const int n_tilde = tilde_psi.size();
  const int J = person_ptr.size() - 1;
  const int N = r_obs0.size();
  const int M = quad_nodes.size();
  const int n_u = K - 2;

  if (J < 1 || N < 1 || M < 1 || K < 2 || R < 1 || n_loc_blocks < 1 ||
      x_obs0.size() != N || quad_weights.size() != M ||
      loc_idx0.size() != N * n_loc_blocks || person_ptr[0] != 0 ||
      person_ptr[J] != N) {
    stop("Invalid dimensions passed to `mfrm_cpp_drift_mml_eval()`.");
  }

  std::vector<double> psi_v(static_cast<size_t>(n_tilde));
  for (int p = 0; p < n_tilde; ++p) {
    psi_v[static_cast<size_t>(p)] = tilde_psi[p];
  }

  // Step cumulative sums g_mat (R x K)
  std::vector<double> g_mat(static_cast<size_t>(R) * static_cast<size_t>(K), 0.0);
  if (n_u > 0) {
    if (step_rater) {
      for (int r = 0; r < R; ++r) {
        double u_sum = 0.0;
        double cum = 0.0;
        double* g_r = &g_mat[static_cast<size_t>(r) * static_cast<size_t>(K)];
        for (int m = 0; m < n_u; ++m) {
          const double u_rm = psi_v[static_cast<size_t>(u_tilde_offset + r * n_u + m)];
          u_sum += u_rm;
          cum += u_rm;
          g_r[static_cast<size_t>(m + 1)] = cum;
        }
        g_r[static_cast<size_t>(K - 1)] = cum - u_sum;
      }
    } else {
      std::vector<double> g_shared(static_cast<size_t>(K), 0.0);
      double u_sum = 0.0;
      double cum = 0.0;
      for (int m = 0; m < n_u; ++m) {
        const double u_m = psi_v[static_cast<size_t>(u_tilde_offset + m)];
        u_sum += u_m;
        cum += u_m;
        g_shared[static_cast<size_t>(m + 1)] = cum;
      }
      g_shared[static_cast<size_t>(K - 1)] = cum - u_sum;
      for (int r = 0; r < R; ++r) {
        std::copy(g_shared.begin(), g_shared.end(),
                  &g_mat[static_cast<size_t>(r) * static_cast<size_t>(K)]);
      }
    }
  }

  // Observation location linear predictors & contiguous index copies
  std::vector<int> r_v(static_cast<size_t>(N));
  std::vector<int> x_v(static_cast<size_t>(N));
  std::vector<int> loc_v(static_cast<size_t>(N) * static_cast<size_t>(n_loc_blocks));
  std::vector<double> eta_loc(static_cast<size_t>(N), 0.0);
  for (int i = 0; i < N; ++i) {
    r_v[static_cast<size_t>(i)] = r_obs0[i];
    x_v[static_cast<size_t>(i)] = x_obs0[i];
    double acc = 0.0;
    const int base_l = i * n_loc_blocks;
    for (int b = 0; b < n_loc_blocks; ++b) {
      const int idx = loc_idx0[base_l + b];
      loc_v[static_cast<size_t>(base_l + b)] = idx;
      acc += psi_v[static_cast<size_t>(idx)];
    }
    eta_loc[static_cast<size_t>(i)] = acc;
  }

  std::vector<double> nodes_v(static_cast<size_t>(M));
  std::vector<double> log_w(static_cast<size_t>(M));
  for (int q = 0; q < M; ++q) {
    nodes_v[static_cast<size_t>(q)] = quad_nodes[q];
    log_w[static_cast<size_t>(q)] = std::log(quad_weights[q]);
  }

  int max_obs_j = 1;
  for (int j = 0; j < J; ++j) {
    const int cnt = person_ptr[j + 1] - person_ptr[j];
    if (cnt > max_obs_j) max_obs_j = cnt;
  }

  double ll_mml = 0.0;
  const double D2 = D * D;
  std::vector<double> grad_tilde_v(static_cast<size_t>(n_tilde), 0.0);
  writable::doubles_matrix<> out_pi_jm(J, M);
  writable::doubles out_eap(J);
  writable::doubles out_psd(J);

  std::vector<double> p_buf(
      static_cast<size_t>(max_obs_j) * static_cast<size_t>(M) * static_cast<size_t>(K), 0.0);
  std::vector<double> log_L_j(static_cast<size_t>(M), 0.0);
  std::vector<double> pi_j(static_cast<size_t>(M), 0.0);
  std::vector<double> eta_k(static_cast<size_t>(K), 0.0);
  std::vector<double> gm_acc(static_cast<size_t>(std::max(1, n_u)), 0.0);

  const size_t mat_sz = need_hessian
                            ? (static_cast<size_t>(n_tilde) * static_cast<size_t>(n_tilde))
                            : 0u;
  std::vector<double> H_comp_v(mat_sz, 0.0);
  std::vector<double> H_miss_v(mat_sz, 0.0);
  std::vector<int> active_map(need_hessian ? static_cast<size_t>(n_tilde) : 0u, -1);
  std::vector<int> active_params;
  if (need_hessian) active_params.reserve(static_cast<size_t>(n_tilde));
  std::vector<double> S_j(
      need_hessian ? (static_cast<size_t>(n_tilde) * static_cast<size_t>(M)) : 0u, 0.0);
  std::vector<double> mu_vec(need_hessian ? static_cast<size_t>(n_tilde) : 0u, 0.0);
  std::vector<double> cbar_q(need_hessian ? static_cast<size_t>(M) : 0u, 0.0);
  std::vector<double> tail_qm(
      need_hessian ? (static_cast<size_t>(std::max(1, n_u)) * static_cast<size_t>(M)) : 0u, 0.0);
  std::vector<double> ctail_qm(
      need_hessian ? (static_cast<size_t>(std::max(1, n_u)) * static_cast<size_t>(M)) : 0u, 0.0);
  std::vector<int> loc_slots(static_cast<size_t>(n_loc_blocks), 0);
  std::vector<int> u_slots(static_cast<size_t>(std::max(1, n_u)), 0);

  auto get_active_slot = [&](int p) -> int {
    int a = active_map[static_cast<size_t>(p)];
    if (a < 0) {
      a = static_cast<int>(active_params.size());
      active_map[static_cast<size_t>(p)] = a;
      active_params.push_back(p);
      std::fill_n(&S_j[static_cast<size_t>(a) * static_cast<size_t>(M)], M, 0.0);
    }
    return a;
  };

  for (int j = 0; j < J; ++j) {
    const int p_start = person_ptr[j];
    const int p_end = person_ptr[j + 1];
    const int n_obs_j = p_end - p_start;

    std::copy(log_w.begin(), log_w.end(), log_L_j.begin());

    for (int l = 0; l < n_obs_j; ++l) {
      const int i = p_start + l;
      const double eta_i = eta_loc[static_cast<size_t>(i)];
      const int r_i = r_v[static_cast<size_t>(i)];
      const int x_i = x_v[static_cast<size_t>(i)];
      const double* g_r = &g_mat[static_cast<size_t>(r_i) * static_cast<size_t>(K)];

      for (int q = 0; q < M; ++q) {
        const double diff = nodes_v[static_cast<size_t>(q)] - eta_i;
        eta_k[0] = 0.0;
        double mx = 0.0;
        for (int c = 1; c < K; ++c) {
          const double ek = D * (static_cast<double>(c) * diff - g_r[static_cast<size_t>(c)]);
          eta_k[static_cast<size_t>(c)] = ek;
          if (ek > mx) mx = ek;
        }
        double denom = 0.0;
        double* p_iq = &p_buf[(static_cast<size_t>(l) * static_cast<size_t>(M) +
                               static_cast<size_t>(q)) * static_cast<size_t>(K)];
        for (int c = 0; c < K; ++c) {
          const double ex = std::exp(eta_k[static_cast<size_t>(c)] - mx);
          p_iq[static_cast<size_t>(c)] = ex;
          denom += ex;
        }
        const double log_denom = mx + std::log(denom);
        const double inv_denom = 1.0 / denom;
        for (int c = 0; c < K; ++c) {
          p_iq[static_cast<size_t>(c)] *= inv_denom;
        }
        log_L_j[static_cast<size_t>(q)] += (eta_k[static_cast<size_t>(x_i)] - log_denom);
      }
    }

    double row_max = log_L_j[0];
    for (int q = 1; q < M; ++q) {
      if (log_L_j[static_cast<size_t>(q)] > row_max) {
        row_max = log_L_j[static_cast<size_t>(q)];
      }
    }
    double row_sum = 0.0;
    for (int q = 0; q < M; ++q) {
      const double ex = std::exp(log_L_j[static_cast<size_t>(q)] - row_max);
      pi_j[static_cast<size_t>(q)] = ex;
      row_sum += ex;
    }
    ll_mml += row_max + std::log(row_sum);

    const double inv_row_sum = 1.0 / row_sum;
    double eap1 = 0.0;
    double eap2 = 0.0;
    for (int q = 0; q < M; ++q) {
      const double w = pi_j[static_cast<size_t>(q)] * inv_row_sum;
      pi_j[static_cast<size_t>(q)] = w;
      out_pi_jm(j, q) = w;
      const double th = nodes_v[static_cast<size_t>(q)];
      eap1 += w * th;
      eap2 += w * th * th;
    }
    out_eap[j] = eap1;
    out_psd[j] = std::sqrt(std::max(0.0, eap2 - eap1 * eap1));

    for (int l = 0; l < n_obs_j; ++l) {
      const int i = p_start + l;
      const int r_i = r_v[static_cast<size_t>(i)];
      const int x_i = x_v[static_cast<size_t>(i)];
      const double c_obs = static_cast<double>(x_i);
      const size_t base_l = static_cast<size_t>(i) * static_cast<size_t>(n_loc_blocks);
      const int u_base_i = step_rater ? (u_tilde_offset + r_i * n_u) : u_tilde_offset;

      if (!need_hessian) {
        double g_loc_i = 0.0;
        if (n_u > 0) {
          std::fill(gm_acc.begin(), gm_acc.begin() + n_u, 0.0);
        }
        for (int q = 0; q < M; ++q) {
          const double* p_iq = &p_buf[(static_cast<size_t>(l) * static_cast<size_t>(M) +
                                       static_cast<size_t>(q)) * static_cast<size_t>(K)];
          const double wq = pi_j[static_cast<size_t>(q)];
          double cbar = 0.0;
          for (int c = 1; c < K; ++c) {
            cbar += static_cast<double>(c) * p_iq[static_cast<size_t>(c)];
          }
          g_loc_i += wq * (-D * (c_obs - cbar));
          if (n_u > 0) {
            double acc_p = 0.0;
            for (int m = n_u - 1; m >= 0; --m) {
              const int c_add = m + 1;
              acc_p += p_iq[static_cast<size_t>(c_add)];
              const double ind = (x_i >= c_add && x_i <= (K - 2)) ? 1.0 : 0.0;
              gm_acc[static_cast<size_t>(m)] += wq * (-D * (ind - acc_p));
            }
          }
        }
        for (int b = 0; b < n_loc_blocks; ++b) {
          grad_tilde_v[static_cast<size_t>(loc_v[base_l + static_cast<size_t>(b)])] += g_loc_i;
        }
        for (int m = 0; m < n_u; ++m) {
          grad_tilde_v[static_cast<size_t>(u_base_i + m)] += gm_acc[static_cast<size_t>(m)];
        }
      } else {
        for (int b = 0; b < n_loc_blocks; ++b) {
          loc_slots[static_cast<size_t>(b)] = get_active_slot(loc_v[base_l + static_cast<size_t>(b)]);
        }
        for (int m = 0; m < n_u; ++m) {
          u_slots[static_cast<size_t>(m)] = get_active_slot(u_base_i + m);
        }

        double g_loc_i = 0.0;
        double info_bb_i = 0.0;
        for (int q = 0; q < M; ++q) {
          const double* p_iq = &p_buf[(static_cast<size_t>(l) * static_cast<size_t>(M) +
                                       static_cast<size_t>(q)) * static_cast<size_t>(K)];
          const double wq = pi_j[static_cast<size_t>(q)];
          double cbar = 0.0;
          double cbar2 = 0.0;
          for (int c = 1; c < K; ++c) {
            const double cd = static_cast<double>(c);
            const double pc = p_iq[static_cast<size_t>(c)];
            cbar += cd * pc;
            cbar2 += cd * cd * pc;
          }
          cbar_q[static_cast<size_t>(q)] = cbar;
          const double var_c = std::max(0.0, cbar2 - cbar * cbar);
          const double s_beta = -D * (c_obs - cbar);
          g_loc_i += wq * s_beta;
          info_bb_i += wq * var_c;

          for (int b = 0; b < n_loc_blocks; ++b) {
            S_j[static_cast<size_t>(loc_slots[static_cast<size_t>(b)]) * static_cast<size_t>(M) +
                static_cast<size_t>(q)] += s_beta;
          }

          if (n_u > 0) {
            double acc_p = 0.0;
            double acc_cp = 0.0;
            for (int m = n_u - 1; m >= 0; --m) {
              const int c_add = m + 1;
              const double pc = p_iq[static_cast<size_t>(c_add)];
              acc_p += pc;
              acc_cp += static_cast<double>(c_add) * pc;
              tail_qm[static_cast<size_t>(m) * static_cast<size_t>(M) + static_cast<size_t>(q)] = acc_p;
              ctail_qm[static_cast<size_t>(m) * static_cast<size_t>(M) + static_cast<size_t>(q)] = acc_cp;
              const double ind = (x_i >= c_add && x_i <= (K - 2)) ? 1.0 : 0.0;
              const double s_u = -D * (ind - acc_p);
              S_j[static_cast<size_t>(u_slots[static_cast<size_t>(m)]) * static_cast<size_t>(M) +
                  static_cast<size_t>(q)] += s_u;
            }
          }
        }
        info_bb_i *= D2;

        const size_t nt = static_cast<size_t>(n_tilde);
        for (int b1 = 0; b1 < n_loc_blocks; ++b1) {
          const size_t g1 = static_cast<size_t>(loc_v[base_l + static_cast<size_t>(b1)]);
          grad_tilde_v[g1] += g_loc_i;
          H_comp_v[g1 + g1 * nt] += info_bb_i;
          for (int b2 = b1 + 1; b2 < n_loc_blocks; ++b2) {
            const size_t g2 = static_cast<size_t>(loc_v[base_l + static_cast<size_t>(b2)]);
            H_comp_v[g1 + g2 * nt] += info_bb_i;
            H_comp_v[g2 + g1 * nt] += info_bb_i;
          }
        }

        for (int m = 0; m < n_u; ++m) {
          const size_t ui = static_cast<size_t>(u_base_i + m);
          const double* t_m = &tail_qm[static_cast<size_t>(m) * static_cast<size_t>(M)];
          const double* ct_m = &ctail_qm[static_cast<size_t>(m) * static_cast<size_t>(M)];
          const double ind_m = (x_i >= (m + 1) && x_i <= (K - 2)) ? 1.0 : 0.0;
          double gm_i = 0.0;
          double cov_bu = 0.0;
          for (int q = 0; q < M; ++q) {
            const double wq = pi_j[static_cast<size_t>(q)];
            gm_i += wq * (-D * (ind_m - t_m[static_cast<size_t>(q)]));
            cov_bu += wq * (ct_m[static_cast<size_t>(q)] -
                            cbar_q[static_cast<size_t>(q)] * t_m[static_cast<size_t>(q)]);
          }
          cov_bu *= D2;
          grad_tilde_v[ui] += gm_i;
          for (int b = 0; b < n_loc_blocks; ++b) {
            const size_t gb = static_cast<size_t>(loc_v[base_l + static_cast<size_t>(b)]);
            H_comp_v[gb + ui * nt] += cov_bu;
            H_comp_v[ui + gb * nt] += cov_bu;
          }
          for (int mp = m; mp < n_u; ++mp) {
            const size_t ui2 = static_cast<size_t>(u_base_i + mp);
            const double* t_mp = &tail_qm[static_cast<size_t>(mp) * static_cast<size_t>(M)];
            double cov_uu = 0.0;
            for (int q = 0; q < M; ++q) {
              const double tmp = t_mp[static_cast<size_t>(q)];
              cov_uu += pi_j[static_cast<size_t>(q)] * (tmp - t_m[static_cast<size_t>(q)] * tmp);
            }
            cov_uu *= D2;
            H_comp_v[ui + ui2 * nt] += cov_uu;
            if (mp != m) {
              H_comp_v[ui2 + ui * nt] += cov_uu;
            }
          }
        }
      }
    }

    if (need_hessian) {
      const size_t nt = static_cast<size_t>(n_tilde);
      const int K_j = static_cast<int>(active_params.size());
      for (int a = 0; a < K_j; ++a) {
        const double* Sa = &S_j[static_cast<size_t>(a) * static_cast<size_t>(M)];
        double mu = 0.0;
        for (int q = 0; q < M; ++q) {
          mu += pi_j[static_cast<size_t>(q)] * Sa[static_cast<size_t>(q)];
        }
        mu_vec[static_cast<size_t>(a)] = mu;
      }
      for (int a = 0; a < K_j; ++a) {
        const size_t pa = static_cast<size_t>(active_params[static_cast<size_t>(a)]);
        const double* Sa = &S_j[static_cast<size_t>(a) * static_cast<size_t>(M)];
        const double mua = mu_vec[static_cast<size_t>(a)];
        for (int b = a; b < K_j; ++b) {
          const size_t pb = static_cast<size_t>(active_params[static_cast<size_t>(b)]);
          const double* Sb = &S_j[static_cast<size_t>(b) * static_cast<size_t>(M)];
          double e2 = 0.0;
          for (int q = 0; q < M; ++q) {
            e2 += pi_j[static_cast<size_t>(q)] * Sa[static_cast<size_t>(q)] * Sb[static_cast<size_t>(q)];
          }
          const double cov_ab = e2 - mua * mu_vec[static_cast<size_t>(b)];
          H_miss_v[pa + pb * nt] += cov_ab;
          if (b != a) {
            H_miss_v[pb + pa * nt] += cov_ab;
          }
        }
      }
      for (int a = 0; a < K_j; ++a) {
        active_map[static_cast<size_t>(active_params[static_cast<size_t>(a)])] = -1;
      }
      active_params.clear();
    }
  }

  writable::doubles out_grad_tilde(n_tilde);
  for (int p = 0; p < n_tilde; ++p) {
    out_grad_tilde[p] = grad_tilde_v[static_cast<size_t>(p)];
  }

  if (!need_hessian) {
    return writable::list({
      "ll"_nm = ll_mml,
      "grad_tilde"_nm = out_grad_tilde,
      "pi_jm"_nm = out_pi_jm,
      "person_eap"_nm = out_eap,
      "person_psd"_nm = out_psd
    });
  }

  writable::doubles_matrix<> out_H_mml_tilde(n_tilde, n_tilde);
  const size_t nt = static_cast<size_t>(n_tilde);
  for (int c = 0; c < n_tilde; ++c) {
    const size_t cs = static_cast<size_t>(c);
    for (int r = 0; r < n_tilde; ++r) {
      const size_t rs = static_cast<size_t>(r);
      const double v_rc = H_comp_v[rs + cs * nt] - H_miss_v[rs + cs * nt];
      const double v_cr = H_comp_v[cs + rs * nt] - H_miss_v[cs + rs * nt];
      out_H_mml_tilde(r, c) = 0.5 * (v_rc + v_cr);
    }
  }

  return writable::list({
    "ll"_nm = ll_mml,
    "grad_tilde"_nm = out_grad_tilde,
    "H_mml_tilde"_nm = out_H_mml_tilde,
    "pi_jm"_nm = out_pi_jm,
    "person_eap"_nm = out_eap,
    "person_psd"_nm = out_psd
  });
}

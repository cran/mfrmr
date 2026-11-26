#ifndef MFRMR_FAST_ENGINE_H
#define MFRMR_FAST_ENGINE_H

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace mfrm_fast {

inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
inline constexpr double kInf = std::numeric_limits<double>::infinity();
inline constexpr double kPi = 3.141592653589793238462643383279502884;

// ============================================================================
// 1. Small Dense Reentrant Linear Algebra (Thread-safe, BLAS-independent)
// ============================================================================

inline bool chol_factor(const std::vector<double>& A, int n, std::vector<double>& L) {
  L.assign(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
  if (n == 0) return true;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum = A[static_cast<size_t>(i) * n + j];
      for (int k = 0; k < j; ++k) {
        sum -= L[static_cast<size_t>(i) * n + k] * L[static_cast<size_t>(j) * n + k];
      }
      if (i == j) {
        if (!std::isfinite(sum) || sum <= 1e-15) return false;
        L[static_cast<size_t>(i) * n + i] = std::sqrt(sum);
      } else {
        L[static_cast<size_t>(i) * n + j] = sum / L[static_cast<size_t>(j) * n + j];
      }
    }
  }
  return true;
}

inline void chol_solve_vec(const std::vector<double>& L, int n,
                           const double* b, double* x) {
  if (n == 0) return;
  for (int i = 0; i < n; ++i) {
    double sum = b[i];
    for (int k = 0; k < i; ++k) sum -= L[static_cast<size_t>(i) * n + k] * x[k];
    x[i] = sum / L[static_cast<size_t>(i) * n + i];
  }
  for (int i = n - 1; i >= 0; --i) {
    double sum = x[i];
    for (int k = i + 1; k < n; ++k) sum -= L[static_cast<size_t>(k) * n + i] * x[k];
    x[i] = sum / L[static_cast<size_t>(i) * n + i];
  }
}

inline bool chol_invert_and_logdet(const std::vector<double>& A, int n,
                                   std::vector<double>& invA, double& logdet) {
  invA.clear();
  logdet = 0.0;
  if (n == 0) return true;
  std::vector<double> inv_sd(static_cast<size_t>(n), 1.0);
  double sum_log_diag = 0.0;
  for (int i = 0; i < n; ++i) {
    double d = A[static_cast<size_t>(i) * n + i];
    if (!std::isfinite(d) || d <= 0.0) return false;
    inv_sd[static_cast<size_t>(i)] = 1.0 / std::sqrt(d);
    sum_log_diag += std::log(d);
  }
  std::vector<double> H(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double v = 0.5 * (A[static_cast<size_t>(i) * n + j] +
                        A[static_cast<size_t>(j) * n + i]) *
                 inv_sd[static_cast<size_t>(i)] * inv_sd[static_cast<size_t>(j)];
      H[static_cast<size_t>(i) * n + j] = v;
      H[static_cast<size_t>(j) * n + i] = v;
    }
  }
  std::vector<double> L;
  if (!chol_factor(H, n, L)) return false;
  invA.assign(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
  double logdet_H = 0.0;
  for (int i = 0; i < n; ++i) {
    logdet_H += 2.0 * std::log(L[static_cast<size_t>(i) * n + i]);
  }
  logdet = logdet_H + sum_log_diag;

  std::vector<double> e(static_cast<size_t>(n), 0.0);
  std::vector<double> col(static_cast<size_t>(n), 0.0);
  for (int j = 0; j < n; ++j) {
    std::fill(e.begin(), e.end(), 0.0);
    e[static_cast<size_t>(j)] = inv_sd[static_cast<size_t>(j)];
    chol_solve_vec(L, n, e.data(), col.data());
    for (int i = 0; i < n; ++i) {
      invA[static_cast<size_t>(i) * n + j] = col[static_cast<size_t>(i)] * inv_sd[static_cast<size_t>(i)];
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double sym = 0.5 * (invA[static_cast<size_t>(i) * n + j] +
                          invA[static_cast<size_t>(j) * n + i]);
      invA[static_cast<size_t>(i) * n + j] = sym;
      invA[static_cast<size_t>(j) * n + i] = sym;
    }
  }
  return true;
}

inline bool damped_spd_solve(const std::vector<double>& H, int n,
                             const std::vector<double>& g,
                             std::vector<double>& step,
                             double init_ridge = 0.0) {
  step.assign(static_cast<size_t>(n), 0.0);
  if (n == 0) return true;
  double max_diag = 1.0;
  for (int i = 0; i < n; ++i) {
    double d = std::abs(H[static_cast<size_t>(i) * n + i]);
    if (std::isfinite(d) && d > max_diag) max_diag = d;
  }
  std::vector<double> Hw = H;
  std::vector<double> L;
  double ridge = init_ridge;
  for (int attempt = 0; attempt < 12; ++attempt) {
    Hw = H;
    if (ridge > 0.0) {
      for (int i = 0; i < n; ++i) Hw[static_cast<size_t>(i) * n + i] += ridge;
    }
    if (chol_factor(Hw, n, L)) {
      chol_solve_vec(L, n, g.data(), step.data());
      bool ok = true;
      for (int i = 0; i < n; ++i) {
        if (!std::isfinite(step[static_cast<size_t>(i)])) { ok = false; break; }
      }
      if (ok) return true;
    }
    ridge = (ridge == 0.0) ? (1e-8 * max_diag) : (ridge * 10.0);
  }
  return false;
}

inline void gauss_hermite_normal(int Q, std::vector<double>& nodes,
                                 std::vector<double>& weights,
                                 std::vector<double>& log_weights) {
  nodes.assign(static_cast<size_t>(Q), 0.0);
  weights.assign(static_cast<size_t>(Q), 1.0);
  log_weights.assign(static_cast<size_t>(Q), 0.0);
  if (Q <= 1) return;

  std::vector<double> A(static_cast<size_t>(Q) * static_cast<size_t>(Q), 0.0);
  std::vector<double> V(static_cast<size_t>(Q) * static_cast<size_t>(Q), 0.0);
  for (int i = 0; i < Q; ++i) V[static_cast<size_t>(i) * Q + i] = 1.0;
  for (int i = 0; i < Q - 1; ++i) {
    double s = std::sqrt(static_cast<double>(i + 1));
    A[static_cast<size_t>(i) * Q + (i + 1)] = s;
    A[static_cast<size_t>(i + 1) * Q + i] = s;
  }

  for (int sweep = 0; sweep < 80; ++sweep) {
    double max_off = 0.0;
    for (int p = 0; p < Q - 1; ++p) {
      for (int q = p + 1; q < Q; ++q) {
        max_off = std::max(max_off, std::abs(A[static_cast<size_t>(p) * Q + q]));
      }
    }
    if (max_off < 1e-15) break;
    for (int p = 0; p < Q - 1; ++p) {
      for (int q = p + 1; q < Q; ++q) {
        double apq = A[static_cast<size_t>(p) * Q + q];
        if (std::abs(apq) < 1e-16) continue;
        double app = A[static_cast<size_t>(p) * Q + p];
        double aqq = A[static_cast<size_t>(q) * Q + q];
        double tau = (aqq - app) / (2.0 * apq);
        double t = (tau >= 0.0) ? (1.0 / (tau + std::sqrt(1.0 + tau * tau)))
                                : (-1.0 / (-tau + std::sqrt(1.0 + tau * tau)));
        double c = 1.0 / std::sqrt(1.0 + t * t);
        double s = t * c;
        A[static_cast<size_t>(p) * Q + p] = app - t * apq;
        A[static_cast<size_t>(q) * Q + q] = aqq + t * apq;
        A[static_cast<size_t>(p) * Q + q] = 0.0;
        A[static_cast<size_t>(q) * Q + p] = 0.0;
        for (int r = 0; r < Q; ++r) {
          if (r == p || r == q) continue;
          double arp = A[static_cast<size_t>(r) * Q + p];
          double arq = A[static_cast<size_t>(r) * Q + q];
          A[static_cast<size_t>(r) * Q + p] = c * arp - s * arq;
          A[static_cast<size_t>(p) * Q + r] = A[static_cast<size_t>(r) * Q + p];
          A[static_cast<size_t>(r) * Q + q] = s * arp + c * arq;
          A[static_cast<size_t>(q) * Q + r] = A[static_cast<size_t>(r) * Q + q];
        }
        for (int r = 0; r < Q; ++r) {
          double vrp = V[static_cast<size_t>(r) * Q + p];
          double vrq = V[static_cast<size_t>(r) * Q + q];
          V[static_cast<size_t>(r) * Q + p] = c * vrp - s * vrq;
          V[static_cast<size_t>(r) * Q + q] = s * vrp + c * vrq;
        }
      }
    }
  }

  std::vector<int> ord(static_cast<size_t>(Q));
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
    return A[static_cast<size_t>(a) * Q + a] < A[static_cast<size_t>(b) * Q + b];
  });
  double wsum = 0.0;
  for (int k = 0; k < Q; ++k) {
    int idx = ord[static_cast<size_t>(k)];
    nodes[static_cast<size_t>(k)] = A[static_cast<size_t>(idx) * Q + idx];
    double v0 = V[idx];
    weights[static_cast<size_t>(k)] = v0 * v0;
    wsum += weights[static_cast<size_t>(k)];
  }
  for (int k = 0; k < Q; ++k) {
    weights[static_cast<size_t>(k)] /= wsum;
    log_weights[static_cast<size_t>(k)] = std::log(weights[static_cast<size_t>(k)]);
  }
}

// ============================================================================
// 2. Category Probability & Moments Helper (Supports Per-Group Category Bounds)
// ============================================================================

struct CatMoments {
  double log_py = 0.0;
  double mean = 0.0;
  double var = 0.0;
  double m3 = 0.0;
  double m4 = 0.0;
};

// Note: step_cum[k] == +kInf marks a structurally unsupported category k!
inline CatMoments eval_cat_probs_and_moments(double eta,
                                             const double* step_cum,
                                             int C, int y_obs,
                                             double* probs_out,
                                             int k_min = 0,
                                             int k_max = -1) {
  if (k_max < 0) k_max = C - 1;
  for (int k = 0; k < k_min; ++k) probs_out[k] = 0.0;
  for (int k = k_max + 1; k < C; ++k) probs_out[k] = 0.0;

  double max_lp = -kInf;
  for (int k = k_min; k <= k_max; ++k) {
    double sck = step_cum[k];
    if (!std::isfinite(sck)) {
      probs_out[k] = -kInf;
    } else {
      double lp = static_cast<double>(k - k_min) * eta - sck;
      probs_out[k] = lp;
      if (lp > max_lp) max_lp = lp;
    }
  }
  if (!std::isfinite(max_lp)) {
    CatMoments bad;
    bad.log_py = -kInf;
    return bad;
  }
  double sum_exp = 0.0;
  for (int k = k_min; k <= k_max; ++k) {
    if (std::isfinite(probs_out[k])) {
      probs_out[k] = std::exp(probs_out[k] - max_lp);
      sum_exp += probs_out[k];
    } else {
      probs_out[k] = 0.0;
    }
  }
  double log_z = max_lp + std::log(sum_exp);
  double inv_sum = 1.0 / sum_exp;
  double mu = 0.0;
  for (int k = k_min; k <= k_max; ++k) {
    probs_out[k] *= inv_sum;
    mu += static_cast<double>(k) * probs_out[k];
  }
  double v = 0.0, m3 = 0.0, m4 = 0.0;
  for (int k = k_min; k <= k_max; ++k) {
    double d = static_cast<double>(k) - mu;
    double d2 = d * d;
    double pk = probs_out[k];
    v += d2 * pk;
    m3 += d2 * d * pk;
    m4 += d2 * d2 * pk;
  }
  CatMoments out;
  if (y_obs >= k_min && y_obs <= k_max && std::isfinite(step_cum[y_obs])) {
    out.log_py = static_cast<double>(y_obs - k_min) * eta - step_cum[y_obs] - log_z;
  } else {
    out.log_py = -kInf;
  }
  out.mean = mu;
  out.var = v;
  out.m3 = m3;
  out.m4 = m4;
  return out;
}

inline double expected_score_from_eta(double eta, const double* step_cum,
                                      int C, double rating_min,
                                      int k_min = 0, int k_max = -1) {
  if (!std::isfinite(eta)) return kNaN;
  if (k_max < 0) k_max = C - 1;
  double max_lp = -kInf;
  double lps[64];
  std::vector<double> heap_lp;
  double* lp = (C <= 64) ? lps : (heap_lp.resize(static_cast<size_t>(C)), heap_lp.data());
  for (int k = k_min; k <= k_max; ++k) {
    if (!std::isfinite(step_cum[k])) {
      lp[k] = -kInf;
    } else {
      lp[k] = static_cast<double>(k - k_min) * eta - step_cum[k];
      if (lp[k] > max_lp) max_lp = lp[k];
    }
  }
  if (!std::isfinite(max_lp)) return kNaN;
  double z = 0.0, num = 0.0;
  for (int k = k_min; k <= k_max; ++k) {
    if (std::isfinite(lp[k])) {
      double w = std::exp(lp[k] - max_lp);
      z += w;
      num += static_cast<double>(k) * w;
    }
  }
  return rating_min + num / z;
}

// Generalized GPCM / GMFRM category probabilities & moments:
//   slope_formula == 0 (Full-Predictor, Uto & Ueno 2020): u_k = alpha * ((k - k_min)*(theta + base_eta) - S_k)
//   slope_formula == 1 (Ability-Only, sirt::rm.facets):   u_k = (k - k_min)*(alpha * theta + base_eta) - S_k
inline CatMoments eval_gpcm_probs_and_moments(double theta,
                                              double base_eta,
                                              double alpha,
                                              int slope_formula,
                                              const double* step_cum,
                                              int C, int y_obs,
                                              double* probs_out,
                                              int k_min = 0,
                                              int k_max = -1) {
  if (alpha == 1.0) {
    return eval_cat_probs_and_moments(theta + base_eta, step_cum, C, y_obs, probs_out, k_min, k_max);
  }
  if (k_max < 0) k_max = C - 1;
  for (int k = 0; k < k_min; ++k) probs_out[k] = 0.0;
  for (int k = k_max + 1; k < C; ++k) probs_out[k] = 0.0;

  double max_lp = -kInf;
  const double lin = (slope_formula == 1) ? (alpha * theta + base_eta) : (alpha * (theta + base_eta));
  for (int k = k_min; k <= k_max; ++k) {
    double sck = step_cum[k];
    if (!std::isfinite(sck)) {
      probs_out[k] = -kInf;
    } else {
      double lp = static_cast<double>(k - k_min) * lin - ((slope_formula == 1) ? sck : (alpha * sck));
      probs_out[k] = lp;
      if (lp > max_lp) max_lp = lp;
    }
  }
  if (!std::isfinite(max_lp)) {
    CatMoments bad;
    bad.log_py = -kInf;
    return bad;
  }
  double sum_exp = 0.0;
  for (int k = k_min; k <= k_max; ++k) {
    if (std::isfinite(probs_out[k])) {
      probs_out[k] = std::exp(probs_out[k] - max_lp);
      sum_exp += probs_out[k];
    } else {
      probs_out[k] = 0.0;
    }
  }
  double log_z = max_lp + std::log(sum_exp);
  double inv_sum = 1.0 / sum_exp;
  double mu = 0.0;
  for (int k = k_min; k <= k_max; ++k) {
    probs_out[k] *= inv_sum;
    mu += static_cast<double>(k) * probs_out[k];
  }
  double v = 0.0, m3 = 0.0, m4 = 0.0;
  for (int k = k_min; k <= k_max; ++k) {
    double d = static_cast<double>(k) - mu;
    double d2 = d * d;
    double pk = probs_out[k];
    v += d2 * pk;
    m3 += d2 * d * pk;
    m4 += d2 * d2 * pk;
  }
  CatMoments out;
  if (y_obs >= k_min && y_obs <= k_max && std::isfinite(step_cum[y_obs])) {
    double sc_y = (slope_formula == 1) ? step_cum[y_obs] : (alpha * step_cum[y_obs]);
    out.log_py = static_cast<double>(y_obs - k_min) * lin - sc_y - log_z;
  } else {
    out.log_py = -kInf;
  }
  out.mean = mu;
  out.var = v;
  out.m3 = m3;
  out.m4 = m4;
  return out;
}

inline double expected_score_from_gpcm(double theta, double base_eta,
                                       double alpha, int slope_formula,
                                       const double* step_cum,
                                       int C, double rating_min,
                                       int k_min = 0, int k_max = -1) {
  if (alpha == 1.0) {
    return expected_score_from_eta(theta + base_eta, step_cum, C, rating_min, k_min, k_max);
  }
  if (!std::isfinite(theta) || !std::isfinite(base_eta) || !std::isfinite(alpha) || alpha <= 0.0) return kNaN;
  if (k_max < 0) k_max = C - 1;
  double max_lp = -kInf;
  double lps[64];
  std::vector<double> heap_lp;
  double* lp = (C <= 64) ? lps : (heap_lp.resize(static_cast<size_t>(C)), heap_lp.data());
  const double lin = (slope_formula == 1) ? (alpha * theta + base_eta) : (alpha * (theta + base_eta));
  for (int k = k_min; k <= k_max; ++k) {
    if (!std::isfinite(step_cum[k])) {
      lp[k] = -kInf;
    } else {
      lp[k] = static_cast<double>(k - k_min) * lin - ((slope_formula == 1) ? step_cum[k] : (alpha * step_cum[k]));
      if (lp[k] > max_lp) max_lp = lp[k];
    }
  }
  if (!std::isfinite(max_lp)) return kNaN;
  double z = 0.0, num = 0.0;
  for (int k = k_min; k <= k_max; ++k) {
    if (std::isfinite(lp[k])) {
      double w = std::exp(lp[k] - max_lp);
      z += w;
      num += static_cast<double>(k) * w;
    }
  }
  return rating_min + num / z;
}

// ============================================================================
// 3. Connectivity & Support Audit (Union-Find in O(N alpha(V)))
// ============================================================================

struct ConnectivityResult {
  bool connected = true;
  int n_subsets = 1;
  std::vector<int> person_subset;
  std::vector<std::vector<int>> facet_subset;
  std::vector<int> person_extreme_flag; // -1 min, 0 non-extreme, +1 max
  std::vector<double> person_obs_count;
  std::vector<double> person_raw_score;
  std::vector<double> person_min_score;
  std::vector<double> person_max_score;
  bool interior_categories_supported = true;
};

struct FastDesignSpec {
  int n_obs = 0;
  int n_person = 0;
  int n_facets = 0;
  int n_cat = 2; // C = K + 1 (categories 0 .. K)
  double rating_min = 0.0;
  int model_code = 0; // 0 = RSM, 1 = PCM
  int step_facet = -1; // -1 for RSM, 0..F-1 for PCM
  std::vector<int> facet_n_levels;
  std::vector<double> facet_signs; // -1 for severity/difficulty, +1 for positive
  std::vector<int> facet_dummy;    // 1 if dummy fixed at 0
  std::vector<int> facet_noncenter;// 1 if unconstrained
  int centering_scope = 0;         // 0 = remainder sum=0, 1 = facet total sum=0
  std::vector<std::vector<int>> anchor_fixed;   // [f][l] -> 0/1
  std::vector<std::vector<double>> anchor_val;  // [f][l]

  // Optional custom affine facet design: beta_f = custom_facet_d0[f] + custom_facet_B[f] * xi_f
  // Enables explanatory RSM/PCM (LLTM / rater covariates), group_anchors, and joint facet_interactions!
  std::vector<int> has_custom_facet_map;              // [F] 0/1
  std::vector<int> custom_facet_q;                    // [F] number of free columns
  std::vector<std::vector<double>> custom_facet_d0;   // [F][Lf]
  std::vector<std::vector<double>> custom_facet_B;    // [F][Lf * q_f] row-major
  std::vector<std::vector<int>> custom_facet_anchored;// [F][Lf] 1 if anchored/group-anchored to origin

  // Optional GPCM / GMFRM log-slope affine design (1-family, 2-family, or K-family):
  // log_alpha_s = slope_d0[s] + slope_B[s] * xi_slope_s,  alpha_i = prod_s exp(log_alpha_s[l_{i, f_s}])
  int n_slope_facets = 0;
  int slope_formula = 0; // 0 = full_predictor (Uto & Ueno 2020), 1 = ability_only (sirt::rm.facets)
  std::vector<int> slope_facets;             // [S] 0-based facet index in 0..F-1
  std::vector<int> slope_q;                  // [S] free columns for family s
  std::vector<std::vector<double>> slope_d0; // [S][L_{f_s}]
  std::vector<std::vector<double>> slope_B;  // [S][L_{f_s} * q_s] row-major
  std::vector<std::vector<int>> slope_fixed; // [S][L_{f_s}] 1 if fixed

  // Optional per-step-group structural category bounds [cat_min[g], cat_max[g]] in 0..K
  // and custom affine step map: tau = custom_step_d0 + custom_step_B * xi_step
  std::vector<int> step_cat_min;       // [n_step_groups], default 0
  std::vector<int> step_cat_max;       // [n_step_groups], default K
  bool has_custom_step_map = false;
  int custom_step_q = 0;
  std::vector<double> custom_step_d0;  // [n_step_groups * K]
  std::vector<double> custom_step_B;   // [(n_step_groups * K) * custom_step_q] row-major
  std::vector<int> custom_step_fixed;  // [n_step_groups * K] 1 if fixed/inactive (SE = NaN)

  // Observation arrays (0-based indices)
  std::vector<int> person;     // [n_obs]
  std::vector<int> score;      // [n_obs] in 0 .. n_cat - 1
  std::vector<double> weight;  // [n_obs]
  std::vector<int> facet_flat; // row-major [n_obs * n_facets + f]
  std::vector<int> step_group; // [n_obs] in 0 .. n_step_groups - 1
  int n_step_groups = 1;

  // Person grouping CSR for cache-friendly person loops
  std::vector<int> p_row_ptr;  // [n_person + 1]
  std::vector<int> p_obs_idx;  // [n_obs]

  // Latent regression design X_person [n_person * pop_dim] (row-major)
  int pop_mode = 0; // 0 = fixed N(0, sigma_fixed^2), 1 = free (mu/beta, sigma)
  int pop_dim = 0;
  double fixed_sigma = 1.0;
  std::vector<double> pop_X; // [n_person * pop_dim]

  void ensure_step_bounds() {
    int K = n_cat - 1;
    if (static_cast<int>(step_cat_min.size()) != n_step_groups) {
      step_cat_min.assign(static_cast<size_t>(n_step_groups), 0);
    }
    if (static_cast<int>(step_cat_max.size()) != n_step_groups) {
      step_cat_max.assign(static_cast<size_t>(n_step_groups), K);
    }
  }

  void build_person_csr() {
    ensure_step_bounds();
    p_row_ptr.assign(static_cast<size_t>(n_person) + 1, 0);
    for (int i = 0; i < n_obs; ++i) {
      p_row_ptr[static_cast<size_t>(person[static_cast<size_t>(i)]) + 1]++;
    }
    for (int n = 0; n < n_person; ++n) {
      p_row_ptr[static_cast<size_t>(n) + 1] += p_row_ptr[static_cast<size_t>(n)];
    }
    p_obs_idx.assign(static_cast<size_t>(n_obs), 0);
    std::vector<int> cur(p_row_ptr.begin(), p_row_ptr.begin() + n_person);
    for (int i = 0; i < n_obs; ++i) {
      int p = person[static_cast<size_t>(i)];
      p_obs_idx[static_cast<size_t>(cur[static_cast<size_t>(p)]++)] = i;
    }
  }
};

inline ConnectivityResult check_connectivity(const FastDesignSpec& spec) {
  ConnectivityResult res;
  const int P = spec.n_person;
  const int F = spec.n_facets;
  const int K = spec.n_cat - 1;

  std::vector<int> s_min = spec.step_cat_min;
  std::vector<int> s_max = spec.step_cat_max;
  if (static_cast<int>(s_min.size()) != spec.n_step_groups) s_min.assign(static_cast<size_t>(spec.n_step_groups), 0);
  if (static_cast<int>(s_max.size()) != spec.n_step_groups) s_max.assign(static_cast<size_t>(spec.n_step_groups), K);

  std::vector<int> f_offset(static_cast<size_t>(F), 0);
  int total_nodes = P;
  for (int f = 0; f < F; ++f) {
    f_offset[static_cast<size_t>(f)] = total_nodes;
    if (!spec.facet_dummy[static_cast<size_t>(f)]) {
      total_nodes += spec.facet_n_levels[static_cast<size_t>(f)];
    }
  }

  std::vector<int> parent(static_cast<size_t>(total_nodes));
  std::vector<int> rnk(static_cast<size_t>(total_nodes), 0);
  std::vector<int> active(static_cast<size_t>(total_nodes), 0);
  std::iota(parent.begin(), parent.end(), 0);

  auto find_root = [&](int x) {
    while (parent[static_cast<size_t>(x)] != x) {
      parent[static_cast<size_t>(x)] = parent[static_cast<size_t>(parent[static_cast<size_t>(x)])];
      x = parent[static_cast<size_t>(x)];
    }
    return x;
  };
  auto unite = [&](int a, int b) {
    int ra = find_root(a), rb = find_root(b);
    if (ra == rb) return;
    if (rnk[static_cast<size_t>(ra)] < rnk[static_cast<size_t>(rb)]) std::swap(ra, rb);
    parent[static_cast<size_t>(rb)] = ra;
    if (rnk[static_cast<size_t>(ra)] == rnk[static_cast<size_t>(rb)]) rnk[static_cast<size_t>(ra)]++;
  };

  int first_anchored = -1;
  for (int f = 0; f < F; ++f) {
    if (spec.facet_dummy[static_cast<size_t>(f)]) continue;
    if (!spec.has_custom_facet_map.empty() && spec.has_custom_facet_map[static_cast<size_t>(f)]) {
      if (!spec.custom_facet_anchored.empty() &&
          static_cast<int>(spec.custom_facet_anchored[static_cast<size_t>(f)].size()) == spec.facet_n_levels[static_cast<size_t>(f)]) {
        // Explicit anchor / group_anchor mask: tie only anchored/grouped elements to fixed origin
        for (int l = 0; l < spec.facet_n_levels[static_cast<size_t>(f)]; ++l) {
          if (spec.custom_facet_anchored[static_cast<size_t>(f)][static_cast<size_t>(l)] == 1) {
            int u = f_offset[static_cast<size_t>(f)] + l;
            if (first_anchored < 0) first_anchored = u;
            else unite(first_anchored, u);
          }
        }
      } else {
        // Explanatory facets (facet_formulas) tie levels through shared covariates
        for (int l = 1; l < spec.facet_n_levels[static_cast<size_t>(f)]; ++l) {
          unite(f_offset[static_cast<size_t>(f)], f_offset[static_cast<size_t>(f)] + l);
        }
      }
      continue;
    }
    for (int l = 0; l < spec.facet_n_levels[static_cast<size_t>(f)]; ++l) {
      if (spec.anchor_fixed[static_cast<size_t>(f)][static_cast<size_t>(l)]) {
        int u = f_offset[static_cast<size_t>(f)] + l;
        if (first_anchored < 0) first_anchored = u;
        else unite(first_anchored, u);
      }
    }
  }

  res.person_obs_count.assign(static_cast<size_t>(P), 0.0);
  res.person_raw_score.assign(static_cast<size_t>(P), 0.0);
  res.person_min_score.assign(static_cast<size_t>(P), 0.0);
  res.person_max_score.assign(static_cast<size_t>(P), 0.0);

  std::vector<int> cat_counts(static_cast<size_t>(spec.n_step_groups) * spec.n_cat, 0);

  for (int i = 0; i < spec.n_obs; ++i) {
    double w = spec.weight[static_cast<size_t>(i)];
    if (w <= 0.0) continue;
    int p = spec.person[static_cast<size_t>(i)];
    int y = spec.score[static_cast<size_t>(i)];
    int sg = spec.step_group[static_cast<size_t>(i)];
    active[static_cast<size_t>(p)] = 1;
    res.person_obs_count[static_cast<size_t>(p)] += w;
    res.person_raw_score[static_cast<size_t>(p)] += w * static_cast<double>(y);
    res.person_min_score[static_cast<size_t>(p)] += w * static_cast<double>(s_min[static_cast<size_t>(sg)]);
    res.person_max_score[static_cast<size_t>(p)] += w * static_cast<double>(s_max[static_cast<size_t>(sg)]);
    cat_counts[static_cast<size_t>(sg) * spec.n_cat + y]++;

    for (int f = 0; f < F; ++f) {
      if (spec.facet_dummy[static_cast<size_t>(f)]) continue;
      int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
      int u = f_offset[static_cast<size_t>(f)] + l;
      active[static_cast<size_t>(u)] = 1;
      unite(p, u);
    }
  }

  std::map<int, int> root_to_sub;
  res.person_subset.assign(static_cast<size_t>(P), 0);
  for (int p = 0; p < P; ++p) {
    if (!active[static_cast<size_t>(p)]) continue;
    int r = find_root(p);
    auto it = root_to_sub.find(r);
    if (it == root_to_sub.end()) {
      int id = static_cast<int>(root_to_sub.size()) + 1;
      root_to_sub[r] = id;
      res.person_subset[static_cast<size_t>(p)] = id;
    } else {
      res.person_subset[static_cast<size_t>(p)] = it->second;
    }
  }
  res.facet_subset.resize(static_cast<size_t>(F));
  for (int f = 0; f < F; ++f) {
    int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
    res.facet_subset[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), 0);
    if (spec.facet_dummy[static_cast<size_t>(f)]) continue;
    for (int l = 0; l < Lf; ++l) {
      int u = f_offset[static_cast<size_t>(f)] + l;
      if (!active[static_cast<size_t>(u)]) continue;
      int r = find_root(u);
      auto it = root_to_sub.find(r);
      int id = (it == root_to_sub.end()) ? (static_cast<int>(root_to_sub.size()) + 1) : it->second;
      root_to_sub[r] = id;
      res.facet_subset[static_cast<size_t>(f)][static_cast<size_t>(l)] = id;
    }
  }

  res.n_subsets = std::max(1, static_cast<int>(root_to_sub.size()));
  res.connected = (res.n_subsets == 1);

  res.person_extreme_flag.assign(static_cast<size_t>(P), 0);
  for (int p = 0; p < P; ++p) {
    if (res.person_obs_count[static_cast<size_t>(p)] <= 0.0) continue;
    if (res.person_raw_score[static_cast<size_t>(p)] <= res.person_min_score[static_cast<size_t>(p)] + 1e-12) {
      res.person_extreme_flag[static_cast<size_t>(p)] = -1;
    } else if (res.person_raw_score[static_cast<size_t>(p)] >= res.person_max_score[static_cast<size_t>(p)] - 1e-12) {
      res.person_extreme_flag[static_cast<size_t>(p)] = 1;
    }
  }

  res.interior_categories_supported = true;
  for (int g = 0; g < spec.n_step_groups; ++g) {
    int a = s_min[static_cast<size_t>(g)];
    int b = s_max[static_cast<size_t>(g)];
    for (int k = a + 1; k < b; ++k) {
      if (cat_counts[static_cast<size_t>(g) * spec.n_cat + k] == 0) {
        res.interior_categories_supported = false;
      }
    }
  }
  return res;
}

// ============================================================================
// 4. Unified Structural Affine Basis Mapping:
//    beta_f = facet_d0[f] + facet_B[f] * xi_f
//    tau    = step_d0     + step_B     * xi_step
// ============================================================================

struct StructuralMap {
  int n_xi = 0;
  int n_facet_xi = 0;
  std::vector<int> facet_xi_start;
  std::vector<int> facet_xi_len;
  std::vector<std::vector<double>> facet_d0; // [F][Lf]
  std::vector<std::vector<double>> facet_B;  // [F][Lf * len_f]
  std::vector<std::vector<int>> facet_fixed_mask; // [F][Lf] -> 1 if fixed/zero-var

  int step_xi_start = 0;
  int step_xi_len = 0;
  int n_step_groups = 1;
  int K = 1;
  std::vector<int> cat_min; // [G]
  std::vector<int> cat_max; // [G]
  bool is_general_step = false;
  std::vector<double> step_d0;   // [G * K] (NaN if transition outside [a_g, b_g-1])
  std::vector<double> step_B;    // [(G * K) * step_xi_len]
  std::vector<int> step_fixed;   // [G * K] 1 if fixed or inactive
  std::vector<double> cum_d0;    // [G * C] (+kInf if category unsupported)
  std::vector<double> cum_B;     // [G * C * step_xi_len]
  std::vector<std::vector<int>> group_step_cols; // [G] active local step col indices (0..step_xi_len-1)

  // GPCM / GMFRM log-slope affine coordinates (appended after step_xi inside xi)
  int slope_xi_start = 0;
  int slope_xi_len = 0;
  int n_slope_facets = 0;
  std::vector<int> slope_facets;
  std::vector<int> slope_xi_offset; // [S] offset inside xi
  std::vector<int> slope_q;         // [S]
  std::vector<std::vector<double>> slope_d0;
  std::vector<std::vector<double>> slope_B;
  std::vector<std::vector<int>> slope_fixed;

  int m_cat() const { return K + 1; }

  static StructuralMap build(const FastDesignSpec& spec) {
    StructuralMap m;
    const int F = spec.n_facets;
    const int C = spec.n_cat;
    m.K = C - 1;
    m.n_step_groups = spec.n_step_groups;
    m.cat_min = spec.step_cat_min.empty() ? std::vector<int>(static_cast<size_t>(m.n_step_groups), 0) : spec.step_cat_min;
    m.cat_max = spec.step_cat_max.empty() ? std::vector<int>(static_cast<size_t>(m.n_step_groups), m.K) : spec.step_cat_max;

    m.facet_xi_start.assign(static_cast<size_t>(F), 0);
    m.facet_xi_len.assign(static_cast<size_t>(F), 0);
    m.facet_d0.resize(static_cast<size_t>(F));
    m.facet_B.resize(static_cast<size_t>(F));
    m.facet_fixed_mask.resize(static_cast<size_t>(F));

    int cur = 0;
    for (int f = 0; f < F; ++f) {
      int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
      m.facet_xi_start[static_cast<size_t>(f)] = cur;
      m.facet_d0[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), 0.0);
      m.facet_fixed_mask[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), 0);

      if (spec.facet_dummy[static_cast<size_t>(f)]) {
        std::fill(m.facet_fixed_mask[static_cast<size_t>(f)].begin(),
                  m.facet_fixed_mask[static_cast<size_t>(f)].end(), 1);
        continue;
      }

      if (!spec.has_custom_facet_map.empty() && spec.has_custom_facet_map[static_cast<size_t>(f)]) {
        int qf = spec.custom_facet_q[static_cast<size_t>(f)];
        m.facet_xi_len[static_cast<size_t>(f)] = qf;
        m.facet_d0[static_cast<size_t>(f)] = spec.custom_facet_d0[static_cast<size_t>(f)];
        m.facet_B[static_cast<size_t>(f)] = spec.custom_facet_B[static_cast<size_t>(f)];
        for (int l = 0; l < Lf; ++l) {
          bool any_nz = false;
          for (int c = 0; c < qf; ++c) {
            if (std::abs(m.facet_B[static_cast<size_t>(f)][static_cast<size_t>(l) * qf + c]) > 1e-15) {
              any_nz = true; break;
            }
          }
          if (!any_nz) m.facet_fixed_mask[static_cast<size_t>(f)][static_cast<size_t>(l)] = 1;
        }
        cur += qf;
        continue;
      }

      double anchored_sum = 0.0;
      std::vector<int> unanc;
      for (int l = 0; l < Lf; ++l) {
        if (spec.anchor_fixed[static_cast<size_t>(f)][static_cast<size_t>(l)]) {
          double av = spec.anchor_val[static_cast<size_t>(f)][static_cast<size_t>(l)];
          m.facet_d0[static_cast<size_t>(f)][static_cast<size_t>(l)] = av;
          m.facet_fixed_mask[static_cast<size_t>(f)][static_cast<size_t>(l)] = 1;
          anchored_sum += av;
        } else {
          unanc.push_back(l);
        }
      }
      double target = (spec.centering_scope == 1) ? (-anchored_sum) : 0.0;
      if (unanc.empty()) {
        m.facet_xi_len[static_cast<size_t>(f)] = 0;
      } else if (spec.facet_noncenter[static_cast<size_t>(f)]) {
        int qf = static_cast<int>(unanc.size());
        m.facet_xi_len[static_cast<size_t>(f)] = qf;
        m.facet_B[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf) * qf, 0.0);
        for (int c = 0; c < qf; ++c) {
          m.facet_B[static_cast<size_t>(f)][static_cast<size_t>(unanc[static_cast<size_t>(c)]) * qf + c] = 1.0;
        }
        cur += qf;
      } else {
        int last_l = unanc.back();
        int qf = static_cast<int>(unanc.size()) - 1;
        m.facet_xi_len[static_cast<size_t>(f)] = qf;
        m.facet_d0[static_cast<size_t>(f)][static_cast<size_t>(last_l)] = target;
        if (qf == 0) {
          m.facet_fixed_mask[static_cast<size_t>(f)][static_cast<size_t>(last_l)] = 1;
        } else {
          m.facet_B[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf) * qf, 0.0);
          for (int c = 0; c < qf; ++c) {
            m.facet_B[static_cast<size_t>(f)][static_cast<size_t>(unanc[static_cast<size_t>(c)]) * qf + c] = 1.0;
            m.facet_B[static_cast<size_t>(f)][static_cast<size_t>(last_l) * qf + c] = -1.0;
          }
        }
        cur += qf;
      }
    }
    m.n_facet_xi = cur;
    m.step_xi_start = cur;

    const int G = m.n_step_groups;
    const int K = m.K;
    bool custom_bounds = false;
    for (int g = 0; g < G; ++g) {
      if (m.cat_min[static_cast<size_t>(g)] != 0 || m.cat_max[static_cast<size_t>(g)] != K) {
        custom_bounds = true;
        break;
      }
    }
    m.is_general_step = spec.has_custom_step_map || custom_bounds;

    if (spec.has_custom_step_map) {
      m.step_xi_len = spec.custom_step_q;
      m.step_d0 = spec.custom_step_d0;
      m.step_B = spec.custom_step_B;
      m.step_fixed = spec.custom_step_fixed;
    } else {
      // Standard sum-to-zero within each step group's active range [a_g, b_g]
      int q_tot = 0;
      for (int g = 0; g < G; ++g) {
        int Kg = m.cat_max[static_cast<size_t>(g)] - m.cat_min[static_cast<size_t>(g)];
        if (Kg >= 2) q_tot += (Kg - 1);
      }
      m.step_xi_len = q_tot;
      m.step_d0.assign(static_cast<size_t>(G) * K, kNaN);
      m.step_B.assign(static_cast<size_t>(G) * K * std::max(1, q_tot), 0.0);
      m.step_fixed.assign(static_cast<size_t>(G) * K, 1);
      int col_cursor = 0;
      for (int g = 0; g < G; ++g) {
        int ag = m.cat_min[static_cast<size_t>(g)];
        int bg = m.cat_max[static_cast<size_t>(g)];
        int Kg = bg - ag;
        for (int t = ag; t < bg; ++t) {
          m.step_d0[static_cast<size_t>(g) * K + t] = 0.0;
        }
        if (Kg >= 2) {
          for (int t = 0; t < Kg - 1; ++t) {
            int row_t = g * K + (ag + t);
            int row_last = g * K + (bg - 1);
            m.step_B[static_cast<size_t>(row_t) * q_tot + (col_cursor + t)] = 1.0;
            m.step_B[static_cast<size_t>(row_last) * q_tot + (col_cursor + t)] = -1.0;
            m.step_fixed[static_cast<size_t>(row_t)] = 0;
          }
          m.step_fixed[static_cast<size_t>(g) * K + (bg - 1)] = 0;
          col_cursor += (Kg - 1);
        }
      }
    }

    cur += m.step_xi_len;
    m.slope_xi_start = cur;
    m.n_slope_facets = spec.n_slope_facets;
    if (m.n_slope_facets > 0) {
      m.slope_facets = spec.slope_facets;
      m.slope_q = spec.slope_q;
      m.slope_d0 = spec.slope_d0;
      m.slope_B = spec.slope_B;
      m.slope_fixed = spec.slope_fixed;
      m.slope_xi_offset.assign(static_cast<size_t>(m.n_slope_facets), 0);
      for (int s = 0; s < m.n_slope_facets; ++s) {
        m.slope_xi_offset[static_cast<size_t>(s)] = cur;
        cur += m.slope_q[static_cast<size_t>(s)];
      }
      m.slope_xi_len = cur - m.slope_xi_start;
    }
    m.n_xi = cur;

    // Precompute cumulative offset S_{0,g,k} and cumulative basis A_cum[g, k, c]
    const int qs = m.step_xi_len;
    m.cum_d0.assign(static_cast<size_t>(G) * C, kInf);
    m.cum_B.assign(static_cast<size_t>(G) * C * std::max(1, qs), 0.0);
    m.group_step_cols.resize(static_cast<size_t>(G));

    for (int g = 0; g < G; ++g) {
      int ag = m.cat_min[static_cast<size_t>(g)];
      int bg = m.cat_max[static_cast<size_t>(g)];
      m.cum_d0[static_cast<size_t>(g) * C + ag] = 0.0;
      double run_d0 = 0.0;
      for (int k = ag + 1; k <= bg; ++k) {
        double d_step = m.step_d0[static_cast<size_t>(g) * K + (k - 1)];
        run_d0 += std::isfinite(d_step) ? d_step : 0.0;
        m.cum_d0[static_cast<size_t>(g) * C + k] = run_d0;
      }
      if (qs > 0) {
        for (int c = 0; c < qs; ++c) {
          double run_b = 0.0;
          bool any_nz = false;
          for (int k = ag + 1; k <= bg; ++k) {
            run_b += m.step_B[static_cast<size_t>(g * K + (k - 1)) * qs + c];
            if (std::abs(run_b) < 1e-15) run_b = 0.0;
            m.cum_B[(static_cast<size_t>(g) * C + k) * qs + c] = run_b;
            if (run_b != 0.0) any_nz = true;
          }
          if (any_nz) m.group_step_cols[static_cast<size_t>(g)].push_back(c);
        }
      }
    }
    return m;
  }

  void expand(const FastDesignSpec& spec, const double* xi,
              std::vector<std::vector<double>>& beta,
              std::vector<double>& tau_mat,     // [G * K]
              std::vector<double>& step_cum_mat // [G * C]
  ) const {
    const int F = spec.n_facets;
    const int C = m_cat();
    beta.resize(static_cast<size_t>(F));
    for (int f = 0; f < F; ++f) {
      int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
      beta[static_cast<size_t>(f)] = facet_d0[static_cast<size_t>(f)];
      int qf = facet_xi_len[static_cast<size_t>(f)];
      if (qf > 0) {
        int st = facet_xi_start[static_cast<size_t>(f)];
        const auto& Bf = facet_B[static_cast<size_t>(f)];
        for (int l = 0; l < Lf; ++l) {
          double s = beta[static_cast<size_t>(f)][static_cast<size_t>(l)];
          const double* row = Bf.data() + static_cast<size_t>(l) * qf;
          for (int c = 0; c < qf; ++c) s += row[c] * xi[st + c];
          beta[static_cast<size_t>(f)][static_cast<size_t>(l)] = s;
        }
      }
    }

    tau_mat = step_d0;
    step_cum_mat = cum_d0;
    if (step_xi_len > 0) {
      const double* xi_s = xi + step_xi_start;
      for (int g = 0; g < n_step_groups; ++g) {
        int ag = cat_min[static_cast<size_t>(g)];
        int bg = cat_max[static_cast<size_t>(g)];
        const auto& cols = group_step_cols[static_cast<size_t>(g)];
        for (int t = ag; t < bg; ++t) {
          int row_idx = g * K + t;
          if (!std::isfinite(tau_mat[static_cast<size_t>(row_idx)])) continue;
          double s = tau_mat[static_cast<size_t>(row_idx)];
          for (int c : cols) {
            s += step_B[static_cast<size_t>(row_idx) * step_xi_len + c] * xi_s[c];
          }
          tau_mat[static_cast<size_t>(row_idx)] = s;
        }
        for (int k = ag + 1; k <= bg; ++k) {
          double sc = step_cum_mat[static_cast<size_t>(g) * C + k];
          const double* cb = cum_B.data() + (static_cast<size_t>(g) * C + k) * step_xi_len;
          for (int c : cols) sc += cb[c] * xi_s[c];
          step_cum_mat[static_cast<size_t>(g) * C + k] = sc;
        }
      }
    }
  }

  // Sparse derivative terms for d(s_f * beta_{f, l}) / d(xi)
  struct LevelDeriv {
    std::vector<int> cols;
    std::vector<double> vals;
  };

  std::vector<std::vector<LevelDeriv>> build_level_derivs(const FastDesignSpec& spec) const {
    const int F = spec.n_facets;
    std::vector<std::vector<LevelDeriv>> out(static_cast<size_t>(F));
    for (int f = 0; f < F; ++f) {
      int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
      out[static_cast<size_t>(f)].resize(static_cast<size_t>(Lf));
      int qf = facet_xi_len[static_cast<size_t>(f)];
      if (qf <= 0 || spec.facet_dummy[static_cast<size_t>(f)]) continue;
      double s_f = spec.facet_signs[static_cast<size_t>(f)];
      int st = facet_xi_start[static_cast<size_t>(f)];
      const auto& Bf = facet_B[static_cast<size_t>(f)];
      for (int l = 0; l < Lf; ++l) {
        auto& ld = out[static_cast<size_t>(f)][static_cast<size_t>(l)];
        const double* row = Bf.data() + static_cast<size_t>(l) * qf;
        for (int c = 0; c < qf; ++c) {
          double v = s_f * row[c];
          if (std::abs(v) > 1e-15) {
            ld.cols.push_back(st + c);
            ld.vals.push_back(v);
          }
        }
      }
    }
    return out;
  }

  void propagate_covariance(const FastDesignSpec& spec,
                            const std::vector<double>& cov_xi, // [n_xi * n_xi]
                            std::vector<std::vector<double>>& facet_se,
                            std::vector<std::vector<double>>& facet_cov,
                            std::vector<double>& step_se,
                            std::vector<double>* step_cov_out = nullptr) const {
    const int F = spec.n_facets;
    facet_se.resize(static_cast<size_t>(F));
    facet_cov.resize(static_cast<size_t>(F));
    bool has_cov = (static_cast<int>(cov_xi.size()) == n_xi * n_xi && n_xi > 0);
    if (has_cov) {
      for (double v : cov_xi) {
        if (!std::isfinite(v)) { has_cov = false; break; }
      }
    }

    for (int f = 0; f < F; ++f) {
      int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
      facet_se[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), kNaN);
      facet_cov[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf) * Lf, kNaN);
      if (spec.facet_dummy[static_cast<size_t>(f)]) {
        std::fill(facet_se[static_cast<size_t>(f)].begin(), facet_se[static_cast<size_t>(f)].end(), 0.0);
        std::fill(facet_cov[static_cast<size_t>(f)].begin(), facet_cov[static_cast<size_t>(f)].end(), 0.0);
        continue;
      }
      int qf = facet_xi_len[static_cast<size_t>(f)];
      if (qf == 0) {
        std::fill(facet_se[static_cast<size_t>(f)].begin(), facet_se[static_cast<size_t>(f)].end(), 0.0);
        std::fill(facet_cov[static_cast<size_t>(f)].begin(), facet_cov[static_cast<size_t>(f)].end(), 0.0);
        continue;
      }
      if (!has_cov) continue;

      int st = facet_xi_start[static_cast<size_t>(f)];
      const auto& Bf = facet_B[static_cast<size_t>(f)];
      std::fill(facet_cov[static_cast<size_t>(f)].begin(), facet_cov[static_cast<size_t>(f)].end(), 0.0);
      for (int l1 = 0; l1 < Lf; ++l1) {
        const double* r1 = Bf.data() + static_cast<size_t>(l1) * qf;
        for (int l2 = 0; l2 <= l1; ++l2) {
          const double* r2 = Bf.data() + static_cast<size_t>(l2) * qf;
          double s = 0.0;
          for (int a = 0; a < qf; ++a) {
            if (r1[a] == 0.0) continue;
            for (int b = 0; b < qf; ++b) {
              if (r2[b] != 0.0) s += r1[a] * cov_xi[static_cast<size_t>(st + a) * n_xi + (st + b)] * r2[b];
            }
          }
          facet_cov[static_cast<size_t>(f)][static_cast<size_t>(l1) * Lf + l2] = s;
          facet_cov[static_cast<size_t>(f)][static_cast<size_t>(l2) * Lf + l1] = s;
        }
      }
      for (int l = 0; l < Lf; ++l) {
        if (facet_fixed_mask[static_cast<size_t>(f)][static_cast<size_t>(l)]) {
          facet_se[static_cast<size_t>(f)][static_cast<size_t>(l)] = 0.0;
        } else {
          double v = facet_cov[static_cast<size_t>(f)][static_cast<size_t>(l) * Lf + l];
          facet_se[static_cast<size_t>(f)][static_cast<size_t>(l)] = (v >= 0.0) ? std::sqrt(v) : kNaN;
        }
      }
    }

    const int n_st_rows = n_step_groups * K;
    step_se.assign(static_cast<size_t>(n_st_rows), kNaN);
    if (step_cov_out) step_cov_out->assign(static_cast<size_t>(n_st_rows) * n_st_rows, kNaN);
    if (step_xi_len == 0) {
      for (int r = 0; r < n_st_rows; ++r) {
        if (!spec.has_custom_step_map && std::isfinite(step_d0[static_cast<size_t>(r)])) {
          step_se[static_cast<size_t>(r)] = 0.0;
        }
      }
      if (step_cov_out) std::fill(step_cov_out->begin(), step_cov_out->end(), 0.0);
    } else if (has_cov) {
      int st = step_xi_start;
      int qs = step_xi_len;
      if (step_cov_out) {
        std::fill(step_cov_out->begin(), step_cov_out->end(), 0.0);
        for (int r1 = 0; r1 < n_st_rows; ++r1) {
          if (step_fixed[static_cast<size_t>(r1)]) continue;
          const double* b1 = step_B.data() + static_cast<size_t>(r1) * qs;
          for (int r2 = 0; r2 <= r1; ++r2) {
            if (step_fixed[static_cast<size_t>(r2)]) continue;
            const double* b2 = step_B.data() + static_cast<size_t>(r2) * qs;
            double s = 0.0;
            for (int a = 0; a < qs; ++a) {
              if (b1[a] == 0.0) continue;
              for (int b = 0; b < qs; ++b) {
                if (b2[b] != 0.0) s += b1[a] * cov_xi[static_cast<size_t>(st + a) * n_xi + (st + b)] * b2[b];
              }
            }
            (*step_cov_out)[static_cast<size_t>(r1) * n_st_rows + r2] = s;
            (*step_cov_out)[static_cast<size_t>(r2) * n_st_rows + r1] = s;
          }
        }
      }
      for (int r = 0; r < n_st_rows; ++r) {
        if (step_fixed[static_cast<size_t>(r)]) {
          step_se[static_cast<size_t>(r)] = kNaN; // Fixed or structurally inactive step: SE is NA
        } else {
          const double* br = step_B.data() + static_cast<size_t>(r) * qs;
          double v = 0.0;
          for (int a = 0; a < qs; ++a) {
            if (br[a] == 0.0) continue;
            for (int b = 0; b < qs; ++b) {
              if (br[b] != 0.0) v += br[a] * cov_xi[static_cast<size_t>(st + a) * n_xi + (st + b)] * br[b];
            }
          }
          step_se[static_cast<size_t>(r)] = (v >= 0.0) ? std::sqrt(v) : kNaN;
        }
      }
    }
  }

  // Expand GPCM log-slopes, slopes, and observation-level effective slope alpha_obs[0..n_obs-1]
  // Returns false if any log-slope overflows [-12, +12].
  bool expand_slopes(const FastDesignSpec& spec, const double* xi,
                     std::vector<std::vector<double>>& log_slope_est,
                     std::vector<std::vector<double>>& slope_est,
                     std::vector<double>& alpha_obs) const {
    if (n_slope_facets <= 0) {
      log_slope_est.clear();
      slope_est.clear();
      alpha_obs.assign(static_cast<size_t>(spec.n_obs), 1.0);
      return true;
    }
    const int F = spec.n_facets;
    log_slope_est.resize(static_cast<size_t>(n_slope_facets));
    slope_est.resize(static_cast<size_t>(n_slope_facets));
    for (int s = 0; s < n_slope_facets; ++s) {
      int f = slope_facets[static_cast<size_t>(s)];
      int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
      int qs = slope_q[static_cast<size_t>(s)];
      int st = slope_xi_offset[static_cast<size_t>(s)];
      log_slope_est[static_cast<size_t>(s)] = slope_d0[static_cast<size_t>(s)];
      slope_est[static_cast<size_t>(s)].assign(static_cast<size_t>(Lf), 1.0);
      const auto& Bs = slope_B[static_cast<size_t>(s)];
      for (int l = 0; l < Lf; ++l) {
        double g = log_slope_est[static_cast<size_t>(s)][static_cast<size_t>(l)];
        if (qs > 0) {
          const double* row = Bs.data() + static_cast<size_t>(l) * qs;
          for (int c = 0; c < qs; ++c) g += row[c] * xi[st + c];
        }
        if (!std::isfinite(g) || std::abs(g) > 12.0) return false;
        log_slope_est[static_cast<size_t>(s)][static_cast<size_t>(l)] = g;
        slope_est[static_cast<size_t>(s)][static_cast<size_t>(l)] = std::exp(g);
      }
    }
    alpha_obs.assign(static_cast<size_t>(spec.n_obs), 1.0);
    for (int i = 0; i < spec.n_obs; ++i) {
      double a = 1.0;
      for (int s = 0; s < n_slope_facets; ++s) {
        int f = slope_facets[static_cast<size_t>(s)];
        int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
        a *= slope_est[static_cast<size_t>(s)][static_cast<size_t>(l)];
      }
      if (!std::isfinite(a) || a <= 1e-8 || a > 1e6) return false;
      alpha_obs[static_cast<size_t>(i)] = a;
    }
    return true;
  }

  std::vector<std::vector<LevelDeriv>> build_slope_derivs(const FastDesignSpec& spec) const {
    std::vector<std::vector<LevelDeriv>> out(static_cast<size_t>(n_slope_facets));
    for (int s = 0; s < n_slope_facets; ++s) {
      int f = slope_facets[static_cast<size_t>(s)];
      int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
      out[static_cast<size_t>(s)].resize(static_cast<size_t>(Lf));
      int qs = slope_q[static_cast<size_t>(s)];
      if (qs <= 0) continue;
      int st = slope_xi_offset[static_cast<size_t>(s)];
      const auto& Bs = slope_B[static_cast<size_t>(s)];
      for (int l = 0; l < Lf; ++l) {
        auto& ld = out[static_cast<size_t>(s)][static_cast<size_t>(l)];
        const double* row = Bs.data() + static_cast<size_t>(l) * qs;
        for (int c = 0; c < qs; ++c) {
          if (std::abs(row[c]) > 1e-15) {
            ld.cols.push_back(st + c);
            ld.vals.push_back(row[c]);
          }
        }
      }
    }
    return out;
  }

  void propagate_slope_covariance(const FastDesignSpec& spec,
                                  const std::vector<double>& cov_xi,
                                  const std::vector<std::vector<double>>& slope_est,
                                  std::vector<std::vector<double>>& slope_se,
                                  std::vector<std::vector<double>>& slope_cov) const {
    slope_se.resize(static_cast<size_t>(n_slope_facets));
    slope_cov.resize(static_cast<size_t>(n_slope_facets));
    if (n_slope_facets <= 0) return;
    bool has_cov = (static_cast<int>(cov_xi.size()) == n_xi * n_xi && n_xi > 0);
    if (has_cov) {
      for (double v : cov_xi) {
        if (!std::isfinite(v)) { has_cov = false; break; }
      }
    }
    for (int s = 0; s < n_slope_facets; ++s) {
      int f = slope_facets[static_cast<size_t>(s)];
      int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
      slope_se[static_cast<size_t>(s)].assign(static_cast<size_t>(Lf), kNaN);
      slope_cov[static_cast<size_t>(s)].assign(static_cast<size_t>(Lf) * Lf, kNaN);
      int qs = slope_q[static_cast<size_t>(s)];
      if (qs == 0) {
        std::fill(slope_se[static_cast<size_t>(s)].begin(), slope_se[static_cast<size_t>(s)].end(), 0.0);
        std::fill(slope_cov[static_cast<size_t>(s)].begin(), slope_cov[static_cast<size_t>(s)].end(), 0.0);
        continue;
      }
      if (!has_cov) continue;
      int st = slope_xi_offset[static_cast<size_t>(s)];
      const auto& Bs = slope_B[static_cast<size_t>(s)];
      const auto& al = slope_est[static_cast<size_t>(s)];
      std::fill(slope_cov[static_cast<size_t>(s)].begin(), slope_cov[static_cast<size_t>(s)].end(), 0.0);
      for (int l1 = 0; l1 < Lf; ++l1) {
        const double* r1 = Bs.data() + static_cast<size_t>(l1) * qs;
        for (int l2 = 0; l2 <= l1; ++l2) {
          const double* r2 = Bs.data() + static_cast<size_t>(l2) * qs;
          double cov_log = 0.0;
          for (int a = 0; a < qs; ++a) {
            if (r1[a] == 0.0) continue;
            for (int b = 0; b < qs; ++b) {
              if (r2[b] != 0.0) cov_log += r1[a] * cov_xi[static_cast<size_t>(st + a) * n_xi + (st + b)] * r2[b];
            }
          }
          double c_alpha = al[static_cast<size_t>(l1)] * al[static_cast<size_t>(l2)] * cov_log;
          slope_cov[static_cast<size_t>(s)][static_cast<size_t>(l1) * Lf + l2] = c_alpha;
          slope_cov[static_cast<size_t>(s)][static_cast<size_t>(l2) * Lf + l1] = c_alpha;
        }
      }
      for (int l = 0; l < Lf; ++l) {
        if (!slope_fixed.empty() && slope_fixed[static_cast<size_t>(s)][static_cast<size_t>(l)]) {
          slope_se[static_cast<size_t>(s)][static_cast<size_t>(l)] = 0.0;
        } else {
          double v = slope_cov[static_cast<size_t>(s)][static_cast<size_t>(l) * Lf + l];
          slope_se[static_cast<size_t>(s)][static_cast<size_t>(l)] = (v >= 0.0) ? std::sqrt(v) : kNaN;
        }
      }
    }
  }
};

inline int obs_log_slope_xi_terms(int i, const FastDesignSpec& spec,
                                  const std::vector<std::vector<StructuralMap::LevelDeriv>>& sderivs,
                                  int* out_cols, double* out_vals) {
  int nnz = 0;
  const int F = spec.n_facets;
  for (int s = 0; s < spec.n_slope_facets; ++s) {
    int f = spec.slope_facets[static_cast<size_t>(s)];
    int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
    const auto& ld = sderivs[static_cast<size_t>(s)][static_cast<size_t>(l)];
    for (size_t k = 0; k < ld.cols.size(); ++k) {
      out_cols[nnz] = ld.cols[k];
      out_vals[nnz] = ld.vals[k];
      nnz++;
    }
  }
  return nnz;
}

inline int obs_eta_xi_terms(int i, const FastDesignSpec& spec,
                            const std::vector<std::vector<StructuralMap::LevelDeriv>>& lderivs,
                            int* out_cols, double* out_vals) {
  int nnz = 0;
  const int F = spec.n_facets;
  for (int f = 0; f < F; ++f) {
    int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
    const auto& ld = lderivs[static_cast<size_t>(f)][static_cast<size_t>(l)];
    for (size_t k = 0; k < ld.cols.size(); ++k) {
      out_cols[nnz] = ld.cols[k];
      out_vals[nnz] = ld.vals[k];
      nnz++;
    }
  }
  return nnz;
}

// ============================================================================
// 5. Polytomous Warm (1989) WLE Scorer (Supports Per-Group Category Bounds)
// ============================================================================

struct WleResult {
  std::vector<double> theta;
  std::vector<double> se;
};

inline WleResult score_wle_persons(const FastDesignSpec& spec,
                                   const std::vector<double>& base_eta_obs,
                                   const std::vector<double>& step_cum_mat,
                                   const std::vector<double>* alpha_obs_ptr = nullptr) {
  const int P = spec.n_person;
  const int C = spec.n_cat;
  const int K = C - 1;
  const int sform = spec.slope_formula;
  std::vector<int> s_min = spec.step_cat_min.empty() ? std::vector<int>(static_cast<size_t>(spec.n_step_groups), 0) : spec.step_cat_min;
  std::vector<int> s_max = spec.step_cat_max.empty() ? std::vector<int>(static_cast<size_t>(spec.n_step_groups), K) : spec.step_cat_max;

  WleResult res;
  res.theta.assign(static_cast<size_t>(P), kNaN);
  res.se.assign(static_cast<size_t>(P), kNaN);

  std::vector<double> probs(static_cast<size_t>(C));
  for (int p = 0; p < P; ++p) {
    int beg = spec.p_row_ptr[static_cast<size_t>(p)];
    int end = spec.p_row_ptr[static_cast<size_t>(p) + 1];
    if (beg == end) continue;

    double raw = 0.0, min_raw = 0.0, max_raw = 0.0, mean_off = 0.0, w_tot = 0.0;
    for (int idx = beg; idx < end; ++idx) {
      int i = spec.p_obs_idx[static_cast<size_t>(idx)];
      double w = spec.weight[static_cast<size_t>(i)];
      if (w <= 0.0) continue;
      int sg = spec.step_group[static_cast<size_t>(i)];
      raw += w * spec.score[static_cast<size_t>(i)];
      min_raw += w * s_min[static_cast<size_t>(sg)];
      max_raw += w * s_max[static_cast<size_t>(sg)];
      mean_off += w * base_eta_obs[static_cast<size_t>(i)];
      w_tot += w;
    }
    if (w_tot <= 0.0 || max_raw <= min_raw) continue;
    double prop = ((raw - min_raw) + 0.5) / ((max_raw - min_raw) + 1.0);
    double th = std::log(prop / (1.0 - prop)) - (mean_off / w_tot);
    th = std::max(-8.0, std::min(8.0, th));

    double lo = -25.0, hi = 25.0;
    auto eval_wle_g = [&](double t, double& g_out, double& dg_out, double& I_out) {
      double S_diff = 0.0, I = 0.0, dI = 0.0, d2I = 0.0;
      for (int idx = beg; idx < end; ++idx) {
        int i = spec.p_obs_idx[static_cast<size_t>(idx)];
        double w = spec.weight[static_cast<size_t>(i)];
        if (w <= 0.0) continue;
        int sg = spec.step_group[static_cast<size_t>(i)];
        const double* sc = step_cum_mat.data() + static_cast<size_t>(sg) * C;
        double al = (alpha_obs_ptr && i < static_cast<int>(alpha_obs_ptr->size()))
                        ? (*alpha_obs_ptr)[static_cast<size_t>(i)]
                        : 1.0;
        CatMoments cm = eval_gpcm_probs_and_moments(
            t, base_eta_obs[static_cast<size_t>(i)], al, sform, sc, C,
            spec.score[static_cast<size_t>(i)], probs.data(),
            s_min[static_cast<size_t>(sg)], s_max[static_cast<size_t>(sg)]);
        double al2 = al * al;
        S_diff += w * al * (static_cast<double>(spec.score[static_cast<size_t>(i)]) - cm.mean);
        I += w * al2 * cm.var;
        dI += w * al2 * al * cm.m3;
        d2I += w * al2 * al2 * (cm.m4 - 3.0 * cm.var * cm.var);
      }
      I_out = I;
      if (I <= 1e-15) {
        g_out = S_diff;
        dg_out = -1e-8;
        return;
      }
      g_out = S_diff + dI / (2.0 * I);
      dg_out = -I + (d2I * I - dI * dI) / (2.0 * I * I);
    };

    double g_lo = 0.0, dg_lo = 0.0, I_lo = 0.0;
    double g_hi = 0.0, dg_hi = 0.0, I_hi = 0.0;
    eval_wle_g(lo, g_lo, dg_lo, I_lo);
    eval_wle_g(hi, g_hi, dg_hi, I_hi);

    for (int iter = 0; iter < 45; ++iter) {
      double g = 0.0, dg = 0.0, I = 0.0;
      eval_wle_g(th, g, dg, I);
      if (g > 0.0) lo = std::max(lo, th);
      else if (g < 0.0) hi = std::min(hi, th);
      if (std::abs(g) < 1e-11) break;
      double step = (dg < -1e-12) ? (-g / dg) : (g > 0.0 ? 0.5 : -0.5);
      double cand = th + std::max(-2.0, std::min(2.0, step));
      if (cand <= lo || cand >= hi) cand = 0.5 * (lo + hi);
      if (std::abs(cand - th) < 1e-12) { th = cand; break; }
      th = cand;
    }
    double g_fin = 0.0, dg_fin = 0.0, I_fin = 0.0;
    eval_wle_g(th, g_fin, dg_fin, I_fin);
    res.theta[static_cast<size_t>(p)] = th;
    res.se[static_cast<size_t>(p)] = (I_fin > 1e-14) ? (1.0 / std::sqrt(I_fin)) : kNaN;
  }
  return res;
}

// ============================================================================
// 6. Single-Pass Diagnostics, Fair(M), Fair(Z), Real SE, & RajuSA-Welch Bias
// ============================================================================

inline double compute_zstd_facets(double mnsq, double df, bool whexact, double cap = 9.0) {
  if (!std::isfinite(mnsq) || !std::isfinite(df) || mnsq < 0.0 || df < 1.0) return kNaN;
  double q = std::sqrt(2.0 / df);
  double z = 3.0 * (std::cbrt(mnsq) - 1.0) / q + (whexact ? 0.0 : (q / 3.0));
  if (cap > 0.0 && std::isfinite(cap)) {
    z = std::max(-cap, std::min(cap, z));
  }
  return z;
}

struct ElementDiagTable {
  std::vector<double> obs_n;
  std::vector<double> obs_score;
  std::vector<double> exp_score;
  std::vector<double> cond_info;
  std::vector<double> infit_mnsq;
  std::vector<double> infit_zstd;
  std::vector<double> outfit_mnsq;
  std::vector<double> outfit_zstd;
  std::vector<double> fair_m;
  std::vector<double> fair_z;
};

struct FastFitOutput {
  bool converged = false;
  int iterations = 0;
  double loglik = kNaN;
  double objective = kNaN;
  double max_abs_grad = kNaN;
  int n_free_par = 0;

  ConnectivityResult connectivity;

  std::vector<double> person_est;
  std::vector<double> person_se;
  std::vector<double> person_real_se;
  std::vector<double> person_cond_se;
  std::vector<double> person_eap;
  std::vector<double> person_psd;
  std::vector<double> person_wle;
  std::vector<double> person_wle_se;
  ElementDiagTable person_diag;

  std::vector<std::vector<double>> facet_est;
  std::vector<std::vector<double>> facet_se;
  std::vector<std::vector<double>> facet_real_se;
  std::vector<std::vector<double>> facet_joint_se;
  std::vector<std::vector<double>> facet_cov_model;
  std::vector<std::vector<double>> facet_cov_sandwich;
  std::vector<ElementDiagTable> facet_diag;

  std::vector<double> step_est;
  std::vector<double> step_se;
  std::vector<double> step_cov_model;
  std::vector<double> step_cov_sandwich;
  std::vector<double> step_cum;

  // Optional GPCM / GMFRM discrimination estimates & delta-method covariances
  std::vector<std::vector<double>> log_slope_est;
  std::vector<std::vector<double>> slope_est;
  std::vector<std::vector<double>> slope_se;
  std::vector<std::vector<double>> slope_joint_se;
  std::vector<std::vector<double>> slope_cov_model;
  std::vector<std::vector<double>> slope_cov_sandwich;
  std::vector<double> alpha_obs;

  // Direct unconstrained structural coordinates xi (for explanatory formulas, step constraints, warm starts, & vcov)
  std::vector<double> xi_est;
  std::vector<double> xi_cov_model;
  std::vector<double> xi_cov_sandwich;

  std::vector<double> pop_beta;
  std::vector<double> pop_se_model;
  std::vector<double> pop_se_sandwich;
  double pop_sigma = 1.0;
  double person_rel_eap = kNaN;

  struct LocalBiasCell {
    int level_a = 0;
    int level_b = 0;
    int raw_obs_n = 0;
    double obs_count = 0.0;
    double obs_score = 0.0;
    double exp_score = 0.0;
    double bias_est = 0.0;
    double bias_se = kNaN;
    double t_stat = kNaN;
    double df = 0.0;
  };
  std::vector<LocalBiasCell> local_bias_cells;
};

inline void compute_fast_diagnostics_and_fair(
    const FastDesignSpec& spec,
    const std::vector<double>& theta,
    const std::vector<bool>& theta_exclude_from_diag,
    const std::vector<std::vector<double>>& beta,
    const std::vector<double>& step_cum_mat,
    bool whexact,
    ElementDiagTable& p_diag,
    std::vector<ElementDiagTable>& f_diag,
    const std::vector<double>* alpha_obs_ptr = nullptr) {
  const int P = spec.n_person;
  const int F = spec.n_facets;
  const int C = spec.n_cat;
  const int K = C - 1;
  const int sform = spec.slope_formula;
  std::vector<int> s_min = spec.step_cat_min.empty() ? std::vector<int>(static_cast<size_t>(spec.n_step_groups), 0) : spec.step_cat_min;
  std::vector<int> s_max = spec.step_cat_max.empty() ? std::vector<int>(static_cast<size_t>(spec.n_step_groups), K) : spec.step_cat_max;

  auto init_tbl = [](ElementDiagTable& t, int len) {
    t.obs_n.assign(static_cast<size_t>(len), 0.0);
    t.obs_score.assign(static_cast<size_t>(len), 0.0);
    t.exp_score.assign(static_cast<size_t>(len), 0.0);
    t.cond_info.assign(static_cast<size_t>(len), 0.0);
    t.infit_mnsq.assign(static_cast<size_t>(len), kNaN);
    t.infit_zstd.assign(static_cast<size_t>(len), kNaN);
    t.outfit_mnsq.assign(static_cast<size_t>(len), kNaN);
    t.outfit_zstd.assign(static_cast<size_t>(len), kNaN);
    t.fair_m.assign(static_cast<size_t>(len), kNaN);
    t.fair_z.assign(static_cast<size_t>(len), kNaN);
  };
  init_tbl(p_diag, P);
  f_diag.resize(static_cast<size_t>(F));
  for (int f = 0; f < F; ++f) init_tbl(f_diag[static_cast<size_t>(f)], spec.facet_n_levels[static_cast<size_t>(f)]);

  std::vector<double> p_score_var(static_cast<size_t>(P), 0.0);
  std::vector<double> p_sq_res(static_cast<size_t>(P), 0.0), p_std_sq(static_cast<size_t>(P), 0.0);
  std::vector<double> p_f4_in(static_cast<size_t>(P), 0.0), p_f4_out(static_cast<size_t>(P), 0.0);
  std::vector<std::vector<double>> f_score_var(static_cast<size_t>(F));
  std::vector<std::vector<double>> f_sq_res(static_cast<size_t>(F)), f_std_sq(static_cast<size_t>(F));
  std::vector<std::vector<double>> f_f4_in(static_cast<size_t>(F)), f_f4_out(static_cast<size_t>(F));
  for (int f = 0; f < F; ++f) {
    int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
    f_score_var[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), 0.0);
    f_sq_res[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), 0.0);
    f_std_sq[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), 0.0);
    f_f4_in[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), 0.0);
    f_f4_out[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), 0.0);
  }

  std::vector<double> probs(static_cast<size_t>(C));
  for (int i = 0; i < spec.n_obs; ++i) {
    double w = spec.weight[static_cast<size_t>(i)];
    if (w <= 0.0) continue;
    int p = spec.person[static_cast<size_t>(i)];
    if (theta_exclude_from_diag[static_cast<size_t>(p)] || !std::isfinite(theta[static_cast<size_t>(p)])) {
      continue;
    }
    double th_p = theta[static_cast<size_t>(p)];
    double be = 0.0;
    for (int f = 0; f < F; ++f) {
      int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
      be += spec.facet_signs[static_cast<size_t>(f)] * beta[static_cast<size_t>(f)][static_cast<size_t>(l)];
    }
    int sg = spec.step_group[static_cast<size_t>(i)];
    int y = spec.score[static_cast<size_t>(i)];
    double al = (alpha_obs_ptr && i < static_cast<int>(alpha_obs_ptr->size()))
                    ? (*alpha_obs_ptr)[static_cast<size_t>(i)]
                    : 1.0;
    CatMoments cm = eval_gpcm_probs_and_moments(
        th_p, be, al, sform, step_cum_mat.data() + static_cast<size_t>(sg) * C, C, y, probs.data(),
        s_min[static_cast<size_t>(sg)], s_max[static_cast<size_t>(sg)]);
    if (cm.var <= 1e-12 || !std::isfinite(cm.var)) continue;

    double res = static_cast<double>(y) - cm.mean;
    double res2 = res * res;
    double std2 = res2 / cm.var;
    double f4_in = w * (cm.m4 - cm.var * cm.var);
    double f4_out = w * (cm.m4 / (cm.var * cm.var) - 1.0);

    p_diag.obs_n[static_cast<size_t>(p)] += w;
    p_diag.obs_score[static_cast<size_t>(p)] += w * (spec.rating_min + y);
    p_diag.exp_score[static_cast<size_t>(p)] += w * (spec.rating_min + cm.mean);
    p_diag.cond_info[static_cast<size_t>(p)] += w * (al * al * cm.var);
    p_score_var[static_cast<size_t>(p)] += w * cm.var;
    p_sq_res[static_cast<size_t>(p)] += w * res2;
    p_std_sq[static_cast<size_t>(p)] += w * std2;
    p_f4_in[static_cast<size_t>(p)] += f4_in;
    p_f4_out[static_cast<size_t>(p)] += f4_out;

    for (int f = 0; f < F; ++f) {
      int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
      f_diag[static_cast<size_t>(f)].obs_n[static_cast<size_t>(l)] += w;
      f_diag[static_cast<size_t>(f)].obs_score[static_cast<size_t>(l)] += w * (spec.rating_min + y);
      f_diag[static_cast<size_t>(f)].exp_score[static_cast<size_t>(l)] += w * (spec.rating_min + cm.mean);
      f_diag[static_cast<size_t>(f)].cond_info[static_cast<size_t>(l)] += w * (al * al * cm.var);
      f_score_var[static_cast<size_t>(f)][static_cast<size_t>(l)] += w * cm.var;
      f_sq_res[static_cast<size_t>(f)][static_cast<size_t>(l)] += w * res2;
      f_std_sq[static_cast<size_t>(f)][static_cast<size_t>(l)] += w * std2;
      f_f4_in[static_cast<size_t>(f)][static_cast<size_t>(l)] += f4_in;
      f_f4_out[static_cast<size_t>(f)][static_cast<size_t>(l)] += f4_out;
    }
  }

  auto finalize_fit = [&](ElementDiagTable& t, const std::vector<double>& sc_v,
                          const std::vector<double>& sq_r,
                          const std::vector<double>& std_sq, const std::vector<double>& f4_i,
                          const std::vector<double>& f4_o) {
    for (size_t j = 0; j < t.obs_n.size(); ++j) {
      if (t.obs_n[j] <= 0.0 || sc_v[j] <= 1e-12) continue;
      double infit = sq_r[j] / sc_v[j];
      double outfit = std_sq[j] / t.obs_n[j];
      t.infit_mnsq[j] = infit;
      t.outfit_mnsq[j] = outfit;
      double df_in = (f4_i[j] > 1e-12) ? (2.0 * sc_v[j] * sc_v[j] / f4_i[j]) : kNaN;
      double df_out = (f4_o[j] > 1e-12) ? (2.0 * t.obs_n[j] * t.obs_n[j] / f4_o[j]) : kNaN;
      t.infit_zstd[j] = compute_zstd_facets(infit, df_in, whexact, 9.0);
      t.outfit_zstd[j] = compute_zstd_facets(outfit, df_out, whexact, 9.0);
    }
  };
  finalize_fit(p_diag, p_score_var, p_sq_res, p_std_sq, p_f4_in, p_f4_out);
  for (int f = 0; f < F; ++f) {
    finalize_fit(f_diag[static_cast<size_t>(f)], f_score_var[static_cast<size_t>(f)],
                 f_sq_res[static_cast<size_t>(f)],
                 f_std_sq[static_cast<size_t>(f)], f_f4_in[static_cast<size_t>(f)],
                 f_f4_out[static_cast<size_t>(f)]);
  }

  double theta_sum = 0.0;
  int theta_cnt = 0;
  for (int p = 0; p < P; ++p) {
    if (!theta_exclude_from_diag[static_cast<size_t>(p)] && std::isfinite(theta[static_cast<size_t>(p)])) {
      theta_sum += theta[static_cast<size_t>(p)];
      theta_cnt++;
    }
  }
  double theta_mean = (theta_cnt > 0) ? (theta_sum / theta_cnt) : kNaN;
  std::vector<double> facet_means(static_cast<size_t>(F), 0.0);
  double all_facet_signed_mean = 0.0;
  for (int f = 0; f < F; ++f) {
    double s = 0.0;
    int cnt = 0;
    for (double v : beta[static_cast<size_t>(f)]) {
      if (std::isfinite(v)) { s += v; cnt++; }
    }
    facet_means[static_cast<size_t>(f)] = (cnt > 0) ? (s / cnt) : 0.0;
    all_facet_signed_mean += spec.facet_signs[static_cast<size_t>(f)] * facet_means[static_cast<size_t>(f)];
  }

  std::vector<double> p_eta_m(static_cast<size_t>(P), kNaN), p_eta_z(static_cast<size_t>(P), kNaN);
  for (int p = 0; p < P; ++p) {
    if (std::isfinite(theta[static_cast<size_t>(p)])) {
      p_eta_m[static_cast<size_t>(p)] = theta[static_cast<size_t>(p)] + all_facet_signed_mean;
      p_eta_z[static_cast<size_t>(p)] = theta[static_cast<size_t>(p)];
    }
  }
  std::vector<std::vector<double>> f_eta_m(static_cast<size_t>(F)), f_eta_z(static_cast<size_t>(F));
  for (int f = 0; f < F; ++f) {
    int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
    f_eta_m[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), kNaN);
    f_eta_z[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), kNaN);
    double s_f = spec.facet_signs[static_cast<size_t>(f)];
    double other_sum = all_facet_signed_mean - s_f * facet_means[static_cast<size_t>(f)];
    for (int l = 0; l < Lf; ++l) {
      double b = beta[static_cast<size_t>(f)][static_cast<size_t>(l)];
      if (std::isfinite(b)) {
        f_eta_m[static_cast<size_t>(f)][static_cast<size_t>(l)] = theta_mean + other_sum + s_f * b;
        f_eta_z[static_cast<size_t>(f)][static_cast<size_t>(l)] = s_f * b;
      }
    }
  }

  if (spec.model_code == 0 && spec.n_step_groups == 1 && spec.n_slope_facets == 0) {
    const double* sc = step_cum_mat.data();
    int a0 = s_min[0], b0 = s_max[0];
    for (int p = 0; p < P; ++p) {
      if (std::isfinite(p_eta_m[static_cast<size_t>(p)])) {
        p_diag.fair_m[static_cast<size_t>(p)] = expected_score_from_eta(p_eta_m[static_cast<size_t>(p)], sc, C, spec.rating_min, a0, b0);
        p_diag.fair_z[static_cast<size_t>(p)] = expected_score_from_eta(p_eta_z[static_cast<size_t>(p)], sc, C, spec.rating_min, a0, b0);
      }
    }
    for (int f = 0; f < F; ++f) {
      for (int l = 0; l < spec.facet_n_levels[static_cast<size_t>(f)]; ++l) {
        if (std::isfinite(f_eta_m[static_cast<size_t>(f)][static_cast<size_t>(l)])) {
          f_diag[static_cast<size_t>(f)].fair_m[static_cast<size_t>(l)] =
              expected_score_from_eta(f_eta_m[static_cast<size_t>(f)][static_cast<size_t>(l)], sc, C, spec.rating_min, a0, b0);
          f_diag[static_cast<size_t>(f)].fair_z[static_cast<size_t>(l)] =
              expected_score_from_eta(f_eta_z[static_cast<size_t>(f)][static_cast<size_t>(l)], sc, C, spec.rating_min, a0, b0);
        }
      }
    }
  } else {
    // Compute geometric mean of observed slopes as reference baseline for Person Fair(M)/Fair(Z)
    double gm_alpha = 1.0;
    if (alpha_obs_ptr && !alpha_obs_ptr->empty()) {
      double ls_sum = 0.0, w_sum = 0.0;
      for (int i = 0; i < spec.n_obs; ++i) {
        double w = spec.weight[static_cast<size_t>(i)];
        double al = (*alpha_obs_ptr)[static_cast<size_t>(i)];
        if (w > 0.0 && std::isfinite(al) && al > 0.0) { ls_sum += w * std::log(al); w_sum += w; }
      }
      if (w_sum > 0.0) gm_alpha = std::exp(ls_sum / w_sum);
    }
    std::vector<double> pw_tot(static_cast<size_t>(P), 0.0), p_fm(static_cast<size_t>(P), 0.0), p_fz(static_cast<size_t>(P), 0.0);
    std::vector<std::vector<double>> fw_tot(static_cast<size_t>(F)), f_fm(static_cast<size_t>(F)), f_fz(static_cast<size_t>(F));
    for (int f = 0; f < F; ++f) {
      int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
      fw_tot[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), 0.0);
      f_fm[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), 0.0);
      f_fz[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), 0.0);
    }
    for (int i = 0; i < spec.n_obs; ++i) {
      double w = spec.weight[static_cast<size_t>(i)];
      if (w <= 0.0) continue;
      int p = spec.person[static_cast<size_t>(i)];
      int sg = spec.step_group[static_cast<size_t>(i)];
      int ag = s_min[static_cast<size_t>(sg)], bg = s_max[static_cast<size_t>(sg)];
      const double* sc = step_cum_mat.data() + static_cast<size_t>(sg) * C;
      double al_i = (alpha_obs_ptr && i < static_cast<int>(alpha_obs_ptr->size()))
                        ? (*alpha_obs_ptr)[static_cast<size_t>(i)]
                        : 1.0;
      if (std::isfinite(theta[static_cast<size_t>(p)])) {
        pw_tot[static_cast<size_t>(p)] += w;
        p_fm[static_cast<size_t>(p)] += w * expected_score_from_gpcm(
            theta[static_cast<size_t>(p)], all_facet_signed_mean, gm_alpha, sform, sc, C, spec.rating_min, ag, bg);
        p_fz[static_cast<size_t>(p)] += w * expected_score_from_gpcm(
            theta[static_cast<size_t>(p)], 0.0, gm_alpha, sform, sc, C, spec.rating_min, ag, bg);
      }
      for (int f = 0; f < F; ++f) {
        int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
        double b = beta[static_cast<size_t>(f)][static_cast<size_t>(l)];
        if (std::isfinite(b)) {
          double s_f = spec.facet_signs[static_cast<size_t>(f)];
          double other_sum = all_facet_signed_mean - s_f * facet_means[static_cast<size_t>(f)];
          fw_tot[static_cast<size_t>(f)][static_cast<size_t>(l)] += w;
          f_fm[static_cast<size_t>(f)][static_cast<size_t>(l)] += w * expected_score_from_gpcm(
              std::isfinite(theta_mean) ? theta_mean : 0.0, other_sum + s_f * b, al_i, sform, sc, C, spec.rating_min, ag, bg);
          f_fz[static_cast<size_t>(f)][static_cast<size_t>(l)] += w * expected_score_from_gpcm(
              0.0, s_f * b, al_i, sform, sc, C, spec.rating_min, ag, bg);
        }
      }
    }
    for (int p = 0; p < P; ++p) {
      if (pw_tot[static_cast<size_t>(p)] > 0.0) {
        p_diag.fair_m[static_cast<size_t>(p)] = p_fm[static_cast<size_t>(p)] / pw_tot[static_cast<size_t>(p)];
        p_diag.fair_z[static_cast<size_t>(p)] = p_fz[static_cast<size_t>(p)] / pw_tot[static_cast<size_t>(p)];
      }
    }
    for (int f = 0; f < F; ++f) {
      for (int l = 0; l < spec.facet_n_levels[static_cast<size_t>(f)]; ++l) {
        if (fw_tot[static_cast<size_t>(f)][static_cast<size_t>(l)] > 0.0) {
          f_diag[static_cast<size_t>(f)].fair_m[static_cast<size_t>(l)] =
              f_fm[static_cast<size_t>(f)][static_cast<size_t>(l)] / fw_tot[static_cast<size_t>(f)][static_cast<size_t>(l)];
          f_diag[static_cast<size_t>(f)].fair_z[static_cast<size_t>(l)] =
              f_fz[static_cast<size_t>(f)][static_cast<size_t>(l)] / fw_tot[static_cast<size_t>(f)][static_cast<size_t>(l)];
        }
      }
    }
  }
}

// ============================================================================
// 7. Two-Way Fixed-Measure Local Bias & Pairwise RajuSA-Welch (FACETS 4.5.3)
// ============================================================================

inline std::vector<FastFitOutput::LocalBiasCell> compute_local_bias_2way(
    const FastDesignSpec& spec,
    const std::vector<double>& theta,
    const std::vector<bool>& theta_exclude,
    const std::vector<std::vector<double>>& beta,
    const std::vector<double>& step_cum_mat,
    int facet_a, int facet_b,
    const std::vector<double>* alpha_obs_ptr = nullptr) {
  std::vector<FastFitOutput::LocalBiasCell> out;
  const int F = spec.n_facets;
  if (facet_a < 0 || facet_a >= F || facet_b < 0 || facet_b >= F || facet_a == facet_b) return out;
  const int La = spec.facet_n_levels[static_cast<size_t>(facet_a)];
  const int Lb = spec.facet_n_levels[static_cast<size_t>(facet_b)];
  const int C = spec.n_cat;
  const int K = C - 1;
  const int sform = spec.slope_formula;
  std::vector<int> s_min = spec.step_cat_min.empty() ? std::vector<int>(static_cast<size_t>(spec.n_step_groups), 0) : spec.step_cat_min;
  std::vector<int> s_max = spec.step_cat_max.empty() ? std::vector<int>(static_cast<size_t>(spec.n_step_groups), K) : spec.step_cat_max;
  const int n_cells = La * Lb;

  std::vector<std::vector<int>> cell_rows(static_cast<size_t>(n_cells));
  std::vector<double> base_off_i(static_cast<size_t>(spec.n_obs), kNaN);
  for (int i = 0; i < spec.n_obs; ++i) {
    double w = spec.weight[static_cast<size_t>(i)];
    if (w <= 0.0) continue;
    int p = spec.person[static_cast<size_t>(i)];
    if (theta_exclude[static_cast<size_t>(p)] || !std::isfinite(theta[static_cast<size_t>(p)])) continue;
    double off = 0.0;
    for (int f = 0; f < F; ++f) {
      int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
      off += spec.facet_signs[static_cast<size_t>(f)] * beta[static_cast<size_t>(f)][static_cast<size_t>(l)];
    }
    base_off_i[static_cast<size_t>(i)] = off;
    int la = spec.facet_flat[static_cast<size_t>(i) * F + facet_a];
    int lb = spec.facet_flat[static_cast<size_t>(i) * F + facet_b];
    cell_rows[static_cast<size_t>(la) * Lb + lb].push_back(i);
  }

  std::vector<double> probs(static_cast<size_t>(C));
  double sign_ab = spec.facet_signs[static_cast<size_t>(facet_a)];
  for (int la = 0; la < La; ++la) {
    for (int lb = 0; lb < Lb; ++lb) {
      const auto& rows = cell_rows[static_cast<size_t>(la) * Lb + lb];
      if (rows.empty()) continue;
      double obs_cnt = 0.0, obs_sum = 0.0, exp0_sum = 0.0, min_sum = 0.0, max_sum = 0.0;
      for (int i : rows) {
        double w = spec.weight[static_cast<size_t>(i)];
        int p = spec.person[static_cast<size_t>(i)];
        int sg = spec.step_group[static_cast<size_t>(i)];
        int ag = s_min[static_cast<size_t>(sg)], bg = s_max[static_cast<size_t>(sg)];
        double al = (alpha_obs_ptr && i < static_cast<int>(alpha_obs_ptr->size()))
                        ? (*alpha_obs_ptr)[static_cast<size_t>(i)] : 1.0;
        CatMoments cm = eval_gpcm_probs_and_moments(
            theta[static_cast<size_t>(p)], base_off_i[static_cast<size_t>(i)], al, sform,
            step_cum_mat.data() + static_cast<size_t>(sg) * C, C,
            spec.score[static_cast<size_t>(i)], probs.data(), ag, bg);
        obs_cnt += w;
        obs_sum += w * spec.score[static_cast<size_t>(i)];
        exp0_sum += w * cm.mean;
        min_sum += w * ag;
        max_sum += w * bg;
      }
      double delta = 0.0;
      bool cell_extreme = (obs_sum <= min_sum + 1e-12 || obs_sum >= max_sum - 1e-12);
      double target_sum = obs_sum;
      if (cell_extreme) {
        target_sum = (obs_sum <= min_sum + 1e-12) ? (min_sum + 0.3) : (max_sum - 0.3);
      }
      double info_fin = 0.0;
      for (int iter = 0; iter < 40; ++iter) {
        double e_sum = 0.0, v_sum = 0.0, deriv_sum = 0.0;
        for (int i : rows) {
          double w = spec.weight[static_cast<size_t>(i)];
          int p = spec.person[static_cast<size_t>(i)];
          int sg = spec.step_group[static_cast<size_t>(i)];
          int ag = s_min[static_cast<size_t>(sg)], bg = s_max[static_cast<size_t>(sg)];
          double al = (alpha_obs_ptr && i < static_cast<int>(alpha_obs_ptr->size()))
                          ? (*alpha_obs_ptr)[static_cast<size_t>(i)] : 1.0;
          CatMoments cm = eval_gpcm_probs_and_moments(
              theta[static_cast<size_t>(p)], base_off_i[static_cast<size_t>(i)] + sign_ab * delta, al, sform,
              step_cum_mat.data() + static_cast<size_t>(sg) * C, C,
              spec.score[static_cast<size_t>(i)], probs.data(), ag, bg);
          double eff_sc = (sform == 1) ? 1.0 : al;
          e_sum += w * cm.mean;
          v_sum += w * eff_sc * eff_sc * cm.var;
          deriv_sum += w * eff_sc * (target_sum / obs_cnt - cm.mean);
        }
        (void)deriv_sum;
        info_fin = v_sum;
        double g = sign_ab * ((sform == 1 || !alpha_obs_ptr) ? (target_sum - e_sum) : ([&]() {
          double acc = 0.0;
          for (int i : rows) {
            double w = spec.weight[static_cast<size_t>(i)];
            int p = spec.person[static_cast<size_t>(i)];
            int sg = spec.step_group[static_cast<size_t>(i)];
            int ag = s_min[static_cast<size_t>(sg)], bg = s_max[static_cast<size_t>(sg)];
            double al = (*alpha_obs_ptr)[static_cast<size_t>(i)];
            CatMoments cm = eval_gpcm_probs_and_moments(
                theta[static_cast<size_t>(p)], base_off_i[static_cast<size_t>(i)] + sign_ab * delta, al, sform,
                step_cum_mat.data() + static_cast<size_t>(sg) * C, C,
                spec.score[static_cast<size_t>(i)], probs.data(), ag, bg);
            double yi_adj = cell_extreme ? (target_sum / obs_cnt) : static_cast<double>(spec.score[static_cast<size_t>(i)]);
            acc += w * al * (yi_adj - cm.mean);
          }
          return acc;
        })());
        if (std::abs(g) < 1e-8 || v_sum <= 1e-12) break;
        double step = std::max(-1.5, std::min(1.5, g / v_sum));
        delta += step;
        if (std::abs(step) < 1e-9) break;
      }
      FastFitOutput::LocalBiasCell cell;
      cell.level_a = la;
      cell.level_b = lb;
      cell.raw_obs_n = static_cast<int>(rows.size());
      cell.obs_count = obs_cnt;
      cell.obs_score = obs_cnt * spec.rating_min + obs_sum;
      cell.exp_score = obs_cnt * spec.rating_min + exp0_sum;
      cell.bias_est = delta;
      cell.bias_se = (info_fin > 1e-12) ? (1.0 / std::sqrt(info_fin)) : kNaN;
      cell.df = (obs_cnt > 1.0) ? (obs_cnt - 1.0) : kNaN;
      cell.t_stat = std::isfinite(cell.bias_se) ? (cell.bias_est / cell.bias_se) : kNaN;
      out.push_back(cell);
    }
  }
  return out;
}

// FACETS 4.5.2/4.5.3 RajuSA-Welch (von Davier 2025 modified df) & Classical Welch-Satterthwaite
inline double compute_raju_welch_df(double se1, double n1, double se2, double n2,
                                    bool facets_von_davier, double& df_unrounded) {
  df_unrounded = kNaN;
  if (!std::isfinite(se1) || !std::isfinite(se2) || se1 <= 0.0 || se2 <= 0.0) return kNaN;
  if (!std::isfinite(n1) || !std::isfinite(n2) || n1 <= 1.0 || n2 <= 1.0) return kNaN;
  double scale = std::max(se1, se2);
  double v1 = (se1 / scale) * (se1 / scale);
  double v2 = (se2 / scale) * (se2 / scale);
  double num = (v1 + v2) * (v1 + v2);
  if (!facets_von_davier) {
    double den = (v1 * v1) / (n1 - 1.0) + (v2 * v2) / (n2 - 1.0);
    if (den <= 0.0 || !std::isfinite(den)) return kNaN;
    df_unrounded = num / den;
    return (df_unrounded > 0.0) ? df_unrounded : kNaN;
  }
  // FACETS 4.5.0-4.5.3 (von Davier, 2025; labeled RajuSA-Welch in 4.5.2+):
  double den = (v1 * v1) / (n1 + 1.0) + (v2 * v2) / (n2 + 1.0);
  if (den <= 0.0 || !std::isfinite(den)) return kNaN;
  df_unrounded = (num / den) * ((n1 + n2 - 2.0) / (n1 + n2 + 2.0));
  double rounded_10 = std::round(df_unrounded * 1e10) * 1e-10;
  double df_int = std::floor(rounded_10);
  return (df_int > 0.0) ? df_int : kNaN;
}

} // namespace mfrm_fast

#endif // MFRMR_FAST_ENGINE_H

#ifndef MFRMR_FAST_SOLVERS_H
#define MFRMR_FAST_SOLVERS_H

#include "fast_engine.h"

namespace mfrm_fast {

// ============================================================================
// Compact Self-Contained L-BFGS Optimizer (m = 12, Strong Backtracking)
// ============================================================================

struct LbfgsResult {
  std::vector<double> x;
  double f = kNaN;
  std::vector<double> g;
  int iterations = 0;
  bool converged = false;
};

template <typename EvalFn>
inline LbfgsResult minimize_lbfgs(const std::vector<double>& x0,
                                  EvalFn&& eval_fg,
                                  int max_iter = 250,
                                  double grad_tol = 1e-8,
                                  double rel_ftol = 1e-12,
                                  int history_m = 12) {
  const int n = static_cast<int>(x0.size());
  LbfgsResult res;
  res.x = x0;
  res.g.assign(static_cast<size_t>(n), 0.0);
  if (n == 0) {
    res.f = eval_fg(res.x, res.g);
    res.converged = std::isfinite(res.f);
    return res;
  }

  res.f = eval_fg(res.x, res.g);
  if (!std::isfinite(res.f)) return res;

  auto max_abs = [](const std::vector<double>& v) {
    double m = 0.0;
    for (double val : v) m = std::max(m, std::abs(val));
    return m;
  };
  auto dot = [&](const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += a[static_cast<size_t>(i)] * b[static_cast<size_t>(i)];
    return s;
  };

  if (max_abs(res.g) <= grad_tol) {
    res.converged = true;
    return res;
  }

  std::vector<std::vector<double>> s_hist, y_hist;
  std::vector<double> rho_hist;
  s_hist.reserve(static_cast<size_t>(history_m));
  y_hist.reserve(static_cast<size_t>(history_m));
  rho_hist.reserve(static_cast<size_t>(history_m));

  std::vector<double> q(static_cast<size_t>(n)), dir(static_cast<size_t>(n));
  std::vector<double> x_cand(static_cast<size_t>(n)), g_cand(static_cast<size_t>(n));

  for (int iter = 1; iter <= max_iter; ++iter) {
    res.iterations = iter;
    q = res.g;
    int k_hist = static_cast<int>(s_hist.size());
    std::vector<double> alpha_hist(static_cast<size_t>(k_hist), 0.0);
    for (int i = k_hist - 1; i >= 0; --i) {
      double a = rho_hist[static_cast<size_t>(i)] * dot(s_hist[static_cast<size_t>(i)], q);
      alpha_hist[static_cast<size_t>(i)] = a;
      for (int j = 0; j < n; ++j) {
        q[static_cast<size_t>(j)] -= a * y_hist[static_cast<size_t>(i)][static_cast<size_t>(j)];
      }
    }
    double gamma = 1.0;
    if (k_hist > 0) {
      double ys = dot(y_hist.back(), s_hist.back());
      double yy = dot(y_hist.back(), y_hist.back());
      if (ys > 1e-16 && yy > 1e-16) gamma = ys / yy;
    } else {
      double gnorm = std::sqrt(dot(res.g, res.g));
      if (gnorm > 1.0) gamma = 1.0 / gnorm;
    }
    for (int j = 0; j < n; ++j) dir[static_cast<size_t>(j)] = gamma * q[static_cast<size_t>(j)];
    for (int i = 0; i < k_hist; ++i) {
      double b = rho_hist[static_cast<size_t>(i)] * dot(y_hist[static_cast<size_t>(i)], dir);
      double diff = alpha_hist[static_cast<size_t>(i)] - b;
      for (int j = 0; j < n; ++j) {
        dir[static_cast<size_t>(j)] += diff * s_hist[static_cast<size_t>(i)][static_cast<size_t>(j)];
      }
    }
    for (int j = 0; j < n; ++j) dir[static_cast<size_t>(j)] = -dir[static_cast<size_t>(j)];

    double dg = dot(dir, res.g);
    if (!std::isfinite(dg) || dg >= 0.0) {
      s_hist.clear(); y_hist.clear(); rho_hist.clear();
      double gnorm = std::max(1.0, std::sqrt(dot(res.g, res.g)));
      for (int j = 0; j < n; ++j) dir[static_cast<size_t>(j)] = -res.g[static_cast<size_t>(j)] / gnorm;
      dg = dot(dir, res.g);
    }

    double max_step = 0.0;
    for (double v : dir) max_step = std::max(max_step, std::abs(v));
    double step_scale = (max_step > 4.0) ? (4.0 / max_step) : 1.0;

    bool accepted = false;
    double f_cand = kNaN;
    for (int ls = 0; ls < 30; ++ls) {
      for (int j = 0; j < n; ++j) {
        x_cand[static_cast<size_t>(j)] = res.x[static_cast<size_t>(j)] + step_scale * dir[static_cast<size_t>(j)];
      }
      f_cand = eval_fg(x_cand, g_cand);
      if (std::isfinite(f_cand) && f_cand <= res.f + 1e-4 * step_scale * dg) {
        bool g_ok = true;
        for (double gv : g_cand) if (!std::isfinite(gv)) { g_ok = false; break; }
        if (g_ok) { accepted = true; break; }
      }
      step_scale *= 0.5;
    }

    if (!accepted) {
      if (max_abs(res.g) <= std::max(grad_tol * 100.0, 1e-5)) res.converged = true;
      break;
    }

    std::vector<double> s_vec(static_cast<size_t>(n)), y_vec(static_cast<size_t>(n));
    double max_dx = 0.0;
    for (int j = 0; j < n; ++j) {
      s_vec[static_cast<size_t>(j)] = x_cand[static_cast<size_t>(j)] - res.x[static_cast<size_t>(j)];
      y_vec[static_cast<size_t>(j)] = g_cand[static_cast<size_t>(j)] - res.g[static_cast<size_t>(j)];
      max_dx = std::max(max_dx, std::abs(s_vec[static_cast<size_t>(j)]));
    }
    double ys = dot(y_vec, s_vec);
    if (ys > 1e-14) {
      if (static_cast<int>(s_hist.size()) >= history_m) {
        s_hist.erase(s_hist.begin());
        y_hist.erase(y_hist.begin());
        rho_hist.erase(rho_hist.begin());
      }
      s_hist.push_back(std::move(s_vec));
      y_hist.push_back(std::move(y_vec));
      rho_hist.push_back(1.0 / ys);
    }

    double f_prev = res.f;
    res.x = x_cand;
    res.f = f_cand;
    res.g = g_cand;

    double gmax = max_abs(res.g);
    if (gmax <= grad_tol) {
      res.converged = true;
      break;
    }
    if (std::abs(f_prev - res.f) <= rel_ftol * (1.0 + std::abs(f_prev)) &&
        max_dx <= 1e-10 && gmax <= std::max(grad_tol * 10.0, 1e-6)) {
      res.converged = true;
      break;
    }
  }
  return res;
}

// ============================================================================
// Control Struct (with Optional Warm-Start Vector for Fast Split-Rater DFF)
// ============================================================================

struct FastControl {
  int method_code = 0; // 0 = JML, 1 = JML_jeffreys, 2 = MML, 3 = CML
  int maxit = 400;
  double reltol = 1e-10;
  int quad_points = 31;
  bool adaptive_quad = false;
  double xtreme = 0.0;
  bool whexact = false;
  bool compute_diagnostics = true;
  bool compute_wle = true;
  bool wle_primary = false;
  int bias_facet_a = -1;
  int bias_facet_b = -1;
  double facet_l2_pen = 0.0;
  double slope_l2_pen = 0.0;
  std::vector<double> warm_init;
};

// Optional L2 Gaussian prior penalty on facet locations and/or log-slopes
inline double add_structural_l2_penalties(
    const FastDesignSpec& spec,
    const StructuralMap& smap,
    const FastControl& ctrl,
    const std::vector<std::vector<double>>& beta,
    const std::vector<std::vector<double>>& log_slope_est,
    std::vector<double>& grad) {
  double pen = 0.0;
  if (ctrl.facet_l2_pen > 0.0) {
    const double lam = ctrl.facet_l2_pen;
    for (int f = 0; f < spec.n_facets; ++f) {
      if (spec.facet_dummy[static_cast<size_t>(f)]) continue;
      int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
      int qf = smap.facet_xi_len[static_cast<size_t>(f)];
      int st = smap.facet_xi_start[static_cast<size_t>(f)];
      const auto& Bf = smap.facet_B[static_cast<size_t>(f)];
      for (int l = 0; l < Lf; ++l) {
        if (smap.facet_fixed_mask[static_cast<size_t>(f)][static_cast<size_t>(l)]) continue;
        double b = beta[static_cast<size_t>(f)][static_cast<size_t>(l)];
        pen += 0.5 * lam * b * b;
        if (qf > 0 && !grad.empty()) {
          const double* row = Bf.data() + static_cast<size_t>(l) * qf;
          for (int c = 0; c < qf; ++c) {
            grad[static_cast<size_t>(st + c)] += lam * b * row[c];
          }
        }
      }
    }
  }
  if (ctrl.slope_l2_pen > 0.0 && smap.n_slope_facets > 0) {
    const double lam_s = ctrl.slope_l2_pen;
    for (int s = 0; s < smap.n_slope_facets; ++s) {
      int f = smap.slope_facets[static_cast<size_t>(s)];
      int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
      int qs = smap.slope_q[static_cast<size_t>(s)];
      int st = smap.slope_xi_offset[static_cast<size_t>(s)];
      const auto& Bs = smap.slope_B[static_cast<size_t>(s)];
      for (int l = 0; l < Lf; ++l) {
        if (!smap.slope_fixed.empty() && smap.slope_fixed[static_cast<size_t>(s)][static_cast<size_t>(l)]) continue;
        double g = log_slope_est[static_cast<size_t>(s)][static_cast<size_t>(l)];
        pen += 0.5 * lam_s * g * g;
        if (qs > 0 && !grad.empty()) {
          const double* row = Bs.data() + static_cast<size_t>(l) * qs;
          for (int c = 0; c < qs; ++c) {
            grad[static_cast<size_t>(st + c)] += lam_s * g * row[c];
          }
        }
      }
    }
  }
  return pen;
}

// Helper: populate z_xi [C x n_xi] for observation i with step group sg and probs[0..C-1]
inline void populate_step_z_xi(const StructuralMap& smap, int sg, int C,
                               const double* probs, std::vector<int>& act_cols,
                               double* z_xi) {
  const int qs = smap.step_xi_len;
  if (qs <= 0) return;
  const int n_xi = smap.n_xi;
  const int ag = smap.cat_min[static_cast<size_t>(sg)];
  const int bg = smap.cat_max[static_cast<size_t>(sg)];
  const auto& scols = smap.group_step_cols[static_cast<size_t>(sg)];
  for (int c_loc : scols) {
    int col = smap.step_xi_start + c_loc;
    act_cols.push_back(col);
    double e_a = 0.0;
    for (int k = ag; k <= bg; ++k) {
      e_a += probs[k] * smap.cum_B[(static_cast<size_t>(sg) * C + k) * qs + c_loc];
    }
    for (int k = 0; k < C; ++k) {
      if (k < ag || k > bg) {
        z_xi[static_cast<size_t>(k) * n_xi + col] = 0.0;
      } else {
        z_xi[static_cast<size_t>(k) * n_xi + col] =
            e_a - smap.cum_B[(static_cast<size_t>(sg) * C + k) * qs + c_loc];
      }
    }
  }
}

// ============================================================================
// Solver 1b: Generalized C++ JML for GPCM Slopes & L2-Penalized Facet/Slope Fits
// ============================================================================

inline void solve_jml_gpcm_or_penalized(const FastDesignSpec& spec,
                                        const StructuralMap& smap,
                                        const FastControl& ctrl,
                                        FastFitOutput& out) {
  const int P = spec.n_person;
  const int F = spec.n_facets;
  const int C = spec.n_cat;
  const int n_xi = smap.n_xi;
  const int qs = smap.step_xi_len;
  const bool has_slopes = (smap.n_slope_facets > 0);
  const int sform = spec.slope_formula;
  const auto lderivs = smap.build_level_derivs(spec);
  const auto sderivs = smap.build_slope_derivs(spec);

  std::vector<int> p_to_act(static_cast<size_t>(P), -1);
  std::vector<int> act_to_p;
  bool use_xtreme = (ctrl.xtreme > 0.0);
  for (int p = 0; p < P; ++p) {
    if (out.connectivity.person_obs_count[static_cast<size_t>(p)] > 0.0 &&
        (use_xtreme || out.connectivity.person_extreme_flag[static_cast<size_t>(p)] == 0)) {
      p_to_act[static_cast<size_t>(p)] = static_cast<int>(act_to_p.size());
      act_to_p.push_back(p);
    }
  }
  const int Pa = static_cast<int>(act_to_p.size());
  const int dim = Pa + n_xi;
  out.n_free_par = dim;

  std::vector<int> obs_f_ptr(static_cast<size_t>(spec.n_obs) + 1, 0);
  std::vector<int> obs_f_cols;
  std::vector<double> obs_f_vals;
  std::vector<int> obs_s_ptr(static_cast<size_t>(spec.n_obs) + 1, 0);
  std::vector<int> obs_s_cols;
  std::vector<double> obs_s_vals;
  {
    std::vector<int> tmp_c(std::max(512, n_xi + 16));
    std::vector<double> tmp_v(std::max(512, n_xi + 16));
    for (int i = 0; i < spec.n_obs; ++i) {
      int nnz = obs_eta_xi_terms(i, spec, lderivs, tmp_c.data(), tmp_v.data());
      for (int h = 0; h < nnz; ++h) {
        obs_f_cols.push_back(tmp_c[static_cast<size_t>(h)]);
        obs_f_vals.push_back(tmp_v[static_cast<size_t>(h)]);
      }
      obs_f_ptr[static_cast<size_t>(i) + 1] = static_cast<int>(obs_f_cols.size());
      if (has_slopes) {
        int snz = obs_log_slope_xi_terms(i, spec, sderivs, tmp_c.data(), tmp_v.data());
        for (int h = 0; h < snz; ++h) {
          obs_s_cols.push_back(tmp_c[static_cast<size_t>(h)]);
          obs_s_vals.push_back(tmp_v[static_cast<size_t>(h)]);
        }
      }
      obs_s_ptr[static_cast<size_t>(i) + 1] = static_cast<int>(obs_s_cols.size());
    }
  }

  // Center person ability when all base facets are centered
  bool center_persons = true;
  for (int f = 0; f < F; ++f) {
    if (spec.facet_noncenter[static_cast<size_t>(f)]) { center_persons = false; break; }
    for (int af : spec.anchor_fixed[static_cast<size_t>(f)]) {
      if (af == 1) { center_persons = false; break; }
    }
  }

  std::vector<double> x0(static_cast<size_t>(dim), 0.0);
  for (int ia = 0; ia < Pa; ++ia) {
    int p = act_to_p[static_cast<size_t>(ia)];
    double mn = out.connectivity.person_min_score[static_cast<size_t>(p)];
    double mx = out.connectivity.person_max_score[static_cast<size_t>(p)];
    double prop = ((out.connectivity.person_raw_score[static_cast<size_t>(p)] - mn) + 0.5) / ((mx - mn) + 1.0);
    x0[static_cast<size_t>(ia)] = std::log(prop / (1.0 - prop));
  }
  if (center_persons && Pa > 0) {
    double m = std::accumulate(x0.begin(), x0.begin() + Pa, 0.0) / Pa;
    for (int ia = 0; ia < Pa; ++ia) x0[static_cast<size_t>(ia)] -= m;
  }
  if (static_cast<int>(ctrl.warm_init.size()) == dim) x0 = ctrl.warm_init;

  std::vector<std::vector<double>> beta, log_slope_est, slope_est;
  std::vector<double> tau_mat, step_cum_mat, alpha_obs(static_cast<size_t>(spec.n_obs), 1.0);
  std::vector<double> base_eta_obs(static_cast<size_t>(spec.n_obs), 0.0);

  auto eval_jml_fg = [&](const std::vector<double>& x, std::vector<double>& g) -> double {
    g.assign(static_cast<size_t>(dim), 0.0);
    std::vector<double> th_act(x.begin(), x.begin() + Pa);
    double th_mean = 0.0;
    if (center_persons && Pa > 0) {
      th_mean = std::accumulate(th_act.begin(), th_act.end(), 0.0) / Pa;
      for (int ia = 0; ia < Pa; ++ia) th_act[static_cast<size_t>(ia)] -= th_mean;
    }
    for (int ia = 0; ia < Pa; ++ia) {
      if (!std::isfinite(th_act[static_cast<size_t>(ia)]) || std::abs(th_act[static_cast<size_t>(ia)]) > 35.0) return kInf;
    }
    const double* xi_ptr = x.data() + Pa;
    smap.expand(spec, xi_ptr, beta, tau_mat, step_cum_mat);
    if (!smap.expand_slopes(spec, xi_ptr, log_slope_est, slope_est, alpha_obs)) return kInf;
    for (int i = 0; i < spec.n_obs; ++i) {
      double e = 0.0;
      for (int f = 0; f < F; ++f) {
        int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
        e += spec.facet_signs[static_cast<size_t>(f)] * beta[static_cast<size_t>(f)][static_cast<size_t>(l)];
      }
      if (!std::isfinite(e) || std::abs(e) > 40.0) return kInf;
      base_eta_obs[static_cast<size_t>(i)] = e;
    }

    std::vector<double> probs(static_cast<size_t>(C));
    double nll = 0.0;
    for (int i = 0; i < spec.n_obs; ++i) {
      double w = spec.weight[static_cast<size_t>(i)];
      if (w <= 0.0) continue;
      int p = spec.person[static_cast<size_t>(i)];
      int ia = p_to_act[static_cast<size_t>(p)];
      if (ia < 0) continue;
      int sg = spec.step_group[static_cast<size_t>(i)];
      int ag = smap.cat_min[static_cast<size_t>(sg)], bg = smap.cat_max[static_cast<size_t>(sg)];
      int y = spec.score[static_cast<size_t>(i)];
      const double* sc = step_cum_mat.data() + static_cast<size_t>(sg) * C;
      double tp = th_act[static_cast<size_t>(ia)];
      double be = base_eta_obs[static_cast<size_t>(i)];
      double al = alpha_obs[static_cast<size_t>(i)];

      CatMoments cm = eval_gpcm_probs_and_moments(tp, be, al, sform, sc, C, y, probs.data(), ag, bg);
      if (!std::isfinite(cm.log_py)) return kInf;
      double y_eff = static_cast<double>(y);
      int ex_p = out.connectivity.person_extreme_flag[static_cast<size_t>(p)];
      if (use_xtreme && ex_p != 0) {
        double n_p = out.connectivity.person_obs_count[static_cast<size_t>(p)];
        y_eff += (ex_p < 0 ? 1.0 : -1.0) * (ctrl.xtreme / std::max(1.0, n_p));
        nll -= w * (cm.log_py + (y_eff - static_cast<double>(y)) * al * tp);
      } else {
        nll -= w * cm.log_py;
      }

      g[static_cast<size_t>(ia)] += w * al * (cm.mean - y_eff);
      double eff_sc = (!has_slopes || sform == 1) ? 1.0 : al;
      double g_eta = w * eff_sc * (cm.mean - y_eff);
      for (int ptr = obs_f_ptr[static_cast<size_t>(i)]; ptr < obs_f_ptr[static_cast<size_t>(i) + 1]; ++ptr) {
        g[static_cast<size_t>(Pa + obs_f_cols[static_cast<size_t>(ptr)])] += g_eta * obs_f_vals[static_cast<size_t>(ptr)];
      }
      const auto& scols = smap.group_step_cols[static_cast<size_t>(sg)];
      for (int c_loc : scols) {
        double ea = 0.0;
        for (int k = ag; k <= bg; ++k) {
          ea += probs[static_cast<size_t>(k)] * smap.cum_B[(static_cast<size_t>(sg) * C + k) * qs + c_loc];
        }
        double a_obs = smap.cum_B[(static_cast<size_t>(sg) * C + y) * qs + c_loc];
        g[static_cast<size_t>(Pa + smap.step_xi_start + c_loc)] += w * eff_sc * (a_obs - ea);
      }
      if (has_slopes) {
        double g_ls = 0.0;
        if (sform == 1) {
          g_ls = w * al * tp * (cm.mean - y_eff);
        } else {
          double e_sc = 0.0;
          for (int k = ag; k <= bg; ++k) {
            if (probs[static_cast<size_t>(k)] > 0.0 && std::isfinite(sc[k])) e_sc += probs[static_cast<size_t>(k)] * sc[k];
          }
          g_ls = w * al * ((cm.mean - y_eff) * (tp + be) + (sc[y] - e_sc));
        }
        for (int ptr = obs_s_ptr[static_cast<size_t>(i)]; ptr < obs_s_ptr[static_cast<size_t>(i) + 1]; ++ptr) {
          g[static_cast<size_t>(Pa + obs_s_cols[static_cast<size_t>(ptr)])] += g_ls * obs_s_vals[static_cast<size_t>(ptr)];
        }
      }
    }
    if (center_persons && Pa > 0) {
      double gm = std::accumulate(g.begin(), g.begin() + Pa, 0.0) / Pa;
      for (int ia = 0; ia < Pa; ++ia) g[static_cast<size_t>(ia)] -= gm;
      // Soft centering tether 0.5 * (mean(x_theta))^2 so L-BFGS removes gauge drift
      nll += 0.5 * Pa * th_mean * th_mean;
      for (int ia = 0; ia < Pa; ++ia) g[static_cast<size_t>(ia)] += th_mean;
    }
    if (ctrl.facet_l2_pen > 0.0 || ctrl.slope_l2_pen > 0.0) {
      std::vector<double> g_xi(static_cast<size_t>(n_xi), 0.0);
      nll += add_structural_l2_penalties(spec, smap, ctrl, beta, log_slope_est, g_xi);
      for (int j = 0; j < n_xi; ++j) g[static_cast<size_t>(Pa + j)] += g_xi[static_cast<size_t>(j)];
    }
    return nll;
  };

  LbfgsResult opt = minimize_lbfgs(x0, eval_jml_fg, std::max(ctrl.maxit, 350), 1e-8, 1e-13, 16);
  std::vector<double> x_hat = opt.x;
  if (center_persons && Pa > 0) {
    double m = std::accumulate(x_hat.begin(), x_hat.begin() + Pa, 0.0) / Pa;
    for (int ia = 0; ia < Pa; ++ia) x_hat[static_cast<size_t>(ia)] -= m;
  }
  std::vector<double> g_hat;
  double obj_hat = eval_jml_fg(x_hat, g_hat);

  const double* xi_hat = x_hat.data() + Pa;
  smap.expand(spec, xi_hat, out.facet_est, out.step_est, out.step_cum);
  smap.expand_slopes(spec, xi_hat, out.log_slope_est, out.slope_est, out.alpha_obs);
  out.xi_est.assign(xi_hat, xi_hat + n_xi);

  // Pure unpenalized log-likelihood at solution
  double pure_ll = 0.0;
  std::vector<double> probs(static_cast<size_t>(C));
  for (int i = 0; i < spec.n_obs; ++i) {
    double e = 0.0;
    for (int f = 0; f < F; ++f) {
      int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
      e += spec.facet_signs[static_cast<size_t>(f)] * out.facet_est[static_cast<size_t>(f)][static_cast<size_t>(l)];
    }
    base_eta_obs[static_cast<size_t>(i)] = e;
    double w = spec.weight[static_cast<size_t>(i)];
    int p = spec.person[static_cast<size_t>(i)];
    int ia = p_to_act[static_cast<size_t>(p)];
    if (w > 0.0 && ia >= 0) {
      int sg = spec.step_group[static_cast<size_t>(i)];
      CatMoments cm = eval_gpcm_probs_and_moments(
          x_hat[static_cast<size_t>(ia)], e, out.alpha_obs[static_cast<size_t>(i)], sform,
          out.step_cum.data() + static_cast<size_t>(sg) * C, C, spec.score[static_cast<size_t>(i)], probs.data(),
          smap.cat_min[static_cast<size_t>(sg)], smap.cat_max[static_cast<size_t>(sg)]);
      pure_ll += w * cm.log_py;
    }
  }

  out.iterations = opt.iterations;
  out.objective = obj_hat;
  out.loglik = pure_ll;
  double gmax = 0.0;
  for (double gv : g_hat) gmax = std::max(gmax, std::abs(gv));
  out.max_abs_grad = gmax;
  out.converged = std::isfinite(obj_hat) && (gmax <= 1e-4);

  // Block-Schur observed information SEs: H_schur = H_{xi,xi} - H_{xi,theta} D_theta^{-1} H_{theta,xi}
  std::vector<double> D_theta(static_cast<size_t>(Pa), 0.0);
  for (int i = 0; i < spec.n_obs; ++i) {
    double w = spec.weight[static_cast<size_t>(i)];
    int p = spec.person[static_cast<size_t>(i)];
    int ia = p_to_act[static_cast<size_t>(p)];
    if (w <= 0.0 || ia < 0) continue;
    int sg = spec.step_group[static_cast<size_t>(i)];
    double al = out.alpha_obs[static_cast<size_t>(i)];
    CatMoments cm = eval_gpcm_probs_and_moments(
        x_hat[static_cast<size_t>(ia)], base_eta_obs[static_cast<size_t>(i)], al, sform,
        out.step_cum.data() + static_cast<size_t>(sg) * C, C, spec.score[static_cast<size_t>(i)], probs.data(),
        smap.cat_min[static_cast<size_t>(sg)], smap.cat_max[static_cast<size_t>(sg)]);
    D_theta[static_cast<size_t>(ia)] += w * al * al * cm.var;
  }

  std::vector<double> H_xx(static_cast<size_t>(n_xi) * n_xi, 0.0);
  std::vector<double> H_tx(static_cast<size_t>(Pa) * n_xi, 0.0);
  if (n_xi > 0) {
    const double eps = 1e-5;
    std::vector<double> xp = x_hat, xm = x_hat, gp, gm;
    for (int j = 0; j < n_xi; ++j) {
      double h = eps * std::max(1.0, std::abs(x_hat[static_cast<size_t>(Pa + j)]));
      xp[static_cast<size_t>(Pa + j)] = x_hat[static_cast<size_t>(Pa + j)] + h;
      xm[static_cast<size_t>(Pa + j)] = x_hat[static_cast<size_t>(Pa + j)] - h;
      eval_jml_fg(xp, gp);
      eval_jml_fg(xm, gm);
      for (int ia = 0; ia < Pa; ++ia) {
        H_tx[static_cast<size_t>(ia) * n_xi + j] = (gp[static_cast<size_t>(ia)] - gm[static_cast<size_t>(ia)]) / (2.0 * h);
      }
      for (int k = 0; k < n_xi; ++k) {
        H_xx[static_cast<size_t>(j) * n_xi + k] = (gp[static_cast<size_t>(Pa + k)] - gm[static_cast<size_t>(Pa + k)]) / (2.0 * h);
      }
      xp[static_cast<size_t>(Pa + j)] = x_hat[static_cast<size_t>(Pa + j)];
      xm[static_cast<size_t>(Pa + j)] = x_hat[static_cast<size_t>(Pa + j)];
    }
    std::vector<double> H_schur = H_xx;
    for (int ia = 0; ia < Pa; ++ia) {
      double d = D_theta[static_cast<size_t>(ia)];
      if (d <= 1e-12) continue;
      double inv_d = 1.0 / d;
      const double* row = H_tx.data() + static_cast<size_t>(ia) * n_xi;
      for (int a = 0; a < n_xi; ++a) {
        for (int b = 0; b < n_xi; ++b) {
          H_schur[static_cast<size_t>(a) * n_xi + b] -= inv_d * row[a] * row[b];
        }
      }
    }
    std::vector<double> Vxi_joint, Vxi_marg;
    double ld1 = 0.0, ld2 = 0.0;
    chol_invert_and_logdet(H_xx, n_xi, Vxi_marg, ld1);
    if (!chol_invert_and_logdet(H_schur, n_xi, Vxi_joint, ld2)) Vxi_joint = Vxi_marg;
    out.xi_cov_model = Vxi_marg;
    out.xi_cov_sandwich = Vxi_joint;
    std::vector<double> dummy_step_se;
    smap.propagate_covariance(spec, Vxi_marg, out.facet_se, out.facet_cov_model, out.step_se, &out.step_cov_model);
    smap.propagate_covariance(spec, Vxi_joint, out.facet_joint_se, out.facet_cov_sandwich, dummy_step_se, &out.step_cov_sandwich);
    if (has_slopes) {
      smap.propagate_slope_covariance(spec, Vxi_marg, out.slope_est, out.slope_se, out.slope_cov_model);
      smap.propagate_slope_covariance(spec, Vxi_joint, out.slope_est, out.slope_joint_se, out.slope_cov_sandwich);
    }
  } else {
    smap.propagate_covariance(spec, {}, out.facet_se, out.facet_cov_model, out.step_se, &out.step_cov_model);
    smap.propagate_covariance(spec, {}, out.facet_joint_se, out.facet_cov_sandwich, out.step_se, &out.step_cov_sandwich);
  }

  out.person_est.assign(static_cast<size_t>(P), kNaN);
  out.person_se.assign(static_cast<size_t>(P), kNaN);
  out.person_cond_se.assign(static_cast<size_t>(P), kNaN);
  for (int ia = 0; ia < Pa; ++ia) {
    int p = act_to_p[static_cast<size_t>(ia)];
    out.person_est[static_cast<size_t>(p)] = x_hat[static_cast<size_t>(ia)];
    double se_p = (D_theta[static_cast<size_t>(ia)] > 1e-12) ? (1.0 / std::sqrt(D_theta[static_cast<size_t>(ia)])) : kNaN;
    out.person_se[static_cast<size_t>(p)] = se_p;
    out.person_cond_se[static_cast<size_t>(p)] = se_p;
  }
  const std::vector<double>* alpha_ptr = has_slopes ? &out.alpha_obs : nullptr;
  if (ctrl.compute_wle) {
    WleResult wle = score_wle_persons(spec, base_eta_obs, out.step_cum, alpha_ptr);
    out.person_wle = std::move(wle.theta);
    out.person_wle_se = std::move(wle.se);
  }
  std::vector<bool> diag_exclude(static_cast<size_t>(P), false);
  for (int p = 0; p < P; ++p) {
    if (p_to_act[static_cast<size_t>(p)] < 0) diag_exclude[static_cast<size_t>(p)] = true;
  }
  compute_fast_diagnostics_and_fair(
      spec, out.person_est, diag_exclude, out.facet_est, out.step_cum,
      ctrl.whexact, out.person_diag, out.facet_diag, alpha_ptr);

  out.person_real_se.assign(static_cast<size_t>(P), kNaN);
  for (int p = 0; p < P; ++p) {
    double se = out.person_se[static_cast<size_t>(p)];
    if (std::isfinite(se)) {
      double infit = out.person_diag.infit_mnsq[static_cast<size_t>(p)];
      out.person_real_se[static_cast<size_t>(p)] = std::isfinite(infit) ? (se * std::sqrt(std::max(1.0, infit))) : se;
    }
  }
  out.facet_real_se.resize(static_cast<size_t>(F));
  for (int f = 0; f < F; ++f) {
    int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
    out.facet_real_se[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), kNaN);
    for (int l = 0; l < Lf; ++l) {
      double se = out.facet_se[static_cast<size_t>(f)][static_cast<size_t>(l)];
      if (std::isfinite(se)) {
        double infit = out.facet_diag[static_cast<size_t>(f)].infit_mnsq[static_cast<size_t>(l)];
        out.facet_real_se[static_cast<size_t>(f)][static_cast<size_t>(l)] =
            std::isfinite(infit) ? (se * std::sqrt(std::max(1.0, infit))) : se;
      }
    }
  }
  if (ctrl.bias_facet_a >= 0 && ctrl.bias_facet_b >= 0) {
    out.local_bias_cells = compute_local_bias_2way(
        spec, out.person_est, diag_exclude, out.facet_est, out.step_cum,
        ctrl.bias_facet_a, ctrl.bias_facet_b, alpha_ptr);
  }
}

// ============================================================================
// Solver 1: Ultra-Fast Block-Schur Joint Newton JML
// ============================================================================

inline void solve_jml(const FastDesignSpec& spec, const StructuralMap& smap,
                      const FastControl& ctrl, FastFitOutput& out) {
  if (smap.n_slope_facets > 0 || ctrl.facet_l2_pen > 0.0 || ctrl.slope_l2_pen > 0.0) {
    solve_jml_gpcm_or_penalized(spec, smap, ctrl, out);
    return;
  }
  const int P = spec.n_person;
  const int F = spec.n_facets;
  const int C = spec.n_cat;
  const int K = C - 1;
  const int n_xi = smap.n_xi;
  const auto lderivs = smap.build_level_derivs(spec);

  std::vector<int> p_to_act(static_cast<size_t>(P), -1);
  std::vector<int> act_to_p;
  for (int p = 0; p < P; ++p) {
    if (out.connectivity.person_obs_count[static_cast<size_t>(p)] > 0.0 &&
        out.connectivity.person_extreme_flag[static_cast<size_t>(p)] == 0) {
      p_to_act[static_cast<size_t>(p)] = static_cast<int>(act_to_p.size());
      act_to_p.push_back(p);
    }
  }
  const int Pa = static_cast<int>(act_to_p.size());
  out.n_free_par = Pa + n_xi;

  std::vector<double> theta_act(static_cast<size_t>(Pa), 0.0);
  std::vector<double> xi(static_cast<size_t>(n_xi), 0.0);
  for (int idx_a = 0; idx_a < Pa; ++idx_a) {
    int p = act_to_p[static_cast<size_t>(idx_a)];
    double mn = out.connectivity.person_min_score[static_cast<size_t>(p)];
    double mx = out.connectivity.person_max_score[static_cast<size_t>(p)];
    double prop = ((out.connectivity.person_raw_score[static_cast<size_t>(p)] - mn) + 0.5) /
                  ((mx - mn) + 1.0);
    theta_act[static_cast<size_t>(idx_a)] = std::log(prop / (1.0 - prop));
  }
  if (static_cast<int>(ctrl.warm_init.size()) == Pa + n_xi) {
    std::copy(ctrl.warm_init.begin(), ctrl.warm_init.begin() + Pa, theta_act.begin());
    std::copy(ctrl.warm_init.begin() + Pa, ctrl.warm_init.end(), xi.begin());
  } else if (static_cast<int>(ctrl.warm_init.size()) == n_xi && n_xi > 0) {
    xi = ctrl.warm_init;
  }

  std::vector<std::vector<double>> beta;
  std::vector<double> tau_mat, step_cum_mat;
  std::vector<double> probs(static_cast<size_t>(C));

  std::vector<int> obs_f_ptr(static_cast<size_t>(spec.n_obs) + 1, 0);
  std::vector<int> obs_f_cols;
  std::vector<double> obs_f_vals;
  {
    std::vector<int> tmp_c(std::max(512, n_xi + 16));
    std::vector<double> tmp_v(std::max(512, n_xi + 16));
    for (int i = 0; i < spec.n_obs; ++i) {
      int nnz = obs_eta_xi_terms(i, spec, lderivs, tmp_c.data(), tmp_v.data());
      for (int h = 0; h < nnz; ++h) {
        obs_f_cols.push_back(tmp_c[static_cast<size_t>(h)]);
        obs_f_vals.push_back(tmp_v[static_cast<size_t>(h)]);
      }
      obs_f_ptr[static_cast<size_t>(i) + 1] = static_cast<int>(obs_f_cols.size());
    }
  }

  auto eval_nll_only = [&](const std::vector<double>& th_a, const std::vector<double>& xi_v) {
    smap.expand(spec, xi_v.data(), beta, tau_mat, step_cum_mat);
    double nll = 0.0;
    for (int i = 0; i < spec.n_obs; ++i) {
      double w = spec.weight[static_cast<size_t>(i)];
      if (w <= 0.0) continue;
      int pa = p_to_act[static_cast<size_t>(spec.person[static_cast<size_t>(i)])];
      if (pa < 0) continue;
      double eta = th_a[static_cast<size_t>(pa)];
      for (int f = 0; f < F; ++f) {
        int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
        eta += spec.facet_signs[static_cast<size_t>(f)] * beta[static_cast<size_t>(f)][static_cast<size_t>(l)];
      }
      int sg = spec.step_group[static_cast<size_t>(i)];
      CatMoments cm = eval_cat_probs_and_moments(
          eta, step_cum_mat.data() + static_cast<size_t>(sg) * C, C,
          spec.score[static_cast<size_t>(i)], probs.data(),
          smap.cat_min[static_cast<size_t>(sg)], smap.cat_max[static_cast<size_t>(sg)]);
      if (!std::isfinite(cm.log_py)) return kInf;
      nll -= w * cm.log_py;
    }
    return nll;
  };

  std::vector<double> g_th(static_cast<size_t>(Pa), 0.0), D_th(static_cast<size_t>(Pa), 0.0);
  std::vector<double> g_xi(static_cast<size_t>(n_xi), 0.0);
  std::vector<double> B_mat(static_cast<size_t>(Pa) * std::max(1, n_xi), 0.0);
  std::vector<double> C_mat(static_cast<size_t>(n_xi) * n_xi, 0.0);
  std::vector<double> z_xi(static_cast<size_t>(C) * std::max(1, n_xi), 0.0);
  std::vector<int> act_xi_cols;
  act_xi_cols.reserve(64);

  double cur_nll = eval_nll_only(theta_act, xi);
  double max_grad = 1.0;
  int iter_used = 0;

  for (int iter = 1; iter <= std::min(ctrl.maxit, 80); ++iter) {
    iter_used = iter;
    smap.expand(spec, xi.data(), beta, tau_mat, step_cum_mat);
    std::fill(g_th.begin(), g_th.end(), 0.0);
    std::fill(D_th.begin(), D_th.end(), 0.0);
    std::fill(g_xi.begin(), g_xi.end(), 0.0);
    std::fill(B_mat.begin(), B_mat.end(), 0.0);
    std::fill(C_mat.begin(), C_mat.end(), 0.0);
    cur_nll = 0.0;

    for (int i = 0; i < spec.n_obs; ++i) {
      double w = spec.weight[static_cast<size_t>(i)];
      if (w <= 0.0) continue;
      int pa = p_to_act[static_cast<size_t>(spec.person[static_cast<size_t>(i)])];
      if (pa < 0) continue;
      double eta = theta_act[static_cast<size_t>(pa)];
      for (int f = 0; f < F; ++f) {
        int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
        eta += spec.facet_signs[static_cast<size_t>(f)] * beta[static_cast<size_t>(f)][static_cast<size_t>(l)];
      }
      int sg = spec.step_group[static_cast<size_t>(i)];
      int ag = smap.cat_min[static_cast<size_t>(sg)], bg = smap.cat_max[static_cast<size_t>(sg)];
      int y = spec.score[static_cast<size_t>(i)];
      CatMoments cm = eval_cat_probs_and_moments(
          eta, step_cum_mat.data() + static_cast<size_t>(sg) * C, C, y, probs.data(), ag, bg);
      cur_nll -= w * cm.log_py;

      g_th[static_cast<size_t>(pa)] += w * (cm.mean - static_cast<double>(y));
      D_th[static_cast<size_t>(pa)] += w * cm.var;

      if (n_xi == 0) continue;
      act_xi_cols.clear();
      for (int ptr = obs_f_ptr[static_cast<size_t>(i)]; ptr < obs_f_ptr[static_cast<size_t>(i) + 1]; ++ptr) {
        int col = obs_f_cols[static_cast<size_t>(ptr)];
        double cf = obs_f_vals[static_cast<size_t>(ptr)];
        act_xi_cols.push_back(col);
        for (int k = ag; k <= bg; ++k) {
          z_xi[static_cast<size_t>(k) * n_xi + col] = (static_cast<double>(k) - cm.mean) * cf;
        }
      }
      populate_step_z_xi(smap, sg, C, probs.data(), act_xi_cols, z_xi.data());

      for (int c1 : act_xi_cols) {
        g_xi[static_cast<size_t>(c1)] -= w * z_xi[static_cast<size_t>(y) * n_xi + c1];
        double b_val = 0.0;
        for (int k = ag; k <= bg; ++k) {
          b_val += probs[static_cast<size_t>(k)] * (static_cast<double>(k) - cm.mean) *
                   z_xi[static_cast<size_t>(k) * n_xi + c1];
        }
        B_mat[static_cast<size_t>(pa) * n_xi + c1] += w * b_val;
        for (int c2 : act_xi_cols) {
          if (c2 > c1) continue;
          double c_val = 0.0;
          for (int k = ag; k <= bg; ++k) {
            c_val += probs[static_cast<size_t>(k)] *
                     z_xi[static_cast<size_t>(k) * n_xi + c1] *
                     z_xi[static_cast<size_t>(k) * n_xi + c2];
          }
          C_mat[static_cast<size_t>(c1) * n_xi + c2] += w * c_val;
        }
      }
    }
    for (int c1 = 0; c1 < n_xi; ++c1) {
      for (int c2 = 0; c2 < c1; ++c2) {
        C_mat[static_cast<size_t>(c2) * n_xi + c1] = C_mat[static_cast<size_t>(c1) * n_xi + c2];
      }
    }

    max_grad = 0.0;
    for (double v : g_th) max_grad = std::max(max_grad, std::abs(v));
    for (double v : g_xi) max_grad = std::max(max_grad, std::abs(v));
    if (max_grad <= std::max(1e-10, ctrl.reltol)) {
      out.converged = true;
      break;
    }

    std::vector<double> S = C_mat;
    std::vector<double> rhs = g_xi;
    for (int pa = 0; pa < Pa; ++pa) {
      double inv_d = 1.0 / std::max(D_th[static_cast<size_t>(pa)], 1e-10);
      double g_scaled = g_th[static_cast<size_t>(pa)] * inv_d;
      const double* b_row = B_mat.data() + static_cast<size_t>(pa) * n_xi;
      for (int c1 = 0; c1 < n_xi; ++c1) {
        if (b_row[c1] == 0.0) continue;
        rhs[static_cast<size_t>(c1)] -= b_row[c1] * g_scaled;
        double b1_d = b_row[c1] * inv_d;
        for (int c2 = 0; c2 <= c1; ++c2) {
          S[static_cast<size_t>(c1) * n_xi + c2] -= b1_d * b_row[c2];
        }
      }
    }
    for (int c1 = 0; c1 < n_xi; ++c1) {
      for (int c2 = 0; c2 < c1; ++c2) {
        S[static_cast<size_t>(c2) * n_xi + c1] = S[static_cast<size_t>(c1) * n_xi + c2];
      }
    }

    std::vector<double> d_xi(static_cast<size_t>(n_xi), 0.0);
    if (!damped_spd_solve(S, n_xi, rhs, d_xi, (iter <= 2) ? 1e-4 : 0.0)) break;

    std::vector<double> d_th(static_cast<size_t>(Pa), 0.0);
    for (int pa = 0; pa < Pa; ++pa) {
      double inv_d = 1.0 / std::max(D_th[static_cast<size_t>(pa)], 1e-10);
      double bd = 0.0;
      const double* b_row = B_mat.data() + static_cast<size_t>(pa) * n_xi;
      for (int c1 = 0; c1 < n_xi; ++c1) bd += b_row[c1] * d_xi[static_cast<size_t>(c1)];
      d_th[static_cast<size_t>(pa)] = inv_d * (g_th[static_cast<size_t>(pa)] - bd);
    }

    double mstep = 0.0;
    for (double v : d_th) mstep = std::max(mstep, std::abs(v));
    for (double v : d_xi) mstep = std::max(mstep, std::abs(v));
    double scale = (mstep > 3.0) ? (3.0 / mstep) : 1.0;

    std::vector<double> cand_th(static_cast<size_t>(Pa)), cand_xi(static_cast<size_t>(n_xi));
    bool step_ok = false;
    for (int ls = 0; ls < 25; ++ls) {
      for (int pa = 0; pa < Pa; ++pa) cand_th[static_cast<size_t>(pa)] = theta_act[static_cast<size_t>(pa)] - scale * d_th[static_cast<size_t>(pa)];
      for (int c1 = 0; c1 < n_xi; ++c1) cand_xi[static_cast<size_t>(c1)] = xi[static_cast<size_t>(c1)] - scale * d_xi[static_cast<size_t>(c1)];
      double cand_nll = eval_nll_only(cand_th, cand_xi);
      if (std::isfinite(cand_nll) && cand_nll <= cur_nll + 1e-12 * (1.0 + std::abs(cur_nll))) {
        theta_act = cand_th;
        xi = cand_xi;
        cur_nll = cand_nll;
        step_ok = true;
        break;
      }
      scale *= 0.5;
    }
    if (!step_ok) {
      if (max_grad <= 1e-6) out.converged = true;
      break;
    }
    if (mstep * scale < 1e-11 && max_grad <= 1e-7) {
      out.converged = true;
      break;
    }
  }

  out.iterations = iter_used;
  out.loglik = -cur_nll;
  out.objective = cur_nll;
  out.max_abs_grad = max_grad;
  out.xi_est = xi;

  smap.expand(spec, xi.data(), out.facet_est, out.step_est, out.step_cum);

  out.person_est.assign(static_cast<size_t>(P), kNaN);
  out.person_se.assign(static_cast<size_t>(P), kNaN);
  out.person_cond_se.assign(static_cast<size_t>(P), kNaN);
  std::vector<bool> diag_exclude(static_cast<size_t>(P), false);

  for (int pa = 0; pa < Pa; ++pa) {
    int p = act_to_p[static_cast<size_t>(pa)];
    out.person_est[static_cast<size_t>(p)] = theta_act[static_cast<size_t>(pa)];
  }
  std::vector<double> ext_se(static_cast<size_t>(P), kNaN);
  for (int p = 0; p < P; ++p) {
    if (out.connectivity.person_obs_count[static_cast<size_t>(p)] <= 0.0) {
      diag_exclude[static_cast<size_t>(p)] = true;
      continue;
    }
    int ext = out.connectivity.person_extreme_flag[static_cast<size_t>(p)];
    if (ext == 0) continue;
    diag_exclude[static_cast<size_t>(p)] = true;
    if (ctrl.xtreme > 0.0) {
      double target = (ext < 0) ? (out.connectivity.person_min_score[static_cast<size_t>(p)] + ctrl.xtreme)
                                : (out.connectivity.person_max_score[static_cast<size_t>(p)] - ctrl.xtreme);
      int beg = spec.p_row_ptr[static_cast<size_t>(p)];
      int end = spec.p_row_ptr[static_cast<size_t>(p) + 1];
      double th = (ext < 0) ? -4.0 : 4.0;
      double lo = -30.0, hi = 30.0;
      double v_fin = 0.0;
      for (int it = 0; it < 50; ++it) {
        double esum = 0.0, vsum = 0.0;
        for (int idx = beg; idx < end; ++idx) {
          int i = spec.p_obs_idx[static_cast<size_t>(idx)];
          double w = spec.weight[static_cast<size_t>(i)];
          if (w <= 0.0) continue;
          double eta = th;
          for (int f = 0; f < F; ++f) {
            int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
            eta += spec.facet_signs[static_cast<size_t>(f)] *
                   out.facet_est[static_cast<size_t>(f)][static_cast<size_t>(l)];
          }
          int sg = spec.step_group[static_cast<size_t>(i)];
          CatMoments cm = eval_cat_probs_and_moments(
              eta, out.step_cum.data() + static_cast<size_t>(sg) * C, C,
              spec.score[static_cast<size_t>(i)], probs.data(),
              smap.cat_min[static_cast<size_t>(sg)], smap.cat_max[static_cast<size_t>(sg)]);
          esum += w * cm.mean;
          vsum += w * cm.var;
        }
        v_fin = vsum;
        double diff = target - esum;
        if (diff > 0.0) lo = std::max(lo, th);
        else hi = std::min(hi, th);
        if (std::abs(diff) < 1e-11) break;
        double step = (vsum > 1e-12) ? (diff / vsum) : ((diff > 0.0) ? 1.0 : -1.0);
        double cand = th + std::max(-3.0, std::min(3.0, step));
        if (cand <= lo || cand >= hi) cand = 0.5 * (lo + hi);
        th = cand;
      }
      out.person_est[static_cast<size_t>(p)] = th;
      if (v_fin > 1e-12) ext_se[static_cast<size_t>(p)] = 1.0 / std::sqrt(v_fin);
    } else {
      out.person_est[static_cast<size_t>(p)] = (ext > 0) ? kInf : -kInf;
    }
  }

  std::vector<double> S = C_mat;
  for (int pa = 0; pa < Pa; ++pa) {
    double inv_d = 1.0 / std::max(D_th[static_cast<size_t>(pa)], 1e-12);
    const double* b_row = B_mat.data() + static_cast<size_t>(pa) * n_xi;
    for (int c1 = 0; c1 < n_xi; ++c1) {
      if (b_row[c1] == 0.0) continue;
      for (int c2 = 0; c2 <= c1; ++c2) {
        S[static_cast<size_t>(c1) * n_xi + c2] -= b_row[c1] * inv_d * b_row[c2];
      }
    }
  }
  for (int c1 = 0; c1 < n_xi; ++c1) {
    for (int c2 = 0; c2 < c1; ++c2) S[static_cast<size_t>(c2) * n_xi + c1] = S[static_cast<size_t>(c1) * n_xi + c2];
  }
  double logdet_S = 0.0;
  chol_invert_and_logdet(S, n_xi, out.xi_cov_model, logdet_S);
  std::vector<double> step_joint_se;
  smap.propagate_covariance(spec, out.xi_cov_model, out.facet_joint_se, out.facet_cov_model, step_joint_se, &out.step_cov_model);

  compute_fast_diagnostics_and_fair(
      spec, out.person_est, diag_exclude, out.facet_est, out.step_cum,
      ctrl.whexact, out.person_diag, out.facet_diag);

  out.person_real_se.assign(static_cast<size_t>(P), kNaN);
  for (int p = 0; p < P; ++p) {
    double ci = out.person_diag.cond_info[static_cast<size_t>(p)];
    if (!diag_exclude[static_cast<size_t>(p)] && ci > 1e-12) {
      double se = 1.0 / std::sqrt(ci);
      out.person_se[static_cast<size_t>(p)] = se;
      out.person_cond_se[static_cast<size_t>(p)] = se;
      double infit = out.person_diag.infit_mnsq[static_cast<size_t>(p)];
      out.person_real_se[static_cast<size_t>(p)] =
          std::isfinite(infit) ? (se * std::sqrt(std::max(1.0, infit))) : se;
    } else if (ctrl.xtreme > 0.0 && std::isfinite(ext_se[static_cast<size_t>(p)])) {
      out.person_se[static_cast<size_t>(p)] = ext_se[static_cast<size_t>(p)];
      out.person_cond_se[static_cast<size_t>(p)] = ext_se[static_cast<size_t>(p)];
      out.person_real_se[static_cast<size_t>(p)] = ext_se[static_cast<size_t>(p)];
    }
  }

  out.facet_se.resize(static_cast<size_t>(F));
  out.facet_real_se.resize(static_cast<size_t>(F));
  for (int f = 0; f < F; ++f) {
    int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
    bool custom_f = (!spec.has_custom_facet_map.empty() && spec.has_custom_facet_map[static_cast<size_t>(f)]);
    out.facet_se[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), kNaN);
    out.facet_real_se[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), kNaN);
    for (int l = 0; l < Lf; ++l) {
      if (smap.facet_fixed_mask[static_cast<size_t>(f)][static_cast<size_t>(l)]) {
        out.facet_se[static_cast<size_t>(f)][static_cast<size_t>(l)] = 0.0;
        out.facet_real_se[static_cast<size_t>(f)][static_cast<size_t>(l)] = 0.0;
      } else if (custom_f) {
        double se = out.facet_joint_se[static_cast<size_t>(f)][static_cast<size_t>(l)];
        out.facet_se[static_cast<size_t>(f)][static_cast<size_t>(l)] = se;
        double infit = out.facet_diag[static_cast<size_t>(f)].infit_mnsq[static_cast<size_t>(l)];
        out.facet_real_se[static_cast<size_t>(f)][static_cast<size_t>(l)] =
            (std::isfinite(se) && std::isfinite(infit)) ? (se * std::sqrt(std::max(1.0, infit))) : se;
      } else {
        double ci = out.facet_diag[static_cast<size_t>(f)].cond_info[static_cast<size_t>(l)];
        if (ci > 1e-12) {
          double se = 1.0 / std::sqrt(ci);
          out.facet_se[static_cast<size_t>(f)][static_cast<size_t>(l)] = se;
          double infit = out.facet_diag[static_cast<size_t>(f)].infit_mnsq[static_cast<size_t>(l)];
          out.facet_real_se[static_cast<size_t>(f)][static_cast<size_t>(l)] =
              std::isfinite(infit) ? (se * std::sqrt(std::max(1.0, infit))) : se;
        }
      }
    }
  }

  if (spec.has_custom_step_map) {
    out.step_se = step_joint_se;
  } else {
    out.step_se.assign(static_cast<size_t>(spec.n_step_groups) * K, kNaN);
    std::vector<double> step_info(static_cast<size_t>(spec.n_step_groups) * K, 0.0);
    for (int i = 0; i < spec.n_obs; ++i) {
      double w = spec.weight[static_cast<size_t>(i)];
      if (w <= 0.0) continue;
      int p = spec.person[static_cast<size_t>(i)];
      if (diag_exclude[static_cast<size_t>(p)] || !std::isfinite(out.person_est[static_cast<size_t>(p)])) continue;
      double eta = out.person_est[static_cast<size_t>(p)];
      for (int f = 0; f < F; ++f) {
        int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
        eta += spec.facet_signs[static_cast<size_t>(f)] * out.facet_est[static_cast<size_t>(f)][static_cast<size_t>(l)];
      }
      int sg = spec.step_group[static_cast<size_t>(i)];
      int ag = smap.cat_min[static_cast<size_t>(sg)], bg = smap.cat_max[static_cast<size_t>(sg)];
      eval_cat_probs_and_moments(eta, out.step_cum.data() + static_cast<size_t>(sg) * C, C, 0, probs.data(), ag, bg);
      double p_geq = 0.0;
      for (int k = bg; k >= ag + 1; --k) {
        p_geq += probs[static_cast<size_t>(k)];
        step_info[static_cast<size_t>(sg) * K + (k - 1)] += w * p_geq * (1.0 - p_geq);
      }
    }
    for (size_t idx = 0; idx < step_info.size(); ++idx) {
      if (!smap.step_fixed[idx] && step_info[idx] > 1e-12) {
        out.step_se[idx] = 1.0 / std::sqrt(step_info[idx]);
      }
    }
  }

  if (ctrl.bias_facet_a >= 0 && ctrl.bias_facet_b >= 0) {
    out.local_bias_cells = compute_local_bias_2way(
        spec, out.person_est, diag_exclude, out.facet_est, out.step_cum,
        ctrl.bias_facet_a, ctrl.bias_facet_b);
  }
}

// ============================================================================
// Solver 2: Ultra-Fast Block-Schur Firth / Jeffreys JML (with Exact Schur SE!)
// ============================================================================

inline void solve_jeffreys(const FastDesignSpec& spec, const StructuralMap& smap,
                           const FastControl& ctrl, FastFitOutput& out) {
  const int P = spec.n_person;
  const int F = spec.n_facets;
  const int C = spec.n_cat;
  const int n_xi = smap.n_xi;
  const auto lderivs = smap.build_level_derivs(spec);

  std::vector<int> p_to_act(static_cast<size_t>(P), -1);
  std::vector<int> act_to_p;
  for (int p = 0; p < P; ++p) {
    if (out.connectivity.person_obs_count[static_cast<size_t>(p)] > 0.0) {
      p_to_act[static_cast<size_t>(p)] = static_cast<int>(act_to_p.size());
      act_to_p.push_back(p);
    }
  }
  const int Pa = static_cast<int>(act_to_p.size());
  const int dim = Pa + n_xi;
  out.n_free_par = dim;

  std::vector<int> obs_f_ptr(static_cast<size_t>(spec.n_obs) + 1, 0);
  std::vector<int> obs_f_cols;
  std::vector<double> obs_f_vals;
  {
    std::vector<int> tmp_c(std::max(512, n_xi + 16));
    std::vector<double> tmp_v(std::max(512, n_xi + 16));
    for (int i = 0; i < spec.n_obs; ++i) {
      int nnz = obs_eta_xi_terms(i, spec, lderivs, tmp_c.data(), tmp_v.data());
      for (int h = 0; h < nnz; ++h) {
        obs_f_cols.push_back(tmp_c[static_cast<size_t>(h)]);
        obs_f_vals.push_back(tmp_v[static_cast<size_t>(h)]);
      }
      obs_f_ptr[static_cast<size_t>(i) + 1] = static_cast<int>(obs_f_cols.size());
    }
  }

  std::vector<std::vector<double>> beta;
  std::vector<double> tau_mat, step_cum_mat;
  std::vector<double> probs(static_cast<size_t>(C));
  std::vector<double> D_th(static_cast<size_t>(Pa), 0.0);
  std::vector<double> B_mat(static_cast<size_t>(Pa) * std::max(1, n_xi), 0.0);
  std::vector<double> C_mat(static_cast<size_t>(n_xi) * n_xi, 0.0);
  std::vector<double> invS, alpha_p(static_cast<size_t>(Pa), 0.0);
  std::vector<double> w_p(static_cast<size_t>(Pa) * std::max(1, n_xi), 0.0);
  std::vector<double> z_xi(static_cast<size_t>(C) * std::max(1, n_xi), 0.0);
  std::vector<int> act_cols;
  act_cols.reserve(64);
  double last_nll = kNaN;

  auto eval_jeffreys_fg = [&](const std::vector<double>& x, std::vector<double>& grad) -> double {
    grad.assign(static_cast<size_t>(dim), 0.0);
    smap.expand(spec, x.data() + Pa, beta, tau_mat, step_cum_mat);
    std::fill(D_th.begin(), D_th.end(), 0.0);
    std::fill(B_mat.begin(), B_mat.end(), 0.0);
    std::fill(C_mat.begin(), C_mat.end(), 0.0);
    double nll = 0.0;

    for (int i = 0; i < spec.n_obs; ++i) {
      double w = spec.weight[static_cast<size_t>(i)];
      if (w <= 0.0) continue;
      int pa = p_to_act[static_cast<size_t>(spec.person[static_cast<size_t>(i)])];
      if (pa < 0) continue;
      double eta = x[static_cast<size_t>(pa)];
      if (!std::isfinite(eta) || std::abs(eta) > 35.0) return kInf;
      for (int f = 0; f < F; ++f) {
        int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
        eta += spec.facet_signs[static_cast<size_t>(f)] * beta[static_cast<size_t>(f)][static_cast<size_t>(l)];
      }
      if (!std::isfinite(eta) || std::abs(eta) > 45.0) return kInf;
      int sg = spec.step_group[static_cast<size_t>(i)];
      int ag = smap.cat_min[static_cast<size_t>(sg)], bg = smap.cat_max[static_cast<size_t>(sg)];
      int y = spec.score[static_cast<size_t>(i)];
      CatMoments cm = eval_cat_probs_and_moments(
          eta, step_cum_mat.data() + static_cast<size_t>(sg) * C, C, y, probs.data(), ag, bg);
      if (!std::isfinite(cm.log_py) || cm.var <= 0.0) return kInf;
      nll -= w * cm.log_py;

      grad[static_cast<size_t>(pa)] += w * (cm.mean - static_cast<double>(y));
      D_th[static_cast<size_t>(pa)] += w * cm.var;

      if (n_xi == 0) continue;
      act_cols.clear();
      for (int ptr = obs_f_ptr[static_cast<size_t>(i)]; ptr < obs_f_ptr[static_cast<size_t>(i) + 1]; ++ptr) {
        int col = obs_f_cols[static_cast<size_t>(ptr)];
        double cf = obs_f_vals[static_cast<size_t>(ptr)];
        act_cols.push_back(col);
        for (int k = ag; k <= bg; ++k) z_xi[static_cast<size_t>(k) * n_xi + col] = (static_cast<double>(k) - cm.mean) * cf;
      }
      populate_step_z_xi(smap, sg, C, probs.data(), act_cols, z_xi.data());

      for (int c1 : act_cols) {
        grad[static_cast<size_t>(Pa + c1)] -= w * z_xi[static_cast<size_t>(y) * n_xi + c1];
        double b_val = 0.0;
        for (int k = ag; k <= bg; ++k) {
          b_val += probs[static_cast<size_t>(k)] * (static_cast<double>(k) - cm.mean) * z_xi[static_cast<size_t>(k) * n_xi + c1];
        }
        B_mat[static_cast<size_t>(pa) * n_xi + c1] += w * b_val;
        for (int c2 : act_cols) {
          if (c2 > c1) continue;
          double c_val = 0.0;
          for (int k = ag; k <= bg; ++k) {
            c_val += probs[static_cast<size_t>(k)] * z_xi[static_cast<size_t>(k) * n_xi + c1] * z_xi[static_cast<size_t>(k) * n_xi + c2];
          }
          C_mat[static_cast<size_t>(c1) * n_xi + c2] += w * c_val;
        }
      }
    }
    for (int c1 = 0; c1 < n_xi; ++c1) {
      for (int c2 = 0; c2 < c1; ++c2) C_mat[static_cast<size_t>(c2) * n_xi + c1] = C_mat[static_cast<size_t>(c1) * n_xi + c2];
    }

    double logdet_D = 0.0;
    std::vector<double> S = C_mat;
    for (int pa = 0; pa < Pa; ++pa) {
      double d = D_th[static_cast<size_t>(pa)];
      if (!std::isfinite(d) || d <= 1e-15) return kInf;
      logdet_D += std::log(d);
      double inv_d = 1.0 / d;
      const double* b_row = B_mat.data() + static_cast<size_t>(pa) * n_xi;
      for (int c1 = 0; c1 < n_xi; ++c1) {
        if (b_row[c1] == 0.0) continue;
        for (int c2 = 0; c2 <= c1; ++c2) {
          S[static_cast<size_t>(c1) * n_xi + c2] -= b_row[c1] * inv_d * b_row[c2];
        }
      }
    }
    for (int c1 = 0; c1 < n_xi; ++c1) {
      for (int c2 = 0; c2 < c1; ++c2) S[static_cast<size_t>(c2) * n_xi + c1] = S[static_cast<size_t>(c1) * n_xi + c2];
    }

    double logdet_S = 0.0;
    if (!chol_invert_and_logdet(S, n_xi, invS, logdet_S)) return kInf;
    double obj = nll - 0.5 * (logdet_D + logdet_S);
    last_nll = nll;

    for (int pa = 0; pa < Pa; ++pa) {
      double inv_d = 1.0 / D_th[static_cast<size_t>(pa)];
      const double* b_row = B_mat.data() + static_cast<size_t>(pa) * n_xi;
      double quad = 0.0;
      for (int c1 = 0; c1 < n_xi; ++c1) {
        double acc = 0.0;
        for (int c2 = 0; c2 < n_xi; ++c2) {
          acc += invS[static_cast<size_t>(c1) * n_xi + c2] * (b_row[c2] * inv_d);
        }
        w_p[static_cast<size_t>(pa) * n_xi + c1] = acc;
        quad += (b_row[c1] * inv_d) * acc;
      }
      alpha_p[static_cast<size_t>(pa)] = inv_d + quad;
    }

    for (int i = 0; i < spec.n_obs; ++i) {
      double w = spec.weight[static_cast<size_t>(i)];
      if (w <= 0.0) continue;
      int pa = p_to_act[static_cast<size_t>(spec.person[static_cast<size_t>(i)])];
      if (pa < 0) continue;
      double eta = x[static_cast<size_t>(pa)];
      for (int f = 0; f < F; ++f) {
        int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
        eta += spec.facet_signs[static_cast<size_t>(f)] * beta[static_cast<size_t>(f)][static_cast<size_t>(l)];
      }
      int sg = spec.step_group[static_cast<size_t>(i)];
      int ag = smap.cat_min[static_cast<size_t>(sg)], bg = smap.cat_max[static_cast<size_t>(sg)];
      CatMoments cm = eval_cat_probs_and_moments(
          eta, step_cum_mat.data() + static_cast<size_t>(sg) * C, C, 0, probs.data(), ag, bg);

      act_cols.clear();
      if (n_xi > 0) {
        for (int ptr = obs_f_ptr[static_cast<size_t>(i)]; ptr < obs_f_ptr[static_cast<size_t>(i) + 1]; ++ptr) {
          int col = obs_f_cols[static_cast<size_t>(ptr)];
          double cf = obs_f_vals[static_cast<size_t>(ptr)];
          act_cols.push_back(col);
          for (int k = ag; k <= bg; ++k) z_xi[static_cast<size_t>(k) * n_xi + col] = (static_cast<double>(k) - cm.mean) * cf;
        }
        populate_step_z_xi(smap, sg, C, probs.data(), act_cols, z_xi.data());
      }

      double al_p = alpha_p[static_cast<size_t>(pa)];
      const double* wp_row = w_p.data() + static_cast<size_t>(pa) * n_xi;
      for (int k = ag; k <= bg; ++k) {
        double pk = probs[static_cast<size_t>(k)];
        if (pk <= 0.0) continue;
        double u = static_cast<double>(k) - cm.mean;
        double wv = 0.0, vSv = 0.0;
        for (int c1 : act_cols) {
          double v1 = z_xi[static_cast<size_t>(k) * n_xi + c1];
          wv += wp_row[c1] * v1;
          for (int c2 : act_cols) {
            vSv += v1 * invS[static_cast<size_t>(c1) * n_xi + c2] * z_xi[static_cast<size_t>(k) * n_xi + c2];
          }
        }
        double q_ik = u * u * al_p - 2.0 * u * wv + vSv;
        double factor = 0.5 * w * pk * q_ik;
        grad[static_cast<size_t>(pa)] -= factor * u;
        for (int c1 : act_cols) {
          grad[static_cast<size_t>(Pa + c1)] -= factor * z_xi[static_cast<size_t>(k) * n_xi + c1];
        }
      }
    }
    return obj;
  };

  std::vector<double> x0(static_cast<size_t>(dim), 0.0);
  for (int pa = 0; pa < Pa; ++pa) {
    int p = act_to_p[static_cast<size_t>(pa)];
    double mn = out.connectivity.person_min_score[static_cast<size_t>(p)];
    double mx = out.connectivity.person_max_score[static_cast<size_t>(p)];
    double prop = ((out.connectivity.person_raw_score[static_cast<size_t>(p)] - mn) + 0.5) /
                  ((mx - mn) + 1.0);
    x0[static_cast<size_t>(pa)] = std::max(-3.5, std::min(3.5, std::log(prop / (1.0 - prop))));
  }
  if (static_cast<int>(ctrl.warm_init.size()) == dim) {
    x0 = ctrl.warm_init;
  }
  LbfgsResult opt = minimize_lbfgs(x0, eval_jeffreys_fg, std::max(ctrl.maxit, 250), 1e-7, 1e-13, 16);

  std::vector<double> g_cur;
  double f_cur = eval_jeffreys_fg(opt.x, g_cur);
  auto max_abs_v = [](const std::vector<double>& v) {
    double m = 0.0;
    for (double val : v) m = std::max(m, std::abs(val));
    return m;
  };
  for (int pol = 0; pol < 15 && max_abs_v(g_cur) > 1e-6; ++pol) {
    std::vector<double> rhs_xi(static_cast<size_t>(n_xi), 0.0);
    for (int c1 = 0; c1 < n_xi; ++c1) rhs_xi[static_cast<size_t>(c1)] = g_cur[static_cast<size_t>(Pa + c1)];
    for (int pa = 0; pa < Pa; ++pa) {
      double inv_d = 1.0 / D_th[static_cast<size_t>(pa)];
      const double* b_row = B_mat.data() + static_cast<size_t>(pa) * n_xi;
      for (int c1 = 0; c1 < n_xi; ++c1) rhs_xi[static_cast<size_t>(c1)] -= b_row[c1] * inv_d * g_cur[static_cast<size_t>(pa)];
    }
    std::vector<double> d_xi(static_cast<size_t>(n_xi), 0.0);
    for (int c1 = 0; c1 < n_xi; ++c1) {
      double s = 0.0;
      for (int c2 = 0; c2 < n_xi; ++c2) s += invS[static_cast<size_t>(c1) * n_xi + c2] * rhs_xi[static_cast<size_t>(c2)];
      d_xi[static_cast<size_t>(c1)] = s;
    }
    std::vector<double> d_th(static_cast<size_t>(Pa), 0.0);
    for (int pa = 0; pa < Pa; ++pa) {
      double inv_d = 1.0 / D_th[static_cast<size_t>(pa)];
      const double* b_row = B_mat.data() + static_cast<size_t>(pa) * n_xi;
      double bd = 0.0;
      for (int c1 = 0; c1 < n_xi; ++c1) bd += b_row[c1] * d_xi[static_cast<size_t>(c1)];
      d_th[static_cast<size_t>(pa)] = inv_d * (g_cur[static_cast<size_t>(pa)] - bd);
    }
    double scale = 1.0;
    bool improved = false;
    std::vector<double> x_try(static_cast<size_t>(dim)), g_try;
    for (int ls = 0; ls < 12; ++ls) {
      for (int pa = 0; pa < Pa; ++pa) x_try[static_cast<size_t>(pa)] = opt.x[static_cast<size_t>(pa)] - scale * d_th[static_cast<size_t>(pa)];
      for (int c1 = 0; c1 < n_xi; ++c1) x_try[static_cast<size_t>(Pa + c1)] = opt.x[static_cast<size_t>(Pa + c1)] - scale * d_xi[static_cast<size_t>(c1)];
      double f_try = eval_jeffreys_fg(x_try, g_try);
      if (std::isfinite(f_try) && (f_try < f_cur || max_abs_v(g_try) < max_abs_v(g_cur))) {
        opt.x = x_try;
        f_cur = f_try;
        g_cur = g_try;
        improved = true;
        break;
      }
      scale *= 0.5;
    }
    if (!improved) break;
  }

  f_cur = eval_jeffreys_fg(opt.x, g_cur);
  double gmax = max_abs_v(g_cur);
  out.converged = std::isfinite(f_cur) && (gmax <= 1e-4);
  out.iterations = opt.iterations;
  out.objective = f_cur;
  out.loglik = -last_nll;
  out.max_abs_grad = gmax;
  out.xi_est.assign(opt.x.begin() + Pa, opt.x.end());
  out.xi_cov_model = invS;

  smap.expand(spec, opt.x.data() + Pa, out.facet_est, out.step_est, out.step_cum);

  out.person_est.assign(static_cast<size_t>(P), kNaN);
  out.person_se.assign(static_cast<size_t>(P), kNaN);
  out.person_cond_se.assign(static_cast<size_t>(P), kNaN);
  std::vector<bool> diag_exclude(static_cast<size_t>(P), false);
  for (int p = 0; p < P; ++p) {
    int pa = p_to_act[static_cast<size_t>(p)];
    if (pa < 0) {
      diag_exclude[static_cast<size_t>(p)] = true;
    } else {
      out.person_est[static_cast<size_t>(p)] = opt.x[static_cast<size_t>(pa)];
      double v_joint = alpha_p[static_cast<size_t>(pa)];
      out.person_se[static_cast<size_t>(p)] = (v_joint > 0.0) ? std::sqrt(v_joint) : kNaN;
      out.person_cond_se[static_cast<size_t>(p)] = (D_th[static_cast<size_t>(pa)] > 0.0) ? (1.0 / std::sqrt(D_th[static_cast<size_t>(pa)])) : kNaN;
    }
  }

  smap.propagate_covariance(spec, invS, out.facet_se, out.facet_cov_model, out.step_se, &out.step_cov_model);
  out.facet_joint_se = out.facet_se;

  compute_fast_diagnostics_and_fair(
      spec, out.person_est, diag_exclude, out.facet_est, out.step_cum,
      ctrl.whexact, out.person_diag, out.facet_diag);

  out.person_real_se.assign(static_cast<size_t>(P), kNaN);
  for (int p = 0; p < P; ++p) {
    double se = out.person_se[static_cast<size_t>(p)];
    if (std::isfinite(se)) {
      double infit = out.person_diag.infit_mnsq[static_cast<size_t>(p)];
      out.person_real_se[static_cast<size_t>(p)] = std::isfinite(infit) ? (se * std::sqrt(std::max(1.0, infit))) : se;
    }
  }
  out.facet_real_se.resize(static_cast<size_t>(F));
  for (int f = 0; f < F; ++f) {
    int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
    out.facet_real_se[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), kNaN);
    for (int l = 0; l < Lf; ++l) {
      double se = out.facet_se[static_cast<size_t>(f)][static_cast<size_t>(l)];
      if (std::isfinite(se)) {
        double infit = out.facet_diag[static_cast<size_t>(f)].infit_mnsq[static_cast<size_t>(l)];
        out.facet_real_se[static_cast<size_t>(f)][static_cast<size_t>(l)] =
            std::isfinite(infit) ? (se * std::sqrt(std::max(1.0, infit))) : se;
      }
    }
  }
  if (ctrl.bias_facet_a >= 0 && ctrl.bias_facet_b >= 0) {
    out.local_bias_cells = compute_local_bias_2way(
        spec, out.person_est, diag_exclude, out.facet_est, out.step_cum,
        ctrl.bias_facet_a, ctrl.bias_facet_b);
  }
}

// ============================================================================
// Solver 3: Ultra-Fast C++ MML (Fixed / Adaptive GH, Latent Reg, EAP + WLE,
//           Observed-Info & Sandwich SEs, Affine Facets/Steps, Category Support)
// ============================================================================

inline void solve_mml(const FastDesignSpec& spec, const StructuralMap& smap,
                      const FastControl& ctrl, FastFitOutput& out) {
  const int P = spec.n_person;
  const int F = spec.n_facets;
  const int C = spec.n_cat;
  const int n_xi = smap.n_xi;
  const int qs = smap.step_xi_len;
  const int Q = std::max(5, ctrl.quad_points);
  const int dX = (spec.pop_mode == 1) ? spec.pop_dim : 0;
  const int has_free_sigma = (spec.pop_mode == 1) ? 1 : 0;
  const int dim = n_xi + dX + has_free_sigma;
  const bool has_slopes = (smap.n_slope_facets > 0);
  const int sform = spec.slope_formula;
  out.n_free_par = dim;

  std::vector<double> z_nodes, w_nodes, log_w_nodes;
  gauss_hermite_normal(Q, z_nodes, w_nodes, log_w_nodes);
  const auto lderivs = smap.build_level_derivs(spec);
  const auto sderivs = smap.build_slope_derivs(spec);

  std::vector<int> obs_f_ptr(static_cast<size_t>(spec.n_obs) + 1, 0);
  std::vector<int> obs_f_cols;
  std::vector<double> obs_f_vals;
  std::vector<int> obs_s_ptr(static_cast<size_t>(spec.n_obs) + 1, 0);
  std::vector<int> obs_s_cols;
  std::vector<double> obs_s_vals;
  {
    std::vector<int> tmp_c(std::max(512, n_xi + 16));
    std::vector<double> tmp_v(std::max(512, n_xi + 16));
    for (int i = 0; i < spec.n_obs; ++i) {
      int nnz = obs_eta_xi_terms(i, spec, lderivs, tmp_c.data(), tmp_v.data());
      for (int h = 0; h < nnz; ++h) {
        obs_f_cols.push_back(tmp_c[static_cast<size_t>(h)]);
        obs_f_vals.push_back(tmp_v[static_cast<size_t>(h)]);
      }
      obs_f_ptr[static_cast<size_t>(i) + 1] = static_cast<int>(obs_f_cols.size());

      if (has_slopes) {
        int snz = obs_log_slope_xi_terms(i, spec, sderivs, tmp_c.data(), tmp_v.data());
        for (int h = 0; h < snz; ++h) {
          obs_s_cols.push_back(tmp_c[static_cast<size_t>(h)]);
          obs_s_vals.push_back(tmp_v[static_cast<size_t>(h)]);
        }
      }
      obs_s_ptr[static_cast<size_t>(i) + 1] = static_cast<int>(obs_s_cols.size());
    }
  }

  int max_group_sc = 0;
  for (int g = 0; g < smap.n_step_groups; ++g) {
    max_group_sc = std::max(max_group_sc, static_cast<int>(smap.group_step_cols[static_cast<size_t>(g)].size()));
  }

  std::vector<std::vector<double>> beta, log_slope_est, slope_est;
  std::vector<double> tau_mat, step_cum_mat, alpha_obs(static_cast<size_t>(spec.n_obs), 1.0);
  std::vector<double> base_eta_obs(static_cast<size_t>(spec.n_obs), 0.0);
  std::vector<double> mu_p(static_cast<size_t>(P), 0.0);

  bool adaptive_active = false;
  std::vector<double> adap_mode(static_cast<size_t>(P), 0.0);
  std::vector<double> adap_scale(static_cast<size_t>(P), 1.0);

  auto unpack_pop = [&](const std::vector<double>& phi, double& sigma_out) {
    if (spec.pop_mode == 0) {
      std::fill(mu_p.begin(), mu_p.end(), 0.0);
      sigma_out = spec.fixed_sigma;
    } else {
      for (int p = 0; p < P; ++p) {
        double m = 0.0;
        for (int c = 0; c < dX; ++c) {
          m += spec.pop_X[static_cast<size_t>(p) * dX + c] * phi[static_cast<size_t>(n_xi + c)];
        }
        mu_p[static_cast<size_t>(p)] = m;
      }
      double log_sig = phi[static_cast<size_t>(n_xi + dX)];
      if (log_sig < -12.0 || log_sig > 6.0) sigma_out = kNaN;
      else sigma_out = std::exp(log_sig);
    }
  };

  auto refresh_adaptive_basis = [&](const std::vector<double>& phi) {
    double sigma = 1.0;
    unpack_pop(phi, sigma);
    if (!std::isfinite(sigma) || sigma <= 0.0) return;
    smap.expand(spec, phi.data(), beta, tau_mat, step_cum_mat);
    if (!smap.expand_slopes(spec, phi.data(), log_slope_est, slope_est, alpha_obs)) return;
    for (int i = 0; i < spec.n_obs; ++i) {
      double e = 0.0;
      for (int f = 0; f < F; ++f) {
        int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
        e += spec.facet_signs[static_cast<size_t>(f)] * beta[static_cast<size_t>(f)][static_cast<size_t>(l)];
      }
      base_eta_obs[static_cast<size_t>(i)] = e;
    }
    std::vector<double> probs(static_cast<size_t>(C));
    double inv_sig2 = 1.0 / (sigma * sigma);
    for (int p = 0; p < P; ++p) {
      int beg = spec.p_row_ptr[static_cast<size_t>(p)];
      int end = spec.p_row_ptr[static_cast<size_t>(p) + 1];
      double th = mu_p[static_cast<size_t>(p)];
      if (beg == end) {
        adap_mode[static_cast<size_t>(p)] = th;
        adap_scale[static_cast<size_t>(p)] = sigma;
        continue;
      }
      for (int it = 0; it < 20; ++it) {
        double g = -(th - mu_p[static_cast<size_t>(p)]) * inv_sig2;
        double h = inv_sig2;
        for (int idx = beg; idx < end; ++idx) {
          int i = spec.p_obs_idx[static_cast<size_t>(idx)];
          double w = spec.weight[static_cast<size_t>(i)];
          if (w <= 0.0) continue;
          int sg = spec.step_group[static_cast<size_t>(i)];
          double al = alpha_obs[static_cast<size_t>(i)];
          CatMoments cm = eval_gpcm_probs_and_moments(
              th, base_eta_obs[static_cast<size_t>(i)], al, sform,
              step_cum_mat.data() + static_cast<size_t>(sg) * C, C,
              spec.score[static_cast<size_t>(i)], probs.data(),
              smap.cat_min[static_cast<size_t>(sg)], smap.cat_max[static_cast<size_t>(sg)]);
          g += w * al * (static_cast<double>(spec.score[static_cast<size_t>(i)]) - cm.mean);
          h += w * al * al * cm.var;
        }
        double step = std::max(-2.0, std::min(2.0, g / std::max(h, 1e-8)));
        th += step;
        if (std::abs(step) < 1e-10) break;
      }
      double h = inv_sig2;
      for (int idx = beg; idx < end; ++idx) {
        int i = spec.p_obs_idx[static_cast<size_t>(idx)];
        double w = spec.weight[static_cast<size_t>(i)];
        if (w <= 0.0) continue;
        int sg = spec.step_group[static_cast<size_t>(i)];
        int ag = smap.cat_min[static_cast<size_t>(sg)], bg = smap.cat_max[static_cast<size_t>(sg)];
        double al = alpha_obs[static_cast<size_t>(i)];
        CatMoments cm = eval_gpcm_probs_and_moments(
            th, base_eta_obs[static_cast<size_t>(i)], al, sform,
            step_cum_mat.data() + static_cast<size_t>(sg) * C, C, ag, probs.data(), ag, bg);
        h += w * al * al * cm.var;
      }
      adap_mode[static_cast<size_t>(p)] = th;
      adap_scale[static_cast<size_t>(p)] = 1.0 / std::sqrt(std::max(h, 1e-8));
    }
  };

  auto eval_mml_core = [&](const std::vector<double>& phi,
                           std::vector<double>& grad,
                           std::vector<double>* person_scores_out,
                           std::vector<double>* eap_out,
                           std::vector<double>* psd_out) -> double {
    grad.assign(static_cast<size_t>(dim), 0.0);
    if (person_scores_out) person_scores_out->assign(static_cast<size_t>(P) * dim, 0.0);
    double sigma = 1.0;
    unpack_pop(phi, sigma);
    if (!std::isfinite(sigma) || sigma <= 1e-8) return kInf;
    smap.expand(spec, phi.data(), beta, tau_mat, step_cum_mat);
    if (!smap.expand_slopes(spec, phi.data(), log_slope_est, slope_est, alpha_obs)) return kInf;
    for (int i = 0; i < spec.n_obs; ++i) {
      double e = 0.0;
      for (int f = 0; f < F; ++f) {
        int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
        e += spec.facet_signs[static_cast<size_t>(f)] * beta[static_cast<size_t>(f)][static_cast<size_t>(l)];
      }
      if (!std::isfinite(e) || std::abs(e) > 40.0) return kInf;
      base_eta_obs[static_cast<size_t>(i)] = e;
    }

    std::vector<double> probs(static_cast<size_t>(C)), step_ea_acc(static_cast<size_t>(std::max(1, max_group_sc)));
    std::vector<double> log_post(static_cast<size_t>(Q)), W_q(static_cast<size_t>(Q)), th_q(static_cast<size_t>(Q));
    std::vector<double> dlogL_dth(static_cast<size_t>(Q), 0.0);
    std::vector<double> p_grad_local(static_cast<size_t>(dim), 0.0);
    double total_nll = 0.0;
    double inv_sig2 = 1.0 / (sigma * sigma);

    if (eap_out) eap_out->assign(static_cast<size_t>(P), kNaN);
    if (psd_out) psd_out->assign(static_cast<size_t>(P), kNaN);

    for (int p = 0; p < P; ++p) {
      int beg = spec.p_row_ptr[static_cast<size_t>(p)];
      int end = spec.p_row_ptr[static_cast<size_t>(p) + 1];
      if (beg == end) {
        if (eap_out) (*eap_out)[static_cast<size_t>(p)] = mu_p[static_cast<size_t>(p)];
        if (psd_out) (*psd_out)[static_cast<size_t>(p)] = sigma;
        continue;
      }

      double mp = mu_p[static_cast<size_t>(p)];
      for (int q = 0; q < Q; ++q) {
        double z = z_nodes[static_cast<size_t>(q)];
        if (!adaptive_active) {
          th_q[static_cast<size_t>(q)] = mp + sigma * z;
          log_post[static_cast<size_t>(q)] = log_w_nodes[static_cast<size_t>(q)];
        } else {
          double sp = adap_scale[static_cast<size_t>(p)];
          double tq = adap_mode[static_cast<size_t>(p)] + sp * z;
          th_q[static_cast<size_t>(q)] = tq;
          double diff = tq - mp;
          log_post[static_cast<size_t>(q)] =
              log_w_nodes[static_cast<size_t>(q)] + 0.5 * z * z +
              std::log(sp / sigma) - 0.5 * diff * diff * inv_sig2;
        }
        dlogL_dth[static_cast<size_t>(q)] = 0.0;
      }

      const int m_obs_p = end - beg;
      const int stride_sc = std::max(1, max_group_sc);
      std::vector<double> cached_mean(static_cast<size_t>(m_obs_p) * Q);
      std::vector<double> cached_step_ea(static_cast<size_t>(m_obs_p) * Q * stride_sc, 0.0);
      std::vector<double> cached_logslope_g(has_slopes ? (static_cast<size_t>(m_obs_p) * Q) : 0, 0.0);

      for (int loc = 0; loc < m_obs_p; ++loc) {
        int i = spec.p_obs_idx[static_cast<size_t>(beg + loc)];
        double w = spec.weight[static_cast<size_t>(i)];
        if (w <= 0.0) continue;
        int sg = spec.step_group[static_cast<size_t>(i)];
        int ag = smap.cat_min[static_cast<size_t>(sg)], bg = smap.cat_max[static_cast<size_t>(sg)];
        int y = spec.score[static_cast<size_t>(i)];
        const double* sc = step_cum_mat.data() + static_cast<size_t>(sg) * C;
        const auto& scols = smap.group_step_cols[static_cast<size_t>(sg)];
        const int n_sc = static_cast<int>(scols.size());
        double be = base_eta_obs[static_cast<size_t>(i)];
        double al = alpha_obs[static_cast<size_t>(i)];

        for (int q = 0; q < Q; ++q) {
          double tq = th_q[static_cast<size_t>(q)];
          CatMoments cm = eval_gpcm_probs_and_moments(tq, be, al, sform, sc, C, y, probs.data(), ag, bg);
          log_post[static_cast<size_t>(q)] += w * cm.log_py;
          dlogL_dth[static_cast<size_t>(q)] += w * al * (static_cast<double>(y) - cm.mean);
          cached_mean[static_cast<size_t>(loc) * Q + q] = cm.mean;
          if (has_slopes) {
            if (sform == 1) {
              // u_k = (k - ag)*(alpha * theta + be) - S_k  =>  d(-log P_y)/d(log alpha) = alpha * theta * (mean - y)
              cached_logslope_g[static_cast<size_t>(loc) * Q + q] = al * tq * (cm.mean - static_cast<double>(y));
            } else {
              // u_k = alpha * ((k - ag)*(theta + be) - S_k) => d(-log P_y)/d(log alpha) = E[u] - u_y
              double e_sc = 0.0;
              for (int k = ag; k <= bg; ++k) {
                if (probs[static_cast<size_t>(k)] > 0.0 && std::isfinite(sc[k])) {
                  e_sc += probs[static_cast<size_t>(k)] * sc[k];
                }
              }
              cached_logslope_g[static_cast<size_t>(loc) * Q + q] =
                  al * ((cm.mean - static_cast<double>(y)) * (tq + be) + (sc[y] - e_sc));
            }
          }
          if (n_sc > 0) {
            double* ea_ptr = cached_step_ea.data() + (static_cast<size_t>(loc) * Q + q) * stride_sc;
            for (int m_idx = 0; m_idx < n_sc; ++m_idx) {
              int c_loc = scols[static_cast<size_t>(m_idx)];
              double ea = 0.0;
              for (int k = ag; k <= bg; ++k) {
                ea += probs[static_cast<size_t>(k)] * smap.cum_B[(static_cast<size_t>(sg) * C + k) * qs + c_loc];
              }
              ea_ptr[m_idx] = ea;
            }
          }
        }
      }

      double max_lp = log_post[0];
      for (int q = 1; q < Q; ++q) if (log_post[static_cast<size_t>(q)] > max_lp) max_lp = log_post[static_cast<size_t>(q)];
      if (!std::isfinite(max_lp)) return kInf;
      double sum_exp = 0.0;
      for (int q = 0; q < Q; ++q) {
        W_q[static_cast<size_t>(q)] = std::exp(log_post[static_cast<size_t>(q)] - max_lp);
        sum_exp += W_q[static_cast<size_t>(q)];
      }
      double log_marg = max_lp + std::log(sum_exp);
      total_nll -= log_marg;
      double inv_sum = 1.0 / sum_exp;
      double eap = 0.0;
      for (int q = 0; q < Q; ++q) {
        W_q[static_cast<size_t>(q)] *= inv_sum;
        eap += W_q[static_cast<size_t>(q)] * th_q[static_cast<size_t>(q)];
      }
      if (eap_out) (*eap_out)[static_cast<size_t>(p)] = eap;
      if (psd_out) {
        double pvar = 0.0;
        for (int q = 0; q < Q; ++q) {
          double dt = th_q[static_cast<size_t>(q)] - eap;
          pvar += W_q[static_cast<size_t>(q)] * dt * dt;
        }
        (*psd_out)[static_cast<size_t>(p)] = std::sqrt(std::max(0.0, pvar));
      }

      std::fill(p_grad_local.begin(), p_grad_local.end(), 0.0);

      if (n_xi > 0) {
        for (int loc = 0; loc < m_obs_p; ++loc) {
          int i = spec.p_obs_idx[static_cast<size_t>(beg + loc)];
          double w = spec.weight[static_cast<size_t>(i)];
          if (w <= 0.0) continue;
          int sg = spec.step_group[static_cast<size_t>(i)];
          int y = spec.score[static_cast<size_t>(i)];
          const auto& scols = smap.group_step_cols[static_cast<size_t>(sg)];
          const int n_sc = static_cast<int>(scols.size());
          double al = alpha_obs[static_cast<size_t>(i)];
          double eff_sc = (!has_slopes || sform == 1) ? 1.0 : al;

          double exp_mean = 0.0;
          const double* cm_row = cached_mean.data() + static_cast<size_t>(loc) * Q;
          for (int q = 0; q < Q; ++q) exp_mean += W_q[static_cast<size_t>(q)] * cm_row[q];

          if (n_sc > 0) {
            std::fill(step_ea_acc.begin(), step_ea_acc.begin() + n_sc, 0.0);
            for (int q = 0; q < Q; ++q) {
              double wq = W_q[static_cast<size_t>(q)];
              const double* ea_ptr = cached_step_ea.data() + (static_cast<size_t>(loc) * Q + q) * stride_sc;
              for (int m_idx = 0; m_idx < n_sc; ++m_idx) step_ea_acc[static_cast<size_t>(m_idx)] += wq * ea_ptr[m_idx];
            }
          }
          double g_eta = w * eff_sc * (exp_mean - static_cast<double>(y));
          for (int ptr = obs_f_ptr[static_cast<size_t>(i)]; ptr < obs_f_ptr[static_cast<size_t>(i) + 1]; ++ptr) {
            p_grad_local[static_cast<size_t>(obs_f_cols[static_cast<size_t>(ptr)])] +=
                g_eta * obs_f_vals[static_cast<size_t>(ptr)];
          }
          for (int m_idx = 0; m_idx < n_sc; ++m_idx) {
            int c_loc = scols[static_cast<size_t>(m_idx)];
            double a_obs = smap.cum_B[(static_cast<size_t>(sg) * C + y) * qs + c_loc];
            p_grad_local[static_cast<size_t>(smap.step_xi_start + c_loc)] +=
                w * eff_sc * (a_obs - step_ea_acc[static_cast<size_t>(m_idx)]);
          }
          if (has_slopes) {
            double g_ls = 0.0;
            const double* ls_row = cached_logslope_g.data() + static_cast<size_t>(loc) * Q;
            for (int q = 0; q < Q; ++q) g_ls += W_q[static_cast<size_t>(q)] * ls_row[q];
            g_ls *= w;
            for (int ptr = obs_s_ptr[static_cast<size_t>(i)]; ptr < obs_s_ptr[static_cast<size_t>(i) + 1]; ++ptr) {
              p_grad_local[static_cast<size_t>(obs_s_cols[static_cast<size_t>(ptr)])] +=
                  g_ls * obs_s_vals[static_cast<size_t>(ptr)];
            }
          }
        }
      }

      if (spec.pop_mode == 1) {
        if (!adaptive_active) {
          double g_mu = 0.0, g_logsig = 0.0;
          for (int q = 0; q < Q; ++q) {
            double wq = W_q[static_cast<size_t>(q)];
            g_mu -= wq * dlogL_dth[static_cast<size_t>(q)];
            g_logsig -= wq * dlogL_dth[static_cast<size_t>(q)] * (sigma * z_nodes[static_cast<size_t>(q)]);
          }
          for (int c = 0; c < dX; ++c) {
            p_grad_local[static_cast<size_t>(n_xi + c)] = g_mu * spec.pop_X[static_cast<size_t>(p) * dX + c];
          }
          p_grad_local[static_cast<size_t>(n_xi + dX)] = g_logsig;
        } else {
          double e_diff = 0.0, e_diff2 = 0.0;
          for (int q = 0; q < Q; ++q) {
            double diff = th_q[static_cast<size_t>(q)] - mp;
            e_diff += W_q[static_cast<size_t>(q)] * diff;
            e_diff2 += W_q[static_cast<size_t>(q)] * diff * diff;
          }
          double g_mu = -e_diff * inv_sig2;
          double g_logsig = 1.0 - e_diff2 * inv_sig2;
          for (int c = 0; c < dX; ++c) {
            p_grad_local[static_cast<size_t>(n_xi + c)] = g_mu * spec.pop_X[static_cast<size_t>(p) * dX + c];
          }
          p_grad_local[static_cast<size_t>(n_xi + dX)] = g_logsig;
        }
      }

      for (int j = 0; j < dim; ++j) {
        grad[static_cast<size_t>(j)] += p_grad_local[static_cast<size_t>(j)];
        if (person_scores_out) {
          (*person_scores_out)[static_cast<size_t>(p) * dim + j] = p_grad_local[static_cast<size_t>(j)];
        }
      }
    }

    if (ctrl.facet_l2_pen > 0.0 || ctrl.slope_l2_pen > 0.0) {
      total_nll += add_structural_l2_penalties(spec, smap, ctrl, beta, log_slope_est, grad);
    }
    return total_nll;
  };

  std::vector<double> phi0(static_cast<size_t>(dim), 0.0);
  if (static_cast<int>(ctrl.warm_init.size()) == dim) {
    phi0 = ctrl.warm_init;
  }
  auto fg_wrap = [&](const std::vector<double>& x, std::vector<double>& g) {
    return eval_mml_core(x, g, nullptr, nullptr, nullptr);
  };
  LbfgsResult opt = minimize_lbfgs(phi0, fg_wrap, std::min(ctrl.maxit, 250), 1e-8, 1e-13, 14);
  std::vector<double> phi = opt.x;
  int total_iters = opt.iterations;

  if (ctrl.adaptive_quad) {
    adaptive_active = true;
    for (int cyc = 0; cyc < 5; ++cyc) {
      std::vector<double> prev_phi = phi;
      refresh_adaptive_basis(phi);
      LbfgsResult sub = minimize_lbfgs(phi, fg_wrap, 80, 1e-8, 1e-13, 12);
      phi = sub.x;
      total_iters += sub.iterations;
      double mdiff = 0.0;
      for (int j = 0; j < dim; ++j) mdiff = std::max(mdiff, std::abs(phi[static_cast<size_t>(j)] - prev_phi[static_cast<size_t>(j)]));
      if (mdiff < 1e-7) break;
    }
    refresh_adaptive_basis(phi);
  }

  auto compute_hessian_and_scores = [&](const std::vector<double>& x_at,
                                        std::vector<double>& g_at,
                                        std::vector<double>& H_out,
                                        std::vector<double>& S_person_out,
                                        double& nll_at) {
    if (ctrl.adaptive_quad) refresh_adaptive_basis(x_at);
    nll_at = eval_mml_core(x_at, g_at, &S_person_out, &out.person_eap, &out.person_psd);
    H_out.assign(static_cast<size_t>(dim) * dim, 0.0);
    if (dim == 0) return;

    const double eps = 1e-5;
    std::vector<double> xp = x_at, xm = x_at, gp, gm;
    if (!ctrl.adaptive_quad) {
      for (int j = 0; j < dim; ++j) {
        double h = eps * std::max(1.0, std::abs(x_at[static_cast<size_t>(j)]));
        xp[static_cast<size_t>(j)] = x_at[static_cast<size_t>(j)] + h;
        xm[static_cast<size_t>(j)] = x_at[static_cast<size_t>(j)] - h;
        eval_mml_core(xp, gp, nullptr, nullptr, nullptr);
        eval_mml_core(xm, gm, nullptr, nullptr, nullptr);
        for (int k = 0; k < dim; ++k) {
          H_out[static_cast<size_t>(j) * dim + k] = (gp[static_cast<size_t>(k)] - gm[static_cast<size_t>(k)]) / (2.0 * h);
        }
        xp[static_cast<size_t>(j)] = x_at[static_cast<size_t>(j)];
        xm[static_cast<size_t>(j)] = x_at[static_cast<size_t>(j)];
      }
    } else {
      S_person_out.assign(static_cast<size_t>(P) * dim, 0.0);
      auto eval_adap_nll_per_person = [&](const std::vector<double>& xv, std::vector<double>& p_nll) -> double {
        refresh_adaptive_basis(xv);
        double sigma = 1.0;
        unpack_pop(xv, sigma);
        if (!std::isfinite(sigma) || sigma <= 1e-8) return kInf;
        smap.expand(spec, xv.data(), beta, tau_mat, step_cum_mat);
        if (!smap.expand_slopes(spec, xv.data(), log_slope_est, slope_est, alpha_obs)) return kInf;
        for (int i = 0; i < spec.n_obs; ++i) {
          double e = 0.0;
          for (int f = 0; f < F; ++f) {
            int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
            e += spec.facet_signs[static_cast<size_t>(f)] * beta[static_cast<size_t>(f)][static_cast<size_t>(l)];
          }
          base_eta_obs[static_cast<size_t>(i)] = e;
        }
        p_nll.assign(static_cast<size_t>(P), 0.0);
        std::vector<double> probs(static_cast<size_t>(C)), lp(static_cast<size_t>(Q));
        double inv_s2 = 1.0 / (sigma * sigma), tot = 0.0;
        for (int p = 0; p < P; ++p) {
          int beg = spec.p_row_ptr[static_cast<size_t>(p)], end = spec.p_row_ptr[static_cast<size_t>(p) + 1];
          if (beg == end) continue;
          double mp = mu_p[static_cast<size_t>(p)], sp = adap_scale[static_cast<size_t>(p)], md = adap_mode[static_cast<size_t>(p)];
          for (int q = 0; q < Q; ++q) {
            double z = z_nodes[static_cast<size_t>(q)];
            double tq = md + sp * z, df = tq - mp;
            lp[static_cast<size_t>(q)] = log_w_nodes[static_cast<size_t>(q)] + 0.5 * z * z + std::log(sp / sigma) - 0.5 * df * df * inv_s2;
            for (int idx = beg; idx < end; ++idx) {
              int i = spec.p_obs_idx[static_cast<size_t>(idx)];
              double w = spec.weight[static_cast<size_t>(i)];
              if (w <= 0.0) continue;
              int sg = spec.step_group[static_cast<size_t>(i)];
              lp[static_cast<size_t>(q)] += w * eval_gpcm_probs_and_moments(
                  tq, base_eta_obs[static_cast<size_t>(i)], alpha_obs[static_cast<size_t>(i)], sform,
                  step_cum_mat.data() + static_cast<size_t>(sg) * C, C, spec.score[static_cast<size_t>(i)], probs.data(),
                  smap.cat_min[static_cast<size_t>(sg)], smap.cat_max[static_cast<size_t>(sg)]).log_py;
            }
          }
          double mx = *std::max_element(lp.begin(), lp.end()), se = 0.0;
          for (int q = 0; q < Q; ++q) se += std::exp(lp[static_cast<size_t>(q)] - mx);
          p_nll[static_cast<size_t>(p)] = -(mx + std::log(se));
          tot += p_nll[static_cast<size_t>(p)];
        }
        if (ctrl.facet_l2_pen > 0.0 || ctrl.slope_l2_pen > 0.0) {
          std::vector<double> dummy_g;
          tot += add_structural_l2_penalties(spec, smap, ctrl, beta, log_slope_est, dummy_g);
        }
        return tot;
      };
      std::vector<double> pnll0, pnll_p, pnll_m;
      nll_at = eval_adap_nll_per_person(x_at, pnll0);
      g_at.assign(static_cast<size_t>(dim), 0.0);
      const double h_fd = 1e-4;
      for (int j = 0; j < dim; ++j) {
        double h = h_fd * std::max(1.0, std::abs(x_at[static_cast<size_t>(j)]));
        xp[static_cast<size_t>(j)] = x_at[static_cast<size_t>(j)] + h;
        xm[static_cast<size_t>(j)] = x_at[static_cast<size_t>(j)] - h;
        double fp = eval_adap_nll_per_person(xp, pnll_p);
        double fm = eval_adap_nll_per_person(xm, pnll_m);
        g_at[static_cast<size_t>(j)] = (fp - fm) / (2.0 * h);
        H_out[static_cast<size_t>(j) * dim + j] = (fp - 2.0 * nll_at + fm) / (h * h);
        for (int p = 0; p < P; ++p) {
          S_person_out[static_cast<size_t>(p) * dim + j] = (pnll_p[static_cast<size_t>(p)] - pnll_m[static_cast<size_t>(p)]) / (2.0 * h);
        }
        xp[static_cast<size_t>(j)] = x_at[static_cast<size_t>(j)];
        xm[static_cast<size_t>(j)] = x_at[static_cast<size_t>(j)];
      }
      for (int j = 0; j < dim; ++j) {
        double hj = h_fd * std::max(1.0, std::abs(x_at[static_cast<size_t>(j)]));
        for (int k = j + 1; k < dim; ++k) {
          double hk = h_fd * std::max(1.0, std::abs(x_at[static_cast<size_t>(k)]));
          std::vector<double> xpp = x_at, xpm = x_at, xmp = x_at, xmm = x_at, tmp;
          xpp[static_cast<size_t>(j)] += hj; xpp[static_cast<size_t>(k)] += hk;
          xpm[static_cast<size_t>(j)] += hj; xpm[static_cast<size_t>(k)] -= hk;
          xmp[static_cast<size_t>(j)] -= hj; xmp[static_cast<size_t>(k)] += hk;
          xmm[static_cast<size_t>(j)] -= hj; xmm[static_cast<size_t>(k)] -= hk;
          double v = (eval_adap_nll_per_person(xpp, tmp) - eval_adap_nll_per_person(xpm, tmp) -
                      eval_adap_nll_per_person(xmp, tmp) + eval_adap_nll_per_person(xmm, tmp)) / (4.0 * hj * hk);
          H_out[static_cast<size_t>(j) * dim + k] = v;
          H_out[static_cast<size_t>(k) * dim + j] = v;
        }
      }
      refresh_adaptive_basis(x_at);
      std::vector<double> dummy_g;
      eval_mml_core(x_at, dummy_g, nullptr, &out.person_eap, &out.person_psd);
    }
    for (int j = 0; j < dim; ++j) {
      for (int k = j + 1; k < dim; ++k) {
        double s = 0.5 * (H_out[static_cast<size_t>(j) * dim + k] + H_out[static_cast<size_t>(k) * dim + j]);
        H_out[static_cast<size_t>(j) * dim + k] = s;
        H_out[static_cast<size_t>(k) * dim + j] = s;
      }
    }
  };

  std::vector<double> g_fin, H_full, S_person;
  double nll_fin = kNaN;
  compute_hessian_and_scores(phi, g_fin, H_full, S_person, nll_fin);
  auto max_abs_v = [](const std::vector<double>& v) {
    double m = 0.0;
    for (double val : v) m = std::max(m, std::abs(val));
    return m;
  };
  if (!ctrl.adaptive_quad) {
    for (int pol = 0; pol < 3 && max_abs_v(g_fin) > 1e-9; ++pol) {
      std::vector<double> d_phi;
      if (!damped_spd_solve(H_full, dim, g_fin, d_phi, 0.0)) break;
      std::vector<double> cand(static_cast<size_t>(dim)), g_c;
      double scale = 1.0;
      bool ok = false;
      for (int ls = 0; ls < 10; ++ls) {
        for (int j = 0; j < dim; ++j) cand[static_cast<size_t>(j)] = phi[static_cast<size_t>(j)] - scale * d_phi[static_cast<size_t>(j)];
        double fc = eval_mml_core(cand, g_c, nullptr, nullptr, nullptr);
        if (std::isfinite(fc) && fc <= nll_fin + 1e-13 * (1.0 + std::abs(nll_fin))) {
          phi = cand;
          ok = true;
          break;
        }
        scale *= 0.5;
      }
      if (!ok) break;
      compute_hessian_and_scores(phi, g_fin, H_full, S_person, nll_fin);
    }
  }

  out.iterations = total_iters;
  out.loglik = -nll_fin;
  out.objective = nll_fin;
  out.max_abs_grad = max_abs_v(g_fin);
  out.converged = std::isfinite(nll_fin) && (out.max_abs_grad <= 1e-4);
  out.xi_est = phi; // Full parameter vector [n_xi + dX + has_free_sigma] for warm starts & vcov!

  double sigma_hat = 1.0;
  unpack_pop(phi, sigma_hat);
  out.pop_sigma = sigma_hat;
  out.pop_beta.assign(static_cast<size_t>(dX), 0.0);
  for (int c = 0; c < dX; ++c) out.pop_beta[static_cast<size_t>(c)] = phi[static_cast<size_t>(n_xi + c)];
  smap.expand(spec, phi.data(), out.facet_est, out.step_est, out.step_cum);
  smap.expand_slopes(spec, phi.data(), out.log_slope_est, out.slope_est, out.alpha_obs);

  // Boundary guard: if free sigma collapses near 0, withhold Wald covariance
  std::vector<double> V_model, V_sw;
  double logdet_H = 0.0;
  bool sigma_interior = (spec.pop_mode == 0) || (std::isfinite(sigma_hat) && sigma_hat > 1e-3);
  bool cov_ok = sigma_interior && (dim > 0) && chol_invert_and_logdet(H_full, dim, V_model, logdet_H);
  if (cov_ok) {
    std::vector<double> J_meat(static_cast<size_t>(dim) * dim, 0.0);
    for (int p = 0; p < P; ++p) {
      const double* sp = S_person.data() + static_cast<size_t>(p) * dim;
      for (int a = 0; a < dim; ++a) {
        for (int b = 0; b <= a; ++b) {
          J_meat[static_cast<size_t>(a) * dim + b] += sp[a] * sp[b];
        }
      }
    }
    for (int a = 0; a < dim; ++a) {
      for (int b = 0; b < a; ++b) J_meat[static_cast<size_t>(b) * dim + a] = J_meat[static_cast<size_t>(a) * dim + b];
    }
    std::vector<double> VJ(static_cast<size_t>(dim) * dim, 0.0);
    V_sw.assign(static_cast<size_t>(dim) * dim, 0.0);
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < dim; ++j) {
        double s = 0.0;
        for (int k = 0; k < dim; ++k) s += V_model[static_cast<size_t>(i) * dim + k] * J_meat[static_cast<size_t>(k) * dim + j];
        VJ[static_cast<size_t>(i) * dim + j] = s;
      }
    }
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j <= i; ++j) {
        double s = 0.0;
        for (int k = 0; k < dim; ++k) s += VJ[static_cast<size_t>(i) * dim + k] * V_model[static_cast<size_t>(k) * dim + j];
        V_sw[static_cast<size_t>(i) * dim + j] = s;
        V_sw[static_cast<size_t>(j) * dim + i] = s;
      }
    }
  }
  out.xi_cov_model = V_model;
  out.xi_cov_sandwich = V_sw;

  out.pop_se_model.assign(static_cast<size_t>(dX + has_free_sigma), kNaN);
  out.pop_se_sandwich.assign(static_cast<size_t>(dX + has_free_sigma), kNaN);
  std::vector<double> Vxi_mod(static_cast<size_t>(n_xi) * n_xi, kNaN);
  std::vector<double> Vxi_sw(static_cast<size_t>(n_xi) * n_xi, kNaN);
  if (cov_ok && static_cast<int>(V_model.size()) == dim * dim &&
      static_cast<int>(V_sw.size()) == dim * dim) {
    for (int a = 0; a < n_xi; ++a) {
      for (int b = 0; b < n_xi; ++b) {
        Vxi_mod[static_cast<size_t>(a) * n_xi + b] = V_model[static_cast<size_t>(a) * dim + b];
        Vxi_sw[static_cast<size_t>(a) * n_xi + b] = V_sw[static_cast<size_t>(a) * dim + b];
      }
    }
    for (int c = 0; c < dX + has_free_sigma; ++c) {
      double vm = V_model[static_cast<size_t>(n_xi + c) * dim + (n_xi + c)];
      double vs = V_sw[static_cast<size_t>(n_xi + c) * dim + (n_xi + c)];
      out.pop_se_model[static_cast<size_t>(c)] = (vm >= 0.0) ? std::sqrt(vm) : kNaN;
      out.pop_se_sandwich[static_cast<size_t>(c)] = (vs >= 0.0) ? std::sqrt(vs) : kNaN;
    }
  }
  std::vector<double> step_sw_se;
  smap.propagate_covariance(spec, Vxi_mod, out.facet_se, out.facet_cov_model, out.step_se, &out.step_cov_model);
  smap.propagate_covariance(spec, Vxi_sw, out.facet_joint_se, out.facet_cov_sandwich, step_sw_se, &out.step_cov_sandwich);
  if (has_slopes) {
    smap.propagate_slope_covariance(spec, Vxi_mod, out.slope_est, out.slope_se, out.slope_cov_model);
    smap.propagate_slope_covariance(spec, Vxi_sw, out.slope_est, out.slope_joint_se, out.slope_cov_sandwich);
  }

  for (int i = 0; i < spec.n_obs; ++i) {
    double e = 0.0;
    for (int f = 0; f < F; ++f) {
      int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
      e += spec.facet_signs[static_cast<size_t>(f)] * out.facet_est[static_cast<size_t>(f)][static_cast<size_t>(l)];
    }
    base_eta_obs[static_cast<size_t>(i)] = e;
  }
  const std::vector<double>* alpha_ptr = has_slopes ? &out.alpha_obs : nullptr;
  if (ctrl.compute_wle) {
    WleResult wle = score_wle_persons(spec, base_eta_obs, out.step_cum, alpha_ptr);
    out.person_wle = std::move(wle.theta);
    out.person_wle_se = std::move(wle.se);
  }

  {
    std::vector<double> eaps, psds2;
    for (int p = 0; p < P; ++p) {
      if (out.connectivity.person_obs_count[static_cast<size_t>(p)] > 0.0 &&
          std::isfinite(out.person_eap[static_cast<size_t>(p)]) &&
          std::isfinite(out.person_psd[static_cast<size_t>(p)])) {
        eaps.push_back(out.person_eap[static_cast<size_t>(p)]);
        psds2.push_back(out.person_psd[static_cast<size_t>(p)] * out.person_psd[static_cast<size_t>(p)]);
      }
    }
    if (eaps.size() >= 2) {
      double m_eap = std::accumulate(eaps.begin(), eaps.end(), 0.0) / eaps.size();
      double v_eap = 0.0;
      for (double v : eaps) v_eap += (v - m_eap) * (v - m_eap);
      v_eap /= eaps.size();
      double s_err = std::accumulate(psds2.begin(), psds2.end(), 0.0) / psds2.size();
      out.person_rel_eap = (v_eap + s_err > 0.0) ? (v_eap / (v_eap + s_err)) : 0.0;
    }
  }

  out.person_est = ctrl.wle_primary ? out.person_wle : out.person_eap;
  out.person_se = ctrl.wle_primary ? out.person_wle_se : out.person_psd;
  std::vector<bool> diag_exclude(static_cast<size_t>(P), false);
  for (int p = 0; p < P; ++p) {
    if (out.connectivity.person_obs_count[static_cast<size_t>(p)] <= 0.0) diag_exclude[static_cast<size_t>(p)] = true;
  }

  compute_fast_diagnostics_and_fair(
      spec, out.person_est, diag_exclude, out.facet_est, out.step_cum,
      ctrl.whexact, out.person_diag, out.facet_diag, alpha_ptr);

  out.person_real_se.assign(static_cast<size_t>(P), kNaN);
  out.person_cond_se.assign(static_cast<size_t>(P), kNaN);
  for (int p = 0; p < P; ++p) {
    double se = out.person_se[static_cast<size_t>(p)];
    double ci = out.person_diag.cond_info[static_cast<size_t>(p)];
    if (ci > 1e-12) out.person_cond_se[static_cast<size_t>(p)] = 1.0 / std::sqrt(ci);
    if (std::isfinite(se)) {
      double infit = out.person_diag.infit_mnsq[static_cast<size_t>(p)];
      out.person_real_se[static_cast<size_t>(p)] = std::isfinite(infit) ? (se * std::sqrt(std::max(1.0, infit))) : se;
    }
  }
  out.facet_real_se.resize(static_cast<size_t>(F));
  for (int f = 0; f < F; ++f) {
    int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
    out.facet_real_se[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), kNaN);
    for (int l = 0; l < Lf; ++l) {
      double se = out.facet_se[static_cast<size_t>(f)][static_cast<size_t>(l)];
      if (std::isfinite(se)) {
        double infit = out.facet_diag[static_cast<size_t>(f)].infit_mnsq[static_cast<size_t>(l)];
        out.facet_real_se[static_cast<size_t>(f)][static_cast<size_t>(l)] =
            std::isfinite(infit) ? (se * std::sqrt(std::max(1.0, infit))) : se;
      }
    }
  }
  if (ctrl.bias_facet_a >= 0 && ctrl.bias_facet_b >= 0) {
    out.local_bias_cells = compute_local_bias_2way(
        spec, out.person_est, diag_exclude, out.facet_est, out.step_cum,
        ctrl.bias_facet_a, ctrl.bias_facet_b, alpha_ptr);
  }
}

// ============================================================================
// Solver 4: Ultra-Fast C++ Strict CMLE + WLE (Affine Facets/Steps & Support)
// ============================================================================

inline void solve_cml(const FastDesignSpec& spec, const StructuralMap& smap,
                      const FastControl& ctrl, FastFitOutput& out) {
  const int P = spec.n_person;
  const int F = spec.n_facets;
  const int C = spec.n_cat;
  const int n_xi = smap.n_xi;
  const int qs = smap.step_xi_len;
  out.n_free_par = n_xi;
  const auto lderivs = smap.build_level_derivs(spec);

  struct CellFeature {
    std::vector<int> facet_levels;
    int step_group = 0;
    int a_g = 0;
    int K_g = 1;
    std::vector<double> const_offset; // [K_g + 1]
    std::vector<double> X;            // [(K_g + 1) * n_xi]
  };
  std::map<std::vector<int>, int> cell_map;
  std::vector<CellFeature> cells;
  std::vector<int> obs_cell(static_cast<size_t>(spec.n_obs), 0);

  std::vector<int> tmp_c(std::max(512, n_xi + 16));
  std::vector<double> tmp_v(std::max(512, n_xi + 16));
  for (int i = 0; i < spec.n_obs; ++i) {
    std::vector<int> key(static_cast<size_t>(F) + 1);
    for (int f = 0; f < F; ++f) key[static_cast<size_t>(f)] = spec.facet_flat[static_cast<size_t>(i) * F + f];
    key[static_cast<size_t>(F)] = spec.step_group[static_cast<size_t>(i)];
    auto it = cell_map.find(key);
    if (it != cell_map.end()) {
      obs_cell[static_cast<size_t>(i)] = it->second;
    } else {
      int cid = static_cast<int>(cells.size());
      cell_map[key] = cid;
      obs_cell[static_cast<size_t>(i)] = cid;
      CellFeature cf;
      cf.facet_levels.assign(key.begin(), key.begin() + F);
      cf.step_group = key[static_cast<size_t>(F)];
      cf.a_g = smap.cat_min[static_cast<size_t>(cf.step_group)];
      cf.K_g = smap.cat_max[static_cast<size_t>(cf.step_group)] - cf.a_g;
      int Cg = cf.K_g + 1;
      cf.const_offset.assign(static_cast<size_t>(Cg), 0.0);
      cf.X.assign(static_cast<size_t>(Cg) * std::max(1, n_xi), 0.0);

      double eta0 = 0.0;
      for (int f = 0; f < F; ++f) {
        int l = cf.facet_levels[static_cast<size_t>(f)];
        eta0 += spec.facet_signs[static_cast<size_t>(f)] * smap.facet_d0[static_cast<size_t>(f)][static_cast<size_t>(l)];
      }
      for (int ks = 0; ks <= cf.K_g; ++ks) {
        int k = cf.a_g + ks;
        double cd0 = smap.cum_d0[static_cast<size_t>(cf.step_group) * C + k];
        cf.const_offset[static_cast<size_t>(ks)] = std::isfinite(cd0) ? (static_cast<double>(ks) * eta0 - cd0) : -kInf;
      }
      if (n_xi > 0) {
        int nnz = obs_eta_xi_terms(i, spec, lderivs, tmp_c.data(), tmp_v.data());
        for (int ks = 1; ks <= cf.K_g; ++ks) {
          for (int h = 0; h < nnz; ++h) {
            cf.X[static_cast<size_t>(ks) * n_xi + tmp_c[static_cast<size_t>(h)]] +=
                static_cast<double>(ks) * tmp_v[static_cast<size_t>(h)];
          }
        }
        const auto& scols = smap.group_step_cols[static_cast<size_t>(cf.step_group)];
        for (int ks = 0; ks <= cf.K_g; ++ks) {
          int k = cf.a_g + ks;
          for (int c_loc : scols) {
            cf.X[static_cast<size_t>(ks) * n_xi + (smap.step_xi_start + c_loc)] -=
                smap.cum_B[(static_cast<size_t>(cf.step_group) * C + k) * qs + c_loc];
          }
        }
      }
      cells.push_back(std::move(cf));
    }
  }

  struct CmlRosterGroup {
    std::vector<int> roster;
    std::vector<double> s_counts;
  };
  std::map<std::vector<int>, int> roster_map;
  std::vector<CmlRosterGroup> groups;
  std::vector<double> t_obs(static_cast<size_t>(n_xi), 0.0);
  double const_obs_offset = 0.0;

  for (int p = 0; p < P; ++p) {
    if (out.connectivity.person_obs_count[static_cast<size_t>(p)] <= 0.0 ||
        out.connectivity.person_extreme_flag[static_cast<size_t>(p)] != 0) {
      continue;
    }
    int beg = spec.p_row_ptr[static_cast<size_t>(p)];
    int end = spec.p_row_ptr[static_cast<size_t>(p) + 1];
    std::vector<int> roster;
    roster.reserve(static_cast<size_t>(end - beg));
    int total_s = 0;
    int max_s = 0;
    for (int idx = beg; idx < end; ++idx) {
      int i = spec.p_obs_idx[static_cast<size_t>(idx)];
      int cid = obs_cell[static_cast<size_t>(i)];
      const auto& cf = cells[static_cast<size_t>(cid)];
      int ys = spec.score[static_cast<size_t>(i)] - cf.a_g;
      roster.push_back(cid);
      total_s += ys;
      max_s += cf.K_g;
      const_obs_offset += cf.const_offset[static_cast<size_t>(ys)];
      for (int j = 0; j < n_xi; ++j) {
        t_obs[static_cast<size_t>(j)] += cf.X[static_cast<size_t>(ys) * n_xi + j];
      }
    }
    std::sort(roster.begin(), roster.end());
    auto it = roster_map.find(roster);
    int gid = 0;
    if (it == roster_map.end()) {
      gid = static_cast<int>(groups.size());
      roster_map[roster] = gid;
      CmlRosterGroup grp;
      grp.roster = roster;
      grp.s_counts.assign(static_cast<size_t>(max_s + 1), 0.0);
      groups.push_back(std::move(grp));
    } else {
      gid = it->second;
    }
    groups[static_cast<size_t>(gid)].s_counts[static_cast<size_t>(total_s)] += 1.0;
  }

  auto eval_cml = [&](const std::vector<double>& xi_v,
                      std::vector<double>& grad,
                      std::vector<double>* H_out) -> double {
    grad.assign(static_cast<size_t>(n_xi), 0.0);
    for (int j = 0; j < n_xi; ++j) grad[static_cast<size_t>(j)] = -t_obs[static_cast<size_t>(j)];
    if (H_out) H_out->assign(static_cast<size_t>(n_xi) * n_xi, 0.0);
    double nll = -const_obs_offset;
    for (int j = 0; j < n_xi; ++j) nll -= t_obs[static_cast<size_t>(j)] * xi_v[static_cast<size_t>(j)];

    const int n_cells = static_cast<int>(cells.size());
    std::vector<std::vector<double>> cell_eta(static_cast<size_t>(n_cells));
    for (int cid = 0; cid < n_cells; ++cid) {
      const auto& cf = cells[static_cast<size_t>(cid)];
      int Cg = cf.K_g + 1;
      cell_eta[static_cast<size_t>(cid)].assign(static_cast<size_t>(Cg), -kInf);
      for (int ks = 0; ks <= cf.K_g; ++ks) {
        double c0 = cf.const_offset[static_cast<size_t>(ks)];
        if (!std::isfinite(c0)) continue;
        double v = c0;
        for (int j = 0; j < n_xi; ++j) v += cf.X[static_cast<size_t>(ks) * n_xi + j] * xi_v[static_cast<size_t>(j)];
        if (!std::isfinite(v) || std::abs(v) > 250.0) return kInf;
        cell_eta[static_cast<size_t>(cid)][static_cast<size_t>(ks)] = v;
      }
    }

    std::vector<double> logz, next_logz, means, next_means, covs, next_covs;
    std::vector<double> logs_k(static_cast<size_t>(C)), w_k(static_cast<size_t>(C));
    std::vector<double> cand_mu(static_cast<size_t>(C) * std::max(1, n_xi));

    for (const auto& grp : groups) {
      logz.assign(1, 0.0);
      means.assign(static_cast<size_t>(std::max(1, n_xi)), 0.0);
      if (H_out) covs.assign(static_cast<size_t>(std::max(1, n_xi)) * std::max(1, n_xi), 0.0);

      for (int cid : grp.roster) {
        const auto& cf = cells[static_cast<size_t>(cid)];
        const auto& eta_c = cell_eta[static_cast<size_t>(cid)];
        int prev_nz = static_cast<int>(logz.size());
        int nz = prev_nz + cf.K_g;
        next_logz.assign(static_cast<size_t>(nz), -kInf);
        next_means.assign(static_cast<size_t>(nz) * std::max(1, n_xi), 0.0);
        if (H_out) next_covs.assign(static_cast<size_t>(nz) * std::max(1, n_xi) * std::max(1, n_xi), 0.0);

        for (int s = 0; s < nz; ++s) {
          int k_min = std::max(0, s - prev_nz + 1);
          int k_max = std::min(cf.K_g, s);
          int len_k = k_max - k_min + 1;
          double peak = -kInf;
          for (int h = 0; h < len_k; ++h) {
            int ks = k_min + h;
            int from = s - ks;
            double lv = (std::isfinite(logz[static_cast<size_t>(from)]) && std::isfinite(eta_c[static_cast<size_t>(ks)]))
                            ? (logz[static_cast<size_t>(from)] + eta_c[static_cast<size_t>(ks)])
                            : -kInf;
            logs_k[static_cast<size_t>(h)] = lv;
            if (lv > peak) peak = lv;
          }
          if (!std::isfinite(peak)) continue;
          double sw = 0.0;
          for (int h = 0; h < len_k; ++h) {
            w_k[static_cast<size_t>(h)] = std::isfinite(logs_k[static_cast<size_t>(h)])
                                              ? std::exp(logs_k[static_cast<size_t>(h)] - peak)
                                              : 0.0;
            sw += w_k[static_cast<size_t>(h)];
          }
          next_logz[static_cast<size_t>(s)] = peak + std::log(sw);
          double inv_sw = 1.0 / sw;
          for (int h = 0; h < len_k; ++h) w_k[static_cast<size_t>(h)] *= inv_sw;

          if (n_xi > 0) {
            double* mu_s = next_means.data() + static_cast<size_t>(s) * n_xi;
            for (int h = 0; h < len_k; ++h) {
              int ks = k_min + h;
              int from = s - ks;
              const double* m_from = means.data() + static_cast<size_t>(from) * n_xi;
              const double* x_row = cf.X.data() + static_cast<size_t>(ks) * n_xi;
              double* cm_row = cand_mu.data() + static_cast<size_t>(h) * n_xi;
              double wh = w_k[static_cast<size_t>(h)];
              for (int j = 0; j < n_xi; ++j) {
                double cval = m_from[j] + x_row[j];
                cm_row[j] = cval;
                mu_s[j] += wh * cval;
              }
            }
            if (H_out) {
              double* cov_s = next_covs.data() + static_cast<size_t>(s) * n_xi * n_xi;
              for (int h = 0; h < len_k; ++h) {
                double wh = w_k[static_cast<size_t>(h)];
                if (wh <= 0.0) continue;
                int ks = k_min + h;
                int from = s - ks;
                const double* cm_row = cand_mu.data() + static_cast<size_t>(h) * n_xi;
                const double* cov_from = covs.data() + static_cast<size_t>(from) * n_xi * n_xi;
                for (int a = 0; a < n_xi; ++a) {
                  double da = cm_row[a] - mu_s[a];
                  for (int b = 0; b <= a; ++b) {
                    cov_s[a * n_xi + b] += wh * (da * (cm_row[b] - mu_s[b]) + cov_from[a * n_xi + b]);
                  }
                }
              }
            }
          }
        }
        logz.swap(next_logz);
        means.swap(next_means);
        if (H_out) covs.swap(next_covs);
      }

      for (size_t s = 0; s < grp.s_counts.size(); ++s) {
        double cnt = grp.s_counts[s];
        if (cnt <= 0.0) continue;
        if (!std::isfinite(logz[s])) return kInf;
        nll += cnt * logz[s];
        for (int j = 0; j < n_xi; ++j) grad[static_cast<size_t>(j)] += cnt * means[s * n_xi + j];
        if (H_out) {
          const double* cov_s = covs.data() + s * n_xi * n_xi;
          for (int a = 0; a < n_xi; ++a) {
            for (int b = 0; b <= a; ++b) {
              (*H_out)[static_cast<size_t>(a) * n_xi + b] += cnt * cov_s[a * n_xi + b];
            }
          }
        }
      }
    }
    if (H_out) {
      for (int a = 0; a < n_xi; ++a) {
        for (int b = 0; b < a; ++b) (*H_out)[static_cast<size_t>(b) * n_xi + a] = (*H_out)[static_cast<size_t>(a) * n_xi + b];
      }
    }
    return nll;
  };

  std::vector<double> xi(static_cast<size_t>(n_xi), 0.0), g_xi, H_xi;
  if (static_cast<int>(ctrl.warm_init.size()) == n_xi) xi = ctrl.warm_init;
  double cur_nll = eval_cml(xi, g_xi, &H_xi);
  auto max_abs_v = [](const std::vector<double>& v) {
    double m = 0.0;
    for (double val : v) m = std::max(m, std::abs(val));
    return m;
  };
  int iter_used = 0;
  for (int iter = 1; iter <= std::min(ctrl.maxit, 60); ++iter) {
    iter_used = iter;
    if (max_abs_v(g_xi) <= std::max(1e-10, ctrl.reltol)) {
      out.converged = true;
      break;
    }
    std::vector<double> d_xi;
    if (!damped_spd_solve(H_xi, n_xi, g_xi, d_xi, 0.0)) break;
    double mstep = max_abs_v(d_xi);
    double scale = (mstep > 3.0) ? (3.0 / mstep) : 1.0;
    std::vector<double> cand(static_cast<size_t>(n_xi)), g_cand;
    bool ok = false;
    for (int ls = 0; ls < 20; ++ls) {
      for (int j = 0; j < n_xi; ++j) cand[static_cast<size_t>(j)] = xi[static_cast<size_t>(j)] - scale * d_xi[static_cast<size_t>(j)];
      double fc = eval_cml(cand, g_cand, nullptr);
      if (std::isfinite(fc) && fc <= cur_nll + 1e-12 * (1.0 + std::abs(cur_nll))) {
        xi = cand;
        ok = true;
        break;
      }
      scale *= 0.5;
    }
    if (!ok) break;
    cur_nll = eval_cml(xi, g_xi, &H_xi);
    if (max_abs_v(xi) > 30.0) {
      out.converged = false;
      break;
    }
  }
  if (max_abs_v(g_xi) <= 1e-6 && max_abs_v(xi) <= 30.0) out.converged = true;

  out.iterations = iter_used;
  out.loglik = -cur_nll;
  out.objective = cur_nll;
  out.max_abs_grad = max_abs_v(g_xi);
  out.xi_est = xi;

  smap.expand(spec, xi.data(), out.facet_est, out.step_est, out.step_cum);

  double logdet_H = 0.0;
  chol_invert_and_logdet(H_xi, n_xi, out.xi_cov_model, logdet_H);
  smap.propagate_covariance(spec, out.xi_cov_model, out.facet_se, out.facet_cov_model, out.step_se, &out.step_cov_model);
  out.facet_joint_se = out.facet_se;

  std::vector<double> base_eta_obs(static_cast<size_t>(spec.n_obs), 0.0);
  for (int i = 0; i < spec.n_obs; ++i) {
    double e = 0.0;
    for (int f = 0; f < F; ++f) {
      int l = spec.facet_flat[static_cast<size_t>(i) * F + f];
      e += spec.facet_signs[static_cast<size_t>(f)] * out.facet_est[static_cast<size_t>(f)][static_cast<size_t>(l)];
    }
    base_eta_obs[static_cast<size_t>(i)] = e;
  }
  WleResult wle = score_wle_persons(spec, base_eta_obs, out.step_cum);
  out.person_wle = wle.theta;
  out.person_wle_se = wle.se;
  out.person_est = wle.theta;
  out.person_se = wle.se;
  out.person_cond_se = wle.se;

  std::vector<bool> diag_exclude(static_cast<size_t>(P), false);
  for (int p = 0; p < P; ++p) {
    if (out.connectivity.person_obs_count[static_cast<size_t>(p)] <= 0.0) diag_exclude[static_cast<size_t>(p)] = true;
  }
  compute_fast_diagnostics_and_fair(
      spec, out.person_est, diag_exclude, out.facet_est, out.step_cum,
      ctrl.whexact, out.person_diag, out.facet_diag);

  out.person_real_se.assign(static_cast<size_t>(P), kNaN);
  for (int p = 0; p < P; ++p) {
    double se = out.person_se[static_cast<size_t>(p)];
    if (std::isfinite(se)) {
      double infit = out.person_diag.infit_mnsq[static_cast<size_t>(p)];
      out.person_real_se[static_cast<size_t>(p)] = std::isfinite(infit) ? (se * std::sqrt(std::max(1.0, infit))) : se;
    }
  }
  out.facet_real_se.resize(static_cast<size_t>(F));
  for (int f = 0; f < F; ++f) {
    int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
    out.facet_real_se[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf), kNaN);
    for (int l = 0; l < Lf; ++l) {
      double se = out.facet_se[static_cast<size_t>(f)][static_cast<size_t>(l)];
      if (std::isfinite(se)) {
        double infit = out.facet_diag[static_cast<size_t>(f)].infit_mnsq[static_cast<size_t>(l)];
        out.facet_real_se[static_cast<size_t>(f)][static_cast<size_t>(l)] =
            std::isfinite(infit) ? (se * std::sqrt(std::max(1.0, infit))) : se;
      }
    }
  }
  if (ctrl.bias_facet_a >= 0 && ctrl.bias_facet_b >= 0) {
    out.local_bias_cells = compute_local_bias_2way(
        spec, out.person_est, diag_exclude, out.facet_est, out.step_cum,
        ctrl.bias_facet_a, ctrl.bias_facet_b);
  }
}

// ============================================================================
// Unified Single-Dataset & Multi-Threaded Batch Entry Points
// ============================================================================

inline FastFitOutput run_fast_fit(FastDesignSpec& spec, const FastControl& ctrl) {
  if (spec.p_row_ptr.empty()) spec.build_person_csr();
  FastFitOutput out;
  out.connectivity = check_connectivity(spec);
  StructuralMap smap = StructuralMap::build(spec);

  if (ctrl.method_code == 0) {
    solve_jml(spec, smap, ctrl, out);
  } else if (ctrl.method_code == 1) {
    solve_jeffreys(spec, smap, ctrl, out);
  } else if (ctrl.method_code == 2) {
    solve_mml(spec, smap, ctrl, out);
  } else if (ctrl.method_code == 3) {
    solve_cml(spec, smap, ctrl, out);
  }
  return out;
}

inline std::vector<FastFitOutput> run_fast_batch_parallel(
    std::vector<FastDesignSpec>& specs,
    const FastControl& ctrl,
    int n_threads = 1) {
  const int R = static_cast<int>(specs.size());
  std::vector<FastFitOutput> results(static_cast<size_t>(R));
  if (R == 0) return results;
  int workers = std::max(1, std::min(n_threads, R));
  if (workers == 1) {
    for (int r = 0; r < R; ++r) {
      results[static_cast<size_t>(r)] = run_fast_fit(specs[static_cast<size_t>(r)], ctrl);
    }
    return results;
  }

  std::atomic<int> next_idx{0};
  std::vector<std::thread> pool;
  pool.reserve(static_cast<size_t>(workers));
  for (int w = 0; w < workers; ++w) {
    pool.emplace_back([&]() {
      while (true) {
        int r = next_idx.fetch_add(1, std::memory_order_relaxed);
        if (r >= R) break;
        results[static_cast<size_t>(r)] = run_fast_fit(specs[static_cast<size_t>(r)], ctrl);
      }
    });
  }
  for (auto& th : pool) th.join();
  return results;
}

// Parallel per-spec warm-start runner (used by analyze_dff for parallel split-rater refits!)
inline std::vector<FastFitOutput> run_fast_batch_warm_parallel(
    std::vector<FastDesignSpec>& specs,
    const std::vector<FastControl>& ctrls,
    int n_threads = 1) {
  const int R = static_cast<int>(specs.size());
  std::vector<FastFitOutput> results(static_cast<size_t>(R));
  if (R == 0) return results;
  int workers = std::max(1, std::min(n_threads, R));
  if (workers == 1) {
    for (int r = 0; r < R; ++r) {
      results[static_cast<size_t>(r)] = run_fast_fit(specs[static_cast<size_t>(r)], ctrls[static_cast<size_t>(r)]);
    }
    return results;
  }

  std::atomic<int> next_idx{0};
  std::vector<std::thread> pool;
  pool.reserve(static_cast<size_t>(workers));
  for (int w = 0; w < workers; ++w) {
    pool.emplace_back([&]() {
      while (true) {
        int r = next_idx.fetch_add(1, std::memory_order_relaxed);
        if (r >= R) break;
        results[static_cast<size_t>(r)] = run_fast_fit(specs[static_cast<size_t>(r)], ctrls[static_cast<size_t>(r)]);
      }
    });
  }
  for (auto& th : pool) th.join();
  return results;
}

} // namespace mfrm_fast

#endif // MFRMR_FAST_SOLVERS_H

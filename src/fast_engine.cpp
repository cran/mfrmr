#define R_NO_REMAP
#include <cpp11.hpp>
#include "fast_solvers.h"

using namespace cpp11;
namespace writable = cpp11::writable;

namespace {

inline double to_r_double(double v) {
  return std::isnan(v) ? NA_REAL : v;
}

inline double from_r_double(double v) {
  return (cpp11::is_na(v) || std::isnan(v)) ? mfrm_fast::kNaN : v;
}

inline bool list_has_nonnull(const list& lst, const char* key) {
  SEXP nms = lst.names();
  if (nms == R_NilValue) return false;
  strings s_nms(nms);
  for (R_xlen_t i = 0; i < s_nms.size(); ++i) {
    if (std::string(s_nms[i]) == key) {
      return lst[i] != R_NilValue;
    }
  }
  return false;
}

writable::doubles vec_to_r(const std::vector<double>& v) {
  writable::doubles out(static_cast<R_xlen_t>(v.size()));
  for (size_t i = 0; i < v.size(); ++i) {
    out[static_cast<R_xlen_t>(i)] = to_r_double(v[i]);
  }
  return out;
}

writable::integers ivec_to_r(const std::vector<int>& v) {
  writable::integers out(static_cast<R_xlen_t>(v.size()));
  for (size_t i = 0; i < v.size(); ++i) {
    out[static_cast<R_xlen_t>(i)] = v[i];
  }
  return out;
}

writable::doubles_matrix<> rowmajor_to_rmat(const std::vector<double>& v, int nr, int nc) {
  if (nr <= 0 || nc <= 0) {
    return writable::doubles_matrix<>(0, 0);
  }
  writable::doubles_matrix<> out(nr, nc);
  if (static_cast<int>(v.size()) != nr * nc) {
    for (int i = 0; i < nr; ++i) {
      for (int j = 0; j < nc; ++j) out(i, j) = NA_REAL;
    }
    return out;
  }
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      out(i, j) = to_r_double(v[static_cast<size_t>(i) * nc + j]);
    }
  }
  return out;
}

mfrm_fast::FastDesignSpec parse_design_spec_list(const list& ds) {
  mfrm_fast::FastDesignSpec spec;
  integers person0(ds["person0"]);
  integers score0(ds["score0"]);
  doubles weight(ds["weight"]);
  integers facet_n_levels(ds["facet_n_levels"]);
  doubles facet_signs(ds["facet_signs"]);
  integers facet_dummy(ds["facet_dummy"]);
  integers facet_noncenter(ds["facet_noncenter"]);
  list anchor_fixed_list(ds["anchor_fixed"]);
  list anchor_val_list(ds["anchor_val"]);
  integers_matrix<> facet_mat0(ds["facet_mat0"]);

  spec.n_obs = person0.size();
  spec.n_person = as_cpp<int>(ds["n_person"]);
  spec.n_facets = facet_n_levels.size();
  spec.n_cat = as_cpp<int>(ds["n_cat"]);
  spec.rating_min = as_cpp<double>(ds["rating_min"]);
  spec.model_code = as_cpp<int>(ds["model_code"]);
  spec.step_facet = as_cpp<int>(ds["step_facet_idx"]);
  spec.centering_scope = list_has_nonnull(ds, "centering_scope")
                             ? as_cpp<int>(ds["centering_scope"])
                             : 0;

  spec.facet_n_levels.resize(static_cast<size_t>(spec.n_facets));
  spec.facet_signs.resize(static_cast<size_t>(spec.n_facets));
  spec.facet_dummy.resize(static_cast<size_t>(spec.n_facets));
  spec.facet_noncenter.resize(static_cast<size_t>(spec.n_facets));
  spec.anchor_fixed.resize(static_cast<size_t>(spec.n_facets));
  spec.anchor_val.resize(static_cast<size_t>(spec.n_facets));

  for (int f = 0; f < spec.n_facets; ++f) {
    int Lf = facet_n_levels[f];
    spec.facet_n_levels[static_cast<size_t>(f)] = Lf;
    spec.facet_signs[static_cast<size_t>(f)] = facet_signs[f];
    spec.facet_dummy[static_cast<size_t>(f)] = facet_dummy[f];
    spec.facet_noncenter[static_cast<size_t>(f)] = facet_noncenter[f];
    integers af(anchor_fixed_list[f]);
    doubles av(anchor_val_list[f]);
    spec.anchor_fixed[static_cast<size_t>(f)].assign(af.begin(), af.end());
    spec.anchor_val[static_cast<size_t>(f)].assign(av.begin(), av.end());
  }

  // Optional custom affine facet design maps (explanatory RSM/PCM & joint interactions)
  if (list_has_nonnull(ds, "has_custom_facet_map")) {
    integers has_cf(ds["has_custom_facet_map"]);
    integers cf_q(ds["custom_facet_q"]);
    list cf_d0_list(ds["custom_facet_d0"]);
    list cf_B_list(ds["custom_facet_B"]);
    spec.has_custom_facet_map.assign(has_cf.begin(), has_cf.end());
    spec.custom_facet_q.assign(cf_q.begin(), cf_q.end());
    spec.custom_facet_d0.resize(static_cast<size_t>(spec.n_facets));
    spec.custom_facet_B.resize(static_cast<size_t>(spec.n_facets));
    for (int f = 0; f < spec.n_facets; ++f) {
      if (spec.has_custom_facet_map[static_cast<size_t>(f)]) {
        int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
        int qf = spec.custom_facet_q[static_cast<size_t>(f)];
        doubles d0_f(cf_d0_list[f]);
        spec.custom_facet_d0[static_cast<size_t>(f)].resize(static_cast<size_t>(Lf));
        for (int l = 0; l < Lf; ++l) {
          spec.custom_facet_d0[static_cast<size_t>(f)][static_cast<size_t>(l)] = from_r_double(d0_f[l]);
        }
        spec.custom_facet_B[static_cast<size_t>(f)].assign(static_cast<size_t>(Lf) * std::max(0, qf), 0.0);
        if (qf > 0) {
          doubles_matrix<> B_mat(cf_B_list[f]);
          for (int l = 0; l < Lf; ++l) {
            for (int c = 0; c < qf; ++c) {
              spec.custom_facet_B[static_cast<size_t>(f)][static_cast<size_t>(l) * qf + c] = B_mat(l, c);
            }
          }
        }
      }
    }
    if (list_has_nonnull(ds, "custom_facet_anchored")) {
      list cf_anc_list(ds["custom_facet_anchored"]);
      spec.custom_facet_anchored.resize(static_cast<size_t>(spec.n_facets));
      for (int f = 0; f < spec.n_facets; ++f) {
        if (cf_anc_list[f] != R_NilValue) {
          integers ca(cf_anc_list[f]);
          spec.custom_facet_anchored[static_cast<size_t>(f)].assign(ca.begin(), ca.end());
        }
      }
    }
  }

  // Optional GPCM / GMFRM log-slope affine design maps
  if (list_has_nonnull(ds, "n_slope_facets") && as_cpp<int>(ds["n_slope_facets"]) > 0) {
    spec.n_slope_facets = as_cpp<int>(ds["n_slope_facets"]);
    spec.slope_formula = list_has_nonnull(ds, "slope_formula") ? as_cpp<int>(ds["slope_formula"]) : 0;
    integers sf_idx(ds["slope_facets0"]);
    integers sq_vec(ds["slope_q"]);
    list sd0_list(ds["slope_d0"]);
    list sB_list(ds["slope_B"]);
    list sfix_list(ds["slope_fixed"]);
    spec.slope_facets.assign(sf_idx.begin(), sf_idx.end());
    spec.slope_q.assign(sq_vec.begin(), sq_vec.end());
    spec.slope_d0.resize(static_cast<size_t>(spec.n_slope_facets));
    spec.slope_B.resize(static_cast<size_t>(spec.n_slope_facets));
    spec.slope_fixed.resize(static_cast<size_t>(spec.n_slope_facets));
    for (int s = 0; s < spec.n_slope_facets; ++s) {
      int f = spec.slope_facets[static_cast<size_t>(s)];
      int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
      int qs = spec.slope_q[static_cast<size_t>(s)];
      doubles d0_s(sd0_list[s]);
      integers fx_s(sfix_list[s]);
      spec.slope_d0[static_cast<size_t>(s)].assign(d0_s.begin(), d0_s.end());
      spec.slope_fixed[static_cast<size_t>(s)].assign(fx_s.begin(), fx_s.end());
      spec.slope_B[static_cast<size_t>(s)].assign(static_cast<size_t>(Lf) * std::max(0, qs), 0.0);
      if (qs > 0) {
        doubles_matrix<> B_s(sB_list[s]);
        for (int l = 0; l < Lf; ++l) {
          for (int c = 0; c < qs; ++c) {
            spec.slope_B[static_cast<size_t>(s)][static_cast<size_t>(l) * qs + c] = B_s(l, c);
          }
        }
      }
    }
  }

  spec.person.assign(person0.begin(), person0.end());
  spec.score.assign(score0.begin(), score0.end());
  spec.weight.assign(weight.begin(), weight.end());
  spec.facet_flat.resize(static_cast<size_t>(spec.n_obs) * spec.n_facets);
  for (int i = 0; i < spec.n_obs; ++i) {
    for (int f = 0; f < spec.n_facets; ++f) {
      spec.facet_flat[static_cast<size_t>(i) * spec.n_facets + f] = facet_mat0(i, f);
    }
  }

  if (list_has_nonnull(ds, "step_group0") && list_has_nonnull(ds, "n_step_groups")) {
    integers sg0(ds["step_group0"]);
    spec.step_group.assign(sg0.begin(), sg0.end());
    spec.n_step_groups = as_cpp<int>(ds["n_step_groups"]);
  } else {
    spec.step_group.assign(static_cast<size_t>(spec.n_obs), 0);
    spec.n_step_groups =
        (spec.model_code == 1 && spec.step_facet >= 0 && spec.step_facet < spec.n_facets)
            ? spec.facet_n_levels[static_cast<size_t>(spec.step_facet)]
            : 1;
    if (spec.model_code == 1 && spec.step_facet >= 0 && spec.step_facet < spec.n_facets) {
      for (int i = 0; i < spec.n_obs; ++i) {
        spec.step_group[static_cast<size_t>(i)] =
            spec.facet_flat[static_cast<size_t>(i) * spec.n_facets + spec.step_facet];
      }
    }
  }

  // Optional per-step-group structural category bounds [step_cat_min, step_cat_max]
  if (list_has_nonnull(ds, "step_cat_min")) {
    integers smin(ds["step_cat_min"]);
    spec.step_cat_min.assign(smin.begin(), smin.end());
  }
  if (list_has_nonnull(ds, "step_cat_max")) {
    integers smax(ds["step_cat_max"]);
    spec.step_cat_max.assign(smax.begin(), smax.end());
  }

  // Optional custom affine step map (step_constraints: fixed / shared / linear step basis)
  if (list_has_nonnull(ds, "has_custom_step_map") && as_cpp<bool>(ds["has_custom_step_map"])) {
    spec.has_custom_step_map = true;
    spec.custom_step_q = as_cpp<int>(ds["custom_step_q"]);
    doubles sd0(ds["custom_step_d0"]);
    integers sfix(ds["custom_step_fixed"]);
    spec.custom_step_d0.resize(static_cast<size_t>(sd0.size()));
    for (R_xlen_t r = 0; r < sd0.size(); ++r) {
      spec.custom_step_d0[static_cast<size_t>(r)] = from_r_double(sd0[r]);
    }
    spec.custom_step_fixed.assign(sfix.begin(), sfix.end());
    int n_rows = static_cast<int>(sd0.size());
    int qs = spec.custom_step_q;
    spec.custom_step_B.assign(static_cast<size_t>(n_rows) * std::max(1, qs), 0.0);
    if (qs > 0 && list_has_nonnull(ds, "custom_step_B")) {
      doubles_matrix<> sB(ds["custom_step_B"]);
      for (int r = 0; r < n_rows; ++r) {
        for (int c = 0; c < qs; ++c) {
          spec.custom_step_B[static_cast<size_t>(r) * qs + c] = sB(r, c);
        }
      }
    }
  }

  spec.pop_mode = list_has_nonnull(ds, "pop_mode") ? as_cpp<int>(ds["pop_mode"]) : 0;
  spec.fixed_sigma = list_has_nonnull(ds, "fixed_sigma") ? as_cpp<double>(ds["fixed_sigma"]) : 1.0;
  spec.pop_dim = 0;
  if (spec.pop_mode == 1 && list_has_nonnull(ds, "pop_X_mat")) {
    doubles_matrix<> pop_X_mat(ds["pop_X_mat"]);
    spec.pop_dim = pop_X_mat.ncol();
    if (spec.pop_dim > 0) {
      spec.pop_X.resize(static_cast<size_t>(spec.n_person) * spec.pop_dim);
      for (int p = 0; p < spec.n_person; ++p) {
        for (int c = 0; c < spec.pop_dim; ++c) {
          spec.pop_X[static_cast<size_t>(p) * spec.pop_dim + c] = pop_X_mat(p, c);
        }
      }
    }
  }

  spec.build_person_csr();
  return spec;
}

mfrm_fast::FastControl parse_control_list(const list& ctrl_r) {
  mfrm_fast::FastControl ctrl;
  if (list_has_nonnull(ctrl_r, "method_code")) ctrl.method_code = as_cpp<int>(ctrl_r["method_code"]);
  if (list_has_nonnull(ctrl_r, "maxit")) ctrl.maxit = as_cpp<int>(ctrl_r["maxit"]);
  if (list_has_nonnull(ctrl_r, "reltol")) ctrl.reltol = as_cpp<double>(ctrl_r["reltol"]);
  if (list_has_nonnull(ctrl_r, "quad_points")) ctrl.quad_points = as_cpp<int>(ctrl_r["quad_points"]);
  if (list_has_nonnull(ctrl_r, "adaptive_quad")) ctrl.adaptive_quad = as_cpp<bool>(ctrl_r["adaptive_quad"]);
  if (list_has_nonnull(ctrl_r, "xtreme")) ctrl.xtreme = as_cpp<double>(ctrl_r["xtreme"]);
  if (list_has_nonnull(ctrl_r, "whexact")) ctrl.whexact = as_cpp<bool>(ctrl_r["whexact"]);
  if (list_has_nonnull(ctrl_r, "compute_diagnostics")) ctrl.compute_diagnostics = as_cpp<bool>(ctrl_r["compute_diagnostics"]);
  if (list_has_nonnull(ctrl_r, "wle_primary")) ctrl.wle_primary = as_cpp<bool>(ctrl_r["wle_primary"]);
  if (list_has_nonnull(ctrl_r, "compute_wle")) ctrl.compute_wle = as_cpp<bool>(ctrl_r["compute_wle"]);
  if (list_has_nonnull(ctrl_r, "bias_facet_a")) ctrl.bias_facet_a = as_cpp<int>(ctrl_r["bias_facet_a"]);
  if (list_has_nonnull(ctrl_r, "bias_facet_b")) ctrl.bias_facet_b = as_cpp<int>(ctrl_r["bias_facet_b"]);
  if (list_has_nonnull(ctrl_r, "facet_l2_pen")) ctrl.facet_l2_pen = as_cpp<double>(ctrl_r["facet_l2_pen"]);
  if (list_has_nonnull(ctrl_r, "slope_l2_pen")) ctrl.slope_l2_pen = as_cpp<double>(ctrl_r["slope_l2_pen"]);
  if (list_has_nonnull(ctrl_r, "warm_init")) {
    doubles wi(ctrl_r["warm_init"]);
    ctrl.warm_init.assign(wi.begin(), wi.end());
  }
  return ctrl;
}

writable::list pack_fit_output(const mfrm_fast::FastDesignSpec& spec,
                               const mfrm_fast::FastFitOutput& res) {
  using namespace cpp11::literals;
  const int F = spec.n_facets;
  const int K = spec.n_cat - 1;
  const int G = spec.n_step_groups;
  const int n_st_rows = G * K;
  const int xi_dim = static_cast<int>(res.xi_est.size());

  writable::list f_est(F), f_se(F), f_rse(F), f_jse(F);
  writable::list f_cov_mod(F), f_cov_sw(F), f_diag_list(F), f_sub(F);

  for (int f = 0; f < F; ++f) {
    int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
    f_est[f] = vec_to_r(res.facet_est[static_cast<size_t>(f)]);
    f_se[f] = vec_to_r(res.facet_se[static_cast<size_t>(f)]);
    f_rse[f] = vec_to_r(res.facet_real_se[static_cast<size_t>(f)]);
    f_jse[f] = (f < static_cast<int>(res.facet_joint_se.size()))
                   ? vec_to_r(res.facet_joint_se[static_cast<size_t>(f)])
                   : vec_to_r(std::vector<double>(static_cast<size_t>(Lf), mfrm_fast::kNaN));
    f_cov_mod[f] = (f < static_cast<int>(res.facet_cov_model.size()))
                       ? rowmajor_to_rmat(res.facet_cov_model[static_cast<size_t>(f)], Lf, Lf)
                       : rowmajor_to_rmat({}, Lf, Lf);
    f_cov_sw[f] = (f < static_cast<int>(res.facet_cov_sandwich.size()))
                      ? rowmajor_to_rmat(res.facet_cov_sandwich[static_cast<size_t>(f)], Lf, Lf)
                      : rowmajor_to_rmat({}, Lf, Lf);
    f_sub[f] = ivec_to_r(res.connectivity.facet_subset[static_cast<size_t>(f)]);
    const auto& fd = res.facet_diag[static_cast<size_t>(f)];
    f_diag_list[f] = writable::list({
      "obs_n"_nm = vec_to_r(fd.obs_n),
      "obs_score"_nm = vec_to_r(fd.obs_score),
      "exp_score"_nm = vec_to_r(fd.exp_score),
      "cond_info"_nm = vec_to_r(fd.cond_info),
      "infit_mnsq"_nm = vec_to_r(fd.infit_mnsq),
      "infit_zstd"_nm = vec_to_r(fd.infit_zstd),
      "outfit_mnsq"_nm = vec_to_r(fd.outfit_mnsq),
      "outfit_zstd"_nm = vec_to_r(fd.outfit_zstd),
      "fair_m"_nm = vec_to_r(fd.fair_m),
      "fair_z"_nm = vec_to_r(fd.fair_z)
    });
  }

  writable::list p_diag_r({
    "obs_n"_nm = vec_to_r(res.person_diag.obs_n),
    "obs_score"_nm = vec_to_r(res.person_diag.obs_score),
    "exp_score"_nm = vec_to_r(res.person_diag.exp_score),
    "cond_info"_nm = vec_to_r(res.person_diag.cond_info),
    "infit_mnsq"_nm = vec_to_r(res.person_diag.infit_mnsq),
    "infit_zstd"_nm = vec_to_r(res.person_diag.infit_zstd),
    "outfit_mnsq"_nm = vec_to_r(res.person_diag.outfit_mnsq),
    "outfit_zstd"_nm = vec_to_r(res.person_diag.outfit_zstd),
    "fair_m"_nm = vec_to_r(res.person_diag.fair_m),
    "fair_z"_nm = vec_to_r(res.person_diag.fair_z)
  });

  int nb = static_cast<int>(res.local_bias_cells.size());
  writable::integers b_la(nb), b_lb(nb), b_raw_n(nb);
  writable::doubles b_cnt(nb), b_obs(nb), b_exp(nb), b_est(nb), b_se(nb), b_t(nb), b_df(nb);
  for (int i = 0; i < nb; ++i) {
    const auto& c = res.local_bias_cells[static_cast<size_t>(i)];
    b_la[i] = c.level_a;
    b_lb[i] = c.level_b;
    b_raw_n[i] = c.raw_obs_n;
    b_cnt[i] = to_r_double(c.obs_count);
    b_obs[i] = to_r_double(c.obs_score);
    b_exp[i] = to_r_double(c.exp_score);
    b_est[i] = to_r_double(c.bias_est);
    b_se[i] = to_r_double(c.bias_se);
    b_t[i] = to_r_double(c.t_stat);
    b_df[i] = to_r_double(c.df);
  }
  writable::list bias_r({
    "level_a"_nm = b_la, "level_b"_nm = b_lb, "raw_obs_n"_nm = b_raw_n,
    "obs_count"_nm = b_cnt, "obs_score"_nm = b_obs, "exp_score"_nm = b_exp,
    "bias_est"_nm = b_est, "bias_se"_nm = b_se, "t_stat"_nm = b_t, "df"_nm = b_df
  });

  const int S_sl = spec.n_slope_facets;
  writable::list sl_log(S_sl), sl_est(S_sl), sl_se(S_sl), sl_jse(S_sl), sl_cov_mod(S_sl), sl_cov_sw(S_sl);
  for (int s = 0; s < S_sl; ++s) {
    int f = spec.slope_facets[static_cast<size_t>(s)];
    int Lf = spec.facet_n_levels[static_cast<size_t>(f)];
    sl_log[s] = (s < static_cast<int>(res.log_slope_est.size()))
                    ? vec_to_r(res.log_slope_est[static_cast<size_t>(s)])
                    : vec_to_r(std::vector<double>(static_cast<size_t>(Lf), 0.0));
    sl_est[s] = (s < static_cast<int>(res.slope_est.size()))
                    ? vec_to_r(res.slope_est[static_cast<size_t>(s)])
                    : vec_to_r(std::vector<double>(static_cast<size_t>(Lf), 1.0));
    sl_se[s] = (s < static_cast<int>(res.slope_se.size()))
                   ? vec_to_r(res.slope_se[static_cast<size_t>(s)])
                   : vec_to_r(std::vector<double>(static_cast<size_t>(Lf), mfrm_fast::kNaN));
    sl_jse[s] = (s < static_cast<int>(res.slope_joint_se.size()))
                    ? vec_to_r(res.slope_joint_se[static_cast<size_t>(s)])
                    : vec_to_r(std::vector<double>(static_cast<size_t>(Lf), mfrm_fast::kNaN));
    sl_cov_mod[s] = (s < static_cast<int>(res.slope_cov_model.size()))
                        ? rowmajor_to_rmat(res.slope_cov_model[static_cast<size_t>(s)], Lf, Lf)
                        : rowmajor_to_rmat({}, Lf, Lf);
    sl_cov_sw[s] = (s < static_cast<int>(res.slope_cov_sandwich.size()))
                       ? rowmajor_to_rmat(res.slope_cov_sandwich[static_cast<size_t>(s)], Lf, Lf)
                       : rowmajor_to_rmat({}, Lf, Lf);
  }

  return writable::list({
    "converged"_nm = res.converged,
    "iterations"_nm = res.iterations,
    "loglik"_nm = to_r_double(res.loglik),
    "objective"_nm = to_r_double(res.objective),
    "max_abs_grad"_nm = to_r_double(res.max_abs_grad),
    "n_free_par"_nm = res.n_free_par,
    "connected"_nm = res.connectivity.connected,
    "n_subsets"_nm = res.connectivity.n_subsets,
    "interior_categories_supported"_nm = res.connectivity.interior_categories_supported,
    "person_subset"_nm = ivec_to_r(res.connectivity.person_subset),
    "facet_subset"_nm = f_sub,
    "person_extreme_flag"_nm = ivec_to_r(res.connectivity.person_extreme_flag),
    "person_est"_nm = vec_to_r(res.person_est),
    "person_se"_nm = vec_to_r(res.person_se),
    "person_real_se"_nm = vec_to_r(res.person_real_se),
    "person_cond_se"_nm = vec_to_r(res.person_cond_se),
    "person_eap"_nm = vec_to_r(res.person_eap),
    "person_psd"_nm = vec_to_r(res.person_psd),
    "person_wle"_nm = vec_to_r(res.person_wle),
    "person_wle_se"_nm = vec_to_r(res.person_wle_se),
    "person_diag"_nm = p_diag_r,
    "facet_est"_nm = f_est,
    "facet_se"_nm = f_se,
    "facet_real_se"_nm = f_rse,
    "facet_joint_se"_nm = f_jse,
    "facet_cov_model"_nm = f_cov_mod,
    "facet_cov_sandwich"_nm = f_cov_sw,
    "facet_diag"_nm = f_diag_list,
    "step_est"_nm = rowmajor_to_rmat(res.step_est, G, K),
    "step_se"_nm = rowmajor_to_rmat(res.step_se, G, K),
    "step_cov_model"_nm = rowmajor_to_rmat(res.step_cov_model, n_st_rows, n_st_rows),
    "step_cov_sandwich"_nm = rowmajor_to_rmat(res.step_cov_sandwich, n_st_rows, n_st_rows),
    "step_cum"_nm = rowmajor_to_rmat(res.step_cum, G, spec.n_cat),
    "log_slope_est"_nm = sl_log,
    "slope_est"_nm = sl_est,
    "slope_se"_nm = sl_se,
    "slope_joint_se"_nm = sl_jse,
    "slope_cov_model"_nm = sl_cov_mod,
    "slope_cov_sandwich"_nm = sl_cov_sw,
    "alpha_obs"_nm = vec_to_r(res.alpha_obs),
    "xi_est"_nm = vec_to_r(res.xi_est),
    "xi_cov_model"_nm = rowmajor_to_rmat(res.xi_cov_model, xi_dim, xi_dim),
    "xi_cov_sandwich"_nm = rowmajor_to_rmat(res.xi_cov_sandwich, xi_dim, xi_dim),
    "pop_beta"_nm = vec_to_r(res.pop_beta),
    "pop_se_model"_nm = vec_to_r(res.pop_se_model),
    "pop_se_sandwich"_nm = vec_to_r(res.pop_se_sandwich),
    "pop_sigma"_nm = to_r_double(res.pop_sigma),
    "person_rel_eap"_nm = to_r_double(res.person_rel_eap),
    "local_bias"_nm = bias_r
  });
}

} // namespace

[[cpp11::register]]
writable::list mfrm_cpp_fast_connectivity_spec(list ds) {
  using namespace cpp11::literals;
  auto spec = parse_design_spec_list(ds);
  auto conn = mfrm_fast::check_connectivity(spec);
  writable::list f_sub(spec.n_facets);
  for (int f = 0; f < spec.n_facets; ++f) {
    f_sub[f] = ivec_to_r(conn.facet_subset[static_cast<size_t>(f)]);
  }
  return writable::list({
    "connected"_nm = conn.connected,
    "n_subsets"_nm = conn.n_subsets,
    "interior_categories_supported"_nm = conn.interior_categories_supported,
    "person_subset"_nm = ivec_to_r(conn.person_subset),
    "facet_subset"_nm = f_sub,
    "person_extreme_flag"_nm = ivec_to_r(conn.person_extreme_flag),
    "person_obs_count"_nm = vec_to_r(conn.person_obs_count),
    "person_raw_score"_nm = vec_to_r(conn.person_raw_score),
    "person_min_score"_nm = vec_to_r(conn.person_min_score),
    "person_max_score"_nm = vec_to_r(conn.person_max_score)
  });
}

[[cpp11::register]]
writable::list mfrm_cpp_fast_fit_spec(list ds, list ctrl_r) {
  auto spec = parse_design_spec_list(ds);
  auto ctrl = parse_control_list(ctrl_r);
  auto res = mfrm_fast::run_fast_fit(spec, ctrl);
  return pack_fit_output(spec, res);
}

[[cpp11::register]]
writable::list mfrm_cpp_fast_wle(
    int n_person,
    int n_cat,
    integers person0,
    integers score0,
    doubles weight,
    doubles base_eta_obs,
    doubles_matrix<> step_cum_mat,
    integers step_group0,
    integers step_cat_min,
    integers step_cat_max) {
  using namespace cpp11::literals;
  mfrm_fast::FastDesignSpec spec;
  spec.n_obs = person0.size();
  spec.n_person = n_person;
  spec.n_cat = n_cat;
  spec.person.assign(person0.begin(), person0.end());
  spec.score.assign(score0.begin(), score0.end());
  spec.weight.assign(weight.begin(), weight.end());
  spec.step_group.assign(step_group0.begin(), step_group0.end());
  spec.n_step_groups = step_cum_mat.nrow();
  spec.step_cat_min.assign(step_cat_min.begin(), step_cat_min.end());
  spec.step_cat_max.assign(step_cat_max.begin(), step_cat_max.end());
  spec.build_person_csr();

  std::vector<double> be(base_eta_obs.begin(), base_eta_obs.end());
  std::vector<double> sc(static_cast<size_t>(step_cum_mat.nrow()) * n_cat);
  for (int g = 0; g < step_cum_mat.nrow(); ++g) {
    for (int k = 0; k < n_cat; ++k) {
      double v = step_cum_mat(g, k);
      sc[static_cast<size_t>(g) * n_cat + k] =
          (cpp11::is_na(v) || std::isnan(v)) ? mfrm_fast::kInf : v;
    }
  }
  auto wle = mfrm_fast::score_wle_persons(spec, be, sc);
  return writable::list({
    "theta"_nm = vec_to_r(wle.theta),
    "se"_nm = vec_to_r(wle.se)
  });
}

[[cpp11::register]]
writable::list mfrm_cpp_fast_batch_list(
    list dataset_specs,
    list ctrl_r,
    int n_threads) {
  const int R = dataset_specs.size();
  std::vector<mfrm_fast::FastDesignSpec> specs;
  specs.reserve(static_cast<size_t>(R));
  for (int r = 0; r < R; ++r) {
    specs.push_back(parse_design_spec_list(list(dataset_specs[r])));
  }
  auto ctrl = parse_control_list(ctrl_r);
  auto outs = mfrm_fast::run_fast_batch_parallel(specs, ctrl, n_threads);
  writable::list r_outs(R);
  for (int r = 0; r < R; ++r) {
    r_outs[r] = pack_fit_output(specs[static_cast<size_t>(r)], outs[static_cast<size_t>(r)]);
  }
  return r_outs;
}

// Parallel per-spec warm-started batch runner (for parallel split-rater DFF & profile CIs)
[[cpp11::register]]
writable::list mfrm_cpp_fast_batch_warm_list(
    list dataset_specs,
    list ctrl_list,
    int n_threads) {
  const int R = dataset_specs.size();
  std::vector<mfrm_fast::FastDesignSpec> specs;
  std::vector<mfrm_fast::FastControl> ctrls;
  specs.reserve(static_cast<size_t>(R));
  ctrls.reserve(static_cast<size_t>(R));
  for (int r = 0; r < R; ++r) {
    specs.push_back(parse_design_spec_list(list(dataset_specs[r])));
    ctrls.push_back(parse_control_list(list(ctrl_list[r])));
  }
  auto outs = mfrm_fast::run_fast_batch_warm_parallel(specs, ctrls, n_threads);
  writable::list r_outs(R);
  for (int r = 0; r < R; ++r) {
    r_outs[r] = pack_fit_output(specs[static_cast<size_t>(r)], outs[static_cast<size_t>(r)]);
  }
  return r_outs;
}

[[cpp11::register]]
writable::list mfrm_cpp_fast_batch_shared_design(
    list base_ds,
    integers_matrix<> score_mat0,
    list ctrl_r,
    int n_threads) {
  const int R = score_mat0.ncol();
  const int N = score_mat0.nrow();
  auto base_spec = parse_design_spec_list(base_ds);
  std::vector<mfrm_fast::FastDesignSpec> specs(static_cast<size_t>(R), base_spec);
  for (int r = 0; r < R; ++r) {
    for (int i = 0; i < N; ++i) {
      specs[static_cast<size_t>(r)].score[static_cast<size_t>(i)] = score_mat0(i, r);
    }
  }
  auto ctrl = parse_control_list(ctrl_r);
  auto outs = mfrm_fast::run_fast_batch_parallel(specs, ctrl, n_threads);
  writable::list r_outs(R);
  for (int r = 0; r < R; ++r) {
    r_outs[r] = pack_fit_output(specs[static_cast<size_t>(r)], outs[static_cast<size_t>(r)]);
  }
  return r_outs;
}

// Ultra-fast multi-pair Table 13 local bias + Table 14 RajuSA-Welch / classical Welch contrasts
[[cpp11::register]]
writable::list mfrm_cpp_fast_bias_pairwise(
    list ds,
    doubles theta,
    logicals theta_exclude,
    list beta_list,
    doubles_matrix<> step_cum_mat,
    integers pair_fa,
    integers pair_fb,
    bool facets_von_davier,
    doubles alpha_obs) {
  using namespace cpp11::literals;
  auto spec = parse_design_spec_list(ds);
  const int P = spec.n_person;
  const int F = spec.n_facets;
  const int C = spec.n_cat;

  std::vector<double> th_v(static_cast<size_t>(P), mfrm_fast::kNaN);
  std::vector<bool> ex_v(static_cast<size_t>(P), false);
  for (int p = 0; p < P; ++p) {
    th_v[static_cast<size_t>(p)] = from_r_double(theta[p]);
    ex_v[static_cast<size_t>(p)] = (theta_exclude[p] == TRUE);
  }
  std::vector<std::vector<double>> beta_v(static_cast<size_t>(F));
  for (int f = 0; f < F; ++f) {
    doubles bf(beta_list[f]);
    beta_v[static_cast<size_t>(f)].resize(static_cast<size_t>(bf.size()));
    for (R_xlen_t l = 0; l < bf.size(); ++l) {
      beta_v[static_cast<size_t>(f)][static_cast<size_t>(l)] = from_r_double(bf[l]);
    }
  }
  std::vector<double> sc_v(static_cast<size_t>(step_cum_mat.nrow()) * C);
  for (int g = 0; g < step_cum_mat.nrow(); ++g) {
    for (int k = 0; k < C; ++k) {
      double v = step_cum_mat(g, k);
      sc_v[static_cast<size_t>(g) * C + k] =
          (cpp11::is_na(v) || std::isnan(v)) ? mfrm_fast::kInf : v;
    }
  }

  std::vector<double> alpha_v(alpha_obs.begin(), alpha_obs.end());
  if (alpha_v.size() != static_cast<size_t>(spec.n_obs)) cpp11::stop("Expected one slope per retained rating.");
  const int n_pairs = pair_fa.size();
  writable::list cell_tables(n_pairs);
  writable::list contrast_tables(n_pairs);

  for (int pr = 0; pr < n_pairs; ++pr) {
    int fa = pair_fa[pr];
    int fb = pair_fb[pr];
    auto cells = mfrm_fast::compute_local_bias_2way(
        spec, th_v, ex_v, beta_v, sc_v, fa, fb, &alpha_v);

    int nc = static_cast<int>(cells.size());
    writable::integers c_la(nc), c_lb(nc), c_rn(nc);
    writable::doubles c_cnt(nc), c_obs(nc), c_exp(nc), c_base(nc), c_targ(nc);
    writable::doubles c_est(nc), c_se(nc), c_t(nc), c_df(nc);

    for (int i = 0; i < nc; ++i) {
      const auto& c = cells[static_cast<size_t>(i)];
      double b0 = beta_v[static_cast<size_t>(fa)][static_cast<size_t>(c.level_a)];
      c_la[i] = c.level_a;
      c_lb[i] = c.level_b;
      c_rn[i] = c.raw_obs_n;
      c_cnt[i] = to_r_double(c.obs_count);
      c_obs[i] = to_r_double(c.obs_score);
      c_exp[i] = to_r_double(c.exp_score);
      c_base[i] = to_r_double(b0);
      c_targ[i] = to_r_double(b0 + c.bias_est);
      c_est[i] = to_r_double(c.bias_est);
      c_se[i] = to_r_double(c.bias_se);
      c_t[i] = to_r_double(c.t_stat);
      c_df[i] = to_r_double(c.df);
    }
    cell_tables[pr] = writable::list({
      "level_a"_nm = c_la, "level_b"_nm = c_lb, "raw_obs_n"_nm = c_rn,
      "obs_count"_nm = c_cnt, "obs_score"_nm = c_obs, "exp_score"_nm = c_exp,
      "base_measure"_nm = c_base, "target_measure"_nm = c_targ,
      "bias_est"_nm = c_est, "bias_se"_nm = c_se, "t_stat"_nm = c_t, "df"_nm = c_df
    });

    // Pairwise contrasts within each level of FacetA across contexts of FacetB (Table 14)
    std::vector<int> p_la, p_lb1, p_lb2;
    std::vector<double> p_m1, p_se1, p_n1, p_m2, p_se2, p_n2;
    std::vector<double> p_diff, p_sed, p_t, p_df, p_df_raw;

    for (int i = 0; i < nc; ++i) {
      for (int j = i + 1; j < nc; ++j) {
        if (cells[static_cast<size_t>(i)].level_a != cells[static_cast<size_t>(j)].level_a) continue;
        const auto& c1 = cells[static_cast<size_t>(i)];
        const auto& c2 = cells[static_cast<size_t>(j)];
        double b0 = beta_v[static_cast<size_t>(fa)][static_cast<size_t>(c1.level_a)];
        double m1 = b0 + c1.bias_est;
        double m2 = b0 + c2.bias_est;
        double diff = c1.bias_est - c2.bias_est;
        double sed = (std::isfinite(c1.bias_se) && std::isfinite(c2.bias_se))
                         ? std::sqrt(c1.bias_se * c1.bias_se + c2.bias_se * c2.bias_se)
                         : mfrm_fast::kNaN;
        double tval = (std::isfinite(sed) && sed > 0.0) ? (diff / sed) : mfrm_fast::kNaN;
        double df_unr = mfrm_fast::kNaN;
        double df_val = mfrm_fast::compute_raju_welch_df(
            c1.bias_se, c1.obs_count, c2.bias_se, c2.obs_count, facets_von_davier, df_unr);

        p_la.push_back(c1.level_a);
        p_lb1.push_back(c1.level_b);
        p_lb2.push_back(c2.level_b);
        p_m1.push_back(m1);
        p_se1.push_back(c1.bias_se);
        p_n1.push_back(c1.obs_count);
        p_m2.push_back(m2);
        p_se2.push_back(c2.bias_se);
        p_n2.push_back(c2.obs_count);
        p_diff.push_back(diff);
        p_sed.push_back(sed);
        p_t.push_back(tval);
        p_df.push_back(df_val);
        p_df_raw.push_back(df_unr);
      }
    }
    contrast_tables[pr] = writable::list({
      "level_a"_nm = ivec_to_r(p_la),
      "level_b1"_nm = ivec_to_r(p_lb1),
      "level_b2"_nm = ivec_to_r(p_lb2),
      "target_1"_nm = vec_to_r(p_m1),
      "se_1"_nm = vec_to_r(p_se1),
      "n_1"_nm = vec_to_r(p_n1),
      "target_2"_nm = vec_to_r(p_m2),
      "se_2"_nm = vec_to_r(p_se2),
      "n_2"_nm = vec_to_r(p_n2),
      "contrast"_nm = vec_to_r(p_diff),
      "contrast_se"_nm = vec_to_r(p_sed),
      "t_stat"_nm = vec_to_r(p_t),
      "df"_nm = vec_to_r(p_df),
      "df_unrounded"_nm = vec_to_r(p_df_raw)
    });
  }

  return writable::list({
    "cells"_nm = cell_tables,
    "contrasts"_nm = contrast_tables
  });
}

[[cpp11::register]]
writable::list mfrm_cpp_raju_welch_df(
    doubles se1,
    doubles n1,
    doubles se2,
    doubles n2,
    bool facets_von_davier) {
  using namespace cpp11::literals;
  R_xlen_t m = se1.size();
  writable::doubles df_out(m), df_raw_out(m);
  for (R_xlen_t i = 0; i < m; ++i) {
    double unr = mfrm_fast::kNaN;
    double d = mfrm_fast::compute_raju_welch_df(
        from_r_double(se1[i]), from_r_double(n1[i]),
        from_r_double(se2[i]), from_r_double(n2[i]),
        facets_von_davier, unr);
    df_out[i] = to_r_double(d);
    df_raw_out[i] = to_r_double(unr);
  }
  return writable::list({
    "df"_nm = df_out,
    "df_unrounded"_nm = df_raw_out
  });
}

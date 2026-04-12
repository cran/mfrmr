#include <cpp11.hpp>

#include <cmath>
#include <limits>
#include <vector>

using namespace cpp11;

namespace {

inline writable::doubles_matrix<> build_p_geq(const writable::doubles_matrix<>& probs) {
  const int n = probs.nrow();
  const int k_cat = probs.ncol();
  const int n_steps = k_cat - 1;
  writable::doubles_matrix<> out(n, std::max(n_steps, 0));
  if (n_steps <= 0) {
    return out;
  }

  for (int i = 0; i < n; ++i) {
    out(i, n_steps - 1) = probs(i, k_cat - 1);
    for (int s = n_steps - 2; s >= 0; --s) {
      out(i, s) = out(i, s + 1) + probs(i, s + 1);
    }
  }

  return out;
}

}  // namespace

[[cpp11::register]]
writable::doubles_matrix<> mfrm_cpp_compute_p_geq(doubles_matrix<> probs) {
  const int n = probs.nrow();
  const int k_cat = probs.ncol();
  const int n_steps = k_cat - 1;
  writable::doubles_matrix<> out(n, std::max(n_steps, 0));
  if (n_steps <= 0) {
    return out;
  }

  for (int i = 0; i < n; ++i) {
    out(i, n_steps - 1) = probs(i, k_cat - 1);
    for (int s = n_steps - 2; s >= 0; --s) {
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

  for (int q = 0; q < n_nodes; ++q) {
    if (include_probs) {
      writable::doubles_matrix<> probs_q(n, k_cat);
      for (int i = 0; i < n; ++i) {
        const int person_row = person_int[i] - 1;
        const int obs_col = score_k[i];
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
        const double log_denom = row_max + std::log(denom_sum);

        double lp = log_num[static_cast<size_t>(obs_col)] - log_denom;
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
        const double log_denom = row_max + std::log(denom_sum);

        double lp = log_num[static_cast<size_t>(obs_col)] - log_denom;
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

  for (int q = 0; q < n_nodes; ++q) {
    if (include_probs) {
      writable::doubles_matrix<> probs_q(n, k_cat);
      for (int i = 0; i < n; ++i) {
        const int person_row = person_int[i] - 1;
        const int obs_col = score_k[i];
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
        const double log_denom = row_max + std::log(denom_sum);

        double lp = log_num[static_cast<size_t>(obs_col)] - log_denom;
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
        const double log_denom = row_max + std::log(denom_sum);

        double lp = log_num[static_cast<size_t>(obs_col)] - log_denom;
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
    double second_moment = 0.0;
    for (int k = 0; k < k_cat; ++k) {
      const double p = posterior_prob(i, k);
      first_moment += p * static_cast<double>(k);
      second_moment += p * static_cast<double>(k * k);
    }
    expected_k[i] = first_moment;
    var_k[i] = second_moment - first_moment * first_moment;
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

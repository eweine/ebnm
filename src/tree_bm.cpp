#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

struct GaussMsg {
  double mean;
  double var;
  double log_const;
};

// Product of Gaussian messages: prod_k N(x; means[k], vars[k]) = const * N(x; m, v)
static GaussMsg combine_messages(const std::vector<double>& means,
                                  const std::vector<double>& vars) {
  double m = means[0];
  double v = vars[0];
  double lc = 0.0;

  for (int k = 1; k < (int)means.size(); ++k) {
    double a = means[k];
    double b = vars[k];
    double vpb = v + b;
    lc += -0.5 * std::log(2.0 * M_PI * vpb) - 0.5 * (m - a) * (m - a) / vpb;
    double v_new = 1.0 / (1.0 / v + 1.0 / b);
    m = v_new * (m / v + a / b);
    v = v_new;
  }

  GaussMsg result;
  result.mean = m;
  result.var = v;
  result.log_const = lc;
  return result;
}

// Upward + downward message passing on the tree (root fixed at 0).
//
// All node indices are 1-based from R; converted to 0-based internally.
//
// Args:
//   ntip          - number of tip nodes (nodes 1..ntip in R indexing)
//   total_nodes   - ntip + nnode
//   root          - root node index (1-based)
//   y             - tip observations (length ntip, ordered 1..ntip)
//   s2            - tip variances
//   edge_len_to_child - edge length from parent to each node (length total_nodes,
//                       NA at root position)
//   children      - R list of length total_nodes; element v is integer vector of
//                   1-based child indices for node v+1
//   internal_post - 1-based internal node indices in postorder
//   internal_pre  - 1-based internal node indices in preorder
//   bm_var        - Brownian motion variance parameter
//
// Returns list with: loglik, post_tip_mean, post_tip_var
// [[Rcpp::export]]
List tree_bm_compute_cpp(
    int ntip,
    int total_nodes,
    IntegerVector roots,
    NumericVector y,
    NumericVector s2,
    NumericVector edge_len_to_child,
    List children,
    IntegerVector internal_post,
    IntegerVector internal_pre,
    double bm_var
) {
  NumericVector up_mean(total_nodes);
  NumericVector up_var(total_nodes);
  NumericVector up_logconst(total_nodes, 0.0);

  // Tip messages: p(y_i | x_i) = N(x_i; y_i, s2_i)
  for (int i = 0; i < ntip; ++i) {
    up_mean[i] = y[i];
    up_var[i] = s2[i];
  }

  // Upward pass (postorder)
  int n_post = internal_post.size();
  for (int ki = 0; ki < n_post; ++ki) {
    int v = internal_post[ki] - 1;
    IntegerVector ch_r = as<IntegerVector>(children[v]);
    int nch = ch_r.size();
    if (nch == 0) continue;

    std::vector<double> msg_means(nch), msg_vars(nch);
    double sum_lc = 0.0;
    for (int j = 0; j < nch; ++j) {
      int c = ch_r[j] - 1;
      msg_means[j] = up_mean[c];
      msg_vars[j] = up_var[c] + bm_var * edge_len_to_child[c];
      sum_lc += up_logconst[c];
    }

    GaussMsg comb = combine_messages(msg_means, msg_vars);
    up_mean[v] = comb.mean;
    up_var[v] = comb.var;
    up_logconst[v] = sum_lc + comb.log_const;
  }

  std::vector<int> roots0(roots.size());
  std::vector<bool> is_root(total_nodes, false);
  for (int k = 0; k < roots.size(); ++k) {
    roots0[k] = roots[k] - 1;
    is_root[roots0[k]] = true;
  }

  double loglik = 0.0;
  for (int k = 0; k < (int)roots0.size(); ++k) {
    int r = roots0[k];
    double m_root = up_mean[r];
    double v_root = up_var[r];
    double lc_root = up_logconst[r];
    loglik += lc_root
      - 0.5 * std::log(2.0 * M_PI * v_root)
      - 0.5 * m_root * m_root / v_root;
  }

  // Downward pass: compute outside (top-down) messages for each node.
  NumericVector out_mean(total_nodes, 0.0);
  NumericVector out_var(total_nodes, 0.0);

  // Root is fixed at 0: its children receive a point-mass-at-0 message
  // propagated through the edge, giving out_var[c] = bm_var * edge_len[c].
  for (int k = 0; k < (int)roots0.size(); ++k) {
    int r = roots0[k];
    IntegerVector root_ch = as<IntegerVector>(children[r]);
    for (int j = 0; j < root_ch.size(); ++j) {
      int c = root_ch[j] - 1;
      out_mean[c] = 0.0;
      out_var[c] = bm_var * edge_len_to_child[c];
    }
  }

  // For all other internal nodes (preorder, skip root)
  int n_pre = internal_pre.size();
  for (int ki = 0; ki < n_pre; ++ki) {
    int p = internal_pre[ki] - 1;
    if (is_root[p]) continue;

    IntegerVector ch_r_vec = as<IntegerVector>(children[p]);
    int nch = ch_r_vec.size();
    if (nch == 0) continue;

    // Upward messages from each child shifted by edge length to p
    std::vector<int> ch_0(nch);
    std::vector<double> cup_means(nch), cup_vars(nch);
    for (int j = 0; j < nch; ++j) {
      int c = ch_r_vec[j] - 1;
      ch_0[j] = c;
      cup_means[j] = up_mean[c];
      cup_vars[j] = up_var[c] + bm_var * edge_len_to_child[c];
    }

    // Outside message for each child = combination of:
    //   - outside message from p (out_mean[p], out_var[p])
    //   - upward messages from all siblings (excluding the child itself)
    // then propagated through child's edge length.
    for (int j = 0; j < nch; ++j) {
      int c = ch_0[j];

      std::vector<double> excl_means, excl_vars;
      excl_means.push_back(out_mean[p]);
      excl_vars.push_back(out_var[p]);
      for (int k = 0; k < nch; ++k) {
        if (k != j) {
          excl_means.push_back(cup_means[k]);
          excl_vars.push_back(cup_vars[k]);
        }
      }

      GaussMsg comb_excl = combine_messages(excl_means, excl_vars);
      out_mean[c] = comb_excl.mean;
      out_var[c] = comb_excl.var + bm_var * edge_len_to_child[c];
    }
  }

  // Posterior marginals for tip nodes:
  //   precision = 1/out_var[i] + 1/up_var[i]
  //   post_mean = post_var * (out_mean[i]/out_var[i] + up_mean[i]/up_var[i])
  NumericVector post_tip_mean(ntip), post_tip_var(ntip);
  for (int i = 0; i < ntip; ++i) {
    double prec = 1.0 / out_var[i] + 1.0 / up_var[i];
    post_tip_var[i] = 1.0 / prec;
    post_tip_mean[i] = post_tip_var[i] * (out_mean[i] / out_var[i] + up_mean[i] / up_var[i]);
  }

  return List::create(
    Named("loglik") = loglik,
    Named("post_tip_mean") = post_tip_mean,
    Named("post_tip_var") = post_tip_var
  );
}

// Upward + downward message passing on the tree for a stationary OU prior.
//
// Transition model along an edge of length ell:
//   x_child | x_parent ~ N(a * x_parent, q)
// where a = exp(-alpha * ell) and q = stationary_var * (1 - a^2).
// Root prior:
//   x_root ~ N(0, stationary_var)
//
// [[Rcpp::export]]
List tree_ou_compute_cpp(
    int ntip,
    int total_nodes,
    IntegerVector roots,
    NumericVector y,
    NumericVector s2,
    NumericVector edge_len_to_child,
    List children,
    IntegerVector internal_post,
    IntegerVector internal_pre,
    double alpha,
    double stationary_var
) {
  NumericVector up_mean(total_nodes);
  NumericVector up_var(total_nodes);
  NumericVector up_logconst(total_nodes, 0.0);

  for (int i = 0; i < ntip; ++i) {
    up_mean[i] = y[i];
    up_var[i] = s2[i];
  }

  int n_post = internal_post.size();
  for (int ki = 0; ki < n_post; ++ki) {
    int v = internal_post[ki] - 1;
    IntegerVector ch_r = as<IntegerVector>(children[v]);
    int nch = ch_r.size();
    if (nch == 0) continue;

    std::vector<double> msg_means(nch), msg_vars(nch);
    double sum_lc = 0.0;
    for (int j = 0; j < nch; ++j) {
      int c = ch_r[j] - 1;
      double ell = edge_len_to_child[c];
      double a = std::exp(-alpha * ell);
      double q = stationary_var * (1.0 - a * a);
      double child_var = up_var[c] + q;
      msg_means[j] = up_mean[c] / a;
      msg_vars[j] = child_var / (a * a);
      sum_lc += up_logconst[c] - std::log(a);
    }

    GaussMsg comb = combine_messages(msg_means, msg_vars);
    up_mean[v] = comb.mean;
    up_var[v] = comb.var;
    up_logconst[v] = sum_lc + comb.log_const;
  }

  std::vector<int> roots0(roots.size());
  double loglik = 0.0;
  for (int k = 0; k < roots.size(); ++k) {
    int r = roots[k] - 1;
    roots0[k] = r;
    std::vector<double> root_means(2), root_vars(2);
    root_means[0] = 0.0;
    root_vars[0] = stationary_var;
    root_means[1] = up_mean[r];
    root_vars[1] = up_var[r];
    GaussMsg root_comb = combine_messages(root_means, root_vars);
    loglik += up_logconst[r] + root_comb.log_const;
  }

  NumericVector out_mean(total_nodes, 0.0);
  NumericVector out_var(total_nodes, 0.0);
  for (int k = 0; k < (int)roots0.size(); ++k) {
    int r = roots0[k];
    out_mean[r] = 0.0;
    out_var[r] = stationary_var;
  }

  int n_pre = internal_pre.size();
  for (int ki = 0; ki < n_pre; ++ki) {
    int p = internal_pre[ki] - 1;
    IntegerVector ch_r_vec = as<IntegerVector>(children[p]);
    int nch = ch_r_vec.size();
    if (nch == 0) continue;

    std::vector<int> ch_0(nch);
    std::vector<double> msg_means(nch), msg_vars(nch);
    for (int j = 0; j < nch; ++j) {
      int c = ch_r_vec[j] - 1;
      double ell = edge_len_to_child[c];
      double a = std::exp(-alpha * ell);
      double q = stationary_var * (1.0 - a * a);
      ch_0[j] = c;
      msg_means[j] = up_mean[c] / a;
      msg_vars[j] = (up_var[c] + q) / (a * a);
    }

    for (int j = 0; j < nch; ++j) {
      int c = ch_0[j];
      std::vector<double> excl_means, excl_vars;
      excl_means.push_back(out_mean[p]);
      excl_vars.push_back(out_var[p]);
      for (int k = 0; k < nch; ++k) {
        if (k != j) {
          excl_means.push_back(msg_means[k]);
          excl_vars.push_back(msg_vars[k]);
        }
      }

      GaussMsg comb_excl = combine_messages(excl_means, excl_vars);
      double ell = edge_len_to_child[c];
      double a = std::exp(-alpha * ell);
      double q = stationary_var * (1.0 - a * a);
      out_mean[c] = a * comb_excl.mean;
      out_var[c] = q + a * a * comb_excl.var;
    }
  }

  NumericVector post_tip_mean(ntip), post_tip_var(ntip);
  for (int i = 0; i < ntip; ++i) {
    double prec = 1.0 / out_var[i] + 1.0 / up_var[i];
    post_tip_var[i] = 1.0 / prec;
    post_tip_mean[i] = post_tip_var[i] *
      (out_mean[i] / out_var[i] + up_mean[i] / up_var[i]);
  }

  return List::create(
    Named("loglik") = loglik,
    Named("post_tip_mean") = post_tip_mean,
    Named("post_tip_var") = post_tip_var
  );
}

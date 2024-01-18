// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "data_gibbs_person_parallel.h"
using namespace Rcpp;

// ----------------------------------------------------------------------------|
// Gibbs sampler for the augmented data that is used in DMH
// ----------------------------------------------------------------------------|
IntegerMatrix data_gibbs(IntegerMatrix data,
                         IntegerVector no_categories,
                         NumericMatrix interactions,
                         NumericMatrix thresholds,
                         int iter = 1,
                         bool parallel = false) {
  int no_states = data.nrow();
  int no_nodes = data.ncol();
  IntegerMatrix augmented_data = Rcpp::clone(data);
  int max_no_categories = 0;
  for(int node = 0; node < no_nodes; node++) {
    if(no_categories[node] > max_no_categories) {
      max_no_categories = no_categories[node];
    }
  }

  //The Gibbs sampler ----------------------------------------------------------
  if (parallel) {
    data_gibbs_person_parallel(interactions,
                               thresholds,
                               data,
                               no_categories,
                               no_states,
                               iter,
                               no_nodes,
                               max_no_categories,
                               augmented_data);
  } else {
    data_gibbs_person_serial(interactions,
                             thresholds,
                             data,
                             no_categories,
                             no_states,
                             iter,
                             no_nodes,
                             max_no_categories,
                             augmented_data);
  }

  return augmented_data;
}


// ----------------------------------------------------------------------------|
// DMH algorithm to sample from the full-conditional of the threshold parameters
// ----------------------------------------------------------------------------|
List dmh_thresholds(NumericMatrix interactions,
                    NumericMatrix thresholds,
                    IntegerMatrix observations,
                    IntegerMatrix O_thresholds,
                    IntegerVector no_categories,
                    int no_persons,
                    int no_nodes,
                    double threshold_alpha,
                    double threshold_beta,
                    NumericMatrix proposal_threshold_sd,
                    int m,
                    int t,
                    bool parallel = false) {

  double log_prob;
  double current_state, proposed_state;
  double U;
  int o_ast;

  IntegerMatrix augmented_data(no_persons, no_nodes);

  //Parameters of adaptive proposals -------------------------------------------
  double phi = .75;
  double target_ar = 0.234;
  double epsilon_lo = 1 / no_persons;
  double epsilon_hi = 2.0;

  for(int node = 0; node < no_nodes; node++) {
    for(int category = 0; category < no_categories[node]; category++) {
      current_state = thresholds(node, category);
      proposed_state = R::rnorm(current_state,
                                proposal_threshold_sd(node, category));

      //We set it to current_state if we reject the proposed state
      thresholds(node, category) = proposed_state;

      augmented_data = data_gibbs(observations,
                                  no_categories,
                                  interactions,
                                  thresholds,
                                  m,
                                  parallel);

      o_ast = 0;
      for(int v = 0; v < no_persons; v++) {
        if(augmented_data(v, node) == category + 1) {
          o_ast++;
        }
      }

      log_prob = (O_thresholds(category + 1, node) - o_ast + threshold_alpha) *
        (proposed_state - current_state);

      //Second, we add the ratio of prior probabilities
      log_prob -= (threshold_alpha + threshold_beta) *
        std::log(1 + std::exp(proposed_state));
      log_prob += (threshold_alpha + threshold_beta) *
        std::log(1 + std::exp(current_state));

      U = std::log(R::unif_rand());
      if(U > log_prob) {
        //We changed it to the proposed state for generating the augmented data
        thresholds(node, category) = current_state;
      }

      if(log_prob > 0) {
        log_prob = 1;
      } else {
        log_prob = std::exp(log_prob);
      }
      proposal_threshold_sd(node, category) =
        proposal_threshold_sd(node, category) +
        (log_prob - target_ar) * std::exp(-log(t) * phi);
      if(proposal_threshold_sd(node, category) < epsilon_lo) {
        proposal_threshold_sd(node, category) = epsilon_lo;
      } else if (proposal_threshold_sd(node, category) > epsilon_hi) {
        proposal_threshold_sd(node, category) = epsilon_hi;
      }
    }
  }

  return List::create(Named("thresholds") = thresholds,
                      Named("proposal_threshold_sd") = proposal_threshold_sd);
}

// ----------------------------------------------------------------------------|
// DMH algorithm to sample from the full-conditional of the active interaction
//  parameters (using a cauchy prior) for Bayesian edge selection
// ----------------------------------------------------------------------------|
List dmh_interactions_cauchy(NumericMatrix interactions,
                             NumericMatrix thresholds,
                             IntegerMatrix gamma,
                             IntegerMatrix observations,
                             IntegerMatrix O_interactions,
                             IntegerVector no_categories,
                             int no_persons,
                             int no_nodes,
                             double cauchy_scale,
                             NumericMatrix proposal_interaction_sd,
                             int m,
                             int t,
                             const bool parallel) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  IntegerMatrix augmented_data(no_persons, no_nodes);

  //Parameters of adaptive proposals -------------------------------------------
  double phi = .75;
  double target_ar = 0.234;
  double epsilon_lo = 1 / no_persons;
  double epsilon_hi = 2.0;

  for(int node1 = 0; node1 <  no_nodes - 1; node1++) {
    for(int node2 = node1 + 1; node2 <  no_nodes; node2++) {
      if(gamma(node1, node2) == 1) {
        current_state = interactions(node1, node2);
        proposed_state = R::rnorm(current_state,
                                  proposal_interaction_sd(node1, node2));

        //We set it to current_state if we reject the proposed state
        interactions(node1, node2) = proposed_state;
        interactions(node2, node1) = proposed_state;

        augmented_data = data_gibbs(observations,
                                    no_categories,
                                    interactions,
                                    thresholds,
                                    m,
                                    parallel);

        int o_ast = 0;
        for(int v = 0; v < no_persons; v++) {
          o_ast += augmented_data(v, node1) * augmented_data(v, node2);
        }

        log_prob = (O_interactions(node1, node2) - o_ast) *
          (proposed_state - current_state);

        //Second, we add the ratio of prior probabilities
        log_prob += R::dcauchy(proposed_state, 0.0, cauchy_scale, true);
        log_prob -= R::dcauchy(current_state, 0.0, cauchy_scale, true);

        U = R::unif_rand();
        if(std::log(U) > log_prob) {
          //We changed it to the proposed state for generating augmented data
          interactions(node1, node2) = current_state;
          interactions(node2, node1) = current_state;
        }

        if(log_prob > 0) {
          log_prob = 1;
        } else {
          log_prob = std::exp(log_prob);
        }
        proposal_interaction_sd(node1, node2) = proposal_interaction_sd(node1, node2) +
          (log_prob - target_ar) * std::exp(-log(t) * phi);
        if(proposal_interaction_sd(node1, node2) < epsilon_lo) {
          proposal_interaction_sd(node1, node2) = epsilon_lo;
        } else if (proposal_interaction_sd(node1, node2) > epsilon_hi) {
          proposal_interaction_sd(node1, node2) = epsilon_hi;
        }
      }
    }
  }

  return List::create(Named("interactions") = interactions,
                      Named("proposal_interaction_sd") = proposal_interaction_sd);
}

// ----------------------------------------------------------------------------|
// DMH algorithm to sample from the full-conditional of an edge + interaction
//  pair (using a cauchy prior) for Bayesian edge selection
// ----------------------------------------------------------------------------|
List dmh_edge_interaction_pair_cauchy(NumericMatrix interactions,
                                      NumericMatrix thresholds,
                                      IntegerMatrix gamma,
                                      IntegerMatrix observations,
                                      IntegerMatrix O_interactions,
                                      IntegerVector no_categories,
                                      int no_interactions,
                                      int no_persons,
                                      IntegerMatrix index,
                                      double cauchy_scale,
                                      NumericMatrix theta,
                                      NumericMatrix proposal_interaction_sd,
                                      int m,
                                      const bool parallel) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  int node1;
  int node2;

  int no_nodes = interactions.nrow();
  IntegerMatrix augmented_data(no_persons, no_nodes);

  for(int cntr = 0; cntr < no_interactions; cntr ++) {
    node1 = index(cntr, 1) - 1;
    node2 = index(cntr, 2) - 1;

    current_state = interactions(node1, node2);

    if(gamma(node1, node2) == 0) {
      proposed_state = R::rnorm(current_state, proposal_interaction_sd(node1, node2));
    } else {
      proposed_state = 0.0;
    }

    //We set it to current_state if we reject the proposed state
    interactions(node1, node2) = proposed_state;
    interactions(node2, node1) = proposed_state;

    augmented_data = data_gibbs(observations,
                                no_categories,
                                interactions,
                                thresholds,
                                m,
                                parallel);

    int o_ast = 0;
    for(int v = 0; v < no_persons; v++) {
      o_ast += augmented_data(v, node1) * augmented_data(v, node2);
    }
    log_prob = (O_interactions(node1, node2) - o_ast) *
      (proposed_state - current_state);

    if(gamma(node1, node2) == 0) {
      log_prob += R::dcauchy(proposed_state, 0.0, cauchy_scale, true);
      log_prob -= R::dnorm(proposed_state,
                           current_state,
                           proposal_interaction_sd(node1, node2),
                           true);

      log_prob += log(theta(node1, node2) / (1 - theta(node1, node2)));
    } else {
      log_prob -= R::dcauchy(current_state, 0.0, cauchy_scale,  true);
      log_prob += R::dnorm(current_state,
                           proposed_state,
                           proposal_interaction_sd(node1, node2),
                           true);

      log_prob -= log(theta(node1, node2) / (1 - theta(node1, node2)));
    }

    U = R::unif_rand();
    if(std::log(U) < log_prob) {
      gamma(node1, node2) = 1 - gamma(node1, node2);
      gamma(node2, node1) = 1 - gamma(node2, node1);
    } else {
      interactions(node1, node2) = current_state;
      interactions(node2, node1) = current_state;
    }
  }
  return List::create(Named("interactions") = interactions,
                      Named("gamma") = gamma);
}

// ----------------------------------------------------------------------------|
// A Gibbs step for graphical model parameters for Bayesian edge selection
// ----------------------------------------------------------------------------|
List dmh_gibbs_step_gm(IntegerMatrix observations,
                       IntegerMatrix O_thresholds,
                       IntegerMatrix O_interactions,
                       IntegerVector no_categories,
                       int no_persons,
                       int no_nodes,
                       int no_interactions,
                       int no_thresholds,
                       IntegerMatrix index,
                       NumericMatrix proposal_threshold_sd,
                       NumericMatrix proposal_interaction_sd,
                       int t,
                       int m,
                       double cauchy_scale,
                       double threshold_alpha,
                       double threshold_beta,
                       IntegerMatrix gamma,
                       NumericMatrix interactions,
                       NumericMatrix thresholds,
                       NumericMatrix theta,
                       const bool parallel) {

  if(true) {
    List out = dmh_edge_interaction_pair_cauchy(interactions,
                                                thresholds,
                                                gamma,
                                                observations,
                                                O_interactions,
                                                no_categories,
                                                no_interactions,
                                                no_persons,
                                                index,
                                                cauchy_scale,
                                                theta,
                                                proposal_interaction_sd,
                                                m,
                                                parallel);
    IntegerMatrix gamma = out["gamma"];
    NumericMatrix interactions = out["interactions"];
  }

  //Update interactions (within model move)
  if(true) {
    List out = dmh_interactions_cauchy(interactions,
                                       thresholds,
                                       gamma,
                                       observations,
                                       O_interactions,
                                       no_categories,
                                       no_persons,
                                       no_nodes,
                                       cauchy_scale,
                                       proposal_interaction_sd,
                                       m,
                                       t,
                                       parallel);

    NumericMatrix interactions = out["interactions"];
    NumericMatrix proposal_interaction_sd = out["proposal_interaction_sd"];
  }

  //Update thresholds
  if(true) {
    List out = dmh_thresholds(interactions,
                              thresholds,
                              observations,
                              O_thresholds,
                              no_categories,
                              no_persons,
                              no_nodes,
                              threshold_alpha,
                              threshold_beta,
                              proposal_threshold_sd,
                              m,
                              t,
                              parallel);
    NumericMatrix thresholds = out["thresholds"];
    NumericMatrix proposal_threshold_sd = out["proposal_threshold_sd"];
  }

  return List::create(Named("gamma") = gamma,
                      Named("interactions") = interactions,
                      Named("thresholds") = thresholds,
                      Named("proposal_threshold_sd") = proposal_threshold_sd,
                      Named("proposal_interaction_sd") = proposal_interaction_sd);
}


// ----------------------------------------------------------------------------|
// The Gibbs sampler for Bayesian edge selection
// ----------------------------------------------------------------------------|
// [[Rcpp::export]]
List dmh_gibbs_sampler(IntegerMatrix observations,
                       IntegerMatrix O_thresholds,
                       IntegerMatrix O_interactions,
                       int m,
                       IntegerVector no_categories,
                       IntegerMatrix Index,
                       NumericMatrix proposal_threshold_sd,
                       NumericMatrix proposal_interaction_sd,
                       double cauchy_scale,
                       double threshold_alpha,
                       double threshold_beta,
                       String edge_prior,
                       NumericMatrix theta,
                       double beta_bernoulli_alpha,
                       double beta_bernoulli_beta,
                       IntegerMatrix gamma,
                       NumericMatrix interactions,
                       NumericMatrix thresholds,
                       int iter,
                       int burnin,
                       bool save = false,
                       bool display_progress = false,
                       bool parallel = false) {
  int cntr;
  int no_nodes = observations.ncol();
  int no_persons = observations.nrow();
  int no_interactions = Index.nrow();
  int no_thresholds = sum(no_categories);
  int max_no_categories = max(no_categories);

  IntegerVector v = seq(0, no_interactions - 1);
  IntegerVector order(no_interactions);
  IntegerMatrix index(no_interactions, 3);

  //The resizing based on ``save'' could probably be prettier ------------------
  int nrow = no_nodes;
  int ncol_edges = no_nodes;
  int ncol_thresholds = max_no_categories;

  if(save == true) {
    nrow = iter;
    ncol_edges= no_interactions;
    ncol_thresholds = no_thresholds;
  }

  NumericMatrix out_gamma(nrow, ncol_edges);
  NumericMatrix out_interactions(nrow, ncol_edges);
  NumericMatrix out_thresholds(nrow, ncol_thresholds);

  //Progress bar
  Progress p(iter + burnin, display_progress);

  //The Gibbs sampler ----------------------------------------------------------
  //First, we do burn-in iterations---------------------------------------------
  for(int iteration = 0; iteration < burnin; iteration++) {
    if (Progress::check_abort()) {
      return List::create(Named("gamma") = out_gamma,
                          Named("interactions") = out_interactions,
                          Named("thresholds") = out_thresholds);
    }
    p.increment();

    //Update interactions and model (between model move) -----------------------
    //Use a random order to update the edge - interaction pairs ----------------
    order = sample(v,
                   no_interactions,
                   false,
                   R_NilValue);

    for(int cntr = 0; cntr < no_interactions; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }

    List out = dmh_gibbs_step_gm(observations,
                                 O_thresholds,
                                 O_interactions,
                                 no_categories,
                                 no_persons,
                                 no_nodes,
                                 no_interactions,
                                 no_thresholds,
                                 index,
                                 proposal_threshold_sd,
                                 proposal_interaction_sd,
                                 iteration + 1,
                                 m,
                                 cauchy_scale,
                                 threshold_alpha,
                                 threshold_beta,
                                 gamma,
                                 interactions,
                                 thresholds,
                                 theta,
                                 parallel);

    IntegerMatrix gamma = out["gamma"];
    NumericMatrix interactions = out["interactions"];
    NumericMatrix thresholds = out["thresholds"];
    NumericMatrix proposal_threshold_sd = out["proposal_threshold_sd"];
    NumericMatrix proposal_interaction_sd = out["proposal_interaction_sd"];

    if(edge_prior == "Beta-Bernoulli") {
      int sumG = 0;
      for(int i = 0; i < no_nodes - 1; i++) {
        for(int j = i + 1; j < no_nodes; j++) {
          sumG += gamma(i, j);
        }
      }
      double probability = R::rbeta(beta_bernoulli_alpha + sumG,
                                    beta_bernoulli_beta + no_interactions - sumG);

      for(int i = 0; i < no_nodes - 1; i++) {
        for(int j = i + 1; j < no_nodes; j++) {
          theta(i, j) = probability;
          theta(j, i) = probability;
        }
      }
    }

  }

  //The post burn-in iterations ------------------------------------------------
  for(int iteration = 0; iteration < iter; iteration++) {
    if (Progress::check_abort()) {
      return List::create(Named("gamma") = out_gamma,
                          Named("interactions") = out_interactions,
                          Named("thresholds") = out_thresholds);
    }
    p.increment();

    //Update interactions and model (between model move) -----------------------
    //Use a random order to update the edge - interaction pairs ----------------
    order = sample(v,
                   no_interactions,
                   false,
                   R_NilValue);

    for(int cntr = 0; cntr < no_interactions; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }

    List out = dmh_gibbs_step_gm(observations,
                                 O_thresholds,
                                 O_interactions,
                                 no_categories,
                                 no_persons,
                                 no_nodes,
                                 no_interactions,
                                 no_thresholds,
                                 index,
                                 proposal_threshold_sd,
                                 proposal_interaction_sd,
                                 iteration + 1,
                                 m,
                                 cauchy_scale,
                                 threshold_alpha,
                                 threshold_beta,
                                 gamma,
                                 interactions,
                                 thresholds,
                                 theta,
                                 parallel);

    IntegerMatrix gamma = out["gamma"];
    NumericMatrix interactions = out["interactions"];
    NumericMatrix thresholds = out["thresholds"];
    NumericMatrix proposal_threshold_sd = out["proposal_threshold_sd"];
    NumericMatrix proposal_interaction_sd = out["proposal_interaction_sd"];

    if(edge_prior == "Beta-Bernoulli") {
      int sumG = 0;
      for(int i = 0; i < no_nodes - 1; i++) {
        for(int j = i + 1; j < no_nodes; j++) {
          sumG += gamma(i, j);
        }
      }
      double probability = R::rbeta(beta_bernoulli_alpha + sumG,
                                    beta_bernoulli_beta + no_interactions - sumG);

      for(int i = 0; i < no_nodes - 1; i++) {
        for(int j = i + 1; j < no_nodes; j++) {
          theta(i, j) = probability;
          theta(j, i) = probability;
        }
      }
    }

    //Output -------------------------------------------------------------------
    if(save == TRUE) {
      //Save raw samples -------------------------------------------------------
      cntr = 0;
      for(int node1 = 0; node1 < no_nodes - 1; node1++) {
        for(int node2 = node1 + 1; node2 < no_nodes;node2++) {
          out_gamma(iteration, cntr) = gamma(node1, node2);
          out_interactions(iteration, cntr) = interactions(node1, node2);
          cntr++;
        }
      }
      cntr = 0;
      for(int node = 0; node < no_nodes; node++) {
        for(int category = 0; category < no_categories[node]; category++) {
          out_thresholds(iteration, cntr) = thresholds(node, category);
          cntr++;
        }
      }
    } else {
      //Compute running averages -----------------------------------------------
      for(int node1 = 0; node1 < no_nodes - 1; node1++) {
        for(int node2 = node1 + 1; node2 < no_nodes; node2++) {
          out_gamma(node1, node2) *= iteration;
          out_gamma(node1, node2) += gamma(node1, node2);
          out_gamma(node1, node2) /= iteration + 1;
          out_gamma(node2, node1) = out_gamma(node1, node2);

          out_interactions(node1, node2) *= iteration;
          out_interactions(node1, node2) += interactions(node1, node2);
          out_interactions(node1, node2) /= iteration + 1;
          out_interactions(node2, node1) = out_interactions(node1, node2);
        }

        for(int category = 0; category < no_categories[node1]; category++) {
          out_thresholds(node1, category) *= iteration;
          out_thresholds(node1, category) += thresholds(node1, category);
          out_thresholds(node1, category) /= iteration + 1;
        }
      }
      for(int category = 0; category < no_categories[no_nodes - 1]; category++) {
        out_thresholds(no_nodes - 1, category) *= iteration;
        out_thresholds(no_nodes - 1, category) +=
          thresholds(no_nodes - 1, category);
        out_thresholds(no_nodes - 1, category) /= iteration + 1;
      }
    }
  }

  return List::create(Named("gamma") = out_gamma,
                      Named("interactions") = out_interactions,
                      Named("thresholds") = out_thresholds);
}


// ----------------------------------------------------------------------------|
// DMH algorithm to sample from the full-conditional of the interaction
//  parameters (using a cauchy prior) for Bayesian estimation
// ----------------------------------------------------------------------------|
List dmh_interactions_cauchy_estimation(NumericMatrix interactions,
                                        NumericMatrix thresholds,
                                        IntegerMatrix observations,
                                        IntegerMatrix O_interactions,
                                        IntegerVector no_categories,
                                        int no_persons,
                                        int no_nodes,
                                        double cauchy_scale,
                                        NumericMatrix proposal_interaction_sd,
                                        int m,
                                        int t,
                                        const bool parallel) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  IntegerMatrix augmented_data(no_persons, no_nodes);

  //Parameters of adaptive proposals -------------------------------------------
  double phi = .75;
  double target_ar = 0.234;
  double epsilon_lo = 1 / no_persons;
  double epsilon_hi = 2.0;

  for(int node1 = 0; node1 <  no_nodes - 1; node1++) {
    for(int node2 = node1 + 1; node2 <  no_nodes; node2++) {
      current_state = interactions(node1, node2);
      proposed_state = R::rnorm(current_state,
                                proposal_interaction_sd(node1, node2));

      //We set it to current_state if we reject the proposed state
      interactions(node1, node2) = proposed_state;
      interactions(node2, node1) = proposed_state;

      augmented_data = data_gibbs(observations,
                                  no_categories,
                                  interactions,
                                  thresholds,
                                  m,
                                  parallel);

      int o_ast = 0;
      for(int v = 0; v < no_persons; v++) {
        o_ast += augmented_data(v, node1) * augmented_data(v, node2);
      }

      log_prob = (O_interactions(node1, node2) - o_ast) *
        (proposed_state - current_state);

      //Second, we add the ratio of prior probabilities
      log_prob += R::dcauchy(proposed_state, 0.0, cauchy_scale, true);
      log_prob -= R::dcauchy(current_state, 0.0, cauchy_scale, true);

      U = R::unif_rand();
      if(std::log(U) > log_prob) {
        //We changed it to the proposed state for generating augmented data
        interactions(node1, node2) = current_state;
        interactions(node2, node1) = current_state;
      }

      if(log_prob > 0) {
        log_prob = 1;
      } else {
        log_prob = std::exp(log_prob);
      }
      proposal_interaction_sd(node1, node2) = proposal_interaction_sd(node1, node2) +
        (log_prob - target_ar) * std::exp(-log(t) * phi);
      if(proposal_interaction_sd(node1, node2) < epsilon_lo) {
        proposal_interaction_sd(node1, node2) = epsilon_lo;
      } else if (proposal_interaction_sd(node1, node2) > epsilon_hi) {
        proposal_interaction_sd(node1, node2) = epsilon_hi;
      }
    }
  }

  return List::create(Named("interactions") = interactions,
                      Named("proposal_interaction_sd") = proposal_interaction_sd);
}

// ----------------------------------------------------------------------------|
// A Gibbs step for graphical model parameters for Bayesian estimation
// ----------------------------------------------------------------------------|
List dmh_gibbs_step_gm_estimation(IntegerMatrix observations,
                                  IntegerMatrix O_thresholds,
                                  IntegerMatrix O_interactions,
                                  IntegerVector no_categories,
                                  int no_persons,
                                  int no_nodes,
                                  int no_interactions,
                                  int no_thresholds,
                                  IntegerMatrix index,
                                  NumericMatrix proposal_threshold_sd,
                                  NumericMatrix proposal_interaction_sd,
                                  int t,
                                  int m,
                                  double cauchy_scale,
                                  double threshold_alpha,
                                  double threshold_beta,
                                  NumericMatrix interactions,
                                  NumericMatrix thresholds,
                                  const bool parallel) {

  //Update interactions (within model move)
  if(true) {
    List out = dmh_interactions_cauchy_estimation(interactions,
                                                  thresholds,
                                                  observations,
                                                  O_interactions,
                                                  no_categories,
                                                  no_persons,
                                                  no_nodes,
                                                  cauchy_scale,
                                                  proposal_interaction_sd,
                                                  m,
                                                  t,
                                                  parallel);

    NumericMatrix interactions = out["interactions"];
    NumericMatrix proposal_interaction_sd = out["proposal_interaction_sd"];
  }

  //Update thresholds
  if(true) {
    List out = dmh_thresholds(interactions,
                              thresholds,
                              observations,
                              O_thresholds,
                              no_categories,
                              no_persons,
                              no_nodes,
                              threshold_alpha,
                              threshold_beta,
                              proposal_threshold_sd,
                              m,
                              t,
                              parallel);
    NumericMatrix thresholds = out["thresholds"];
    NumericMatrix proposal_threshold_sd = out["proposal_threshold_sd"];
  }

  return List::create(Named("interactions") = interactions,
                      Named("thresholds") = thresholds,
                      Named("proposal_threshold_sd") = proposal_threshold_sd,
                      Named("proposal_interaction_sd") = proposal_interaction_sd);
}


// ----------------------------------------------------------------------------|
// The Gibbs sampler for Bayesian estimation
// ----------------------------------------------------------------------------|
// [[Rcpp::export]]
List dmh_gibbs_sampler_estimation(IntegerMatrix observations,
                                  IntegerMatrix O_thresholds,
                                  IntegerMatrix O_interactions,
                                  int m,
                                  IntegerVector no_categories,
                                  IntegerMatrix Index,
                                  NumericMatrix proposal_threshold_sd,
                                  NumericMatrix proposal_interaction_sd,
                                  double cauchy_scale,
                                  double threshold_alpha,
                                  double threshold_beta,
                                  NumericMatrix interactions,
                                  NumericMatrix thresholds,
                                  int iter,
                                  int burnin,
                                  bool save = false,
                                  bool display_progress = false,
                                  bool parallel = false) {
  int cntr;
  int no_nodes = observations.ncol();
  int no_persons = observations.nrow();
  int no_interactions = Index.nrow();
  int no_thresholds = sum(no_categories);
  int max_no_categories = max(no_categories);

  IntegerVector v = seq(0, no_interactions - 1);
  IntegerVector order(no_interactions);
  IntegerMatrix index(no_interactions, 3);

  //The resizing based on ``save'' could probably be prettier ------------------
  int nrow = no_nodes;
  int ncol_edges = no_nodes;
  int ncol_thresholds = max_no_categories;

  if(save == true) {
    nrow = iter;
    ncol_edges= no_interactions;
    ncol_thresholds = no_thresholds;
  }

  NumericMatrix out_interactions(nrow, ncol_edges);
  NumericMatrix out_thresholds(nrow, ncol_thresholds);

  //Progress bar
  Progress p(iter + burnin, display_progress);

  //The Gibbs sampler ----------------------------------------------------------
  //First, we do burn-in iterations---------------------------------------------
  for(int iteration = 0; iteration < burnin; iteration++) {
    if (Progress::check_abort()) {
      return List::create(Named("interactions") = out_interactions,
                          Named("thresholds") = out_thresholds);
    }
    p.increment();

    //Update interactions and model (between model move) -----------------------
    //Use a random order to update the edge - interaction pairs ----------------
    order = sample(v,
                   no_interactions,
                   false,
                   R_NilValue);

    for(int cntr = 0; cntr < no_interactions; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }

    List out = dmh_gibbs_step_gm_estimation(observations,
                                            O_thresholds,
                                            O_interactions,
                                            no_categories,
                                            no_persons,
                                            no_nodes,
                                            no_interactions,
                                            no_thresholds,
                                            index,
                                            proposal_threshold_sd,
                                            proposal_interaction_sd,
                                            iteration + 1,
                                            m,
                                            cauchy_scale,
                                            threshold_alpha,
                                            threshold_beta,
                                            interactions,
                                            thresholds,
                                            parallel);

    NumericMatrix interactions = out["interactions"];
    NumericMatrix thresholds = out["thresholds"];
    NumericMatrix proposal_threshold_sd = out["proposal_threshold_sd"];
    NumericMatrix proposal_interaction_sd = out["proposal_interaction_sd"];
  }

  //The post burn-in iterations ------------------------------------------------
  for(int iteration = 0; iteration < iter; iteration++) {
    if (Progress::check_abort()) {
      return List::create(Named("interactions") = out_interactions,
                          Named("thresholds") = out_thresholds);
    }
    p.increment();

    //Update interactions and model (between model move) -----------------------
    //Use a random order to update the edge - interaction pairs ----------------
    order = sample(v,
                   no_interactions,
                   false,
                   R_NilValue);

    for(int cntr = 0; cntr < no_interactions; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }

    List out = dmh_gibbs_step_gm_estimation(observations,
                                            O_thresholds,
                                            O_interactions,
                                            no_categories,
                                            no_persons,
                                            no_nodes,
                                            no_interactions,
                                            no_thresholds,
                                            index,
                                            proposal_threshold_sd,
                                            proposal_interaction_sd,
                                            iteration + 1,
                                            m,
                                            cauchy_scale,
                                            threshold_alpha,
                                            threshold_beta,
                                            interactions,
                                            thresholds,
                                            parallel);

    NumericMatrix interactions = out["interactions"];
    NumericMatrix thresholds = out["thresholds"];
    NumericMatrix proposal_threshold_sd = out["proposal_threshold_sd"];
    NumericMatrix proposal_interaction_sd = out["proposal_interaction_sd"];


    //Output -------------------------------------------------------------------
    if(save == TRUE) {
      //Save raw samples -------------------------------------------------------
      cntr = 0;
      for(int node1 = 0; node1 < no_nodes - 1; node1++) {
        for(int node2 = node1 + 1; node2 < no_nodes;node2++) {
          out_interactions(iteration, cntr) = interactions(node1, node2);
          cntr++;
        }
      }
      cntr = 0;
      for(int node = 0; node < no_nodes; node++) {
        for(int category = 0; category < no_categories[node]; category++) {
          out_thresholds(iteration, cntr) = thresholds(node, category);
          cntr++;
        }
      }
    } else {
      //Compute running averages -----------------------------------------------
      for(int node1 = 0; node1 < no_nodes - 1; node1++) {
        for(int node2 = node1 + 1; node2 < no_nodes; node2++) {
          out_interactions(node1, node2) *= iteration;
          out_interactions(node1, node2) += interactions(node1, node2);
          out_interactions(node1, node2) /= iteration + 1;
          out_interactions(node2, node1) = out_interactions(node1, node2);
        }

        for(int category = 0; category < no_categories[node1]; category++) {
          out_thresholds(node1, category) *= iteration;
          out_thresholds(node1, category) += thresholds(node1, category);
          out_thresholds(node1, category) /= iteration + 1;
        }
      }
      for(int category = 0; category < no_categories[no_nodes - 1]; category++) {
        out_thresholds(no_nodes - 1, category) *= iteration;
        out_thresholds(no_nodes - 1, category) +=
          thresholds(no_nodes - 1, category);
        out_thresholds(no_nodes - 1, category) /= iteration + 1;
      }
    }
  }

  return List::create(Named("interactions") = out_interactions,
                      Named("thresholds") = out_thresholds);
}
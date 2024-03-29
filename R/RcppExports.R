# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

data_gibbs_person_parallel <- function(interactions, thresholds, data, no_categories, no_persons, iter, no_nodes, max_no_categories, augmented_data) {
    invisible(.Call(`_dmhBGM_data_gibbs_person_parallel`, interactions, thresholds, data, no_categories, no_persons, iter, no_nodes, max_no_categories, augmented_data))
}

dmh_gibbs_sampler <- function(observations, O_thresholds, O_interactions, m, no_categories, Index, proposal_threshold_sd, proposal_interaction_sd, cauchy_scale, threshold_alpha, threshold_beta, edge_prior, theta, beta_bernoulli_alpha, beta_bernoulli_beta, gamma, interactions, thresholds, iter, burnin, save = FALSE, display_progress = FALSE, parallel = FALSE) {
    .Call(`_dmhBGM_dmh_gibbs_sampler`, observations, O_thresholds, O_interactions, m, no_categories, Index, proposal_threshold_sd, proposal_interaction_sd, cauchy_scale, threshold_alpha, threshold_beta, edge_prior, theta, beta_bernoulli_alpha, beta_bernoulli_beta, gamma, interactions, thresholds, iter, burnin, save, display_progress, parallel)
}

dmh_gibbs_sampler_estimation <- function(observations, O_thresholds, O_interactions, m, no_categories, Index, proposal_threshold_sd, proposal_interaction_sd, cauchy_scale, threshold_alpha, threshold_beta, interactions, thresholds, iter, burnin, save = FALSE, display_progress = FALSE, parallel = FALSE) {
    .Call(`_dmhBGM_dmh_gibbs_sampler_estimation`, observations, O_thresholds, O_interactions, m, no_categories, Index, proposal_threshold_sd, proposal_interaction_sd, cauchy_scale, threshold_alpha, threshold_beta, interactions, thresholds, iter, burnin, save, display_progress, parallel)
}


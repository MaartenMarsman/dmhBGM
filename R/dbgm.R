#' Bayesian structure learning in Markov Random Fields of mixed binary and
#' ordinal variables using MCMC.
#'
#' The function \code{dbgm} explores the joint posterior distribution of
#' structures and parameters in a Markov Random Field for mixed binary and
#' ordinal variables.
#'
#' A discrete spike and slab prior distribution is stipulated on the pairwise
#' interactions. By formulating it as a mixture of mutually singular
#' distributions, the function can use a combination of Metropolis-Hastings and
#' Gibbs sampling to create a Markov chain that has the joint posterior
#' distribution as invariant. This package uses Double Metropolis Hastings
#' \insertCite{Liang_2010}{dmhBGM} to circumvent having to compute the
#' intractable normalizing constant of the Markov Random Field, and uses an
#' adaptive Metropolis to "learn" an optimal proposal distribution: Adjusting
#' the proposal variance to match the acceptance probability of the
#' random walk Metropolis algorithm to be close to the optimum of \code{.234}
#' using a Robbins-Monro type algorithm.
#'
#' The slab distribution is a Cauchy with an optional scaling parameter. A
#' Beta-prime distribution is used for the exponent of the category parameters.
#' Two prior distributions are implemented for edge inclusion variables (i.e.,
#' the prior probability that an edge is included); the Bernoulli prior and the
#' Beta-Bernoulli prior.
#'
#' @param x A data frame or matrix with \code{n} rows and \code{p} columns
#' containing binary and ordinal variables for \code{n} independent observations
#' and \code{p} variables in the network. Variables are recoded as non-negative
#' integers \code{(0, 1, ..., m)} if not already done. Unobserved categories are
#' collapsed into other categories after recoding (i.e., if category 1 is
#' unobserved, the data will be recoded from (0, 2) to (0, 1)).
#' @param iter The number of iterations of the Gibbs sampler. The default of
#' \code{1e4} is for illustrative purposes. For stable estimates, it is
#' recommended to run the Gibbs sampler for at least \code{1e5} iterations.
#' @param burnin The number of iterations of the Gibbs sampler before its output
#' is saved. Since it may take some time for the Gibbs sampler to converge to
#' the posterior distribution, it is recommended not to set this number too low.
#' @param dmhsamples The DMH approach generates a new, augmented data set to
#' update model parameters with Metropolis Hastings. It must generate a new
#' augmented data set for each parameter of the model in each iteration of the
#' Gibbs sampler. Since we cannot sample these data directly from the MRF, we
#' must use a Gibbs sampler to generate them.  The dmhsamples argument controls how
#' many iterations are used for this "inner" Gibbs sampler. The default of
#' \code{1} is for illustrative purposes. For stable estimates,
#' \insertCite{ParkEtAl_2020}{dmhBGM} recommend running the DMH-Gibbs sampler
#' for at least \code{n * p} iterations.
#' @param cauchy_scale The scale of the Cauchy prior for interactions. Defaults
#' to \code{2.5}.
#' @param edge_prior The prior distribution for the edges or structure of the
#' network. Two prior distributions are currently implemented: The Bernoulli
#' model \code{edge_prior = "Bernoulli"} assumes that the probability that an
#' edge between two variables is included is equal to
#' \code{inclusion_probability} and independent of other edges or variables.
#' When \code{inclusion_probability = 0.5}, this implies that each network
#' structure receives the same prior weight. The Beta-Bernoulli model
#' \code{edge_prior = "Beta-Bernoulli"} assumes a beta prior for the unknown
#' inclusion probability with shape parameters \code{beta_bernoulli_alpha} and
#' \code{beta_bernoulli_beta}. If \code{beta_bernoulli_alpha = 1} and
#' \code{beta_bernoulli_beta = 1}, this means that networks with the same
#' complexity (number of edges) receive the same prior weight. Defaults to
#' \code{edge_prior = "Bernoulli"}.
#' @param inclusion_probability The prior edge inclusion probability for the
#' Bernoulli model. Can be a single probability, or a matrix of \code{p} rows
#' and \code{p} columns specifying an inclusion probability for each edge pair.
#' Defaults to \code{inclusion_probability = 0.5}.
#' @param beta_bernoulli_alpha,beta_bernoulli_beta The two shape parameters of
#' the Beta prior density for the Bernoulli inclusion probability. Must be
#' positive numbers. Defaults to \code{beta_bernoulli_alpha = 1} and
#' \code{beta_bernoulli_beta = 1}.
#' @param threshold_alpha,threshold_beta The shape parameters of the beta-prime
#' prior density for the threshold parameters. Must be positive values. If the
#' two values are equal, the prior density is symmetric about zero. If
#' \code{threshold_beta} is greater than \code{threshold_alpha}, the
#' distribution is skewed to the left, and if \code{threshold_beta} is less than
#' \code{threshold_alpha}, it is skewed to the right. Smaller values tend to
#' lead to more diffuse prior distributions.
#' @param save Should the function collect and return all samples from the Gibbs
#' sampler (\code{save = TRUE})? Or should it only return the (model-averaged)
#' posterior means (\code{save = FALSE})? Defaults to \code{FALSE}.
#' @param display_progress Should the function show a progress bar
#' (\code{display_progress = TRUE})? Or not (\code{display_progress = FALSE})?
#' Defaults to \code{TRUE}.
#' @param parallel Should the DMH iterations be run in parallel? Defaults to \code{FALSE}.
#' Defaults to \code{TRUE}.
#' @param no_cores If the iterations for DMH are run in parallel, how many cores should be used?
#' Defaults to \link{\code{RcppParallel::defaultNumThreads() - 1}}.

#'
#' @return If \code{save = FALSE} (the default), the result is a list of class
#' ``bgms'' containing the following matrices:
#' \itemize{
#' \item \code{gamma}: A matrix with \code{p} rows and \code{p} columns,
#' containing posterior inclusion probabilities of individual edges.
#' \item \code{interactions}: A matrix with \code{p} rows and \code{p} columns,
#' containing model-averaged posterior means of the pairwise associations.
#' \item \code{thresholds}: A matrix with \code{p} rows and \code{max(m)}
#' columns, containing model-averaged category thresholds.
#' }
#'
#' If \code{save = TRUE}, the result is a list of class ``dmhBGM'' containing:
#' \itemize{
#' \item \code{gamma}: A matrix with \code{iter} rows and
#' \code{p * (p - 1) / 2} columns, containing the edge inclusion indicators from
#' every iteration of the Gibbs sampler.
#' \item \code{interactions}: A matrix with \code{iter} rows and
#' \code{p * (p - 1) / 2} columns, containing parameter states from every
#' iteration of the Gibbs sampler for the pairwise associations.
#' \item \code{thresholds}: A matrix with \code{iter} rows and
#' \code{sum(m)} columns, containing parameter states from every iteration of
#' the Gibbs sampler for the category thresholds.
#' }
#' Column averages of these matrices provide the model-averaged posterior means.
#'
#' @references
#' \insertAllCited{}
#'
#' @importFrom RcppParallel RcppParallelLibs
#' @export
dmhbgm = function(x,
               iter = 1e4,
               burnin = 1e3,
               dmhsamples = 1,
               cauchy_scale = 2.5,
               edge_prior = c("Bernoulli", "Beta-Bernoulli"),
               inclusion_probability = 0.5,
               beta_bernoulli_alpha = 1,
               beta_bernoulli_beta = 1,
               threshold_alpha = 0.5,
               threshold_beta = 0.5,
               save = FALSE,
               display_progress = TRUE,
               parallel = FALSE,
               no_cores = RcppParallel::defaultNumThreads() - 1) {

  #Check data input ------------------------------------------------------------
  if(!inherits(x, what = "matrix") && !inherits(x, what = "data.frame"))
    stop("The input x needs to be a matrix or dataframe.")
  if(inherits(x, what = "data.frame"))
    x = data.matrix(x)
  if(ncol(x) < 2)
    stop("The matrix x should have more than one variable (columns).")
  if(nrow(x) < 2)
    stop("The matrix x should have more than one observation (rows).")

  #Check Gibbs input -----------------------------------------------------------
  if(abs(iter - round(iter)) > sqrt(.Machine$double.eps))
    stop("Parameter ``iter'' needs to be a positive integer.")
  if(iter <= 0)
    stop("Parameter ``iter'' needs to be a positive integer.")
  if(abs(burnin - round(burnin)) > sqrt(.Machine$double.eps) || burnin < 0)
    stop("Parameter ``burnin'' needs to be a non-negative integer.")
  if(burnin <= 0)
    stop("Parameter ``burnin'' needs to be a positive integer.")

  #Check DMH input -----------------------------------------------------------
  if(abs(dmhsamples - round(dmhsamples)) > sqrt(.Machine$double.eps))
    stop("Parameter ``dmhsamples'' needs to be a positive integer.")
  if(dmhsamples <= 0)
    stop("Parameter ``dmhsamples'' needs to be a positive integer.")

  #Check prior set-up for the interaction parameters ---------------------------
  interaction_prior = "Cauchy"
  if(cauchy_scale <= 0 || is.na(cauchy_scale) || is.infinite(cauchy_scale))
    stop("The scale of the Cauchy prior needs to be positive.")

  #Check prior set-up for the edge indicators ----------------------------------
  edge_prior = match.arg(edge_prior)
  if(edge_prior == "Bernoulli") {
    if(length(inclusion_probability) == 1) {
      theta = inclusion_probability[1]
      if(is.na(theta) || is.null(theta))
        stop("There is no value specified for the inclusion probability.")
      if(theta <= 0)
        stop("The inclusion probability needs to be positive.")
      if(theta >= 1)
        stop("The inclusion probability cannot exceed the value one.")
      theta = matrix(theta, nrow = ncol(x), ncol = ncol(x))
    } else {
      if(!inherits(inclusion_probability, what = "matrix") &&
         !inherits(inclusion_probability, what = "data.frame"))
        stop("The input for the inclusion probability argument needs to be a single number, matrix, or dataframe.")

      if(inherits(inclusion_probability, what = "data.frame")) {
        theta = data.matrix(inclusion_probability)
      } else {
        theta = inclusion_probability
      }
      if(!isSymmetric(theta))
        stop("The inclusion probability matrix needs to be symmetric.")
      if(ncol(theta) != ncol(x))
        stop("The inclusion probability matrix needs to have as many rows (columns) as there are variables in the data.")

      if(any(is.na(theta[lower.tri(theta)])) ||
         any(is.null(theta[lower.tri(theta)])))
        stop("One or more elements of the elements in inclusion probability matrix are not specified.")
      if(any(theta[lower.tri(theta)] <= 0))
        stop(paste0("The inclusion probability matrix contains negative or zero values;\n",
                    "inclusion probabilities need to be positive."))
      if(any(theta[lower.tri(theta)] >= 1))
        stop(paste0("The inclusion probability matrix contains values greater than or equal to one;\n",
                    "inclusion probabilities cannot exceed or equal the value one."))
    }
  }
  if(edge_prior == "Beta-Bernoulli") {
    theta = matrix(0.5, nrow = ncol(x), ncol = ncol(x))
    if(beta_bernoulli_alpha <= 0 || beta_bernoulli_beta <= 0)
      stop("The scale parameters of the beta distribution need to be positive.")
    if(!is.finite(beta_bernoulli_alpha) || !is.finite(beta_bernoulli_beta))
      stop("The scale parameters of the beta distribution need to be finite.")
    if(is.na(beta_bernoulli_alpha) || is.na(beta_bernoulli_beta) ||
       is.null(beta_bernoulli_alpha) || is.null(beta_bernoulli_beta))
      stop("Values for both scale parameters of the beta distribution need to be specified.")
  }

  #Check prior set-up for the threshold parameters -----------------------------
  if(threshold_alpha <= 0  | !is.finite(threshold_alpha))
    stop("Parameter ``threshold_alpha'' needs to be positive.")
  if(threshold_beta <= 0  | !is.finite(threshold_beta))
    stop("Parameter ``threshold_beta'' needs to be positive.")

  #Check parallel arguments
  if(!isTRUE(parallel) || isFALSE(parallel))
    stop("Parameter ``parallel'' must be TRUE or FALSE.")
  if(isTRUE(parallel)) {
    if(no_cores <= 0)
      stop("Parameter ``no_cores'' must be larger than 0.")

    # ensure we nicely clean up any value set for RcppParallel::setThreadOptions
    old_no_cores <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA)
    if (!is.na(old_no_cores))
      on.exit(RcppParallel::setThreadOptions(numThreads = old_no_cores))
    RcppParallel::setThreadOptions(numThreads = no_cores)
  }

  #Check na.action -------------------------------------------------------------
  na.action = "listwise"

  #Format the data input -------------------------------------------------------
  data = reformat_data(x = x, na.action)
  x = data$x
  no_categories = data$no_categories
  missing_index = data$missing_index
  na.impute = data$na.impute

  no_nodes = ncol(x)
  no_interactions = no_nodes * (no_nodes - 1) / 2
  no_thresholds = sum(no_categories)

  #Proposal set-up for the interaction parameters ------------------------------

  #Set up the variance of the (normal) proposal distribution
  proposal_threshold_sd = matrix(1, nrow = no_nodes, ncol = max(no_categories))
  proposal_interaction_sd = matrix(1, nrow = no_nodes, ncol = no_nodes)

  # Starting value of model matrix:
  gamma = matrix(1,
                 nrow = no_nodes,
                 ncol = no_nodes)

  #Starting values of interactions and thresholds (posterior mode)
  interactions = matrix(0, nrow = no_nodes, ncol = no_nodes)
  thresholds = matrix(0, nrow = no_nodes, ncol = max(no_categories))

  #Precomputing number of observations per category for each node.
  O_thresholds = matrix(0,
                        nrow = max(no_categories) + 1,
                        ncol = no_nodes)
  for(node in 1:no_nodes) {
    for(category in 0:no_categories[node]) {
      O_thresholds[category + 1, node] = sum(x[, node] == category)
    }
  }
  O_interactions = t(x) %*% x

  # Index vector used to sample interactions in a random order.
  Index = matrix(0,
                 nrow = no_nodes * (no_nodes - 1) / 2,
                 ncol = 3)
  cntr = 0
  for(node1 in 1:(no_nodes - 1)) {
    for(node2 in (node1 + 1):no_nodes) {
      cntr =  cntr + 1
      Index[cntr, 1] = cntr
      Index[cntr, 2] = node1
      Index[cntr, 3] = node2
    }
  }

  #The Metropolis within Gibbs sampler -----------------------------------------
  out = dmh_gibbs_sampler(observations = x,
                          O_thresholds,
                          O_interactions,
                          dmhsamples,
                          no_categories  = no_categories,
                          Index = Index,
                          proposal_threshold_sd = proposal_threshold_sd,
                          proposal_interaction_sd = proposal_interaction_sd,
                          cauchy_scale = cauchy_scale,
                          threshold_alpha = threshold_alpha,
                          threshold_beta = threshold_beta,
                          edge_prior = edge_prior,
                          theta = theta,
                          beta_bernoulli_alpha = beta_bernoulli_alpha,
                          beta_bernoulli_beta = beta_bernoulli_beta,
                          gamma = gamma,
                          interactions = interactions,
                          thresholds = thresholds,
                          iter = iter,
                          burnin = burnin,
                          save = save,
                          display_progress = display_progress,
                          parallel = parallel)

  #Preparing the output --------------------------------------------------------
  if(save == FALSE) {
    gamma = out$gamma
    interactions = out$interactions
    tresholds = out$thresholds

    if(is.null(colnames(x))){
      data_columnnames = paste0("node ", 1:no_nodes)
      colnames(interactions) = data_columnnames
      rownames(interactions) = data_columnnames
      colnames(gamma) = data_columnnames
      rownames(gamma) = data_columnnames
      rownames(thresholds) = data_columnnames
    } else {
      data_columnnames <- colnames(x)
      colnames(interactions) = data_columnnames
      rownames(interactions) = data_columnnames
      colnames(gamma) = data_columnnames
      rownames(gamma) = data_columnnames
      rownames(thresholds) = data_columnnames
    }

    colnames(tresholds) = paste0("category ", 1:max(no_categories))

    output = list(gamma = gamma,
                  interactions = interactions,
                  thresholds = thresholds,
                  edge_prior = edge_prior,
                  inclusion_probability = inclusion_probability,
                  beta_bernoulli_alpha = beta_bernoulli_alpha,
                  beta_bernoulli_beta = beta_bernoulli_beta,
                  save = save,
                  colnames = data_columnnames)
    class(output) = "dmhBGM"
    return(output)
  } else {
    gamma = out$gamma
    interactions = out$interactions
    thresholds = out$thresholds

    if(is.null(colnames(x))){
      data_columnnames <- 1:ncol(x)
    } else {
      data_columnnames <- colnames(x)
    }
    p <- ncol(x)
    names_bycol <- matrix(rep(data_columnnames, each = p), ncol = p)
    names_byrow <- matrix(rep(data_columnnames, each = p), ncol = p, byrow = T)
    names_comb <- matrix(paste0(names_byrow, "-", names_bycol), ncol = p)
    names_vec <- names_comb[lower.tri(names_comb)]

    colnames(gamma) = colnames(interactions) = names_vec

    names = character(length = sum(no_categories))
    cntr = 0
    for(node in 1:no_nodes) {
      for(category in 1:no_categories[node]) {
        cntr = cntr + 1
        names[cntr] = paste0("threshold(",node, ", ",category,")")
      }
    }
    colnames(thresholds) = names

    rownames(gamma) = paste0("Iter. ", 1:iter)
    rownames(interactions) = paste0("Iter. ", 1:iter)
    rownames(thresholds) = paste0("Iter. ", 1:iter)

    output = list(gamma = gamma,
                  interactions = interactions,
                  thresholds = thresholds,
                  edge_prior = edge_prior,
                  inclusion_probability = inclusion_probability,
                  beta_bernoulli_alpha = beta_bernoulli_alpha,
                  beta_bernoulli_beta = beta_bernoulli_beta,
                  save = save,
                  colnames = data_columnnames)
    class(output) = "dmhBGM"
    return(output)
  }
}


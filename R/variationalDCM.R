#' @title Variational Bayesian estimation for DCMs
#' @description \code{variationalDCM()} fits DCMs by VB algorithms
#'
#' @section {DINA model}:
#' The DINA model has two types of model parameters: slip \eqn{s_j} and guessing \eqn{g_j} for \eqn{j=1,\cdots,J}.
#' The ideal response of the DINA model is given as
#' \deqn{\eta_{lj}=\prod_k\alpha_{lk}^{q_{jk}}.}
#' This ideal response equals 1 if the attribute mastery pattern \eqn{\boldsymbol{\alpha}_l} of the \eqn{l}-th class satisfies \eqn{\alpha_{lk} \geq q_{jk}} for all \eqn{k} but 0 otherwise. Using a class indicator vector \eqn{\boldsymbol{z}_i}, the ideal response for respndent \eqn{i} is written as
#' \deqn{\eta_{ij} = \prod_l\eta_{lj}z_{il}.}
#' Item response function of the DINA model can be written as
#' \deqn{P(X_{ij}=1|\eta_{ij},s_j,g_j)=(1-s_j)^{\eta_{ij}}g_j^{1-\eta_{ij}}.}
#' The guessing parameter represents the probability that a respondent obtains a correct response to an item when the ideal response is zero; the slip parameter represents the probability of a respondent obtaining an incorrect response when the ideal response is one; that is,
#' \deqn{P(X_{ij}=0|\eta_{ij}=1)=s_j,\\ P(X_{ij}=1|\eta_{ij}=0)=g_j.}
#'
#' @section {DINO model}:
#' The DINO model can be represented by slightly modifying the DINA model. The ideal response of the DINO model is given as
#' \deqn{\eta_{lj}=1-\prod_k(1-\alpha_{lk})^{q_{jk}}.}
#' The ideal response takes the value of 1 when a respondent masters at least one of the relevant attributes with an item; otherwise, it is 0.
#'
#' @section {MC-DINA model}:
#' The MC-DINA model is a extension of the DINA model to capture nominal item responses. Assume each item has \eqn{C_j} options, and we represent an item response as \eqn{X_{ijc}\in\{0,1\}}.
#'
#' @section {Saturated DCM}:
#' The saturated DCM is a generalized model such as the G-DINA and GDM.
#'
#' @param X  N by J item response data for the DINA, DINO, MC-DINA, and
#'   saturated DCM models. Alternatively, T-length list or 3-dim array whose
#'   elements are N by J/T binary item response data matrices for the HM-DCM
#' @param Q  J by K binary Q-matrix for the DINA, DINO, and saturated DCM models. For the MC-DINA model, its size should be J by (K+2). Alternatively, T-length list or 3-dim array whose elements
#'   are J/T by K Q-matrices for the HM-DCM
#' @param model specify one of "dina", "dino", "mc_dina", "satu_dcm", and
#'   "hm_dcm"
#' @param max_it Maximum number of iterations (default: 500)
#' @param epsilon convergence tolerance for iterations (default: 1e-4)
#' @param verbose logical, controls whether to print progress (default: TRUE)
#' @param ... additional arguments such as hyperparameter values
#'
#' @return \code{variationalDCM} returns an object of class
#'   \code{variationalDCM}. We provide \code{summary} function to summarize a
#'   result and users can check the following information:
#' \describe{
#'   \item{Attribute.mastery.patterns}{scored attribute mastery patterns}
#'   \item{Estimated.model.parameters}{estimates and posterior standard deviations of model parameters}
#'   \item{ELBO}{resulting value of evidence lower bound}
#'   \item{time}{time spent in computation}
#' }
#'
#' @references Yamaguchi, K., & Okada, K. (2020). Variational Bayes inference
#'   for the DINA model. \emph{Journal of Educational and Behavioral
#'   Statistics}, 45(5), 569-597. \doi{10.3102/1076998620911934}
#'
#'   Yamaguchi, K. (2020). Variational Bayesian inference for the
#'   multiple-choice DINA model. \emph{Behaviormetrika}, 47(1), 159-187.
#'   \doi{10.1007/s41237-020-00104-w}
#'
#'   Yamaguchi, K., Okada, K. (2020). Variational Bayes Inference Algorithm for
#'   the Saturated Diagnostic Classification Model. \emph{Psychometrika}, 85(4),
#'   973–995. \doi{10.1007/s11336-020-09739-w}
#'
#'   Yamaguchi, K., & Martinez, A. J. (2023). Variational Bayes inference for
#'   hidden Markov diagnostic classification models. \emph{British Journal of
#'   Mathematical and Statistical Psychology}, 00, 1– 25.
#'   \doi{10.1111/bmsp.12308}
#'
#' @examples
#' # fit the DINA model
#' Q = sim_Q_J80K5
#' sim_data = dina_data_gen(Q=Q,I=200)
#' res = variationalDCM(X=sim_data$X, Q=Q, model="dina")
#' summary(res)
#'
#' @export

variationalDCM = function(
    X,
    Q,
    model,
    max_it  = 500,
    epsilon = 1e-04,
    verbose = TRUE,...
){

  t1 = Sys.time()
  variationalDCMcall = match.call()

  if(model == "dina"){
    res = dina(X=X, Q=Q, max_it=max_it, epsilon=epsilon,...)
  }
  else if(model == "dino"){
    res = dino(X=X, Q=Q, max_it=max_it, epsilon=epsilon,...)
  }
  else if(model == "satu_dcm"){
    res = satu_dcm(X=X, Q=Q, max_it=max_it, epsilon=epsilon,...)
  }
  else if(model == "mc_dina"){
    res = mc_dina(X=X, Q=Q, max_it=max_it, epsilon=epsilon,...)
  }
  else if(model == "hm_dcm"){
    res = hm_dcm(X=X, Q=Q, max_it=max_it, epsilon=epsilon,...)
  }

  t2 = Sys.time()
  res$time = t2-t1
  res$call = variationalDCMcall
  class(res) = "variationalDCM"
  return(res)
}

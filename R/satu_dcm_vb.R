#' for the saturated DCM.
#'
#' \code{satu_dcm()} returns variational Bayesian estimates for the saturated DCM.
#'
#' @param X I by J binary matrix, item response data
#' @param Q J by K binary matrix, Q-matrix
#' @param max_it The maximum number of iterations (default: 500)
#' @param epsilon The convergence tolerance for iterations (default: 1e-5)
#' @param seed The seed value (default: 123)
#' @param verbose Logical, controls whether to print progress (default: TRUE)
#'
#' @return A list including:
#' \describe{
#'   \item{theta_est}{the estimates of the conditional response probabilities}
#'   \item{theta_sd}{the standard deviations of theta}
#'   \item{pi_est}{the estimates of categorical distribution's parameters, which is the prior for the attribute mastery indicator}
#'   \item{pi_sd}{the standard deviations of pi}
#'   \item{r_il}{the expectation of the attribute mastery pattern conditioned by the categorical distribution}
#'   \item{a_ast}{the estimated variational Bayesian posterior for theta}
#'   \item{delta_ast}{the estimated variational Bayesian posterior for pi}
#'   \item{A_0}{the estimates of beta distribution's parameter, which is the prior for theta}
#'   \item{B_0}{the estimates of beta distribution's parameter, which is the prior for theta}
#'   \item{delta_0}{the estimates of Dirichlet distribution's parameter, which is the prior for pi}
#'   \item{l_lb}{the computed lower bound of each iteration}
#'   \item{att_pat_est}{the estimated attribute mastery patterns}
#   \item{A}{all of the possible attribute mastery patterns}
#   \item{Q}{the entered Q-matrix}
#   \item{X}{the entered data matrix}
#'   \item{G_j}{the computed G-matrix}
#'   \item{m}{the number of the performed iterations}
#'   \item{seed}{the entered seed number}
#' }
#'
#' @references Yamaguchi, K., Okada, K. (2020). Variational Bayes Inference
#'   Algorithm for the Saturated Diagnostic Classification Model.
#'   \emph{Psychometrika}, 85(4), 973–995.
#'   \url{https://doi.org/10.1007/s11336-020-09739-w}.
#'
#' @export

satu_dcm = function(X,
                    Q,
                    max_it  = 500,
                    epsilon = 1e-05,
                    seed = 123,
                    verbose = TRUE
){
  if(!inherits(X, "matrix")){
    X <- as.matrix(X)
  }
  if(!inherits(Q, "matrix")){
    Q <- as.matrix(Q)
  }

  set.seed(seed);  I <- nrow(X);  J <- nrow(Q);  K <- ncol(Q);  L <- 2^K
  not_zero_q <- apply(Q, 1, function(x) which(x != 0))

  # Attribute pattern matrix

  A <- as.matrix(expand.grid(lapply(1:K, function(x)rep(0:1))))
  A_j <- lapply(not_zero_q, function(x) A[,x,drop=F])
  A_red <- lapply(A_j , function(x) apply(x,1,function(y) paste0(y,collapse = "")))
  A_red_uni <- lapply(A_j , function(x) unique(apply(x,1,function(y) paste0(y,collapse = ""))))

  # Make G-matrix

  G_j <- lapply(1:J, function(j) t(sapply(A_red_uni[[j]], function(x) x == A_red[[j]] ))*1 )
  att_pat <- apply(A, 1, function(x) paste0(x, collapse=""))
  for(j in 1:J) {
    colnames(G_j[[j]]) <- att_pat
  }

  # Hyper parameters

  delta_0 = rep(1,L)# For π

  # Weekly informative prior

  number_of_attributes <- lapply(A_red_uni, function(x) sapply(strsplit(x, ""), function(y)sum(as.numeric(y))) )
  A_0_hyperparam <- seq(from = 1+epsilon, to = 2, length.out = max(unlist( number_of_attributes))+1)
  B_0_hyperparam <- seq(from = 2, to = 1+epsilon, length.out = max(unlist( number_of_attributes))+1)
  A_0 <- lapply( number_of_attributes , function(x){A_0_hyperparam[x + 1] })
  B_0 <- lapply( number_of_attributes , function(x){B_0_hyperparam[x + 1] })

  # Initialization

  r_il <- matrix(1/L ,ncol=L, nrow = I)

  one_vec  = matrix(1,nrow=I,ncol=1)

  # lower bound of log marginal likelihood

  llb_fun <- function(X,G_j,delta_ast,delta_0, A_ast, A_0, B_ast, B_0, r_il){
    tmp1 <- 0
    for(j in 1:length(G_j)){
      tmp1 <- tmp1+ sum(((X[,j]%*%(digamma(A_ast[[j]]) - digamma(A_ast[[j]]+B_ast[[j]])) + (1-X[,j])%*%(digamma(B_ast[[j]]) - digamma(A_ast[[j]]+B_ast[[j]]))) %*% G_j[[j]]) * r_il)
    }
    tmp2 <- sum(r_il * ( one_vec%*%(digamma(delta_ast) - digamma(sum(delta_ast))) - log(r_il)))
    tmp3 <- sum(lgamma(delta_ast))- lgamma(sum(delta_ast)) - (sum(lgamma(delta_0))- lgamma(sum(delta_0))) + sum((delta_0 - delta_ast)*(digamma(delta_ast) - digamma(sum(delta_ast))))

    A_ast_unlist <- unlist(A_ast)
    B_ast_unlist <- unlist(B_ast)
    A_0_unlist <- unlist(A_0)
    B_0_unlist <- unlist(B_0)
    tmp4 <- sum(lbeta(a =A_ast_unlist,  b=B_ast_unlist) - lbeta(a =A_0_unlist,  b=B_0_unlist) +  (A_0_unlist - A_ast_unlist)*(digamma(A_ast_unlist) - digamma(A_ast_unlist + B_ast_unlist ) ) + (B_0_unlist - B_ast_unlist)*(digamma(B_ast_unlist) - digamma(A_ast_unlist + B_ast_unlist ) ) )

    tmp1 + tmp2 + tmp3 + tmp4  }
  l_lb = rep(NA_real_, max_it+1)
  l_lb[1] = 100

  m = 1
  for(m in 1:max_it){

    # VM-step
    delta_ast <- colSums(r_il) + delta_0
    A_ast  <- lapply(1:J, function(j) t(G_j[[j]] %*% t(r_il) %*% X[,j]     + A_0[[j]] ))
    B_ast  <- lapply(1:J, function(j) t(G_j[[j]] %*% t(r_il) %*% (1-X[,j]) + B_0[[j]] ))

    E_log_pi      = digamma(delta_ast) - digamma(sum(delta_ast))
    E_log_theta   = lapply(1:J, function(j) digamma(A_ast[[j]]) - digamma(A_ast[[j]] + B_ast[[j]]))
    E_log_1_theta = lapply(1:J, function(j) digamma(B_ast[[j]]) - digamma(A_ast[[j]] + B_ast[[j]]))

    # VE-step: r_il
    temp <- matrix(0, ncol=I, nrow = L)
    for(j in 1:J) temp <- temp + t(G_j[[j]] ) %*% (t(E_log_theta[[j]]) %*% t(X[,j,drop=F]) +  t(E_log_1_theta[[j]]) %*% t(1 - X[,j,drop=F]))
    log_rho_il <- t(temp + E_log_pi)
    temp <- exp(log_rho_il)
    r_il =  temp / rowSums(temp)

    l_lb[m+1] <- llb_fun(X,G_j,delta_ast,delta_0, A_ast, A_0, B_ast,B_0, r_il)

    if(verbose){
      cat("\riteration = ", m+1, sprintf(": l_lb = %.05f", l_lb[m+1]), sprintf(": last change = %.05f", abs(l_lb[m] - l_lb[m+1])))
    }

    if(abs(l_lb[m] - l_lb[m+1]) < epsilon){
      if(verbose){
        cat("\nreached convergence.")
      }
      break()
    }
  }
  l_lb <- l_lb[-1]

  delta_sum <- sum(delta_ast)
  pi_est <-  delta_ast/delta_sum
  pi_sd <-sqrt(delta_ast*(delta_sum - delta_ast)/(delta_sum^2*(delta_sum+1)) )
  theta_est <- mapply(function(x,y) x/(x+y), A_ast, B_ast)
  theta_sd <- mapply(function(x,y) sqrt((x*y)/(((x+y)^2) *(x+y+1)) ),  A_ast, B_ast)

  list(theta_est = theta_est,
       theta_sd = theta_sd,
       pi_est = pi_est,
       pi_sd = pi_sd,
       r_il  = r_il,
       A_ast = A_ast,
       delta_ast   = delta_ast,
       A_0 = A_0,
       B_0 = B_0,
       delta_0 = delta_0,
       l_lb = l_lb,
       att_pat_est = A[apply(r_il, 1, which.max),],
       #A = A,
       #Q = Q,
       #X = X,
       G_j = G_j,
       m = m,
       seed = seed)
}

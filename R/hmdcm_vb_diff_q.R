#' for the hidden Markov DCM.
#'
#' \code{hmdcm_diff_q()} returns variational Bayesian estimates for the hidden
#' Markov DCM.
#'
#' @param X I by J by T binary 3-dimension array, item response data
#' @param Q J by k by T binary 3-dimension array, Q-matrix
#' @param max_it Maximum number of iterations (default: 500)
#' @param epsilon convergence tolerance for iterations (default: 10E-6)
#' @param verbose logical (default: TRUE)
#' @param A_0 the value of hyperparameter A_0
#' @param B_0 the value of hyperparameter B_0
#' @param delta_0 the value of hyperparameter delta_0
#' @param ommega_0 the value of hyperparameter ommega_0
#' @param Test_versions test version
#' @param test_order test order
#' @param model "general" or "DINA" (default: "gengeral")
#' @param random_start Logical, controls whether to print progress (default: "FALSE")
#'
#' @return A list including:
#' \describe{
#'   \item{theta_est}{the estimated conditional probability of taking a wrong response}
#'   \item{theta_sd}{the estimated standard deviation of theta}
#'   \item{pi_est}{the estimates of variational parameter for class indicator}
#'   \item{pi_sd}{the estimatd standard deviation of pi}
#'   \item{Tau_est}{the estimated transition probability}
#'   \item{Tau_sd}{the estimated standard deviation of the transition probability}
#'   \item{post_max_class}{the result of class analysis}
#'   \item{MAP_att_pat}{the MAP estimates of attribute mastery patterns}
#'   \item{att_master_prob}{the estimated attribute mastery probabilities}
#'   \item{EAP_att_pat}{the EAP estimates of attribute mastery patterns}
#'   \item{A_ast}{the estimated A_ast}
#'   \item{delta_ast}{the estimated delta_ast}
#'   \item{ommega_ast}{the estimated ommega_ast}
#'   \item{E_z_itl}{the resulted expectations of  attribute mastery pattern}
#'   \item{E_z_itl_z_itm1l}{the resulted expectations of  attribute mastery pattern at the time point t and t-1}
#'   \item{A_0}{the entered value of A_0}
#'   \item{B_0}{the entered value of B_0}
#'   \item{delta_0}{the entered value of delta_0}
#'   \item{ommega_0}{the entered value of ommega_0}
#'   \item{l_lb}{the computed lower bound of each iteration}
#'   \item{gamma_t_x_it}{the computed value of the normalizing constant for calculating zeta}
#'   \item{log_zeta_sum}{the computed value of the normalizing constant for calculating variational posterior for class indicator}
#'   \item{A}{all the possible attribute mastery patterns}
#'   \item{Q}{the entered Q-matrix}
#'   \item{X}{the entered data matrix}
#'   \item{G_jt}{the computed G-matrix}
#'   \item{m}{the number of performed iterations}
#' }
#'
#' @references Yamaguchi, K., & Martinez, A. J. (2021). Variational Bayesian
#'   Inference Posterior Approximation Algorithm for Hidden Markov Diagnostic
#'   Classification Models. \url{https://doi.org/10.31234/osf.io/28jkf}.
#'
#' @export

hmdcm_diff_q = function(X,Q,
                        max_it  = 500,
                        epsilon = 1e-05,
                        verbose = TRUE,
                        #
                        # Hyper parameters
                        #
                        A_0 = NULL,
                        B_0 = NULL,
                        delta_0 = NULL,
                        ommega_0 = NULL,
                        #
                        # Other settings
                        #
                        Test_versions,
                        test_order,
                        model="General",
                        random_start = FALSE
){

  indI <- sapply(X, nrow)[1] # Assume all individuals take all tests.
  indK <- ncol(Q[[1]]) # Assume attributes are all same across time and individual.
  indT <- length(Q) # Assume time points is all same across individual.
  indJt <- sapply(Q,nrow) # Assume items presented at a time is just same across time.
  indL <- 2^indK

  #
  # All attribute pattern matrix
  #
  not_zero_q_t <- lapply(Q,function(y)apply(y, 1, function(x) which(x != 0)))
  A <- as.matrix(expand.grid(lapply(1:indK, function(x)rep(0:1))))
  A_jt <- lapply(not_zero_q_t, function(y)lapply(y, function(x) A[,x,drop=F]))

  A_red <- lapply(  A_jt, function(z)lapply(z , function(x) apply(x,1,function(y) paste0(y,collapse = ""))))

  #
  # Unique correct item response probability label for each time point and item.
  #
  A_red_uni <- lapply(  A_jt,function(z)lapply(z , function(x) unique(apply(x,1,function(y) paste0(y,collapse = "")))))

  #
  # Make G-matrix
  #
  # G_j <- lapply(1:J, function(j) t(sapply(A_red_uni[[j]], function(x) x == A_red[[j]] ))*1 )

  if(model == "General"){

    G_jt <- lapply(1:indT,function(time_t)lapply(1:indJt[[time_t]], function(j) outer(A_red_uni[[time_t]][[j]], A_red[[time_t]][[j]], function(x,y) (x == y)*1  )))


    att_pat <- apply(A, 1, function(x) paste0(x, collapse=""))
    for(time_t in 1:indT){
      for(j in 1:indJt[time_t]) {
        colnames(G_jt[[time_t]][[j]]) <- att_pat
        row.names(G_jt[[time_t]][[j]]) <- A_red_uni[[time_t]][[j]]
      }
    }

  }else if(model == "DINA"){

    G_jt <- lapply(1:indT,function(time_t)lapply(1:indJt[time_t], function(j)matrix(0,ncol=indL,nrow=2)))



    att_pat <- apply(A, 1, function(x) paste0(x, collapse=""))
    for(time_t in 1:indT){
      for(j in 1:indJt[time_t]) {
        temp_eta <- apply(t(t(A) ^ Q[[time_t]][j,]),1, prod)
        G_jt[[time_t]][[j]][1,] <- 1 - temp_eta
        G_jt[[time_t]][[j]][2,] <- temp_eta

        colnames(G_jt[[time_t]][[j]]) <- att_pat
        row.names(G_jt[[time_t]][[j]]) <- c("0","1")
      }
    }


  } else {
    cat("Error: Specify model General or DINA.\n")
    #break()
  }





  #
  # Hyper parameter
  #
  if(is.null(delta_0) ){
    delta_0 = rep(1,indL)# For π
    #delta_0 = rep(1/indL,indL)# For π

  }

  if(is.null(ommega_0) ){
    ommega_0 = matrix(1,indL,indL)# For Tau matrix
    #ommega_0 = matrix(1/indL,indL,indL)# For Tau matrix

    for(l in 1:indL){
      for(ld in 1:indL){
        dif_pat <- A[l,] - A[ld,]
        ommega_0[l,ld] <- ifelse(any(dif_pat > 0), 0, 1)
      }
    }


  }

  #
  # Weekly Monotonicity constraint
  #
  number_of_attributes <- lapply(A_red_uni,function(y)lapply(y, function(x) sapply(strsplit(x, ""), function(y)sum(as.numeric(y))) ) )


  if(model == "DINA") {number_of_attributes <- lapply(1:indT,function(time_t){lapply(1:indJt[time_t],function(j)c(0,1))})}

  if(is.null(A_0)){

    A_0_hyperparam <- lapply(number_of_attributes, function(x)seq(from = 1+epsilon, to = 2, length.out = max(unlist( x))+1) )
    A_0 <- vector("list", length = indT)
    for(time_t in 1:indT){
      A_0[[time_t]] <-lapply( number_of_attributes[[time_t]] , function(x){A_0_hyperparam[[time_t]][x + 1] })
    }

  }

  if(is.null(B_0)){
    B_0_hyperparam <- lapply(number_of_attributes, function(x)seq(from = 2, to = 1+epsilon, length.out = max(unlist( x))+1) )
    B_0 <- vector("list", length = indT)
    for(time_t in 1:indT){
      B_0[[time_t]] <-lapply( number_of_attributes[[time_t]] , function(x){B_0_hyperparam[[time_t]][x + 1] })
    }

  }




  #
  # Initialization
  #
  if(random_start == TRUE){
    E_z_itl_temp <- lapply(1:indT, function(time_t)matrix(stats::runif(indI*indL) ,ncol=indL, nrow = indI))
    E_z_itl_temp <- lapply(E_z_itl_temp, function(x)diag(1/rowSums(x)) %*% x)

    E_z_itl <-  array(0, dim=c(indI, indL, indT))
    for(time_t in 1:indT){
      E_z_itl[,,time_t] <- E_z_itl_temp[[time_t]]
    }


    E_z_itl_z_itm1l <- array(0, dim=c(indI, indL, indL, indT-1))

    for(i in 1:indI){
      for(time_t in 1:(indT-1)){

        E_z_itl_z_itm1l_temp <- matrix(stats::runif(indL*indL) ,ncol=indL, nrow = indL)
        E_z_itl_z_itm1l_temp[ommega_0==0] <-0
        E_z_itl_z_itm1l_temp <- E_z_itl_z_itm1l_temp/ sum(E_z_itl_z_itm1l_temp)

        E_z_itl_z_itm1l[i,,,time_t] <- E_z_itl_z_itm1l_temp
      }
    }

  }else{

    E_z_itl <-  array(0, dim=c(indI, indL, indT))
    for(time_t in 1:indT){
      E_z_itl[,,time_t] <- matrix(1/indL, nrow=indI, ncol=indL)

    }

    E_z_itl_z_itm1l <- array(0, dim=c(indI, indL, indL, indT-1))

    for(i in 1:indI){
      for(time_t in 1:(indT-1)){

        E_z_itl_z_itm1l_temp <- matrix(1/(indL*indL),indL,indL)
        E_z_itl_z_itm1l_temp[ommega_0==0] <-0
        E_z_itl_z_itm1l[i,,,time_t] <- E_z_itl_z_itm1l_temp/ sum(E_z_itl_z_itm1l_temp)

      }
    }

  }

  #
  # Evidence of Lower Bound
  #
  llb_fun <- function(delta_ast,delta_0,ommega_ast,ommega_0, A_ast, A_0, B_ast, B_0, log_zeta_sum){

    A_0_unlist <- unlist(A_0)
    B_0_unlist <- unlist(B_0)
    A_ast_unlist <- unlist(A_ast)
    B_ast_unlist <- unlist(B_ast)

    tmp1 <- sum( lbeta(A_ast_unlist, B_ast_unlist) - lbeta(A_0_unlist, B_0_unlist) + (A_0_unlist - A_ast_unlist)*(digamma(A_ast_unlist) - digamma(A_ast_unlist+B_ast_unlist)) + (B_0_unlist - B_ast_unlist)*( digamma(B_ast_unlist)-digamma(A_ast_unlist+B_ast_unlist)) )

    tmp2 <- (sum(lgamma(delta_ast)) - lgamma(sum(delta_ast)))  - (sum(lgamma(delta_0)) - lgamma(sum(delta_0))) + sum((delta_0 - delta_ast)*(digamma(delta_ast) - digamma(sum(delta_ast))) )

    tmp3 <- 0
    for(l in 1:indL){
      ommega_not_0   <- ommega_0[l,]!=0
      tmp3 <- tmp3 + (sum(lgamma(ommega_ast[l,ommega_not_0])) - lgamma(sum(ommega_ast[l,ommega_not_0])))  - (sum(lgamma(ommega_0[l,ommega_not_0])) - lgamma(sum(ommega_0[l,ommega_not_0]))) + sum((ommega_0[l,ommega_not_0] - ommega_ast[l,ommega_not_0])*(digamma(ommega_ast[l,ommega_not_0]) - digamma(sum(ommega_ast[l,ommega_not_0]))) )
    }


    tmp1 + tmp2 + tmp3 + log_zeta_sum
  }


  #
  # Make objects for variational parameters
  #
  E_log_theta <- E_log_1_theta <- B_ast <- A_ast <- A_0
  delta_ast <- delta_0
  ommega_ast <- ommega_0
  ommega_zero_elem <- ommega_0 == 0


  b_z_it <- f_z_it <- array(0, dim=c(indI, indL ,indT) )
  gamma_t_x_it <- matrix(0, nrow=indI,ncol=indT)
  P_til_x_it_z_it <- array(0, dim=c(indI,indL,indT))



  X_reord <- X
  for(i in 1:indI){
    for(time_t in 1:indT){
      X_reord[[test_order[Test_versions[i],time_t ]]][i,] <- X[[time_t]][i,]
    }
  }

  m = 1

  l_lb = rep(NA_real_, max_it+1)
  l_lb[1] = 100

  for(m in 1:max_it){
    if(verbose){
      #cat("m = ",m,": l_lb = ",l_lb[m],"\n")
      cat("\riteration = ", m, sprintf(": l_lb = %.05f", l_lb[m]))
    }

    #
    # M-step and Calculation of Expectations
    #
    delta_ast <- colSums(E_z_itl[,,1]) + delta_0
    ommega_ast <- apply(E_z_itl_z_itm1l, c(2,3),sum) + ommega_0 #Check this point


    E_log_pi      = digamma(delta_ast) - digamma(sum(delta_ast))
    #digamma_omega <- try(digamma(ommega_ast))
    #digamma_omega[ommega_zero_elem]  <- 0
    E_log_tau     = try(digamma(ommega_ast), silent = T)  - digamma(rowSums(ommega_ast))
    E_log_tau[ommega_zero_elem]  <- 0


    #
    # Reorder
    #
    E_z_itl_reord <- E_z_itl
    for(i in 1:indI){
      for(time_t in 1:indT){
        E_z_itl_reord[i,,test_order[Test_versions[i],time_t]] <-  E_z_itl[i,,time_t]

      }
    }


    for(time_t in 1:indT){

      A_ast[[time_t]]  <- lapply(1:indJt[[time_t]], function(j) t(G_jt[[time_t]][[j]] %*% (t(E_z_itl_reord[,,time_t]) %*% X_reord[[time_t]][,j])     + A_0[[time_t]][[j]] ))
      B_ast[[time_t]]  <- lapply(1:indJt[[time_t]], function(j) t(G_jt[[time_t]][[j]] %*% (t(E_z_itl_reord[,,time_t]) %*% (1-X_reord[[time_t]][,j])) + B_0[[time_t]][[j]] ))

      E_log_theta[[time_t]]   = lapply(1:indJt[[time_t]], function(j) digamma(A_ast[[time_t]][[j]]) - digamma(A_ast[[time_t]][[j]] + B_ast[[time_t]][[j]]))
      E_log_1_theta[[time_t]] = lapply(1:indJt[[time_t]], function(j) digamma(B_ast[[time_t]][[j]]) - digamma(A_ast[[time_t]][[j]] + B_ast[[time_t]][[j]]))
    }



    #
    # E-step
    #

    for(time_t in 1:indT){
      temp <- matrix(0, nrow = indI, ncol=indL)
      for(i in 1:indI){

        for(j in 1:indJt[time_t]){
          #          temp <- temp +  ( X[[time_t]][,j]%*% E_log_theta[[time_t]][[j]] + (1-X[[time_t]][,j]) %*% E_log_1_theta[[time_t]][[j]]) %*% G_jt[[time_t]][[j]]
          temp[i,] <- temp[i,] +  ( X[[time_t]][i,j]* E_log_theta[[test_order[Test_versions[i],time_t]]][[j]] + (1-X[[time_t]][i,j]) * E_log_1_theta[[test_order[Test_versions[i],time_t]]][[j]]) %*% G_jt[[test_order[Test_versions[i],time_t]]][[j]]
        }
      }
      P_til_x_it_z_it[,,time_t] <-  exp(temp)
    }


    f_z_it[,,1] <- exp(t(log(t(P_til_x_it_z_it[,,1])) + E_log_pi))
    gamma_t_x_it[,1] <- rowSums(f_z_it[,,1])
    f_z_it[,,1] <- f_z_it[,,1]/gamma_t_x_it[,1] # Normarize

    b_z_it[,,indT] <- 1

    #
    # Recursive calculation
    #

    exp_E_log_tau <- exp(E_log_tau)
    exp_E_log_tau[ommega_zero_elem] <- 0

    for(time_t in 2:indT){
      #
      # calc f
      #

      f_z_it[,,time_t] <-  P_til_x_it_z_it[,,time_t] * (f_z_it[,,time_t-1] %*% exp_E_log_tau)
      gamma_t_x_it[,time_t] <- rowSums(f_z_it[,,time_t])
      f_z_it[,,time_t] <- f_z_it[,,time_t]/gamma_t_x_it[,time_t] # Normarize

      #
      # calc b
      #
      b_z_it[,,indT - time_t + 1] <-  (P_til_x_it_z_it[,,indT - time_t + 2] * b_z_it[,,indT - time_t + 2]) %*% t(exp_E_log_tau)

      b_z_it[,,indT - time_t + 1] <- b_z_it[,,indT - time_t + 1] / rowSums(b_z_it[,,indT - time_t + 1])

    }



    E_z_itl <- f_z_it*b_z_it
    E_z_itl_temp <-  apply(E_z_itl, c(1,3),sum)

    for(time_t in 1:indT){
      E_z_itl[,,time_t] <- E_z_itl[,,time_t]/E_z_itl_temp[,time_t]
    }


    for(l in 1:indL){
      for(time_t in 2:indT){
        E_z_itl_z_itm1l[,l,,time_t-1] <- t(t(P_til_x_it_z_it[,,time_t]*b_z_it[,,time_t]*f_z_it[,l,time_t-1]) * exp_E_log_tau[l,])
      }
    }

    E_z_itl_z_itm1l_temp <- apply(E_z_itl_z_itm1l, c(1,4),sum)


    for(i in 1:indI){
      for(time_t in 1:(indT-1)){
        E_z_itl_z_itm1l[i,,,time_t] <- E_z_itl_z_itm1l[i,,,time_t]/E_z_itl_z_itm1l_temp[i,time_t]
      }
    }

    log_zeta_sum <- sum(log(gamma_t_x_it))

    l_lb[m+1] <- llb_fun(delta_ast,delta_0,ommega_ast,ommega_0, A_ast, A_0, B_ast, B_0, log_zeta_sum)

    if(abs(l_lb[m] - l_lb[m+1]) < epsilon){
      if(verbose){
        cat("\nreached convergence.")
      }
      break()
    }

  }
  l_lb <- l_lb[-1]
  #  plot(l_lb,type="l")


  #
  # Calculation of mean and sd of VB posteriors
  #
  delta_sum <- sum(delta_ast)
  pi_est <-  delta_ast/delta_sum
  pi_sd <-sqrt(delta_ast*(delta_sum - delta_ast)/(delta_sum^2*(delta_sum+1)) )
  names(pi_est) <- att_pat
  names(pi_sd)  <- att_pat

  ommega_sum <- rowSums(ommega_ast)
  Tau_est <-  ommega_ast/ommega_sum
  Tau_sd <- matrix(0, indL, indL)
  for(l in 1:indL) Tau_sd[,l] <- sqrt(ommega_ast[,l]*(ommega_sum - ommega_ast[,l])/(ommega_sum^2*(ommega_sum+1)) )

  colnames(Tau_est) <- att_pat
  row.names(Tau_est) <- att_pat
  colnames(Tau_sd) <- att_pat
  row.names(Tau_sd) <- att_pat


  theta_sd <- theta_est <- vector("list",indT)
  for(time_t in 1:indT){
    theta_est[[time_t]] <- mapply(function(x,y) x/(x+y), A_ast[[time_t]], B_ast[[time_t]])
    theta_sd[[time_t]] <- mapply(function(x,y) sqrt((x*y)/(((x+y)^2) *(x+y+1)) ),  A_ast[[time_t]], B_ast[[time_t]])

  }

  #
  # MAP and EAP of attribute mastery.
  #
  post_max_class <- matrix(0, nrow=indI, ncol=indT)
  EAP_att_pat  <- att_master_prob  <- MAP_att_pat  <- lapply(1:indT, function(time_t) matrix(0, nrow=indI, ncol=indK))
  for(time_t in 1:indT){
    post_max_class[,time_t] <- apply(E_z_itl[,,time_t], 1, function(x)which.max(x) )
    MAP_att_pat[[time_t]] <- A[post_max_class[,time_t],]

    att_master_prob[[time_t]] <- E_z_itl[,,time_t] %*% A
    EAP_att_pat[[time_t]] <- (att_master_prob[[time_t]] > 0.5)*1

  }




  list(theta_est = theta_est,
       theta_sd = theta_sd,
       pi_est = pi_est,
       pi_sd = pi_sd,
       Tau_est = Tau_est,
       Tau_sd = Tau_sd,
       post_max_class = post_max_class,
       MAP_att_pat = MAP_att_pat,
       att_master_prob = att_master_prob,
       EAP_att_pat = EAP_att_pat,
       A_ast = A_ast,
       B_ast = B_ast,
       delta_ast   = delta_ast,
       ommega_ast   = ommega_ast,
       E_z_itl = E_z_itl,
       E_z_itl_z_itm1l = E_z_itl_z_itm1l,
       A_0 = A_0,
       B_0 = B_0,
       delta_0 = delta_0,
       ommega_0 = ommega_0,
       l_lb = l_lb,
       gamma_t_x_it = gamma_t_x_it,
       log_zeta_sum = log_zeta_sum,
       E_z_itl = E_z_itl,
       E_z_itl_z_itm1l = E_z_itl_z_itm1l,
       A = A,
       Q = Q,
       X = X,
       G_jt = G_jt,
       m = m)
}

#' @include variationalDCM.R
#' @export
#' @describeIn variationalDCM print summary information

summary.variationalDCM = function(object){

  output = list(
    Estimated.model.parameters = object$model_params,
    Attribute.mastery.patterns = object$att_pat_est,
    ELBO = tail(object$l_lb, n=1),
    time = object$time
    )
  class(output) = "summary.variationalDCM"
  output
}

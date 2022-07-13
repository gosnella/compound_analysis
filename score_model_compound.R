score_model_compound <- function (Observed, Predicted, score = c("log", "Brier", "spherical"),
                                  weights)
{
  score <- match.arg(score)
  if (!is.factor(Observed))
    stop("Observed must be a factor.")
  l <- nlevels(Observed)
  n <- length(Observed)
  Predicted <- as.matrix(Predicted)
  if (nrow(Predicted) != l)
    stop("nrow(Predicted) not= nlevels(Observed).")
  if (ncol(Predicted) != n)
    stop("ncol(Predicted) not= length(Observed).")
  if (max(abs(colSums(Predicted) - 1)) > sqrt(.Machine$double.eps))
    stop("Some colSums(Predicted) not= 1.")
  y <- as.integer(Observed)
  py <- Predicted[cbind(y, 1:n)]
  if (missing(weights)) {
    avg <- function(x) mean(x)
  }
  else {
    avg <- function(x) weighted.mean(x, weights)
  }
  out <- switch(score, log = -avg(log(py)),
                Brier = 1 - avg(py + py - colSums(Predicted * Predicted)),
                spherical = -avg(py/sqrt(colSums(Predicted * Predicted)))
                                                                                                                                       Predicted))))
  out
}
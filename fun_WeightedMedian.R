#' Weighted quantile
#'
#' Function copied from **spatstat** package.
#'
#' @param x Vector of values
#' @param w Vector of weights
#' @param probs Vector of probabilities
#' @param na.rm Ignore missing data?
#' @export
WeightedQuantile <- function(x, w, probs=seq(0,1,0.25), na.rm=TRUE) {
  x <- as.numeric(as.vector(x))
  w <- as.numeric(as.vector(w))
  if(anyNA(x) || anyNA(w)) {
    ok <- !(is.na(x) | is.na(w))
    x <- x[ok]
    w <- w[ok]
  }
  stopifnot(all(w >= 0))
  if(all(w == 0)) stop("All weights are zero", call.=FALSE)
  #'
  oo <- order(x)
  x <- x[oo]
  w <- w[oo]
  Fx <- cumsum(w)/sum(w)
  #'
  result <- numeric(length(probs))
  for(i in seq_along(result)) {
    p <- probs[i]
    lefties <- which(Fx <= p)
    if(length(lefties) == 0) {
      result[i] <- x[1]
    } else {
      left <- max(lefties)
      result[i] <- x[left]
      if(Fx[left] < p && left < length(x)) {
        right <- left+1
        y <- x[left] + (x[right]-x[left]) * (p-Fx[left])/(Fx[right]-Fx[left])
        if(is.finite(y)) result[i] <- y
      }
    }
  }
  names(result) <- paste0(format(100 * probs, trim = TRUE), "%")
  return(result)
}


#' Weighted median
#'
#' Function copied from **spatstat** package.
#'
#' @param x Vector of values
#' @param w Vector of weights
#' @param na.rm Ignore missing data?
#' @export
WeightedMedian <- function(x, w, na.rm=TRUE) {
  unname(WeightedQuantile(x, probs=0.5, w=w, na.rm=na.rm))
}
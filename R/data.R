require(xts)

#' Bitcoin Price Data 2016-19
#'
#' Time series of Bitcoin closing prices from 31 December 2015 to 31 December 2019 (1044 values).
#' This permits the calculation of 4 calendar years of returns.
#'
#' @docType data
#'
#' @usage data(bitcoin)
#'
#' @format An object of class \code{"xts"}.
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(bitcoin)
#' plot(bitcoin)
#' X <- (diff(log(bitcoin))[-1]) * 100
#' plot(X)
"bitcoin"

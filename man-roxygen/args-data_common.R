#' @param data a data frame
#' @param trt column of `data` specifying treatments, coded using integers,
#'   strings, or factors
#' @param trt_ref reference treatment for the network, as a single integer,
#'   string, or factor. By default takes the first treatment in a "natural" sort
#'   order.
#' @param study column of `data` specifying the studies, coded using integers,
#'   strings, or factors
#' @param y column of `data` specifying a continuous outcome

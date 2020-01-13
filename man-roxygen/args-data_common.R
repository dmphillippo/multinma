#' @param data a data frame
#' @param trt column of `data` specifying treatments, coded using integers,
#'   strings, or factors
#' @param trt_ref reference treatment for the network, as a single integer,
#'   string, or factor. By default takes the first treatment in a "natural" sort
#'   order.
#' @param trt_class column of `data` specifying treatment classes, coded using
#'   integers, strings, or factors. By default, no classes are specified.
#' @param study column of `data` specifying the studies, coded using integers,
#'   strings, or factors
#' @param y column of `data` specifying a continuous outcome

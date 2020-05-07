#' @examples
#' # Set up network of smoking cessation data
#' head(smoking)
#'
#' smk_net <- set_agd_arm(smoking,
#'                        study = studyn,
#'                        trt = trtc,
#'                        r = r,
#'                        n = n,
#'                        trt_ref = "No intervention")
#'
#' # Print details
#' smk_net
#'

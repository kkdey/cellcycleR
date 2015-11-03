#' @title tan inverse of ratio of two numbers, mapping from 0 to 2pi.
#'
#' @param beta2 The numerator for the tan inverse function
#' @param beta1 The denominator for the tan inverse function
#'
#' @description The function computes the tan inverse of ratio of two numbers and maps it to the cycle
#'              from 0 to 2pi. A wrapper function for the atan function in R
#'
#' @author  Kushal K Dey
#' @export
#' @examples
#'  atan(4,1)
#'  atan(2,3)
#'  atan(0,0)
#'

atan3 <- function(beta2, beta1)
{
  if (beta1 ==0 & beta2 ==0)
    stop("encountered a 0/0 scenario")
  if (beta1 == Inf | beta1 == -Inf)
    stop("the denominator value is Inf/-Inf")
  if (beta2 == Inf | beta2 == -Inf)
    stop("the numerator value is Inf/-Inf")
  if (beta1 > 0)
    v <- atan(beta2/beta1);
  if(beta2 >=0 & beta1 <0)
    v <- pi + atan(beta2/beta1);
  if(beta2 <0 & beta1 <0)
    v <- -pi + atan(beta2/beta1);
  if(beta2 >0 & beta1==0)
    v <- pi/2;
  if(beta2 <0 & beta1==0)
    v <- - (pi/2);
  if (v < 0)
    v <- v + 2*pi;
  return(v)
}

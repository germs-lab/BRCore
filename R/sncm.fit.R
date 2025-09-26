#' Fit Sloan Neutral Community Model (SNCM)
#'
#' This function fits the Sloan Neutral Community Model to community data to assess
#' neutral processes in community assembly. It compares the neutral model against
#' binomial and Poisson null models using various statistical measures.
#'
#' @param spp A matrix or data frame where rows represent communities (samples) and
#'   columns represent species/taxa. Values should be abundances or counts.
#' @param pool Optional. A matrix or data frame representing the source pool for
#'   calculating relative abundances. If NULL, the source pool is calculated from
#'   the spp matrix. Default is NULL.
#' @param stats Logical. If TRUE, returns fit statistics including AIC, BIC, R-squared,
#'   and RMSE for model comparison. If FALSE, returns predicted vs observed frequencies
#'   with confidence intervals. Default is TRUE.
#' @param taxon Optional. A data frame containing taxonomic information to merge with
#'   results when stats=FALSE. Should have row names matching species names. Default is NULL.
#'
#' @return If stats=TRUE, returns a data frame with model fit statistics including:
#'   \itemize{
#'     \item m: Migration rate parameter from NLS fit
#'     \item m.ci: Confidence interval for m parameter
#'     \item m.mle: Migration rate from maximum likelihood estimation
#'     \item maxLL, binoLL, poisLL: Log-likelihood values for SNCM, binomial, and Poisson models
#'     \item Rsqr, Rsqr.bino, Rsqr.pois: R-squared values for each model
#'     \item RMSE, RMSE.bino, RMSE.pois: Root mean squared error for each model
#'     \item AIC, BIC: Information criteria for model selection
#'     \item N: Average number of individuals per community
#'     \item Samples: Number of samples/communities
#'     \item Richness: Number of taxa analyzed
#'     \item Detect: Detection limit (1/N)
#'   }
#'   If stats=FALSE, returns a data frame with observed and predicted frequencies
#'   along with confidence intervals for visualization.
#'
#' @details The function implements the Sloan Neutral Community Model which predicts
#'   the occurrence frequency of taxa based on their relative abundance in the source
#'   pool and a migration parameter (m). The model assumes neutral processes govern
#'   community assembly. Three models are compared: SNCM (beta distribution),
#'   binomial null model, and Poisson null model.
#'
#' @references
#'   Sloan, W.T. et al. (2006) Quantifying the roles of immigration and chance in
#'   shaping prokaryote community structure. Environmental Microbiology, 8, 732-740.
#'
#' @examples
#' \dontrun{
#' # Generate example community data
#' spp_data <- matrix(rpois(100*50, 5), nrow=100, ncol=50)
#' colnames(spp_data) <- paste0("Species_", 1:50)
#' 
#' # Fit SNCM and get statistics
#' fit_stats <- sncm.fit(spp_data, stats=TRUE)
#' print(fit_stats)
#' 
#' # Get predictions for plotting
#' predictions <- sncm.fit(spp_data, stats=FALSE)
#' }
#'
#' @importFrom minpack.lm nlsLM
#' @importFrom Hmisc binconf
#' @importFrom stats4 mle
#' @importFrom stats pbeta
#' 
#' @export
sncm.fit <- function(spp, pool=NULL, stats=TRUE, taxon=NULL){
    
    base::options(warn=-1)
    
    #Calculate the number of individuals per community
    N <- base::mean(base::apply(spp, 1, sum))
    
    #Calculate the average relative abundance of each taxa across communities
    if(base::is.null(pool)){
        p.m <- base::apply(spp, 2, mean)
        p.m <- p.m[p.m != 0]
        p <- p.m/N
    } else {
        p.m <- base::apply(pool, 2, mean)
        p.m <- p.m[p.m != 0]
        p <- p.m/N
    }
    
    #Calculate the occurrence frequency of each taxa across communities
    spp.bi <- 1*(spp>0)
    freq <- base::apply(spp.bi, 2, mean)
    freq <- freq[freq != 0]
    
    #Combine
    C <- base::merge(p, freq, by=0)
    C <- C[base::order(C[,2]),]
    C <- base::as.data.frame(C)
    C.0 <- C[!(base::apply(C, 1, function(y) base::any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
    p <- C.0[,2]
    freq <- C.0[,3]
    base::names(p) <- C.0[,1]
    base::names(freq) <- C.0[,1]
    
    #Calculate the limit of detection
    d = 1/N
    
    ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
    m.fit <- minpack.lm::nlsLM(freq ~ stats::pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
    m.ci <- stats::confint(m.fit, 'm', level=0.95)
    
    ##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
    sncm.LL <- function(m, sigma){
        R = freq - stats::pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
        R = stats::dnorm(R, 0, sigma)
        -base::sum(base::log(R))
    }
    m.mle <- stats4::mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=base::length(p))
    
    ##Calculate Akaike's Information Criterion (AIC)
    aic.fit <- stats::AIC(m.mle, k=2)
    bic.fit <- stats::BIC(m.mle)
    
    ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
    freq.pred <- stats::pbeta(d, N*stats::coef(m.fit)*p, N*stats::coef(m.fit)*(1-p), lower.tail=FALSE)
    Rsqr <- 1 - (base::sum((freq - freq.pred)^2))/(base::sum((freq - base::mean(freq))^2))
    RMSE <- base::sqrt(base::sum((freq-freq.pred)^2)/(base::length(freq)-1))
    
    pred.ci <- Hmisc::binconf(freq.pred*base::nrow(spp), base::nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
    
    ##Calculate AIC for binomial model
    bino.LL <- function(mu, sigma){
        R = freq - stats::pbinom(d, N, p, lower.tail=FALSE)
        R = stats::dnorm(R, mu, sigma)
        -base::sum(base::log(R))
    }
    bino.mle <- stats4::mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=base::length(p))
    
    aic.bino <- stats::AIC(bino.mle, k=2)
    bic.bino <- stats::BIC(bino.mle)
    
    ##Goodness of fit for binomial model
    bino.pred <- stats::pbinom(d, N, p, lower.tail=FALSE)
    Rsqr.bino <- 1 - (base::sum((freq - bino.pred)^2))/(base::sum((freq - base::mean(freq))^2))
    RMSE.bino <- base::sqrt(base::sum((freq - bino.pred)^2)/(base::length(freq) - 1))
    
    bino.pred.ci <- Hmisc::binconf(bino.pred*base::nrow(spp), base::nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
    
    ##Calculate AIC for Poisson model
    pois.LL <- function(mu, sigma){
        R = freq - stats::ppois(d, N*p, lower.tail=FALSE)
        R = stats::dnorm(R, mu, sigma)
        -base::sum(base::log(R))
    }
    pois.mle <- stats4::mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=base::length(p))
    
    aic.pois <- stats::AIC(pois.mle, k=2)
    bic.pois <- stats::BIC(pois.mle)
    
    ##Goodness of fit for Poisson model
    pois.pred <- stats::ppois(d, N*p, lower.tail=FALSE)
    Rsqr.pois <- 1 - (base::sum((freq - pois.pred)^2))/(base::sum((freq - base::mean(freq))^2))
    RMSE.pois <- base::sqrt(base::sum((freq - pois.pred)^2)/(base::length(freq) - 1))
    
    pois.pred.ci <- Hmisc::binconf(pois.pred*base::nrow(spp), base::nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
    
    ##Results
    if(stats==TRUE){
        fitstats <- base::data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
        fitstats[1,] <- c(stats::coef(m.fit), stats::coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, base::nrow(spp), base::length(p), d)
        return(fitstats)
    } else {
        A <- base::cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
        A <- base::as.data.frame(A)
        base::colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
        if(base::is.null(taxon)){
            B <- A[base::order(A[,1]),]
        } else {
            B <- base::merge(A, taxon, by=0, all=TRUE)
            base::row.names(B) <- B[,1]
            B <- B[,-1]
            B <- B[base::order(B[,1]),]
        }
        return(B)
    }
}
#' Inverse-Gamma DLM MCMC Simulation
#'
#' MCMC simulation of the reduced Fourier-form Bayesian DLM with conjugate Inverse-Gamma priors.
#' 
#' @param dat vector of length \emph{T} containing a time series of VI data
#' @param psi10 initial value for observation equation precision \eqn{1/\sigma_e^2}
#' @param psi20 initial value for state equation precision \eqn{1/\sigma_w^2}
#' @param s value specifying the period of the Fourier harmonics
#' @param seas the number of Fourier harmonics to include in the model (q)
#' @param mc integer value specifying the number of MCMC iterations to run
#' @param thin degree of thinning of the MCMC chain. Retain every \code{thin}-th iteration, thin=1 retains all iterations.
#' @param priors list containing values of inverse-gamma hyperparameters (\eqn{a_1}, \eqn{a_2}, \eqn{b_1}, \eqn{b_2})
#' @param SVD boolean specifying whether the singular value decomposition (TRUE) or the Cholesky decomposition (FALSE) should be used for matrix operations. Default is TRUE.
#' @param save.theta boolean specifying whether samples of the state vector should be retained. Default is TRUE.
#' @details 
#' blah blah blah
#' @return List containing MCMC samples of \eqn{\sigma_e^2}, \eqn{\sigma_w^2}, and \eqn{\theta_{1:T}} (if \code{save.theta = TRUE}) 
#' \describe{
#' \item{sigma2s:}{mc x 2 matrix containing MCMC samples of \eqn{\sigma_e^2} (first column) and \eqn{\sigma_w^2} (second column)}
#' \item{thetas:}{if \code{save_value = TRUE}, a T x p x mc array containing MCMC samples of the latent states}
#' }
#' @examples
#' data(mtci_data) # attach MTCI data
#' 
#' c.mtci <- scale(mtci_data, scale=FALSE) # center MTCI data
#' 
#' ps <- c(10,200) # initial starting precision values
#' 
#' y <- c.mtci[,5] # pull out 5th location
#' 
#' runMCMC <- Rmcmc_DLM_IG(y, psi10=ps[1], psi20=ps[2], s=46, seas=3, mc=100, thin=1) 
#' # period s = 46
#' # include first 3 harmonics (seas = 3)
#' # run 100 MCMC iterations (mc = 100)
#' # no thinning, (thin = 1)
#' 
#' ### Run multiple chains
#' twochains <- lapply(1:2, function(j) Rmcmc_DLM_IG(y, psi10=ps[j], psi20=ps[j], 
#'                                                   s=46, seas=3, mc=1000, thin=1))  
#' fits <- par_to_match(twochains) # organize output from two chains
Rmcmc_DLM_IG <- function(dat, psi10, psi20, s = 46, 
    seas = 2, mc = 1000, thin = 1, priors = list(a1 = 0.001, a2 = 0.001, b1 = 0.001, b2 = 0.001), SVD = TRUE, save.theta = TRUE) {
    
    # require(abind)
    # require(dlm)
    if (seas > s/2) {
        stop("number of frequencies cannot be more than s/2")
    }
    
    nout <- length(seq(1, mc, thin))
    if (save.theta) {
        output <- rep(list(NA), 2)
    } else {
        output <- list(NA)
    }
    output[[1]] <- matrix(NA, ncol = 2, nrow = nout)
    
    initialvalues <- list(psi1 = psi10, psi2 = psi20)
    
    dlm_mod <- dlm::dlmModTrig(s = s, dV = 1/psi10, dW = 1/psi20, q = seas)
    
    modela <- list(F = dlm::FF(dlm_mod), G = dlm::GG(dlm_mod), m0 = dlm::m0(dlm_mod), C0 = dlm::C0(dlm_mod), V = dlm::V(dlm_mod), W = dlm::W(dlm_mod))
    
    out <- mcmc_local_seas_NA(mc, as.matrix(dat), initialvalues, priors, modela, T_star = sum(!is.na(dat)), SVD)
    
    output[[1]][, 1] <- out[[2]][seq(1, mc, thin)]
    output[[1]][, 2] <- out[[3]][seq(1, mc, thin)]
    colnames(output[[1]]) <- c("sigma2_e", "sigma2_w")
    
    if (save.theta) {
        output[[2]] <- unname(do.call(abind::abind, c(out[[1]], along = 3))[, , seq(1, mc, thin)])
    }
    
    if (save.theta) {
        names(output) <- c("sigma2s", "theta")
    } else {
        names(output) <- c("sigma2s")
    }
    return(output)
}

#' Half-Cauchy DLM MCMC Simulation
#' 
#' MCMC simulation of the reduced Fourier-form Bayesian DLM with uninformative half-Cauchy priors.
#' 
#' @param dat vector of length \emph{T} containing a time series of VI data
#' @param sig2_e0 initial value for observation equation variance \eqn{\sigma_e^2}
#' @param sig2_w0 initial value for state equation variance \eqn{\sigma_w^2}
#' @param s value specifying the period of the Fourier harmonics
#' @param seas the number of Fourier harmonics to include in the model (q)
#' @param mc integer value specifying the number of MCMC iterations to run
#' @param c0 2-dimension vector specifying the scale parameters for the half-Cauchy priors
#' @param thin degree of thinning of the MCMC chain. Retain every \code{thin}-th iteration, thin=1 retains all iterations.
#' @param SVD boolean specifying whether the singular value decomposition (TRUE) or the Cholesky decomposition (FALSE) should be used for matrix operations. Default is TRUE.
#' @param save.theta boolean specifying whether samples of the state vector should be retained. Default is TRUE.
#' @return List containing MCMC samples of \eqn{\sigma_e^2}, \eqn{\sigma_w^2}, and \eqn{\theta_{1:T}} (if \code{save.theta = TRUE}) 
#' \describe{
#' \item{sigma2s:}{mc x 2 matrix containing MCMC samples of \eqn{\sigma_e^2} (first column) and \eqn{\sigma_w^2} (second column)}
#' \item{thetas:}{if \code{save_value = TRUE}, a T x p x mc array containing MCMC samples of the latent states}
#' }
#' @examples
#' data(mtci_data) # attach MTCI data
#' 
#' c.mtci <- scale(mtci_data, scale=FALSE) # center MTCI data
#' 
#' ps <- c(10,200) # initial starting precision values
#' 
#' y <- c.mtci[,5] # pull out 5th location
#' 
#' runMCMC <- Rmcmc_DLM_half_Cauchy(y, sig2_e0=1/ps[1], sig2_w0=1/ps[2], s=46, seas=3, mc=100, thin=1) 
#' # period s = 46
#' # include first 3 harmonics (seas = 3)
#' # run 100 MCMC iterations (mc = 100)
#' # no thinning, (thin = 1)
#' 
#' ### Run multiple chains
#' twochains <- lapply(1:2, function(j) Rmcmc_DLM_half_Cauchy(y, sig2_e0=1/ps[j], 
#'                                            sig2_w0=1/ps[j], s=46, seas=3, mc=1000, thin=1))  
#' fits <- par_to_match(twochains) # organize output from two chains
Rmcmc_DLM_half_Cauchy <- function(dat, sig2_e0, sig2_w0, s = 46, seas = 2, mc = 1000, thin = 1, c0 = c(0.1,0.1),  SVD = FALSE, save.theta = TRUE) {
    # dat: centered time series sig2_e0: initial value for sigma_e^2 sig2_w0: initial value for sigma_e^2 c: 2-dim vector of cauchy
    # scale params for sigma_e, sigma_w priors, respectively mc: number of iterations s: period of seasonality seas: number of
    # harmonics thin: return every 'thin' iteration (if thin > 1, number of iterations returned is floor(mc/thin)) SVD:
    # (TRUE/FALSE) TRUE: use SVD decomp for MVN in FFBS, FALSE: use Cholesky decomp. Better to use SVD to avoid singularities from
    # small variances in Sigma save.theta: (TRUE/FALSE) TRUE: return draws of theta, sigma_e^2, sigma_w^2. FALSE: return draws of
    # sigma_e^2, sigma_w^2
    
    # require(abind)
    # require(dlm)
    
    if (length(sig2_e0) != 1 | length(sig2_w0) != 1) {
        stop("Supply initial values for only one chain.")
    }
    
    if (seas > s/2) {
        stop("number of frequencies cannot be more than s/2")
    }
    
    nout <- length(seq(1, mc, thin))
    if (save.theta) {
        output <- rep(list(NA), 2)
    } else {
        output <- list(NA)
    }
    output[[1]] <- matrix(NA, ncol = 2, nrow = nout)
    
    initialvalues <- list(sig2_e = sig2_e0, sig2_w = sig2_w0)
    
    dlm_mod <- dlm::dlmModTrig(s = s, dV = sig2_e0, dW = sig2_w0, q = seas)
    
    modela <- list(F = dlm::FF(dlm_mod), G = dlm::GG(dlm_mod), m0 = dlm::m0(dlm_mod), C0 = dlm::C0(dlm_mod), V = dlm::V(dlm_mod), W = dlm::W(dlm_mod))
    
    out <- mcmc_seas_IG_NA(mc, as.matrix(dat), initialvalues, c0, modela, T_star = sum(!is.na(dat)), SVD, save.theta)
    
    output[[1]][, 1] <- out[[1]][seq(1, mc, thin)]
    output[[1]][, 2] <- out[[2]][seq(1, mc, thin)]
    colnames(output[[1]]) <- c("sigma2_e", "sigma2_w")
    if (save.theta) {
        output[[2]] <- unname(do.call(abind::abind, c(out[[3]], along = 3))[, , seq(1, mc, thin)])
    }
    
    if (save.theta) {
        names(output) <- c("sigma2s", "theta")
    } else {
        names(output) <- c("sigma2s")
    }
    return(output)
}

#### Set up functions #### faster than apply to compute sums on a 3d array
sums.along <- function(a, i) {
    n <- length(dim(a))
    b <- aperm(a, c(seq_len(n)[-i], i))
    rowSums(b, dims = n - 1)
}

## smooth curves (F theta_t) & uncertainty function
qfits <- function(theta, burn = 1000, thin = 1, q = 2) {
    # require(abind)
    it <- dim(theta[[1]])[3]
    df <- unname(do.call(abind::abind, c(theta, along = 3)))
    keep <- df[, , -c(1:burn, (it + 1):(it + burn))]
    keep <- keep[, , seq(1, dim(keep)[3], thin)]
    if (q == 1) {
        St <- keep[, seq(1, 2 * q, 2), ]
    } else {
        St <- sums.along(keep[, seq(1, 2 * q, 2), ], 2)
    }
    t(apply(St, 1, quantile, c(0.025, 0.5, 0.975)))[-1, ]
}

#' Convert to \code{mcmc.list()}
#' 
#' Function to convert list of samples of \eqn{\sigma_e^2} and \eqn{\sigma_w^2} from multiple MCMC chains from \code{Rmcmc_DLM_IG} or \code{Rmcmc_DLM_half_Cauchy} to \code{coda} object \code{mcmc.list()}. Helpful in order to use MCMC diagnostic functions from the \code{coda} package on chains of \eqn{\sigma_e^2} and \eqn{\sigma_w^2}.
#' 
#' @param mod_out List of model output as obtained from \code{par_to_match}.
#' @param thin Number of iterations to thin the MCMC chains
#' @param burn Number of iterations to use as burn-in.
#' @return \code{mcmc.list()} \code{coda} object containing MCMC samples
#' @examples
#' data(mtci_data) # attach MTCI data
#' 
#' c.mtci <- scale(mtci_data, scale=FALSE) # center MTCI data
#' 
#' ps <- c(10,200) # initial starting precision values
#' 
#' y <- c.mtci[,5] # pull out 5th location
#' 
#' ### Run multiple chains
#' twochains <- lapply(1:2, function(j) Rmcmc_DLM_half_Cauchy(y, sig2_e0=1/ps[j], sig2_w0=1/ps[j], 
#'                                                            s=46, seas=3, mc=1000, thin=1))  
#' mcmc_smp <- par_to_match(twochains) # organize output from two chains
#' 
#' fit_coda <- to_mcmc_list(mcmc_smp) # convert to mcmc.list() object
#' 
#' coda::gelman.diag(fit_coda) # compute Gelman-Rubin diagnostic on chains of sigma_e and sigma_w
to_mcmc_list <- function(mod_out, thin = 1, burn = 0) {
    # require(coda)
    nr <- nrow(mod_out$sigma2_e)
    d <- list()
    for (i in 1:ncol(mod_out$sigma2_e)) {
        d[[i]] <- coda::as.mcmc(data.frame(V.5 = sqrt(mod_out$sigma2_e[seq(burn + 1, nr, thin), i]), W.5 = sqrt(mod_out$sigma2_w[seq(burn + 
            1, nr, thin), i])))
    }
    d <- coda::as.mcmc.list(d)
    return(d)
}

#' Combine MCMC Output from chains
#' 
#' Function to rearrange output from multiple MCMC chains.
#' 
#' @param mod_out List of MCMC output from running multiple MCMC chains using \code{lapply()}, or \code{doParallel::mclapply()}
#' @return List of organized MCMC samples.
#' \describe{
#' \item{sigma2_e}{mc x nchains matrix containing samples of \eqn{\sigma_e^2}}
#' \item{sigma2_w}{mc x nchains matrix containing samples of \eqn{\sigma_w^2}}
#' \item{theta}{list of length nchains containing samples of \eqn{\theta}}
#' }
#' @examples
#' data(mtci_data) # attach MTCI data
#' 
#' c.mtci <- scale(mtci_data, scale=FALSE) # center MTCI data
#' 
#' ps <- c(10,200) # initial starting precision values
#' 
#' y <- c.mtci[,5] # pull out 5th location
#' 
#' ### Run multiple chains
#' twochains <- lapply(1:2, function(j) Rmcmc_DLM_half_Cauchy(y, sig2_e0=1/ps[j], sig2_w0=1/ps[j], 
#'                                                            s=46, seas=3, mc=1000, thin=1))  
#' fits <- par_to_match(twochains) # organize output from two chains
par_to_match <- function(mod_out) {
    sigma2_e <- sapply(mod_out, function(u) u$sigma2s[, 1])
    sigma2_w <- sapply(mod_out, function(u) u$sigma2s[, 2])
    if (length(mod_out[[1]]) == 2) {
        theta <- lapply(mod_out, function(u) u$theta)
        return(list(sigma2_e = sigma2_e, sigma2_w = sigma2_w, theta = theta))
    } else {
        return(list(sigma2_e = sigma2_e, sigma2_w = sigma2_w))
    }
}

#' Importance of Individual Harmonics
#' 
#' Function to compute credible intervals for each state to assess importance of individual harmonics
#' 
#' @param theta T x p x mc array of MCMC samples of theta
#' @param burn number of MCMC samples to use as burn-in
#' @param thin number of iterations to thin MCMC samples
#' @param alpha number specifying level of credible interval, e.g. alpha = 0.05 for 95\% credible interval
#' @return List containing credible intervals for each harmonic state at each timepoint.
#' \describe{
#' \item{ints}{array containing credible interval bounds for each state.}
#' \item{props}{proportion of credible intervals which did not contain 0 for each harmonic.}
#' }
indiv_harm_fits <- function(theta, burn = 1000, thin = 1, alpha = 0.01) {
    q <- dim(theta)[2]/2
    # require(abind) it <- dim(theta[[1]])[3] df <- unname(do.call(abind, c(theta, along=3))) keep <- df[,,-c(1:burn,
    # (it+1):(it+burn))] keep <- keep[,,seq(1,dim(keep)[3],thin)]
    
    St <- theta[, seq(1, 2 * q, 2), ]
    
    ints <- apply(St, c(1, 2), quantile, c(alpha/2, 0.5, 1 - alpha/2))[, -1, ]
    df.m <- reshape2::melt(ints)
    wd1 <- apply(ints, c(2, 3), function(x) (x[1] > 0 & x[3] > 0) | (x[1] < 0 & x[3] < 0))
    return(list(ints = ints, props = colMeans(wd1)))
}

#' Simulate Seasonal DLM
#' 
#' Function to simulate time series from a univariate dynamic linear model with static \eqn{F}, \eqn{G}, \eqn{V}, \eqn{W}.
#' 
#' @param T number specifying the length of the simulated time series
#' @param theta_0 p-dim vector of initial states
#' @param F observation equation matrix
#' @param G state equation transition matrix
#' @param V value specifying observation equation variance
#' @param W value specifying state equation variance
#' @return List of simulated values:
#' \describe{
#' \item{y}{vector of simulated observations}
#' \item{v}{vector of simulated observation errors}
#' \item{w}{vector of simulated innovation errors}
#' \item{theta}{matrix of simulated state vectors}
#' }
simulate_dlm_seas <- function(T, theta_0, F, G, V, W) {
    y <- numeric(length = T)
    v <- numeric(length = T)
    w <- matrix(NA, nrow = T, ncol = nrow(G))
    theta <- matrix(NA, nrow = T + 1, ncol = length(theta_0))
    theta[1, ] <- theta_0
    for (i in 1:T) {
        w[i, ] <- rnorm(nrow(G), 0, sqrt(W))
        v[i] <- rnorm(1, 0, sqrt(V))
        theta[i + 1, ] <- G %*% theta[i, ] + w[i, ]
        y[i] <- F %*% theta[i + 1, ] + v[i]
    }
    return(list(y = y, v = v, w = w, theta = theta))
}

dvnc <- function(y, sig2_e, theta, F) {
    # -2 log(f(y|theta))
    T = sum(!is.na(y))
    as.numeric(T * log(2 * pi) + T * log(sig2_e) + 1/sig2_e * sum((y - c(F %*% t(theta)))^2, na.rm = T))
}

#' Deviance Information Criterion
#' 
#' Function to compute DIC from MCMC samples 
#' 
#' @param y \emph{T}-dim vector of observations 
#' @param sig2es vector containing mc samples of \eqn{\sigma_e^2}
#' @param thetas T x p x mc array of samples of \eqn{\theta}
#' @param F observation equation matrix
#' @return List containing:
#' \describe{
#' \item{DIC}{computed DIC}
#' \item{ENP}{effective number of parameters}
#' }
DIC <- function(y, sig2es, thetas, F) {
    # computer posterior mean deviance
    Ds <- sapply(1:length(sig2es), function(i) dvnc(y, sig2es[i], thetas[, , i], F))
    pmd <- mean(Ds)
    # compute deviance of posterior mean
    sig2e_m <- mean(sig2es)
    theta_m <- apply(thetas, c(1, 2), mean)
    dpm <- dvnc(y, sig2e_m, theta_m, F)
    # compute effective number of parameters
    pd <- pmd - dpm
    # return DIC
    return(list(DIC = pmd + pd, ENP = pd))
}

#' Bayesian Predictive Information Criterion
#' 
#' Compute BPIC from MCMC samples.
#' 
#' @param y \emph{T}-dim vector of observations 
#' @param sig2es vector containing mc samples of \eqn{\sigma_e^2}
#' @param thetas T x p x mc array of samples of \eqn{\theta}
#' @param F observation equation matrix
#' @return List containing:
#' \describe{
#' \item{BPIC}{computed BPIC}
#' \item{ENP}{effective number of parameters}
#' }
BPIC <- function(y, sig2es, thetas, F) {
    # computer posterior mean deviance
    Ds <- sapply(1:length(sig2es), function(i) dvnc(y, sig2es[i], thetas[, , i], F))
    pmd <- mean(Ds)
    # compute deviance of posterior mean
    sig2e_m <- mean(sig2es)
    theta_m <- apply(thetas, c(1, 2), mean)
    dpm <- dvnc(y, sig2e_m, theta_m, F)
    # compute effective number of parameters
    pd <- pmd - dpm
    # return DIC
    return(list(BPIC = pmd + 2 * pd, ENP = pd))
}

#' RMSE
#' 
#' Function to compute root mean squared error.
#' 
#' @param y \emph{T}-dim vector of observations.
#' @param fit object containing posterior summaries of \eqn{SVI_{1:T}} as obtained from \code{MCMC_to_pheno_events()}
#' @return List containing:
#' \describe{
#' \item{rmse}{overall rmse}
#' \item{rmse_year}{rmse computed separately for each year}
#' }
rmse <- function(y, fit) {
    rmse <- sqrt(sum((y - fit$fits[, 2])^2, na.rm = T))
    rmse_year <- sapply(1:5, function(i) sqrt(sum((y[((i - 1) * 46 + 1):(i * 46)] - fit$fits[((i - 1) * 46 + 1):(i * 46), 2])^2, 
        na.rm = T)))
    return(list(rmse, rmse_year))
}

#' Phenological Event Estimation
#' 
#' Function to estimate phenological events from posterior samples obtained from multiple MCMC chains
#' 
#' @param output MCMC ouput from multiple chains, use \code{par_to_match()} around list of multiple chains
#' @param q number of harmonics included in the model
#' @param burnin number of iterations used for burnin
#' @param thin number of iterations to use to thin the chains (\code{thin = 1} keeps all samples)
#' @param cf value specifying the minimum difference between sos and pog in identification rule. If \code{cf = NULL} and \code{cf_calc=TRUE}, \code{cf} is s_cf*posterior mode of \eqn{\sigma_e}, else \code{cf} \eqn{\neq} NULL allows user to specify an alternate value for the cutoff.
#' @param cf_calc boolean where TRUE calculates \code{cf} as s_cf*posterior mode of \eqn{\sigma_e} and FALSE requires a specified cutoff value for \code{cf}.
#' @param s_cf value specifying the multiple of the posterior mode of \eqn{\sigma_e} to use in the computation of \code{cf} if \code{cf_calc = TRUE}
#' @param shift number of composites past January 1 to consider for the previous year (protects against POG that occur around Jan 1st).
#' @param d value specifying the minimum duration of a growing season
#' @param alpha value specifying the confidence level for (1 - alpha)*100\% intervals for phenological event timings
#' @param s period of harmonic time series 
#' @param year_start value specifying the first year of the time range
#' @param year_end value specifying the last year of the time range
#' @param return_points (TRUE/FALSE) should events identified from all MCMC draws be returned, or just summaries?
#' @return List containing:
#' \describe{
#' \item{intervals}{estimates and intervals for all identified phenological events}
#' \item{probs}{estimated probabilities of multiple growing seasons in each year}
#' \item{nseas_year}{highest probability number of growing seasons each year}
#' \item{nused}{number of MCMC samples used to estimate each phenological events}
#' \item{fits}{posterior median estimates and credible intervals of \eqn{SVI_{1:T}}}
#' \item{pts}{if \code{return_points = TRUE}, contains estimated events from all MCMC samples}
#' \item{var_summaries}{posterior summaries for observation and state equation standard deviations (\eqn{\sigma_e}, \eqn{\sigma_w})}
#' }
#' @examples
#' data(mtci_data) # attach MTCI data
#' 
#' c.mtci <- scale(mtci_data, scale=FALSE) # center MTCI data
#' 
#' ps <- c(10,200) # initial starting precision values
#' 
#' y <- c.mtci[,5] # pull out 5th location
#' 
#' ### Run multiple chains
#' twochains <- lapply(1:2, function(j) Rmcmc_DLM_half_Cauchy(y, sig2_e0=1/ps[j], sig2_w0=1/ps[j], 
#'                                                            s=46, seas=3, mc=1000, thin=1))  
#' # period s = 46
#' # include first 3 harmonics (seas = 3)
#' # run 1000 MCMC iterations (mc = 1000)
#' # no thinning, (thin = 1)
#'   
#' mcmc_smp <- par_to_match(twochains) # organize output from two chains
#' 
#' ### Estimate phenological events
#' fit_sum = MCMC_to_pheno_events(mcmc_smp, shift=8, thin=2, d=8, s_cf=3, burnin=100, q=3)
#' 
MCMC_to_pheno_events <- function(output, q = 2, burnin = 1000, thin = 2, cf = NULL, cf_calc = TRUE, s_cf = 2, shift = 8, 
    d = 5, cut2 = 0.8, alpha = 0.1, s = 46, year_start = 2003, year_end = 2007, return_points = FALSE) {
    # output: MCMC ouput from more than one chain of model (use par_to_match() around list of multiple chains) q: number of
    # harmonics in model burnin: number of iterations used for burnin thin: take ever 'thin' iteration, if want to thin (thin=1 is
    # no thinning) cf: min difference between sos and pog in identification rule, if NULL and cf_calc=TRUE, cf = s_cf*posterior
    # mode of sigma_e, else if cf !=NULL allows user to specify an alternate value for cutoff cf_calc: TRUE calculates cf as
    # s_cf*posterior mode of sigma_e, FALSE requires cf!=NULL specified cutoff value s_cf: multiple of posterior mode of sigma_e
    # for cutoff (i.e. s_cf is k, k sigma_e in the paper) shift: number of composites past January 1 to consider for the previous
    # year (protects against POG that occur around Jan 1st) d: This is L in the paper. Min duration of a growing season cut2:
    # cutoff to decide more than one growing season (don't need this anymore) alpha: (1 - alpha)*100% intervals for pheno event
    # timings s: period of time series year_start: first year of time range year_end: last year of time range return_points:
    # (TRUE/FALSE) should points identified from all draws be returned?
    
    # require(MCMCglmm)
    # require(abind)
    
    if (cf_calc) {
        cf <- s_cf * MCMCglmm::posterior.mode(coda::as.mcmc(sqrt(c(output$sigma2_e[seq(burnin, nrow(output$sigma2_e), thin), ]))))
    } else {
        if (is.null(cf)) {
            stop("Must supply value for cf if cf_calc = FALSE!")
        }
    }
    
    it <- dim(output$theta[[1]])[3]
    df <- unname(do.call(abind::abind, c(output$theta, along = 3)))
    keep <- df[, , -c(1:burnin, (it + 1):(it + burnin))]
    keep <- keep[, , seq(1, dim(keep)[3], thin)]
    if (q == 1) {
        St <- as.data.frame(keep[-1, seq(1, 2 * q, 2), ])
    } else {
        St <- as.data.frame(sums.along(keep[, seq(1, 2 * q, 2), ], 2)[-1, ])
    }
    mseas2 <- t(apply(St, 1, quantile, c(0.025, 0.5, 0.975)))
    
    get_ints <- mcmc_to_intervals_cpp_fixed(as.matrix(St), cf, d, shift, s, year_start, year_end, cut2, alpha, return_points)
    
    ints = as.data.frame(do.call(rbind, get_ints$ints$intervals))
    names(ints) = c("Type", "lT", "mT", "uT", "mY")
    ints <- ints[ints$Type != 0, ]
    
    ints$Type = factor(ints$Type)
    levels(ints$Type) = c("POG", "SOS")
    
    probs = as.data.frame(get_ints$ints$probs)
    names(probs) = c("P0", "P1", "P2", "P3")
    
    # sigma_e
    sum_v <- c(MCMCglmm::posterior.mode(coda::mcmc(sqrt(c(output$sigma2_e[-c(1:burnin), ])))), quantile(sqrt(c(output$sigma2_e[-c(1:burnin), ])), 
        c(0.5, 0.025, 0.975)))
    
    # sigma_w
    sum_w <- c(MCMCglmm::posterior.mode(coda::mcmc(sqrt(c(output$sigma2_w[-c(1:burnin), ])))), quantile(sqrt(c(output$sigma2_w[-c(1:burnin), ])), 
        c(0.5, 0.025, 0.975)))
    
    var_summaries <- data.frame(rbind(sum_v, sum_w))
    
    if (return_points) {
        return(list(intervals = ints, probs = probs, nseas_year = get_ints$ints$nseas_year, nused = get_ints$ints$nused, fits = mseas2, 
            pts = get_ints$pts, var_summaries = var_summaries))
    } else {
        return(list(intervals = ints, probs = probs, nseas_year = get_ints$ints$nseas_year, nused = get_ints$ints$nused, fits = mseas2, 
            var_summaries = var_summaries))
    }
}











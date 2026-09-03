##' Relabel the components of a mixture to order them by mean.
##'
##' The sampler does not label the components in any consistent
##' manner, and labels may switch from draw to draw.  This function
##' relabels the components with each draw so that they are ordered in
##' order of their mean.
##'
##' @title Component Relabelling
##' @param s a list of `mcarray` objects returned by [JAGSsample()]
##' @return a list of `mcarray` objects.
##' @export
## ---- relabelComponents
relabelComponents <- function(s) {
  ## Apply a permutation to a mcarray
  applyPerm <- function(x,p) {
    ## If dimensions match
    if(identical(dim(x),dim(p))) {
      att <- attributes(x)
      x <- x[p]
      attributes(x) <- att
    }
    x
  }

  ## Determine the permutation that orders mixture components by mu
  mu <- s$mu
  dim(mu) <- c(dim(mu)[1],prod(dim(mu)[-1]))
  perm <- array(seq_len(length(mu)),dim(mu))
  for(j in seq_len(ncol(mu))) perm[,j] <- perm[order(mu[,j]),j]
  dim(perm) <- dim(s$mu)

  ## Permute components
  s$mu <- applyPerm(s$mu,perm)
  s$w <- applyPerm(s$w,perm)
  s$sigma <- applyPerm(s$sigma,perm)
  s
}
## ----


##' Generate the data to pass to JAGS.
##'
##' This function generates the data list required by the JAGS
##' scripts.  Here `Y` is the matrix of length bin densities for each
##' length bin and replicate haul. Each column of this matrix
##' corresponds to a length bin and each row a replicate haul. The
##' length bins are assumed contiguous and `breaks` must be and
##' ordered vector of bin boundaries of length one greater than the
##' number of bins.
##'
##' Length bins with no non-zero densities in any replicate are not
##' informative, and are removed. Within each bin the densities are
##' sorted so that the non-zero densities are the first `NPos`
##' elements.
##'
##' @title Generate data for JAGS
##' @param Y a matrix of length bin densities.
##' @param breaks the length bin boundaries.
##' @param Ncomp the number of mixture components
##' @return a list with elements
##' * `Y` - the reduced density matrix
##' * `lwr` - the lower limit of the used length bins
##' * `upr` - the upper limit of the used length bins
##' * `Npos` - the number of non-zero densities for each bin
##' * `NHaul` - the number of replicate hauls
##' * `Ncomp` - the number of components to fit
##' @export
## ---- jagsData
jagsData <- function(Y,breaks,Ncomp) {

  ## Retain bins with > 1 non-zero density
  keep <- which(colSums(Y>0)>1)
  Y <- Y[,keep]
  ## Move non-zero densities to the front
  for(k in seq_len(ncol(Y))) Y[,k] <- sort(Y[,k],decreasing=TRUE)

  list(Y=Y,
      lwr=breaks[keep],
      upr=breaks[keep+1],
      Npos=colSums(Y>0),
      NBin=ncol(Y),
      NHaul=nrow(Y),
      NComp=Ncomp)
}
## ----

##' Generate JAGS code for the model priors
##'
##' Translates a list of formula of the form `par ~ dist(a,b,...)` to
##' JAGS code for the model priors.
##'
##' @title Prior JAGS code
##' @param Ncomp the number of mixture components
##' @param prior a list of formula specifying priors
##' @return a list with components
##' * `par` - the parameters for which priors are defined
##' * `code` - JAGS code
##' @export
## ---- jagsPriors
jagsPriors <- function(Ncomp,prior) {

  ## Extract parameter name
  getPar <- function(f) {
    if(length(f) != 3 || as.character(f[[1]])!="~") stop("Invalid prior:",f)
    as.character(f[[2]])
  }

  ## Generate JAGS code
  genCode <- function(f) {
    ## Extract components
    var <- as.character(f[[2]])
    rhs <- f[[3]]
    dist <- as.character(rhs[[1]])

    ## Create parameter list
    rhs[[1]] <- as.name("list")
    pars <- eval(rhs, envir = environment(f))
    pars <- lapply(pars, function(x) format(rep_len(x,Ncomp), trim=TRUE))

    ## Generate JAGS statements
    fmt <- paste0("  %s[%d] ~ %s(", paste(rep_len("%s",length(pars)),collapse =","), ")")
    paste0("\n  ## Prior for ",var,"\n",
           paste(do.call(sprintf,c(list(fmt, var, seq_len(Ncomp), dist), pars)),
                 collapse="\n"))
  }

  list(par=sapply(prior, getPar),
       code=paste(lapply(prior,genCode),collapse="\n"))
}
## ----


##' Generate JAGS code for a Normal mixture.
##'
##' Generate the JAGS code for the mixture of Normals and the model
##' parameters.  Three variant models are provided for the standard
##' deviations sigma of the components
##'
##' * `ind` - sigma is independent of mu, and priors must be provided
##' for either `sigma` or the precision `tau`.
##'
##' * `cv` - sigma is proportional to mu, `sigma = kappa*mu` and a
##' prior must be provided for coefficient of proportionality `kappa`
##'
##' * `lin` - sigma is a linear function of the mean, `sigma = a +
##' b*mu` and priors must be provided for `a` and `b`.
##'
##' @title Normal Mixture JAGS code
##' @param Ncomp the number of mixture components
##' @param prior a list of formula specifying priors
##' @param sigma the model for sigma
##' @return a list with components
##' * `par` - the parameters to sample
##' * `code` - JAGS code
##' @export
## ---- normalMixture
normalMixture <- function(Ncomp,prior,
                          sigma = c("ind","cv","lin")){
  sigma <- match.arg(sigma)
  pr <- jagsPriors(Ncomp, prior)

  code <- "
  for(j in 1:NBin) {
    ## Mixture components
    for(k in 1:NComp) {
      comp[j,k] <- w[k]*(pnorm(upr[j],mu[k],tau[k])-pnorm(lwr[j],mu[k],tau[k]))
    }
    lambda[j] <- sum(comp[j,])
  }
"

  switch(sigma,

         ## sigma is independent of mu
         ind={
           par <- c("w","mu","sigma")
           if("tau" %in% pr$par) {
             mp <- setdiff(c("w","mu","tau"),pr$par)
             if(length(mp)>0) stop("Missing priors: ",paste(mp,collapse=", "))
             code <- paste("  sigma <- 1/sqrt(tau)",code,sep="\n")
           } else {
             mp <- setdiff(c("w","mu","sigma"),pr$par)
             if(length(mp)>0) stop("Missing priors: ",paste(mp,collapse=", "))
             code <- paste("  tau <- 1/sigma^2",code,sep="\n")
           }
         },

         ## sigma is proportional to mu
         cv={
           par <- c("w","mu","sigma","kappa")
           mp <- setdiff(c("w","mu","kappa"),pr$par)
           if(length(mp)>0) stop("Missing priors: ",paste(mp,collapse=", "))
           code <- paste("  sigma <- kappa*mu\n  tau <- 1/sigma^2",code,sep="\n")
         },

         ## sigma a linear function of mu
         lin={
           par <- c("w","mu","sigma","a","b")
           mp <- setdiff(c("w","mu","a","b"),pr$par)
           if(length(mp)>0) stop("Missing priors: ",paste(mp,collapse=", "))
           code <- paste("  sigma <- a+b*mu\n  tau <- 1/sigma^2",code,sep="\n")
         })

  ## Add code for priors
  code <- paste(code, pr$code, sep="\n")
  list(par=par, code=code)
}
## ----

##' Generate JAGS code for delta log Normal bin densities.
##'
##' Generate the JAGS code for modelling delta log Normal bin
##' densities.  The delta log Normal model is parameterized in term of
##' parameters p - the probability of a non-zero response, and m and s
##' - the mean and standard deviation of the log of the non-zero
##' densities.  Three variant models are provided for the p and s parameters
##'
##' * `const` - p and s are constant across bins
##'
##' * `fixed` - p and s vary across bins
##'
##' * `random` - log p and logit s across bins are modelled as Normal
##'   random effects.
##'
##' @title DLN Response JAGS code
##' @param DLN the model for the p and s parameters
##' @return a list with two components
##' * `par` - the parameters to sample
##' * `code` - JAGS code
##' @export
## ---- DLNResponse
DLNResponse <- function(DLN = c("fixed","const","random")) {


  DLN <- match.arg(DLN)
  switch(DLN,

         ## p and s are constant across bins
         const={
           par <- c("p","s")
           code <- "
  ## Priors for nuisance parameters
  p ~ dbeta(1,1)
  t ~ dgamma(0.01,0.01)
  s <- 1/sqrt(t);

  for(j in 1:NBin) {

    ## Log Normal mean
    m[j] <- log(lambda[j]/p)-1/(2*t)

    ## Likelihood
    Npos[j] ~ dbin(p[j],NHaul)
    for(i in 1:Npos[j]) {
      Y[i,j] ~ dlnorm(m[j],t)
    }
  }"
         },


         ## p and s are constant across bins
         fixed={
           par <- c("p","s")
           code <-"
  for(j in 1:NBin) {

    ## Priors for nuisance parameters
    p[j] ~ dbeta(1,1)
    t[j] ~ dgamma(0.01,0.01)
    s[j] <- 1/sqrt(t[j])

    ## Log Normal mean
    m[j] <- log(lambda[j]/p[j])-1/(2*t[j])

    ## Likelihood
    Npos[j] ~ dbin(p[j],NHaul)
    for(i in 1:Npos[j]) {
      Y[i,j] ~ dlnorm(m[j],t[j])
    }
  }"
         },

         ## p and s are random effects
         random={
           par <- c("p","s","p.mu","p.tau","s.mu","s.tau")
           code <-"
  ## Priors for random effect parameters
  p.mu <- dnorm(0,0.01)
  p.tau <- dgamma(0.01,0.01)
  s.mu <- dnorm(0,0.01)
  s.tau <- dgamma(0.01,0.01)

  for(j in 1:NBin) {

    ## Priors for nuisance parameters
    logit(p[j]) ~ dnorm(p.mu,p.tau)
    log(s[j]) ~ dnorm(s.mu,s.tau)
    t[j] <- 1/s[j]^2

    ## Log Normal mean
    m[j] <- log(lambda[j]/p[j])-1/(2*t[j])

    ## Likelihood
    Npos[j] ~ dbin(p[j],NHaul)
    for(i in 1:Npos[j]) {
      Y[i,j] ~ dlnorm(m[j],t[j])
    }
  }"
         })

  list(par=par,code=code)

}
## ----




##' Create a model object for a Normal mixture model with delta log
##' Normal bin densities.
##'
##' Creats a model object for use with [JAGSsample()].  The model fits
##' a Normal mixture to binned length density data, assuming the
##' length bin densities have a delta log Normal distribution.
##'
##' The data are supplied as a matrix `Y` of the density in each
##' length bin and replicate haul where each column of this matrix
##' corresponds to a length bin and each row a replicate haul. The
##' length bins are assumed contiguous and `breaks` must be and
##' ordered vector of bin boundaries of length one greater than the
##' number of bins.
##'
##' The mean density in each bin are assumed to follow a Normal
##' mixture model with `Ncomp` components.  Three variant models are
##' provided for the standard deviations sigma of the components
##'
##' * `ind` - sigma is independent of mu, and priors must be provided
##' for either `sigma` or the precision `tau`.
##'
##' * `cv` - sigma is proportional to mu, `sigma = kappa*mu` and a
##' prior must be provided for coefficient of proportionality `kappa`
##'
##' * `lin` - sigma is a linear function of the mean, `sigma = a +
##' b*mu` and priors must be provided for `a` and `b`.
##'
##' The densities in each bin are assumed delta log Normal distributed
##' parameterized in term of parameters p - the probability of a
##' non-zero response, and m and s - the mean and standard deviation
##' of the log of the non-zero densities.  Three variant models are
##' provided for the p and s parameters
##'
##' * `const` - p and s are constant across bins
##'
##' * `fixed` - p and s vary across bins
##'
##' * `random` - log p and logit s across bins are modelled as Normal
##'   random effects.
##'
##' @title Delta LogNormal Mixture Model
##' @param Y a matrix of length bin densities.
##' @param breaks the length bin boundaries.
##' @param Ncomp the number of mixture components
##' @param prior a list of formula specifying priors
##' @param sigma the model for sigma
##' @param DLN the model for the DLN parameters
##' @return a `DMIX` object with components
##' * `call` - the model call
##' * `sigma` - the model for the sigma parameters
##' * `DLN` - the model for the delta log Normal parameters
##' * `jags.par` - the parameters to sample
##' * `jags.model` - the JAGS script
##' * `jags.data` - the JAGS data
##' @export
## ---- DLNMix
DLNMix <- function(Y,breaks,Ncomp,prior,
                   sigma = c("ind","cv","lin"),
                   DLN = c("fixed","const","random")) {

  ## Determine model for sigma and the delta log Normal
  sigma <- match.arg(sigma)
  DLN <- match.arg(DLN)

  ## Mixture model
  mix <- normalMixture(Ncomp,prior,sigma)
  ## Response model
  resp <- DLNResponse(DLN)

  model <- list(
    call = match.call(),
    sigma=sigma,
    DLN=DLN,
    jags.par=c(mix$par,resp$par),
    jags.model=paste("model {", mix$code, resp$code,"}",sep="\n"),
    jags.data=jagsData(Y,breaks,Ncomp))
  class(model) <- c("DLNMIX","DMIX")
  model
}
## ----


##' Print a DMIX object
##'
##' @title Print DMIX model objects
##' @param x a `DMIX` model object.
##' @param ... currently ignored.
##' @keywords internal
##' @export
## ---- print.DMIX
print.DMIX <- function(x,...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat(x$jags.model,sep="\n")
}
## ----



##' Call JAGS to sample from a DMIX model.
##'
##' This function calls JAGS to sample from a DMIX model generated.
##' `JAGSsample` communicates with JAGS through the rjags functions
##' [jags.model()][rjags::jags.model],
##' [update.jags()][rjags::update.jags],
##' [jags.samples()][rjags::jags.samples] and
##' [load.module()][rjags::load.module] and the manual pages for these
##' functions should be consulted for more details.
##'
##' @title Fit DMIX model
##' @param model a DMIX model object.
##' @param n.iter the number of iterations to run for sampling.
##' @param n.thin the thinning interval for sampling.
##' @param n.chains the number of chains to run.
##' @param n.adapt the number of iterations to run for adaptation.
##' @param n.update the number of iterations to run for initial update.
##' @param quiet whether to suppress JAGS reporting.
##' @return A list of `mcarray` objects.
##' @importFrom stats update
##' @importFrom rjags load.module jags.model jags.samples
##' @export
## ---- JAGSsample
JAGSsample <- function(model,n.iter=5000,n.thin=1,n.chains=4,
                       n.adapt=1000,n.update=1000,
                       quiet=!interactive()) {

  ## Ensure the glm module is loaded
  load.module("glm",quiet=quiet)
  pb <- if(quiet) "none" else getOption("jags.pb","text")
  ## Create the JAGS model
  jags <- jags.model(file=textConnection(model$jags.model),
                     data=model$jags.data,
                     inits=model$jags.init,
                     n.chains=n.chains,
                     n.adapt=n.adapt,
                     quiet=quiet)
  ## Burnin
  update(jags,n.iter=n.update,progress.bar=pb)
  ## Sample
  s <- jags.samples(jags,model$jags.par,n.iter=n.iter,thin=n.thin,progress.bar=pb)
  s
}
## ----



##' Convert a list of `mcarray` objects to one `mcmc.list` object.
##'
##' This function converts a list of `mcarray` objects to a single
##' `mcmc.list` object.  Note that individual `mcarray` objects may be
##' converted with `as.mcmc.list()`.  This function converts each
##' `mcarray` object to an `mcmc.list` object and then merges these
##' into a single `mcmc.list` object.
##'
##' @title Convert mcarray to mcmc.list
##' @param s a list of `mcarray` objects.
##' @return An `mcmc.list` object.
##' @seealso [coda::as.mcmc.list()]
##' @importFrom coda mcmc mcmc.list as.mcmc.list
##' @export
## ---- as.coda
as.coda <- function(s,name=TRUE) {

  ## Merge a list of mcmc.list objects
  merge <- function(...) {
    dots <- list(...)
    ## Extract sampling parameters
    start <- unique(sapply(dots,start))
    end <- unique(sapply(dots,end))
    thin <- unique(sapply(dots,thin))
    if(length(start)!=1L || length(end)!=1L || length(thin)!=1L)
      warning("All mcmc.list objects must have the same start, end and thin")
    ## Combine mcmc.list objects
    do.call(mcmc.list,
            lapply(do.call(mapply,c(cbind,dots,SIMPLIFY=FALSE)),
                   mcmc,start=start,end=end,thin=thin))
  }


  ## Convert each mcarray to an mcmc.list and merge them
  do.call(merge,lapply(s,as.mcmc.list))
}
## ----


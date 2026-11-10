## ---------------------------------------------------------------------------
## coloc 6.0.3, for the tests to compare against.
##
## Copied verbatim from commit 8f20f0b (coloc 6.0.3): the functions coloc()
## replaces and the ones of theirs that have changed since; the arithmetic they
## share with the new engine; and the dataset-processing and SuSiE-fitting
## functions, process.dataset() and runsusie() with what each calls, and
## check_alignment(). Never edit them: they are what "the same as 6.0.3" means.
##
## The copies call each other, not the package's versions. What they still
## take from outside is in the packages coloc imports: susieR's susie_rss(),
## data.table and stats. Both sides of a comparison use the installed susieR,
## so a change to susieR shows only in test-v603.R's tests pinned as numbers.
##
## get_match_6_0_3(), at the end, is how a test runs them.
## ---------------------------------------------------------------------------

v603 <- new.env()
local({

## R/claudia.R:9-11
Var.data <- function(f, N) {
    1 / (2 * N * f * (1 - f))
}

## R/claudia.R:21-23
Var.data.cc <- function(f, N, s) {
    1 / (2 * N * f * (1 - f) * s * (1 - s))
}

## R/claudia.R:33-37
logsum <- function(x) {
    my.max <- max(x)                              ##take out the maximum value in log form
    my.res <- my.max + log(sum(exp(x - my.max )))
    return(my.res)
}

## R/claudia.R:48-52
logdiff <- function(x,y) {
    my.max <- max(x,y)                              ##take out the maximum value in log form
    my.res <- my.max + log(max(exp(x - my.max) - exp(y-my.max), 0))
    return(my.res)
}

## R/claudia.R:67-84
approx.bf.p <- function(p,f,type, N, s, suffix=NULL) {
    if(type=="quant") {
        sd.prior <- 0.15
        V <- Var.data(f, N)
    } else {
        sd.prior <- 0.2
        V <- Var.data.cc(f, N, s)
    }
    z <- qnorm(0.5 * p, lower.tail = FALSE)
    ## Shrinkage factor: ratio of the prior variance to the total variance
    r <- sd.prior^2 / (sd.prior^2 + V)
    ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
    lABF = 0.5 * (log(1-r) + (r * z^2))
    ret <- data.frame(V,z,r,lABF)
    if(!is.null(suffix))
        colnames(ret) <- paste(colnames(ret), suffix, sep=".")
    return(ret)
}

## R/claudia.R:102-111
approx.bf.estimates <- function (z, V, type, suffix=NULL, sdY=1,
                                 effect_priors=c(quant=0.15,cc=0.2)) {
    sd.prior <- if (type == "quant") { effect_priors["quant"] * sdY } else { effect_priors["cc"] }
    r <- sd.prior^2/(sd.prior^2 + V)
    lABF = 0.5 * (log(1 - r) + (r * z^2))
    ret <- data.frame(V, z, r, lABF)
    if(!is.null(suffix))
        colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
    return(ret)
}

## R/claudia.R:123-141
combine.abf <- function(l1, l2, p1, p2, p12, quiet=FALSE) {
    stopifnot(length(l1)==length(l2))
    lsum <- l1 + l2
    lH0.abf <- 0
    lH1.abf <- log(p1) + logsum(l1)
    lH2.abf <- log(p2) + logsum(l2)
    lH3.abf <- log(p1) + log(p2) + logdiff(logsum(l1) + logsum(l2), logsum(lsum))
    lH4.abf <- log(p12) + logsum(lsum)

    all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
    my.denom.log.abf <- logsum(all.abf)
    pp.abf <- exp(all.abf - my.denom.log.abf)
    names(pp.abf) <- paste("PP.H", (1:length(pp.abf)) - 1, ".abf", sep = "")
    if(!quiet) {
        print(signif(pp.abf,3))
        print(paste("PP abf for shared variant: ", signif(pp.abf["PP.H4.abf"],3)*100 , '%', sep=''))
    }
    return(pp.abf)
}

## R/claudia.R:143-190
combine_abf_weighted <- function(l1, l2, p1, p2, p12,
                                prior_weights1, prior_weights2,
                                quiet = FALSE) {

  stopifnot(length(l1) == length(l2))

  Q <- length(l1)

  if (is.null(prior_weights1)) {
    p1_vec <- rep(p1, Q)
  } else {
    prior_weights1 <- prior_weights1 / sum(prior_weights1)
    p1_vec <- Q * p1 * prior_weights1
  }

  if (is.null(prior_weights2)) {
    p2_vec <- rep(p2, Q)
  } else {
    prior_weights2 <- prior_weights2 / sum(prior_weights2)
    p2_vec <- Q * p2 * prior_weights2
  }

  stopifnot(length(p1_vec) == length(p2_vec))
  stopifnot(length(p1_vec) == length(l1))

  p12_vec <- p1_vec * p2_vec * (p12 / (p1 * p2))

  lsum <- l1 + l2
  lH0_abf <- 0
  lH1_abf <- logsum(log(p1_vec) + l1)
  lH2_abf <- logsum(log(p2_vec) + l2)
  lH3_abf <- logdiff(logsum(log(p1_vec) + l1) + logsum(log(p2_vec) + l2),
                     logsum(log(p1_vec) + log(p2_vec) + lsum))
  lH4_abf <- logsum(log(p12_vec) + lsum)

  all_abf <- c(lH0_abf, lH1_abf, lH2_abf, lH3_abf, lH4_abf)
  denom_log_abf <- logsum(all_abf)
  pp_abf <- exp(all_abf - denom_log_abf)
  names(pp_abf) <- paste0("PP.H", (1:length(pp_abf)) - 1, ".abf")

  if(!quiet) {
    print(signif(pp_abf,3))
    print(paste("PP abf for shared variant: ", signif(pp_abf["PP.H4.abf"],3) * 100 ,
                "%", sep = ""))
  }

  pp_abf
}

## R/claudia.R:205-214
sdY.est <- function(vbeta, maf, n) {
    warning("estimating sdY from maf and varbeta, please directly supply sdY if known")
    oneover <- 1/vbeta
    nvx <- 2 * n * maf * (1-maf)
    m <- lm(nvx ~ oneover - 1)
    cf <- coef(m)[['oneover']]
    if(cf < 0)
        stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
    return(sqrt(cf))
}

## R/claudia.R:228-291
process.dataset <- function(d, suffix, ...) {
                                        #message('Processing dataset')

    nd <- names(d)
    ## if (! 'type' %in% nd)
    ##   stop("dataset ",suffix,": ",'The variable type must be set, otherwise the Bayes factors cannot be computed')

    ## if(!(d$type %in% c("quant","cc")))
    ##     stop("dataset ",suffix,": ","type must be quant or cc")

    ## if(d$type=="cc" & "pvalues" %in% nd) {
    ## if(!( "s" %in% nd))
    ##     stop("dataset ",suffix,": ","please give s, proportion of samples who are cases, if using p values")
    ## if(!("MAF" %in% nd))
    ##     stop("dataset ",suffix,": ","please give MAF if using p values")
    ## if(d$s<=0 || d$s>=1)
    ##     stop("dataset ",suffix,": ","s must be between 0 and 1")
    ## }

    ## if(d$type=="quant") {
    ##     if(!("sdY" %in% nd || ("MAF" %in% nd && "N" %in% nd )))
    ##         stop("dataset ",suffix,": ","must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated")
    ## }

    if("beta" %in% nd && "varbeta" %in% nd) {  ## use beta/varbeta.  sdY should be estimated by now for quant
        ## if(length(d$beta) != length(d$varbeta))
        ##   stop("dataset ",suffix,": ","Length of the beta vectors and variance vectors must match")
        ## if(!("snp" %in% nd))
        ##   d$snp <- sprintf("SNP.%s",1:length(d$beta))
        ## if(length(d$snp) != length(d$beta))
        ##   stop("dataset ",suffix,": ","Length of snp names and beta vectors must match")

        if(d$type=="quant" && !('sdY' %in% nd))
            d$sdY <- sdY.est(d$varbeta, d$MAF, d$N)
        df <- approx.bf.estimates(z=d$beta/sqrt(d$varbeta),
                                  V=d$varbeta, type=d$type, suffix=suffix, sdY=d$sdY)
        df$snp <- as.character(d$snp)
        if("position" %in% nd)
            df <- cbind(df,position=d$position)
        return(df)
    }

    if("pvalues" %in% nd & "MAF" %in% nd & "N" %in% nd) { ## no beta/varbeta: use p value / MAF approximation
        ## if (length(d$pvalues) != length(d$MAF))
        ##   stop('Length of the P-value vectors and MAF vector must match')
        ## if(!("snp" %in% nd))
        ##   d$snp <- sprintf("SNP.%s",1:length(d$pvalues))
        df <- data.frame(pvalues = d$pvalues,
                         MAF = d$MAF,
                         N=d$N,
                         snp=as.character(d$snp))
        snp.index <- which(colnames(df)=="snp")
        colnames(df)[-snp.index] <- paste(colnames(df)[-snp.index], suffix, sep=".")
        ## keep <- which(df$MAF>0 & df$pvalues > 0) # all p values and MAF > 0
        ## df <- df[keep,]
        abf <- approx.bf.p(p=df$pvalues, f=df$MAF, type=d$type, N=df$N, s=d$s, suffix=suffix)
        df <- cbind(df, abf)
        if("position" %in% nd)
            df <- cbind(df,position=d$position)
        return(df)
    }

    stop("Must give, as a minimum, one of:\n(beta, varbeta, type, sdY)\n(beta, varbeta, type, MAF)\n(pvalues, MAF, N, type)")
}

## R/claudia.R:356-363
adjust_prior=function(p,nsnps,suffix="") {
    if(nsnps * p >= 1) { ## for very large regions
        warning(paste0("p",suffix," * nsnps >= 1, setting p",suffix,"=1/(nsnps + 1)"))
        1/(nsnps + 1)
    } else {
        p
    }
}

## R/claudia.R:411-485
coloc.abf <- function(dataset1, dataset2, MAF=NULL,
                      p1=1e-4, p2=1e-4, p12=1e-5,
                      prior_weights1 = NULL, prior_weights2 = NULL, ...) {
    if(!("MAF" %in% names(dataset1)) & !is.null(MAF))
        dataset1$MAF <- MAF
    if(!("MAF" %in% names(dataset2)) & !is.null(MAF))
        dataset2$MAF <- MAF
    check_dataset(d=dataset1,1)
    check_dataset(d=dataset2,2)

    df1 <- process.dataset(d=dataset1, suffix="df1")
    df2 <- process.dataset(d=dataset2, suffix="df2")

    if (!is.null(prior_weights1) && (nrow(df1) != length(prior_weights1))) {
      stop("Length of prior_weights1 must match size of dataset 1")
    }
    if (!is.null(prior_weights2) && (nrow(df2) != length(prior_weights2))) {
      stop("Length of prior_weights2 must match size of dataset 2")
    }

    if (!is.null(prior_weights1) && any(is.na(prior_weights1) | (prior_weights1 < 0))) {
      stop("prior_weights1 must contain non-negative weights.")
    }

    if (!is.null(prior_weights2) && any(is.na(prior_weights2) | (prior_weights2 < 0))) {
      stop("prior_weights2 must contain non-negative weights.")
    }

    if (!is.null(prior_weights1) && !is.null(prior_weights2)) {
      warning("Prior weights specified for both traits.\n",
              "  The two weight vectors should be derived from independent sources of\n",
              "  information, i.e. not computed using the same method.")
    }

    p1=adjust_prior(p1,nrow(df1),"1")
    p2=adjust_prior(p2,nrow(df2),"2")

    merged.df <- merge(df1,df2)
    p12=adjust_prior(p12,nrow(merged.df),"12")

    prior_weights1 <- prior_weights1[which(df1$snp %in% merged.df$snp)]
    prior_weights2 <- prior_weights2[which(df2$snp %in% merged.df$snp)]

    if(!nrow(merged.df))
        stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")

    merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + lABF.df2)
    ## add SNP.PP.H4 - post prob that each SNP is THE causal variant for a shared signal
    my.denom.log.abf <- logsum(merged.df$internal.sum.lABF)
    merged.df$SNP.PP.H4 <- exp(merged.df$internal.sum.lABF - my.denom.log.abf)

    is_weighted <- !is.null(prior_weights1) || !is.null(prior_weights2)
    if (is_weighted) {
      pp.abf <- combine_abf_weighted(merged.df$lABF.df1, merged.df$lABF.df2,
                                     p1, p2, p12, prior_weights1, prior_weights2)
    } else {
      pp.abf <- combine.abf(merged.df$lABF.df1, merged.df$lABF.df2, p1, p2, p12)
    }

    common.snps <- nrow(merged.df)
    results <- c(nsnps=common.snps, pp.abf)
    output<-list(summary=results,
                 results=merged.df,
                 priors=c(p1=p1,p2=p2,p12=p12))

    if (is_weighted) {
      output$weights <- list(
        prior_weights1 = prior_weights1,
        prior_weights2 = prior_weights2
      )
    }

    class(output) <- c("coloc_abf",class(output))
    return(output)
}

## R/susie.R:16-41
logbf_to_pp=function(bf,pi, last_is_null) {
  n=if(last_is_null) {
      ncol(bf)-1 # number of snps, null is at n+1
    } else {
      n=ncol(bf) # number of snps
    }
  if(length(pi)==1) { # scalar pi
    if(pi > 1/n)
      pi=1/n
    pi=if(last_is_null) {
         c(rep(pi,n),1-n*pi)
       } else {
         rep(pi,n)
       }
  }
  if(any(pi == 0)) { # avoid division by zero
    pi[ pi==0 ] = 1e-16
    pi=pi/sum(pi)
  }
  if(last_is_null) {
    bf=bf - bf[,ncol(bf)] ## scale bf in case needed, so log BF for null = 0
  }
  priors=matrix(log(pi),nrow(bf),ncol(bf),byrow=TRUE)
  denom=matrix(apply(bf + priors,1,logsum), nrow(bf), ncol(bf))
  exp(bf + priors - denom)
}

## R/susie.R:70-100
coloc.susie=function(dataset1, dataset2, back_calculate_lbf=FALSE, susie.args=list(), ...) {
  ## if(!requireNamespace("susieR", quietly = TRUE)) {
  ##   message("please install susieR https://github.com/stephenslab/susieR")
  ##   return(NULL)
  ## }
  if("susie" %in% class(dataset1))
    s1=dataset1
  else
    s1=do.call("runsusie", c(list(d=dataset1,suffix=1),susie.args))
  if("susie" %in% class(dataset2))
    s2=dataset2
  else
    s2=do.call("runsusie", c(list(d=dataset2,suffix=2),susie.args))
  cs1=s1$sets
  cs2=s2$sets
  if(is.null(cs1$cs) || is.null(cs2$cs) || length(cs1$cs)==0 || length(cs2$cs)==0 ) {
	  warning("at least one dataset has no credible sets, nothing to colocalise")
    return(data.table(nsnps=NA))
	  }

  idx1=cs1$cs_index
  idx2=cs2$cs_index
  bf1=s1$lbf_variable[idx1,,drop=FALSE]
  bf2=s2$lbf_variable[idx2,,drop=FALSE]

  ret=coloc.bf_bf(bf1,bf2, ...)
  ## renumber index to match
  ret$summary[,idx1:=cs1$cs_index[idx1]]
  ret$summary[,idx2:=cs2$cs_index[idx2]]
  ret
}

## R/susie.R:132-148
coloc.susie_bf=function(dataset1,bf2, p1=1e-4, p2=1e-4, p12=5e-6, susie.args=list(), ...) {
  if("susie" %in% class(dataset1))
    s1=dataset1
  else
    s1=do.call("runsusie", c(list(d=dataset1,suffix=1),susie.args))
  cs1=s1$sets
  if(is.null(cs1$cs) || length(cs1$cs)==0 )
    return(data.table(nsnps=NA))
  idx1=cs1$cs_index
  ## alpha: an L by p matrix of posterior inclusion probabilites
  bf1=s1$lbf_variable[idx1,,drop=FALSE][,setdiff(colnames(s1$lbf_variable),""),drop=FALSE]

  ret=coloc.bf_bf(bf1,bf2, ...)
  ## renumber index to match
  ret$summary[,idx1:=cs1$cs_index[idx1]]
  ret
}

## R/susie.R:231-384
coloc.bf_bf=function(bf1,bf2, p1=1e-4, p2=1e-4, p12=5e-6, overlap.min=0.5,
                     trim_by_posterior=TRUE, prior_weights1 = NULL, prior_weights2 = NULL) {

  if ((length(p1) != 1 || length(p2) != 1 || length(p12) != 1) &&
      (!is.null(prior_weights1) || !is.null(prior_weights2))) {
      stop("Prior weights can only be used with single values of the prior parameters.")
  }

  if (!is.null(prior_weights1) && any(is.na(prior_weights1) | (prior_weights1 <= 0))) {
    stop("prior_weights1 must contain non-NA/non-negative weights.")
  }

  if (!is.null(prior_weights2) && any(is.na(prior_weights2) | (prior_weights2 <= 0))) {
    stop("prior_weights2 must contain non-NA/non-negative weights.")
  }
  
  if(is.vector(bf1))
    bf1=matrix(bf1,nrow=1,dimnames=list(NULL,names(bf1)))
  if(is.vector(bf2))
    bf2=matrix(bf2,nrow=1,dimnames=list(NULL,names(bf2)))
  
  if (!is.null(prior_weights1) && ncol(bf1) != length(prior_weights1)) {
    stop("Length of prior_weights1 must match size of dataset 1")
  }
  if (!is.null(prior_weights2) && ncol(bf2) != length(prior_weights2)) {
    stop("Length of prior_weights2 must match size of dataset 2")
  } 
  
  todo <- as.data.table(expand.grid(i=1:nrow(bf1),j=1:nrow(bf2)))
  todo[,pp4:=0]
  isnps=setdiff(intersect(colnames(bf1),colnames(bf2)),
                "null")
  if(!length(isnps))
    return(data.table(nsnps=NA))
    ## scale bf in case needed
    nulllast <- function(m) m[,c(setdiff(colnames(m),"null"), "null"), drop=FALSE]
  if("null" %in% colnames(bf1))
    bf1=nulllast(bf1) - matrix(bf1[,"null"],nrow(bf1),ncol(bf1))
  if("null" %in% colnames(bf2))
    bf2=nulllast(bf2) - matrix(bf2[,"null"],nrow(bf2),ncol(bf2))

  ## check whether isnps covers the signal for each trait
  pp1=logbf_to_pp(bf1,p1, last_is_null="null" %in% colnames(bf1))
  pp2=logbf_to_pp(bf2,p2, last_is_null="null" %in% colnames(bf2))
  ## pp1=logbf_to_pp(bf1,p1, last_is_null=TRUE)
  ## pp2=logbf_to_pp(bf2,p2, last_is_null=TRUE)
  ph0.1=if("null" %in% colnames(pp1)) { pp1[,"null"] } else { 1 - rowSums(pp1) }
  ph0.2=if("null" %in% colnames(pp2)) { pp2[,"null"] } else { 1 - rowSums(pp2) }
  prop1=rowSums(pp1[,c(isnps),drop=FALSE]) / rowSums(pp1[,setdiff(colnames(pp1),"null"),drop=FALSE])
  prop2=rowSums(pp2[,c(isnps),drop=FALSE]) / rowSums(pp2[,setdiff(colnames(pp2),"null"),drop=FALSE])
  if(trim_by_posterior==TRUE) {
    drop=sapply(1:nrow(todo), function(k) {
      prop1[todo$i[k]] < overlap.min | prop2[todo$j[k]] < overlap.min
    })
    if(all(drop)) {
      warning("snp overlap too small between datasets: too few snps with high posterior in one trait represented in other")
      return(
        list(summary=cbind(
               data.table(nsnps=length(isnps),
                          hit1=colnames(pp1)[apply(pp1,1,which.max)][todo$i],
                          hit2=colnames(pp2)[apply(pp2,1,which.max)][todo$j],
                          PP.H0.abf=pmin(ph0.1[todo$i],ph0.2[todo$j]),
                          PP.H1.abf=NA, PP.H2.abf=NA, PP.H3.abf=NA, PP.H4.abf=NA),
               todo[,.(idx1=i,idx2=j)]))
      )
    }
    if(any(drop))
      todo=todo[!drop,,drop=FALSE]
  }
  ## check p12
  if(length(p12)>1) { # only proceed if didn't need to trim snps
    if(length(isnps)!=ncol(bf1) || length(isnps)!=ncol(bf2))
      stop("p12 must have length 1 or equal length to p1 and p2, but different numbers of snps between datasets")
  }

  ## restrict bf1/2, p1, p2 and prior weights for incomplete snp overlap
  if(length(isnps)!=ncol(bf1)) {
    keep=match(isnps,colnames(bf1))
    prior_weights1 <- prior_weights1[keep]
    bf1=bf1[,keep,drop=FALSE]
    if(length(p1)>1)
      p1=p1[keep]
  }
  if(length(isnps)!=ncol(bf2)) {
    keep=match(isnps,colnames(bf2))
    prior_weights2 <- prior_weights2[keep]
    bf2=bf2[,keep,drop=FALSE]
    if(length(p2)>1)
      p2=p2[keep]
  }
  ## sort p12 if length(p1)>1 || length(p2)>1
  if(length(p12)==1) {
    if(length(p1)>1 || length(p2)>1) {
      p1_at_dist=min(p1)
      p2_at_dist=min(p2)
      c12=p12 / p1_at_dist / p2_at_dist
      p12=p1*p2*c12
    }
  }

  results <- PP <- vector("list",nrow(todo))
  ## results=lapply(1:nrow(todo), function(k) {
  for(k in 1:nrow(todo)) {
    df <- data.frame(snp=isnps, bf1=bf1[todo$i[k], isnps], bf2=bf2[todo$j[k], isnps])
    df$internal.sum.lABF <- with(df, bf1 + bf2)
    my.denom.log.abf <- logsum(df$internal.sum.lABF)
    df$SNP.PP.H4 <- exp(df$internal.sum.lABF - my.denom.log.abf)

    if (!is.null(prior_weights1) || !is.null(prior_weights2)) {
      pp.abf <- combine_abf_weighted(df$bf1, df$bf2,
                                     p1, p2, p12,
                                     prior_weights1,
                                     prior_weights2,
                                     quiet = TRUE)
    } else {
      pp.abf <- combine.abf(df$bf1, df$bf2, p1, p2, p12, quiet=TRUE)
    }

    PP[[k]]=df$SNP.PP.H4
    if(all(is.na(df$SNP.PP.H4))) {
      df$SNP.PP.H4=0
      pp.abf[1:5]=c(1,0,0,0,0)
    }
    common.snps <- nrow(df)
    hit1=names(which.max(bf1[todo$i[k],])) #df$snp[ which.max(abs(df$bf1)) ]
    if(is.null(hit1)) {
      hit1="-"
      pp.abf[c(1,3)]=c(0,1)
    }
    hit2=names(which.max(bf2[todo$j[k],])) #df$snp[ which.max(abs(df$bf2)) ]
    if(is.null(hit2)) {
      hit2="-"
      pp.abf[c(1,2)]=c(0,1)
    }
    results[[k]]=do.call("data.frame",c(list(nsnps=common.snps, hit1=hit1, hit2=hit2),
                                        as.list(pp.abf)))
  }
  results <- as.data.table(do.call("rbind",results))
  PP <- as.data.table(do.call("cbind",PP))
  if(nrow(todo)>1)
    colnames(PP)=paste0("SNP.PP.H4.row",1:nrow(todo))
  else
    colnames(PP)="SNP.PP.H4.abf"

  results=cbind(results,todo[,.(idx1=i,idx2=j)])
  PP=cbind(data.table(snp=isnps),PP)
  list(summary=results,
       results=PP,
       priors=if(length(p1)==1 && length(p2)==1 && length(p12)==1) {
                c(p1=p1,p2=p2,p12=p12)
              } else {
                list(p1=p1,p2=p2,p12=p12)
              })
}

## R/susie.R:421-479
runsusie=function(d,suffix=1,
                  maxit=100,repeat_until_convergence=TRUE,
                  s_init=NULL, ...) {
  check_dataset(d,suffix,req=c("beta","varbeta","LD","snp","N"))
  check_ld(d,d$LD)
  ##make copies of z and LD so we can subset if needed
  if(!("z" %in% names(d)))
    z=d$beta/sqrt(d$varbeta)
  else
    z=d$z
  LD=d$LD[d$snp,d$snp] # just in case
  names(z)=d$snp
  snp=d$snp

  converged=FALSE;
  ## set some defaults for susie arguments
  susie_args=list(...)
  if("max_iter" %in% names(susie_args)) {
    maxit=susie_args$max_iter
    susie_args = susie_args[ setdiff(names(susie_args), "max_iter") ]
  }
  ## at 0.12.6 susieR introduced need for n = sample size
  if(!("n" %in% names(susie_args)))
	  susie_args=c(list(n=d$N), susie_args)

  while(!converged) {
    message("running max iterations: ",maxit)
    res=do.call(susie_rss,
                c(list(z=z, R=LD, max_iter=maxit), susie_args))
    converged=res$converged; #s_init=res; maxit=maxit*2
    message("\tconverged: ",converged)
    if(!converged && repeat_until_convergence==FALSE)
      stop("susie_rss() did not converge in ",maxit," iterations. Try running with run_until_convergence=TRUE")
    if(!converged)
      maxit=maxit * 100 # no point in half measures!
  }
 res=annotate_susie(res, snp, LD)
  ## ## prune sets in high LD
  ## if(!is.null(r2.prune))
  ##   res=.susie_prune(res,r2.prune)

  ## if(!all(keep)) {
  ##   res$trimz=trimz
  ##   res$trimpeaks=trimpeaks
  ##   snps_to_add=setdiff(dbak$snp,snp)
  ##   res$pip=c(structure(rep(0,length(snps_to_add)), names=snps_to_add),
  ##             res$pip)
  ##   res$alpha=cbind(matrix(apply(res$alpha,1,min,na.rm=TRUE),
  ##                          nrow(res$alpha),length(snps_to_add),
  ##                          dimnames=list(rownames(res$alpha),snps_to_add)),
  ##                   res$alpha)
  ##   res$alpha=res$alpha/matrix(rowSums(res$alpha),nrow(res$alpha),ncol(res$alph))
  ##   res$lbf_variable=cbind(matrix(apply(res$lbf_variable,1,min,na.rm=TRUE),
  ##                                 nrow(res$lbf_variable),length(snps_to_add),
  ##                                 dimnames=list(rownames(res$lbf_variable),snps_to_add)),
  ##                          res$lbf_variable)
  ## }
  res
}

## R/susie.R:499-510
annotate_susie=function(res,snp, LD) {
    ## if(ncol(res$lbf_variable) != length(snp)+1)
    ##     stop("length of snp vector should be 1 less than ncol(res$lbf_variable)")
    colnames(res$lbf_variable) = c(snp,"null")[1:ncol(res$lbf_variable)]
    colnames(res$alpha) = c(snp,"null")[1:ncol(res$alpha)]
    names(res$pip)=snp
    if(length(res$sets$cs))
        res$sets$cs = lapply(res$sets$cs, function(x) { names(x) = snp[x]; x })
    res$sld=.susie_setld(res$sets$cs,LD)
    res$pruned=FALSE
    res
}

## R/susie.R:525-543
.susie_setld=function(s,ld) {
  ## stmp=lapply(s, setdiff, ncol(ld)+1)
  stmp=lapply(s, names)
  if(!length(stmp))
    return(0)
  if(length(stmp)==1)
    return(matrix(0,1,1,dimnames=list(names(stmp),names(stmp))))
  sld=matrix(0,length(stmp),length(stmp),dimnames=list(names(stmp),names(stmp)))
  for(i in 2:length(stmp)) {
    for(j in 1:(i-1)) {
      if(length(stmp[[i]]) && length(stmp[[j]])) {
        sld[i,j]=max(ld[stmp[[i]],stmp[[j]]]^2, na.rm=TRUE)
      } else {
        sld[i,j]=0
      }
    }
  }
  sld
}

## R/check.R:65-169
check_dataset <- function(d,
                          suffix="",
                          req=c("type","snp"),
                          warn.minp=1e-6) {
    if(!is.list(d) )
        stop("dataset ",suffix,": is not a list")
    recognised_items=c("beta","varbeta","pvalues","MAF","snp","position","N","type","s","sdY","LD")
    nd <- intersect(names(d), recognised_items)

    ## no missing values - make people clean their own data rather
    ## than make assumptions here for datasets I don't know
    n <- 0
    
    if(length(setdiff(req,nd)))
        stop("dataset ",suffix,": missing required element(s) ",paste(setdiff(req,nd),collapse=", "))

    for(v in nd) {
        if(any(is.na(d[[v]])))
            stop("dataset ",suffix,": ",v," contains missing values")
    }

    if (!(d$type %in% c("quant", "cc")))
        stop("dataset ", suffix, ": ", "type must be quant or cc")

    ## snps should be unique
    if("snp" %in% nd && any(duplicated(d$snp)))
        stop("dataset ",suffix,": duplicated snps found")
    if("snp" %in% nd && is.factor(d$snp))
        stop("dataset ",suffix,": snp should be a character vector but is a factor")

    ## MAF should be > 0, < 1
    if("MAF" %in% nd && (!is.numeric(d$MAF) || any(is.na(d$MAF)) ||
                         any(d$MAF<=0) || any(d$MAF>=1)))
        stop("dataset ",suffix,": MAF should be a numeric, strictly >0 & <1")

    ## lengths of vector arguments should match
    l <- -1 # impossible length
    shouldmatch <- intersect(nd, c("pvalues","MAF","beta","varbeta","snp","position"))
    for(v in shouldmatch)
        if(l<0) { ## update
            l <- length(d[[v]])
        } else { ## check
            if(length(d[[v]])!=l) {
                stop("dataset ",suffix,": lengths of inputs don't match: ")
                print(shouldmatch)
            }
        }
    ## N should be scalar if present
    if("N" %in% nd && (!is.numeric(d[["N"]]) || length(d[["N"]])!=1))
       stop("N should be a single numeric. coloc does not accomodate per-snp sample sizes.")

    ## type of data
    if (! ('type' %in% nd))
        stop("dataset ",suffix,": variable type not set")
    if(!(d$type %in% c("quant","cc")))
        stop("dataset ",suffix,": ","type must be quant or cc")

    ## no beta/varbeta
    if(("s" %in% nd) && (!is.numeric(d$s) || d$s<=0 || d$s>=1))
        stop("dataset ",suffix,": ","s must be between 0 and 1")
    if(!("beta" %in% nd) || !("varbeta" %in% nd)) { # need to estimate var (Y)
        if(!("pvalues" %in% nd) || !( "MAF" %in% nd))
            stop("dataset ",suffix,": ","require p values and MAF if beta, varbeta are unavailable")
        if(any(d$pvalues<=0))
            stop("pvalues should not be negative or exactly 0")
        if(d$type=="cc" && !("s" %in% nd))
            stop("dataset ",suffix,": ","require, s, proportion of samples who are cases, if beta, varbeta are unavailable")
        if (!('N' %in% nd) || is.null(d$N) || any(d$N<=0) )
            stop("dataset ",suffix,": sample size N <=0 or not set")
        p=d$pvalues
    } else {
        p=pnorm( -abs( d$beta/sqrt(d$varbeta) ) ) * 2
    }

    # Infinite values.
    if (any(is.infinite(d$beta) | is.infinite(d$varbeta))) {
      stop("dataset ", suffix,": Infinite values in beta and/or varbeta")
    }

    # zero variance.
    if (any(d$varbeta==0)) {
      stop("dataset ", suffix,": Zero values in varbeta")
    }

    ## minp
    if(min(p) > warn.minp)
        warning("minimum p value is: ",format.pval(min(p)),"\nIf this is what you expected, this is not a problem.\nIf this is not as small as you expected, please check you supplied var(beta) and not sd(beta) for the varbeta argument. If that's not the explanation, please check the 02_data vignette.")

    ## sdY
    if(d$type=="quant" && !("sdY" %in% nd))
        if(!("MAF" %in% nd && "N" %in% nd ))
            stop("dataset ",suffix,": ","must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated")

    if("LD" %in% nd) {
        if(nrow(d$LD)!=ncol(d$LD))
            stop("LD not square")
        if(!identical(colnames(d$LD),rownames(d$LD)))
            stop("LD rownames != colnames")
        if(length(setdiff(d$snp,colnames(d$LD))))
            stop("colnames in LD do not contain all SNPs")
    }

    ## if we reach here, no badness detected
    NULL
}

## R/check.R:180-193
check_ld <- function(D,LD) {
    if(is.null(LD))
      stop("LD required")
    if(!is.matrix(LD))
      stop("LD must be of class matrix")
    ## if(any(LD[upper.tri(LD)]==1))
    ##   stop("LD includes SNPs in perfect LD")
    if(nrow(LD)!=ncol(LD))
        stop("LD not square")
    if(!identical(colnames(LD),rownames(LD)))
      stop("LD rownames != colnames")
    if(length(setdiff(D$snp,colnames(LD))))
        stop("colnames in LD do not contain all SNPs")
}

## R/check.R:204-218
check_alignment <- function(D,thr=0.2,do_plot=TRUE) {
  check_dataset(D)
  bprod=outer(D$beta/sqrt(D$varbeta),D$beta/sqrt(D$varbeta),"*")
  ## plot(bprod,D$LD[D$snp,D$snp],xlab="product of z scores",ylab="LD")
  tmp=(bprod/D$LD)[abs(D$LD) > 0.2]
  if(do_plot) {
  hist(tmp,
       xlab="ratio of product of Z scores to LD",
       main="alignment check plot",
       sub="expect most values to be positive\nsymmetry is a warning sign of potentially poor alignment")
  legend("topright",legend=paste0("% positive = ",100*round(mean( tmp > 0 ),3)))
  abline(v=0,col="red")
  }
  return(mean(tmp > 0))
}

}, envir = v603)

## What 6.0.3 returns for a call, written exactly as it was made in 6.0.3:
## get_match_6_0_3(coloc.abf(D1, D2)) runs 6.0.3's coloc.abf() on the test's
## D1 and D2. Warnings and errors come through; what 6.0.3 prints does not.
get_match_6_0_3 <- function(call) {
    value <- NULL
    ## all.names, or a copy named with a leading dot, such as .susie_setld(),
    ## would be missed and the package's version run instead
    utils::capture.output(
        value <- eval(substitute(call), as.list(v603, all.names = TRUE), parent.frame()))
    value
}

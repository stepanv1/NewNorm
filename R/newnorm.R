###  TODO there should be a default Annot matrix with required information for all probes.
###  ###  sigA and sigB have had log transforamtions prior this should be incoporated (log (1+x))

## if you use GenomeStudio to extract control probes, then they would be in a certain order
## that we could use

newnorm <- function(sigA, sigB, Annot=default.Annot, quantiledat=NULL,
                    controlred, controlgrn, cp.types, cell_type, ncmp=3,
                    save.quant=TRUE, save.loess=TRUE, apply.loess=FALSE, logit.quant=FALSE)
{
    # functions
    logitfn <- function(x) { log(x/(1-x)) }
    ##
    mult1 <- function(x, vec1) return(x*vec1)
    ##
    sum1 <- function(x, v2) { return(x + v2) }
    ##
    extractqnt <- function(x, i, AB, ncmp)  { 
        return(x$fitted.values[i,AB,ncmp])  }
    ##
    ## just the usual lm functgion
    lmfunction1<- function(yvec, xmat)  {
        if (sum(is.na(yvec))==0)  {
            fit1 <- lm(yvec ~ cell_type/Disease + Sex + age_quant, data = xmat)
            cfit1 <- summary(fit1)$coefficients
            fit2 <- anova(fit1)
            outp <- c(cfit1[2:nrow(cfit1),4],fit2[4,5])
            names(outp) <- c(rownames(cfit1[2:nrow(cfit1)]),rownames(fit2)[4])
            return(outp)   
        }}
    ##
    qqfn <- function(vect1)  {
        nnonmiss <- sum(!is.na(vect1))
        xaxs <- sort(-log10( (1:nnonmiss - 0.5) / nnonmiss))
        plot(c(0, max(xaxs)), c(0, max(-log10(vect1), na.rm=TRUE)), type='n', 
             xlab='Expected', ylab='Observed')
        abline(0,1)
        points(xaxs, sort(-log10(vect1[!is.na(vect1)])), pch=16, cex=0.4)
    }
    ##
    qqfn2 <- function(vect1, vect2)  {
        nnonmiss <- sum(!is.na(vect1))
        xaxs <- sort(-log10( (1:nnonmiss - 0.5) / nnonmiss))
        plot(c(0, max(xaxs)), c(0, max(-log10(vect1), na.rm=TRUE)), type='n', 
             xlab='Expected', ylab='Observed')
        abline(0,1)
        points(xaxs, sort(-log10(vect1[!is.na(vect1)])), pch=16, cex=0.4, col=1)
        nnonmiss <- sum(!is.na(vect2))
        xaxs <- sort(-log10( (1:nnonmiss - 0.5) / nnonmiss))
        points(xaxs, sort(-log10(vect2[!is.na(vect2)])), pch=16, cex=0.4, col=2)
    }
    
    NotNumeric<-function(a){
            return(any(!is.finite(as.matrix(a)))) 
        }
    
    #checking sanity of the data
    print("Checking sanity of the data...")
    if (NotNumeric(sigA)){stop("There are non-numeric values in the matrix", '\n')}
    if (NotNumeric(sigB)){stop("There are non-numeric values in the matrix", '\n')}
    if (NotNumeric(controlred)){stop("There are non-numeric values in the matrix", '\n')}
    if (NotNumeric(controlgrn)){stop("There are non-numeric values in the matrix", '\n')}
    if (any(cell_type == '' | typeof(cell_type) != 'character' | is.na(cell_type))) {{stop("There are non-character values in cell_type", '\n')}}
    if (any(cp.types == '' | typeof(cp.types) != 'character' | is.na(cp.types))) {{stop("There are non-character values in cp.types", '\n')}}
    
    #check that ID's match between sigA, sigB and control probe data
    if (!(identical(colnames(controlgrn), colnames(sigA)) &
              identical(colnames(controlred), colnames(sigA)) &
              identical(colnames(sigB), colnames(sigA)))) {
        if ((ncol(controlgrn) == ncol(sigA)) & (ncol(controlred) == ncol(sigA)) 
            & (ncol(controlred) == ncol(sigB)) & (length(cp.types) == ncol(sigB)) & (length(cell_type) == ncol(sigB))){
            ix <- sort(colnames(controlred), index.return = T)$ix
            controlgrn<-controlgrn[,ix]; sigA<-sigA[,ix]; sigB<-sigB[,ix];
            if (!(identical(colnames(controlgrn), colnames(sigA)) &
                      identical(colnames(controlred), colnames(sigA)) &
                      identical(colnames(sigB), colnames(sigA)))){
                stop("Sample names do not match.", '\n')}
            
        } else{stop("Data dimensions or samples names do not match.", '\n')}
    }
    
    ### TODO add a check that all probes in sigA and sigB also exist in Annot (Annot can be bigger but cannot be smaller)
    ### then extract relevant rows from Annot
    
    print("Data is ok.")
    
    
    nr <- nrow(sigA)
    qntllist <- c((0.5)/nr, seq(0.001, 0.009, 0.001), seq(0.01,0.05,0.01), 
                  seq(0.07,0.93,by=0.02), seq(0.95,0.99,0.01), seq(0.991,0.999,0.001),
                  (nr-0.5)/nr)
    nqnt <- length(qntllist)  # number of desired quantiles
    
    ####! TODO good to add special treatment of X and Y chromosomes.  See Fortin et al. particularly for Y
    if (save.quant)  {
        quantilesA.red <- matrix(NA, ncol(sigA), nqnt)   
        quantilesA.grn <- quantilesA.red;  quantilesA.II <- quantilesA.red
        quantilesB.red <- matrix(NA, ncol(sigA), nqnt)  
        quantilesB.grn <- quantilesA.red; quantilesB.II <- quantilesA.red
        #  need separate estimation for X  but the X chromosome is currently removed
        #quantiles.Xred <- quantiles.red
        #quantiles.Xgrn <- quantiles.red
        #quantiles.XII <- quantiles.red
        for (j in (1:nqnt))  {
            quantilesA.red[,j] <- apply(sigA[Annot$Type=="I" & Annot$Color=='Red',], 
                                        2, quantile, qntllist[j])  
            quantilesA.grn[,j] <- apply(sigA[Annot$Type=="I" & Annot$Color=='Grn',], 
                                        2, quantile, qntllist[j])  
            quantilesA.II[,j] <- apply(sigA[Annot$Type=="II",], 
                                       2, quantile, qntllist[j])  
            quantilesB.red[,j] <- apply(sigB[Annot$Type=="I" & Annot$Color=='Red',], 
                                        2, quantile, qntllist[j])  
            quantilesB.grn[,j] <- apply(sigB[Annot$Type=="I" & Annot$Color=='Grn',], 
                                        2, quantile, qntllist[j])  
            quantilesB.II[,j] <- apply(sigB[Annot$Type=="II",], 
                                       2, quantile, qntllist[j])  
            #  quantiles.Xred[,j] <- apply(JPNorm_ABeta[Annot$Type=="I" & Annot$Color=='Red' &
            #               Annot$chr=="X",], 2, quantile, j/nqnt)  
            #  quantiles.Xgrn[,j] <- apply(JPNorm_ABeta[Annot$Type=="I" & Annot$Color=='Grn' &
            #               Annot$chr=="X",], 2, quantile, j/nqnt)  
            #  quantiles.XII[,j] <- apply(JPNorm_ABeta[Annot$Type=="II" & Annot$chr=="X",], 
            #                             2, quantile, j/nqnt)  
        } # for j
        save(quantilesA.red, file="quantilesA.red.RData")
        save(quantilesA.grn, file="quantilesA.grn.RData")
        save(quantilesA.II, file="quantilesA.II.RData")
        save(quantilesB.red, file="quantilesB.red.RData")
        save(quantilesB.grn, file="quantilesB.grn.RData")
        save(quantilesB.II, file="quantilesB.II.RData")
    }  # if save.quant
    if (!save.quant)  {
        load("quantilesA.red.RData")
        load("quantilesA.grn.RData")
        load("quantilesA.II.RData")
        load("quantilesB.red.RData")
        load("quantilesB.grn.RData")
        load("quantilesB.II.RData")
    }
    
    # TODO get cp.types from control probe data automatically
    # do log(x+1) of the control probes prior to code below
    
    # assume log transformation has already been done
    # construct control probe summaries, averages by type of control probe
    # then specifically create columns by cell type
    cp.type.tab <- table(cp.types)
    for (k in (1:length(cp.type.tab))) {
        temp1 <- apply(controlred[cp.types==names(cp.type.tab[k]),],2,mean) 
        temp2 <- apply(controlgrn[cp.types==names(cp.type.tab[k]),],2,mean)
        if (k==1) {
            mat.by.ct <- cbind(temp1,temp2)  }
        if (k>1) {  mat.by.ct <- cbind(mat.by.ct, cbind(temp1,temp2))  }
    } # end for   mat.by.ct is 15 columns - means for each type of control probe
    
    ind=matrix(FALSE,ncol=length(cell_type),nrow=length(unique(cell_type)))
    for(j in 1:length(unique(cell_type))){
        ind[j,]<- (cell_type==names(table(cell_type))[j])
    }
    # now add additional sets of 15 columns, each specific to one cell type
    ctl.covmat <- cbind(mat.by.ct, apply(mat.by.ct, 2, mult1, ind[1,]))
    for(j in 2:length(unique(cell_type))){
        ctl.covmat <- cbind(ctl.covmat, apply(mat.by.ct, 2, mult1, ind[j,]))
    }
    
    
    if (save.quant)  {
        save(ctl.covmat, file="ctl.covmat.RData")
    }
    if (!save.quant)  {
    load("ctl.covmat.RData")
    }
    
    # here are the actual fits, one for each quantile in qntllist
    fit2.red <- list() ; fit2.grn <- list();  fit2.II <- list()
    for (j in (1:nqnt))  {
        fit2.red <- c(fit2.red, list(
            plsr(cbind(quantilesA.red[,j], quantilesB.red[,j]) ~ ctl.covmat, ncomp=ncmp)))
        fit2.grn <- c(fit2.grn, list(
            plsr(cbind(quantilesA.grn[,j], quantilesB.grn[,j]) ~ ctl.covmat, ncomp=ncmp)))
        fit2.II <- c(fit2.II, list( 
            plsr(cbind(quantilesA.II[,j], quantilesB.II[,j]) ~ ctl.covmat, ncomp=ncmp)))
    }
    
    # smooth the predictions so can get a value for any quantile
    loessfitsA.red <- list()
    loessfitsA.grn <- list()
    loessfitsA.II <- list()
    loessfitsB.red <- list()
    loessfitsB.grn <- list()
    loessfitsB.II <- list()
    if (!logit.quant) {xq <- qntllist }
    if (logit.quant)  {xq <- logitfn(qntllist)  }
    for (i in (1:nrow(quantilesA.red))) {    # i: over samples
        # assemble the fitted values into a vector
        fittedvals <-  unlist(lapply(fit2.red, extractqnt, i=i, AB=1, ncmp=ncmp))
        coef.fit <- loess(fittedvals ~ xq, span = 0.25)
        loessfitsA.red <- c(loessfitsA.red, list(coef.fit))  
        fittedvals <-  unlist(lapply(fit2.grn, extractqnt, i=i, AB=1, ncmp=ncmp))
        coef.fit <- loess(fittedvals ~ xq, span = 0.25)
        loessfitsA.grn <- c(loessfitsA.grn, list(coef.fit))  
        fittedvals <-  unlist(lapply(fit2.II, extractqnt, i=i, AB=1, ncmp=ncmp))
        coef.fit <- loess(fittedvals ~ xq, span = 0.25)
        loessfitsA.II <- c(loessfitsA.II, list(coef.fit))  
        
        fittedvals <-  unlist(lapply(fit2.red, extractqnt, i=i, AB=2, ncmp=ncmp))
        coef.fit <- loess(fittedvals ~ xq, span = 0.25)
        loessfitsB.red <- c(loessfitsB.red, list(coef.fit))  
        fittedvals <-  unlist(lapply(fit2.grn, extractqnt, i=i, AB=2, ncmp=ncmp))
        coef.fit <- loess(fittedvals ~ xq, span = 0.25)
        loessfitsB.grn <- c(loessfitsB.grn, list(coef.fit))  
        fittedvals <-  unlist(lapply(fit2.II, extractqnt, i=i, AB=2, ncmp=ncmp))
        coef.fit <- loess(fittedvals ~ xq, span = 0.25)
        loessfitsB.II <- c(loessfitsB.II, list(coef.fit))  
    }
    if (save.loess)  {
        save(loessfitsA.red, file="loessfitsA.red")
        save(loessfitsA.grn, file="loessfitsA.grn")
        save(loessfitsA.II, file="loessfitsA.II")
        save(loessfitsB.red, file="loessfitsB.red")
        save(loessfitsB.grn, file="loessfitsB.grn")
        save(loessfitsB.II, file="loessfitsB.II")
    }
    loessfits <- list(loessfitsA.red, loessfitsA.grn, loessfitsA.II,
                      loessfitsB.red, loessfitsB.grn, loessfitsB.II)
    if (apply.loess)  {
        betamatrices <- applynewnorm(loessfits, sigA, sigB, Annot)
        return(betamatrices)
    } else {   return(loessfits)  }
}  # end function



### TODO add in a new function that does cross-validation with different numbers of components and plots results.  
### I have code to help with this.

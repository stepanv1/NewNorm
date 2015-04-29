
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
#####################################################
#  ifquant=TRUE quantiles will be provided in quantiledat as a list
#  ifquant=FALSE, calculate quantiles from sigA and sigB and Annot
#  controldata has red, green, pval, for each of the control probes in consecutive cols
#  column order must match in sigA, sigB, and control data
newnorm <- function(ifquant, sigA, sigB=NULL, Annot=NULL, quantiledat=NULL,
                    controlred, controlgrn, cp.types, cell_type, ncmp=3,
                    save.loess=TRUE, applyloess=FALSE, logit.quant=FALSE)
{
    nr <- nrow(sigA)
    qntllist <- c((0.5)/nr, seq(0.001, 0.009, 0.001), seq(0.01,0.05,0.01), 
                  seq(0.07,0.93,by=0.02), seq(0.95,0.99,0.01), seq(0.991,0.999,0.001),
                  (nr-0.5)/nr)
    nqnt <- length(qntllist)  # number of desired quantiles
    
    if (!(ifquant))  {
        # (1) quantiles of data
        # Could include more quantiles especially at lower and upper ends
        quantilesA.red <- matrix(NA, ncol(sigA), nqnt)   # nqnt + 1 is important
        quantilesA.grn <- quantilesA.red;  quantilesA.II <- quantilesA.red
        quantilesB.red <- matrix(NA, ncol(sigA), nqnt)   # nqnt + 1 is important
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
    }  # if ifquant
    if (ifquant)  {
        load("quantilesA.red.RData")
        load("quantilesA.grn.RData")
        load("quantilesA.II.RData")
        load("quantilesB.red.RData")
        load("quantilesB.grn.RData")
        load("quantilesB.II.RData")
    }
    
    # assume log transformation has already bee done
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
    # now add 3 additional sets of 15 columns, each specific to one cell type
    ctl.covmat <- cbind(mat.by.ct, apply(mat.by.ct, 2, mult1, ind[1,]))
    for(j in 2:length(unique(cell_type))){
        ctl.covmat <- cbind(ctl.covmat, apply(mat.by.ct, 2, mult1, ind[j,]))
    }
    
    # here are the actual fits
    fit2.red <- plsr(cbind(quantilesA.red, quantilesB.red) ~ ctl.covmat, ncomp=ncmp)
    fit2.grn <- plsr(cbind(quantilesA.grn, quantilesB.grn) ~ ctl.covmat, ncomp=ncmp)
    fit2.II <- plsr(cbind(quantilesA.II, quantilesB.II) ~ ctl.covmat, ncomp=ncmp)
    
    # smooth the predictions so can get a value for any quantile
    if (save.loess)  {
        loessfitsA.red <- list()
        loessfitsA.grn <- list()
        loessfitsA.II <- list()
        loessfitsB.red <- list()
        loessfitsB.grn <- list()
        loessfitsB.II <- list()
        if (!logit.quant) {xq <- qntllist }
        if (logit.quant)  {xq <- logitfn(qntllist)  }
        for (i in (1:dim(fit2.red$fitted.values)[1]))  {
            coef.fit <- loess(fit2.red$fitted.values[i,1:nqnt,ncmp] ~ xq, span = 0.2)
            loessfitsA.red <- c(loessfitsA.red, list(coef.fit))  
            coef.fit <- loess(fit2.grn$fitted.values[i,1:nqnt,ncmp] ~ xq, span = 0.2)
            loessfitsA.grn <- c(loessfitsA.grn, list(coef.fit))  
            coef.fit <- loess(fit2.II$fitted.values[i,1:nqnt,ncmp] ~ xq, span = 0.2)
            loessfitsA.II <- c(loessfitsA.II, list(coef.fit))  
            coef.fit <- loess(fit2.red$fitted.values[i,(nqnt+1):(2*nqnt),ncmp] ~ xq, span = 0.2)
            loessfitsB.red <- c(loessfitsB.red, list(coef.fit))  
            coef.fit <- loess(fit2.grn$fitted.values[i,(nqnt+1):(2*nqnt),ncmp] ~ xq, span = 0.2)
            loessfitsB.grn <- c(loessfitsB.grn, list(coef.fit))  
            coef.fit <- loess(fit2.II$fitted.values[i,(nqnt+1):(2*nqnt),ncmp] ~ xq, span = 0.2)
            loessfitsB.II <- c(loessfitsB.II, list(coef.fit))  
        }
        save(loessfitsA.red, file="loessfitsA.red.RData")
        save(loessfitsA.grn, file="loessfitsA.grn.RData")
        save(loessfitsA.II, file="loessfitsA.II.RData")
        save(loessfitsB.red, file="loessfitsB.red.RData")
        save(loessfitsB.grn, file="loessfitsB.grn.RData")
        save(loessfitsB.II, file="loessfitsB.II.RData")
    }  # if save.loess
    if (!save.loess)  {
        load(file="loessfitsA.red.RData")
        load(file="loessfitsA.grn.RData")
        load(file="loessfitsA.II.RData")
        load(file="loessfitsB.red.RData")
        load(file="loessfitsB.grn.RData")
        load(file="loessfitsB.II.RData")
    }
    loessfits <- list(loessfitsA.red, loessfitsA.grn, loessfitsA.II,
                      loessfitsB.red, loessfitsB.grn, loessfitsB.II)
    if (applyloess)  {
        betamatrices <- applynewnorm(loessfits, sigA, sigB, Annot)
        return(betamatrices)
    } else {   return(loessfits)  }
}  # end function
###############################################################
# separate pls fits for each quantile
newnorm2 <- function(ifquant, sigA, sigB=NULL, Annot=NULL, quantiledat=NULL,
                     controlred, controlgrn, cp.types, cell_type, ncmp=3,
                     save.loess=TRUE, applyloess=FALSE, logit.quant=FALSE)
{
    nr <- nrow(sigA)
    qntllist <- c((0.5)/nr, seq(0.001, 0.009, 0.001), seq(0.01,0.05,0.01), 
                  seq(0.07,0.93,by=0.02), seq(0.95,0.99,0.01), seq(0.991,0.999,0.001),
                  (nr-0.5)/nr)
    nqnt <- length(qntllist)  # number of desired quantiles
    
    if (ifquant==FALSE)  {
        # (1) quantiles of data
        # Could include more quantiles especially at lower and upper ends
        quantilesA.red <- matrix(NA, ncol(sigA), nqnt)   # nqnt + 1 is important
        quantilesA.grn <- quantilesA.red;  quantilesA.II <- quantilesA.red
        quantilesB.red <- matrix(NA, ncol(sigA), nqnt)   # nqnt + 1 is important
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
    }  # if ifquant
    if (ifquant==TRUE)  {
        load("quantilesA.red.RData")
        load("quantilesA.grn.RData")
        load("quantilesA.II.RData")
        load("quantilesB.red.RData")
        load("quantilesB.grn.RData")
        load("quantilesB.II.RData")
    }
    
    # assume log transformation has already bee done
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
    
    ind1 <- (cell_type=="MONO"); ind1[is.na(ind1) ] <- FALSE
    ind2 <- (cell_type=="TCELL"); ind2[is.na(ind2) ] <- FALSE
    ind3 <- (cell_type=="BCELL"); ind3[is.na(ind3) ] <- TRUE
    # now add 3 additional sets of 15 columns, each specific to one cell type
    ctl.covmat <- cbind(mat.by.ct, apply(mat.by.ct, 2, mult1, ind1))
    ctl.covmat <- cbind(ctl.covmat, apply(mat.by.ct, 2, mult1, ind2))
    ctl.covmat <- cbind(ctl.covmat, apply(mat.by.ct, 2, mult1, ind3))
    
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
    if (applyloess)  {
        betamatrices <- applynewnorm(loessfits, sigA, sigB, Annot)
        return(betamatrices)
    } else {   return(loessfits)  }
}  # end function

######################
# now obtain predictions when all the control covariate data are the same 
# by cell type)
#   ctl.covmat.mono <- apply(ctl.covmat[ind1,],  2, mean)
#   ctl.covmat.tcell <- apply(ctl.covmat[ind2,], 2, mean)
#   ctl.covmat.bcell <- apply(ctl.covmat[ind3,], 2, mean)
#   xmono <- matrix(ctl.covmat.mono, nrow=1, ncol=length(ctl.covmat.mono), 
#		    byrow=TRUE);  colnames(xmono) <- colnames(ctl.covmat)
#   xtcell <- matrix(ctl.covmat.tcell, nrow=1, ncol=length(ctl.covmat.mono), 
#		    byrow=TRUE);  colnames(xtcell) <- colnames(ctl.covmat)
#   xbcell <- matrix(ctl.covmat.bcell, nrow=1, ncol=length(ctl.covmat.mono), 
#		    byrow=TRUE);  colnames(xbcell) <- colnames(ctl.covmat)

# predictions for the  mean value of the control probes, by cell type
#   preds.red.mono <- predict(fit2.red, newdata = xmono, ncomp=ncmp)
#   preds.grn.mono <- predict(fit2.grn, newdata = xmono, ncomp=ncmp)
#   pred.II.mono <- predict(fit2.II, newdata = xmono, ncomp=ncmp)
#   preds.red.tcell <- predict(fit2.red, newdata = xtcell, ncomp=ncmp)
#   preds.grn.tcell <- predict(fit2.grn, newdata = xtcell, ncomp=ncmp)
#   pred.II.tcell <- predict(fit2.II, newdata = xtcell, ncomp=ncmp)
#   preds.red.bcell <- predict(fit2.red, newdata = xbcell, ncomp=ncmp)
#   preds.grn.bcell <- predict(fit2.grn, newdata = xbcell, ncomp=ncmp)
#   pred.II.bcell <- predict(fit2.II, newdata = xbcell, ncomp=ncmp)
##############################3
applynewnorm <- function(loessfits, sigA, sigB, Annot)  {
    
    wh.red <- which(Annot$Type=='I' & Annot$Color=="Red")
    wh.grn <- which(Annot$Type=='I' & Annot$Color=="Grn")
    wh.II <- which(Annot$Type=='II')
    
    rankmatA.red <- (apply(sigA[wh.red,],2,  rank) - 0.5)/length(wh.red)
    rankmatA.grn <- (apply(sigA[wh.grn,],2,  rank) - 0.5)/length(wh.grn)
    rankmatA.II <- (apply(sigA[wh.II,],2,  rank) - 0.5)/length(wh.II)
    
    rankmatB.red <- (apply(sigB[wh.red,],2,  rank) - 0.5)/length(wh.red)
    rankmatB.grn <- (apply(sigB[wh.grn,],2,  rank) - 0.5)/length(wh.grn)
    rankmatB.II <- (apply(sigB[wh.II,],2,  rank) - 0.5)/length(wh.II)
    
    predmatA.red <- matrix(NA, nrow=nrow(rankmatA.red), ncol = ncol(rankmatA.red))
    predmatA.grn <- matrix(NA, nrow=nrow(rankmatA.grn), ncol = ncol(rankmatA.grn))
    predmatA.II <- matrix(NA, nrow=nrow(rankmatA.II), ncol = ncol(rankmatA.II))
    predmatB.red <- matrix(NA, nrow=nrow(rankmatB.red), ncol = ncol(rankmatB.red))
    predmatB.grn <- matrix(NA, nrow=nrow(rankmatB.grn), ncol = ncol(rankmatB.grn))
    predmatB.II <- matrix(NA, nrow=nrow(rankmatB.II), ncol = ncol(rankmatB.II))
    for (i in (1:ncol(predmatA.red)))  {
        predmatA.red[,i] <- predict(loessfits[[1]][[i]], newdata = logitfn(rankmatA.red[,i]))
        predmatA.grn[,i] <- predict(loessfits[[2]][[i]], newdata = logitfn(rankmatA.grn[,i]))
        predmatA.II[,i] <- predict(loessfits[[3]][[i]], newdata = logitfn(rankmatA.II[,i]))
        predmatB.red[,i] <- predict(loessfits[[4]][[i]], newdata = logitfn(rankmatB.red[,i]))
        predmatB.grn[,i] <- predict(loessfits[[5]][[i]], newdata = logitfn(rankmatB.grn[,i]))
        predmatB.II[,i] <- predict(loessfits[[6]][[i]], newdata = logitfn(rankmatB.II[,i]))
    }
    predmatA <- matrix(NA, nrow(sigA), ncol(predmatA.red))
    predmatB <- predmatA
    predmatA[wh.red,] <- predmatA.red
    predmatA[wh.grn,] <- predmatA.grn
    predmatA[wh.II,] <- predmatA.II
    predmatB[wh.red,] <- predmatB.red
    predmatB[wh.grn,] <- predmatB.grn
    predmatB[wh.II,] <- predmatB.II
    
    newBeta <- (exp(predmatB)-1)/(exp(predmatA) + exp(predmatB)-2)
    
    #rm(predmatA); rm(predmatB); rm(predmatA.red); rm(predmatA.grn); rm(predmatA.II)
    #rm(predmatB.red); rm(predmatB.grn); rm(predmatB.II)
    
    # now for residuals
    #   residmat <- (exp(sigB)-1)/(exp(sigA) + exp(sigB) - 2 + 0.1)  - predBeta
    origBeta <- (exp(sigB)-1)/(exp(sigA) + exp(sigB) - 2  + 0.1)
    origBeta <- (origBeta * 999 + 0.5)/1000 
    origBetamn <- apply(origBeta,1,mean)
    
    return(list(origBeta, newBeta))
}
###############################

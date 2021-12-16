## Required Libraries to run these functions:
#library(survival)
#library(tidyverse)
#library(ggplot2)
#library(patchwork)
#library(GSVA)
#library(survMisc)
#library(dplyr)

#' provide cutoff threshold and gene expression values for specific gene
#'
#' @param gene the gene we want to use to classify samples based on expression
#' @param mat tibble with samples and gene expression (FPKM)
#' @param os tibble with samples, overall survival (OS) values, and Vital.Status values
#' @param cut method to use to separate into high and low groups, can be either median (defualt) or cutP
#' @return the cutoff threshold and expression values (log2(FPKM+1)) as a list
#' @examples
#' data = data.frame(data = matrix(rnorm(100),ncol=10))
#' expr = dplyr::tibble(symbol = paste0(rep("gene",10),1:10), data)
#' colnames(expr) = c("symbol", letters[1:10])
#' os = dplyr::tibble(sample = letters[1:10], OS = 10:19, Vital.Status = c(rep("Alive",8),rep("Dead",2)))
#' mCut("gene1", mat = expr, os = os, cut = "cutP")
mCut = function(gene, mat, os, cut = "median") {
    geneExp = mat %>% dplyr::filter(symbol == gene) %>% dplyr::select(!"symbol") %>% unlist(., use.names = FALSE) %>% as.numeric() # extract all values for gene    
    geneExp = log2(geneExp + 1) # retrieve log2(FPKM+1)
    if (cut == "median") {
        cutoff = median(geneExp[2:length(geneExp)]) # find median for gene (removing NA at beginning due to symbol column)
    } else if (cut == "cutP") {
        os["Exp"] = geneExp
        cox.os = suppressWarnings(survival::coxph(survival::Surv(time=os$OS, event=os$Vital.Status=="Dead")~Exp, data = os)) # get survival based on Expression
        c = survMisc::cutp(cox.os)$Exp # get cutp output
        data.table::setorder(c, "Exp") # create table with cut points
        percentile <- ecdf(c$Exp) # find cdf of Expression points
        cutoff <- as.numeric(c[order(c$Q), ][nrow(c), 1]) # find optimal expression CutOff
        print(paste0("Cutoff = ", signif(cutoff,3)))
    } else {
        print("invalid cut method, please insert either median or cutP ... returning os without separation")
        os$group = "Low"
        return(os)
    }
    out = list(signif(cutoff,3), geneExp)
    names(out) = c("cutoff","geneExp")
    return(out)
}

# function mSet 
    # inputs:
    # gs -- gene set we want to use to classify samples based on median expression
    # mat -- matrix with rownames as Ensembl IDs and colnames as sample names
    # os -- tibble with samples, OS values, empty column "group"
    # ref -- reference tibble to use with helper function
    # cut -- method to use to separate into high and low (either median (default) or cutP)
    # outputs:
    # os -- cutoff value and GSVA scores in list
mSet = function(gs, mat, os, ref, cut = "median") {
    gs = findENS(gs, ref) # convert to ensmbl ids
    mat = mat[rownames(mat) %in% unlist(gs),] # extract only genes of interest
    mat = log2(mat + 1) # retrieve log2(FPKM+1)
    scores = callGSVA(mat,gs)
    #scores = suppressWarnings(gsva(mat, gs, verbose=FALSE, method="zscore", parallel.sz=2)) # calculate ssGSEA enrichment scores
    if (cut == "median") {
        cutoff = median(scores) # find median score for gene set
        #highSamples = colnames(mat)[scores > med]
    } else if (cut == "cutP") {
        os["Scores"] = 0
        # add scores into os
        for (i in 1:length(scores)) {
            #os[os$sample == colnames(scores)[i],"Scores"] = scores[i]
            os[os$sample == rownames(scores)[i],"Scores"] = scores[i]
        }
        cox.os = coxph(Surv(time=os$OS, event=os$Vital.Status=="Dead")~Scores, data = os) # get survival based on GSEA scores
        c = cutp(cox.os)$Scores # get cutp output
        data.table::setorder(c, Scores) # create table with cut points
        percentile <- ecdf(c$Scores) # find cdf of GSEA score points
        cutoff <- as.numeric(c[order(c$Q), ][nrow(c), 1]) # find optimal GSEA score CutOff
        #highSamples = colnames(mat[,-1])[scores > signif(low,3)] # divide samples
        #lowSamples = colnames(mat[,-1])[scores <= signif(low,3)] # divide samples
        print(paste0("Cutoff = ", signif(cutoff,3)))
    } else {
        print("invalid cut method, please insert either median or cutP ... returning os without separation")
        os$group = "Low"
        return(os)
    }
    #os$group[os$sample %in% lowSamples] = "Low"
    #os$group[os$sample %in% highSamples] = "High"
    #os$group = factor(os$group, levels = c("Low","High"))# change levels so that the reference group is Low expression
    #return(os)
    out = list(signif(cutoff,3),scores)
    names(out) = c("cutoff","scores")
    return(out) # return cutoff point and Scores
}

# function findENS
    # inputs:
    # genes -- gene names we want to find ENSMBL ids for
    # ref -- reference matrix to search
    # outputs:
    # ens -- list of ENSMBL ids for the input genes
findENS = function(genes, ref) {
    ens = "" # create empty vector to add onto
    for (i in 1:length(genes)) {
        ind = which(ref$symbol == genes[i]) # find index or indices of gene symbol
        if (length(ind) > 1) {
            ens = c(ens, ref$ensg[ind[1]]) # take first index
        } else if (length(ind) == 1) {
            ens = c(ens, ref$ensg[ind])
        } else {
            print(paste(genes[i], "not found in ENSEMBL reference."))
        }
    }
    ens = ens[2:length(ens)] # remove "" at beginning
    return(ens)
}

# function plotSurv -- not working
    # inputs: 
    # gene -- gene we want to use to split expression, and then calculate survival from
    # mat -- tibble with samples and gene expression
    # os -- tibble with samples, OS values, empty column "group"
    # subtype -- type of MPAL we are analyzing (BM or TM)
    # outputs: 
    # save plot in plots folder - GENENAMEsurvivalSUBTYPE.pdf
    # if gene not found, print "gene not found in dataset"
plotSurv = function(gene, mat, ENS, os, subtype, subfolder, ref) {
    if (length(gene) > 1) {
        os = suppressWarnings(mSet(gene, ENS, os, ref)) # calculate median, classify as high/low expression
    } else {
        os = suppressWarnings(mCut(gene, mat, os)) # calculate median, classify as high/low expression
    }
    # check if gene is not found in dataset
    if (length(unique(os$group)) == 1) {
        print("Gene not found in dataset!")
    } else {
        pdf(file = paste0("plots/", subfolder, "/",gene,"survival",subtype,".pdf"))
        plot(survfit(Surv(time=os$OS, event=os$Vital.Status=="Dead")~group, data = os), 
            main = paste0("Plot of Survival Curves by ",gene," Expression in ",subtype," MPAL"), 
            xlab = "Length of Survival (days)",ylab="Probability of Survival",col=c("red","blue"))
        legend("topright", legend=c("High","Low"),fill=c("red","blue"),bty="n")
        dev.off()
    }
}

# function survStats
    # inputs:
    # os -- tibble with samples, OS values, and filled in column "group"
    # outputs:
    # stats -- tibble filled in with test statistics (Comparison, survdiffP, coxHR, coxP)
survStats = function(os) {
    stats = data.frame("Comparison" = "", "survdiffP" = "", "coxHR" = "", "coxP" = "")
    surv = Surv(time=os$OS, event=os$Vital.Status=="Dead")
    sdOut = survdiff(surv~group, data = os) # save survdiff output
    cpOut = coxph(surv~group, data = os) # save coxph output
    stats$survdiffP = 1 - pchisq(sdOut$chisq, length(sdOut$n) - 1) # fill in survdiff p-value
    stats$coxHR = summary(cpOut)$coefficients[2] # fill in coxph Hazard Ratio
    stats$coxP = summary(cpOut)$coefficients[5] # fill in coxph p-value
    return(stats)
}

# function checkSurv
checkSurv = function(gene, ENS, mat, os, ref) {
    if (length(gene) > 1) {
        os = mSet(gs, ENS, os, ref)
    } else {
        os = mCut(gene, mat, os)
    }
    
    check = data.frame("OS" = os$OS, "group" = os$group, "geneExp" = 0)
    for (i in 1:nrow(check)) {
        sample_i = os$sample[i] # extract sample id from os
        check$geneExp[i] = unlist(mat[mat$symbol == gene, colnames(mat) == sample_i]) # add gene exp for gene/sample into check
    }
    print(check %>% group_by(group) %>% summarize(meanExp = mean(geneExp), meanOS = mean(OS)))

    return(check)
}

# FROM SURVIVALGENIE
#'@name callGSVA
#'@aliases callGSVA
#'@title GSVA enrichment analysis
#'@description Estimates GSVA enrichment zscores.
#'@usage callGSVA(x,y)
#'@param x A data frame or matrix of gene or probe expression values where rows corrospond to genes and columns corrospond to samples
#'@param y A list of genes as data frame or vector
#'@details This function uses "zscore" gene-set enrichment method in the estimation of gene-set enrichment scores per sample.
#'@return A gene-set by sample matrix of GSVA enrichment zscores.
#'@import GSVA
#'@examples 
#'g <- 10 ## number of genes
#'s <- 30 ## number of samples
#'## sample data matrix with values ranging from 1 to 10
#'rnames <- paste("g", 1:g, sep="")
#'cnames <- paste("s", 1:s, sep="")
#'expr <- matrix(sample.int(10, size = g*s, replace = TRUE), nrow=g, ncol=s, dimnames=list(rnames, cnames))
#'## genes of interest
#'genes <- paste("g", 1:g, sep="")
#'## Estimates GSVA enrichment zscores.
#'callGSVA(expr,genes)
#'@seealso GSVA
#'@export
callGSVA = function(x,y) {
    if(missing(x)){
    stop("input expression data missing!")
    }
    if(missing(y)){
    stop("input gene set missing!")
    }
    #genes <- list(set1=y)
    genes = list(y)
    gsva.results <- gsva(x, genes, method="zscore",verbose=FALSE, parallel.sz=2)
    tr_gsva.results <- t(gsva.results)
    colnames(tr_gsva.results) <- c("GSVA score")
    return (tr_gsva.results)
}

# function cut 
    # inputs:
    # os -- tibble with samples, OS values
    # out -- (output from mCut or mSet), list with cutoff value and either GSVA scores or Gene Expression values
    # outputs:
    # os -- input tibble now with group column separated into High and Low values based on cutoff and scores/expression
cut = function(os, out) {
    if (names(out)[2] == "geneExp") {
        os["Exp"] = out[[2]]
        os[os$Exp > out[[1]],"group"] = "High"
        os[os$Exp <= out[[1]],"group"] = "Low"
    } else {
        os["Scores"] = out[[2]]
        os[os$Scores > out[[1]],"group"] = "High"
        os[os$Scores <= out[[1]],"group"] = "Low"
    }
    os$group = factor(os$group, levels = c("Low","High"))# change levels so that the reference group is Low expression
    return(os$group)
}
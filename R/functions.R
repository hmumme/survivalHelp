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
#' @param mat expression values (FPKM) tibble with samples as columns and symbol column containing gene names
#' @param os tibble with samples, overall survival (OS) values, and Vital.Status values
#' @param cut method to use to separate into high and low groups, can be either median (defualt) or cutP
#' @return the cutoff threshold and expression values (log2(FPKM+1)) as a list
#' @export
#' @examples
#' data = abs(data.frame(data = matrix(rnorm(100),ncol=10)))
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

#' provide cutoff threshold and GSEA values for gene set of two or more genes
#'
#' @param gs the genes we want to use to classify samples based on expression
#' @param ENS expression values (FPKM) tibble with samples as columns and ENSMBL ids as rows
#' @param os tibble with samples, overall survival (OS) values, and Vital.Status values
#' @param ref reference tibble with ensmble ids and gene symbols as columns, converts to ENS using this
#' @param cut method to use to separate into high and low groups, can be either median (defualt) or cutP
#' @return the cutoff threshold and GSEA scores as a list
#' @export
#' @examples
#' expr = abs(matrix(rnorm(100),ncol=10))
#' rownames(expr) = paste0(rep("ENSG",10),1:10)
#' colnames(expr) = letters[1:10]
#' ref = dplyr::tibble(ensg = paste0(rep("ENSG",10), 1:10), symbol = paste0(rep("gene",10), 10:1))
#' os = dplyr::tibble(sample = letters[1:10], OS = 10:19, Vital.Status = c(rep("Alive",8),rep("Dead",2)))
#' mSet(c("gene1","gene4","gene5"), ENS = expr, os = os, ref = ref, cut = "cutP")
mSet = function(gs, ENS, os, ref, cut = "median") {
    gs = findENS(gs, ref) # convert to ensmbl ids
    ENS = ENS[rownames(ENS) %in% unlist(gs),] # extract only genes of interest
    ENS = log2(ENS + 1) # retrieve log2(FPKM+1)
    scores = callGSVA(ENS,gs)
    if (cut == "median") {
        cutoff = median(scores) # find median score for gene set
    } else if (cut == "cutP") {
        os["Scores"] = 0
        # add scores into os
        for (i in 1:length(scores)) {
            os[os$sample == rownames(scores)[i],"Scores"] = scores[i]
        }
        cox.os = suppressWarnings(survival::coxph(survival::Surv(time=os$OS, event=os$Vital.Status=="Dead")~Scores, data = os)) # get survival based on GSEA scores
        c = survMisc::cutp(cox.os)$Scores # get cutp output
        data.table::setorder(c, Scores) # create table with cut points
        percentile <- ecdf(c$Scores) # find cdf of GSEA score points
        cutoff <- as.numeric(c[order(c$Q), ][nrow(c), 1]) # find optimal GSEA score CutOff
        print(paste0("Cutoff = ", signif(cutoff,3)))
    } else {
        print("invalid cut method, please insert either median or cutP ... returning os without separation")
        os$group = "Low"
        return(os)
    }
    out = list(signif(cutoff,3),scores)
    names(out) = c("cutoff","scores")
    return(out) # return cutoff point and Scores
}

#' convert gene symbols to ENSMBL ids
#'
#' @param genes gene symbols
#' @param ref reference tibble with ensmble ids and gene symbols as columns
#' @return character vector of ENSMBL ids for genes
#' @examples
#' genes = c("gene3","gene4","gene10")
#' ref = dplyr::tibble(ensg = paste0(rep("ENSG",10), 1:10), symbol = paste0(rep("gene",10), 10:1))
#' findENS(genes, ref)
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

# function survStats
    # inputs:
    # os -- tibble with samples, OS values, and filled in column "group"
    # outputs:
    # stats -- tibble filled in with test statistics (Comparison, survdiffP, coxHR, coxP)
#' calculate survival statistics for split survival data
#' 
#' @param os tibble with sample, OS, Vital.Status, and group (High/Low) columns
#' @return stats tibble filled with test statistics for survival data: Comparison, survdiffP, coxHR, coxP
#' @export
survStats = function(os) {
    stats = data.frame("Comparison" = "", "survdiffP" = "", "coxHR" = "", "coxP" = "")
    surv = survival::Surv(time=os$OS, event=os$Vital.Status=="Dead")
    sdOut = survival::survdiff(surv~group, data = os) # save survdiff output
    cpOut = survival::coxph(surv~group, data = os) # save coxph output
    stats$survdiffP = 1 - pchisq(sdOut$chisq, length(sdOut$n) - 1) # fill in survdiff p-value
    stats$coxHR = summary(cpOut)$coefficients[2] # fill in coxph Hazard Ratio
    stats$coxP = summary(cpOut)$coefficients[5] # fill in coxph p-value
    return(stats)
}

# FROM SURVIVALGENIE
#'@name callGSVA
#'@aliases callGSVA
#'@title GSVA enrichment analysis
#'@description Estimates GSVA enrichment zscores (from SurvivalGenie).
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
#' split survival data into High and Low groups based on cutoff value and scores/expression values
#'
#' @param os tibble with columns: samples, OS values, and Vital.Status
#' @param out output from mCut or mSet. It contains the cutoff threshold and either GSEA scores ("scores") or Gene Expression ("geneExp") values for each sample
#' @return the os["group"] column with samples split into High and Low groups
#' @export
#' @examples
#' os = dplyr::tibble(sample = letters[1:10], OS = 10:19, Vital.Status = c(rep("Alive",8),rep("Dead",2)))
#' out = list("cutoff" = 0.44, "scores" = rnorm(10))
#' os["group"] = cut(os,out)
#' os
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
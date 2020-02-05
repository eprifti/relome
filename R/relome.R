################################################################
#  _____       _                   ____            _           #
# |_   _|     | |                 / __ \          (_)          #
#   | |  _ __ | |_ ___  __ _ _ __| |  | |_ __ ___  _  ___ ___  #
#   | | | '_ \| __/ _ \/ _` | '__| |  | | '_ ` _ \| |/ __/ __| #
#   | |_| | | | ||  __| (_| | |  | |__| | | | | | | | (__\__ \ #
# |_____|_| |_|\__\___|\__, |_|   \____/|_| |_| |_|_|\___|___/ #
#                       __/ |                                  #
#                      |___/                                   #
################################################################

################################################################
# @script: relome.R
# @note:explore relations between variables
# @author: Edi Prifti
# @date: August 2016
################################################################



# This function plots tests, p and q values for relations using the phenoPairwiseRelations object
plotRelations <- function(data, ppr, rel.list, var.interest, threshold=0.05, 
                          verbose=FALSE, 
                          bg=list(scatter="gold2",
                                   category=gray.colors(2)),
                          col=list(scatter="black",
                                  category=gray.colors(2)),
                          pch=21
                          )
{
  y <- data[,var.interest]
  if(is.numeric(y)){
    for(i in 1:length(rel.list)){
      if (verbose) print(paste(i,names(rel.list)[i]))
      x <- data[,names(rel.list)[i]]
      p <- extractSignificant(ppr$p, interest = names(rel.list)[i], threshold = threshold)[[1]][var.interest]
      q <- extractSignificant(ppr$q, interest = names(rel.list)[i], threshold = threshold)[[1]][var.interest]
      test <- ppr$test[var.interest,names(rel.list)[i]]
      if(is.numeric(x)){
        plot(y~x, xlab=names(rel.list)[i], ylab=var.interest,
             bg=bg$scatter, col=col$scatter, pch=pch, main=paste(names(rel.list)[i],
                                                "\np =", signif(p,digits=2),
                                                "; q =", signif(q,digits=2),
                                                "\ntest =", test))
        abline(lm(y~x),col="red")
      }else{
        boxplot(y~x, xlab=names(rel.list)[i], ylab=var.interest, notch=T, col=gray.colors(length(levels(x))),
                pch=pch, main=paste(names(rel.list)[i],
                                   "\np =", signif(p,digits=2),
                                   "; q =", signif(q,digits=2),
                                   "\ntest =", test))
      }
    }
  }else if(is.factor(y)){ # if a factor
    for(i in 1:length(rel.list)){
      x <- data[,names(rel.list)[i]]
      p <- extractSignificant(ppr$p, interest = names(rel.list)[i], threshold = threshold)[[1]][var.interest]
      q <- extractSignificant(ppr$q, interest = names(rel.list)[i], threshold = threshold)[[1]][var.interest]
      test <- ppr$test[var.interest,names(rel.list)[i]]
      if(is.numeric(x)){
        plot(x~y, ylab=names(rel.list)[i], xlab=var.interest, notch=T, col=col$category,
             pch=pch, main=paste(names(rel.list)[i],
                                "\np =", signif(p,digits=2),
                                "; q =", signif(q,digits=2),
                                "\ntest =", test))
      }else{
        plot(y~x, xlab=names(rel.list)[i], ylab=var.interest, col=col$category, pch=pch, main=paste(names(rel.list)[i],
                                                                                                   "\np =", signif(p,digits=2),
                                                                                                   "; q =", signif(q,digits=2),
                                                                                                   "\ntest =", test))
      }
    }
  }else{ warning("Unknown type for the variable of interest!")}
}


plotPhenoRelations <- function(data, 
                               ppr, 
                               rel.list, 
                               var.interest, 
                               threshold=0.05, 
                               verbose=FALSE, 
                               bg=list(scatter="gold2",
                                        category=gray.colors(2)),
                               col=list(scatter="black",
                                        category=gray.colors(2)),
                               pch=21
                               )
{
  y <- data[,var.interest]
  if(is.numeric(y))
  {
    for(i in 1:length(rel.list))
    {
      if (verbose) print(paste(i,names(rel.list)[i]))
      x <- data[,names(rel.list)[i]]
      # p <- extractSignificant(ppr$p, interest = names(rel.list)[i], threshold = threshold)[[1]][var.interest]
      # q <- extractSignificant(ppr$q, interest = names(rel.list)[i], threshold = threshold)[[1]][var.interest]
      # get the test results
      p <- ppr$p[var.interest, names(rel.list)[i]]
      q <- ppr$q[var.interest, names(rel.list)[i]]
      test <- ppr$test[var.interest, names(rel.list)[i]]
      
      if(is.numeric(x))
      {
        plot(y~x, xlab=names(rel.list)[i], ylab=var.interest,
             bg=bg$scatter, col=col$scatter, pch=pch, main=paste(names(rel.list)[i],
                                                "\np =", signif(p,digits=2),
                                                "; q =", signif(q,digits=2),
                                                "\ntest =", test))
        abline(lm(y~x), col="red")
      }else
      {
        boxplot(y~x, xlab=names(rel.list)[i], ylab=var.interest, notch=T, col=gray.colors(length(levels(x))),
                pch=pch, main=paste(names(rel.list)[i],
                                   "\np =", signif(p,digits=2),
                                   "; q =", signif(q,digits=2),
                                   "\ntest =", test))
      }
    }
  }else if(is.factor(y))
  { # if a factor
    for(i in 1:length(rel.list))
    {
      x <- data[,names(rel.list)[i]]
      # get the test results
      p <- ppr$p[var.interest, names(rel.list)[i]]
      q <- ppr$q[var.interest, names(rel.list)[i]]
      test <- ppr$test[var.interest, names(rel.list)[i]]
      
      if(is.numeric(x))
      {
        plot(x~y, ylab=names(rel.list)[i], xlab=var.interest, notch=T, col=col$category, 
             pch=pch, main=paste(names(rel.list)[i],
                                "\np =", signif(p,digits=2),
                                "; q =", signif(q,digits=2),
                                "\ntest =", test))
      }else
      {
        plot(y~x, xlab=names(rel.list)[i], ylab=var.interest, col=col$category, 
             pch=pch, main=paste(names(rel.list)[i],
                                "\np =", signif(p,digits=2),
                                "; q =", signif(q,digits=2),
                                "\ntest =", test))
      }
    }
  }else{ warning("Unknown type for the variable of interest!")}
}

# modified multiple test adjustment line by line instead of a whole 20151128 EP
relome <- function (data, adjust = "BH", adjust.by.var = TRUE, verbose = FALSE) {
  # create an empty square matrix
  data.relations <- matrix(NA, ncol(data), ncol(data))
  colnames(data.relations) <- colnames(data)
  rownames(data.relations) <- colnames(data)
  data.relations.test <- data.relations
  
  for (i in 1:ncol(data)) {
    if (verbose) print(paste(i, colnames(data)[i])) # log
    for (j in 1:ncol(data)) {
      if (verbose) cat("\t",paste(j, colnames(data)[j]),";\t") # log
      x <- data[, i]
      y <- data[, j]
      
      
      if(is.numeric(x))
      {
        if(var(x, na.rm = TRUE)==0)
        {
          if (verbose) print("This test was not performed due to low variance")
          data.relations[i, j] <- NA
          data.relations.test[i, j] <- NA
          next
        }
      }
      
      if(is.numeric(y))
      {
        if(var(y, na.rm = TRUE)==0)
        {
          if (verbose) print("This test was not performed due to low variance")
          data.relations[i, j] <- NA
          data.relations.test[i, j] <- NA
          next
        }
      }
      
      if(is.factor(x))
      {
        if(length(levels(x))<=1)
        {
          if (verbose) print("This test was not performed due to low number of levels")
          data.relations[i, j] <- NA
          data.relations.test[i, j] <- NA
          next
        }
      }
      
      if(is.factor(y))
      {
        if(length(levels(y))<=1)
        {
          if (verbose) print("This test was not performed due to low number of levels")
          data.relations[i, j] <- NA
          data.relations.test[i, j] <- NA
          next
        }
      }    
        
      res <- all.test(x,y)
      if (verbose) cat(res$test,"\n") # log
      data.relations[i, j] <- res$p
      data.relations.test[i, j] <- res$test
      
    }
  }
  if(adjust.by.var) # if we are testing only one variable with the rest we can adjust only for one variable for multiple testing
  {
    data.relations.q <- apply(data.relations, 2, p.adjust, method=adjust) # this is on the columns. 
    # q-vqlues here tend to be smaller (more significant) than in the else when all tests are considered
  }else # otherwise, we make a long vector with all the tests and run adjustment before transforming into a matrix again
  {
    # adjust all p-values
    tmp <- p.adjust(as.vector(data.relations), method = adjust)
    data.relations.q <- matrix(tmp, nrow = nrow(data.relations), byrow = FALSE)
    colnames(data.relations.q) <- colnames(data.relations)
    rownames(data.relations.q) <- rownames(data.relations)
  }
  
  return(list(p = data.relations, q = data.relations.q, test=data.relations.test))
}


# modified multiple test adjustment line by line instead of a whole 20151128 EP
all.test <- function (x, y, adjust = "BH", verbose = FALSE) {
  p <- NA
  test <- "null"
  if (is.factor(x)) {
    if (is.factor(y)) { # if both variables are factors
      if (nrow(table(x, y)) == 1 | ncol(table(x, y)) == 1 | sum(table(x, y)) == 0 | sum(rowSums(table(x,y))!=0)<2 | sum(rowSums(table(y,x))!=0)<2) {
        warning("This is particular case 1")
      }else {
        if (verbose) print("chisq.test") # log
        #p <- chisq.test(table(x, y))$p.value
        
        p <- chisq.test(x, y)$p.value
        test <- "chisq.test"
      }
    } else { # if y is a numerical
      if(length(na.omit(y))>4){ # sanity check
        if (nortest::lillie.test(y)$p.value > 0.05) { # if non normal distribution
          if (length(table(x)) > 1) {
            if (length(table(x)) == 2) {
              if (all(colSums(table(y, x)) > 1)) {
                if (verbose) print("t.test") # log
                p <- t.test(y ~ x)$p.value
                test <- "t.test"
              }
            }else { # if more than two categories use a linear model
              if (sum(table(y, x)) > 1 & length(unique(x[!is.na(x)])) !=  1 & sum(colSums(table(y, x)) > 0) > 1) {
                if (verbose) print("linear model") # log
                tmp <- lmp(lm(y ~ x))
                if (!is.null(tmp)) {
                  p <- tmp
                  test <- "linear model"
                } # end is null result from linear model
              }
            }
          }
        } # end nortest
        else { # if non normal
          if (length(table(x)) > 1) {
            if (length(table(x)) == 2) {
              if (all(colSums(table(y, x)) > 1)) {
                if (verbose) print("wilcox.test") # log
                p <- wilcox.test(y ~ x)$p.value
                test <- "wilcox.test"
              }
            }else {
              if (sum(table(y, x)) > 1 & length(unique(x[!is.na(x)])) != 1 & sum(colSums(table(y, x)) > 0) > 1) {
                if (verbose) print("kruskal.test") # log
                tmp <- kruskal.test(y ~ x)
                if (!is.null(tmp)) {
                  p <- tmp$p.value
                  test <- "kruskal.test"
                } # end is null result from linear model
              } # end validity of the test
            } # end else
          } # end at least two categories
        } # end else normality
      } # sanity check
    } # end else numerical
  } # end x factor
  else { # If x is numerical
    if (is.factor(y)) { # and y is a factor
      if(length(na.omit(x))>4){ # sanity check
        if (nortest::lillie.test(x)$p.value > 0.05) { # if normal
          if (length(table(y)) > 1) { # if at least two factors
            if (length(table(y)) == 2) { # if exactly 2
              if (all(colSums(table(x, y)) > 1)) {
                if (verbose) print("t.test") # log
                p <- t.test(x ~ y)$p.value
                test <- "t.test"
              }
            }
            else {
              if (sum(table(x, y)) > 1 & length(unique(y[!is.na(y)])) != 1 & sum(colSums(table(x, y)) > 0) > 1) {
                if (verbose) print("linear model") # log
                tmp <- lmp(lm(x ~ y))
                if (!is.null(tmp)) {
                  p <- tmp
                  test <- "lm"
                }
              }
            }
          }
        } # end normality
        else { # if non normal
          if (length(table(y)) > 1) {
            if (length(table(y)) == 2) { # if only two categories
              if (all(colSums(table(x, y)) > 1)) {
                if (verbose) print("wilcox.test") # log
                p <- wilcox.test(x ~ y)$p.value
                test <- "wilcox.test"
              }
            } else {
              if (sum(table(x, y)) > 1 & length(unique(y[!is.na(y)])) != 1 & sum(colSums(table(x, y)) > 0) > 1) {
                if (verbose) print("kruskal.test") # log
                tmp <- kruskal.test(x ~ y)
                if (!is.null(tmp)) {
                  p <- tmp$p.value
                  test <- "kruskal.test"
                } # end is null
              } # end validity kurskal
            } 
          } # end of testing at least one category
        } # end else normality
      } # end sanity check
    } else { # if both variables are numeric
      if (length(x[!is.na(x)]) > 4 & length(x[!is.na(y)]) > 4) {
        if (nortest::lillie.test(x)$p.value < 0.05 | nortest::lillie.test(y)$p.value < 0.05) {
          if (verbose) print("spearman correlation") # log
          p <- Hmisc::rcorr(x, y, type = "spearman")$P[1, 2]
          test <- "cor.spearman"
        } else {
          if (verbose) print("pearson correlation") # log
          p <- Hmisc::rcorr(x, y, type = "pearson")$P[1, 2]
          test <- "cor.pearson"
        }
      } # at least 4 values not na
    } # end block numeric variables
  }
  return(list(p = p, test = test))
}


# A novel method that tests pairwisely a given test
all.test.fixed <- function (x, y, adjust = "BH", verbose = FALSE) 
{
  p <- NA
  test <- "null"
  if (is.factor(x)) {
    if (is.factor(y)) { # if both variables are factors
      if (nrow(table(x, y)) == 1 | ncol(table(x, y)) == 1 | sum(table(x, y)) == 0 | sum(rowSums(table(x,y))!=0)<2 | sum(rowSums(table(y,x))!=0)<2) {
        warning("This is particular case 1")
      }else {
        if (verbose) print("chisq.test") # log
        #p <- chisq.test(table(x, y))$p.value
        
        p <- chisq.test(x, y)$p.value
        test <- "chisq.test"
      }
    } else { # if y is a numerical
      if(length(na.omit(y))>4){ # sanity check
        if (nortest::lillie.test(y)$p.value > 0.05) { # if non normal distribution
          if (length(table(x)) > 1) {
            if (length(table(x)) == 2) {
              if (all(colSums(table(y, x)) > 1)) {
                if (verbose) print("t.test") # log
                p <- t.test(y ~ x)$p.value
                test <- "t.test"
              }
            }else { # if more than two categories use a linear model
              if (sum(table(y, x)) > 1 & length(unique(x[!is.na(x)])) !=  1 & sum(colSums(table(y, x)) > 0) > 1) {
                if (verbose) print("linear model") # log
                tmp <- lmp(lm(y ~ x))
                if (!is.null(tmp)) {
                  p <- tmp
                  test <- "linear model"
                } # end is null result from linear model
              }
            }
          }
        } # end nortest
        else { # if non normal
          if (length(table(x)) > 1) {
            if (length(table(x)) == 2) {
              if (all(colSums(table(y, x)) > 1)) {
                if (verbose) print("wilcox.test") # log
                p <- wilcox.test(y ~ x)$p.value
                test <- "wilcox.test"
              }
            }else {
              if (sum(table(y, x)) > 1 & length(unique(x[!is.na(x)])) != 1 & sum(colSums(table(y, x)) > 0) > 1) {
                if (verbose) print("kruskal.test") # log
                tmp <- kruskal.test(y ~ x)
                if (!is.null(tmp)) {
                  p <- tmp$p.value
                  test <- "kruskal.test"
                } # end is null result from linear model
              } # end validity of the test
            } # end else
          } # end at least two categories
        } # end else normality
      } # sanity check
    } # end else numerical
  } # end x factor
  else { # If x is numerical
    if (is.factor(y)) { # and y is a factor
      if(length(na.omit(x))>4){ # sanity check
        if (nortest::lillie.test(x)$p.value > 0.05) { # if normal
          if (length(table(y)) > 1) { # if at least two factors
            if (length(table(y)) == 2) { # if exactly 2
              if (all(colSums(table(x, y)) > 1)) {
                if (verbose) print("t.test") # log
                p <- t.test(x ~ y)$p.value
                test <- "t.test"
              }
            }
            else {
              if (sum(table(x, y)) > 1 & length(unique(y[!is.na(y)])) != 1 & sum(colSums(table(x, y)) > 0) > 1) {
                if (verbose) print("linear model") # log
                tmp <- lmp(lm(x ~ y))
                if (!is.null(tmp)) {
                  p <- tmp
                  test <- "lm"
                }
              }
            }
          }
        } # end normality
        else { # if non normal
          if (length(table(y)) > 1) {
            if (length(table(y)) == 2) { # if only two categories
              if (all(colSums(table(x, y)) > 1)) {
                if (verbose) print("wilcox.test") # log
                p <- wilcox.test(x ~ y)$p.value
                test <- "wilcox.test"
              }
            } else {
              if (sum(table(x, y)) > 1 & length(unique(y[!is.na(y)])) != 1 & sum(colSums(table(x, y)) > 0) > 1) {
                if (verbose) print("kruskal.test") # log
                tmp <- kruskal.test(x ~ y)
                if (!is.null(tmp)) {
                  p <- tmp$p.value
                  test <- "kruskal.test"
                } # end is null
              } # end validity kurskal
            } 
          } # end of testing at least one category
        } # end else normality
      } # end sanity check
    } else { # if both variables are numeric
      if (length(x[!is.na(x)]) > 4 & length(x[!is.na(y)]) > 4) {
        if (nortest::lillie.test(x)$p.value < 0.05 | nortest::lillie.test(y)$p.value < 0.05) {
          if (verbose) print("spearman correlation") # log
          p <- Hmisc::rcorr(x, y, type = "spearman")$P[1, 2]
          test <- "cor.spearman"
        } else {
          if (verbose) print("pearson correlation") # log
          p <- Hmisc::rcorr(x, y, type = "pearson")$P[1, 2]
          test <- "cor.pearson"
        }
      } # at least 4 values not na
    } # end block numeric variables
  }
  return(list(p = p, test = test))
}

extractSignificant <- function (relation.matrix, interest, threshold = 0.05) {
  res <- list()
  for (i in 1:length(interest)) {
    tmp.val <- relation.matrix[, interest[i]]
    tmp <- (tmp.val < threshold)
    tmp[is.na(tmp)] <- FALSE
    res[[i]] <- tmp.val[tmp]
  }
  names(res) <- interest
  return(res)
}

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") 
    stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = F)
  attributes(p) <- NULL
  return(p)
}

#' resetFactors function.
#'
#' @title resetFactors
#' @description It will recreate the factor variable to relevel
#' @param data: a data frame containing the variables to explore
#' @param side: 2 for the columns, If something else the dataset will be transposed
#' @param character.2.factor: should the character vectors be transformed to factor
#' @param verbose: print information on the variables that are treated (default:FALSE)
#' @return the updated data frame
#' @export
resetFactors <- function(data, side = 2, character.2.factor = TRUE, verbose = FALSE)
{
  
  if(length(class(data)) > 1)
  {
    warning("resetFactors: class data type is not data.frame")
    data <- as.data.frame(data)
  }
  
  if(side!=2)
  {
    res <- data.frame(t(data))
  }else
  {
    res <- data.frame(data)  
  }
  
  types <- sapply(data.frame(res), class)
  
  for(i in 1:ncol(res))
  {
    if(verbose)
    {
      print(paste(i,"   =>   ",colnames(data)[i]))
    }
    if(types[i]=="factor")
    {
      res[,i] <- as.factor(as.character(data[,i]))  
    }
    if(character.2.factor)
    {
      if(types[i]=="character")
      {
        res[,i] <- as.factor(as.character(data[,i]))  
      }
    }
  }
  return(res)  
}


#' plotClinHistograms function.
#'
#' @title plotClinHistograms
#' @description function that plots histograms and barplots of a given dataset for descriptive analysis.
#' @param data: a data frame containing the variables to explore
#' @param side: 2 for the columns, If something else the dataset will be transposed
#' @param width: pdf width
#' @param height: pdf height
#' @param filename: pdf file name
#' @param mfrow: pdf mfrow
#' @param verbose: print information on the execution process
#' @param col: colors for the histograms
#' @return nothing
#' @export
plotClinHistograms <- function(data, side = 2, 
                               width=15, height=10, filename ="res.pdf", 
                               mfrow=c(4,7), verbose=TRUE, 
                               col=c("blueviolet","cornflowerblue","black")){
  # side in columns
  if(side!=2){
    print("Transposing data!")
    data <- as.data.frame(t(data))
  }
  
  if(length(class(data)) > 1)
  {
    warning("plotClinHistograms: class data type is not data.frame")
    data <- as.data.frame(data)
  }
  
  # make the plot
  pdf(file=filename, width=width, height=height)
  par(mfrow=mfrow)
  for(i in 1:ncol(data))
  {
    name <- colnames(data)[i]
    x <- data[,i]
    
    if(verbose) print(paste(i, name))
    
    if(all(is.na(x))) # if everything is empty plot empty
    {
      plot.new()
    }else
    {
      if(is.factor(x))
      {
        plot(x, main = name, xlab = "", col=col[1], las=2)  
      }else
      {
        if(is.numeric(x))
        {
          hist(x, main = name, xlab = "", col=col[2])  
        }else
        {
          warning("printing this as a date")
          plot(as.Date(x), main = name, xlab = "",col=col[3]) 
        }
      }
    } # end else empty
  } # end loop variable
  dev.off() # close pdf
}

#' main runRelome function.
#'
#' @title runRelome
#' @description Runs the different functions that are needed to compute and visualize relations beween variables.
#' @param data: a data frame containing the variables to explore
#' @param interest: a subset of variables to be extracted and used to focus the analysis
#' @param threshold: significance threshold of the variables to extract.
#' @param adjust: multiple testing adjustment method
#' @param adjust.by.var: weather to adjust only for one variable of interest or for all the tests performed.
#' @param zoom.p: use the p value as significance selector. If false it will be the adjusted p-value.
#' @param verbose: print running information.
#' @param plot: If TRUE produce the co-relation graphics.
#' @param mfrow: grid disposition of the graph (default:c(4,6))
#' @param width: PDF width in inches (default:15)
#' @param height: PDF heaight in inches (default:10)
#' @param save.all: weather to save intermediary files (default:TRUE)
#' @param return.all: weather to return all the results (default:TRUE)
#' @param bg: a list containg background color information for the graphs, scatterplots and tables.
#' @param col: a list containg border color information for the graphs, scatterplots and tables.
#' @param pch: the point shape for the scatterplot
#' @param inpdf: if TRUE the plots will be saved in a pdf, otherwise  they will be sent to the default canvas (default:TRUE).
#' @return an object containing the different results
#' @export
runRelome <- function(data, interest = "", threshold=0.05, 
                      adjust = "BH", adjust.by.var = TRUE, zoom.p = TRUE, 
                      verbose=TRUE,  plot=TRUE, 
                      mfrow = c(4,6), width = 15, height = 10, save.all=TRUE, rerun.all=FALSE,
                      bg = list(scatter="gold2",
                                 category=gray.colors(2)),
                      col = list(scatter="black",
                                 category=gray.colors(2)),
                      pch=21,
                      inpdf=TRUE
                      )
{
  
  if(file.exists("ppr.rda")){
    if(rerun.all){
      ppr <- relome(data = data, verbose = verbose, adjust = adjust)
      if(save.all){
        save(ppr, file="ppr.rda", compress=TRUE)
      }
    }else{
      load("ppr.rda")
    }
  }else{
    ppr <- relome(data = data, verbose = verbose, adjust = adjust, adjust.by.var = adjust.by.var)
    if(save.all){
      save(ppr, file="ppr.rda", compress=TRUE)
    }
  }
  
  # Richness linked phenotypes
  if(zoom.p)
  {
    if(verbose) print("Selecting variables based on p-values")
    rel <- lapply(extractSignificant(relation.matrix = ppr$p, interest = interest, threshold = threshold), sort, na.last=FALSE)
  }else
  {
    if(verbose) print("Selecting variables based on q-values")
    rel <- lapply(extractSignificant(relation.matrix = ppr$q, interest = interest, threshold = threshold), sort, na.last=TRUE)
  }
  if(plot){
    if(verbose) print("Making plots")
    for(i in 1:length(rel))
    {
      if (length(rel)==1)
      {
        file.name <- paste("relome_res_",i,"_", interest,".pdf",sep="")
      }else
      {
        file.name <- paste("relome_res_",i,"_", names(rel[i]),".pdf",sep="")
      }
      if(verbose) print(paste("Making plot", file.name))
      if(length(rel[[i]])==0)
      {
        print("No significant relations here ...")
      }else 
      {
        if(inpdf)
        {
          pdf(file=file.name, width=width, height=height)
          par(mfrow=mfrow)
          plotPhenoRelations(data = data, ppr = ppr, rel.list = rel[[i]], var.interest=interest[i], threshold = threshold, bg=bg, col=col, pch=pch)
          dev.off()
        }else
        {
          par(mfrow=mfrow)
          plotPhenoRelations(data = data, ppr = ppr, rel.list = rel[[i]], var.interest=interest[i], threshold = threshold, bg=bg, col=col, pch=pch)
        }
      }
    }
  }
  return(list(ppr=ppr, rel=rel))
}







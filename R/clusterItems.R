##
#' @title Item cluster approach to partial MI
#'
#' @description Find subsets of items (clusters) which function homogeneously
#' and may be possible anchor candidates.  Can currently deal with two-group models
#' for metric items or dichotomous items (Rasch, 2PL, or probit model).
#'
#' @param res_testMI Result object of \code{\link{testMI}}. Alternatively, a MI model
#' can be described using the arguments of \code{\link{testMI}} (see ...).
#' @param clusterWhat String, either \code{"loadings"}, \code{"difficulties"}, or
#' \code{c("loadings", "difficulties")}. Indicates, for which parameter type MI does not hold.
#' Usually, this will be indicated per factor in a single string
#' (e.g. \code{"Factor1 configural Factor2 weak"})
#' @param method Criteria used for clustering.
#' One can either specify the maximum difference in item parameters between any
#' two items within a cluster \code{"threshold"} or base clustering
#' on nonsignificance of item parameter differences of any two items within
#' a cluster \code{"sigTest"}. The threshold approach is based upon the approach described
#' in Pohl, Schulze, & Stets (in press); Pohl & Schulze (2020) and Schulze & Pohl (2021).
#' The significance test is used as described in Bechger & Maris (2015).
#' @param loadThreshold numeric value for the loading threshold. Only used if method is
#' \code{"threshold"} and if applicable in the model.
#' @param intThreshold numeric value for the intercept/difficulty threshold. Only used if
#' method is \code{"threshold"}.
#' @param alphaValue Type I error used for clustering via significance testing.
#' This argument only needs to be specified if method is \code{"sigtest"}.

#'
#' @param ... Arguments of \code{\link{testMI}}, if \code{\link{testMI}} has not
#' been called before and a MI model is to be described.
#'
#' @references
#' Bechger, T. M., & Maris, G. (2015). A statistical test for diﬀerential item pair functioning.
#' Psychometrika, 80(2), 317-340. \cr \cr
#' Pohl, S., Schulze, D., & Stets, E. (in press). Partial measurement invariance: Extending and
#' evaluating the cluster approach for identifying anchor items. Applied Psychological Assessment. \cr \cr
#' Pohl, S., & Schulze, D. (2020). Assessing group comparisons or change over
#' time under measurement non-invariance: The cluster approach for nonuniform DIF.
#' Psychological Test and Assessment Modeling, 62(2), 281-303. \cr \cr
#' Schulze, D., & Pohl, S. (2021). Finding clusters of measurement invariant items for
#' continuous covariates. Structural Equation Modeling: A Multidisciplinary Journal, 28(2), 219-228.
#'
#'
#' @return A list containing
#' \itemize{
#'  \item Data (filtered for potential missings in the covariate)
#'  \item Model specification(s)
#'  \item Clustering results organized per factor. If inherited from
#'  \code{\link{testMI}}, MI tests are also present.
#'  }
#' Usually accessed via summary().
#'
#' @usage res_clusterItems <- clusterItems(res_testMI = NULL,
#'                                         clusterWhat = NULL,
#'                                         method = NULL,
#'                                         alphaValue = NULL,
#'                                         loadThreshold = NULL,
#'                                         intThreshold = NULL,
#'                                         ...)
#'
#' @importFrom msm "deltamethod"
#'
#' @export


clusterItems <- function(res_testMI = NULL,
                         clusterWhat,      # loadings, difficulties
                         method,              # sigTest, threshold
                         alphaValue = NULL,   # p value for clustering by sigtest
                         loadThreshold = NULL,# threshold for loading DIF clustering
                         intThreshold = NULL, # threshold for intercept/difficulty DIF clustering # X moved to line
                         # inherited arguments from testMI
                         MIlevel,       # weak, strong, strict
                         items = NULL,
                         data = NULL,
                         group = NULL,
                         stdItems = TRUE,
                         itemType = "metric",
                         dichModel = "factor") {

  pValue <- alphaValue

  ## checks
  stopifnot(method %in% c("sigTest", "threshold"))

  if (method == "sigTest" && any(c(!is.null(loadThreshold),
                                   !is.null(intThreshold)))) {
    warning(paste0("Will ignore the following:",
                   c("loadThreshold",
                     "intThreshold")[c(!is.null(loadThreshold),
                                       !is.null(intThreshold))]))
  }
  if (method == "threshold" && any(c(!is.null(pValue)))) {
    warning(paste0("Will ignore the following: pValue"))
  }
  #if (anyDuplicated(unlist(res$items)) != 0) {
  #  stop("Cross loadings not allowed for clusterItems().", call. = FALSE)   # checking cross loadings
  #}

  ## initial model setup, if not inherited from testMI
  if (!is.null(res_testMI)) {
    res <- res_testMI
  } else {
    MItargetLevel <- MIlevel
    dich <- ifelse(itemType == "metric", FALSE, TRUE)
    res <- setUpModel(items = items,
                      data = data,
                      group = group,
                      MItargetLevel = MItargetLevel,
                      stdItems = stdItems,
                      dich = dich,
                      dichModel = dichModel)
  }
  items <- res$model$items
  lvs <- names(items)
  data <- res$data
  missings <- res$model$missings
  estim <- res$model$estim

  if (length(unique(res$data[, res$group])) > 2) {
    stop("More than 2 groups are not supported yet.", call. = FALSE)
  }

  ## parsing holding MI levels
  partialMI <- list()  # internal representation of input MI levels per factor
  #lvsPlaces <- lapply(lvs, FUN = regexpr, MIholdingLevel) # parse input string
  #lvsPlaces <- append(lvsPlaces, nchar(MIholdingLevel) + 1) # add final stop in order to make loop below work
  #levels <- c("configural",
  #            "weak",
  #            "strong",
  #            "strict")
  #if (any(lvsPlaces < 0) && !(MIholdingLevel %in% levels)) {
  #  stop("Error while reading MIholdingLevel. Couldn't retrieve factor names.", call. = FALSE)   # check for misspelled factor names
  #}
  for (lv in lvs) {
  #  if (MIholdingLevel %in% levels) {
  #    partialMI[[lv]] <- MIholdingLevel # catch cases, where a single MI level is stated for all factors
  #  } else {
  #    partialMI[[lv]] <- levels[which(lapply(levels,   # get desired MI level per factor
  #                                           FUN = regexpr,
  #                                           substring(MIholdingLevel,
  #                                                     lvsPlaces[[which(lvs %in% lv)]] +
  #                                                       attr(lvsPlaces[[which(lvs %in% lv)]], "match.length"),
  #                                                     lvsPlaces[[which(lvs %in% lv) + 1]] - 1)) > 0)]
  #    if (length(partialMI[[lv]]) == 0) {
  #      stop("Error while reading MIholdingLevel. Check your statements.", call. = FALSE)
  #    }
    partialMI[[lv]] <- clusterWhat
    partialMI[[lv]][partialMI[[lv]] == "loadings"] <- "weak"
    partialMI[[lv]][partialMI[[lv]] == "difficulties"] <- "strong"
  }
    res[[lv]]$clusterWhat <- clusterWhat
    #nonLevels <- c(1:which(levels %in% partialMI[[lv]]),  # delete MI levels <= holding level from list
    #               (which(levels %in% res[[lv]]$MItargetLevel) + 1):4)  # delete MI > desired level from list
    #partialMI[[lv]] <- levels[-nonLevels]  # delete
    if (res$model$dichModel == "Rasch") {
      partialMI[[lv]][partialMI[[lv]] == "weak"] <- NA
    }
  #}

  if (any(unlist(partialMI) == "strict", na.rm = TRUE)) {
    stop("Clustering of item residuals is not implemented yet. Please choose 'strong' instead of 'strict'.",
         call. = FALSE)
  }
  # switches
  catItems <- NULL
  int <- "intercepts"

  if (res$model$dich) {
    catItems <- items[[lv]]
    int <- "thresholds"
  }

  for (lv in lvs) { # loop over factors
    if (length(partialMI[[lv]]) == 0) {
      message(paste0("Skipping clustering for ", lv, " as MI level already holds."))
    }
    ##  multiple step approach (loadings -> intercepts)
    #  1st step: loadings (L)
    pkg <- ifelse(res$model$dichModel == "factor",
                  "lavaan", "mirt")
    drids <- NULL

    if ("weak" %in% partialMI[[lv]]) {
      # check
      loadings <- NULL
      for (lv0 in lvs) {
        loadings <- rbind(loadings, getItemParams(res[[lv0]]$configural[[1]],
                                                  type = "a1",
                                                  package = pkg))
      }
      if (min(loadings) < 0) {
        stop(paste0("Negative loadings occurred with ",
                    unique(rownames(which(loadings[, c(1, 2)] < 0, arr.ind = TRUE))),
                    ". Cannot continue."), call. = FALSE)
      }

      drids <- getDrids(res[[lv]]$configural[[1]], package = pkg) # getting relative DIFs for loadings

      if (method == "sigTest") {
        clusterLStep <- kMeansBMtest(res[[lv]]$configural[[1]],
                                     items[[lv]],
                                     type = "L",
                                     package = pkg,
                                     pValue = pValue)
      }
      if (method == "threshold") {
        clusterLStep <- kMeansThresh(drids$L[, 1], loadThreshold)
      }
    } else {
      clusterLStep <- NULL    # if partial MI is not needed for loadings this gives the following section something to chew on
      clusterLStep$cluster <- rep(1, length(items[[lv]]))
      names(clusterLStep$cluster) <- items[[lv]]
    }
    res[[lv]][["itemClustering"]]$finalClustering <- clusterLStep$cluster # set up list entry for final clustering
    res[[lv]][["itemClustering"]]$clusterLStep <- clusterLStep  # save Lstep result
    res[[lv]][["itemClustering"]]$clusterLStep$drids <- drids
    res[[lv]][["itemClustering"]]$notClustered <- FALSE
    res[[lv]][["itemClustering"]]$clusterMethod = c(method, pValue, loadThreshold, intThreshold)

    # 2nd step: thresholds / intercepts (A)
    if ("strong" %in% partialMI[[lv]]) {
      clusterLStep_1 <- unname(which(table(clusterLStep$cluster) == 1))       #  check for 1-item alpha clusters...
      clusterLStep_forAStep <- unname(which(table(clusterLStep$cluster) > 1)) #  ...and here all clusters with more than 1 item

      for (currLstep in clusterLStep_forAStep) {
        currCluster <- which(clusterLStep$cluster == currLstep)

        if (pkg == "lavaan") {
          mod1 <- paste0(res[[lv]]$configuralModel, "\n",
                         paste0("l", currCluster, "1 == l", currCluster, "2",
                                collapse = "\n"),  # add equality constraints for loadings of current anchor
                         "\n Factor~~ c(1, NA)*Factor")
          output <- cfa(mod1,
                        data = data,
                        estimator = estim,
                        group = res$model$group,
                        meanstructure = TRUE,
                        std.lv = TRUE,
                        ordered = catItems,
                        missing = missings)
        }
        if (pkg == "mirt") {
          dat0 <- res$data[, -which(colnames(res$data) %in% res$model$group)] # prepare data for mirt
          dat0 <- dat0[, res$model$items[[lv]]]
          anchorMod <- paste0("F1 = 1-", length(clusterLStep$cluster),"
                              CONSTRAINB = (", paste(currCluster, collapse = ", "),
                              ", a1)")
          constr <- "free_variances"
          if (res$model$dichModel == "Rasch") {
            anchorMod <- paste0("F1 = 1-", length(clusterLStep$cluster))
            constr <- ""
          }
          output <- multipleGroup(data = dat0,
                                  model = anchorMod,
                                  itemtype = res$model$dichModel,
                                  method = res$model$estim,
                                  group = as.factor(res$data[, res$model$group]),
                                  invariance = constr,
                                  SE = TRUE,
                                  verbose = FALSE)
        }
        drids <- getDrids(output,
                          package = pkg)

        if (method == "sigTest") {
          currAstep <- kMeansBMtest(output,
                                    names(currCluster),
                                    type = "A",
                                    package = pkg,
                                    pValue = pValue) # clustering
        }

        if (method == "threshold") {
          currAstep <- kMeansThresh(drids$A[currCluster, 1], intThreshold)
        }
        res[[lv]][["itemClustering"]][["modelAStep"]][[paste0("Loading cluster=", currLstep)]] <- list(drids = drids)
        res[[lv]][["itemClustering"]]$finalClustering[currCluster] <- currAstep$cluster + currLstep*1000   #  *1000 produces unique cluster labels
        historySteps <- list(clusteringHistory = NA,
                             cluster = rep(NA, length(clusterLStep$cluster)))  #  save history
        historySteps <- currAstep$clusteringHistory
        historySteps[["cluster"]][currCluster] <- currAstep$cluster
        res[[lv]][["itemClustering"]][[paste0("history: Load cluster=", currLstep)]] <- historySteps
      }
      for (i in clusterLStep_1) {
        res[[lv]][["itemClustering"]][["clusterAStep"]][clusterLStep$cluster == i] <- rep(i, 1)
      } #  for 1-item clusters: give cluster codes that will not occur otherwise (below 1000)
      res[[lv]][["itemClustering"]]$finalClustering <- as.integer(factor(res[[lv]][["itemClustering"]]$finalClustering))   # unify cluster labels
      res[[lv]][["itemClustering"]]$notClustered <- FALSE
    } else {
      if (length(partialMI[[lv]]) == 0) {
        res[[lv]][["itemClustering"]]$finalClustering <- rep(1, length(clusterLStep$cluster)) # if clustering of intercepts is not done: put dummy result here
        res[[lv]][["itemClustering"]]$notClustered <- TRUE
      }
    }
  }

  class(res) <- "clusterItems"   # give the result list a class attribute in order to be able to apply a customized summary-function
  return(res)
}

#' @title Summary function for clusterItems
#' @rdname summary-clusterItems
#' @description Giving a summary of the found item clusters per factor.
#' @param order Choose wether to order cluster solution by \code{"clusters"} (default)
#' or \code{"items"}.
#' @export

summary.clusterItems <- function(res,
                                 order = "clusters") { # items, clusters
  stopifnot(inherits(res, "clusterItems"))

  lvs <- names(res$model$items)

  if (res[[lvs[1]]]$itemClustering$clusterMethod[1] == "sigTest") {
    cat(paste0("Clustering by sign. test with p of ",
               res[[lvs[1]]]$itemClustering$clusterMethod[2], ".\n", "\n"))
  }
  if (res[[lvs[1]]]$itemClustering$clusterMethod[1] == "threshold") {
    cat(paste0("Clustering by threshold criterion with load threshold ",
               res[[lvs[1]]]$itemClustering$clusterMethod[3], " and intercept threshold ",
               res[[lvs[1]]]$itemClustering$clusterMethod[4], ".\n", "\n"))
  }


  for (lv in lvs) {
    clustering <- as.data.frame(res[[lv]]$itemClustering$finalClustering)
    cat(paste0(#"Factor: ", lv, "   ",
               if (res[[lv]]$itemClustering$notClustered) {"(not clustered)"
               } else {
                 paste0(max(clustering), " cluster",
                        if (max(clustering) > 1) "s", " found after clustering ",
                        paste(res[[lv]]$clusterWhat, collapse = " and ")#res[[lv]]$MIholdingLevel, " -> ",
                        #res[[lv]]$MItargetLevel, ")."
                        )},
               "\n"))
    rownames(clustering) <- res$model$items[[lv]]
    colnames(clustering) <- "cluster"
    if (order == "clusters") {
      clustering <- clustering[order(clustering$cluster), , drop = FALSE]
    }
    print(clustering)
  }
}


###################
# Level 2 functions


### calculate relative DIF values (aka DRIDs)
# is steered by parameter labels and thus rather general
# USED IN: clusterItems, kMeansBMtest
getDrids <- function(mod,
                     package = "lavaan"){ #lavaan, mirt, mplus
  if (package == "lavaan") {
    pars <- parameterestimates(mod)

    intcpts1 <- pars[grep("a.+1", pars$label), "est"]
    intcpts2 <- pars[grep("a.+2", pars$label), "est"]
    dridA <- outer(intcpts1, intcpts1, "-") - outer(intcpts2, intcpts2, "-")
    colnames(dridA) <- pars[grep("a.+1", pars$label), "lhs"]

    loads1 <- pars[grep("l.+1", pars$label), "est"]
    loads2 <- pars[grep("l.+2", pars$label), "est"]
    dridL <- outer(log(loads1), log(loads1), "-") - outer(log(loads2),
                                                          log(loads2), "-")
    colnames(dridL) <- pars[grep("l.+1", pars$label), "rhs"]
  } else { # mirt
    itemPars <- getItemParams(mod, SE = FALSE, package = "mirt")
    r1 <- outer(itemPars[, 1], itemPars[, 1], "-")
    r2 <- outer(itemPars[, 3], itemPars[, 3], "-")
    dridA <- r1 - r2
    colnames(dridA) <-  extract.mirt(mod, "itemnames")

    r1 <- outer(log(itemPars[, 2]), log(itemPars[, 2]), "-")
    r2 <- outer(log(itemPars[, 4]), log(itemPars[, 4]), "-")
    dridL <- r1 - r2
    colnames(dridL) <- extract.mirt(mod, "itemnames")
  }

  res <- list(L = dridL, A = dridA)
  return(res)
}


### kmeans clustering controlled by Bechger & Maris' test  # Sat May 22 17:00:17 2021 ------------------------------
# Can cluster for alpha or lambda. 2 group only.
# USED IN: clusterItems
kMeansBMtest <- function(output, # fitted lavaan/mirt model
                         items,  # item names
                         type,   # "A" or "L" for clustering alpha or lambda
                         package,
                         pValue) {  # cut-off p value to be used
  k <- 1
  resL <- list()
  drids <- getDrids(output, package)

  history <- list(list(clusCode = rep(1, length(items)),
                       test = DIFtestBM(type,
                                        which(colnames(drids$L) %in% items),  # filter items
                                        drids[[type]][, which(colnames(drids$L) %in% items)[1]][which(colnames(drids$L) %in% items)],
                                        output,
                                        package)))
  currPs <- history[[1]][["test"]][3]

  while (min(currPs, na.rm = TRUE) < pValue) {
    currPs <- NULL
    k <- k + 1
    if (k == length(items)) {
      break
    }
    clusCode <- Ckmeans.1d.dp::Ckmeans.1d.dp(drids[[type]][, which(colnames(drids$L) %in% items)[1]][which(colnames(drids$L) %in% items)], k = k)$cluster
    multipleItemClusters <- unname(which(table(clusCode) != 1))        # filter for clusters with more than 1 item
    history[[k]] <- list(clusCode = clusCode)
    for (i in multipleItemClusters) {
      history[[k]]$test[[i]] <- DIFtestBM(type,
                                          testItems = which(clusCode %in% i),    # take items of current cluster
                                          drids = drids[[type]][, which(clusCode %in% i)[1]][which(clusCode %in% i)], # take column of first item of the cluster
                                          output,
                                          package)
      currPs[i] <- history[[k]]$test[[i]][3]
    }
  }
  if (k == length(items)) {
    resL[["cluster"]] <- seq(1:length(items))
  } else {
    resL[["cluster"]] <- Ckmeans.1d.dp:::Ckmeans.1d.dp(drids[[type]][, which(colnames(drids$L) %in% items)[1]][which(colnames(drids$L) %in% items)], k = k)$cluster
  }
  names(resL[["cluster"]]) <- items
  resL[["clusteringHistory"]] <- history
  return(resL)
}


### Threshold clustering for any parameter type
# USED IN: clusterItems
kMeansThresh <- function(drids,         # numeric[n]: DRID-values for the items (taken from delta-R matrix)
                         thresh = NULL) {  # numeric[n]: Threshold for the stopping rule. Threshold is maximum cluster width on logit scale # X added whitespace and moved to line
  k <- 1
  currRange <- max(drids) - min(drids)
  while (max(currRange) > thresh && k < length(drids)) {
    k <- k + 1
    clusCode <- Ckmeans.1d.dp(drids, k = k)$cluster                             # X added whitespace arund =
    currRange <- NULL
    for (i in 1:max(clusCode)) {
      clusCode2 <- clusCode     # helper copy
      clusCode2[clusCode2 != i] <- NA         # make picking vector
      clusCode2[clusCode2 == i] <- 1
      currRange[i] <- max(clusCode2*drids, na.rm = TRUE) - min(clusCode2*drids,
                                                               na.rm = TRUE)   # get current range # X added whitespace around X, T -> TRUE
    }
  }
  res <- list(cluster = Ckmeans.1d.dp(drids, k = k)$cluster)                    # X added whitespace around =
  return(res)
}



###################
# Level 3 functions

### Do a single Bechger & Maris-style signifiance test for significant DIF # Sat May 22 15:39:44 2021 ------------------------------
# USED IN: kMeansBMtest()
DIFtestBM <- function(type,      # either "A" or "L" for alpha or lambda parameter
                      testItems, # item numbers of of items to test
                      drids,     # relative DIFs
                      output,    # lavaan/mirt output
                      package) { # lavaan /mirt
  if (type == "A") {
    if (package == "lavaan") {
      paramNames <- list(paste0("a", outer(testItems, 1, FUN = paste0)),
                         paste0("a", outer(testItems, 2, FUN = paste0)))
    }
    if (package == "mirt") {
      mirtIndices <- rownames(vcov(output))
      indNumb <- grep("d.", mirtIndices)
      paramNames <- list(mirtIndices[indNumb[testItems]],
                         mirtIndices[indNumb[length(indNumb)/2 + testItems]])
    }
    sigma <- calCov(vcov(output), paramNames)
  }

  if (type == "L") {
    if (package == "lavaan") {
      paramNames <- list(paste0("l", outer(testItems, 1, FUN = paste0)),
                         paste0("l", outer(testItems, 2, FUN = paste0)))
      pars <- parameterestimates(output)
      pars <- pars[pars$label %in% unlist(paramNames), "est"]
    }
    if (package == "mirt") {
      mirtIndices <- rownames(vcov(output))
      indNumb <- grep("a1.", mirtIndices)
      paramNames <- list(mirtIndices[indNumb[testItems]],
                         mirtIndices[indNumb[length(indNumb)/2 + testItems]])
      pars <- getItemParams(output, "mirt", "a1", SE = FALSE)
      pars <- c(pars[testItems, 1], pars[testItems, 2])
    }
    sigma <- calCov(vcov(output), paramNames, pars, type = "L")
  }
  itemNumber <- length(testItems)
  DIFtest <- mahalanobis(drids[-1], rep(0,(itemNumber - 1)), sigma[-1, -1])  # as taken from Bechger & Maris (2015)
  DIFp <- 1 - pchisq(DIFtest, itemNumber - 1)
  res <- c(DIFtest, itemNumber - 1, DIFp)
  return(res)
}


#### retrieve item params from models in a digestable way.
## USED IN: getDrids()
getItemParams <- function(mod,
                          package,               # lavaan, mirt, mplus
                          type = c("d", "a1"),   # character[i]: names of the item parameters to be extracted (mirt terminology: d/a1)
                          SE = TRUE              # logical: should parameter standard errors be reported as well?
){
  if (package == "mirt") {
    containsSE <- !is.na(vcov(mod)[1,1]) # does the model include SEs? (estimated with SE = TRUE?)
    if (!containsSE && SE) {
      SE <- FALSE
    }
    itNames <- extract.mirt(mod, "itemnames")
    grpNames <- extract.mirt(mod, "groupNames")
    out <- data.frame(matrix(NA, nrow = length(itNames), ncol = 1))
    cfg <- coef(mod, printSE = TRUE)
    res <- data.frame(t(data.frame(cfg)))
    for (g in grpNames) {
      if (length(grpNames) > 1) {
        grpIdx <- grep(g, rownames(res))
      } else {
        grpIdx <- 1:nrow(res)
      }
      s <- strsplit(rownames(res),".", fixed = TRUE)
      parVec <- sapply(s, function(x) x[length(x)])
      for (p in type) {
        parIdx <- which(parVec == p)
        curr <- res[intersect(parIdx, grpIdx),]
        if (!is.numeric(curr)) curr <- curr[,"par"]
        if (length(grpNames) == 1) {
          colName <- p
        } else {
          colName <- paste(g, p, sep = ".")
        }
        out[colName] <- curr
        if (SE && containsSE) {
          curr <- res[intersect(parIdx,grpIdx),]
          out[paste(colName, "SE",sep = ".")] <- curr[, "SE"]
        }
      }
    }
    out <- out[-1] # this is a really dirty hack but necessary because the preallocation of out above creates a stupidly named first variable
    rownames(out) <- itNames
  }
  if (package == "lavaan") {
    type <- "=~"
    pars <- parTable(mod)
    itemNames <- unique(pars[pars$op == "=~", "rhs"])
    out <- matrix(nrow = length(itemNames), ncol = 4)
    rownames(out) <- itemNames
    colnames(out) <- c("par1", "par2", "se1", "se2")

    out[, 1] <- pars[pars$op == type & pars$group == 1, "est"]
    out[, 2] <- pars[pars$op == type & pars$group == 2, "est"]
    out[, 3] <- pars[pars$op == type & pars$group == 1, "se"]
    out[, 4] <- pars[pars$op == type & pars$group == 2, "se"]
  }
  return(as.data.frame(out))
}



###################
# Level 4 functions

### Calculates transformation of parameter vcov matrix as needed for BM sig test # Sat May 22 15:36:23 2021 ------------------------------
calCov <- function(vcov,
                   paramNames,
                   params,
                   type = "A") {     # for L or A (loadings or intercepts)
  if (type == "L") {
    vcov <- vcov[unlist(paramNames), unlist(paramNames)]
    formulas <- list()
    for (f in 1:length(unlist(params))) {
      formulas <- append(formulas, as.formula(paste0("~log(x", f, ")")))  # looping here is hacky but couldnt think of another easy way
    }
    vcov <- as.matrix(msm::deltamethod(formulas,
                                       unlist(params),
                                       vcov,
                                       ses = FALSE),
                      length(unlist(params)),
                      length(unlist(params))) # convert cov matrix to log scale via the delta method

    colnames(vcov) <- rownames(vcov) <- unlist(paramNames)
  }

  sigma <- matrix(0,
                  nrow = length(paramNames[[1]]),
                  ncol = length(paramNames[[1]]))
  for (g in 1:2) {
    for (i in paramNames[[g]]) {
      for (j in paramNames[[g]]) {
        sigma[which(paramNames[[g]] %in% i), which(paramNames[[g]] %in% j)]  <-
          sigma[which(paramNames[[g]] %in% i), which(paramNames[[g]] %in% j)] +
          vcov[i, j] +    # formula from Bechger & Maris, footnote 6
          vcov[paramNames[[g]][1], paramNames[[g]][1]] -
          vcov[i, paramNames[[g]][1]] -
          vcov[j, paramNames[[g]][1]]
      }
    }
  }
  return(sigma)
}

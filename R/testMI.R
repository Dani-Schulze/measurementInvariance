# # Fri Jun 11 09:38:18 2021 ------------------------------

#' @title Global MI test
#'
#' @description Global measurement invariance (MI) tests for the common MI levels
#' (configural, weak, strong, and strict).  Can currently deal with multiple-group models (two or more groups)
#' for metric items or dichotomous items (Rasch, 2PL, or probit model).
#'
# @param model A SEM written in \code{\link[lavaan:lavaan]{lavaan}}'s syntax. Multidimensional models will be
# broken apart internally into single-factor models for the purpose of MI testing.
#' @param items Describes the measurement model. It can either be:
#' \itemize{
#'  \item A vector of item names, e.g. \code{c("item1", "item2", "item3", "item4")}
#'  \item A data frame section containing the items e.g. \code{data[1:4]}
#'  }
#' This creates a single factor model out of the item names.
# The latent variable
# is internally called "Factor". \cr  \cr
# For testing named factors or multiple factors sequentially it can either be:
# \itemize{
#  \item A named list with vectors of item names, e.g.
#        \code{list(myFactor = c("item1", "item2", "item3", "item4")}
#  \item A named list with data frame sections containing the items , e.g.
#       \code{list(myFactor = data[1:4]}
# \item or a combination thereof.
#  }
#' @param data A data frame containing all the used variables. Missing values have
#' to be declared as \code{NA}.
#' @param group String, indicating the grouping variable.
#' @param MIlevel String, indicating the level of MI to be tested. Can be either \code{configural}, \code{weak},
#' \code{strong}, or \code{strict}.
# Alternatively,
# can be indicated per factor of interest in a single string
# (e.g. \code{"Factor1 weak Factor2 strong"}).
# If MI is not to be tested for all factors mentioned in \code{model}
#  the user can simply drop
# those factors from this list (e.g. \code{"Factor1 weak"} ).
#' @param stdItems Should continuous items be standardized? Defaults to TRUE.
#' Useful for making DIF sizes comparable.
#' @param itemType Can be \code{"dichotomous"} or \code{"metric"} (default).
#' @param dichModel String indicating which measurement model is applied in the case of dichotomous items.
#' Values are \code{"factor"} (default) for probit link via WLSMV in lavaan.
#' \code{"Rasch"} or  \code{"2PL"} for the respective IRT models to be estimated
#' with the EM estimator in mirt.
#'
#' @return A list containing
#' \itemize{
#'  \item Data (filtered for potential missings in the covariate)
#'  \item Model specification(s)
#'  \item Results of the MI tests
# organized per factor
#'  }
#' Usually accessed via summary() or getModel().
#'
#' @usage res_testMI <- testMI(items = NULL,
#'                      data = NULL,
#'                      group = NULL,
#'                      MIlevel = NULL,
#'                      stdItems = TRUE,
#'                      itemType = "metric",
#'                      dichModel = "factor")
#Test several factors sequentially:
#res_testMI <- testMI(items = list(Factor1 = c("item1", "item2", "item3", "item4"),
#                                  Factor2 = c("item5", "item6", "item7", "item8")),
#                     data = data,
#                     group = "country",
#                     MItargetLevel = "strong")
#'
#' @export

testMI <- function(items = NULL,
                   data,
                   group,
                   MIlevel,
                   stdItems = TRUE,
                   itemType = "metric",
                   dichModel = "factor") {

  MItargetLevel <- MIlevel
  if (!exists('MItargetLevel')) {
    stop("Measurement invariance (MIlevel) type must be specified.")
  }

  dich <- ifelse(itemType == "metric", FALSE, TRUE)

  res <- setUpModel(items = items,
                    data = data,
                    group = group,
                    MItargetLevel = MItargetLevel,
                    stdItems = stdItems,
                    dich = dich,
                    dichModel = dichModel)

  items <- res$model$items
  lvs <- names(items)
  data <- res$data
  missings <- res$model$missings
  estim <- res$model$estim

  ### fit models separately per factor
  for (lv in lvs) {

    ## lavaan models (metric items, dich items with non-IRT model)
    if (dichModel == "factor") {
      # switches
      catItems <- NULL
      int <- "intercepts"
      if (dich) {
        catItems <- items[[lv]]
        int <- "thresholds"
      }
      mod <- paste0(lv, "=~", paste(items[[lv]], collapse = " + "))

      tryCatch({
        if (res[[lv]]$MItargetLevel %in% c("weak", "strong", "strict")) {
          if (!dich) {
            cfaweak <- cfa(mod,
                           data = data,
                           estimator = estim,
                           group = group,
                           meanstructure = TRUE,
                           std.lv = TRUE,
                           ordered = catItems,
                           group.equal = c("loadings"),
                           missing = res$model$missings)
            res[[lv]][["weak"]] <- list(cfaweak)
          }
          if (res[[lv]]$MItargetLevel %in% c("strong", "strict")) {
            cfastrong <- cfa(mod,
                             data = data,
                             estimator = estim,
                             group = group,
                             meanstructure = TRUE,
                             std.lv = TRUE,
                             ordered = catItems,
                             group.equal = c("loadings", int),
                             missing = res$model$missings)
            res[[lv]][["strong"]] <- list(cfastrong)

            if (res[[lv]]$MItargetLevel == "strict") {
              cfastrict <- cfa(mod,
                               data = data,
                               estimator = estim,
                               group = group,
                               meanstructure = TRUE,
                               std.lv = TRUE,
                               ordered = catItems,
                               group.equal = c("loadings", int, "residuals"),
                               missing = res$model$missings)
              res[[lv]][["strict"]] <- list(cfastrict)
            }
          }
        }
      })
    }

    ## IRT models in mirt
    if (dichModel != "factor") {
      # checks
      if (res[[lv]]$MItargetLevel == "strict") {
        stop("Strict MI is not relevant in IRT models.", call. = FALSE)
      }

      if (length(lvs) > 2) estim <- "MHRM"
      dat0 <- data[, unlist(items[[lv]])]  # prepare data for mirt

      if (res[[lv]]$MItargetLevel %in% c("weak", "strong")) {
        if (dichModel == "Rasch") {
          message("Skipping weak MI for Rasch model.")
        } else {
          cfaweak <- multipleGroup(dat0,
                                   model = 1,
                                   itemtype = dichModel,
                                   group = as.factor(data[, group]),
                                   invariance = c("slopes", "free_variances"),
                                   method = estim,
                                   SE = TRUE,
                                   verbose = FALSE)
          res[[lv]][["weak"]] <- list(cfaweak)
        }
        if (res[[lv]]$MItargetLevel %in% c("strong")) {
          cfastrong <- multipleGroup(dat0,
                                     model = 1,
                                     itemtype = dichModel,
                                     group = as.factor(data[, group]),
                                     invariance = c("intercepts", "free_means", "slopes", "free_variances"),
                                     method = estim,
                                     SE = TRUE,
                                     verbose = FALSE)
          res[[lv]][["strong"]] <- list(cfastrong)
        }
      }
    }
  }

  class(res) <- "testMI"
  return(res)
}



#' @title Summary function for testMI
#' @rdname summary-testMI
#' @description Giving a summary of fit indices of an estimated testMI object per factor.
#' @export

summary.testMI <- function(res) {  # Tue May 25 08:57:08 2021 ------------------------------
  stopifnot(inherits(res, "testMI"))

  lvs <- names(res$model$items)
  row <- 0

  cat(paste0("Two group model", if (length(res$lvs) > 1) "s", " with:\n"))
  tab <- t(data.frame(unclass(table(res$data[[res$model$group]]))))
  tab <- cbind(tab, sum(tab[1, ]))
  rownames(tab) <- ""
  colnames(tab)[ncol(tab)] <- "total"
  print(tab)
  cat("\n")

  if (res$model$dichModel == "factor") {
## lavaan
    Res <- data.frame(matrix(0, ncol = 15))

    inds <- c("chisq.scaled",
              "df.scaled",
              "pvalue.scaled",
              "cfi.scaled",
              "rmsea.scaled",
              "rmsea.ci.lower.scaled" ,
              "srmr")

    Tech <- data.frame(c("      lavaan",  # tech summary
                         res$model$estim,
                         if (!res$model$dich) "continuous" else "dichotomous",
                         res$model$missingsOut,
                         if (!res$model$dich) if (res$model$stdItems) "yes" else "no"))
    rownames(Tech) <- c("package",
                        "estimator",
                        "item type",
                        "item missings",
                        if (!res$model$dich) "standardized")
    colnames(Tech) <- ""

    for (lv in lvs) {
      row <- row + 1
      drids <- getDrids(res[[lv]]$configural[[1]], "lavaan")
      cnfg <- fitMeasures(res[[lv]]$configural[[1]])
      Res[row, 2:8] <- cnfg[inds]

      if (res[[lv]]$MItargetLevel %in% c("weak", "strong", "strict") && !res$model$dich) {
        weak <- fitMeasures(res[[lv]]$weak[[1]])
        row <- row + 1
        Res[row, 2:15] <- c(weak[inds],
                            anova(res[[lv]]$configural[[1]],
                                  res[[lv]]$weak[[1]])[2, 5:7],
                            weak["cfi.scaled"] - cnfg["cfi.scaled"],
                            weak["rmsea.scaled"] - cnfg["rmsea.scaled"],
                            weak["srmr"] - cnfg["srmr"],
                            abs(max(drids$L[, 1] -  min(drids$L[, 1]))))
      }
      if (res[[lv]]$MItargetLevel %in% c("weak", "strong", "strict") && res$model$dich) {
        row <- row + 1
        Res[row, 2:15] <- 0
      }
      if (res[[lv]]$MItargetLevel %in% c("strong", "strict") && !res$model$dich) {
        strong <- fitMeasures(res[[lv]]$strong[[1]])
        row <- row + 1
        Res[row, 2:15] <- c(strong[inds],
                            anova(res[[lv]]$weak[[1]],
                                  res[[lv]]$strong[[1]])[2, 5:7],
                            strong["cfi.scaled"] - weak["cfi.scaled"],
                            strong["rmsea.scaled"] - weak["rmsea.scaled"],
                            strong["srmr"] - weak["srmr"],
                            abs(max(drids$A[, 1] -  min(drids$A[, 1]))))
      }
      if (res[[lv]]$MItargetLevel %in% c("strong", "strict") && res$model$dich) {
        strong <- fitMeasures(res[[lv]]$strong[[1]])
        row <- row + 1
        Res[row, 2:15] <- c(strong[inds],
                            anova(res[[lv]]$configural[[1]],
                                  res[[lv]]$strong[[1]])[2, 5:7],
                            strong["cfi.scaled"] - cnfg["cfi.scaled"],
                            strong["rmsea.scaled"] - cnfg["rmsea.scaled"],
                            strong["srmr"] - cnfg["srmr"],
                            abs(max(drids$A[, 1] -  min(drids$A[, 1]))))
      }
      if (res[[lv]]$MItargetLevel %in% "strict") {
        strict <- fitMeasures(res[[lv]]$strict[[1]])
        row <- row + 1
        Res[row, 2:15] <- c(strict[inds],
                            anova(res[[lv]]$strong[[1]],
                                  res[[lv]]$strict[[1]])[2, 5:7],
                            strict["cfi.scaled"] - strong["cfi.scaled"],
                            strict["rmsea.scaled"] - strong["rmsea.scaled"],
                            strict["srmr"] - strong["srmr"], "-")
      }
    }
    colnames(Res) <- c("MI level",
                       "chi2",
                       "df",
                       "p",
                       "CFI",
                       "RMSEA",
                       "RMSEA 90% lower",
                       "SRMR",
                       "diff chi2",
                       "diff df",
                       "diff p",
                       "diff CFI",
                       "diff RMSEA",
                       "diff SRMR",
                       "Effect (DIF range)")
    Res$chi2 <- round(Res$chi2, 2)
    Res[, 4:15] <- round(Res[, 4:15], 3)

    MIlevels <- c("configural", "weak", "strong", "strict")
    levels <- NULL
    for (lv in lvs) {
      levels <- append(levels,
                       MIlevels[1:which(c("configural", "weak", "strong", "strict") %in% res[[lv]]$MItargetLevel)])
    }
    Res[, 1] <- levels

    Res[Res$'MI level' == "configural", 9:15] <- rep("", 7)
    if (res$model$dich) {
      Res[Res$'MI level' == "weak", 2:15] <- rep("-", 14)
    }
    Res <- as.data.frame(t(Res))

    # nam <- NULL
    # for (lv in lvs) nam <- append(nam,
    #                              c(lv, rep("",
    #                                        which(c("configural", "weak", "strong", "strict") %in% res[[lv]]$MItargetLevel) - 1)))
  colnames(Res) <- NULL


  } else {
## mirt
    Res <- data.frame(matrix(0, ncol = 16))
    Tech <- data.frame(c("     mirt",  # tech summary
                         res$model$estim,
                         res$model$dichModel,
                         res$model$missingsOut))
    rownames(Tech) <- c("package",
                        "estimator",
                        "item type",
                        "item missings")
    colnames(Tech) <- ""

    for (lv in lvs) {
      row <- row + 1
      drids <- getDrids(res[[lv]]$configural[[1]], "mirt")
      if (res$model$missings == "none") {
        M2conf <- M2(res[[lv]]$configural[[1]])
      } else {
        M2conf <- rep(0, 10)
      }

      Res[row, 2:10] <- c(anova(res[[lv]]$configural[[1]])[c(1, 3, 5)],
                                      M2conf[c(1:5, 10)])

      if (res[[lv]]$MItargetLevel %in% c("weak", "strong") &&
          res$model$dichModel != "Rasch") {

        if (res$model$missings == "none") {
          M2weak <- M2(res[[lv]]$weak[[1]])
        } else {
          M2weak <- rep(0, 10)
        }
        row <- row + 1
        Res[row, 2:16] <- c(anova(res[[lv]]$weak[[1]],
                                              res[[lv]]$configural[[1]],
                                              verbose = FALSE)[1, c(1, 3, 5)],
                                        M2weak[c(1:5, 10)],
                                        anova(res[[lv]]$weak[[1]],
                                              res[[lv]]$configural[[1]],
                                              verbose = FALSE)[2, 7:9],
                                        -(M2conf - M2weak)[c(4, 10)],
                                      abs(max(drids$L[, 1] -  min(drids$L[, 1]))))
      }

      if (res[[lv]]$MItargetLevel %in% c("weak", "strong") &&
          res$model$dichModel == "Rasch") {
        row <- row + 1
        Res[row, 2:16] <- 0
      }
      if (res[[lv]]$MItargetLevel %in% "strong" &&
          res$model$dichModel != "Rasch") {

        if (res$model$missings == "none") {
          M2strong <- M2(res[[lv]]$strong[[1]])
        } else {
          M2strong <- rep(0, 10)
        }
        row <- row + 1
        Res[row, 2:16] <- c(anova(res[[lv]]$strong[[1]],
                                              res[[lv]]$weak[[1]],
                                              verbose = FALSE)[1, c(1, 3, 5)],
                                        M2strong[c(1:5, 10)],
                                        anova(res[[lv]]$strong[[1]],
                                              res[[lv]]$weak[[1]],
                                              verbose = FALSE)[2, 7:9],
                                        -(M2weak - M2strong)[c(4, 10)],
                            abs(max(drids$A[, 1] -  min(drids$A[, 1]))))
      }
      if (res[[lv]]$MItargetLevel %in% "strong" &&
          res$model$dichModel == "Rasch") {

        if (res$model$missings == "none") {
          M2strong <- M2(res[[lv]]$strong[[1]])
        } else {
          M2strong <- rep(0, 10)
        }
        row <- row + 1
        Res[row, 2:16] <- c(anova(res[[lv]]$strong[[1]],
                                              res[[lv]]$configural[[1]],
                                              verbose = FALSE)[1, c(1, 3, 5)],
                                        M2strong[c(1:5, 10)],
                                        anova(res[[lv]]$strong[[1]],
                                              res[[lv]]$configural[[1]],
                                              verbose = FALSE)[2, 7:9],
                                        -(M2conf - M2strong)[c(4, 10)],
                            abs(max(drids$A[, 1] -  min(drids$A[, 1]))))
      }
    }
    colnames(Res) <- c("MI level",
                       "AIC",
                       "SABIC",
                       "BIC",
                       "M2",
                       "df",
                       "p",
                       "RMSEA",
                       "RMSEA 90% lower",
                       "CFI",
                       "diff chi2",
                       "diff df",
                       "diff p",
                       "diff RMSEA",
                       "diff CFI",
                       "Effect (DIF range)")
    Res[, 2:16] <- round(Res[, 2:16], 3)

    MIlevels <- c("configural", "weak", "strong")
    levels <- NULL
    for (lv in lvs) {
      levels <- append(levels,
                       MIlevels[1:which(c("configural", "weak", "strong") %in% res[[lv]]$MItargetLevel)])
    }
    Res[, 1] <- levels

    if (!(res$model$missings == "none")) {
      Res[, c(4:10, 14:15)] <- "*"
    }
    Res[Res$'MI level' == "configural", 11:16] <- rep("", 6)
    if (res$model$dichModel == "Rasch") {
      Res[Res$'MI level' == "weak", 2:16] <- rep("-", 15)
    }
    Res <- as.data.frame(t(Res))

    #nam <- NULL
    #for (lv in lvs) nam <- append(nam,
    #                              c(lv, rep("",
    #                                        which(c("configural", "weak", "strong") %in% res[[lv]]$MItargetLevel) - 1)))
    colnames(Res) <- NULL
  }

  print(Res)
  print(Tech)

  if (!(res$model$missings == "none") &&
      res$model$dichModel != "factor") {
    cat("\n* Some fit statistics could not be calculated due to missing values. (Maydeu-Olivares & Joe, 2006)\n")
  }

  if (any(res$model$MI %in% c("weak", "strong", "strict")) &&
      res$model$dichModel == "factor" && res$model$dich) {
    cat("\nWeak MI not testable with dichotomous items. (Wu & Estabrook, 2016)")
  }
  cat("\nUse getModel() to access parameter estimates or to further process the results.")
}

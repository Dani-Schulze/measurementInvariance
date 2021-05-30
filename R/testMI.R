# Sun May 23 14:19:57 2021 ------------------------------

#' @title testMI
#'
#' @description Global measurement invariance (MI) tests for the common MI levels
#' (configural, weak, strong, and strict).
#'
#' @param model A SEM written in \code{\link[lavaan:lavaan]{lavaan}}'s syntax. Multidimensional models will be
#' broken apart internally into single-factor models for the purpose of MI testing.
#' @param items Alternative to the model argument or lazy users with many items,
#' creates a single factor model out of the item names. The latent variable is internally
#' called "Factor".
#' @param data A data frame containing all variables of relevance. Missings have
#' to be declared \code{NA}.
#' @param group String, indicating the grouping variable
#' @param MIlevel String, either \code{configural}
#' , \code{weak}
#' , \code{strong}
#' , or \code{strict}
#' . Alternatively,
#' can be indicated per factor of interest in a single string
#' (e.g. \code{"Factor1 weak Factor2 strong"}
#' ).
#' If MI shall not be tested for all factors mentioned in \code{model}
#'  the user can simply drop
#' those factors from this list (e.g. \code{"Factor1 weak"}
#' ).
#' @param stdItems Should continuous items be standardized? Defaults to TRUE.
#' @param dich Are the items dichotomous? Defaults to FALSE.
#' @param dichModel For dichotomous items, which model type should be applied?
#' Can be \code{"factor"}
#' (default) for a categorical SEM in lavaan, \code{"Rasch"}
#'  or  \code{"2PL"}
#' for the respective IRT models to be estimated in mirt.
#'
#' @return A list with all estimated models (and some internally used info).
#' Usually accessed via summary() or getModel().
#'
#' @usage model = "Factor1 =~ item1 + item2 + item3 + item4
#'                 Factor2 =~ item5 + item6 + item7 + item8"
#' res_testMI <- testMI(model,
#'                      data,
#'                      group = "country",
#'                      MIlevel = "strong")
#'
#' @export

testMI <- function(model = NULL,
                   items = NULL,
                   data,
                   group,
                   MIlevel,
                   stdItems = TRUE,
                   dich = FALSE,
                   dichModel = "factor") {

  if (!exists('MIlevel')) {
    stop("Measurement invariance (MIlevel) type must be specified.")
  }
  stopifnot(dichModel %in% c("factor", "2PL", "Rasch"))

  ## parse model structure
  if (!is.null(model)) {
    pars0 <- lavParTable(model) # get basic model structure
    lvs <- unique(pars0[pars0$op == "=~", "lhs"]) # get factor names in model
  } else {
    lvs <- "Factor"
  }

  levels <- c("configural",
              "weak",
              "strong",
              "strict")
  if (!(MIlevel %in% levels)) {
    lvs <- lvs[sapply(lvs, FUN = grepl, MIlevel)] # filter for factors that are actually mentioned in MIlevel statement (leaving out factors not of interest)
  }

  if (!is.null(model)) {
    items <- list()
    for (i in lvs) { # get item assignments
      items[[i]] <- unique(pars0[pars0$op == "=~" & pars0$lhs == i, "rhs"])
    }
  } else {
    items0 <- list(Factor = items)
    items <- items0
  }

  ## parsing desired MI levels
  MI <- list()  # internal representation of input MI levels per factor
  lvsPlaces <- lapply(lvs, FUN = regexpr, MIlevel) # parse input string
  lvsPlaces <- append(lvsPlaces, nchar(MIlevel) + 1) # add final stop in order to make loop below work
  if (any(lvsPlaces < 0) && !(MIlevel %in% levels) &&
      nchar(MIlevel) > (sum(nchar(lvs)) + 5*length(lvs))) {
    stop("Error while reading MIlevel. Couldn't retrieve factor names.", call. = FALSE)   # check for misspelled factor names
  }
  for (lv in lvs) {
    if (MIlevel %in% levels) {
      MI[[lv]] <- MIlevel # catch cases, where a single MI level is stated for all factors
    } else {
      MI[[lv]] <- levels[which(lapply(levels,   # get desired MI level per factor
                                      FUN = regexpr,
                                      substring(MIlevel,
                                                lvsPlaces[[which(lvs %in% lv)]] +
                                                  attr(lvsPlaces[[which(lvs %in% lv)]], "match.length"),
                                                lvsPlaces[[which(lvs %in% lv) + 1]] - 1)) > 0)]
    }
  }
  if (length(MI) == 0) {
    stop("Error while reading MIlevel. Check your statements.", call. = FALSE)
    if (any(apply(MI, FUN = length) == 0)) {
      stop("Error while reading MIlevel. Check your statements.", call. = FALSE)
    }
  }

  nGroups <- length(unique(na.omit(data[, group])))
  message(paste0("Input is a cross-sectional model with ", nGroups,
                 " group", if (length(unique(data[, group])) > 1) "s", " and ",
                 length(lvs), " factor", if (length(lvs) > 1) "s fitted seperately", "."))

  ## checks
  if (any(unlist(lapply(items, FUN = length)) < 3)) {
    stop("Factors with less than 3 items are not identified.", call. = FALSE)
  }
  if (any(unlist(lapply(items, FUN = length)) == 3)) {
    warning("Configural MI for factors with 3 items cannot be separately tested.", call. = FALSE)
  }
  if (!(all(c(unlist(items), group) %in% colnames(data)))) {
    stop("Failed to find all variables in the data frame.", call. = FALSE)
  }
  if (dich) {
    dichCheck <- apply(na.omit(data[, unlist(items)]), 2, unique)
    if (is.list(dichCheck)) {
      stop("Items are not all dichotomous as declared. Problem with item ",
           names(which(lapply(dichCheck, length) > 2)), ".", call. = FALSE)
    }
  }

  ## missing handling
  if (any(is.na(data[, group]))) {
    message(paste0(sum(is.na(data[, group])), " cases deleted due to missing data in the covariate ",
                   group, ". Final sample: ", sum(!is.na(data[, group])), " subjects."))
    data <- data[!is.na(data[, group]), ]
  }
  missings <- "none"
  if (any(is.na(data[, unlist(items)]))) {
    if (!dich) {
      message("Missing data present in the items. Using FIML.")
      missings <- "fiml"
    }
    if (dich && dichModel == "factor") {
      message("Missing data present in the items. Using pairwise deletion under WLSMV. (Assumes missing completely at random! Alternatively, IRT models can be used.)")
      missings <- "pairwise"
    }
    if (dich && dichModel != "factor") {
      message("Missing data present in the items. Using FIML-type approach in IRT analysis.")
      missings <- "fiml"
    }
  }
  missingsOut <- missings # for printing information to the user in summary()

  data[, unlist(items)] <- data.matrix(data[, unlist(items)]) # assure numericly coded items

  if (stdItems && !dich && dichModel == "factor") {
    data[, unlist(items)] <- scale(data[, unlist(items)])   # standardization (cont items)
  }

  res <- list(model = model,
              data = data,
              group = group,
              lvs = lvs,
              items = items,
              settings = list(MI = MI,
                              dich = dich,
                              dichModel = dichModel,
                              stdItems = stdItems))

  ### fit models separately per factor
  for (lv in lvs) {

    ## lavaan models (metric items, dich items with non-IRT model)
    if (dichModel == "factor") {
      # switches
      estim <- "MLR"
      catItems <- NULL
      int <- "intercepts"
      if (dich) {
        estim <- "WLSMV"
        catItems <- items[[lv]]
        int <- "thresholds"
      }

      tryCatch({
        mod <- paste0(lv, "=~", paste(items[[lv]], collapse = " + "))
        modconfig <- paste0(paste0("F=~",
                                   " c(l", 1:length(items[[lv]]),
                                   "1, l",
                                   1:length(items[[lv]]),
                                   "2)*",
                                   items[[lv]],
                                   collapse = "\n"),
                            "\n") # adding labels for parameters to find them later more easily
        if (dich) {
          modconfig <- paste0(modconfig,
                              paste0(items[[lv]],
                                     " | c(a", 1:length(items[[lv]]),
                                     "1, a",
                                     1:length(items[[lv]]),
                                     "2)*t1",
                                     collapse = "\n"))
        } else {
          modconfig <- paste0(modconfig, paste0(items[[lv]],
                                                " ~ c(a", 1:length(items[[lv]]),
                                                "1, a", 1:length(items[[lv]]),
                                                "2)*1", collapse = "\n"))
        }

        if (missings == "none") {
          missings <- "listwise"   # set lavaan's default for missing data (only used if there are none!)
          missingsOut <- "none"
        }
        res[[lv]][["configuralModel"]] <- modconfig
        cfaconfig <- cfa(modconfig,
                         data = data,
                         estimator = estim,
                         group = group,
                         meanstructure = TRUE,
                         std.lv = TRUE,
                         ordered = catItems,
                         missing = missings)
        res[[lv]][["configural"]] <- list(cfaconfig)  # containing results in an extra list as older R version cannot handle it otherwise

        if (MI[[lv]] %in% c("weak", "strong", "strict")) {
          if (!dich) {
            cfaweak <- cfa(mod,
                           data = data,
                           estimator = estim,
                           group = group,
                           meanstructure = TRUE,
                           std.lv = TRUE,
                           ordered = catItems,
                           group.equal = c("loadings"),
                           missing = missings)
            res[[lv]][["weak"]] <- list(cfaweak)
          }
          if (MI[[lv]] %in% c("strong", "strict")) {
            cfastrong <- cfa(mod,
                             data = data,
                             estimator = estim,
                             group = group,
                             meanstructure = TRUE,
                             std.lv = TRUE,
                             ordered = catItems,
                             group.equal = c("loadings", int),
                             missing = missings)
            res[[lv]][["strong"]] <- list(cfastrong)

            if (MI[[lv]] == "strict") {
              cfastrict <- cfa(mod,
                               data = data,
                               estimator = estim,
                               group = group,
                               meanstructure = TRUE,
                               std.lv = TRUE,
                               ordered = catItems,
                               group.equal = c("loadings", int, "residuals"),
                               missing = missings)
              res[[lv]][["strict"]] <- list(cfastrict)
            }
          }
        }
      })
    }

    ## IRT models in mirt
    if (dichModel != "factor") {
      # checks
      if (MI[[lv]] == "strict") {
        stop("Strict MI is not relevant in IRT models.", call. = FALSE)
      }

      estim <- "EM"
      if (length(lvs) > 2) estim <- "MHRM"
      dat0 <- data[, unlist(items[[lv]])]  # prepare data for mirt
      if (is.numeric(data[, group])) { # deal with numerically coded groups for mirt
        data[, group] <- factor(data[, group],
                                levels = unique(data[, group]),
                                labels = paste0("group", unique(data[, group])))
        res$data <- data
      }

      cfaconfig <- multipleGroup(dat0,
                                 model = 1,
                                 itemtype = dichModel,
                                 group = data[, group],
                                 method = estim,
                                 SE = TRUE,
                                 verbose = FALSE)
      res[[lv]][["configural"]] <- list(cfaconfig)

      if (MI[[lv]] %in% c("weak", "strong")) {
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
        if (MI[[lv]] %in% c("strong")) {
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

  if (is.null(res[[lvs[1]]]$configural)) {
    stop("Something went wrong. Check your input.", call. = FALSE)
  }

  res$settings$missingsOut <- missingsOut
  res$settings$missings <- missings
  res$settings$estim <- estim
  class(res) <- "testMI"
  return(res)
}



#' @title Summary function for testMI
#' @rdname summary-testMI
#' @description Giving a summary of fit indices of an estimated testMI object per factor.
#' @export

summary.testMI <- function(res) {  # Tue May 25 08:57:08 2021 ------------------------------
  stopifnot(inherits(res, "testMI"))

  cat("Two group model with:\n")
  tab <- t(data.frame(unclass(table(res$data[[res$group]]))))
  tab <- cbind(tab, sum(tab[1, ]))
  rownames(tab) <- ""
  colnames(tab)[ncol(tab)] <- "total"
  print(tab)
  cat("\n")

  if (res$settings$dichModel == "factor") {
    ## lavaan
    Res <- data.frame(matrix(0, ncol = 14))

    inds <- c("chisq.scaled",
              "df.scaled",
              "pvalue.scaled",
              "cfi.scaled",
              "rmsea.scaled",
              "rmsea.ci.lower.scaled" ,
              "srmr")

    Tech <- data.frame(c("      lavaan",  # tech summary
                         res$settings$estim,
                         if (!res$settings$dich) "continuous" else "dichotomous",
                         res$settings$missingsOut,
                         if (!res$settings$dich) if (res$settings$stdItems) "yes" else "no"))
    rownames(Tech) <- c("package",
                        "estimator",
                        "item type",
                        "item missings",
                        if (!res$settings$dich) "standardized")
    colnames(Tech) <- ""

    for (lv in res$lvs) {
      le <- which(c("configural", "weak", "strong", "strict") %in% res$settings$MI[[lv]]) # remember row depth
      f <- which(res$lvs %in% lv)

      cnfg <- fitMeasures(res[[lv]]$configural[[1]])
      Res[(le*f - le + 1), 2:8] <- cnfg[inds]

      if (res$settings$MI[[lv]] %in% c("weak", "strong", "strict") && !res$settings$dich) {
        weak <- fitMeasures(res[[lv]]$weak[[1]])
        Res[(le*f - le + 2), 2:14] <- c(weak[inds],
                                        anova(res[[lv]]$configural[[1]],
                                              res[[lv]]$weak[[1]])[2, 5:7],
                                        weak["cfi.scaled"] - cnfg["cfi.scaled"],
                                        weak["rmsea.scaled"] - cnfg["rmsea.scaled"],
                                        weak["srmr"] - cnfg["srmr"])
      }
      if (res$settings$MI[[lv]] %in% c("weak", "strong", "strict") && res$settings$dich) {
        Res[(le*f - le + 2), 2:14] <- 0
      }
      if (res$settings$MI[[lv]] %in% c("strong", "strict") && !res$settings$dich) {
        strong <- fitMeasures(res[[lv]]$strong[[1]])
        Res[(le*f - le + 3), 2:14] <- c(strong[inds],
                                        anova(res[[lv]]$weak[[1]],
                                              res[[lv]]$strong[[1]])[2, 5:7],
                                        strong["cfi.scaled"] - weak["cfi.scaled"],
                                        strong["rmsea.scaled"] - weak["rmsea.scaled"],
                                        strong["srmr"] - weak["srmr"])
      }
      if (res$settings$MI[[lv]] %in% c("strong", "strict") && res$settings$dich) {
        strong <- fitMeasures(res[[lv]]$strong[[1]])
        Res[(le*f - le + 3), 2:14] <- c(strong[inds],
                                        anova(res[[lv]]$configural[[1]],
                                              res[[lv]]$strong[[1]])[2, 5:7],
                                        strong["cfi.scaled"] - cnfg["cfi.scaled"],
                                        strong["rmsea.scaled"] - cnfg["rmsea.scaled"],
                                        strong["srmr"] - cnfg["srmr"])
      }
      if (res$settings$MI[[lv]] %in% "strict") {
        strict <- fitMeasures(res[[lv]]$strict[[1]])
        Res[(le*f - le + 4), 2:14] <- c(strict[inds],
                                        anova(res[[lv]]$strong[[1]],
                                              res[[lv]]$strict[[1]])[2, 5:7],
                                        strict["cfi.scaled"] - strong["cfi.scaled"],
                                        strict["rmsea.scaled"] - strong["rmsea.scaled"],
                                        strict["srmr"] - strong["srmr"])
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
                       "diff SRMR")
    Res$chi2 <- round(Res$chi2, 2)
    Res[, 4:14] <- round(Res[, 4:14], 3)
    Res[, 1] <- rep(c("configural", "weak", "strong", "strict")[1:le],
                    length(res$lvs))
    Res[Res$'MI level' == "configural", 9:14] <- rep("", 6)
    if (res$settings$dich) {
      Res[Res$'MI level' == "weak", 2:14] <- rep("-", 13)
    }
    Res <- as.data.frame(t(Res))

    nam <- NULL
    for (lv in res$lvs) nam <- append(nam, c(lv, rep("", le - 1)))
    colnames(Res) <- nam
  } else {
    ## mirt
    Res <- data.frame(matrix(0, ncol = 15))
    Tech <- data.frame(c("     mirt",  # tech summary
                         res$settings$estim,
                         res$settings$dichModel,
                         res$settings$missingsOut))
    rownames(Tech) <- c("package",
                        "estimator",
                        "item type",
                        "item missings")
    colnames(Tech) <- ""

    for (lv in res$lvs) {
      le <- which(c("configural", "weak", "strong") %in% res$settings$MI[[lv]]) # remember row depth
      f <- which(res$lvs %in% lv)

      if (res$settings$missings == "none") {
        M2conf <- M2(res[[lv]]$configural[[1]])
      } else {
        M2conf <- rep(0, 10)
      }

      Res[(le*f - le + 1), 2:10] <- c(anova(res[[lv]]$configural[[1]])[c(1, 3, 5)],
                                      M2conf[c(1:5, 10)])

      if (res$settings$MI[[lv]] %in% c("weak", "strong") &&
          res$settings$dichModel != "Rasch") {

        if (res$settings$missings == "none") {
          M2weak <- M2(res[[lv]]$weak[[1]])
        } else {
          M2weak <- rep(0, 10)
        }

        Res[(le*f - le + 2), 2:15] <- c(anova(res[[lv]]$weak[[1]],
                                              res[[lv]]$configural[[1]],
                                              verbose = FALSE)[1, c(1, 3, 5)],
                                        M2weak[c(1:5, 10)],
                                        anova(res[[lv]]$weak[[1]],
                                              res[[lv]]$configural[[1]],
                                              verbose = FALSE)[2, 7:9],
                                        -(M2conf - M2weak)[c(4, 10)])
      }
      if (res$settings$MI[[lv]] %in% c("weak", "strong") &&
          res$settings$dichModel == "Rasch") {
        Res[(le*f - le + 2), 2:15] <- 0
      }
      if (res$settings$MI[[lv]] %in% "strong" &&
          res$settings$dichModel != "Rasch") {

        if (res$settings$missings == "none") {
          M2strong <- M2(res[[lv]]$strong[[1]])
        } else {
          M2strong <- rep(0, 10)
        }

        Res[(le*f - le + 3), 2:15] <- c(anova(res[[lv]]$strong[[1]],
                                              res[[lv]]$weak[[1]],
                                              verbose = FALSE)[1, c(1, 3, 5)],
                                        M2strong[c(1:5, 10)],
                                        anova(res[[lv]]$strong[[1]],
                                              res[[lv]]$weak[[1]],
                                              verbose = FALSE)[2, 7:9],
                                        -(M2weak - M2strong)[c(4, 10)])
      }
      if (res$settings$MI[[lv]] %in% "strong" &&
          res$settings$dichModel == "Rasch") {

        if (res$settings$missings == "none") {
          M2strong <- M2(res[[lv]]$strong[[1]])
        } else {
          M2strong <- rep(0, 10)
        }
        Res[(le*f - le + 3), 2:15] <- c(anova(res[[lv]]$strong[[1]],
                                              res[[lv]]$configural[[1]],
                                              verbose = FALSE)[1, c(1, 3, 5)],
                                        M2strong[c(1:5, 10)],
                                        anova(res[[lv]]$strong[[1]],
                                              res[[lv]]$configural[[1]],
                                              verbose = FALSE)[2, 7:9],
                                        -(M2conf - M2strong)[c(4, 10)])
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
                       "diff CFI")
    Res[, 2:15] <- round(Res[, 2:15], 3)
    Res[, 1] <- rep(c("configural", "weak", "strong")[1:le],
                    length(res$lvs))
    if (!(res$settings$missings == "none")) {
      Res[, c(4:10, 14:15)] <- "*"
    }
    Res[Res$'MI level' == "configural", 11:15] <- rep("", 5)
    if (res$settings$dichModel == "Rasch") {
      Res[Res$'MI level' == "weak", 2:15] <- rep("-", 14)
    }
    Res <- as.data.frame(t(Res))

    nam <- NULL
    for (lv in res$lvs) nam <- append(nam, c(lv, rep("", le - 1)))
    colnames(Res) <- nam
  }

  print(Res)
  print(Tech)

  if (!(res$settings$missings == "none") &&
      res$settings$dichModel != "factor") {
    cat("\n* Some fit statistics could not be calculated due to missing values. (Maydeu-Olivares & Joe, 2006)\n")
  }

  if (any(res$settings$MI %in% c("weak", "strong", "strict")) &&
      res$settings$dichModel == "factor" && res$settings$dich) {
    cat("\nWeak MI not testable with dichotomous items. (Wu & Estabrook, 2016)")
  }
  cat("\nUse getModel() to access parameter estimates or to further process the results.")
}

# Fri Jun 11 13:44:39 2021 ------------------------------

# internal function for parsing model input and delivering it
# to testMI() and clusterItems()

setUpModel <- function(items = NULL,
                       data,
                       group,
                       MItargetLevel = NULL,
                       stdItems = TRUE,
                       dich = FALSE,
                       dichModel = "factor") {

  stopifnot(dichModel %in% c("factor", "2PL", "Rasch"))

  ## parse model structure from lavaan - currently defunct
  #if (!is.null(model)) {
  #  pars0 <- lavParTable(model) # get basic model structure
  #  lvs <- unique(pars0[pars0$op == "=~", "lhs"]) # get factor names in model
  #} else {
  #  lvs <- "Factor"
  #}

  if (is.list(items)) {
    if (all(!sapply(items, is.character))) { # change listing of items for data frame entries
      for (i in names(items)[!sapply(items, is.character)]) {
        items[[i]] <- colnames(items[[i]])
      }
    }
  } else { # if items are given as vector for a single factor
    if (!is.character(items)) { # change listing of items for data frame entries
      items <- list(Factor = colnames(items))
    } else {
      items <- list(Factor = items)
    }
  }

  lvs <- names(items)

  #if (!(MIlevel %in% levels)) {  # only needed for old model command
  #  lvs <- lvs[sapply(lvs, FUN = grepl, MIlevel)] # filter for factors that are actually mentioned in MIlevel statement (leaving out factors not of interest)
  #}

  #if (any(sapply(items, length) == 1)) { # check for items given by the shortcut (e.g. item1-item10) . UNDER DEV
  #  for (i in names(items)[sapply(items, length) == 1]) {
  #itemSections <- unlist(strsplit(items$Speed, "[ ]+"))
  # itemLists <- unlist(strsplit(itemSections, "[-]+"))
  # regmatches(ert, regexpr("^[^0-9]*", ert))
  # regmatches(ert, regexpr("[0-9]*$", ert), invert = TRUE)
  # }
  # }

  #nGroups <- length(unique(na.omit(data[, group])))         # messaging about model input. Currently defunct.
  #message(paste0("Input is a cross-sectional model with ", nGroups,
  #               " group", if (length(unique(data[, group])) > 1) "s", " and ",
  #               length(lvs), " factor", if (length(lvs) > 1) "s fitted seperately", "."))

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

  ## parsing target MI levels
  levels <- c("configural",
              "weak",
              "strong",
              "strict")
  MI <- list()  # internal representation of input MI levels per factor
  lvsPlaces <- lapply(lvs, FUN = regexpr, MItargetLevel) # parse input string
  lvsPlaces <- append(lvsPlaces, nchar(MItargetLevel) + 1) # add final stop in order to make loop below work
  if (any(lvsPlaces < 0) && !(MItargetLevel %in% levels) &&
      nchar(MItargetLevel) > (sum(nchar(lvs)) + 5*length(lvs))) {
    stop("Error while reading MItargetLevel. Couldn't retrieve factor names.", call. = FALSE)   # check for misspelled factor names
  }
  for (lv in lvs) {
    if (MItargetLevel %in% levels) {
      MI[[lv]] <- MItargetLevel # catch cases, where a single MI level is stated for all factors
    } else {
      MI[[lv]] <- levels[which(lapply(levels,   # get desired MI level per factor
                                      FUN = regexpr,
                                      substring(MItargetLevel,
                                                lvsPlaces[[which(lvs %in% lv)]] +
                                                  attr(lvsPlaces[[which(lvs %in% lv)]], "match.length"),
                                                lvsPlaces[[which(lvs %in% lv) + 1]] - 1)) > 0)]
    }
  }
  if (length(MI) == 0) {
    stop("Error while reading MItargetLevel Check your statements.", call. = FALSE)
    if (any(apply(MI, FUN = length) == 0)) {
      stop("Error while reading MItargetLevel Check your statements.", call. = FALSE)
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
      #  message("Missing data present in the items. Using FIML.")
      missings <- "fiml"
    }
    if (dich && dichModel == "factor") {
      message("Missing data present in the items. Using pairwise deletion under WLSMV. (Assumes missing completely at random! Alternatively, IRT models can be used.)")
      missings <- "pairwise"
    }
    if (dich && dichModel != "factor") {
      # message("Missing data present in the items. Using FIML-type approach in IRT analysis.")
      missings <- "fiml"
    }
  }
  missingsOut <- missings # for printing information to the user in summary()

  data[, unlist(items)] <- data.matrix(data[, unlist(items)]) # assure numericly coded items

  if (stdItems && !dich && dichModel == "factor") {
    data[, unlist(items)] <- scale(data[, unlist(items)])   # standardization (cont items)
  }

  res <- list(data = data,
              model = list(items = items,
                           group = group,
                           dich = dich,
                           dichModel = dichModel,
                           stdItems = stdItems))

  ### fit models separately per factor
  for (lv in lvs) {
    res[[lv]] <- list(MItargetLevel = MI[[lv]])

# lavaan models (metric items, dich items with non-IRT model)
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
        modconfig <- paste0(paste0("Factor=~",
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


      })
    }
# IRT models in mirt
    if (dichModel != "factor") {
      if (MItargetLevel == "strict") {
        stop("Strict MI is not sensible in IRT models.", call. = FALSE)
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

    }
  }

  if (is.null(res[[lvs[1]]]$configural)) {
    stop("Something went wrong. Check your input.", call. = FALSE)
  }

  res$model$missingsOut <- missingsOut
  res$model$missings <- missings
  res$model$estim <- estim

  return(res)
}


## Loading
setwd("X:/4$e/FU/# R Package")
library(testthat)
library(MBESS)
data(HS)
DataCont <- HS[, c("t10_addition", "t11_code", "t12_counting_groups_of_dots", # grab the "speed" items (see help file in MBESS package)
                   "t13_straight_and_curved_capitals",
                   "t20_deduction", "t21_numerical_puzzles", "t22_problem_reasoning", # math items
                   "t23_series_completion", "t24_woody_mccall", "sex")]
colnames(DataCont) <- sub("_.*", "", colnames(DataCont)) # shorten variable names for MPlus' convenience

Daten <- foreign::read.spss("data 3.sav", use.value.labels = F, to.data.frame = TRUE)
scales <- list(Agressive_Impulses = c("YSC001T0","YSC002T0","YSC003T0","YSC004T0","YSC005T0","YSC006T0"),
               Cleanliness = c("YSC009T0","YSC010T0","YSC011T0","YSC012T0","YSC013T0","YSC014T0","YSC015T0","YSC016T0","YSC036T0","YSC037T0","YSC038T0","YSC039T0"),
               Responsibility = c("YSC007T0","YSC008T0","YSC040T0","YSC041T0","YSC042T0","YSC043T0","YSC044T0","YBOCS_T0_565758"),
               Keeping_Order = c("YSC024T0","YSC030T0","YSC049T0","YSC052T0"),
               Magical_Thinking = c("YSC021T0","YSC023T0","YSC033T0","YSC034T0","YSC035T0","YSC060T0"),
               Pure_Repetitions = c("YSC046T0","YSC047T0","YSC048T0","YSC054T0","YSC055T0","YSC059T0"),
               Mental_Urges = c("YSC027T0","YSC031T0","YSC032T0","YSC051T0","YSC053T0"),
               small_scales = c("YSC017T0","YSC018T0","YSC019T0", "YSC025T0","YSC026T0","YSC045T0",
                                "YSC022T0","YSC028T0","YSC029T0"))
data <- Daten[, c(scales[[1]], scales[[7]])]
data[data <= 1] <- 0  # current
data[data != 0] <- 1
gender <- Daten$GENDER
gender[gender == 0] <- NA
gender[gender == 1] <- "m"
gender[gender == 2] <- "w"
data$gender <- gender
DataDich <- data


library(TAM)
data("data.fims.Aus.Jpn.scored")
dataDich2 <- data.fims.Aus.Jpn.scored[, -13]


#############
test_that("testMI, Cont, single MIlevel", {
  res <- testMI(items = c("t20", "t21", "t22", "t23", "t24"),
                group = "sex",
                data = DataCont,
                MIlevel = "strong")
  expect_equal(unname(fitmeasures(res$Factor$configural[[1]])[6]),
               17.26,
               tolerance = 0.1)
  expect_output(summary(res),
                regexp = "total")
  summary(res)
})

DataCont[1, 1] <- NA

#############
test_that("testMI, Missings, Cont, single MIlevel", {
  res <- testMI(items = c("t10", "t11", "t12", "t13"),
                group = "sex",
                data = DataCont,
                MIlevel = "strong")
  expect_equal(unname(fitmeasures(res$Factor$configural[[1]])[6]),
               16.83,
               tolerance = 0.1)
  expect_output(summary(res),
                regexp = "fiml")
  summary(res)
})


res <- testMI(items = list(Speed = c("t10", "t11", "t12", "t13"),
                           Math = c("t20", "t21", "t22", "t23", "t24")),
              group = "sex",
              data = DataCont,
              MItargetLevel = "Speed strong Math weak")

################
test_that("testMI, Cont, Missings, normal different MIlevel", {
  expect_equal(unname(fitmeasures(res$Speed$configural[[1]])[6]),
               16.83,
               tolerance = 0.1)
  expect_output(summary(res),
                regexp = "strong")
  expect_output(summary(res),
                regexp = "fiml")
  summary(res)
})


#############
test_that("clusterItems,Cont, res_testMI, Missings, single MIholding, threshold", {
  res2 <- clusterItems(res_testMI = res,
                       MIholdingLevel = "configural",
                       method = "threshold",
                       loadThreshold = 0.3,
                       intThreshold = 0.2)
  expect_equal(max(res2$Speed$itemClustering$finalClustering),
                   4)
  expect_output(summary(res2),
                regexp = "configural -> weak")
  expect_output(summary(res2),
                regexp = "t13")
  expect_output(summary(res2),
                regexp = "t20")
  expect_output(summary(res2),
                regexp = "Clustering by threshold criterion with load threshold 0.3 and intercept threshold 0.2.")
  summary(res2)
})


#############
test_that("clusterItems, Cont, res_testMI, Missings, normal MIholding, threshold", {
  res2 <- clusterItems(res_testMI = res,
                       MIholdingLevel = "Speed weak Math configural",
                       method = "threshold",
                       loadThreshold = 0.3,
                       intThreshold = 0.2)
  expect_equal(max(res2$Speed$itemClustering$finalClustering),
               3)
  expect_output(summary(res2,
                        order = "items"),
               regexp = "t10       2
t11       1")
  summary(res2,
          order = "items")
})


#############
test_that("clusterItems, Cont, res_testMI, Missings, MIholding with MI factor, threshold", {
  expect_message({res2 <- clusterItems(res_testMI = res,
                       MIholdingLevel = "Speed configural Math weak",
                       method = "threshold",
                       loadThreshold = 0.3,
                       intThreshold = 0.2)},
                 "Skipping clustering for Math as MI level already holds.")
  expect_output(summary(res2),
                regexp = "not clustered")
  expect_output(summary(res2),
                regexp = "t20       1
t21       1
t22       1
t23       1
t24       1")
  summary(res2)
})


#############
test_that("clusterItems, Cont, new model, Missings, MIholding with MI factor, sigTest", {
  res2 <- clusterItems(MIholdingLevel = "Speed configural Math weak",
                       method = "sigTest",
                       pValue = 0.05,
                       items = list(Speed = c("t10", "t11", "t12", "t13"),
                                    Math = c("t20", "t21", "t22", "t23", "t24")),
                       group = "sex",
                       data = Data,
                       MItargetLevel = "strong")
  expect_equal(res2$Speed$itemClustering$clusterLStep$clusteringHistory[[1]]$test[3],
               0.029,
               tolerance = 0.001)
  expect_output(summary(res2),
                regexp = "t24       1
t20       2
t21       2
t22       2
t23       2")
  summary(res2)
})



###################
test_that("testMI, Dich, catSEM, single MIlevel, giving df section", {
  expect_message({res <- testMI(items = list(AI = DataDich[1:6],
                                             MU = DataDich[7:11]),
                                data = DataDich,
                                group = "gender",
                                MItargetLevel = "strong",
                                dich = T,
                                dichModel = "factor")},
                 "18 cases deleted due to missing data in the covariate gender. Final sample: 758 subjects.")
  expect_equal(unname(fitmeasures(res$MU$strong[[1]])[6]),
               48.63,
               tolerance = 0.1)

  summary(res)
})

res <- testMI(items = list(AI = DataDich[1:6],
                           MU = DataDich[7:11]),
              data = DataDich,
              group = "gender",
              MItargetLevel = "strong",
              dich = T,
              dichModel = "factor")


###############
test_that("clusterItems, Dich, catSEM, threshold, single MIlevel", {
  res2 <- clusterItems(res,
                       MIholding = "configural",
                       method = "threshold",
                       loadThreshold = 0.3,
                       intThreshold = 0.4)
  summary(res2)
})


###############
cat("\n\nclusterItems, Dich, catSEM, sigTest, single MIlevel")
res2 <- clusterItems(res,
                     MIholding = "configural",
                     method = "sigTest",
                     pValue = 0.05)
summary(res2)

################
cat("\n\ntestMI, Dich, Missings, catSEM, single MIlevel")
DataDich[1, 1] <- NA
res <- testMI(model,
              data = DataDich,
              group = "gender",
              MIlevel = "strong",
              dich = T,
              dichModel = "factor")
summary(res)

#################
cat("\n\ntestMI, Dich, Missings, 2PL, single MIlevel")
res <- testMI(model,
              data = DataDich,
              group = "gender",
              MIlevel = "strong",
              dich = T,
              dichModel = "Rasch")
summary(res)

#################
cat("\n\ntestMI, Dich, Missings, 2PL, single MIlevel")
res <- testMI(model,
              data = DataDich,
              group = "gender",
              MIlevel = "strong",
              dich = T,
              dichModel = "2PL")
summary(res)

###############
cat("\n\nclusterItems, Dich, Missings, 2PL, sigTest, single MIlevel")
res2 <- clusterItems(res,
                     MIholding = "configural",
                     method = "sigTest",
                     pValue = 0.05)
summary(res2)

###############
cat("\n\nclusterItems, Dich, Missings, Rasch, sigTest, single MIlevel")
res2 <- clusterItems(res,
                     MIholding = "configural",
                     method = "sigTest",
                     pValue = 0.05)
summary(res2)

###############
cat("\n\nclusterItems, Dich, Missings, 2PL, thresh, single MIlevel")
res2 <- clusterItems(res,
                     MIholding = "configural",
                     method = "threshold",
                     loadThreshold = 0.3,
                     intThreshold = 0.4)
summary(res2)


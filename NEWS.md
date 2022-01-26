# measurementInvariance 0.3.1-180921

- Modifying Bayesian model averaging to fit procedure as described in 2022 SEM Paper

# measurementInvariance 0.3.1-180921

- overhauling the internal structure, thus changing some function arguments
- functions work now stand alone
- adding a plot function for the DIF clustering
- adding a print function for the DIF values


# measurementInvariance 0.3.0-300521

- added modelAveraging function for 2PL, Rasch and metric items [for IRT models, these Bayesian estimations are not yet done via precompiled STAN code. Thus they are currently rather slow and crash-prone. Should STAN crash R, simply try it again, it will probably work out.]
- added plotAverage function for displaying results of the above
- extended getModel in order to process Bayesian models
- included the summary functions into the help files
- some bugfixes in partialMI


# measurementInvariance 0.2.0-250521

First dev version on gitHub.

Features of the core functions:

- testMI: Global test of MI. Multidimensional measurement models for continuous & dichotomous items. The latter can be treated by categorical SEM or IRT models (Rasch or 2PL). MI covariate is a multiple group (2 or more) setting. MIlevels supported: configural, weak, strong, strict.
- clusterItems: Finding subsets of items. Multidimensional measurement models for continuous & dichotomous items. The latter can be treated by categorical SEM or IRT models (Rasch or 2PL). MI covariate is a 2 group setting. MI levels supported: weak, strong. Clustering techniques: kMeans via threshold or significance test.
- partialMI: Estimate single partial MI models. Multidimensional measurement models for continuous & dichotomous items. The latter can be treated by categorical SEM or IRT models (Rasch or 2PL). MI covariate is a 2 group setting. MI levels supported: weak, strong.
- bayModAveraging: Apply Bayesian model averaging. Unidimensional 2PL model. MI covariate is a 2 group setting. MI levels supported: weak, strong.


* Added a `NEWS.md` file to track changes to the package.

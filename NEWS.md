# measurementInvariance 0.2.0-250521

First dev version on gitHub.

Features of the core functions:

- testMI: Global test of MI. Multidimensional measurement models for continuous & dichotomous items. The latter can be treated by categorical SEM or IRT models (Rasch or 2PL). MI covariate is a multiple group (2 or more) setting. MIlevels supported: configural, weak, strong, strict.
- clusterItems: Finding subsets of items. Multidimensional measurement models for continuous & dichotomous items. The latter can be treated by categorical SEM or IRT models (Rasch or 2PL). MI covariate is a 2 group setting. MI levels supported: weak, strong. Clustering techniques: kMeans via threshold or significance test.
- partialMI: Estimate single partial MI models. Multidimensional measurement models for continuous & dichotomous items. The latter can be treated by categorical SEM or IRT models (Rasch or 2PL). MI covariate is a 2 group setting. MI levels supported: weak, strong.
- bayModAveraging: Apply Bayesian model averaging. Unidimensional 2PL model. MI covariate is a 2 group setting. MI levels supported: weak, strong.


* Added a `NEWS.md` file to track changes to the package.

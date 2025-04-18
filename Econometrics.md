---
name: Econometrics
topic: Econometrics
maintainer: Achim Zeileis, Grant McDermott, Kevin Tappe
email: Achim.Zeileis@R-project.org
version: 2025-04-03
source: https://github.com/cran-task-views/Econometrics/
---

Base R ships with a lot of functionality useful for (computational) econometrics,
in particular in the stats package. This functionality is complemented by many
packages on CRAN, a brief overview is given below. There is also a certain
overlap between the tools for econometrics in this view and those in the task
views on `r view("Finance")`, `r view("TimeSeries")`, and `r view("CausalInference")`.

The packages in this view can be roughly structured into the following topics.
If you think that some package is missing from the list, please file an issue in
the GitHub repository or contact the maintainer.


### Basic linear regression

- *Estimation and standard inference:* Ordinary least squares (OLS)
  estimation for linear models is provided by `lm()` (from stats) and standard
  tests for model comparisons are available in various methods such as
  `summary()` and `anova()`.
- *Further inference and nested model comparisons:* Functions analogous to
  the basic `summary()` and `anova()` methods that also support asymptotic
  tests (*z* instead of *t* tests, and Chi-squared instead of *F* tests) and
  plug-in of other covariance matrices are `coeftest()` and `waldtest()` in
  `r pkg("lmtest", priority = "core")`. (Non)linear hypothesis testing for a
  wide range of R packages can implemented through the `deltamethod()` 
  function of `r pkg("marginaleffects", priority = "core")`. This expands on
  older (non)linear hypothesis test functions like `linearHypothesis()` and
  `deltaMethod()` from `r pkg("car", priority = "core")`.
- *Robust standard errors:* HC, HAC, clustered, and bootstrap covariance
  matrices are available in `r pkg("sandwich", priority = "core")` and can be
  plugged into the inference functions mentioned above.
- *Nonnested model comparisons:* Various tests for comparing non-nested
  linear models are available in `r pkg("lmtest")` (encompassing test, J test,
  Cox test). The Vuong test for comparing other non-nested models is provided
  by `r pkg("nonnest2")` (and specifically for count data regression in
  `r pkg("pscl")`).
- *Diagnostic checking:* The packages `r pkg("car")` and `r pkg("lmtest")`
  provide a large collection of regression diagnostics and diagnostic tests.
- *Miscellaneous:* Much of the above functionality is bundled together in
  `r pkg("fixest", priority = "core")`, which provides a number of in-built
  convenience features that users may find attractive. This includes 
  robust standard error specification, multi-model estimation, custom 
  hypothesis testing, etc.


### Microeconometrics

- *Generalized linear models (GLMs):* Many standard microeconometric models
  belong to the family of generalized linear models and can be fitted by
  `glm()` from package stats. This includes in particular logit and probit
  models for modeling choice data and Poisson models for count data.
- *Effects and marginal effects:*  Effects for typical values of regressors
  in GLMs and various other probabilistic regression models can be obtained
  and visualized using `r pkg("effects")`. Marginal effect tables and
  corresponding visualizations for a wide range of models can be be produced
  with `r pkg("marginaleffects", priority = "core")`. Other implementations
  of marginal effects for certain models are in `r pkg("margins")` and
  `r pkg("mfx")`. Interactive visualizations of both effects and marginal
  effects are possible in `r pkg("LinRegInteractive")`.
- *Binary responses:* The standard logit and probit models (among many
  others) for binary responses are GLMs that can be estimated by `glm()` with
  `family = binomial`. Bias-reduced GLMs that are robust to complete and
  quasi-complete separation are provided by `r pkg("brglm")`. Discrete choice
  models estimated by simulated maximum likelihood are implemented in
  `r pkg("Rchoice")`. `r pkg("bife")` provides binary choice models with fixed
  effects. Heteroscedastic probit models (and other heteroscedastic GLMs) are
  implemented in `r pkg("glmx")` along with parametric link functions and
  goodness-of-link tests for GLMs.
- *Count responses:* The basic Poisson regression is a GLM that can be
  estimated by `glm()` with `family = poisson` as explained above. Negative
  binomial GLMs are available via `glm.nb()` in package `r pkg("MASS")`.
  Another implementation of negative binomial models is provided by
  `r pkg("aod")`, which also contains other models for overdispersed data.
  Zero-inflated and hurdle count models are provided in package `r pkg("pscl")`.
  A reimplementation by the same authors is currently under
  development in `r rforge("countreg")` on R-Forge which also encompasses
  separate functions for zero-truncated regression, finite mixture models etc.
- *Multinomial responses:* Multinomial models with individual-specific
  covariates only are available in `multinom()` from package `r pkg("nnet")`.
  An implementation with both individual- and choice-specific variables is
  `r pkg("mlogit")`. Generalized multinomial logit models
  (e.g., with random effects etc.) are in `r pkg("gmnl")`. A flexible
  framework of various customizable choice models (including multinomial logit
  and nested logit among many others) is implemented in the `r pkg("apollo")`
  package. The newer `r pkg("logitr")` package combines many of the features
  from these preceding packages and also offers some meaningful performance
  improvements for fast estimation of multinomial and mixed logit models.
  Simulated maximum likelihood estimation of mixed logit models, especially
  for large data sets, is available in `r pkg("mixl")`.
  Generalized additive models (GAMs) for multinomial responses can be
  fitted with the `r pkg("VGAM")` package. A Bayesian approach to multinomial
  probit models is provided by `r pkg("MNP")`. Various Bayesian multinomial
  models (including logit and probit) are available in `r pkg("bayesm")`.
  The package `r pkg("RSGHB")` fits various hierarchical Bayesian
  specifications based on direct specification of the likelihood function.
  Furthermore, the `r pkg ("RprobitB")` package implements latent class mixed
  multinomial probit models for approximations of the true underlying mixing
  distribution.
- *Ordered responses:* Proportional-odds regression for ordered responses is
  implemented in `polr()` from package `r pkg("MASS")`. The package
  `r pkg("ordinal")` provides cumulative link models for ordered data which
  encompasses proportional odds models but also includes more general
  specifications. Bayesian ordered probit models are provided by
  `r pkg("bayesm")` and `r pkg("RprobitB")`.
- *Censored responses:* Basic censored regression models (e.g., tobit models)
  can be fitted by `survreg()` in `r pkg("survival")`, a convenience interface
  `tobit()` is in package `r pkg("AER", priority = "core")`. Further censored
  regression models, including models for panel data, are provided in
  `r pkg("censReg")`. Censored regression models with conditional
  heteroscedasticity are in `r pkg("crch")`. Furthermore, hurdle models for
  left-censored data at zero can be estimated with `r pkg("mhurdle")`. Models
  for sample selection are available in `r pkg("sampleSelection")` and
  `r pkg("ssmrob")` using classical and robust inference, respectively.
  Package `r pkg("matchingMarkets")` corrects for selection bias when the sample
  is the result of a stable matching process (e.g., a group formation or college
  admissions problem).
- *Truncated responses:* `r pkg("crch")` for truncated (and potentially
  heteroscedastic) Gaussian, logistic, and t responses. Homoscedastic Gaussian
  responses are also available in `r pkg("truncreg")`.
- *Fraction and proportion responses:* Beta regression for responses in (0, 1) is in
  `r pkg("betareg")` and `r pkg("gamlss")`.
- *Duration responses:* Many classical duration models can be fitted with
  `r pkg("survival")`, e.g., Cox proportional hazard models with `coxph()` or
  Weibull models with `survreg()`. Many more refined models can be found in
  the `r view("Survival")` task view.
- *High-dimensional fixed effects:* Linear and generalized linear models with
  potentially high-dimensional fixed effects, also for multiple groups, can be
  fitted with `r pkg("fixest", priority = "core")`, using optimized parallel
  C++ code. Other implementations of high-dimensional fixed effects are in
  `r pkg("lfe")` and `r pkg("alpaca")` for linear and generalized linear models,
  respectively.
- *Miscellaneous:* Further more refined tools for microeconometrics are
  provided in the `r pkg("micEcon")` family of packages: Analysis with
  Cobb-Douglas, translog, and quadratic functions is in `r pkg("micEcon")`;
  the constant elasticity of scale (CES) function is in `r pkg("micEconCES")`;
  the symmetric normalized quadratic profit (SNQP) function is in
  `r pkg("micEconSNQP")`. The almost ideal demand system (AIDS) is in
  `r pkg("micEconAids")`. Stochastic frontier analysis (SFA) is in
  `r pkg("frontier")`.
  Semiparametric SFA in is available in `r pkg("semsfa")` and spatial SFA in
  `r pkg("ssfa")`. The package `r pkg("bayesm")`
  implements a Bayesian approach to microeconometrics and marketing. Inference
  for relative distributions is contained in package `r pkg("reldist")`.


### Common research designs for causal inference

We review packages related to some common research designs for causal
inference below. This section is necessarily brief and should be paired with
the `r view("CausalInference")` task view, since is there a high degree of
overlap.

#### Difference-in-differences and synthetic control

- *Basic difference-in-differences (DiD):* The canonical 2x2 DiD model (two 
  units, two periods) can be estimated as a simple interaction between two
  factor variables in `lm()` or `glm()`, etc. Similarly, the equivalent 
  two-way fixed effects (TWFE) design can be obtained using factors to control 
  for unit and time fixed effects. However, for high-dimensional datasets TWFE 
  is more conveniently estimated using a dedicated panel data package like 
  `r pkg("fixest")` or `r pkg("plm")`. The former even provides a convenience 
  `i()` operator for constructing and interacting factors in TWFE settings.
- *Advanced DiD and TWFE corrections:* Despite its long-standing popularity,
  recent research has uncovered various problems with (naive) TWFE; for example,
  severe bias in the presence of staggered treatment rollout. A [cottage
  industry](https://asjadnaqvi.github.io/DiD/docs/02_R/) of workarounds and 
  alternative estimators now exists to address these problems. R package 
  implementations include: `r pkg("bacondecomp")`, `r pkg("did")`,
  `r pkg("did2s")`, `r pkg("DRDID")`, `r pkg("etwfe")`, `r pkg("fixest")` (via
  the `sunab()` function), and `r pkg("gsynth")`.
- *Synthetic control:* The original synthetic control (SC) implementation is
  available through `r pkg("Synth")`, while `r pkg("tidysynth")` offers a newer
  SC implementation with various enhancements (speed, inspection, etc.) 
  Similarly, `r pkg("gsynth")` generalizes the original SC implementation to
  multiple treated units and variable treatment periods, and also supports
  additional estimation methods like the EM algorithm and matrix completion.

#### Instrumental variables

- *Basic instrumental variables (IV) regression:* Two-stage least squares
  (2SLS) is provided by `r pkg("ivreg", priority = "core")`, which separates
  out the dedicated 2SLS routines previously found in `r pkg("AER")`. Another
  implementation is available as `tsls()` in package `r pkg("sem")`.
- *Binary responses:* The `r pkg("LARF")` package estimates local average
  response functions for binary treatments and binary instruments.
- *Panel data:* Several panel data model packages (see below) provide their own
  dedicated IV routines for efficient estimation in the presence of
  high-dimensional data. These include `r pkg("fixest")` and `r pkg("lfe")` for
  fixed effects, and `r pkg("plm")` for first-difference, between, and multiple
  random effects methods.
- *Miscellaneous:* `r pkg("REndo")` fits linear models with endogenous
  regressor using various latent instrumental variable approaches. 
  `r pkg("SteinIV")` provides semi-parametric IV estimators, including JIVE and
  SPS.

#### Regression discontinuity design

- Regression discontinuity design (RDD) methods are implemented in 
  `r pkg("rdrobust")` (offering robust confidence interval construction and
  bandwidth selection), `r pkg("rddensity")` (density discontinuity testing
  (also known as manipulation testing)), `r pkg("rdlocrand")` (inference under 
  local randomization), and `r pkg("rdmulti")` (analysis with multiple cutoffs 
  or scores).
- Tools to perform power, sample size and minimum detectable effects (MDE) 
  calculations are available in `r pkg("rdpower")`, while `r pkg("RATest")` 
  provides a collection of randomization tests, including a permutation test
  for the continuity assumption of the baseline covariates in the sharp RDD.


### Panel data models

- *Panel standard errors:* A simple approach for panel data is to fit the
  pooling (or independence) model (e.g., via `lm()` or `glm()`) and only
  correct the standard errors. Different types of clustered, panel, and
  panel-corrected standard errors are available in `r pkg("sandwich")`
  (incorporating prior work from `r pkg("multiwayvcov")`),
  `r pkg("clusterSEs")`, `r pkg("pcse")`, `r pkg("clubSandwich")`, 
  `r pkg("plm", priority = "core")`, and `r pkg("geepack")`, respectively. 
  The latter two require estimation of the pooling/independence models via
  `plm()` and `geeglm()` from the respective packages (which also provide 
  other types of models, see below).
- *Linear panel models:* `r pkg("fixest", priority = "core")` provides very
  efficient fixed-effect routines that scale to high-dimensional data and
  multiple fixed-effects.  `r pkg("plm")`, providing a wide range of within,
  between, and random-effect methods (among others) along with corrected
  standard errors, tests, etc. Various dynamic panel models are
  available in `r pkg("plm")`, with estimation based on moment conditions in
  `r pkg("pdynmc")`, and dynamic panel models with fixed effects in
  `r pkg("OrthoPanels")`. `r pkg("feisr")` provides fixed effects individual
  slope (FEIS) models. Panel vector autoregressions are implemented in
  `r pkg("panelvar")`.
- *GLMs and generalized estimation equations*. The aformentioned `r pkg("fixest")`
  supports a variety of GLM-like models in addition to linear panel models. 
  This includes efficient fixed-effect estimation of logit, probit, Poisson,
  and negative binomial models. Similar functionality is provided by 
  `r pkg("alpaca")` (which also accounts for incidental parameter problems) 
  and `r pkg("pglm")`. `r pkg("penppml")` further extends the high-dimensional
  case through penalized Poisson Pseudo Maximum Likelihood (PPML) regressions,
  using lasso or ridge penalties. GEE models for panel data (or longitudinal
  data in statistical jargon) are available in in `r pkg("geepack")`.
- *Mixed effects models:* Linear and nonlinear models for panel data (and
  more general multi-level data) are available in `r pkg("lme4")` and `r pkg("nlme")`.
- *Instrumental variables:* `r pkg("fixest")`. See also above.
- *Miscellaneous:* Threshold regression
  and unit root tests are in `r pkg("pdR")`. The panel data approach method
  for program evaluation is available in `r pkg("pampe")`. Dedicated fast data
  preprocessing for panel data econometrics is provided by `r pkg("collapse")`.


### Further regression models

- *Nonlinear least squares modeling:* `nls()` in package stats.
- *Quantile regression:* `r pkg("quantreg")` (including linear, nonlinear,
  censored, locally polynomial and additive quantile regressions).
- *Generalized method of moments (GMM) and generalized empirical likelihood
  (GEL):* `r pkg("gmm")`.
- *Spatial econometric models:* The `r view("Spatial")` view gives details
  about handling spatial data, along with information about (regression)
  modeling. In particular, spatial regression models can be fitted using
  `r pkg("spatialreg")` and `r pkg("sphet")` (the latter using a GMM approach).
  `r pkg("splm")` is a package for spatial panel models. Spatial probit models
  are available in `r pkg("spatialprobit")` and spatial seemingly unrelated
  regression (SUR) models in `r pkg("spsur")`.
- *Bayesian model averaging (BMA):* A comprehensive toolbox for BMA is
  provided by `r pkg("BMS")` including flexible prior selection, sampling,
  etc. A different implementation is in `r pkg("BMA")` for linear models,
  generalizable linear models and survival models (Cox regression).
- *Linear structural equation models:* `r pkg("lavaan")` and `r pkg("sem")`.
  See also the `r view("Psychometrics")` task view for more details.
- *Machine learning:* There are several packages that combine machine
  learning techniques with econometric inference (especially for identifying
  causal effects). These include `r pkg("grf")` for causal random forests
  and estimation of heterogeneous treatment effects, `r pkg("DoubleML")`
  for double machine learning of a wide range of models from the mlr3 family,
  and `r pkg("hdm")` for selected high-dimensional econometric models.
  For a more general overview see the `r view("MachineLearning")` task view.
- *Simultaneous equation estimation:* `r pkg("systemfit")`.
- *Nonparametric methods:* `r pkg("np")` using kernel smoothing and
  `r pkg("NNS")` using partial moments.
- *Linear and nonlinear mixed-effect models:* `r pkg("nlme")` and
  `r pkg("lme4")`.
- *Generalized additive models (GAMs):* `r pkg("mgcv")`, `r pkg("gam")`,
  `r pkg("gamlss")` and `r pkg("VGAM")`.
- *Design-based inference:* `r pkg("estimatr")` contains fast procedures for
  several design-appropriate estimators with robust standard errors and
  confidence intervals including linear regression, instrumental variables
  regression, difference-in-means, among others.
- *Extreme bounds analysis:* `r pkg("ExtremeBounds")`.
- *Miscellaneous:* The packages `r pkg("VGAM")`, `r pkg("rms")` and
  `r pkg("Hmisc")` provide several tools for extended handling of
  (generalized) linear regression models.


### Time series data and models

- The `r view("TimeSeries")` task view provides much more detailed information
  about both basic time series infrastructure and time series models. Here,
  only the most important aspects relating to econometrics are briefly
  mentioned. Time series models for financial econometrics (e.g., GARCH,
  stochastic volatility models, or stochastic differential equations, etc.)
  are described in the `r view("Finance")` task view.
- *Infrastructure for regularly spaced time series:* The class `"ts"` in
  package stats is R's standard class for regularly spaced time series
  (especially annual, quarterly, and monthly data). It can be coerced back and
  forth without loss of information to `"zooreg"` from package `r pkg("zoo",
  priority = "core")`.
- *Infrastructure for irregularly spaced time series:* `r pkg("zoo")`
  provides infrastructure for both regularly and irregularly spaced time
  series (the latter via the class `"zoo"`) where the time information can be
  of arbitrary class. This includes daily series (typically with `"Date"` time
  index) or intra-day series (e.g., with `"POSIXct"` time index). An extension
  based on `r pkg("zoo")` geared towards time series with different kinds of
  time index is `r pkg("xts")`. Further packages aimed particularly at finance
  applications are discussed in the `r view("Finance")` task view.
- *Classical time series models:* Simple autoregressive models can be
  estimated with `ar()` and ARIMA modeling and Box-Jenkins-type analysis can
  be carried out with `arima()` (both in the stats package). An enhanced
  version of `arima()` is in `r pkg("forecast", priority = "core")`.
- *Linear regression models:* A convenience interface to `lm()` for
  estimating OLS and 2SLS models based on time series data is `r pkg("dynlm")`.
  Linear regression models with AR error terms via GLS is
  possible using `gls()` from `r pkg("nlme")`.
- *Structural time series models:* Standard models can be fitted with
  `StructTS()` in stats. Further packages are discussed in the
  `r view("TimeSeries")` task view.
- *Filtering and decomposition:* `decompose()` and `HoltWinters()` in stats.
  The basic function for computing filters (both rolling and autoregressive)
  is `filter()` in stats. Many extensions to these methods, in particular for
  forecasting and model selection, are provided in the `r pkg("forecast")`
  package.
- *Vector autoregression:* Simple models can be fitted by `ar()` in stats,
  more elaborate models are provided in package `r pkg("vars")` along with
  suitable diagnostics, visualizations etc. Structural smooth transition vector
  autoregressive models are in `r pkg("sstvars")` and panel vector
  autoregressions in `r pkg("panelvar")`. 
- *Unit root and cointegration tests:* `r pkg("urca", priority = "core")`,
  `r pkg("tseries", priority = "core")`, `r pkg("CADFtest")`. See also
  `r pkg("plm", priority = "core")` for panel unit root tests.
- *Miscellaneous:*
  - `r pkg("tsDyn")` - Threshold and smooth transition models.
  - `r pkg("midasr")` - *MIDAS regression* and other econometric methods for
    mixed frequency time series data analysis.
  - `r pkg("gets")` - GEneral-To-Specific (GETS) model selection for either
    ARX models with log-ARCH-X errors, or a log-ARCH-X model of the log
    variance.
  - `r pkg("bimets")` - Econometric modeling of time series data using
    flexible specifications of simultaneous equation models.
  - `r pkg("dlsem")` - Distributed-lag linear structural equation models.
  - `r pkg("lpirfs")` - Local projections impulse response functions.
  - `r pkg("apt")` - Asymmetric price transmission models.


### Data sets

- *Textbooks and journals:* Packages `r pkg("AER")`, `r pkg("Ecdat")`, and
`r pkg("wooldridge")` contain a comprehensive collections of data sets from
  various standard econometric textbooks (including Greene, Stock & Watson,
  Wooldridge, Baltagi, among others) as well as several data sets from the
  Journal of Applied Econometrics and the Journal of Business & Economic
  Statistics data archives. `r pkg("AER")` and `r pkg("wooldridge")`
  additionally provide extensive sets of examples reproducing analyses from
  the textbooks/papers, illustrating various econometric methods. In
  `r pkg("pder")` a wide collection of data sets for "Panel Data Econometrics
  with R" (Croissant & Millo 2018) is available. The
  `r github("ccolonescu/PoEdata")` package on GitHub provides the data sets from
  "Principles of Econometrics" (4th ed, by Hill, Griffiths, and Lim 2011).
- *Penn World Table:* `r pkg("pwt")` provides versions 5.6, 6.x, 7.x. Version
  8.x and 9.x data are available in `r pkg("pwt8")` and `r pkg("pwt9")`,
  respectively.
- *Time series and forecasting data:* The packages `r pkg("expsmooth")`,
`r pkg("fma")`, and `r pkg("Mcomp")` are data packages with time series data
  from the books "Forecasting with Exponential Smoothing: The State Space
  Approach" (Hyndman, Koehler, Ord, Snyder, 2008, Springer) and
  "Forecasting: Methods and Applications" (Makridakis, Wheelwright, Hyndman,
  3rd ed., 1998, Wiley) and the M-competitions, respectively.
- *Empirical Research in Economics:* Package `r pkg("erer")` contains
  functions and datasets for the book of "Empirical Research in Economics:
  Growing up with R" (Sun 2015).
- *Panel Study of Income Dynamics (PSID):* `r pkg("psidR")` can build panel
  data sets from the Panel Study of Income Dynamics (PSID).
- World Bank data and statistics: The `r pkg("wbstats")` package provides
  programmatic access to the World Bank API.


### Miscellaneous

- *Model tables:* A flexible implementation of side-by-side summary tables for
  a wide range of statistical models along with corresponding visualizations
  and data summary tables is provided in `r pkg("modelsummary")`. Other
  implementations as well as further utilities for integrating econometric
  and statistical results in scientific papers etc. are discussed in the
  `r view("ReproducibleResearch")` task view.
- *Matrix manipulations:* As a vector- and matrix-based language, base R
  ships with many powerful tools for doing matrix manipulations, which are
  complemented by the packages `r pkg("Matrix")` and `r pkg("SparseM")`.
- *Optimization and mathematical programming:* R and many of its contributed
  packages provide many specialized functions for solving particular
  optimization problems, e.g., in regression as discussed above. Further
  functionality for solving more general optimization problems, e.g.,
  likelihood maximization, is discussed in the the `r view("Optimization")`
  task view.
- *Bootstrap:* In addition to the recommended `r pkg("boot")` package, there
  are some other general bootstrapping techniques available in
  `r pkg("bootstrap")` or `r pkg("simpleboot")` as well some bootstrap techniques
  designed for time-series data, such as the maximum entropy bootstrap in
  `r pkg("meboot")` or the `tsbootstrap()` from `r pkg("tseries")`.
  The `r pkg("fwildclusterboot")` package provides a fast wild cluster
  bootstrap implementation for linear regression models, especially when
  the number of clusters is low.
- *Inequality:* For measuring inequality, concentration and poverty the
  package `r pkg("ineq")` provides some basic tools such as Lorenz curves,
  Pen's parade, the Gini coefficient, Herfindahl-Hirschman index and many more.
  `r pkg("wINEQ")` provides these and other inequality measures for weighted
  data along with bootstrapping methods.
- *Structural change:* R is particularly strong when dealing with structural
  changes and changepoints in parametric models, see `r pkg("strucchange")`
  and `r pkg("segmented")`.
- *Exchange rate regimes:* Methods for inference about exchange rate regimes,
  in particular in a structural change setting, are provided by `r pkg("fxregime")`.
- *Global value chains:* Tools and decompositions for global value chains are
  in `r pkg("gvc")` and `r pkg("decompr")`.
- *Regression discontinuity design:* A variety of methods are provided in the
  `r pkg("rdd")`, `r pkg("rdrobust")`, and `r pkg("rdlocrand")` packages. The
  `r pkg("rdpower")` package offers power calculations for regression
  discontinuity designs. And `r pkg("rdmulti")` implements analysis with
  multiple cutoffs or scores.
- *Gravity models:* Estimation of log-log and multiplicative gravity models
  is available in `r pkg("gravity")`.
- *z-Tree:* `r pkg("zTree")` can import data from the z-Tree software for
  developing and carrying out economic experiments.
- *Numerical standard errors:* `r pkg("nse")` implements various numerical
  standard errors for time series data, especially in simulation experiments
  with correlated outcome sequences.


### Links
- Articles: [Special Volume on "Econometrics in R" in JSS (2008)](http://www.jstatsoft.org/v27/)
- Book: [Applied Econometrics with R (Kleiber & Zeileis; 2008)](https://www.zeileis.org/teaching/AER/)
- Book: [Introduction to Econometrics with R (Hanck, Arnold, Gerber, & Schmelzer; 2021)](https://www.Econometrics-with-R.org/)
- Book: [Introduction to Econometrics with R (Oswald, Robin, & Viers; 2020)](https://scpoecon.github.io/ScPoEconometrics/)
- Book: [Causal Inference: The Mixtape (Cunningham; 2021)](https://mixtape.scunning.com/)
- Book: [Hands-On Intermediate Econometrics Using R (Vinod; 2008)](https://doi.org/10.1142/6895)
- Book: [Learning Microeconometrics with R (Adams; 2021)](https://sites.google.com/view/microeconometricswithr)
- Book: [Panel Data Econometrics with R (Croissant & Millo; 2018)](https://doi.org/10.1002/9781119504641)
- Book: [Principles of Econometrics with R (Colonescu; 2016)](https://bookdown.org/ccolonescu/RPoE4/)
- Book: [Spatial Econometrics (Kelejian & Piras; 2017)](https://doi.org/10.1016/C2016-0-04332-2)
- Book: [Statistical Inference via Data Science (Ismay & Kim; 2022)](https://moderndive.com/)
- Book: [The Effect (Huntington-Klein; 2022)](https://theeffectbook.net/)
- Book: [Using R for Introductory Econometrics (Heiss; 2019)](http://www.urfie.net/)
- Book: [Quantitative Economics with R (Dayal; 2020)](https://doi.org/10.1007/978-981-15-2035-8)
- Course: [Applied Empirical Methods (Goldsmith-Pinkham; 2021)](https://github.com/paulgp/applied-methods-phd)
- Course: [Data Science for Economists (McDermott; 2021)](https://github.com/uo-ec607/lectures)
- Course: [Econometrics In-Class Labs (Ransom; 2021)](https://tyleransom.github.io/econometricslabs.html)
- Course: [Introduction to Econometrics (Rubin; 2021)](https://github.com/edrubin/EC421W22)
- Course: [PhD Econometrics (Rubin; 2022)](https://github.com/edrubin/EC607S21)
- Course: [Program Evaluation for Public Service (Heiss, 2022)](https://evalsp22.classes.andrewheiss.com/)
- Course: [Statistical Rethinking (McElreath; 2022)](https://github.com/rmcelreath/stat_rethinking_2022)
- Website: [Stata2R](https://stata2r.github.io/)

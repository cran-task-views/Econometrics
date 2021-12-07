---
name: Econometrics
topic: Econometrics
maintainer: Achim Zeileis, Grant McDermott
email: Achim.Zeileis@R-project.org
version: 2021-12-07
source: https://github.com/cran-task-views/Econometrics
---

Base R ships with a lot of functionality useful for computational econometrics,
in particular in the stats package. This functionality is complemented by many
packages on CRAN, a brief overview is given below. There is also a considerable
overlap between the tools for econometrics in this view and those in the task
views on `r view("Finance")`, `r view("SocialSciences")`, and
`r view("TimeSeries")`.

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
  tests ( *z* instead of *t* tests, and Chi-squared instead of *F* tests) and
  plug-in of other covariance matrices are `coeftest()` and `waldtest()` in
  `r pkg("lmtest", priority = "core")`. Tests of more general linear hypotheses
  are implemented in `linearHypothesis()` and for nonlinear hypotheses in
  `deltaMethod()` in `r pkg("car", priority = "core")`.
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


### Microeconometrics

- *Generalized linear models (GLMs):* Many standard microeconometric models
  belong to the family of generalized linear models and can be fitted by
  `glm()` from package stats. This includes in particular logit and probit
  models for modeling choice data and Poisson models for count data. Effects
  for typical values of regressors in these models can be obtained and
  visualized using `r pkg("effects")`. Marginal effects tables for certain
  GLMs can be obtained using the `r pkg("margins")` and `r pkg("mfx")`
  packages. Interactive visualizations of both effects and marginal effects
  are possible in `r pkg("LinRegInteractive")`.
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
  package. Generalized additive models (GAMs) for multinomial responses can be
  fitted with the `r pkg("VGAM")` package. A Bayesian approach to multinomial
  probit models is provided by `r pkg("MNP")`. Various Bayesian multinomial
  models (including logit and probit) are available in `r pkg("bayesm")`.
  Furthermore, the package `r pkg("RSGHB")` fits various hierarchical Bayesian
  specifications based on direct specification of the likelihood function.
- *Ordered responses:* Proportional-odds regression for ordered responses is
  implemented in `polr()` from package `r pkg("MASS")`. The package
  `r pkg("ordinal")` provides cumulative link models for ordered data which
  encompasses proportional odds models but also includes more general
  specifications. Bayesian ordered probit models are provided by
  `r pkg("bayesm")`.
- *Censored responses:* Basic censored regression models (e.g., tobit models)
  can be fitted by `survreg()` in `r pkg("survival")`, a convenience interface
  `tobit()` is in package `r pkg("AER", priority = "core")`. Further censored
  regression models, including models for panel data, are provided in
  `r pkg("censReg")`. Censored regression models with conditional
  heteroscedasticity are in `r pkg("crch")`. Furthermore, hurdle models for
  left-censored data at zero can be estimated with `r pkg("mhurdle")`. Models
  for sample selection are available in `r pkg("sampleSelection")` and
  `r pkg("ssmrob")` using classical and robust inference, respectively. Package
  `r pkg("matchingMarkets")` corrects for selection bias when the sample is
  the result of a stable matching process (e.g., a group formation or college
  admissions problem).
- *Truncated responses:* `r pkg("crch")` for truncated (and potentially
  heteroscedastic) Gaussian, logistic, and t responses. Homoscedastic Gaussian
  responses are also available in `r pkg("truncreg")`.
- *Fraction and proportion responses:* Fractional response models are in
  `r pkg("frm")`. Beta regression for responses in (0, 1) is in
  `r pkg("betareg")` and `r pkg("gamlss")`.
- *Duration responses:* Many classical duration models can be fitted with
  `r pkg("survival")`, e.g., Cox proportional hazard models with `coxph()` or
  Weibull models with `survreg()`. Many more refined models can be found in
  the `r view("Survival")` task view. The Heckman and Singer mixed
  proportional hazard competing risk model is available in `r pkg("durmod")`.
- *High-dimensional fixed effects:* Linear models with potentially
  high-dimensional fixed effects, also for multiple groups, can be fitted by
  `r pkg("lfe")`. The corresponding GLMs are covered in `r pkg("alpaca")`.
  Another implementation, based on C++ code covering both OLS and GLMs is in
  `r pkg("fixest")`.
- *Miscellaneous:* Further more refined tools for microeconometrics are
  provided in the `r pkg("micEcon")` family of packages: Analysis with
  Cobb-Douglas, translog, and quadratic functions is in `r pkg("micEcon")`;
  the constant elasticity of scale (CES) function is in `r pkg("micEconCES")`;
  the symmetric normalized quadratic profit (SNQP) function is in
  `r pkg("micEconSNQP")`. The almost ideal demand system (AIDS) is in
  `r pkg("micEconAids")`. Stochastic frontier analysis (SFA) is in
  `r pkg("frontier")` and certain special cases also in `r pkg("sfa")`.
  Semiparametric SFA in is available in `r pkg("semsfa")` and spatial SFA in
  `r pkg("spfrontier")` and `r pkg("ssfa")`. The package `r pkg("bayesm")`
  implements a Bayesian approach to microeconometrics and marketing. Inference
  for relative distributions is contained in package `r pkg("reldist")`.


### Instrumental variables

- *Basic instrumental variables (IV) regression:* Two-stage least squares
  (2SLS) is provided by `r pkg("ivreg")` (previously in `r pkg("AER")`). Other
  implementations are in `tsls()` in package `r pkg("sem")`, in
  `r pkg("ivpack")`, and `r pkg("lfe")` (with particular focus on multiple group
  fixed effects).
- *Binary responses:* An IV probit model via GLS estimation is available in
  `r pkg("ivprobit")`. The `r pkg("LARF")` package estimates local average
  response functions for binary treatments and binary instruments.
- *Panel data:* Certain basic IV models for panel data can also be estimated
  with standard 2SLS functions (see above). Dedicated IV panel data models are
  provided by `r pkg("ivfixed")` (fixed effects) and `r pkg("ivpanel")`
  (between and random effects).
- *Miscellaneous:* `r pkg("REndo")` fits linear models with endogenous
  regressor using various latent instrumental variable approaches.


### Panel data models

- *Panel standard errors:* A simple approach for panel data is to fit the
  pooling (or independence) model (e.g., via `lm()` or `glm()`) and only
  correct the standard errors. Different types of clustered, panel, and
  panel-corrected standard errors are available in `r pkg("sandwich")`
  (incorporating prior work from `r pkg("multiwayvcov")`),
  `r pkg("clusterSEs")`, `r pkg("pcse")`, `r pkg("clubSandwich")`, `r pkg("plm",
  priority = "core")`, and `r pkg("geepack")`, respectively. The latter two
  require estimation of the pooling/independence models via `plm()` and
  `geeglm()` from the respective packages (which also provide other types of
  models, see below).
- *Linear panel models:* `r pkg("plm")`, providing a wide range of within,
  between, and random-effect methods (among others) along with corrected
  standard errors, tests, etc. Another implementation of several of these
  models is in `r pkg("Paneldata")`. Various dynamic panel models are
  available in `r pkg("plm")`, with estimation based on moment conditions in
  `r pkg("pdynmc")`, and dynamic panel models with fixed effects in
  `r pkg("OrthoPanels")`. `r pkg("feisr")` provides fixed effects individual
  slope (FEIS) models. Panel vector autoregressions are implemented in
  `r pkg("panelvar")`.
- *Generalized estimation equations and GLMs:* GEE models for panel data (or
  longitudinal data in statistical jargon) are in `r pkg("geepack")`. The
  `r pkg("pglm")` package provides estimation of GLM-like models for panel data.
- *Mixed effects models:* Linear and nonlinear models for panel data (and
  more general multi-level data) are available in `r pkg("lme4")` and `r pkg("nlme")`.
- *Instrumental variables:* `r pkg("ivfixed")` and `r pkg("ivpanel")`, see
  also above.
- *Miscellaneous:* Autocorrelation and heteroscedasticity correction are
  available in `r pkg("wahc")`. Threshold regression
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
  are available in `r pkg("spatialprobit")`.
- *Bayesian model averaging (BMA):* A comprehensive toolbox for BMA is
  provided by `r pkg("BMS")` including flexible prior selection, sampling,
  etc. A different implementation is in `r pkg("BMA")` for linear models,
  generalizable linear models and survival models (Cox regression).
- *Linear structural equation models:* `r pkg("lavaan")` and `r pkg("sem")`.
  See also the `r view("Psychometrics")` task view for more details.
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
  package stats is R\'s standard class for regularly spaced time series
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
  suitable diagnostics, visualizations etc. Panel vector autoregressions are
  available in `r pkg("panelvar")`.
- *Unit root and cointegration tests:* `r pkg("urca", priority = "core")`,
  `r pkg("tseries", priority = "core")`, `r pkg("CADFtest")`. See also
  `r pkg("pco")` for panel cointegration tests.
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
- *Canadian monetary aggregates:* `r pkg("CDNmoney")`.
- *Penn World Table:* `r pkg("pwt")` provides versions 5.6, 6.x, 7.x. Version
  8.x and 9.x data are available in `r pkg("pwt8")` and `r pkg("pwt9")`,
  respectively.
- *Time series and forecasting data:* The packages `r pkg("expsmooth")`,
`r pkg("fma")`, and `r pkg("Mcomp")` are data packages with time series data
  from the books \'Forecasting with Exponential Smoothing: The State Space
  Approach\' (Hyndman, Koehler, Ord, Snyder, 2008, Springer) and
  \'Forecasting: Methods and Applications\' (Makridakis, Wheelwright, Hyndman,
  3rd ed., 1998, Wiley) and the M-competitions, respectively.
- *Empirical Research in Economics:* Package `r pkg("erer")` contains
  functions and datasets for the book of \'Empirical Research in Economics:
  Growing up with R\' (Sun, forthcoming).
- *Panel Study of Income Dynamics (PSID):* `r pkg("psidR")` can build panel
  data sets from the Panel Study of Income Dynamics (PSID).
- World Bank data and statistics: The `r pkg("wbstats")` package provides
  programmatic access to the World Bank API.


### Miscellaneous

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
- *Inequality:* For measuring inequality, concentration and poverty the
  package `r pkg("ineq")` provides some basic tools such as Lorenz curves,
  Pen\'s parade, the Gini coefficient and many more.
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
- Journal of Statistical Software: [Special Volume on \'Econometrics in R\' (2008)](http://www.jstatsoft.org/v27/)
- Book: [Applied Econometrics with R (Kleiber and Zeileis)](https://eeecon.uibk.ac.at/~zeileis/teaching/AER/)
- Book: [Using R for Introductory Econometrics (Heiss)](http://www.urfie.net/)
- Book: [Introduction to Econometrics with R (Hanck, Arnold, Gerber, Schmelzer)](https://www.Econometrics-with-R.org/)
- Book: [Hands-On Intermediate Econometrics Using R (Vinod)](https://doi.org/10.1142/6895)
- Book: [Panel Data Econometrics with R (Croissant & Millo)](https://doi.org/10.1002/9781119504641)
- Book: [Spatial Econometrics (Kelejian and Piras)](https://doi.org/10.1016/C2016-0-04332-2)
- Manual: [Principles of Econometrics with R (Colonescu)](https://bookdown.org/ccolonescu/RPoE4/)
- Manual: [Introduction to Econometrics with R (Oswald, Robin, Viers)](https://scpoecon.github.io/ScPoEconometrics/)
- Manual: [Econometrics In-Class Labs (Ransom)](https://tyleransom.github.io/econometricslabs.html)
- Manual: [Data Science for Economists (McDermott)](https://github.com/uo-ec607/lectures)
- [A Brief Guide to R for Beginners in Econometrics](https://mondo.su.se/access/content/user/ma@su.se/Public/)
- [R for Economists](http://www.mayin.org/ajayshah/KB/R/R_for_economists.html)

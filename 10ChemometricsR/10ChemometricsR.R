#' ---
#' pagetitle: "R for Chemistry Data Analysis and Chemometrics. Chemometrics with R"
#' output:
#'   revealjs::revealjs_presentation:
#'     theme: simple
#'     slide_level: 2
#'     highlight: pygments
#'     center: false
#'     self_contained: true
#'     css: "../css/styles.css"
#'     reveal_options:
#'       slideNumber: true
#'       previewLinks: false
#'       transition: 0
#'       background_transition: 0
#' ---
#' 
## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, dev="svg")

#' ##
#' 
#' <img src="../images/IQSlogo.png" style="border-style:none;box-shadow:none;
#' position:absolute;margin:0;top:20px;left:20px;max-width:200px;height:auto;" />
#' 
#' <div style="font-size:1.5em;font-weight:700;margin-top:200px;">R for Chemistry Data Analysis and Chemometrics</div>
#' <div style="font-size:1.4em;font-weight:500;color:#333333;">Chemometrics with R</div>
#' <div style="font-size:1.2em;margin-top:40px;color:#333333;">Jordi Cuadros</div>
#' <div style="margin-top:80px;color:#333333;">January 2022</div>
#' 
#' # Introduction and References
#' 
#' ## Chemometrics
#' 
#' 
#' 
#' ## R packages for Chemistry and Chemometrics
#' 
#' - `chemometrics`, https://cran.r-project.org/web/packages/chemometrics/index.html
#' - `ChemometricsWithR`, https://github.com/rwehrens/ChemometricsWithR/
#' - `ChemoSpec`, https://cloud.r-project.org/web/packages/ChemoSpec/index.html
#' - `mdatools`, https://cran.r-project.org/web/packages/mdatools/index.html
#' 
#' Updated information can be found following the CRAN Task View: Chemometrics and Computational Physics, (https://cran.r-project.org/web/views/ChemPhys.html).
#' 
#' 
#' ## Chemometrics with R -- References
#' 
#' - Harvey, D. (2021). Chemometrics Using R (Harvey). Libretexts. https://chem.libretexts.org/Bookshelves/Analytical_Chemistry/Chemometrics_Using_R_(Harvey)
#' - Kucheryavskiy, s. (2021). mda.tools: make chemometrics easy. https://mdatools.com/
#' - Varmuza, K., & Filzmoser, P. (2016). Introduction to multivariate statistical analysis in chemometrics. CRC press.
#' - Wehrens, R. (2011). Chemometrics with R. Heidelberg, Germany: Springer.
#' 
#' 
#' # Chemical Data Preprocessing
#' 
#' ## Noise
#' 
#' ## Baseline
#' 
#' ## Peak Alignment
#' 
#' ## Scalling
#' 
#' ## Missing Data
#' 
#' ## Outliers
#' 
#' 
#' # Multivariant Exploratory Analysis
#' 
#' ## Correlations
#' 
#' ## Factor Analysis
#' 
#' ## Self-Organizing Maps
#' 
#' ## Clustering
#' 
#' 
#' # Classification
#' 
#' ## Discriminant Analysis
#' 
#' ## kNN
#' 
#' ## Classification Trees & Forests
#' 
#' ## SVN
#' 
#' ## Neural Networks 
#' 
#' 
#' # Modelling
#' 
#' ## Errors & Regression
#' 
#' ## Multiple Regression
#' 
#' ## Step-wise Regression
#' 
#' ## Regression Trees
#' 
#' ## PCR/PLS
#' 
#' ## Ridge-Regression
#' 
#' ## Non-linear Regression
#' 
#' 
#' 
#' # Multiple Curve Resolution
#' 

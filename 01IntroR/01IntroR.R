#' ---
#' pagetitle: "R for Chemistry Data Analysis and Chemometrics. Introduction to R"
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
## ----setup, include=FALSE---------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, dev="svg")

#' 
#' ## 
#' 
#' <img src="../images/IQSlogo.png" style="border-style:none;box-shadow:none;
#' position:absolute;margin:0;top:20px;left:20px;max-width:200px;height:auto;" />
#' 
#' <div style="font-size:1.5em;font-weight:700;margin-top:200px;">R for Chemistry Data Analysis and Chemometrics</div>
#' <div style="font-size:1.4em;font-weight:500;color:#333333;">Introduction to R</div>
#' <div style="font-size:1.2em;margin-top:40px;color:#333333;">Jordi Cuadros, Vanessa Serrano</div>
#' <div style="margin-top:80px;color:#333333;">January 2022</div>
#' 
#' 
#' # R Basics
#' 
#' ## R
#' 
#' R is a software environment for statistical computing and graphics. It is open source, free to use, and has a really large and active community. There more than 10000 packages with additional functions, tools and data sets. It is among the most used tools in data analysis and data science.
#' 
#' <https://www.r-project.org/about.html><br/>
#' <https://www.rdocumentation.org/><br/>
#' <https://stackoverflow.com/questions/tagged/r>
#' 
#' 
#' ## Installing R & RStudio
#' 
#' -   **R**
#' 
#'     R is available for multiple platforms at <https://cran.r-project.org/>.
#' 
#' -   **RStudio**
#' 
#'     RStudio is an Integrated Development Environment (IDE) for R.
#' 
#'     It can be downloaded and installed from <https://www.rstudio.com/>. RStudio Desktop (Open Source License) version is available at no cost.
#' 
#' ## RStudio Interface
#' 
#' ![](../images/RStudio.png)
#' 
#' ## Calculations in R
#' 
#' Calculations can be done in the RStudio console. Enter the expression, press `Intro` to execute.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
2 + 3 
8 ^ 10
log(100) # Natural logarithm

#' 
#' [\# introduces a single-line comment]{.bibref}
#' 
#' ## Using Scripts
#' 
#' A script is a text file that contains a set or R instructions. This allows saving and restoring workflows as well as running them in unmanned mode. In R, they are commonly saved with the `.R` extension.
#' 
#' **RStudio** is a highly recommended editor for R. It allows creating, executing and debugging these scripts. Scripts are edited on the top-left pane.
#' 
#' Current selection (or instruction if nothing is selected) is executed by pressing `Ctrl + Return` or hitting the `Run` command button on the interface.
#' 
#' ## Help
#' 
#' -   To look up a function
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## ? "mean"
## help("mean")

#' 
#' -   To search for a text in the help files
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## ?? "anova"
## help.search("anova")

#' 
#' In **RStudio**, `F1` can be used to search for the selected text in help. This works only in the pane for editing scripts.
#' 
#' # Basic Data Types
#' 
#' ## Basic Data Types
#' 
#' Basic or atomic data type are
#' 
#' -   `numeric`: 15-digit decimal number
#' -   `integer`: 32-bit integer (to \~2·10^9^)
#' -   `logical`: `TRUE` or `FALSE`
#' -   `character`: character string of undefined length
#' 
#' Although less common in science, there also `complex` and `raw` data types.
#' 
#' ## Basic Data Types -- `numeric`
#' 
#' To store numerical quantities which are continuous in nature.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
a <- 2
a
class(a)
b <- 13.6788956789
b
print(b, digits = 10)

#' 
#' ------------------------------------------------------------------------
#' 
#' ### Common operators for `numeric` data
#' 
#' |                     |                   |
#' |:--------------------|:-----------------:|
#' | Addition            |        `+`        |
#' | Subtraction         |        `-`        |
#' | Product             |        `*`        |
#' | Division            |        `/`        |
#' | Power               |    `^` o `**`     |
#' | Modulus (remainder) |       `%%`        |
#' | Integer division    |       `%/%`       |
#' | Comparisons         | `== > < >= <= !=` |
#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
a+b
a-b
a*b
a^b
a/b


#' 
#' ## Basic Data Types -- `integer`
#' 
#' To store integer numbers, such as counts and indices.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
n <- as.integer(340000)
class(n)
n <- 2L
class(n)

#' 
#' ## Basic Data Types -- `character`
#' 
#' To store text (character strings) of any length.
#' 
#' Character-string literals are quoted (using double or single quotes).
#' 
## ---- echo = TRUE-----------------------------------------------------------------
a <- "aaa"
b <- "bbb"

paste(a, b, "hola", sep = ", ") 

#' 
#' ------------------------------------------------------------------------
#' 
#' ### Common operators for `character` data
#' 
#' |             |                   |
#' |:------------|:-----------------:|
#' | Comparisons | `== > < >= <= !=` |
#' 
#' ## Basic Data Types -- `logical`
#' 
#' To store values that are either `r T` or `r F`. This is the result of a comparison.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
3 == 2
b <- 3 != 2 
b

#' 
#' ------------------------------------------------------------------------
#' 
#' ### Common operators for `logical` data
#' 
#' |                      |         |
#' |:---------------------|:-------:|
#' | Logical intersection |   `&`   |
#' | Logical union        |   `|`   |
#' | Logical negation     |   `!`   |
#' | Comparisons          | `== !=` |
#' 
#' For multiple logical values, `any` and `all` can be used for logical union and intersection respectively.
#' 
#' ------------------------------------------------------------------------
#' 
#' For example...
#' 
## ---- echo = TRUE-----------------------------------------------------------------
a <- TRUE
b <- F
a & b # Operator AND
a | b # Operator OR
!b  # Operator NOT
all(a, !b, T)
any(a, b, F)

#' 
#' ## Conversion between Basic Data Types
#' 
#' Data can be converted to a different data type by using a type conversion function. These are named `as.` followed by the destination data type in R.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
as.character(2)
as.numeric(TRUE)

#' 
#' ------------------------------------------------------------------------
#' 
#' Be aware that not all conversions are possible.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
as.numeric("two")
as.logical("2+2==4")

#' 
#' To check for a specific data type, functions `is.` followed by the data type are also available.
#' 
#' ## Constants
#' 
#' Some constants are readily available in R.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
pi
LETTERS  # This is vector. We'll come back to that.

#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## Inf
## NaN
## NA
## NULL

#' 
#' ------------------------------------------------------------------------
#' 
#' Functions to check if a value matches one of these constants are also available.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
is.na(3/0)
is.null(NULL)
is.infinite(-3e999)

#' 
#' # Data Structures
#' 
#' ## Data Structures
#' 
#' -   Vector
#' -   Factor
#' -   Data frame
#' -   Other data structures
#'     -   Ordered Factor
#'     -   List
#'     -   Matrix
#' -   Objects: S3, S4, R6... Check <https://adv-r.hadley.nz/oo.html> for more information.
#' 
#' 
#' ## Data Structures -- vector
#' 
## ---- echo = TRUE-----------------------------------------------------------------
a <- c(2, 3, 4)
str(a)
a[2]

#' 
#' ------------------------------------------------------------------------
#' 
#' Most operations are applied to vectors element-wise.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
a
a + a 
a > 2.5

#' 
#' ------------------------------------------------------------------------
#' 
#' Be cautious with data recycling...
#' 
## ---- echo = TRUE-----------------------------------------------------------------
a <- c(2, 3, 4)
b <- c(10, 20)
a * b

#' 
#' ------------------------------------------------------------------------
#' 
#' Most R functions take vectors and return vectors.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
a <- c(2, 3, 4)
abs(sin(a))
exp(a)

#' 
#' ------------------------------------------------------------------------
#' 
#' Other returns aggregated values.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
length(a)
sum(a)
mean(a)

#' 
#' Also `sd`, `max`, `min`...
#' 
#' ------------------------------------------------------------------------
#' 
#' The presence of an element in vector can be checked with `%in%`.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
a <- c(2, 3, 4)
3 %in% a
c(2,5,3,4) %in% a

#' 
#' ------------------------------------------------------------------------
#' 
#' ### Sequence vectors
#' 
#' Sequence vectors are created with `:` or `seq`.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
1:10
seq(1, 10, by = 2)
seq(1, 3, length.out = 5) 

#' 
#' ------------------------------------------------------------------------
#' 
#' ### Sorting vectors
#' 
#' Vectors are ordered with `sort` or `order`.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
a <- c(8, 2, 5, 3)
sort(a)
a[order(a)]
a[order(a,decreasing = T)]

#' 
#' ------------------------------------------------------------------------
#' 
#' ### Selecting elements
#' 
## ---- echo = TRUE-----------------------------------------------------------------
a <- 10:40
a[3]
a[4:6]
a[7] <- 0
b <- a[a!=0]
b
length(b)

#' 
#' 
#' ## YOUR TURN {data-background=#eeffcc}
#' 
#' 1.    Create a numerical vector that includes all multiples of seven up to 1000.
#' 
#' 2.    Exclude the numbers that have a last digit equal to 3. Here is a hint (`%%` is the modulus/remainder operation): try `1004 %% 10`
#' 
#' 3.    How many numbers are left in the vector?
#' 
#' 4.    How many of these numbers have a 5 in their representation?
#' 
#' 
#' ## Data Structures -- factor
#' 
#' A factor is an indexed character vector. It is usually created from a character vector.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
a <- c("hola", "adeu","hola", "adeu", "adeu", "bye")
b <- as.factor(a)
as.character(b)   # returns a character vector
as.numeric(b)     # returns an integer vector of level indices
str(b)

#' 
#' ------------------------------------------------------------------------
#' 
#' The `factor` function allows manually coding or re-coding the factor.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
a <- factor(c(3, 1, 3, 1, 1, 2), labels = c("adeu", "bye", "hola"))
a
levels(a)

#' 
#' ## Data Structures -- data frame
#' 
#' A data frame is a rectangular structure of data, organized such as each column is a vector. Different columns may be of different data types.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
dfA <- data.frame(int = 1:10,
                  let = sample(letters, 10, replace = TRUE), 
                  ran = rnorm(10))
dfA

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
dim(dfA) # Dimensions
nrow(dfA) # Row count
ncol(dfA) # Column count

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
str(dfA)
head(dfA, 3) # First rows, 6 by default
tail(dfA, 2) # Last rows

#' 
#' ------------------------------------------------------------------------
#' 
#' To access the data in a data frame, we use indices for row and column (starting at 1). Variables can also be selected by column name by using `$`.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
dfA[2,3]
dfA[,1]
dfA$let

#' 
#' Additional options exist to select rows and columns. We'll go into that later on.
#' 
#' 
#' ## YOUR TURN {data-background=#eeffcc}
#' 
#' 1.    Create a data frame that includes five different properties for eight *n*-alkanes. Data can be obtained from <https://chem.libretexts.org/Bookshelves/Organic_Chemistry/Book%3A_Basic_Principles_of_Organic_Chemistry_(Roberts_and_Caserio)/04%3A_Alkanes/4.02%3A_Physical_Properties_of_Alkanes_and_The_Concept_of_Homology> or <https://doi.org/10.1007/s40747-020-00262-0> for example. Include, at least, IUPAC name, formula number of carbons and boiling point. 
#' 
#' 
#' ## Other Data Structures -- list
#' 
#' Lists are one-index structures that can hold data of different types.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
a <- list(2, "2", FALSE)
b <- list(3, "hola", c(2, 3, 4))

#' 
#' ::: {style="column-count:2;"}
## ---- echo = TRUE-----------------------------------------------------------------
a

#' 
#'  
#' 
#' <p style="display:block;break-after:column;">
#' 
#' </p>
#' 
## ---- echo = TRUE-----------------------------------------------------------------
b

#' :::
#' 
#' ------------------------------------------------------------------------
#' 
#' ::: {style="column-count:2;"}
## ---- echo = TRUE-----------------------------------------------------------------
length(a)
a[[3]]
b[[3]][1]

#' 
#'  
#' 
#' <p style="display:block;break-after:column;">
#' 
#' </p>
#' 
## ---- echo = TRUE-----------------------------------------------------------------
str(b)

#' :::
#' 
#' ## Other Data Structures -- ordered factor
#' 
## ---- echo = TRUE-----------------------------------------------------------------
grades <- c("Pass", "Fail", "Good", "Fail",
           "Good", "Excellent", "Pass")
grades <- factor(grades,
        levels = c("Fail", "Pass",
                   "Good", "Excellent"),
        ordered = TRUE)
str(grades)
levels(grades)

#' 
#' ## Other Data Structures -- matrix
#' 
#' A matrix is a rectangular data structure where all data share the same data type.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
a <- matrix(c(2, 4, -1, 5), ncol = 2)
a
a[2,2]

#' 
#' ------------------------------------------------------------------------
#' 
#' ::: {style="column-count:2;"}
#' 
## ---- echo = TRUE-----------------------------------------------------------------
a * a # Element-wise product
a %*% a # Matrix product
t(a) # Transposition
det(a) # Determinant

#' 
#'  
#' 
#' <p style="display:block;break-after:column;">
#' 
#' </p>
#' 
## ---- echo = TRUE-----------------------------------------------------------------
solve(a) # Inverse
a %*% solve(a)

#' :::
#' 
#' 
#' # Basic Graphics in R
#' 
#' ## ¿Why Using Graphics?
#' 
#' Main uses of graphics and visualization in data analysis are
#' 
#' -   data exploration and interpretation,
#' -   non-evident pattern discovery, and
#' -   communication of analysis results.
#' 
#' ## Graphics in R
#' 
#' R includes many paradigms for graphics development. These include
#' 
#' -   R `base` plots,
#' -   Grammar Of Graphics (*GoG*) based representations, with the `ggplot2` package,
#' -   Lattice charts, using the `lattice` package,
#' -   A formula-based grammar of graphics, `ggformula`,
#' -   `ggvis`, an interactive grammar of graphics framework, among many others.
#' 
#' We will here discuss the simplest charts, those from `base` R.
#' 
#' 
#' ## Graphics in `base` R
#' 
#' Graphics in `base` R are usually constructed from vectors by using specific functions depending on the desired chart.
#' 
#' To exemplify some of these functions we will use the `airquality` data set, one of many data sets included in `base`R.
#' 
#' Access help for details on the data set..
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## help("airquality")

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
str(airquality)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
head(airquality, 10)

#' 
#' ------------------------------------------------------------------------
#' 
#' We will show how to make
#' 
#' -   scatterplots,
#' -   histograms,
#' -   bar plots, and
#' -   boxplots
#' 
#' 
#' ## Graphics in `base` R -- scatterplot
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## plot(airquality$Solar.R,airquality$Ozone)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE, eval = TRUE---------------------------------------------------
plot(airquality$Solar.R,airquality$Ozone)

#' 
#' ------------------------------------------------------------------------
#' 
#' ## Graphics in `base` R -- histogram
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## hist(airquality$Ozone)

#' 
#' With a density curve...
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## hist(airquality$Ozone, breaks = 23, freq = FALSE)
## lines(density(airquality$Ozone, na.rm=TRUE))

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE, eval = TRUE---------------------------------------------------
hist(airquality$Ozone)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE, eval = TRUE---------------------------------------------------
hist(airquality$Ozone, breaks = 23, freq = FALSE)
lines(density(airquality$Ozone, na.rm=TRUE))

#' 
#' ------------------------------------------------------------------------
#' 
#' ## Graphics in `base` R -- bar plot
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## barplot(table(cut(airquality$Wind, breaks=seq(0,22,by=2))))

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE, eval = TRUE---------------------------------------------------
barplot(table(cut(airquality$Wind, breaks=seq(0,22,by=2))))

#' 
#' ## Graphics in `base` R -- boxplot
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## boxplot(airquality$Ozone)
## points(mean(airquality$Ozone, na.rm=TRUE),pch="+")

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE, eval = TRUE---------------------------------------------------
boxplot(airquality$Ozone)
points(mean(airquality$Ozone, na.rm=TRUE),pch="+")

#' 
#' 
#' # Basic Stats with R
#' 
#' ## Statistics Goals
#' 
#' **Statistics** is the mathematics sub-discipline which covers the collection, analysis and interpretation of data. It has two main goals:
#' 
#' 1.  To describe a set of data (**descriptive statistics**), y
#' 2.  To extract conclusions about the population from the available data (**inferential statisitcs**)
#' 
#' In this statistics review, we will again use the `airquality` data set. Details on the data set are avalable at...
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## help("airquality")

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
str(airquality)

#' 
#' 
#' ## Descriptive statistics for one variable
#' 
#' Often, we need to start by exploring and describing each variable data. Let's start here. 
#' 
#' Common aspects to consider when analyzing a single variable are 
#' 
#' -   distribution: frequencies and quantiles,
#' -   central tendency (also referred as localization o position),
#' -   spread, and
#' -   analysis of outliers.
#' 
#' 
#' ## Descriptive statistics for one variable -- distribution
#' 
#' ### Summary 
## ---- echo = TRUE-----------------------------------------------------------------
summary(airquality)

#' 
#' ------------------------------------------------------------------------
#' 
#' ### Sample size or number of data points
## ---- echo = TRUE-----------------------------------------------------------------
length(airquality$Ozone)
n <- sum(!is.na(airquality$Ozone))
n

#' 
#' ------------------------------------------------------------------------
#' 
#' ### Absolute and relative frequencies
#' 
#' For a quantitative variable (`numeric`)...
#' 
## ---- echo = TRUE-----------------------------------------------------------------
frec_abs <- table(cut(airquality$Solar.R,
      breaks=c(0,50,100,150,200,
               250,300,350)))
frec_abs

frec_rel <- frec_abs / sum(frec_abs)
frec_rel

#' 
#' ------------------------------------------------------------------------
#' 
#' For a qualitative variable (usually `factor`)...
#' 
## ---- echo = TRUE-----------------------------------------------------------------
frec_abs <- table(as.factor(airquality$Month))
frec_abs

frec_rel <- frec_abs / sum(frec_abs)
frec_rel

#' 
#' ------------------------------------------------------------------------
#' 
#' ### Quantiles
#' 
## ---- echo = TRUE-----------------------------------------------------------------
quantile(airquality$Wind,.1)
quantile(airquality$Wind,0:5*.2)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
# Excel PERCENTILE and R default are type 7
c(min = min(airquality$Wind), quantile(airquality$Wind,0:5*.2),
  max = max(airquality$Wind))    
# Excel PERCENTILE.EXC, SPSS and Minitab default are type 6
c(min = min(airquality$Wind), quantile(airquality$Wind,0:5*.2, type=6),
  max = max(airquality$Wind))


#' 
#' 
#' ## Descriptive statistics for one variable -- central tendency
#' 
## ---- echo = TRUE-----------------------------------------------------------------
mean(airquality$Ozone)
mean(airquality$Ozone, na.rm = TRUE)
median(airquality$Ozone, na.rm = TRUE)

#' 
#' ## Descriptive statistics for one variable -- spread
#' 
## ---- echo = TRUE-----------------------------------------------------------------
var(airquality$Ozone, na.rm = TRUE) # denominator is n-1
sd(airquality$Ozone, na.rm = TRUE)
range(airquality$Ozone, na.rm = TRUE)
diff(range(airquality$Ozone, na.rm = TRUE))

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
IQR(airquality$Ozone, na.rm = TRUE)
mad(airquality$Ozone, na.rm = TRUE) # mean absolute deviation

#' 
#' ## Descriptive statistics for one variable -- outliers
#' 
#' There different methods to analyze outliers in R. A nice description can be found at <https://statsandr.com/blog/outliers-detection-in-r/>
#' 
## ---- echo = TRUE-----------------------------------------------------------------
# Hampel filter  
lmin <- median(airquality$Wind,.25) - 3 * mad(airquality$Wind)
lmax <- median(airquality$Wind,.75) + 3 * mad(airquality$Wind) 
c(as.numeric(lmin),as.numeric(lmax))
airquality$Wind[airquality$Wind > lmax | airquality$Wind < lmin]

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
# IQR method
lmin <- quantile(airquality$Wind,.25) - 1.5 * IQR(airquality$Wind)
lmax <- quantile(airquality$Wind,.75) + 1.5 * IQR(airquality$Wind) 
c(as.numeric(lmin),as.numeric(lmax))
airquality$Wind[airquality$Wind > lmax | airquality$Wind < lmin]
boxplot.stats(airquality$Wind)$out

#' 
#' 
#' ## Descriptive statistics for one variable -- plots
#' 
#' Plots are also useful to describe a data set.
#' 
#' -   For quantitative variables: boxplot and histogram
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## boxplot(airquality$Wind)
## points(mean(airquality$Wind),pch=3)
## 
## hist(airquality$Wind)

#' 
#' -   For qualitative variables: bar plot
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## barplot(table(airquality$Month))

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
boxplot(airquality$Wind)
points(mean(airquality$Wind),pch=3)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
hist(airquality$Wind)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
barplot(table(airquality$Month))

#' 
#' 
#' ## Descriptive statistics for two variables
#' 
#' To study and describe the relation between two variables, common techniques include
#' 
#' -   contingency tables,
#' -   correlation coefficients, and
#' -   scatterplots.
#' 
#' 
#' ## Descriptive statistics for two variables -- cotingency table
#' 
#' The contingency table is a two-way frequency table.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
table(cut(airquality$Wind,breaks = seq(1,21,by=5)),
      cut(airquality$Temp,breaks = seq(50,100,by=10)))

#' 
#' ## Descriptive statistics for two variables -- correlation coeffcient
#' 
#' ### Pearson product-moment correlation coefficent
#' 
## ---- echo = TRUE-----------------------------------------------------------------
cor(airquality$Temp,airquality$Wind)

#' 
#' ### Spearman correlation coefficent
#' 
## ---- echo = TRUE-----------------------------------------------------------------
cor(airquality$Temp,airquality$Wind,method = "spearman")

#' 
#' ## Descriptive statistics for two variables -- scatterplot
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## plot(airquality$Temp,airquality$Wind)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
plot(airquality$Temp,airquality$Wind)

#' 
#' 
#' ## Probability Distributions
#' 
#' A probability distribution of a random number corresponds to the abstraction of the density distribution of a set of infinite number of data with the same origin.
#' 
#' A large number of experimental distributions have related theoretical models. The most common theoretical distributions are the normal distribution (for continuous variables), the uniform distribution (for continuous or discrete variables), and the binomial distribution (for the number of successes of a discrete event with a defined probability).
#' 
#' ------------------------------------------------------------------------
#' 
#' In **R**, all theoretical distributions share the same system of functions
#' 
#' -   `r<dist>`: to generate random numbers according to the distribution
#' -   `d<dist>`: to calculate the density for a value of the variable
#' -   `q<dist>`: to obtain the value of the variable (quantile) for an accumulated probability
#' -   `p<dist>`: to obtain the accumulated probability for a value of the variable
#' 
#' ------------------------------------------------------------------------
#' 
#' For example, for the normal distribution...
#' 
## ---- echo = TRUE-----------------------------------------------------------------
rnorm(10)
dnorm(0)
qnorm(.95)
pnorm(1.64)

#' 
#' ## Probability Distributions -- normal
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## df <- data.frame(x = rnorm(1000, mean = 3, sd = 1))
## dfT <-data.frame(x = seq(0,6,length.out=101),
##       y = dnorm(seq(0,6,length.out=101),mean=3,sd=1))
## hist(df$x,breaks=2*ceiling(max(df$x)-min(df$x)))
## lines(dfT$x,dfT$y*1000*0.5)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
df <- data.frame(x = rnorm(1000, mean = 3, sd = 1))
dfT <-data.frame(x = seq(0,6,length.out=101),
      y = dnorm(seq(0,6,length.out=101),mean=3,sd=1))
hist(df$x,breaks=2*ceiling(max(df$x)-min(df$x)))
lines(dfT$x,dfT$y*1000*0.5)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
pnorm(5,mean = 3,sd = 1)
qnorm(.98,mean = 3,sd = 1)
pnorm(1:5,mean = 0,sd = 1)
qnorm(c(0.95,0.975,.99,.995,.999),mean = 0,sd = 1)


#' 
#' ## Probability Distributions -- uniforme
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## df <- data.frame(x = runif(1000, min = 10, max = 20))
## dfT <-data.frame(x = seq(10,20,length.out=101),
##       y = dunif(seq(10,20,length.out=101), min=10, max=20))
## hist(df$x,breaks=2*ceiling(max(df$x)-min(df$x)))
## lines(dfT$x,dfT$y*1000*0.5)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
df <- data.frame(x = runif(1000, min = 10, max = 20))
dfT <-data.frame(x = seq(10,20,length.out=101),
      y = dunif(seq(10,20,length.out=101), min=10, max=20))
hist(df$x,breaks=2*ceiling(max(df$x)-min(df$x)))
lines(dfT$x,dfT$y*1000*0.5)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
punif(12, min=10, max=20)
qunif(.90, min=10, max=20)

#' 
#' ## Probability Distributions -- binomial
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## df <- data.frame(x = rbinom(100,5,prob=0.5))
## barplot(table(df$x))

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
df <- data.frame(x = rbinom(100,5,prob=0.5))
barplot(table(df$x))

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
pbinom(2,5,prob=0.5)
qbinom(.5,5,.5)

#' 
#' ## Probability Distributions -- other
#' 
#' **R** includes many other distributions that can be looked up in the help files.
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## ? "Distributions"

#' 
#' 
#' ## Inference
#' 
#' En términos estadísticos, la **inferencia** consiste en la obtención de información sobre la población --el conjunto de los valores existentes-- a partir de una muestra representativa de la misma.
#' 
#' Si los datos que tenemos no constituyen una muestra *representativa*, cualquier información que queramos extrapolar a la población estará sesgada, es decir, desplazada de su valor real.
#' 
#' ------------------------------------------------------------------------
#' 
#' Los resultados de la inferencia se suelen mostrar mediante una de las dos formas siguientes:
#' 
#' -   el valor de un parámetro para la población, asociado normalmente a un error o un intervalo de confianza,
#' -   la probabilidad (*p*) que unos datos, o cualquiera más extremo que ellos, se produzcan siendo cierta una determinada hipótesis, denominada hipótesis nula.
#' 
#' ------------------------------------------------------------------------
#' 
#' Presentaremos a continuación, las técnicas más comúnmente usadas para
#' 
#' -   analizar el ajuste a una **distribución**,
#' -   estudiar la **dispersión** de la población, y
#' -   estimar o comparar valores de **tendencia central** --localización o posición--
#' 
#' ## Inferencia -- distribución
#' 
#' Para contrastar la hipótesis que un conjunto de datos puede representarse mediante un determinada distribución teórica, las pruebas más comunes son la **prueba de bondad de ajuste de chi cuadrado** (`chisq.test`) y la **prueba de Kolmogorov-Smirnov** (`ks.test`). En ambos casos, los argumentos son las frecuencias absolutas de cada clase para la distribución experimental y las probabilidades esperadas de acuerdo con la distribución teórica o de referencia.
#' 
#' ------------------------------------------------------------------------
#' 
#' **Ejemplo**
#' 
#' Se ha tirado cien veces un dado de 10 caras, y se han obtenido cada una de las caras las veces que se muestran en la tabla siguiente.
#' 
#' |     |     |     |     |     |     |     |     |     |     |
#' |:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
#' |  1  |  2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  | 10  |
#' |  9  |  6  |  7  | 15  | 11  | 11  |  9  | 13  |  9  | 10  |
#' 
#' Si el dado fuese normal, todas las caras tendrían la misma probabilidad de aparecer. ¿Lo es?
#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
obs <- c(9, 6, 7, 15, 11, 11, 9, 13, 9, 10)  
exp <- rep(.1, 10)
 
chisq.test(x = obs, p = exp)

#' 
#' ------------------------------------------------------------------------
#' 
#' Gráficamente, y aunque sin la misma objetividad, la comparación de distribuciones también suele hacerse mediante el gráfico de cuantiles --en R, es la función `qqplot`.
#' 
#' Si ambos conjuntos de datos corresponden a la misma distribución, sus cuantiles deben quedar alineados junto a la diagonal del gráfico (especialmente en la región central del mismo).
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## exp <- rep(1:10, obs)
## teo <- rep(1:10, 10)
## qqplot(x = exp, y = quantile(teo, (1:length(exp))/length(exp)))
## qqline(y = exp, distribution = function(x)  quantile(teo,x),
##        probs= c(1/length(exp),1))

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
exp <- rep(1:10, obs)
teo <- rep(1:10, 10)

#' 
## ---- echo = FALSE----------------------------------------------------------------
qqplot(x = exp, y = quantile(teo, (1:length(exp))/length(exp)))
qqline(y = exp, distribution = function(x)  quantile(teo,x),
       probs= c(1/length(exp),1))

#' 
#' ------------------------------------------------------------------------
#' 
#' Para comprobar si un conjunto de datos puede proceder de una población aleatoria de distribución normal, existen pruebas estadísticas más específicas. Entre las muchas pruebas existentes, son de uso común las pruebas de Shapiro-Wilk y Anderson-Darling. Solo la primera está disponible R *base*.
#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
x1 <- rnorm(20, mean = 4, sd = 5)
x2 <- rbeta(20, shape1 = 5, shape2 = .5, ncp = 4)

shapiro.test(x1)
shapiro.test(x2)


#' 
#' ------------------------------------------------------------------------
#' 
#' Gráficamente se puede usar el gráfico de cuantiles comentado anteriormente.
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## qqnorm(x1)
## qqline(x1)
## 
## qqnorm(x2)
## qqline(x2)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
qqnorm(x1)
qqline(x1)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
qqnorm(x2)
qqline(x2)

#' 
#' ## Inferencia -- dispersión
#' 
#' Las inferencias más habituales respecto a la dispersión de una población constituyen los siguientes casos
#' 
#' -   determinar el intervalo de confianza de la varianza o la desviación típica de la distribución,
#' -   comparar la dispersión de dos conjuntos de datos,
#' -   comparar la dispersión de tres o más conjuntos de datos --o comprobar la homocedasticidad de distintos conjuntos de datos--.
#' 
#' ------------------------------------------------------------------------
#' 
#' Cuando los datos estan normalmente distribuidos, el **intervalo de confianza de la varianza** se calcula a partir de la distribución de chi cuadrado.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
x <- rnorm(50, mean = 0, sd = 2)
df <- length(x) - 1
lower <- var(x) * df / qchisq(1 - 0.05/2, df)
upper <- var(x) * df / qchisq(0.05/2, df)
c(lower = lower, var = var(x), upper = upper)

#' 
#' ------------------------------------------------------------------------
#' 
#' Para obtener el **intervalo de confianza para la desviación estándar**, cuando los datos estan normalmente distribuidos, se calcula la desviación estándar a partir de las varainzas calculadas más arriba.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
c(lower = sqrt(lower), sd = sd(x), upper = sqrt(upper))

#' 
#' ------------------------------------------------------------------------
#' 
#' Para la **comparación de la dispersión de las poblaciones para dos conjuntos de datos**, cuando los datos estan normalmente distribuidos, se lleva a cabo con un prueba de F --función `var.test` en R--.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
x <- rnorm(50, mean = 0, sd = 2)
y <- rnorm(30, mean = 1, sd = 1)
var.test(x, y)

#' 
#' ------------------------------------------------------------------------
#' 
#' Si los datos no estan normalmente distribuidos, la comparación de dispersiones puede hacerse mediante el test de Ansari.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
x <- rlnorm(50, meanlog = 2, sdlog = 1)
y <- rlnorm(30, meanlog = 2, sdlog = .2)
ansari.test(x, y)

#' 
#' ------------------------------------------------------------------------
#' 
#' Para comparar las dispersiones de más de un conjunto de datos para los cuáles se pueda asumir normalidad, se usa la prueba de Bartlett --`bartlett.test` en R--.
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## x1 <- round(rnorm(20, mean = 1, sd = 2),1)
## x2 <- round(rnorm(20, mean = 3, sd = 2),1)
## x3 <- round(rnorm(20, mean = 5, sd = 2),1)
## list(x1,x2,x3)
## bartlett.test(list(x1,x2,x3))

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
x1 <- round(rnorm(20, mean = 1, sd = 2),1)
x2 <- round(rnorm(20, mean = 3, sd = 2),1)
x3 <- round(rnorm(20, mean = 5, sd = 2),1)
list(x1,x2,x3)
bartlett.test(list(x1,x2,x3))

#' 
#' ------------------------------------------------------------------------
#' 
#' Si los datos no cumplen con la hipótesis de normalidad, en lugar de la prueba de Bartlett se usa la prueba de Levene o la de Brown--Forsythe --`leveneTest` en el paquete `car` de R--.
#' 
#' ------------------------------------------------------------------------
#' 
#' Gráficamente y aunque de forma menos objetiva, la comparación de dispersiones puede hacerse usando gráficos de caja, bien directamente, bien centrando los datos.
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## x1 <- round(rnorm(20, mean = 1, sd = 2),1)
## x2 <- round(rnorm(20, mean = 3, sd = 2),1)
## x3 <- round(rnorm(20, mean = 5, sd = 2),1)
## 
## xs <- data.frame(
##   group = factor(c(rep(1, length(x1)),rep(2, length(x2)),
##             rep(3, length(x3)))),
##   y = c(x1,x2,x3),
##   centered = c(x1 - median(x1),x2 - median(x2),x3 - median(x3)))
## 
## boxplot(y~group,data=xs)
## boxplot(centered~group,data=xs)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
x1 <- round(rnorm(20, mean = 1, sd = 2),1)
x2 <- round(rnorm(20, mean = 3, sd = 2),1)
x3 <- round(rnorm(20, mean = 5, sd = 2),1)

xs <- data.frame(
  group = factor(c(rep(1, length(x1)),rep(2, length(x2)),
            rep(3, length(x3)))),
  y = c(x1,x2,x3),
  centered = c(x1 - median(x1),x2 - median(x2),x3 - median(x3)))

boxplot(y~group,data=xs)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
boxplot(centered~group,data=xs)

#' 
#' ## Inferencia -- localización
#' 
#' Las principales inferencias respecto a la localización --o tendencia central-- de los datos corresponden a
#' 
#' -   establecer un intervalo de confianza para la tendencia central de la distribución,
#' -   comparar la tendencia central de un conjunto de datos con un valor preestablecido,
#' -   comparar la tendencia central de dos conjuntos de datos,
#' -   comparar la tendencia central de tres o más conjuntos de datos.
#' 
#' ------------------------------------------------------------------------
#' 
#' Para establecer el **intervalo de confianza de la media**, cuando los datos están normalmente distribuidos, se usa la distribución t. El intervalo puede calcularse a partir de la función `qt` o visualizarlo como resultado de la función `t.test`.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
x <- rnorm(10, mean = 2, sd = .5)
x

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
t.test(x)

#' 
#' ------------------------------------------------------------------------
#' 
#' Si los datos no están normalmente distribuidos, la función `wilcox.test` ofrece entre sus resultados el **intervalo de confianza para la mediana**.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
x <- rnorm(10, mean = 3, sd = .5) ^ 3
x

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
wilcox.test(x, conf.int = TRUE)

#' 
#' ------------------------------------------------------------------------
#' 
#' Las dos funciones usadas para obtener los intervalos de confianza permiten **comparar la tendencia central de un conjunto de datos con un valor preestablecido**. Como se ha indicado antes, la prueba de t require normalidad de los datos; la prueba de Wilcoxon, no.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
x <- rnorm(10, mean = 2, sd = .5)
x

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
t.test(x, mu = 1.5)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
wilcox.test(x ^ 3, mu = 8)

#' 
#' ------------------------------------------------------------------------
#' 
#' Para **comparar la localización** de dos conjuntos de datos, deben tenerse en cuenta los siguientes aspectos
#' 
#' -   ¿los datos de ambas muestras están asociados?
#' -   ¿los datos estan normalmente distribuidos?
#' -   ¿tienen ambos conjuntos de datos la misma dispersión?
#' 
#' ------------------------------------------------------------------------
#' 
#' Si los datos están **asociados**, es decir que existen pares de datos en ambas muestras tales que comparten fuentes de variación, entonces
#' 
#' -   si ambas **muestras están normalmente distribuidas**, o la diferencia entre conjuntos de datos está normalmente distribuido, se usa la prueba de t sobre la diferencia de los valores entre ambos conjuntos de datos. En R, puede calcularse la diferencia y luego aplicar la función `t.test` o bien usar directamente la función `t.test` con la opción `paired = TRUE`,
#' -   si los datos **no están normalmente distribuidos** se usa la prueba de Wilcoxon --`wilcox.test` sobre la resta o con la opción `paired = TRUE` en R--.
#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
x <- rnorm(10, mean = 3, sd = 1)
y <- rnorm(10, mean = 3.5, sd = 1)
list(x = x,y = y)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
t.test(x - y)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
t.test(x, y, paired = TRUE)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
x <- rlnorm(10, meanlog = 3, sdlog = 1)
y <- rlnorm(10, meanlog = 3.5, sdlog = 1)
list(x = x,y = y)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
wilcox.test(x - y)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
wilcox.test(x, y, paired = TRUE)

#' 
#' ------------------------------------------------------------------------
#' 
#' Si los datos **no están asociados**, entonces
#' 
#' -   si ambas **muestras están normalmente distribuidas**, o el conjunto de datos --escalados y centrados en cada muestra-- está normalmente distribuido, se usa la prueba de t. En R, se usa la función `t.test`; la opción `var.equal` permite escoger entre la prueba calculada con la varianza de los datos, y la aproximación de Welch, que se usa cuando la varianzas no son iguales.
#' -   si los datos **no están normalmente distribuidos** se usa la prueba U de Mann-Whitney --`wilcox.test` en R--.
#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
x <- rnorm(10, mean = 3, sd = 1)
y <- rnorm(14, mean = 3.5, sd = 2)
list(x = x,y = y)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
t.test(x, y, var.equal = FALSE)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
x <- rlnorm(10, meanlog = 3, sdlog = 1)
y <- rlnorm(14, meanlog = 3.5, sdlog = 2)
list(x = x,y = y)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
wilcox.test(x, y)

#' 
#' ------------------------------------------------------------------------
#' 
#' Para la **comparación de tres o más grupos**, las técnicas más habituales son
#' 
#' -   si los datos de los distintos conjuntos de datos están **normalmente distribuidos** y presentan dispersiones comparables --homocedasticidad--, la técnica recomendada es el análisis de la varianza --`aov` en R--;
#' -   si los datos de los distintos conjuntos de datos están **normalmente distribuidos** y no presentan dispersiones comparables, la técnica recomendada es el análisis de la varianza de Welch --`oneway.test` en R--;
#' -   si los **datos no están normalmente distribuidos** pero las distribuciones tienen formas parecidas y se mantiene la homocedaticidad, se usa la prueba de Kruskal-Wallis --en R, `krukal.test`--.
#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
x <- c(rnorm(10, mean = 3, sd = 1),
       rnorm(8, mean = 3.2, sd = 1),
       rnorm(10, mean = 4, sd = 1))
g <- factor(c(rep(1,10), rep(2,8), rep(3,10)))
test <- aov(x ~ g)
summary(test)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
x <- c(rnorm(10, mean = 3, sd = 1),
       rnorm(8, mean = 3.2, sd = 4),
       rnorm(10, mean = 4, sd = 2))
g <- factor(c(rep(1,10), rep(2,8), rep(3,10)))
oneway.test(x ~ g)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
x <- c(runif(10, min = 2, max = 4),
       runif(8, min = 1.5, max = 3.5),
       runif(10, min = 3, max = 5))
g <- factor(c(rep(1,10), rep(2,8), rep(3,10)))
kruskal.test(x ~ g)

#' 
#' ------------------------------------------------------------------------
#' 
#' Tanto el análisis de la varianza como la prueba de Kruskal-Wallis solo permiten contrastar si algún par de grupos presentan valores de localización distintos, pero no permiten identificar cuál o cuáles.
#' 
#' Para identificar qué pares de grupos presentan diferencias estadísticamente significativos se usan pruebas *post-hoc*
#' 
#' -   En el caso del análisis de varianza (o Welch ANOVA), se puede usar la prueba de Tukey --`TukeyHSD` en R,
#' -   Para la prueba de Kruskal-Wallis, se puede usar la prueba de Dunn --`dunn.test` en el paquete `dunn.test` de R--
#' 
#' ------------------------------------------------------------------------
#' 
#' Gráficamente y aunque con menos objetividad, la comparación de la loclización puede hacerse usando gráficos de caja.
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## x <- c(rnorm(10, mean = 3, sd = 1),
##        rnorm(8, mean = 3.2, sd = 1),
##        rnorm(10, mean = 4, sd = 1))
## g <- factor(c(rep(1,10), rep(2,8), rep(3,10)))
## 
## xs <- data.frame(y = x, group = g)
## 
## boxplot(y~group, data=xs)
## 

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
x <- c(rnorm(10, mean = 3, sd = 1),
       rnorm(8, mean = 3.2, sd = 1),
       rnorm(10, mean = 4, sd = 1))
g <- factor(c(rep(1,10), rep(2,8), rep(3,10)))

xs <- data.frame(y = x, group = g)
  
boxplot(y~group, data=xs)

#' 
#' ## Ajuste de modelos
#' 
#' El ajuste de modelos se expresa de forma habitual en R mediante un objeto `formula`.
#' 
#' En esta notación, se indica la relación entre variables mediante una expresión de tres términos donde la variable dependiente se indica a la izquierda y las variables dependientes a la derecha.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
form1 <- y ~ x              # recta
class(form1)

#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## form1 <- y ~ log(x)         # logaritmo
## form1 <- y ~ poly(x,4)      # polinomio de grado 4
## form1 <- y ~ x + 0          # recta que pasa por el origen
## form1 <- y ~ I(x^.5)        # raíz cuadrada

#' 
#' ## Ajuste de modelos -- regresión lineal
#' 
#' El ajuste de un modelo lineal por mínimos cuadrados ordinarios se hace en R mediante la función `lm`.
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## fit1 <- lm(mpg ~ wt, data=mtcars)
## summary(fit1)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
fit1 <- lm(mpg ~ wt, data=mtcars)
summary(fit1)

#' 
#' ------------------------------------------------------------------------
#' 
#' Una vez ajustado el modelo, este puede usarse para predecir nuevos datos (o calcular los valores ajustados) mediante la función `predict`.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
df <- data.frame(x=mtcars$wt, y=mtcars$mpg,
  predict(fit1, newdata=mtcars, interval="prediction"))
head(df)

#' 
#' ------------------------------------------------------------------------
#' 
#' El modelo se puede visualizar gráficamente a partir de estas predicciones.
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## df2 <- data.frame(
##   wt = seq(min(mtcars$wt), max(mtcars$wt), length.out = 201),
##   predict(fit1,
##   newdata = data.frame(wt = seq(min(mtcars$wt), max(mtcars$wt),
##                              length.out = 201)),
##   interval="prediction"))
## 
## plot(range(df2$wt),range(c(df2$lwr,df2$upr)),type="n")
## lines(df2$wt,df2$fit,type="l",lwd=2)
## lines(df2$wt,df2$upr,col="grey",lwd=1)
## lines(df2$wt,df2$lwr,col="grey",lwd=1)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
df2 <- data.frame(
  wt = seq(min(mtcars$wt), max(mtcars$wt), length.out = 201),
  predict(fit1, 
  newdata = data.frame(wt = seq(min(mtcars$wt), max(mtcars$wt),
                             length.out = 201)), 
  interval="prediction"))

plot(range(df2$wt),range(c(df2$lwr,df2$upr)),type="n")
lines(df2$wt,df2$fit,type="l",lwd=2)
lines(df2$wt,df2$upr,col="grey",lwd=1)
lines(df2$wt,df2$lwr,col="grey",lwd=1)

#' 
#' ## Ajuste de modelos -- comprobación
#' 
#' Para comprobar la adecuación del modelo y del método de ajuste usado (OLS) deben analizarse los residuales.
#' 
#' Estos deben ser aleatorios, estar normalmente distribuidos, mostrar homocedasticidad y no mostrar evidencia de puntos especialmente influyentes.
#' 
#' Estas comprobaciones pueden hacerse a partir de los datos de los residuales y las pruebas estadísticas oportunas.
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## fit1$residuals

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
fit1$residuals

#' 
#' ------------------------------------------------------------------------
#' 
#' Gráficamente, los mismas comprobaciones pueden hacerse a partir de la representación gráfica del modelo.
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## plot(fit1, which=1:6)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
plot(fit1, which=1)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
plot(fit1, which=2)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
plot(fit1, which=3)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
plot(fit1, which=4)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
plot(fit1, which=5)

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = FALSE----------------------------------------------------------------
plot(fit1, which=6)

#' 
#' ------------------------------------------------------------------------
#' 
#' Estas mismas ideas son la base de la aplicación en R de las técnicas de ajuste multilineal (`lm`), ajuste lineal generalizado (`glm`) y ajuste no lineal por mínimos cuadrados (`nls`).
#' 
#' 
#' # Manipulación avanzada de tablas de datos
#' 
#' ## Tabla de datos o *data frame*
#' 
#' Com se ha viso anteriormente, el *data frame* es el tipo de dato más usado para almacenar tablas de datos.
#' 
#' A menudo, para poder realizar gráficos y/o aplicar distintos procedimientos estadísticos, es necesario manipular la tabla de datos. A esto dedicaremos este apartado.
#' 
#' ------------------------------------------------------------------------
#' 
#' Entre las operaciones habituales que haremos sobre una tabla de datos, hay
#' 
#' -   renombrar filas o columnas,
#' -   añadir columnas o filas,
#' -   segmentar,
#' -   eliminar filas o columnas,
#' -   cambiar el formato de un conjunto de datos,
#' -   crear tablas de datos de resumen, y
#' -   unir tablas
#' 
#' Terminaremos este bloque introduciendo los funciones básicas del paquete `dplyr` que facilita la realización de algunas de estas operaciones.
#' 
#' ------------------------------------------------------------------------
#' 
#' Partimos de una tabla de datos sintética...
#' 
## ---- echo = TRUE-----------------------------------------------------------------
df <- data.frame(1:5, letters[1:5], c(rep("a", 3), rep("b", 2)))
df

#' 
#' ## Renombrar filas o columnas
#' 
## ---- echo = TRUE-----------------------------------------------------------------
colnames(df) <- c("var1", "var2", "var3") 
rownames(df) <- paste("subject00", 1:5, sep = "")
df

#' 
#' ## Añadir columnas o filas
#' 
## ---- echo = TRUE-----------------------------------------------------------------
df2 <- cbind(df, rnorm(5)) # añadir un vector al data frame
df2$var5 <- 5:1 # assignando valores a una nueva variable
df2

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
df2 <- rbind(df, list(6, "e", "b"))
df2

#' 
#' ## Segmentar
#' 
#' Existen tres formas básicas de seleccionar filas o columnas de una tablas de datos
#' 
#' -   mediante índices numéricos,
#' -   usando los nombres de filas y columnas, y
#' -   mediante vectores lógicos
#' 
#' ## Segmentar -- mediante índices
#' 
## ---- echo = TRUE-----------------------------------------------------------------
df[1:3,]
df[,c(1,3)]

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
df[-3,-2]

#' 
#' ## Segmentar -- usando nombres
#' 
## ---- echo = TRUE-----------------------------------------------------------------
df[,"var2"]
df$var3
df[,c("var2","var3")]

#' 
#' ## Segmentar -- mediante vectores lógicos
#' 
## ---- echo = TRUE-----------------------------------------------------------------
df[c(T,T,F,T,F), c(T,F,T)]
df[df[,1] == 3 | df[,3] == "b",]

#' 
#' 
#' 
#' 
#' ## Cambiar el formato de un conjunto de datos
#' 
#' Un mismo conjunto de datos puede representar en una tabla de datos de acuerdo con distintas oragnizaciones o formatos.
#' 
#' Se denomina formato ordenado, *tidy*, aquella organización en la que las filas representan individuos, las columnas, variables de medidas, y las intersecciones los valores de dichas medidas.
#' 
#' ------------------------------------------------------------------------
#' 
#' Cuando existen más de una variable que corresponde al mismo tipo de medida, entonces los datos pueden organizarse de dos formas
#' 
#' -   usando una columna distinta para cada variable, formato ancho y *tidy*, o
#' -   usando una columna para los valores y otra para indicar cuál es la variable representada en esta fila.
#' 
#' El formato más adecuado dependerá de los análisis y visualizaciones que quieran hacerse.
#' 
#' ------------------------------------------------------------------------
#' 
#' Partiremos de un subconjunto de `flights`...
#' 
## ---- echo = TRUE, eval = FALSE---------------------------------------------------
## if(!require("nycflights13")) {
##   install.packages("nycflights13")
##   library("nycflights13")
## }

#' 
## ---- echo = FALSE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("nycflights13")) {
  install.packages("nycflights13", repos="https://cloud.r-project.org/",
         quiet=TRUE, type="binary")
  library("nycflights13")
}

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
fl_ny2ws_W <- flights[flights$dest %in% c("IAD","BWI"),
                    c("origin","dest","carrier","arr_delay","dep_delay")]
fl_ny2ws_W <- cbind(key = 1:nrow(fl_ny2ws_W), fl_ny2ws_W)
head(fl_ny2ws_W)

#' 
#' ## Cambiar el formato de un conjunto de datos -- ancho a largo
#' 
## ---- echo = TRUE-----------------------------------------------------------------
fl_ny2ws_L <- rbind(
  cbind(edge = rep("origin", nrow(fl_ny2ws_W)), fl_ny2ws_W[,c(1,4)],
        airport = fl_ny2ws_W[,2], delay = fl_ny2ws_W[,6]),
  cbind(edge = rep("dest", nrow(fl_ny2ws_W)), fl_ny2ws_W[,c(1,4)],
        airport = fl_ny2ws_W[,3], delay = fl_ny2ws_W[,5]))
colnames(fl_ny2ws_L) <- c("edge","key","carrier","airport","delay")

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
head(fl_ny2ws_L)
tail(fl_ny2ws_L)

#' 
#' ## Cambiar el formato de un conjunto de datos -- largo a ancho
#' 
## ---- echo = TRUE-----------------------------------------------------------------
fl_ny2ws_W2p1 <- fl_ny2ws_L[fl_ny2ws_L$edge=="origin",]
fl_ny2ws_W2p2 <- fl_ny2ws_L[fl_ny2ws_L$edge=="dest",]
fl_ny2ws_W2p1 <- fl_ny2ws_W2p1[order(fl_ny2ws_W2p1$key),-1]
fl_ny2ws_W2p2 <- fl_ny2ws_W2p2[order(fl_ny2ws_W2p2$key),-c(1,2)]
fl_ny2ws_W2 <- cbind(fl_ny2ws_W2p1,fl_ny2ws_W2p2[,-1])
colnames(fl_ny2ws_W2) <- c("key", "carrier", "origin", "dep_delay",
                           "dest","arr_delay")


#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
head(fl_ny2ws_W2)
tail(fl_ny2ws_W2)

#' 
#' ## Eliminar filas o columnas
#' 
#' La forma más habitual de eliminar filas o columnas es segmentando la tabla de datos. Sin embargo, una columna también puede eliminarse asignando la misma a `NULL`.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
df <- data.frame(1:5, letters[1:5], c(rep("a", 3), rep("b", 2)))
colnames(df) <- c("var1", "var2", "var3") 
rownames(df) <- paste("subject00", 1:5, sep = "")
df <- df[-2,]
df$var2 <- NULL

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
df

#' 
#' ------------------------------------------------------------------------
#' 
#' Si lo que se desea es eliminar una variable del entorno de trabajo entonces se usa la función `rm`.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
rm(df)

#' 
#' ## Creación de resúmenes a partir de datos en formato ancho
#' 
## ---- echo = TRUE-----------------------------------------------------------------
sum_fl_ny2ws <- data.frame(edge=c("origin","dest"))

sum_fl_ny2ws$mean_delay <- apply(fl_ny2ws_W[,c("dep_delay","arr_delay")],2,
                                 mean,na.rm=TRUE)
sum_fl_ny2ws$r_delay <- apply(fl_ny2ws_W[,c("dep_delay","arr_delay")],2,
                      function(x) diff(range(x,na.rm=TRUE)))
sum_fl_ny2ws$s_delay <- apply(fl_ny2ws_W[,c("dep_delay","arr_delay")],2,
                              sd,na.rm=TRUE)
sum_fl_ny2ws

#' 
#' ## Creación de una tabla de resumen a partir de datos en formato largo
#' 
## ---- echo = TRUE-----------------------------------------------------------------
sum_fl_ny2ws <- data.frame(edge=c("origin","dest"))
sum_fl_ny2ws$mean_delay <- by(fl_ny2ws_L$delay,fl_ny2ws_L$edge,
                              mean,na.rm=TRUE)
sum_fl_ny2ws$r_delay <- by(fl_ny2ws_L$delay,fl_ny2ws_L$edge,
                           function(x) diff(range(x,na.rm=TRUE)))
sum_fl_ny2ws$s_delay <- by(fl_ny2ws_L$delay,fl_ny2ws_L$edge,sd,na.rm=TRUE)
sum_fl_ny2ws

#' 
#' ## `dplyr`
#' 
#' Estas mismas operaciones que se han comentado más arriba, se pueden llevar a cabo también a partir de un sistema coherente de métodos creado como una gramática para la manipulación de datos, el paquete `dplyr` --que forma parte de `tidyverse`--.
#' 
#' <http://dplyr.tidyverse.org/>
#' 
## ---- echo = FALSE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("tidyverse")) {
  install.packages("tidyverse", repos="https://cloud.r-project.org/",
         quiet=TRUE, type="binary")
  library("tidyverse")
}

#' 
#' ------------------------------------------------------------------------
#' 
#' En `dplyr` las operaciones básicas, se realizan mediante cuatro métodos:
#' 
#' -   `select`: segmentar columnas,
#' -   `filter`: segmentar filas,
#' -   `mutate`: crear nuevas columnas, y
#' -   `arrange`: ordenar.
#' 
#' Estas se encadenan usando el operador pipe, `%>%`.
#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE-----------------------------------------------------------------
df <- flights %>% dplyr::select(origin, dest, arr_delay) %>% 
  filter(origin == "LGA" & (dest == "IAD" | dest == "BWI")) %>%
  mutate(arr_delay_h=arr_delay/60) %>% 
  arrange(-arr_delay_h)
df

#' 
#' ------------------------------------------------------------------------
#' 
#' Para la creación de resúmenes a partir de tablas en formato largo, es muy útil y cómoda la combinación `group_by` y `summarize`.
#' 
## ---- echo = TRUE-----------------------------------------------------------------
df %>% group_by(dest) %>% summarise(mean_delay = mean(arr_delay, na.rm=TRUE))

#' 
#' ------------------------------------------------------------------------
#' 
#' Por último, las conversiones entre formatos también pueden hacerse a partir de las funciones `spread` y `gather` del paquete `tidyr` --incluido en `tidyverse` y que se integra de forma natural en la gramática propuestas por `dplyr`--.

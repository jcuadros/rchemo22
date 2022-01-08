#' ---
#' pagetitle: "R for Chemistry Data Analysis and Chemometrics. Visual Analytics in R"
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
#' editor_options: 
#'   chunk_output_type: console
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
#' <div style="font-size:1.4em;font-weight:500;color:#333333;">Visual Analytics in R</div>
#' <div style="font-size:1.2em;margin-top:40px;color:#333333;">Jordi Cuadros, Vanessa Serrano</div>
#' <div style="margin-top:80px;color:#333333;">January 2022</div>
#' 
#' 
#' # Package Management
#' 
#' ## Installation and Loading of R Packages
#' 
#' R has many contributed to solve specific needs and problems in several disciplines. These are usually hosted in repositories, such as CRAN (<https://cran.r-project.org/>), Bioconductor (<https://www.bioconductor.org/>), R-Forge (<https://r-forge.r-project.org/>), ROpenSci (<https://ropensci.org/>), Stan packages (<http://mc-stan.org/r-packages/>)... Yet other packages are hosted in version control websites as GitHub or GitLab.
#' 
#' ----
#' 
#' In all cases, to make use of a package, it has to be installed and loaded in memory. 
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## installed.packages()[,1] # Lists the installed packages
## (.packages())            # Lists the packages currently loaded

#' 
#' Some relevant resources to find useful packages are
#' <https://www.rdocumentation.org/> and <https://cran.rstudio.com/web/views/>.
#' 
#' ----
#' 
#' To install a package, we use
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## install.packages("dplyr")

#' To load a package in memory, we use
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## library("dplyr")

#' 
#' 
#' Alternatively, RStudio offers a menu-based package management feature.
#' 
#' ----
#' 
#' In a script, to avoid repeated installations of a package and to ensure its availability, we can use
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## if(!require("dplyr")) {
##   install.packages("dplyr",
##                    repos="https://cloud.r-project.org/")
##   library("dplyr")
## }

#' 
#' ----
#' 
#' If the package is installed but not loaded, we can still access its functions by using the `::`operator.
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## dplyr::band_instruments

#' 
#' ----
#' 
#' To browse the help files for a packages, run
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## help(package="dplyr")

#' 
#' Additionnaly some packages have a vignette and some demo code that can be inspected with the`vignette()` and `demo()` functions.
#' 
#' 
#' Last, to learn more about packages and package development, you may want to read <https://r-pkgs.org/index.html>.
#' 
#' 
#' # Data Import and Export
#' 
#' ## Data Files {.small}
#' 
#' It is often useful to save and load data for its transfer among different application or just to keep a copy of it for back-up or later use. 
#' 
#' Although can read and write many different types of files, we will focus on the most relevant format for statistics, chemistry and chemometrics. We will discuss here how to handle
#' 
#' - text files,
#' - structured data files, such as XML and JSON, 
#' - files from other computational programs (i.e. Excel and Matlab), 
#' - RData files, the default format for R, and
#' - some data files which are specific for chemistry and chemical analysis.
#' 
#' 
#' ## A Previous: Set the Working Directory
#' 
#' Before working with external files, it is important to set the working directory. `setwd` and `getwd` functions can be used to set the working directory from the R code, although its use is not recommended in scripts that may need to be shared.
#' 
#' RStudio allows setting the working directory from the *Session* menu and form the *Files* pane.
#' 
#' 
#' ## Text Files
#' 
#' Let's start from the `wine` data set in the `FactoMineR` package. 
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## if(!require("FactoMineR")) {
##   install.packages("FactoMineR",
##                    repos="https://cloud.r-project.org/")
##   library("FactoMineR")
## }
## data("wine")
## help("wine")
## str(wine)
## 

#' 
## ---- echo = FALSE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("FactoMineR")) {
  install.packages("FactoMineR",
                   repos="https://cloud.r-project.org/")
  library("FactoMineR")
}
data("wine")

#' 
#' ----
#' 
## ---- echo = FALSE---------------------------------------------------------
str(wine)

#' 
#' ----
#' 
#' To write a delimited text file from a data frame, we use the `write.table` function or any of its derivative functions. To read, `read.table` is to be used.
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## write.table(wine, file="wine.csv", sep=",", dec=".",
##             quote=TRUE, fileEncoding="UTF-8", row.names=FALSE)

#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## wine2 <- read.table("wine.csv", sep=",", dec=".",
##            quote="\"", fileEncoding="UTF-8", header=TRUE)

#' 
#' In case we want a tab-delimited text file, we will use `\t` for the argument `sep`.
#' 
#' ----
#' 
#' Other options exist for managing text-files.
#' 
#' In case the format is not well-defined, we may want to use `readLines` and `writeLines` to read or write the file line by line.
#' 
#' If the file is large, `read_table` from the `readr` package or `vroom` from the `vroom` package are usually faster.
#' 
#' 
#' ## Structured Text Files -- XML
#' 
#' An example of an XML file containing the structure of caffeine can be found at <https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/2519/record/XML>. It can be read in R with `read_xml` from the `xml2` package.
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## if(!require("xml2")) {
##   install.packages("xml2",
##                    repos="https://cloud.r-project.org/")
##   library("xml2")
## }

#' 
## ---- echo = FALSE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("xml2")) {
  install.packages("xml2",
                   repos="https://cloud.r-project.org/")
  library("xml2")
}

#' 
#' ## {.small}
#' 
#' Information can be extracted using XPath selectors, <https://www.w3schools.com/xml/xpath_intro.asp>.
#' 
## ---- echo = TRUE----------------------------------------------------------
xmlObj <- read_xml(
  "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/2519/record/XML",
  options = "RECOVER")
xml_ns_strip(xmlObj)  # Important!
xmlObj_Name <- xml_find_all(xmlObj,xpath=
  ".//PC-InfoData//*[text()='Traditional']/../../..//PC-InfoData_value_sval/text()")
as.character(xmlObj_Name)

#' 
#' ----
#' 
#' NMR spectra can be obtained in CML (Chemistry Markup Language, <http://www.xml-cml.org/spec/>) format from <https://www.nmrshiftdb.org>.
#' 
#' For example, the 13C-NMR spectra for butane is available at <https://www.nmrshiftdb.org/NmrshiftdbServlet/nmrshiftdbaction/searchorpredict/smiles/CCCC/spectrumtype/13C>.
#' 
#' 
#' ## Structured Text Files -- JSON
#' 
#' JSON files are another common format for data storage and transmission. In R, the `jsonlite`package offers functions to import and convert this files to list or data frames.
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## if(!require("jsonlite")) {
##   install.packages("jsonlite",
##                    repos="https://cloud.r-project.org/")
##   library("jsonlite")
## }

#' 
## ---- echo = FALSE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("jsonlite")) {
  install.packages("jsonlite",
                   repos="https://cloud.r-project.org/")
  library("jsonlite")
}

#' 
#' ## {.small}
#' 
## ---- echo = TRUE----------------------------------------------------------
jsBz <- fromJSON(
  "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/241/JSON")
jsBz_NI <- jsBz$Record$Section[
  jsBz$Record$Section$TOCHeading=="Names and Identifiers",]
jsBz_NI_OI <- jsBz_NI$Section[[1]][
  jsBz_NI$Section[[1]]$TOCHeading=="Other Identifiers",]
jsBz_NI_OI$Section[[1]]$TOCHeading
CAS <- jsBz_NI_OI$Section[[1]][jsBz_NI_OI$Section[[1]]$TOCHeading=="CAS",]
table(unlist(CAS$Information[[1]]$Value$StringWithMarkup))

#' 
#' 
#' ## Computational Programs -- Excel
#' 
#' Let's use the `wine` data set again.
#' 
## ---- echo = TRUE----------------------------------------------------------
str(wine)

#' 
#' ----
#' 
#' To write a data.frame to Excel, we can use the `write_xlsx` from the `writexl` package. To read an Excel file, `read_xlsx` from the `readxl` package can be used.
#' 
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## if(!require("writexl")) {
##   install.packages("writexl",
##                    repos="https://cloud.r-project.org/")
##   library("writexl")
## }
## 
## if(!require("readxl")) {
##   install.packages("readxl",
##                    repos="https://cloud.r-project.org/")
##   library("readxl")
## }

#' 
## ---- echo = FALSE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("writexl")) {
  install.packages("writexl",
                   repos="https://cloud.r-project.org/")
  library("writexl")
}

if(!require("readxl")) {
  install.packages("readxl",
                   repos="https://cloud.r-project.org/")
  library("readxl")
}

#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## write_xlsx(list(wineSheet = wine),
##                   path = "wine.xlsx")
## wineExcel <- readxl::read_xlsx(path = "wine.xlsx")

#' 
#' ## Computational Programs -- Matlab
#' 
#' To open and write Matlab data files, we can use the `R.matlab` package. Relevant functions are called `readMat` and `writeMat`.
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## if(!require("R.matlab")) {
##   install.packages("R.matlab",
##                    repos="https://cloud.r-project.org/")
##   library("R.matlab")
## }

#' 
## ---- echo = FALSE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("R.matlab")) {
  install.packages("R.matlab",
                   repos="https://cloud.r-project.org/")
  library("R.matlab")
}

#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## writeMat("wine.mat", wine = wine)
## wineMat <-readMat("wine.mat")
## wineMat <- data.frame(as.character(unlist(wineMat[[1]][[1]])),
##                       as.character(unlist(wineMat[[1]][2])),
##                       as.data.frame(wineMat[[1]][3:31]))
## 

#' 
#' Reading a Matlab data file returns a list that will need *ad hoc* manipulation.
#' 
#' 
#' ## RData Files
#' 
#' To save and read data in the R data format, we use `save` and `load`. These allow storing and restoring any set of variables of the environment. When loaded, variables are recovered with the same names they had on saving.  
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## save(wine, file="wine.rda")

#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## print(load("wine.rda"))          # print shows the name of the loaded objects
## wineR <- get(load("wine.rda"))   # get allows storing the loaded information
##                                  # with a different name

#' 
#' ## Chemical Data Formats {.small}
#' 
#' Besides the general data file formats already discussed, there are several formats used specifically to store chemical information. Some significant ones are 
#' 
#' - Chemical table files, as MOL or SDF, <https://en.wikipedia.org/wiki/Chemical_table_file>,
#' - JCAMP-DX, <http://jcamp-dx.org/>,
#' - NIST MSP format, <https://chemdata.nist.gov/mass-spc/ms-search/docs/Ver20Man_11.pdf>,
#' - AnIML, <https://www.animl.org/>,
#' - ICARTT, <https://www-air.larc.nasa.gov/missions/etc/IcarttDataFormat.htm>, and
#' - Allotrope Data Format, <https://www.allotrope.org/>.
#' 
#' More information on chemical data files can be found at <https://en.wikipedia.org/wiki/Chemical_file_format>, <https://en.wikipedia.org/wiki/Mass_spectrometry_data_format>, <http://unichrom.com/chrom/uc-ffe.shtml> or <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3518119/>.
#' 
#' 
#' ## Chemical Data Formats -- Chemical Table File
#' 
#' Chemical Table Files are a set of text-based file types designed to store molecular information, *e.g.* positions of the atoms and connection tables. Some of this format allow storing additional information such as conformers, additional properties and identifiers, or spectra.
#' 
#' Currently the simplest way to read a MOL or SDF file in R is the `read.SDFset` function in the `ChemmineR` Bioconductor package.
#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## if (!require("BiocManager", quietly = TRUE))
##   install.packages("BiocManager",
##                    repos="https://cloud.r-project.org/")
## 
## if (!require("ChemmineR", quietly = TRUE)) {
##   BiocManager::install("ChemmineR", update=FALSE)
##   library("ChemmineR")
## }

#' 
## ---- echo = FALSE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager",
                   repos="https://cloud.r-project.org/")

if (!require("ChemmineR", quietly = TRUE)) {
  BiocManager::install("ChemmineR", update=FALSE)
  library("ChemmineR")
}

#' 
## ---- echo = TRUE, eval=FALSE----------------------------------------------
## sdf1 <- read.SDFset(
##   "https://cactus.nci.nih.gov/chemical/structure/CCCC/file?format=sdf")
## draw_sdf(sdf1[[1]], filename=NULL)

#' 
#' ----
#' 
## ---- echo = FALSE---------------------------------------------------------
sdf1 <- read.SDFset("https://cactus.nci.nih.gov/chemical/structure/CCCC/file?format=sdf")
draw_sdf(sdf1[[1]], filename=NULL)

#' 
#' 
#' ----
#' 
## ---- echo = TRUE----------------------------------------------------------
header(sdf1)
MW(sdf1[[1]])


#' 
#' ## Chemical Data Formats -- JCAMP-DX
#' 
#' JCAMP-DX is an text-based open standard to store and distribute spectral data. Currently the best way to read these files in R is the `readJDX` function in the `readJDX` package.
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## if (!require("readJDX")) {
##   install.packages("readJDX",
##                    repos="https://cloud.r-project.org/")
##   library("readJDX")
## }

#' 
## ---- echo = FALSE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if (!require("readJDX")) {
  install.packages("readJDX",
                   repos="https://cloud.r-project.org/")
  library("readJDX")
}

#' 
#' ----
#' 
#' ### Acetone IR
#' 
## ---- echo = TRUE, eval=FALSE----------------------------------------------
## jdx1 <- readJDX ("../data/67-64-1-IR.jdx")
## plot(jdx1$Acetone$x,jdx1$Acetone$y, type="l",
##      xlab="wave number", ylab="T")

#' 
#' ::: {.bibref}
#' JCAMP-DX downloaded from <https://webbook.nist.gov/cgi/cbook.cgi?ID=C67641>, in NIST Chemistry WebBook, NIST Standard Reference Database Number 69, Eds. P.J. Linstrom and W.G. Mallard, National Institute of Standards and Technology, Gaithersburg MD, 20899, https://doi.org/10.18434/T4D303, (on January 5, 2022).  
#' :::
#' 
#' ----
#' 
## ---- echo = FALSE---------------------------------------------------------
jdx1 <- readJDX ("../data/67-64-1-IR.jdx")
plot(jdx1$Acetone$x,jdx1$Acetone$y, type="l",
     xlab="wave number", ylab="T")

#' 
#' 
#' ## Chemical Data Formats -- NIST MSP
#' 
#' MSP is a text-based NIST-promoted format to store collections of spectra. It's a common format in the metabolomics community. 
#' 
#' Downloadable spectral databases can be found at
#' 
#' - <http://prime.psc.riken.jp/compms/msdial/main.html#MSP>
#' - <https://mona.fiehnlab.ucdavis.edu/downloads>
#' - <https://chemdata.nist.gov/dokuwiki/doku.php?id=peptidew:cdownload>
#' 
#' ----
#' 
#' While a functional package to read MSP files doesn't seem to be available, these files can be read with standard functions for text files.
#' 
#' 
## ---- echo = TRUE----------------------------------------------------------
    msp1 <- readLines("../data/MSMS-Neg-MassBankEU.msp")
    msp1 <- paste(msp1, collapse="\n")
    msp1 <- unlist(strsplit(msp1, "\n\n"))

    msp1_1 <- unlist(strsplit(msp1[1],"\n"))
    spectrum <- msp1_1[1:length(msp1_1) > which(substr(msp1_1,1,10)=="Num Peaks:")]
    spectrum <- read.table(text=spectrum, sep="\t")
    colnames(spectrum) <- c("m_z","int")
    spectrum

#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
##     spectrum <- data.frame(m_z=rep(spectrum$m_z,each=3),
##                           int=rep(spectrum$int,each=3),
##                           i=1:3)
##     spectrum$int[spectrum$i!=2] <- 0
## 
##     plot(x=spectrum$m_z, y=spectrum$int, xlim=c(40,300),
##          type="l", xlab="m/z", ylab="")

#' 
#' ----
#' 
## ---- echo = FALSE---------------------------------------------------------
    spectrum <- data.frame(m_z=rep(spectrum$m_z,each=3),
                          int=rep(spectrum$int,each=3),
                          i=1:3)
    spectrum$int[spectrum$i!=2] <- 0

    plot(x=spectrum$m_z, y=spectrum$int, xlim=c(40,300),
         type="l", xlab="m/z", ylab="")

#' 
#' ::: {.bibref}
#' `MSMS-Neg-MassBankEU.msp` downloaded from http://prime.psc.riken.jp/compms/msdial/main.html#MSP 
#' 
#' 
#' ## YOUR TURN {data-background="#eeffcc"}
#' 
#' 1.    Import the classical wine data of Forina et al. It is available at <http://archive.ics.uci.edu/ml/datasets/wine>.
#' 2.    Load the data from the NIR spectra for diesel fuels which is described and available at <https://eigenvector.com/resources/data-sets/>.
#' 3.    Import the data from mango aroma components available at <https://figshare.com/articles/dataset/Mango_Mangifera_indica_Aroma_Discriminate_Cultivars_and_Ripeness_Stages/14303913>. Warning: it needs some cleaning!
#' 4.    Extract the MS of lactic acid from the Golm database, at <http://gmd.mpimp-golm.mpg.de/download/> and plot it. It is available in the MDN35 ALK data set.
#' 
#' 
#' # Advanced manipulation of data frames (part 2 - the tidy way)
#' 
#' ## Data frame {.small}
#' 
#' As already presented, we usually have to arrange data we have in a data frame by
#' 
#' -   renaming columns or rows,
#' -   adding columns or rows,
#' -   segmenting (selecting columns or rows),
#' -   removing columns or rows,
#' -   creating a data frame from combining vectors, 
#' -   joining data frames,
#' -   reshaping data tables, and
#' -   summarizing or aggregating the data.
#' 
#' Many of these operations can be performed by using a set of functions that operate in a single common paradigm, sometimes defined as grammar for data manipulation, implemented in the `dplyr` and `tidyr` R packages.
#' 
#' ------------------------------------------------------------------------
#' 
#' We usually install and load both packages jointly with some other packages that work nicely together. These are grouped in the `tidyverse` package.
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## if (!require("tidyverse")) {
##   install.packages("tidyverse",
##                    repos="https://cloud.r-project.org/")
##   library("tidyverse")
## }

#' 
## ---- echo = FALSE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if (!require("tidyverse")) {
  install.packages("tidyverse",
                   repos="https://cloud.r-project.org/")
  library("tidyverse")
}

#' 
#' 
#' ## Functions
#' 
#' R supports a functional approach to programming where functions can be defined and transferred as arguments. Many calculations in R are integrated as functions. Functions are a critical concept in a tidy data management.
#' 
#' Functions are defined with the `function` reserved word (or `\` since R 4.1.0).
#' 
## ---- echo = TRUE----------------------------------------------------------
sq <- function(x) {x ^ 2}
pow <- \(x,y) {x ^ y}

sq(4)
pow(3,4)

#' 
#' 
#' ## Pipe Operators, `%>%` and `|>`
#' 
#' Pipes are a syntax variant for transferring arguments to functions. It allows for a more process oriented organization of the code. R includes two types of pipes: `%>%` form the `magrittr` package and `|>` included in R 4.1.0.
#' 
## ---- echo = TRUE----------------------------------------------------------
sq(4)
4 %>% sq    # With %>%, parenthesis are not required when there is no arguments
4 |> sq()

#' 
#' ----
#' 
## ---- echo = TRUE----------------------------------------------------------
pow(sq(sq(2)),3)
2 %>% sq %>% sq %>% pow(3)

#' 
#' ----
#' 
#' By default, the left operand of the pipe is sent as the first argument of the function on the right. It can be sent to a different position in the calculation by including a dot (`.`) in the right-hand side of the pipe.
#' 
## ---- echo = TRUE----------------------------------------------------------
pow(sq(sq(2)),sq(sq(2))/100)
2 %>% sq %>% sq(.) %>% pow(.,./100)

#' 
#' 
#' ## `dplyr` and `tidyr`
#' 
#' `dplyr` and `tidyr` are two packages of the `tidyverse` designd for manipulating and tidying data. A full description of both packages can be found at <https://dplyr.tidyverse.org/> and <https://tidyr.tidyverse.org/>.
#' 
#' Current cheatsheets for both packages can be found at <https://www.rstudio.com/resources/cheatsheets/>.
#' 
#' ----
#' 
#' Data manipulation in these packages in commonly done using a limited set of functions. The most important are
#' 
#' - `select` and `filter`, for subsetting,
#' - `mutate`, for creating or updating calculated columns,
#' - `arrange`, for sorting,
#' - `pivot_wider` and `pivot_longer`, for reshaping, and
#' - `group by` and `summarize`, for aggregating. 
#' 
#' 
#' ----
#' 
#' Let's see how we operate in this logic, using as an example the `wine` data set.
#' 
## ---- echo = TRUE----------------------------------------------------------
str(wine) 

#' 
#' ## Renaming columns or rows
#' 
## ---- echo = TRUE----------------------------------------------------------
newColnames <- c("label", "soil", "odorBf", "aromaBf",
                 "fruityBf", "flowerBf", "spiceBf",
                 "visual", "nuance", "surface", "odorI",
                 "odorQ", "fruity", "flower", "spice",
                 "plante", "phenolic", "aromaI", "aromaP",
                 "aromaQ","attackI", "acidity", "astringency",
                 "alcohol", "balance", "smooth", "bitterness",
                 "intensity", "harmony", "overall", "typical")
wine2 <- wine %>% rename_with(function(x) newColnames)
colnames(wine2)


#' 
#' ----
#' 
## ---- echo = TRUE----------------------------------------------------------
wine2 <- wine2 %>% rownames_to_column() %>%
  mutate(rowname = paste0("W",sprintf("%0.2d",1:21)))
wine2$rowname

#' 
#' 
#' ## Adding columns or rows
#' 
## ---- echo = TRUE----------------------------------------------------------
wine2 <- wine2 %>% mutate(newCol=1:21)
head(wine2[,c(1,30:33)])

#' 
#' ------------------------------------------------------------------------
#' 
## ---- echo = TRUE----------------------------------------------------------
wine2 <- wine2 %>% add_row(wine2[21,])
tail(wine2[,c(1,30:33)])

#' 
#' ## Segmenting
#' 
## ---- echo = TRUE----------------------------------------------------------
wine2 %>% select(rowname,label,soil,overall) %>% 
  filter(overall>3.6)

#' 
#' ----
#' 
## ---- echo = TRUE----------------------------------------------------------
wine2 %>% rowwise() %>% mutate(worse=min(c_across(4:30))) %>% 
  filter(worse<1.5) %>% select(rowname, worse)

#' 
#' ----
#' 
## ---- echo = TRUE----------------------------------------------------------
wine2 %>% select(overall) %>% str
wine2 %>% pull(overall) %>% str

#' 
#' 
#' ## Removing columns or rows
#' 
## ---- echo = TRUE----------------------------------------------------------
wine2 <- wine2 %>% select(-newCol)
tail(wine2[,30:ncol(wine2)])

#' 
## ---- echo = TRUE----------------------------------------------------------
wine3 <- wine2 %>% filter(overall<3.5)
wine3 %>% pull(rowname)


#' 
#' 
#' ## Joining data frames
#' 
## ---- echo = TRUE----------------------------------------------------------
qualityLevel <- data.frame(level=1:4,
  qualityNom = c("Horrendous","Bad", "Good", "Excellent"))
wine2 <- wine2 %>% mutate(level=round(overall,0)) %>% 
  inner_join(qualityLevel)
head(wine2[,30:ncol(wine2)])

#' 
#' ## Reshaping data tables
#' 
#' The same information can be presented in many different ways, even in a data table. These are commonly referred as formats o shapes. 
#' 
#' A data frame is considered **tidy** when every row refers to an individual or observational unit and the columns to the variables being meaured.  
#' 
#' ----
#' 
#' When we have many variables that correspond to the same type of measurement, a data frame can be organized in two formats: a long format and a wide format.
#' 
#' In **wide** format, rows are observational units and there are different columns for the different variables being measured.
#' 
#' In **long** format, rows are observations and there is a column for the values and a column for for the variable being measured.
#' 
#' Different analyses and visualizations will require having the data organized in one or other format.
#' 
#' ----
#' 
#' The data set `wine` is available in a wide format. Each row is a type of wine.
#' 
#' Let's convert it to long...
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## wineL <- wine %>% rownames_to_column() %>%
##   pivot_longer(4:32,
##                names_to = "determination",
##                values_to = "value")
## str(wineL)

#' 
#' ----
#' 
## ---- echo = FALSE---------------------------------------------------------
wineL <- wine %>% rownames_to_column() %>% 
  pivot_longer(4:32,
               names_to = "determination",
               values_to = "value")
str(wineL)

#' 
## ---- echo = TRUE----------------------------------------------------------
head(wineL)

#' 
#' ----
#' 
#' Let's go back to wide...
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## wineW <- pivot_wider(wineL, id_cols = 1:3,
##                       names_from = "determination",
##                       values_from = "value")
## str(wineW)

#' 
#' ----
#' 
## ---- echo = FALSE---------------------------------------------------------
wineW <- pivot_wider(wineL, id_cols = 1:3,
                      names_from = "determination",
                      values_from = "value")
str(wineW)

#' 
#' ----
#' 
## ---- echo = TRUE----------------------------------------------------------
head(wineW,2)

#' 
#' 
#' ## Summarizing or Aggregating the Data.
#' 
#' ### For wide data frames
#' 
## ---- echo = TRUE----------------------------------------------------------
wineW %>% select(4:32) %>% summarize(across(
  everything(),list(mean=mean,sd=sd))) %>% t

#' 
#' ----
#' 
#' ### For long data frames
#' 
## ---- echo = TRUE----------------------------------------------------------
wineL %>% group_by(determination) %>% 
  summarize(mean=mean(value), sd=sd(value))

#' 
#' 
#' # Advanced Graphics in R
#' 
#' ## Grammar of Graphics (GoG)
#' 
#' The **grammar of graphics** is a theoretical approximation to the study of the graphics components. According to this, a graphical representation can be built from layers, by defining different elements
#' 
#' - the identification of the data to be represented,
#' - the association of the data variables to the graphical elements, and
#' - the definition of the graphical geometry. 
#' 
#' Additional elements of the representation can also be adjusted like scales, formats, annotations..., after the basic structure of the graphics is established.
#' 
#' ::: {.bibref}
#' Wilkinson, L. (2006). The grammar of graphics. Springer Science & Business Media.
#' :::
#' 
#' ## `ggplot2`
#' `ggplot2` is an R package which implements the grammar of graphics. A full reference can be found at <http://ggplot2.tidyverse.org>.
#' 
#' A current cheatsheet is available at <https://www.rstudio.com/resources/cheatsheets/>.
#' 
#' ::: {.bibref}
#' Wickham, H. (2010). A layered grammar of graphics. *Journal of Computational and Graphical Statistics*, 19(1), 3-28.
#' :::
#' 
#' ----
#' 
#' `ggplot2` is part of `tidyverse` package (although it can be installed and loaded independently).
#' 
## ---- echo = TRUE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("tidyverse")) {
  install.packages("tidyverse", repos="https://cloud.r-project.org/")
  library("tidyverse")
}

#' 
#' ## `ggplot2` -- Layers and Graphical Elements
#' 
#' In `ggplot2`, each graphical element constituting the presentation of a data set is a layer. One or more layers make a graph.
#' 
#' Each layer is defined by the specification of its elements. The main ones are the data, the aesthetic mapping --the relation between the data and the graphical elements--, and the geometry.
#' 
#' ----
#' 
#' In `ggplot2`, a plot is an R object a it can be stored. It is constructed additively, by summing different elements to the `gg` object constructed by the `ggplot2` function.
#' 
#' For example,
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## plot1 <- ggplot(data = anscombe,
##         mapping = aes(x = x1, y = y1))  # Data and aesthetic mapping
## plot1 <- plot1 + geom_point()       # Geometry
## plot1

#' 
#' ----
#' 
## ---- echo = FALSE---------------------------------------------------------
plot1 <- ggplot(data = anscombe,
        mapping = aes(x = x1, y = y1))  # Data and aesthetic mapping 
plot1 <- plot1 + geom_point()       # Geometry
plot1

#' 
#' ## `ggplot2` -- Data
#' In `ggplot2`, the data is the first argument of the `ggplot` function. It is usually a data frame --or an object that can be converted to one--.
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## plot1 <- ggplot(data = anscombe,

#' 
#' ## `ggplot2` -- Aesthetic mapping
#' 
#' The aesthetic mapping establishes the relationships between data variables and graphical elements. It is the second argument of the  `ggplot`function and it must be created form the `aes` support function.
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
##         mapping = aes(x = x1, y = y1))

#' 
#' ----
#' 
#' The most adequate mappings will depend on the type of variable and the type of plot --geometry--.
#' 
#' For quantitative variables, the most usual mappings are positions --`x`, `y`...--, `size`, and colors --`color`, `fill`--.
#' 
#' For qualitative variables cualitativas, we opt for positions --`x`, `y`...--, colors --`color`, `fill`--, and `shape` as the mappings of preference.
#' 
#' 
#' ## `ggplot2` -- Geometry
#' 
#' The geometries (`geom_...`) indicate the format --type-- of the plot, 
#' meaning how the graphical variables are put together. One or more geometries (each will be a layer) can be included in a plot by adding them to the `gg` object. 
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## plot1 <- plot1 + geom_point()

#' 
#' ----
#' 
#' Common geometries are
#' 
#' - `geom_point`
#' - `geom_line`, `geom_vline`, `geom_hline`
#' - `geom_bar`
#' - `geom_histogram`
#' - `geom_boxplot`
#' 
#' For each geometry, it is important to check what are the required and available mappings. This is available in the help pages and in the cheatsheet available at <https://www.rstudio.com/resources/cheatsheets/>.
#' 
#' 
#' ## `ggplot2` Graphs
#' 
#' We will discuss here how to create the most common graphs in `ggplot2`. We will be adding specific considerations where required.
#' We will explore
#' 
#' - scatter plots,
#' - histograms,
#' - bar plots, and  
#' - boxplots.
#' 
#' By default, the plots in `ggplot2` are constructed from a data frame from which the variables to be represented are selected.
#' 
#' ----
#' 
#' We will use a sample of 1000 points from the `diamonds` data set --included in the `ggplot2` package-- to create the examples.
#' 
## ---- echo = TRUE----------------------------------------------------------
diaM <- diamonds[sample(1:nrow(diamonds),1000),]
str(diaM)

#' 
#' 
#' ## `ggplot2` Graphs -- Scatter Plot
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## ggplot(diaM, aes(x=carat,y=price)) + geom_point()

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE--------------------------------------------
ggplot(diaM, aes(x=carat,y=price)) + geom_point()

#' 
#' ----
#' 
#' Adding a third variable and adjusting formats...
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## ggplot(diaM, aes(x=carat,y=price,color=cut)) +
##   geom_point(alpha=.8,shape=21,size=3)

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE--------------------------------------------
ggplot(diaM, aes(x=carat,y=price,color=cut)) + 
  geom_point(alpha=.8,shape=21,size=3)

#' 
#' -----
#' 
#' With a trend line...
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## ggplot(diaM, aes(x=carat,y=price,color=cut)) +
##   geom_point(alpha=.8,shape=21,size=3) +
##   geom_smooth(method="lm",se=FALSE)

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE--------------------------------------------
ggplot(diaM, aes(x=carat,y=price,color=cut)) + 
  geom_point(alpha=.8,shape=21,size=3) +
  geom_smooth(method="lm",se=FALSE)

#' 
#' ## `ggplot2` Graphs -- Histogram
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## ggplot(diaM, aes(x=price)) + geom_histogram(binwidth=1000)

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE--------------------------------------------
ggplot(diaM, aes(x=price)) + geom_histogram(binwidth=1000)

#' 
#' ----
#' 
#' Adding the cut variable...
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## ggplot(diaM, aes(x=price,fill=cut)) +
##   geom_histogram(position='dodge',binwidth=1000)

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE--------------------------------------------
ggplot(diaM, aes(x=price,fill=cut)) +
  geom_histogram(position='dodge',binwidth=1000)

#' 
#' ----
#' 
#' In relative frequencies, per group...
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## ggplot(diaM, aes(x=price,y=..density..,fill=cut)) +
##   geom_histogram(position='dodge',binwidth=1000)

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE--------------------------------------------
ggplot(diaM, aes(x=price,y=..density..,fill=cut)) +
  geom_histogram(position='dodge',binwidth=1000)

#' 
#' ----
#' 
#' Let's try with a density plot...
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## ggplot(diaM, aes(x=price,fill=cut)) +
##   geom_density(alpha=.3)

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE--------------------------------------------
ggplot(diaM, aes(x=price,fill=cut)) +
  geom_density(alpha=.3)

#' 
#' 
#' ## `ggplot2` Graphs -- Bar Plot
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## ggplot(diaM, aes(x=clarity)) +
##   geom_bar()

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE--------------------------------------------
ggplot(diaM, aes(x=clarity)) +
  geom_bar()

#' 
#' ----
#' 
#' Adding the cut variable...
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## ggplot(diaM, aes(x=clarity, fill=cut)) +
##   geom_bar()

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE--------------------------------------------
ggplot(diaM, aes(x=clarity, fill=cut)) +
  geom_bar()

#' 
#' ----
#' 
#' To compare absolute frequencies, dodged bars work better...
#' 
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## ggplot(diaM, aes(x=clarity, fill=cut)) +
##   geom_bar(position="dodge")

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE--------------------------------------------
ggplot(diaM, aes(x=clarity, fill=cut)) +
  geom_bar(position="dodge")

#' 
#' ----
#' 
#' To compare cumulative relative frequencies, per-class stacked bars are more useful.
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## ggplot(diaM, aes(x=clarity, fill=cut)) +
##   geom_bar(position="fill")

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE--------------------------------------------
ggplot(diaM, aes(x=clarity, fill=cut)) +
  geom_bar(position="fill")

#' 
#' ----
#' 
#' Bar plots can also be constructed from frequency counts. In this case, the `y` variable has to be stated and the argument `stat="identity"` has to be added to the `geom_bar`.
#' 
#' 
#' ## `ggplot2` Graphs -- Boxplot
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## ggplot(diaM, aes(x=1, y=price)) +
##   geom_boxplot()

#' 
#' *NOTE*: For the `geom_boxplot`, both `x` and `y` are required. If there isn't an independent variable `x=1` should be included.
#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE--------------------------------------------
ggplot(diaM, aes(x=1, y=price)) +
  geom_boxplot()

#' 
#' ----
#' 
#' Adding the cut variable...
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## ggplot(diaM, aes(x=cut, y=price)) +
##   geom_boxplot()

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE--------------------------------------------
ggplot(diaM, aes(x=cut, y=price)) +
  geom_boxplot()

#' 
#' ----
#' 
#' The chart can be improved by showing all the data points with the position modified with some random noise --jittered--.
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## ggplot(diaM, aes(x=cut, y=price)) +
##   geom_boxplot(outlier.shape = NA) +
##   geom_jitter(shape = 21, alpha=.5,height=0,width=.2)

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE--------------------------------------------
ggplot(diaM, aes(x=cut, y=price)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 21, alpha=.5,height=0,width=.2)

#' 
#' ----
#' 
#' Or including a violin plot, `geom_violin`, and an addtional point with the the mean.
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## medias <- diaM %>% group_by(cut) %>%
##   summarise(price=mean(price))
## 
## ggplot(diaM, aes(x=cut, y=price)) +
##   geom_violin() +
##   geom_boxplot(outlier.shape = NA, width = 0.1) +
##   geom_point(data=medias,shape=3)

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE--------------------------------------------
medias <- diaM %>% group_by(cut) %>%
  summarise(price=mean(price))

ggplot(diaM, aes(x=cut, y=price)) + 
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.1) +
  geom_point(data=medias,shape=3)

#' 
#' 
#' ## `ggplot2` -- Additional Elements
#' 
#' Beyond the already presented elements --data, aesthetic mapping, and geometries--, other `ggplot2` elements help adjust many other aspects of the plots, for example
#' 
#' - `scale_...` controls aspects related to the presentation of the graphical variables,
#' - `coord_...` adjustes the coordinate system used in the plot,
#' - `labs` sets plot and axes titles,
#' - `theme`, `theme_...` control the non-data ink elements of the plot, and
#' - `facet_` facilitates the creation of multiple plots.
#' 
#' ----
#' 
#' One last example...
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## ggplot(diaM, aes(x=carat, y = price,
##                  shape = cut, col = clarity)) +
##   geom_point(alpha=.6)
## 

#' 
#' ----
#' 
## ---- echo = FALSE---------------------------------------------------------
ggplot(diaM, aes(x=carat, y = price,
                 shape = cut, col = clarity)) +
  geom_point(alpha=.6)


#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE--------------------------------------------
## ggplot(diaM, aes(x=carat, y = price,
##                  shape = cut, col = clarity)) +
##   geom_point(alpha=.6) +
##   scale_x_continuous(breaks=1:3) +
##   scale_y_continuous(trans="log10") +
##   scale_color_brewer(palette="Spectral")+
##   facet_grid(cut ~ clarity) +
##   theme_bw() +
##   theme(legend.position = "none",
##         text = element_text(size=10))
## 

#' 
#' ----
#' 
## ---- echo = FALSE---------------------------------------------------------
ggplot(diaM, aes(x=carat, y = price,
                 shape = cut, col = clarity)) +
  geom_point(alpha=.8) +
  scale_x_continuous(breaks=1:3) +
  scale_y_continuous(trans="log10") +
  scale_color_brewer(palette="Spectral")+
  facet_grid(cut ~ clarity) + 
  theme_bw() + 
  theme(legend.position = "none",
        text = element_text(size=10))


#' 
#' ## YOUR TURN {data-background="#eeffcc"}
#' 
#' 1.    Revisit the `OrchardSprays` data set and make two plots that highlight the conclusions that can be reached from analyzing it.
#' 

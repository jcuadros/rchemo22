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
#' editor_options: 
#'   chunk_output_type: console
#' ---
#' 
## ----setup, include=FALSE----------------------------------------------------
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
#' 
#' # Introduction and References
#' 
#' ## Chemometrics
#' 
#' "Chemometrics is the chemical discipline that uses mathematical, statistical, and other methods employing formal logic to design or select optimal measurement procedures and experiments, and to provide maximum relevant chemical information by analyzing chemical data."
#' 
#' ::: {.bibref}
#' Vékey, K., Telekes, A., & Vertes, A. (Eds.). (2011). Medical applications of mass spectrometry. Elsevier.
#' :::
#' 
#' ## Main Chemometrics Applications
#' 
#' -   Design of experiments (before data acquisition, mind works better than software here)
#' -   Chemical data preprocessing
#' -   Multivariate exploratory analysis
#' -   Classification
#' -   Modelling
#' -   Multiple Curve Resolution
#'     
#' 
#' ## R packages for Chemistry and Chemometrics {.small}
#' 
#' - `chemometrics`, https://cran.r-project.org/web/packages/chemometrics/index.html
#' - `ChemometricsWithR`, https://github.com/rwehrens/ChemometricsWithR/
#' - `ChemoSpec`, https://cloud.r-project.org/web/packages/ChemoSpec/index.html
#' - `mdatools`, https://cran.r-project.org/web/packages/mdatools/index.html
#' 
#' Updated information can be found following the CRAN Task View: Chemometrics and Computational Physics, (https://cran.r-project.org/web/views/ChemPhys.html).
#' 
#' There are many other packages that will useful for specific tasks. These will appear along the session.
#' 
#' 
#' ## Chemometrics with R -- References {.small}
#' 
#' -   Harvey, D. (2021). Chemometrics Using R (Harvey). Libretexts. <https://chem.libretexts.org/Bookshelves/Analytical_Chemistry/Chemometrics_Using_R_(Harvey)>
#' -   Kucheryavskiy, s. (2021). mda.tools: make chemometrics easy. <https://mdatools.com/>
#' -   Stanstrup, J., Broeckling, C. D., Helmus, R., Hoffmann, N., Mathé, E., Naake, T., ... & Neumann, S. (2019). The metaRbolomics Toolbox in Bioconductor and beyond. Metabolites, 9(10), 200. <https://www.mdpi.com/2218-1989/9/10/200/htm>
#' -   Varmuza, K., & Filzmoser, P. (2016). Introduction to multivariate statistical analysis in chemometrics. CRC press.
#' -   Wehrens, R. (2011). Chemometrics with R. Heidelberg, Germany: Springer.
#' 
#' 
#' # Chemical Data Preprocessing
#' 
#' ## Before Starting, some Package Management
#' 
## ---- echo = TRUE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
pckgs<-c("tidyverse","ggthemes","GGally",
         "ggalt","gtools","gridExtra")
pckgs2Install<-pckgs[!(pckgs %in% library()$results[,1])]
pckgs2Load<-pckgs[!(pckgs %in% (.packages()))]
for(pckg in pckgs2Install) {
  install.packages(pckg,
    repos="https://cloud.r-project.org/",
    quiet=TRUE, type="binary")
}
for(pckg in pckgs2Load) {
  library(pckg, character.only = TRUE)
}

#' 
#' ## From Raw Data to Processed Data {.small}
#' 
#' Many times raw data needs to be processed to be analyzed and compared to preexisting data.
#' 
#' Common processing includes
#' 
#' -   Noise removal
#' -   Baseline adjustment
#' -   Peak alignment
#' -   Binning
#' -   Peak picking
#' -   Scaling
#' -   Managing missing data
#' -   Removing outliers
#' -   Comparing complex data
#' 
#' In the case of chemical data, this is specially true is for spectral data.
#' 
#' ----
#' 
#' ### Spectral Data
#' 
#' There are different types of spectral data
#' 
#' -   Raw/Profile Data
#' -   Clean Data (with some adjustments)
#' -   Binned Data
#' -   Centroided Data
#' -   Peak Lists
#'     -   With Intensity
#'     -   Without Intensity
#' 
#' ----
#' 
#' ### Data Sets
#' 
#' Let's explore this process with several data sets:
#' 
#' -   a gas chromatogram from `ptw::gaschrom`,
#' -   a NIR spectra from `pls::gasoline`,
#' -   a set of NIR of wood samples from the RNIR web book, <https://guifh.github.io/RNIR/index.html>, 
#' -   an extract from a LC-MS from the `msdata` package, and
#' -   a MALDI-TOF spectra from a milk mixture available in the `baseline` package.
#' 
#' ----
#' 
## ---- echo = TRUE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----

if(!require("ptw")) {
  install.packages("ptw", repos="https://cloud.r-project.org/")
  library("ptw")
}
data("gaschrom")


#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## help("gaschrom", package="ptw")
## 

#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## gchrom <- data.frame(rt=1:5000,
##                      I=gaschrom[1,])
## plot(gchrom, type="l")
## 

#' 
#' ----
#' 
## ---- echo = FALSE-----------------------------------------------------------
gchrom <- data.frame(rt=1:5000,
                     I=gaschrom[1,])
plot(gchrom, type="l")


#' 
#' ----
#' 
## ---- echo = TRUE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("pls")) {
  install.packages("pls", repos="https://cloud.r-project.org/")
  library("pls")
}
data(gasoline)

#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## help("gasoline", package="pls")
## 

#' 
#' ----
#' 
## ---- echo = TRUE,eval = FALSE-----------------------------------------------
## gasolineNIR <- data.frame(wl=seq(900,1700,by=2),
##                           logRef=gasoline$NIR[4,])
## plot(gasolineNIR, type="l")
## 

#' 
#' ----
#' 
## ---- echo = FALSE-----------------------------------------------------------
gasolineNIR <- data.frame(wl=seq(900,1700,by=2),
                          logRef=gasoline$NIR[4,])
plot(gasolineNIR, type="l")


#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## woodNIR <- get(load("../data/mydataPHAZIR.Rdata"))
## 
## spectra <- as.matrix(woodNIR$NIR)
## dfSpectrum <- data.frame(wl = colnames(spectra),
##                          Abs = as.numeric(spectra[1,])) %>%
##   mutate(wl=as.numeric(substr(wl,2,nchar(wl))))
## 
## plot(dfSpectrum, type="l")

#' 
#' ----
#' 
## ---- echo = FALSE-----------------------------------------------------------
woodNIR <- get(load("../data/mydataPHAZIR.Rdata"))

spectra <- as.matrix(woodNIR$NIR)
dfSpectrum <- data.frame(wl = colnames(spectra),
                         Abs = as.numeric(spectra[1,])) %>% 
  mutate(wl=as.numeric(substr(wl,2,nchar(wl))))

plot(dfSpectrum, type="l")

#' 
#' 
#' ----
#' 
## ---- echo = TRUE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager",
                   repos="https://cloud.r-project.org/")
if (!require("msdata", quietly = TRUE)) {
  BiocManager::install("msdata", update=FALSE)
  library("msdata")
}

dir(system.file("sciex", package = "msdata"))


#' 
#' ----
#' 
## ---- echo = TRUE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if (!require("mzR", quietly = TRUE)) {
  BiocManager::install("mzR", update=FALSE)
  library("mzR")
}

sciexMS <- mzR::openMSfile(
  system.file("sciex/20171016_POOL_POS_1_105-134.mzML",
              package = "msdata"))

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
sciexSpec <- mzR::spectra(sciexMS)
length(sciexSpec)
dim(sciexSpec[[1]])

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
plot(sciexSpec[[1]], type="l")


#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
plot(sciexSpec[[931]], type="l")


#' 
#' ---
#' 
## ---- echo = TRUE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("baseline")) {
  install.packages("baseline",
                   repos="https://cloud.r-project.org/")
  library("baseline")
}

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
data("milk")
str(milk)

#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## help("milk")
## 
## milkMS <- data.frame(mz = as.numeric(colnames(milk$spectra)),
##                      I = milk$spectra[20,])
## plot(milkMS, type="l")

#' 
#' ----
#' 
## ---- echo = FALSE-----------------------------------------------------------
milkMS <- data.frame(mz = as.numeric(colnames(milk$spectra)),
                     I = milk$spectra[20,])
plot(milkMS, type="l")


#' 
#' ## Noise Removal
#' 
#' Main methods for noise removal are:
#' 
#' -   Moving averages,
#' -   Savitzky-Golay filters, <https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter>, and
#' -   FFT filters.
#' 
#' Running medians, threshold, aggregating repetitions and consensus approaches are also used.
#' 
#' 
#' ## Noise Removal -- Moving Average
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## milkMS_mave <- milkMS
## milkMS_mave$I <- as.numeric(stats::filter(milkMS_mave$I,
##               rep(1,13)/13))
## 
## g1 <- ggplot(milkMS, aes(x=mz, y=I)) +
##   geom_line() +
##   theme_classic()
## g2 <- ggplot(milkMS_mave, aes(x=mz, y=I)) +
##   geom_line() +
##   theme_classic()
## grid.arrange(g1,g2)
## 

#' 
#' ----
#' 
## ---- echo = FALSE, warning = FALSE------------------------------------------
milkMS_mave <- milkMS
milkMS_mave$I <- as.numeric(stats::filter(milkMS_mave$I, 
              rep(1,13)/13))

g1 <- ggplot(milkMS, aes(x=mz, y=I)) + 
  geom_line() +
  theme_classic()
g2 <- ggplot(milkMS_mave, aes(x=mz, y=I)) + 
  geom_line() +
  theme_classic()
grid.arrange(g1,g2)


#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## g1 <- ggplot(milkMS, aes(x=mz, y=I)) +
##   geom_line() +
##   coord_cartesian(xlim=c(10000,11000),ylim=c(0,500)) +
##   theme_classic()
## g2 <- ggplot(milkMS_mave, aes(x=mz, y=I)) +
##   geom_line() +
##   coord_cartesian(xlim=c(10000,11000),ylim=c(0,500)) +
##   theme_classic()
## grid.arrange(g1,g2)
## 

#' 
#' ----
#' 
## ---- echo = FALSE, warning = FALSE------------------------------------------
g1 <- ggplot(milkMS, aes(x=mz, y=I)) + 
  geom_line() +
  coord_cartesian(xlim=c(10000,11000),ylim=c(0,500)) +
  theme_classic()
g2 <- ggplot(milkMS_mave, aes(x=mz, y=I)) + 
  geom_line() +
  coord_cartesian(xlim=c(10000,11000),ylim=c(0,500)) +
  theme_classic()
grid.arrange(g1,g2)


#' 
#' ## Noise Removal -- Savitzky-Golay
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## milkMS_sg7 <- milkMS
## milkMS_sg7$I <- as.numeric(stats::filter(milkMS_mave$I,
##               c(-2,3,6,7,6,3,-2)/21))
## 
## g1 <- ggplot(milkMS, aes(x=mz, y=I)) +
##   geom_line() +
##   coord_cartesian(xlim=c(10000,11000),ylim=c(0,500)) +
##   theme_classic()
## g2 <- ggplot(milkMS_sg7, aes(x=mz, y=I)) +
##   geom_line() +
##   coord_cartesian(xlim=c(10000,11000),ylim=c(0,500)) +
##   theme_classic()
## grid.arrange(g1,g2)
## 

#' 
#' ----
#' 
## ---- echo = FALSE, warning = FALSE------------------------------------------
milkMS_sg7 <- milkMS
milkMS_sg7$I <- as.numeric(stats::filter(milkMS_mave$I, 
              c(-2,3,6,7,6,3,-2)/21))

g1 <- ggplot(milkMS, aes(x=mz, y=I)) + 
  geom_line() +
  coord_cartesian(xlim=c(10000,11000),ylim=c(0,500)) +
  theme_classic()
g2 <- ggplot(milkMS_sg7, aes(x=mz, y=I)) + 
  geom_line() +
  coord_cartesian(xlim=c(10000,11000),ylim=c(0,500)) +
  theme_classic()
grid.arrange(g1,g2)


#' 
#' ## Noise Removal -- FFT
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## I_fft <- fft(milkMS$I)
## I_fft[513:(length(I_fft)-512)] <- 0
## milkMS_fft <- data.frame(mz = milkMS$mz,
##         I=Re(fft(I_fft, inverse=TRUE)/length(I_fft)))
## 
## g1 <- ggplot(milkMS, aes(x=mz, y=I)) +
##   geom_line() +
##   theme_classic()
## g2 <- ggplot(milkMS_fft, aes(x=mz, y=I)) +
##   geom_line() +
##   theme_classic()
## grid.arrange(g1,g2)
## 

#' 
#' ----
#' 
## ---- echo = FALSE, warning = FALSE------------------------------------------
I_fft <- fft(milkMS$I)
I_fft[513:(length(I_fft)-512)] <- 0
milkMS_fft <- data.frame(mz = milkMS$mz,
        I=Re(fft(I_fft, inverse=TRUE)/length(I_fft)))

g1 <- ggplot(milkMS, aes(x=mz, y=I)) + 
  geom_line() +
  theme_classic()
g2 <- ggplot(milkMS_fft, aes(x=mz, y=I)) + 
  geom_line() +
  theme_classic()
grid.arrange(g1,g2)


#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## g1 <- ggplot(milkMS, aes(x=mz, y=I)) +
##   geom_line() +
##   coord_cartesian(xlim=c(10000,15000),ylim=c(0,4000)) +
##   theme_classic()
## g2 <- ggplot(milkMS_fft, aes(x=mz, y=I)) +
##   geom_line() +
##   coord_cartesian(xlim=c(10000,15000),ylim=c(0,4000)) +
##   theme_classic()
## grid.arrange(g1,g2)
## 

#' 
#' ----
#' 
## ---- echo = FALSE, warning = FALSE------------------------------------------
g1 <- ggplot(milkMS, aes(x=mz, y=I)) + 
  geom_line() +
  coord_cartesian(xlim=c(10000,15000),ylim=c(0,4000)) +
  theme_classic()
g2 <- ggplot(milkMS_fft, aes(x=mz, y=I)) + 
  geom_line() +
  coord_cartesian(xlim=c(10000,15000),ylim=c(0,4000)) +
  theme_classic()
grid.arrange(g1,g2)


#' 
#' 
#' ## Baseline Adjustment
#' 
#' Baseline adjustment is commonly done by finding local minima and fitting a curve to them. This curve is then substracted to the spectrum.
#' 
#' Common methods for baseline adjustment are:
#' 
#' -   Multiplicative Scatter Correction (MSC), for NIR, and
#' -   Asymmetric Least Squares (ALS).
#' 
#' MSC is available in the `pls` package. ALS and many other methods, such as LOWESS or a polynomial fit, are implemented in the `baseline` package.
#' 
#' 
#' ## Baseline Adjustment - MSC {.small}
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## originalspectra <- as.matrix(woodNIR$NIR)
## newspectra<-msc(originalspectra)
## mydataOrig <- data.frame(id=woodNIR[,5], NIR = originalspectra) %>%
##   pivot_longer(cols=-1, names_to = "wl", values_to = "Abs") %>%
##   mutate(wl=as.numeric(substr(wl,6,nchar(wl))))
## mydataMSC <- data.frame(id=woodNIR[,5], NIR = newspectra) %>%
##   pivot_longer(cols=-1, names_to = "wl", values_to = "Abs") %>%
##   mutate(wl=as.numeric(substr(wl,6,nchar(wl))))
## 
## g1 <- ggplot(mydataOrig, aes(x=wl, y=Abs, color=as.factor(id))) +
##   geom_line() +
##   theme_classic() + theme(legend.position = "none")
## g2 <- ggplot(mydataMSC, aes(x=wl, y=Abs, color=as.factor(id))) +
##   geom_line() +
##   theme_classic() + theme(legend.position = "none")
## grid.arrange(g1,g2)

#' 
#' ----
#' 
## ---- echo = FALSE-----------------------------------------------------------
originalspectra <- as.matrix(woodNIR$NIR)
newspectra<-msc(originalspectra)            
mydataOrig <- data.frame(id=woodNIR[,5], NIR = originalspectra) %>% 
  pivot_longer(cols=-1, names_to = "wl", values_to = "Abs") %>% 
  mutate(wl=as.numeric(substr(wl,6,nchar(wl))))
mydataMSC <- data.frame(id=woodNIR[,5], NIR = newspectra) %>% 
  pivot_longer(cols=-1, names_to = "wl", values_to = "Abs") %>% 
  mutate(wl=as.numeric(substr(wl,6,nchar(wl))))             

g1 <- ggplot(mydataOrig, aes(x=wl, y=Abs, color=as.factor(id))) + 
  geom_line() +
  theme_classic() + theme(legend.position = "none") 
g2 <- ggplot(mydataMSC, aes(x=wl, y=Abs, color=as.factor(id))) + 
  geom_line() +
  theme_classic() + theme(legend.position = "none") 
grid.arrange(g1,g2)

#' 
#' ## Baseline Adjustment - ALS
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## milkMS_als <- milkMS
## sp <- baseline.als(t(milkMS_als$I), lambda=10)
## milkMS_als$I <- sp$corrected[1,]
## 
## g1 <- ggplot(milkMS, aes(x=mz, y=I)) +
##   geom_line() +
##   theme_classic()
## g2 <- ggplot(milkMS_als, aes(x=mz, y=I)) +
##   geom_line() +
##   theme_classic()
## grid.arrange(g1,g2)
## 

#' 
#' ----
#' 
## ---- echo = FALSE, warning = FALSE------------------------------------------
milkMS_als <- milkMS
sp <- baseline.als(t(milkMS_als$I), lambda=10)
milkMS_als$I <- sp$corrected[1,]

g1 <- ggplot(milkMS, aes(x=mz, y=I)) + 
  geom_line() +
  theme_classic()
g2 <- ggplot(milkMS_als, aes(x=mz, y=I)) + 
  geom_line() +
  theme_classic()
grid.arrange(g1,g2)


#' 
#' 
#' ## Peak Alignment
#' 
#' Different phenomena can provoke misalignments among spectra. In order to compare them peaks must be aligned. In many cases, this is experimentally done with standards or solved through binning; but in some cases, data processing is required, for example in LC chromatograms.
#' 
#' Shifts can be corrected using time warping functions, for example the `ptw` or the `dtw` functions in the package with the same names.
#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## gchrom1 <- data.frame(rt=1:5000,
##                      I=gaschrom[1,])
## gchrom2 <- data.frame(rt=1:5000,
##                      I=gaschrom[16,])
## 
## ggplot(gchrom1,aes(x=rt, y=I)) +
##   geom_line(aes(y=I+50), data=gchrom2, color="red") +
##   geom_line() +
##   theme_classic()
## 

#' 
#' ----
#' 
## ---- echo = FALSE-----------------------------------------------------------
gchrom1 <- data.frame(rt=1:5000,
                     I=gaschrom[1,])
gchrom2 <- data.frame(rt=1:5000,
                     I=gaschrom[16,])

ggplot(gchrom1,aes(x=rt, y=I)) +
  geom_line(aes(y=I+50), data=gchrom2, color="red") +
  geom_line() +
  theme_classic()


#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
wcc(gchrom1[,2], gchrom2[,2], length(gchrom1[,2]))
res <- ptw(gchrom1[,2], gchrom2[,2], init.coef=c(0,1,0))
summary(res)    # WCC is a distance

#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## gchrom2c <- data.frame(rt=gchrom1$rt,
##                        I=res$warped.sample[1,])
## ggplot(gchrom1,aes(x=rt, y=I)) +
##   geom_line(aes(y=I+50), data=gchrom2c, color="green") +
##   geom_line() +
##   theme_classic()

#' 
#' ----
#' 
## ---- echo = FALSE, warning = FALSE------------------------------------------
gchrom2c <- data.frame(rt=gchrom1$rt,
                       I=res$warped.sample[1,])
ggplot(gchrom1,aes(x=rt, y=I)) +
  geom_line(aes(y=I+50), data=gchrom2c, color="green") +
  geom_line() +
  theme_classic()

#' 
#' 
#' ## Binning
#' 
#' Especially in MS, binning is a required process to allow for spectra comparison.
#' 
#' ## {.small}
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## dfMS1 <- data.frame(mz = sciexSpec[[1]][,1],
##                     I = sciexSpec[[1]][,2])
## dfMS2 <- dfMS1 %>% mutate(mz=round(2*mz)/2) %>%
##   group_by(mz) %>% summarize(I=sum(I),.groups="drop")
## dfMS3 <- dfMS1 %>% mutate(mz=round(mz/2)*2) %>%
##   group_by(mz) %>% summarize(I=sum(I),.groups="drop")
## 
## g1 <- ggplot(dfMS1, aes(x=mz, xend=mz,y=0, yend=I))+
##   geom_segment() + theme_classic()
## g2 <- ggplot(dfMS2, aes(x=mz, xend=mz,y=0, yend=I))+
##   geom_segment() + theme_classic()
## g3 <- ggplot(dfMS3, aes(x=mz, xend=mz,y=0, yend=I))+
##   geom_segment() + theme_classic()
## 
## grid.arrange(g1,g2,g3)
## 

#' 
#' ----
#' 
## ---- echo = FALSE-----------------------------------------------------------
dfMS1 <- data.frame(mz = sciexSpec[[1]][,1],
                    I = sciexSpec[[1]][,2])
dfMS2 <- dfMS1 %>% mutate(mz=round(2*mz)/2) %>% 
  group_by(mz) %>% summarize(I=sum(I),.groups="drop")
dfMS3 <- dfMS1 %>% mutate(mz=round(mz/2)*2) %>% 
  group_by(mz) %>% summarize(I=sum(I),.groups="drop")

g1 <- ggplot(dfMS1, aes(x=mz, xend=mz,y=0, yend=I))+
  geom_segment() + theme_classic()
g2 <- ggplot(dfMS2, aes(x=mz, xend=mz,y=0, yend=I))+
  geom_segment() + theme_classic()
g3 <- ggplot(dfMS3, aes(x=mz, xend=mz,y=0, yend=I))+
  geom_segment() + theme_classic()

grid.arrange(g1,g2,g3)


#' 
#' ## Peak Picking
#' 
#' It is also common in spectra processing to convert the profile data to a peak list.
#' 
#' This is usually done looking for local maxima and/or applying some threshold consideration.
#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
str(dfMS2)
dfMS2[dfMS2$I/max(dfMS2$I)>.05,]

span <- 100
dfMS1_g <- embed(dfMS1$I,span)
dfMS1_gmax <- span + 1 - apply(dfMS1_g,1,which.max)
dfMS1_gmax[dfMS1_gmax==1 | dfMS1_gmax==span] <- NA
dfMS1_peaksPos <- dfMS1_gmax + 0:(length(dfMS1_gmax)-1)
dfMS1_peaksPos <- unique(na.omit(dfMS1_peaksPos))
dfMS1$mz[dfMS1_peaksPos]

#' 
#' ## Scaling
#' 
#' In many cases, scaling can change the outcomes of an analysis, so it is important to decide what scaling makes sense. Some examples can be
#' 
#' -   Scaling towards an internal reference,
#' -   Scaling towards the maximum signal,
#' -   Scaling towards the total area, or
#' -   Not scaling at all.
#' 
#' 
#' ## Outliers and Missing Data Imputation
#' 
#' Although analysis of leverage points, dtection of outliers and imputation of missing data are always relevant, we will not discuss it further here. Just keep in mind there is always a necessary step in any data analysis.
#' 
#' 
#' ## Data Matching
#' 
#' All this process is many times done to compare chemical entities: spectra, molecules, mixtures... In any of this cases, decisions will need to be made in order to decide when to data match and/or how distance must be calculated between two data instances.
#' 
#' Some examples of strategies for matching data are
#' 
#' -   Fingerprints
#' -   Layered notations, like InChI,
#' -   Hashed codes, like InChIKey or SPLASH
#' -   Correlation coefficients, 
#' -   Distances on the principal components space.
#' 
#' 
#' # Multivariant Exploratory Analysis
#' 
#' ## Introduction and Data
#' 
#' Whether it comes from spectral data or not, it is common to have multivariant data in chemical analysis. In fact, chemometrics is commonly restricted to the studya and analysis of multivariant data.
#' 
#' Common multivariant techniques in chemometrics are
#' 
#' -   Exploratory anlysis
#' -   Classification
#' -   Modelling
#' -   Multiple Curve Resolution
#' 
#' We will discuss here some aspects of multivariant exploratory analysis.
#' 
#' ----
#' 
#' We'll use the classical Forina's `wines` dataset, which is included in the `kohonen` package.
#' 
## ---- echo = TRUE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("kohonen")) {
  install.packages("kohonen",
                   repos="https://cloud.r-project.org/")
  library("kohonen")
}

data("wines")

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
dfWines <- data.frame(wines)
str(dfWines)

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
summary(dfWines[,1:5])

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
summary(dfWines[,6:10])

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
summary(dfWines[,11:13])

#' 
#' 
#' ## Correlations
#' 
## ---- echo = TRUE------------------------------------------------------------
corWines <- data.frame(cor(dfWines)) %>% 
  rownames_to_column("var1") %>% 
  pivot_longer(2:14, names_to="var2", values_to = "cor") %>% 
  filter(var1<var2) %>% 
  arrange(desc(abs(cor)))
str(corWines)

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
head(corWines)

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
tail(corWines)

#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## if(!require("corrplot")) {
##   install.packages("corrplot",
##                    repos="https://cloud.r-project.org/")
##   library("corrplot")
## }
## 
## corrplot(cor(wines))

#' 
#' ----
#' 
## ---- echo = FALSE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("corrplot")) {
  install.packages("corrplot",
                   repos="https://cloud.r-project.org/")
  library("corrplot")
}
corrplot(cor(wines))

#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## corrplot(cor(wines), order="hclust")

#' 
#' ----
#' 
## ---- echo = FALSE-----------------------------------------------------------
corrplot(cor(wines), order="hclust")

#' 
#' ----
#' 
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## ggpairs(dfWines, axisLabels = "none",
##         lower = list(continuous = wrap("points",
##               alpha = 0.3, shape = 21, size = 1))) +
##   theme_classic()

#' 
#' ----
#' 
## ---- echo = FALSE, message = FALSE------------------------------------------
ggpairs(dfWines, axisLabels = "none",
        lower = list(continuous = wrap("points",
              alpha = 0.3, shape = 21, size = 1))) + 
  theme_classic()

#' 
#' 
#' ## Factor Analysis
#' 
#' It consists in searching for a limited number of factor that can conserve most of the data variability.
#' 
#' The most common way to perform this analysis is
#' 
#' -   Principal Components Analysis (PCA), followed by
#' -   Rotation
#' 
#' 
#' ## Factor Analysis -- PCA Assumptions
#' 
## ---- echo = TRUE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("REdaS")) {
  install.packages("REdaS",
                   repos="https://cloud.r-project.org/")
  library("REdaS")
}

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
REdaS::KMOS(dfWines)  # Kaiser-Meyer-Olkin Sample Adequacy

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
REdaS::bart_spher(cor(dfWines))  # Bartlett's Sphericty Test

#' 
#' 
#' ## Factor Analysis -- PCA
#' 
## ---- echo = TRUE------------------------------------------------------------
pcaWines <- princomp(dfWines, cor=TRUE)
summary(pcaWines)

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
str(pcaWines)

#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## biplot(pcaWines)

#' 
#' ----
#' 
## ---- echo = FALSE-----------------------------------------------------------
biplot(pcaWines)

#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## biplot(pcaWines, choices=c(1,3))

#' 
#' ----
#' 
## ---- echo = FALSE-----------------------------------------------------------
biplot(pcaWines, choices=c(1,3))

#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## screeplot(pcaWines)

#' 
#' ----
#' 
## ---- echo = FALSE-----------------------------------------------------------
screeplot(pcaWines)

#' 
#' 
#' ## Factor Analysis -- Rotation
#' 
## ---- echo = TRUE------------------------------------------------------------
rfa <- varimax(pcaWines$loadings[,1:3], eps=1e-14)
rfa

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
tfa <- promax(pcaWines$loadings[,1:3], m=4)
tfa

#' 
#' ----
#' 
## ---- echo = TRUE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("GPArotation")) {
  install.packages("GPArotation",
                   repos="https://cloud.r-project.org/")
  library("GPArotation")
}

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
obfa <- GPArotation::oblimin(pcaWines$loadings[,1:3])
obfa

#' 
#' 
#' ## Clustering
#' 
#' The most common clustering strategies are
#' 
#' -   Hierarchical clustering, especially for small data sets, and
#' -   K-means, for large data sets but consider selecting the features before clustering.
#' 
#' 
#' ## Clustering - Hierarchical
#' 
## ---- echo = TRUE------------------------------------------------------------
dfWinesScaled <- data.frame(scale(dfWines))
clWines_hca <- hclust(stats::dist(dfWinesScaled,
                                  method="euclidean"),
                      method="ward.D2")
str(clWines_hca)

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
winesGroups_hca <- cutree(clWines_hca, k=3)
table(winesGroups_hca)

#' 
#' ----
#' 
## ---- echo = TRUE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("ggdendro")) {
  install.packages("ggdendro",
                   repos="https://cloud.r-project.org/")
  library("ggdendro")
}

#' 
#' ----
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## ggdendrogram(clWines_hca, labels=FALSE) +
##   geom_hline(yintercept=9.2, color="red") +
##   theme(axis.text.x = element_blank())

#' 
#' ----
#' 
## ---- echo = FALSE-----------------------------------------------------------
ggdendrogram(clWines_hca, labels=FALSE) +
  geom_hline(yintercept=9.2, color="red") +
  theme(axis.text.x = element_blank())

#' 
#' 
#' ## Clustering - K-means
#' 
## ---- echo = TRUE------------------------------------------------------------
dfWinesScaled <- data.frame(scale(dfWines))
clWines_km <- kmeans(dfWinesScaled, centers=3)
str(clWines_km)

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
winesGroups_km <- clWines_km$cluster
table(winesGroups_km)


#' 
#' ----
#' 
## ---- echo = TRUE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("mclust")) {
  install.packages("mclust",
                   repos="https://cloud.r-project.org/")
  library("mclust")
}

#' 
#' ----
#' 
## ---- echo = TRUE------------------------------------------------------------
adjustedRandIndex(winesGroups_hca,winesGroups_km)

#' 
#' ## {.small}
#' 
## ---- echo = TRUE, eval = FALSE----------------------------------------------
## dfWinesClust <- cbind(dfWines,
##                       hca=factor(winesGroups_hca),
##                       km=factor(winesGroups_km))
## g1 <- ggplot(dfWinesClust,
##     aes(x=flavonoids, y=alcohol, color=hca)) +
##   geom_point(shape=3) + geom_encircle() + theme_classic() +
##   theme(legend.position = "top")
## g2 <- ggplot(dfWinesClust,
##     aes(x=flavonoids, y=alcohol, color=km)) +
##   geom_point(shape=3) + geom_encircle() + theme_classic() +
##   theme(legend.position = "top")
## 
## grid.arrange(g1,g2,ncol=2)
## 

#' 
#' ----
#' 
## ---- echo = FALSE-----------------------------------------------------------
dfWinesClust <- cbind(dfWines,
                      hca=factor(winesGroups_hca),
                      km=factor(winesGroups_km))
g1 <- ggplot(dfWinesClust,
    aes(x=flavonoids, y=alcohol, color=hca)) + 
  geom_point(shape=3) + geom_encircle() + theme_classic() +
  theme(legend.position = "top")
g2 <- ggplot(dfWinesClust,
    aes(x=flavonoids, y=alcohol, color=km)) + 
  geom_point(shape=3) + geom_encircle() + theme_classic() +
  theme(legend.position = "top")

grid.arrange(g1,g2,ncol=2)
  

#' 
#' ## Self-Organizing Maps
#' 
#' ## YOUR TURN {data-background="#eeffcc"}
#' 
#' 1.    Explore the data from mango aroma components available at <https://figshare.com/articles/dataset/Mango_Mangifera_indica_Aroma_Discriminate_Cultivars_and_Ripeness_Stages/14303913>. Warning: it needs some cleaning!
#' 
#' 
#' # Classification
#' 
#' ## MORE TO COME
#' 
#' 
#' 
#' ## Discriminant Analysis
#' 
#' ## kNN
#' 
#' ## Classification -- Other ML Approaches
#' 
#' ### Classification Trees & Forests
#' 
#' ### SVN
#' 
#' ### Neural Networks 
#' 
#' 
#' # Modelling
#' 
#' ## Calibration Models
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
#' ## Structural Equation Models
#' 
#' 
#' # Multiple Curve Resolution
#' 

############
#
#  Script to compile all Rmd files to HTML and R
#  @author: Jordi Cuadros
#
############

if(!require("knitr")) {
  install.packages("knitr")
  library("knitr")
}

if(!require("rmarkdown")) {
  install.packages("rmarkdown")
  library("rmarkdown")
}

if(!require("revealjs")) {
  install.packages("revealjs")
  library("revealjs")
}

RmdFiles <- dir(pattern="*.Rmd", recursive=TRUE)
for(i in seq(RmdFiles)) {
  render(RmdFiles[i],encoding ="UTF-8")
  file.remove(gsub(".Rmd",".R",RmdFiles[i],fixed=TRUE))
  purl(RmdFiles[i],gsub(".Rmd",".R",RmdFiles[i],fixed=TRUE),
       encoding ="UTF-8", documentation = 2)
}


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
#' ---
#' 
## ----setup, include=FALSE-------------------------------------------------------
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
#' ## Instalación y carga de paquetes en R
#' R tiene muchos paquetes para resolver problemas específicos. Para usar un paquete adicional este debe instalarse y cargarse en memoria.
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## installed.packages()[,1] # Lista los paquetes instalados
## (.packages()) # Lista los paquetes en memoria

#' 
#' Dos recursos importantes para buscar y identificar paquetes relevantes son
#' https://www.rdocumentation.org/
#' https://cran.rstudio.com/web/views/
#' 
#' ----
#' 
#' Para instalar un paquete, por ejemplo "nycflights13"
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## install.packages("nycflights13")

#' Para cargar un paquete en memoria, se usa
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## library("nycflights13")

#' 
#' En RStudio la gestión de paquetes también puede hacerse desde la interfaz del programa
#' 
#' ----
#' 
#' En un script y para garantizar la disponibilidad de un paquete
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## if(!require("nycflights13")) {
##   install.packages("nycflights13")
##   library("nycflights13")
## }

#' 
#' ----
#' 
#' Para acceder a la documentación de un paquete
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## help(package="nycflights13")

#' 
#' Algunos paquetes tienen información adicional a la que se puede acceder con las funciones `vignette()` y `demo()`.
#' 
#' 
#' 
#' # Importación y exportación de datos
#' 
#' ## Archivos de datos
#' A menudo es necesario leer y guardar los datos en algún formato tal que permita la importación y exportación de los mismos y su intercambio con otros programas o entornos.
#' 
#' Aunque puede importar y exportar archivos de muy distintos formatos, solo comentaremos como trabajar con
#' 
#' - archivos de texto, y
#' - archivos RData, el formato de almacenamiento de datos de R.
#' 
#' ----
#' 
#' Para los distintos ejemplos usaremos el conjunto de datos `flights` del paquete `nycflights13`.
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
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
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## str(flights)

#' 
#' ----
#' 
## ---- echo = FALSE--------------------------------------------------------------
str(flights)

#' 
#' ----
#' 
#' Partiremos de un subconjunto de `flights`, para ello empezamos segmentando el `data.frame`.
#' 
## ---- echo = TRUE---------------------------------------------------------------
fl_ny2ws <- flights[flights$dest %in% c("IAD","BWI"),
                    c("origin","dest","carrier","arr_delay","air_time")]
head(fl_ny2ws)

#' 
#' ## Escribir y leer archivos de texto
#' Para escribir un archivo de texto (delimitado) desde un `data.frame` se usa la función `write.table` o cualquiera de sus derivadas. Para leer, la función es `read.table`.
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## write.table(fl_ny2ws, file="fl_ny2ws.csv", sep=",", dec=".",
##             quote=TRUE, fileEncoding="UTF-8", row.names=FALSE)

#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## fl_ny2ws <- read.table("fl_ny2ws.csv", sep=",", dec=".",
##            quote="\"", fileEncoding="UTF-8", header=TRUE)

#' 
#' Si se desea o dispone de un archivo delimitado por tabulaciones, el carácter tabulador se indica como `\t`.
#' 
#' 
#' ## Escribir y leer archivos de datos de R
#' Para leer y guardar datos en el formato propio de R, se usan las funciones `save` y `load`. Estas permiten almacenar y recuperar cualquier conjunto de variables del entorno de trabajo. Al recuperarlas se recuperan con el mismo nombre con el que se almacenaron.  
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## save(fl_ny2ws, file="fl_ny2ws.rda")

#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## print(load("fl_ny2ws.rda"))   #print muestra el nombre de los objetos

#' 
#' 
#' # Manipulación avanzada de tablas de datos
#' 
#' ## Tabla de datos o *data frame*
#' Com se ha viso anteriormente, el *data frame* es el tipo de dato más usado para almacenar tablas de datos.
#' 
#' A menudo, para poder realizar gráficos y/o aplicar distintos procedimientos estadísticos, es necesario manipular la tabla de datos. A esto dedicaremos este apartado.
#' 
#' ----
#' 
#' Entre las operaciones habituales que haremos sobre una tabla de datos, hay
#' 
#' - renombrar filas o columnas,
#' - añadir columnas o filas,
#' - segmentar,
#' - cambiar el formato de un conjunto de datos,
#' - eliminar filas o columnas,
#' - crear tablas de datos de resumen, y
#' - unir tablas
#' 
#' Terminaremos este bloque introduciendo los funciones básicas del paquete `dplyr` que facilita la realización de algunas de estas operaciones.
#' 
#' ----
#' 
#' Partimos de una tabla de datos sintética...
#' 
## ---- echo = TRUE---------------------------------------------------------------
df <- data.frame(1:5, letters[1:5], c(rep("a", 3), rep("b", 2)))
df

#' 
#' ## Renombrar filas o columnas
#' 
## ---- echo = TRUE---------------------------------------------------------------
colnames(df) <- c("var1", "var2", "var3") 
rownames(df) <- paste("subject00", 1:5, sep = "")
df

#' 
#' ## Añadir columnas o filas
#' 
## ---- echo = TRUE---------------------------------------------------------------
df2 <- cbind(df, rnorm(5)) # añadir un vector al data frame
df2$var5 <- 5:1 # assignando valores a una nueva variable
df2

#' 
#' ----
#' 
## ---- echo = TRUE---------------------------------------------------------------
df2 <- rbind(df, list(6, "e", "b"))
df2

#' 
#' ## Segmentar
#' 
#' Existen tres formas básicas de seleccionar filas o columnas de una tablas de datos
#' 
#' - mediante índices numéricos,
#' - usando los nombres de filas y columnas, y
#' - mediante vectores lógicos
#' 
#' ## Segmentar -- mediante índices
#' 
## ---- echo = TRUE---------------------------------------------------------------
df[1:3,]
df[,c(1,3)]

#' 
#' ----
#' 
## ---- echo = TRUE---------------------------------------------------------------
df[-3,-2]

#' 
#' ## Segmentar -- usando nombres
#' 
## ---- echo = TRUE---------------------------------------------------------------
df[,"var2"]
df$var3
df[,c("var2","var3")]

#' 
#' ## Segmentar -- mediante vectores lógicos
#' 
## ---- echo = TRUE---------------------------------------------------------------
df[c(T,T,F,T,F), c(T,F,T)]
df[df[,1] == 3 | df[,3] == "b",]

#' 
#' ## Cambiar el formato de un conjunto de datos
#' 
#' Un mismo conjunto de datos puede representar en una tabla de datos de acuerdo con distintas oragnizaciones o formatos.
#' 
#' Se denomina formato ordenado, *tidy*, aquella organización en la que las filas representan individuos, las columnas, variables de medidas, y las intersecciones los valores de dichas medidas.  
#' 
#' ----
#' 
#' Cuando existen más de una variable que corresponde al mismo tipo de medida, entonces los datos pueden organizarse de dos formas
#' 
#' - usando una columna distinta para cada variable, formato ancho y *tidy*, o
#' - usando una columna para los valores y otra para indicar cuál es la variable representada en esta fila.
#' 
#' El formato más adecuado dependerá de los análisis y visualizaciones que quieran hacerse.
#' 
#' ----
#' 
#' Partiremos de un subconjunto de `flights`...
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
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
#' ----
#' 
## ---- echo = TRUE---------------------------------------------------------------
fl_ny2ws_W <- flights[flights$dest %in% c("IAD","BWI"),
                    c("origin","dest","carrier","arr_delay","dep_delay")]
fl_ny2ws_W <- cbind(key = 1:nrow(fl_ny2ws_W), fl_ny2ws_W)
head(fl_ny2ws_W)

#' 
#' ## Cambiar el formato de un conjunto de datos -- ancho a largo
## ---- echo = TRUE---------------------------------------------------------------
fl_ny2ws_L <- rbind(
  cbind(edge = rep("origin", nrow(fl_ny2ws_W)), fl_ny2ws_W[,c(1,4)],
        airport = fl_ny2ws_W[,2], delay = fl_ny2ws_W[,6]),
  cbind(edge = rep("dest", nrow(fl_ny2ws_W)), fl_ny2ws_W[,c(1,4)],
        airport = fl_ny2ws_W[,3], delay = fl_ny2ws_W[,5]))
colnames(fl_ny2ws_L) <- c("edge","key","carrier","airport","delay")

#' 
#' ----
#' 
## ---- echo = TRUE---------------------------------------------------------------
head(fl_ny2ws_L)
tail(fl_ny2ws_L)

#' 
#' ## Cambiar el formato de un conjunto de datos -- largo a ancho
## ---- echo = TRUE---------------------------------------------------------------
fl_ny2ws_W2p1 <- fl_ny2ws_L[fl_ny2ws_L$edge=="origin",]
fl_ny2ws_W2p2 <- fl_ny2ws_L[fl_ny2ws_L$edge=="dest",]
fl_ny2ws_W2p1 <- fl_ny2ws_W2p1[order(fl_ny2ws_W2p1$key),-1]
fl_ny2ws_W2p2 <- fl_ny2ws_W2p2[order(fl_ny2ws_W2p2$key),-c(1,2)]
fl_ny2ws_W2 <- cbind(fl_ny2ws_W2p1,fl_ny2ws_W2p2[,-1])
colnames(fl_ny2ws_W2) <- c("key", "carrier", "origin", "dep_delay",
                           "dest","arr_delay")


#' 
#' ----
#' 
## ---- echo = TRUE---------------------------------------------------------------
head(fl_ny2ws_W2)
tail(fl_ny2ws_W2)

#' 
#' ## Eliminar filas o columnas
#' La forma más habitual de eliminar filas o columnas es segmentando la tabla de datos. Sin embargo, una columna también puede eliminarse asignando la misma a `NULL`.
#' 
## ---- echo = TRUE---------------------------------------------------------------
df <- data.frame(1:5, letters[1:5], c(rep("a", 3), rep("b", 2)))
colnames(df) <- c("var1", "var2", "var3") 
rownames(df) <- paste("subject00", 1:5, sep = "")
df <- df[-2,]
df$var2 <- NULL

#' 
#' ----
#' 
## ---- echo = TRUE---------------------------------------------------------------
df

#' 
#' ----
#' 
#' Si lo que se desea es eliminar una variable del entorno de trabajo entonces se usa la función `rm`.
#' 
## ---- echo = TRUE---------------------------------------------------------------
rm(df)

#' 
#' ## Creación de resúmenes a partir de datos en formato ancho
## ---- echo = TRUE---------------------------------------------------------------
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
## ---- echo = TRUE---------------------------------------------------------------
sum_fl_ny2ws <- data.frame(edge=c("origin","dest"))
sum_fl_ny2ws$mean_delay <- by(fl_ny2ws_L$delay,fl_ny2ws_L$edge,
                              mean,na.rm=TRUE)
sum_fl_ny2ws$r_delay <- by(fl_ny2ws_L$delay,fl_ny2ws_L$edge,
                           function(x) diff(range(x,na.rm=TRUE)))
sum_fl_ny2ws$s_delay <- by(fl_ny2ws_L$delay,fl_ny2ws_L$edge,sd,na.rm=TRUE)
sum_fl_ny2ws

#' 
#' ## `dplyr`
#' Estas mismas operaciones que se han comentado más arriba, se pueden llevar a cabo también a partir de un sistema coherente de métodos creado como una gramática para la manipulación de datos, el paquete `dplyr` --que forma parte de `tidyverse`--.
#' 
#' http://dplyr.tidyverse.org/
#' 
#' 
## ---- echo = FALSE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("tidyverse")) {
  install.packages("tidyverse", repos="https://cloud.r-project.org/",
         quiet=TRUE, type="binary")
  library("tidyverse")
}

#' 
#' ----
#' 
#' En `dplyr` las operaciones básicas, se realizan mediante cuatro métodos:
#' 
#' - `select`: segmentar columnas,
#' - `filter`: segmentar filas,
#' - `mutate`: crear nuevas columnas, y
#' - `arrange`: ordenar.
#' 
#' Estas se encadenan usando el operador pipe, `%>%`.
#' 
#' ----
#' 
## ---- echo = TRUE---------------------------------------------------------------
df <- flights %>% dplyr::select(origin, dest, arr_delay) %>% 
  filter(origin == "LGA" & (dest == "IAD" | dest == "BWI")) %>%
  mutate(arr_delay_h=arr_delay/60) %>% 
  arrange(-arr_delay_h)
df

#' 
#' ----
#' 
#' Para la creación de resúmenes a partir de tablas en formato largo, es muy útil y cómoda la combinación `group_by` y `summarize`.
#' 
## ---- echo = TRUE---------------------------------------------------------------
df %>% group_by(dest) %>% summarise(mean_delay = mean(arr_delay, na.rm=TRUE))

#' 
#' ----
#' 
#' Por último, las conversiones entre formatos también pueden hacerse a partir de las funciones `spread` y `gather` del paquete `tidyr` --incluido en `tidyverse` y que se integra de forma natural en la gramática propuestas por `dplyr`--.
#' 
#' 
#' 
#' # Advanced Graphics in R
#' 
#' ## Gramática de gráficos (*GoG*)
#' La **gramática de gráficos** es una aproximación teórica al estudio de los componentes de un gráfico. De acuerdo con este análisis, un gráfico se puede construir mediante la especificación de un conjunto de capas y componentes que definen los datos, la asociación de los mismos a aspectos del gráfico, la especificación de la relación entre los valores de las variables de los datos con las del gráfico, la estructura geométrica del gráfico...
#' 
#' <p class="bibref">Wilkinson, L. (2006). The grammar of graphics. Springer Science & Business Media.</p>
#' 
#' ## `ggplot2`
#' `ggplot2` es un paquete de R que implementa de la gramática de gráficos.
#' 
#' <p class="bibref">Wickham, H. (2010). A layered grammar of graphics. *Journal of Computational and Graphical Statistics*, 19(1), 3-28.</p>
#' 
#' Referencias:
#' 
#' - http://ggplot2.tidyverse.org/reference/
#' - https://github.com/rstudio/cheatsheets/raw/master/data-visualization-2.1.pdf
#' - http://www.ling.upenn.edu/~joseff/avml2012/
#' - http://r-statistics.co/Complete-Ggplot2-Tutorial-Part1-With-R-Code.html
#' 
#' ----
#' 
#' `ggplot2` forma parte del paquete `tidyverse` (aunque también puede instalarse y cargarse autónomamente).
## ---- echo = TRUE, results = 'hide', message = FALSE, warning = FALSE, error = FALSE----
if(!require("tidyverse")) {
  install.packages("tidyverse", repos="https://cloud.r-project.org/",
         quiet=TRUE, type="binary")
  library("tidyverse")
}

#' 
#' ## `ggplot2` -- capas y elementos del gráfico
#' En `ggplot2` cada elemento gráfico que representa un conjunto de datos constituye una capa. Una o más capas constituyen un gráfico.
#' 
#' Cada capa queda definida mediante la especificación de sus elementos. Los principales son 
#' 
#' - datos,
#' - mapeado estético, y
#' - geometrías
#' 
#' ----
#' 
#' En `ggplot2`, los gráficos constituyen un objeto de R y se construyen de forma aditiva.
#' 
#' Por ejemplo,
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## grafico <- ggplot(data = anscombe,
##         mapping = aes(x = x1, y = y1))  # Datos y mapeado estético
## grafico <- grafico + geom_point()       # Geometría
## 
## grafico

#' 
#' ----
#' 
## ---- echo = FALSE--------------------------------------------------------------
grafico <- ggplot(data = anscombe,
        mapping = aes(x = x1, y = y1))  # Datos y mapeado estético 
grafico <- grafico + geom_point()       # Geometría

grafico

#' 
#' ## `ggplot2` -- datos
#' En `ggplot2`, el elemento `data` (datos) se introduce como primer argumento de la función `ggplot`. Debe corresponder a una tabla de datos o un tipo de datos convertible a tabla de datos.
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## grafico <- ggplot(data = anscombe,

#' 
#' ## `ggplot2` -- mapeado estético
#' El `mapping` (mapeado estético) corresponde al establecimiento de relaciones entre variables de los datos y variables del gráfico. Es el segundo argumento de la función `ggplot`y debe crearse con al función de apoyo `aes`. 
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
##         mapping = aes(x = x1, y = y1))

#' 
#' ----
#' 
#' Para variables cuantitativas, los mapeados más comunes corresponden a
#' 
#' - posiciones: `x`, `y`...
#' - tamaño: `size`
#' - color: `color`, `fill`
#' 
#' ----
#' 
#' Para variable cualitativas, los mapeados más frecuentes son
#' 
#' - posiciones: `x`, `y`...
#' - color: `color`, `fill`
#' - forma: `shape`
#' 
#' 
#' ## `ggplot2` -- geometrías
#' Las geometrías (`geom_`) indican la forma que debe tener el gráfico, es decir, cómo se articulan las variables del gráfico. Se añaden al gráfico sumándose al objeto creado por `ggplot`.
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## grafico <- grafico + geom_point()

#' 
#' ----
#' 
#' Son geometrías de uso común
#' 
#' - `geom_point`
#' - `geom_line`, `geom_vline`, `geom_hline`
#' - `geom_bar`
#' - `geom_histogram`
#' - `geom_boxplot`
#' 
#' Un resumen de las geometrías y su relación con las variables del gráfico que reconoce cada una de ellas figura en  https://github.com/rstudio/cheatsheets/raw/master/data-visualization-2.1.pdf.
#' 
#' 
#' ## Gráficos en `ggplot2`
#' Veremos cómo crear en `ggplot2` los gráficos más habituales, añadiendo algunas consideraciones para aquellos casos donde los gráficos realizados con `base` tienen prestaciones insuficientes.
#' 
#' - gráficos de dispersión,
#' - histogramas,
#' - diagramas de barras, y 
#' - diagramas de caja (*boxplot*).
#' 
#' De forma general, los gráficos en `ggplot2` se construyen a partir de tablas de datos (*data frames*), de los cuales se seleccionan las variables a representar.
#' 
#' ----
#' 
#' Usaremos 1000 datos del conjunto de datos `diamonds` para crear los distintos ejemplos.
#' 
## ---- echo = TRUE---------------------------------------------------------------
diaM <- diamonds[sample(1:nrow(diamonds),1000),]
str(diaM)

#' 
#' ## Gráficos en `ggplot2` -- gráfico de dispersión
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## ggplot(diaM, aes(x=carat,y=price)) + geom_point()

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE-------------------------------------------------
ggplot(diaM, aes(x=carat,y=price)) + geom_point()

#' 
#' ----
#' 
#' Añadiendo una tercera variable (`cut`) y modificando algunos aspectos de formato...
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## ggplot(diaM, aes(x=carat,y=price,color=cut)) +
##   geom_point(alpha=.8,shape=21,size=3)

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE-------------------------------------------------
ggplot(diaM, aes(x=carat,y=price,color=cut)) + 
  geom_point(alpha=.8,shape=21,size=3)

#' 
#' -----
#' 
#' Añadiendo líneas de tendencia...
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## ggplot(diaM, aes(x=carat,y=price,color=cut)) +
##   geom_point(alpha=.8,shape=21,size=3) +
##   geom_smooth(method="lm",se=FALSE)

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE-------------------------------------------------
ggplot(diaM, aes(x=carat,y=price,color=cut)) + 
  geom_point(alpha=.8,shape=21,size=3) +
  geom_smooth(method="lm",se=FALSE)

#' 
#' ## Gráficos en `ggplot2` -- histograma
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## ggplot(diaM, aes(x=price)) + geom_histogram(binwidth=1000)

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE-------------------------------------------------
ggplot(diaM, aes(x=price)) + geom_histogram(binwidth=1000)

#' 
#' ----
#' 
#' Y en función del corte...
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## ggplot(diaM, aes(x=price,fill=cut)) +
##   geom_histogram(position='dodge',binwidth=1000)

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE-------------------------------------------------
ggplot(diaM, aes(x=price,fill=cut)) +
  geom_histogram(position='dodge',binwidth=1000)

#' 
#' ----
#' 
#' En frecuencias relativas (por grupo)...
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## ggplot(diaM, aes(x=price,y=..density..,fill=cut)) +
##   geom_histogram(position='dodge',binwidth=1000)

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE-------------------------------------------------
ggplot(diaM, aes(x=price,y=..density..,fill=cut)) +
  geom_histogram(position='dodge',binwidth=1000)

#' 
#' ----
#' 
#' Quizás funcione mejor un gráfico de densidades...
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## ggplot(diaM, aes(x=price,fill=cut)) +
##   geom_density(alpha=.3)

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE-------------------------------------------------
ggplot(diaM, aes(x=price,fill=cut)) +
  geom_density(alpha=.3)

#' 
#' ## Gráficos en `ggplot2` -- diagrama de barras
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## ggplot(diaM, aes(x=clarity)) + geom_bar()

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE-------------------------------------------------
ggplot(diaM, aes(x=clarity)) + geom_bar()

#' 
#' ----
#' 
#' En función de la claridad...
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## ggplot(diaM, aes(x=clarity, fill=cut)) + geom_bar()

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE-------------------------------------------------
ggplot(diaM, aes(x=clarity, fill=cut)) + geom_bar()

#' 
#' ----
#' 
#' Para comparar entre frecuencias absolutas, funcionan mejor las barras separadas.
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## ggplot(diaM, aes(x=clarity, fill=cut)) + geom_bar(position="dodge")

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE-------------------------------------------------
ggplot(diaM, aes(x=clarity, fill=cut)) + geom_bar(position="dodge")

#' 
#' ----
#' 
#' Para comparar entre frecuencias relativas acumuladas, son mejores las barras apiladas en frecuencia relativa (para cada clase).
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## ggplot(diaM, aes(x=clarity, fill=cut)) + geom_bar(position="fill")

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE-------------------------------------------------
ggplot(diaM, aes(x=clarity, fill=cut)) + geom_bar(position="fill")

#' 
#' ----
#' 
#' Los diagramas de barras también pueden crearse a partir de tablas de datos agrupados. En este caso, debe indicarse qué variable es la `y` e incluir `stat="identity"` en el `geom_bar`.
#' 
#' ## Gráficos en `ggplot2` -- diagrama de caja
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## ggplot(diaM, aes(x=1, y=price)) + geom_boxplot()

#' 
#' *NOTA: Para el `geom_boxplot` se requieren dos variables. Si no hay variable independiente se puede incluir `x=1`.*
#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE-------------------------------------------------
ggplot(diaM, aes(x=1, y=price)) + geom_boxplot()

#' 
#' ----
#' 
#' Y en función del corte...
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## ggplot(diaM, aes(x=cut, y=price)) + geom_boxplot()

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE-------------------------------------------------
ggplot(diaM, aes(x=cut, y=price)) + geom_boxplot()

#' 
#' ----
#' 
#' El gráfico se puede mejorar mostrando todos los puntos, con una posición aleatorizada.
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## ggplot(diaM, aes(x=cut, y=price)) +
##   geom_boxplot(outlier.shape = NA) +
##   geom_jitter(shape = 21, alpha=.5,height=0,width=.2)

#' 
#' ----
#' 
## ---- echo = FALSE, eval = TRUE-------------------------------------------------
ggplot(diaM, aes(x=cut, y=price)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 21, alpha=.5,height=0,width=.2)

#' 
#' ----
#' 
#' O incluyendo un *violin plot* y un punto para la media...
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
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
## ---- echo = FALSE, eval = TRUE-------------------------------------------------
medias <- diaM %>% group_by(cut) %>%
  summarise(price=mean(price))

ggplot(diaM, aes(x=cut, y=price)) + 
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.1) +
  geom_point(data=medias,shape=3)

#' 
#' ## `ggplot2` -- elementos adicionales
#' 
#' Además de los elementos ya presentados (datos, mapeado estético y geometrías), otros elementos de `ggplot2` permiten controlar aspectos adicionales del gráfico, por ejemplo
#' 
#' - `scale_...`: controlan los aspectos relativos a la presentación de las variables del gráfico,
#' - `coord_...`: establecen el sistema de coordenadas usado en la geometría,
#' - `labs`: establece los títulos del gráfico,
#' - `theme`, `theme_...`: controlan la part del gráfico que no corresponde a datos (*non-data ink*), y
#' - `facet_`: permite la creació de secuencia de gráficos en función de una o dos variables
#' 
#' ----
#' 
#' Un ejemplo para terminar...
#' 
## ---- echo = TRUE, eval = FALSE-------------------------------------------------
## ggplot(diaM, aes(x=carat, y = price, shape = cut, col = clarity)) +
##   geom_point(alpha=.6) +
##   scale_x_continuous(breaks=1:3) +
##   scale_y_continuous(trans="log10") +
##   scale_color_brewer(palette="Spectral")+
##   facet_grid(cut ~ clarity) +
##   theme_bw() +
##   theme(legend.position = "none",text = element_text(size=10))
## 

#' 
#' ----
#' 
#' 
## ---- echo = FALSE--------------------------------------------------------------
ggplot(diaM, aes(x=carat, y = price, shape = cut, col = clarity)) +
  geom_point(alpha=.8) +
  scale_x_continuous(breaks=1:3) +
  scale_y_continuous(trans="log10") +
  scale_color_brewer(palette="Spectral")+
  facet_grid(cut ~ clarity) + 
  theme_bw() + 
  theme(legend.position = "none",text = element_text(size=10))


#' 

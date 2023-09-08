#Librerias generales usadas durante el TFM
library("ggplot2")
library("writexl")
library("readxl")
library("dplyr")
library(tidyverse)
library(reshape2)
library(gridExtra)
library(ggtext)

#Cargamos la tabla presencia/ausencia general sin anotación
tabla <- read_excel("/home/paablofdeez/TFM/presencia_ausencia_hb8_hb27.xlsx")
tabla <- subset(tabla,select = -2)

#Cargamos la anotación de la estirpe HB8 y eliminamos columnas innecesarias
ttha <- read.csv("/home/paablofdeez/TFM/tthafinlimpio.csv",header = FALSE,sep="\t")
ttha <- subset(ttha,select=-1)

#Cargamos la anotación de la estirpe HB27 y eliminamos valores NA
ttc <- read.csv("/home/paablofdeez/TFM/hb27limpio.csv",header=FALSE,sep="\t")
ttc <- na.omit(ttc)

#Renombramos las columans de ambos datasets para que sean iguales, con el fin de poder ejecutar el merge 
names(ttc) <- c("Start","Stop","anno")
names(ttha) <- c("Start","Stop","anno")

#Cargamos la matriz de anotación general
thermus <- read_excel("/home/paablofdeez/TFM/matriz_general/Thermus/matrizglobalthermus.xlsx")

#Fusionamos la matriz general con la anotación HB27 en función de las posiciones de inicio-fin de las secuencias (deberían ser identicas)
#con esto conseguimos añadir la columna "Gene" a la anotación HB27, que es la misma que presenta la matriz de presencia-ausencia, obteniedo
#una columna identica a través de la cual hacer el merge
thermus_anno_ttc <- merge(x=thermus,y=ttc,by=c("Start","Stop"),all.x=TRUE)

#Eliminamos columnas innecesarias y reanotamos
thermus_anno_ttc <- subset(thermus_anno_ttc, select=-4)
names(thermus_anno_ttc) <- c("Start","Stop","Gene","Annotation_hb27")

#Borramos valores NA
thermus_anno_ttc <- thermus_anno_ttc[!is.na(thermus_anno_ttc$Anno), ]

#Repetimos el proceso en la estirpe Hb8
thermus_anno_ttha <- merge(x=thermus,y=ttha,by=c("Start","Stop"),all.x=TRUE)

thermus_anno_ttha <- subset(thermus_anno_ttha, select=-4)
names(thermus_anno_ttha) <- c("Start","Stop","Gene","Annotation_HB8")

thermus_anno_ttha <- thermus_anno_ttha[!is.na(thermus_anno_ttha$Gene), ]

#Finalmente fusionamos ambas tablas para tener las anotaciones de ambas estirpes en la tabla de presencia-ausencia
tabla_ttha <- merge(x=tabla,y=thermus_anno_ttha,by="Gene",all.x=TRUE)
tabla_ttha <- subset(tabla_ttha,select=-219)

tabla_final <- merge(x=tabla_ttha,y=thermus_anno_ttc, by="Gene",all.x = TRUE)
tabla_final <- subset(tabla_final,select=-220)

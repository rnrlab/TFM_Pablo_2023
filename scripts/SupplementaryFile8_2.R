library(readr)
library("ggplot2")
library("writexl")
library("readxl")
library("dplyr")
library(tidyverse)
library(reshape2)

#Datos obtenidos a través de SAPPHIRE
datos <- read_csv("/home/paablofdeez/TFM/sapphire_resultados/deinococcus/finales.csv")
names(datos) <- c("ID proteina","Termofilia","Valor_Termofilia")

#SupplementaryFile1_2
deinococcus <- read_xlsx('/home/paablofdeez/TFM/matriz_general/Deinococcus/nuevamatriz/archivo_fasta_bakta/fastabaktaxlsx.xlsx')
ppanggolin <- read_xlsx('/home/paablofdeez/TFM/matriz_general/Deinococcus/nuevamatriz/archivos_ppanggolin/fastas_ppanggolin.xlsx')
ppanggolin <- subset(ppanggolin,select=c(3,4,8))
names(ppanggolin) <- c("Start","Stop","PPanGGolin_Type")

deinococcus <- merge(deinococcus, ppanggolin, by.x=c("Start","Stop"),by.y = c("Start","Stop"),all.x = TRUE)


coincidencias <- merge(datos,deinococcus, by.x="ID proteina",by.y="Locus Tag",all.x = TRUE)
coincidencias <- subset(coincidencias, select=c(1,2,3,12))

filtrados90p <- coincidencias %>%
  filter(Valor_Termofilia >= 0.9 & PPanGGolin_Type == "persistent")
n_filtrados90p <- (nrow(filtrados90p)/2934)*100

filtrados90s <- coincidencias %>%
  filter(Valor_Termofilia >= 0.9 & PPanGGolin_Type == "shell")
n_filtrados90s <- (nrow(filtrados90s)/2934)*100

filtrados90c <- coincidencias %>%
  filter(Valor_Termofilia >= 0.9 & PPanGGolin_Type == "cloud")
n_filtrados90c <- (nrow(filtrados90c)/2934)*100

filtrados90 <- data.frame(n_filtrados90c,n_filtrados90s,n_filtrados90p)
names(filtrados90) <- c("Cloud","Shell","Persistent")



filtrados80p <- coincidencias %>%
  filter(Valor_Termofilia >= 0.8 & Valor_Termofilia < 0.9 & PPanGGolin_Type == "persistent")
n_filtrados80p <- (nrow(filtrados80p)/2934)*100

filtrados80s <- coincidencias %>%
  filter(Valor_Termofilia >= 0.8 & Valor_Termofilia < 0.9 & PPanGGolin_Type == "shell")
n_filtrados80s <- (nrow(filtrados80s)/2934)*100

filtrados80c <- coincidencias %>%
  filter(Valor_Termofilia >= 0.8 & Valor_Termofilia < 0.9 & PPanGGolin_Type == "cloud")
n_filtrados80c <- (nrow(filtrados80c)/2934)*100

filtrados80 <- data.frame(n_filtrados80c,n_filtrados80s,n_filtrados80p)
names(filtrados80) <- c("Cloud","Shell","Persistent")






filtrados70p <- coincidencias %>%
  filter(Valor_Termofilia >= 0.7 & Valor_Termofilia < 0.8 & PPanGGolin_Type == "persistent")
n_filtrados70p <- (nrow(filtrados70p)/2934)*100

filtrados70s <- coincidencias %>%
  filter(Valor_Termofilia >= 0.7 & Valor_Termofilia < 0.8 & PPanGGolin_Type == "shell")
n_filtrados70s <- (nrow(filtrados70s)/2934)*100

filtrados70c <- coincidencias %>%
  filter(Valor_Termofilia >= 0.7 & Valor_Termofilia < 0.8 & PPanGGolin_Type == "cloud")
n_filtrados70c <- (nrow(filtrados70c)/2934)*100

filtrados70 <- data.frame(n_filtrados70c,n_filtrados70s,n_filtrados70p)
names(filtrados70) <- c("Cloud","Shell","Persistent")





filtrados60p <- coincidencias %>%
  filter(Valor_Termofilia >= 0.6 & Valor_Termofilia < 0.7 & PPanGGolin_Type == "persistent")
n_filtrados60p <- (nrow(filtrados60p)/2934)*100

filtrados60s <- coincidencias %>%
  filter(Valor_Termofilia >= 0.6 & Valor_Termofilia < 0.7 & PPanGGolin_Type == "shell")
n_filtrados60s <- (nrow(filtrados60s)/2934)*100

filtrados60c <- coincidencias %>%
  filter(Valor_Termofilia >= 0.6 & Valor_Termofilia < 0.7 & PPanGGolin_Type == "cloud")
n_filtrados60c <- (nrow(filtrados60c)/2934)*100

filtrados60 <- data.frame(n_filtrados60c,n_filtrados60s,n_filtrados60p)
names(filtrados60) <- c("Cloud","Shell","Persistent")


filtrados50p <- coincidencias %>%
  filter(Valor_Termofilia >= 0.5 & Valor_Termofilia < 0.6 & PPanGGolin_Type == "persistent")
n_filtrados50p <- (nrow(filtrados50p)/2934)*100

filtrados50s <- coincidencias %>%
  filter(Valor_Termofilia >= 0.5 & Valor_Termofilia < 0.6 & PPanGGolin_Type == "shell")
n_filtrados50s <- (nrow(filtrados50s)/2934)*100

filtrados50c <- coincidencias %>%
  filter(Valor_Termofilia >= 0.5 & Valor_Termofilia < 0.6 & PPanGGolin_Type == "cloud")
n_filtrados50c <- (nrow(filtrados50c)/2934)*100

filtrados40p <- coincidencias %>%
  filter(Valor_Termofilia >= 0.4 & Valor_Termofilia < 0.5 & PPanGGolin_Type == "persistent")
n_filtrados40p <- (nrow(filtrados40p)/2934)*100

filtrados40s <- coincidencias %>%
  filter(Valor_Termofilia >= 0.4 & Valor_Termofilia < 0.5 & PPanGGolin_Type == "shell")
n_filtrados40s <- (nrow(filtrados40s)/2934)*100

filtrados40c <- coincidencias %>%
  filter(Valor_Termofilia >= 0.4 & Valor_Termofilia < 0.5 & PPanGGolin_Type == "cloud")
n_filtrados40c <- (nrow(filtrados40c)/2934)*100


filtrados30p <- coincidencias %>%
  filter(Valor_Termofilia >= 0.3 & Valor_Termofilia < 0.4 & PPanGGolin_Type == "persistent")
n_filtrados30p <- (nrow(filtrados30p)/2934)*100

filtrados30s <- coincidencias %>%
  filter(Valor_Termofilia >= 0.3 & Valor_Termofilia < 0.4 & PPanGGolin_Type == "shell")
n_filtrados30s <- (nrow(filtrados30s)/2934)*100

filtrados30c <- coincidencias %>%
  filter(Valor_Termofilia >= 0.3 & Valor_Termofilia < 0.4 & PPanGGolin_Type == "cloud")
n_filtrados30c <- (nrow(filtrados30c)/2934)*100


filtrados20p <- coincidencias %>%
  filter(Valor_Termofilia >= 0.2 & Valor_Termofilia < 0.3 & PPanGGolin_Type == "persistent")
n_filtrados20p <- (nrow(filtrados20p)/2934)*100

filtrados20s <- coincidencias %>%
  filter(Valor_Termofilia >= 0.2 & Valor_Termofilia < 0.3 & PPanGGolin_Type == "shell")
n_filtrados20s <- (nrow(filtrados20s)/2934)*100

filtrados20c <- coincidencias %>%
  filter(Valor_Termofilia >= 0.2 & Valor_Termofilia < 0.3 & PPanGGolin_Type == "cloud")
n_filtrados20c <- (nrow(filtrados20c)/2934)*100


filtrados10p <- coincidencias %>%
  filter(Valor_Termofilia >= 0.1 & Valor_Termofilia < 0.2 & PPanGGolin_Type == "persistent")
n_filtrados10p <- (nrow(filtrados10p)/2934)*100

filtrados10s <- coincidencias %>%
  filter(Valor_Termofilia >= 0.1 & Valor_Termofilia < 0.2 & PPanGGolin_Type == "shell")
n_filtrados10s <- (nrow(filtrados10s)/2934)*100

filtrados10c <- coincidencias %>%
  filter(Valor_Termofilia >= 0.1 & Valor_Termofilia < 0.2 & PPanGGolin_Type == "cloud")
n_filtrados10c <- (nrow(filtrados10c)/2934)*100


filtrados0p <- coincidencias %>%
  filter(Valor_Termofilia >= 0 & Valor_Termofilia < 0.1 & PPanGGolin_Type == "persistent")
n_filtrados0p <- (nrow(filtrados0p)/2934)*100

filtrados0s <- coincidencias %>%
  filter(Valor_Termofilia >= 0 & Valor_Termofilia < 0.1 & PPanGGolin_Type == "shell")
n_filtrados0s <- (nrow(filtrados0s)/2934)*100

filtrados0c <- coincidencias %>%
  filter(Valor_Termofilia >= 0 & Valor_Termofilia < 0.1 & PPanGGolin_Type == "cloud")
n_filtrados0c <- (nrow(filtrados0c)/2934)*100



filtrados50 <- data.frame(n_filtrados50c,n_filtrados50s,n_filtrados50p)
names(filtrados50) <- c("Cloud","Shell","Persistent")

filtrados40 <- data.frame(n_filtrados40c,n_filtrados40s,n_filtrados40p)
names(filtrados40) <- c("Cloud","Shell","Persistent")

filtrados30 <- data.frame(n_filtrados30c,n_filtrados30s,n_filtrados30p)
names(filtrados30) <- c("Cloud","Shell","Persistent")

filtrados20 <- data.frame(n_filtrados20c,n_filtrados20s,n_filtrados20p)
names(filtrados20) <- c("Cloud","Shell","Persistent")

filtrados10 <- data.frame(n_filtrados10c,n_filtrados10s,n_filtrados10p)
names(filtrados10) <- c("Cloud","Shell","Persistent")

filtrados0 <- data.frame(n_filtrados0c,n_filtrados0s,n_filtrados0p)
names(filtrados0) <- c("Cloud","Shell","Persistent")

df <- rbind(filtrados90,filtrados80,filtrados70,filtrados60,filtrados50,filtrados40,filtrados30,filtrados20,filtrados10,filtrados0)
rownames(df) <- c("90-100 %","80-90 %","70-80 %", "60-70 %", "50-60 %","40-50 %","30-40 %","20-30% ","10-20 %","0-10 %")



df$Uvr <- rownames(df)

# Convertir el dataframe en formato largo
df_long <- reshape2::melt(df, id.vars = "Uvr", variable.name = "Subgrupo", value.name = "Valor")

# Definir los colores para cada subgrupo
colores <-c("Persistent" = "royalblue", "Shell" = "lightblue", "Cloud" = "lightgreen")

plot_grupos <- ggplot(df_long, aes(x = Uvr, y = Valor, fill = Subgrupo)) +
  labs(title = "Thermophilic character of proteins. 5044 Proteins") +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colores) +
  labs(x = " % Termophilic  ", y = "% Clusters ", fill = "Genome\n type") +
  scale_y_continuous(labels = scales::number_format(scale = 1e-0)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, face = "bold",angle=45,hjust=1), 
        axis.text.y = element_text(size = 12, face = "bold")) +
  geom_vline(xintercept = 4.5, colour = "black") +
  theme(axis.title.y = element_text(size=14, face= "bold")) +
  theme(axis.title.x = element_text(size=14, face= "bold")) +
  theme(legend.position = "none") +
  coord_cartesian(ylim=c(0,35))


print(plot_grupos)


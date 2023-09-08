library("ggplot2")
library("writexl")
library("readxl")
library("dplyr")
library(tidyverse)
library(reshape2)

#Lectura del data frame y preprocesamiento de los datos

#SupplementaryFile1_1
thermus <- read_excel('/home/paablofdeez/TFM/matriz_general/Thermus/matrizglobalthermus.xlsx')
thermus <- subset(thermus, select=-14)
thermus <- subset(thermus, select=-6)
names(thermus) <- c("Numero","Start","Stop","SequenceID","Type","Locus_Tag","Gene","Product","DbXrefs","Family","Contig","Strand","Nb_copy_in_org","Genome_class","Persistent_neighbors","Shell_neighbors","Cloud_neighbors")

#Filtrado por clusters
thermus_cluster <- subset(thermus, !duplicated(Family))

#Lectura del data frame y preprocesamiento de los datos

#SupplementaryFile1_2
deinococcus <- read_excel('/home/paablofdeez/TFM/matriz_general/Deinococcus/matrizdeino.xlsx')
deinococcus <- subset(deinococcus, select=-6)
deinococcus <- subset(deinococcus, select=-14)
deinococcus <- subset(deinococcus, select=-20)
deinococcus <- subset(deinococcus, select=-19)
deinococcus <- subset(deinococcus, select=-18)
names(deinococcus)<- c("Numero","Start","Stop","SequenceID","Type","Locus_Tag","Gene","Product","DbXrefs","Family","Contig","Strand","Nb_copy_in_org","Genome_class","Persistent_neighbors","Shell_neighbors","Cloud_neighbors")

#Filtrado por clusters
deinococcus_cluster <- subset(deinococcus,!duplicated(Family))

#Conteo de apariciones para los distintos grupos
thermus_cluster_persistent <- length(grep("persistent",thermus_cluster$Genome_class))
thermus_cluster_shell <- length(grep("shell",thermus_cluster$Genome_class))
thermus_cluster_cloud <- length(grep("cloud",thermus_cluster$Genome_class))


deinococcus_cluster_persistent <- length(grep("persistent",deinococcus_cluster$Genome_class))
deinococcus_cluster_shell <- length(grep("shell",deinococcus_cluster$Genome_class))
deinococcus_cluster_cloud <- length(grep("cloud",deinococcus_cluster$Genome_class))

#Construccion dataframes con los resultados
persistent_df <- data.frame(deinococcus_cluster_persistent,thermus_cluster_persistent)
names(persistent_df)<- c("Deinococcus","Thermus")

shell_df <- data.frame(deinococcus_cluster_shell,thermus_cluster_shell)
names(shell_df)<- c("Deinococcus","Thermus")

cloud_df <- data.frame(deinococcus_cluster_cloud,thermus_cluster_cloud)
names(cloud_df)<- c("Deinococcus","Thermus")

#Agrupacion de resultados
grupo <- rbind(cloud_df,shell_df,persistent_df)                    
rownames(grupo)<- c("Cloud","Shell","Persistent")

grupo$Deinococcus <- 100*(grupo$Deinococcus/94570)
grupo$Deinococcus <- grupo$Deinococcus + 0.38631
grupo$Thermus <- 100*(grupo$Thermus/255160)

grupo_transposed <- t(grupo)
df <- as.data.frame(grupo_transposed)

#Se añade una columna con los rownames
df$Uvr <- rownames(df)

# Convertimos el df en formato largo para poder representarlo
df_long <- reshape2::melt(df, id.vars = "Uvr", variable.name = "Subgrupo", value.name = "Valor")

# Se definen los colores utilizados para los Genomy_Type
colores <- c("Persistent" = "royalblue", "Shell" = "lightblue", "Cloud" = "lightgreen")


positions= c("Thermus","Deinococcus")

plot_grupos <- ggplot(df_long, aes(x = Uvr, y = Valor, fill = Subgrupo)) +
  labs(title = "Genome classes by family ") +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = colores) +
  labs( x="",y="",fill = "Genome \ntype") +
  scale_y_continuous(labels = scales::number_format(scale = 1e-0)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, face = "bold"), 
        axis.text.y = element_text(size = 12, face = "bold"))+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.9,hjust=0.4)) +
  scale_x_discrete(limits = positions)

print(plot_grupos)





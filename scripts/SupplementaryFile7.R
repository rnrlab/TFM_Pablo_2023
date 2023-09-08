library("ggplot2")
library("writexl")
library("readxl")
library("dplyr")
library(tidyverse)
library(reshape2)
library(gridExtra)
library(ggtext)

#Cargamos los datos de Thermaceae, eliminamos columnas innecesarias y renombramos

#SupplementaryFile4_1
thermus <- read_excel("/home/paablofdeez/TFM/matriz_general/Thermus/ThermusMatrizGlobal_AnotacionKEGG.xlsx")
thermus <- subset(thermus, select=-2)
colnames(thermus)[colnames(thermus) == "Gene.y"] <- "Family"
#Filtramos por clusters
thermus_cluster <- subset(thermus,!duplicated(Family))

#Se realizan los conteos de clusters en base a la función
#Persistent
thermus_persistent <- thermus_cluster %>%
  filter(PPanGGolin_Type == "persistent")

thermus_persistent_dna <- length(grep("DNA replication",thermus_persistent$Function_2))
thermus_persistent_BER <- length(grep("Base excision repair",thermus_persistent$Function_2))
thermus_persistent_NER <- length(grep("Nucleotide excision repair",thermus_persistent$Function_2))
thermus_persistent_mismatch <- length(grep("Mismatch repair",thermus_persistent$Function_2))
thermus_persistent_homologous <- length(grep("Homologous recombination",thermus_persistent$Function_2))
thermus_persistent_brite <- length(grep("DNA repair and recombination proteins",thermus_persistent$Function_2))

#Shell
thermus_shell <- thermus_cluster %>%
  filter(PPanGGolin_Type == "shell")

thermus_shell_dna <- length(grep("DNA replication",thermus_shell$Function_2))
thermus_shell_BER <- length(grep("Base excision repair",thermus_shell$Function_2))
thermus_shell_NER <- length(grep("Nucleotide excision repair",thermus_shell$Function_2))
thermus_shell_mismatch <- length(grep("Mismatch repair",thermus_shell$Function_2))
thermus_shell_homologous <- length(grep("Homologous recombination",thermus_shell$Function_2))
thermus_shell_brite <- length(grep("DNA repair and recombination proteins",thermus_shell$Function_2))

#Cloud
thermus_cloud <- thermus_cluster %>%
  filter(PPanGGolin_Type == "cloud")

thermus_cloud_dna <- length(grep("DNA replication",thermus_cloud$Function_2))
thermus_cloud_BER <- length(grep("Base excision repair",thermus_cloud$Function_2))
thermus_cloud_NER <- length(grep("Nucleotide excision repair",thermus_cloud$Function_2))
thermus_cloud_mismatch <- length(grep("Mismatch repair",thermus_cloud$Function_2))
thermus_cloud_homologous <- length(grep("Homologous recombination",thermus_cloud$Function_2))
thermus_cloud_brite <- length(grep("DNA repair and recombination proteins",thermus_cloud$Function_2))

#creacion de los dataframes
thermus_persistent_df <- data.frame(thermus_persistent_dna,thermus_persistent_BER,thermus_persistent_NER,thermus_persistent_mismatch,thermus_persistent_homologous,thermus_persistent_brite)
thermus_shell_df <- data.frame(thermus_shell_dna, thermus_shell_BER, thermus_shell_NER, thermus_shell_mismatch, thermus_shell_homologous, thermus_shell_brite)
thermus_cloud_df <- data.frame(thermus_cloud_dna, thermus_cloud_BER, thermus_cloud_NER, thermus_cloud_mismatch, thermus_cloud_homologous, thermus_cloud_brite)
names(thermus_persistent_df) <- c("uno","dos","tres","cuatro","cinco","seis")
names(thermus_shell_df) <- c("uno","dos","tres","cuatro","cinco","seis")
names(thermus_cloud_df) <- c("uno","dos","tres","cuatro","cinco","seis")
rownames(thermus_persistent_df) <- c("Persistent")
rownames(thermus_shell_df) <- c("Shell")
rownames(thermus_cloud_df)<- c("Cloud")
thermus_df <- rbind(thermus_cloud_df,thermus_shell_df,thermus_persistent_df)
total_cloud= sum(thermus_cloud_df)
total_persistent=sum(thermus_persistent_df)
total_shell=sum(thermus_shell_df)
total= total_cloud+total_shell+total_persistent

thermus_df$uno <- 100 * (thermus_df$uno /total)
thermus_df$dos <- 100 * (thermus_df$dos /total)
thermus_df$tres <- 100 * (thermus_df$tres /total)
thermus_df$cuatro<- 100 * (thermus_df$cuatro /total)
thermus_df$cinco<- 100 * (thermus_df$cinco/total)
thermus_df$seis <- 100 * (thermus_df$seis /total)
names(thermus_df) <- c("DNA \nreplication","Base excision \nrepair","Nucleotide excision \nrepair","Mismatch \nrepair","Homologous \nrecombination","DNA repair and \nrecombination")
thermus_family_transposed <- t(thermus_df)


# Convertir el dataframe transpuesto en un nuevo dataframe
df <- as.data.frame(thermus_family_transposed)

# Añadir una columna con los nombres de las columnas originales
df$Uvr <- rownames(df)

# Convertir el dataframe en formato largo
df_long <- reshape2::melt(df, id.vars = "Uvr", variable.name = "Subgrupo", value.name = "Valor")

# Definir los colores para cada subgrupo
colores <- c("Persistent" = "royalblue", "Shell" = "lightblue", "Cloud" = "green")

# Generar el gráfico de barras con apilamiento
plot_thermus <- ggplot(df_long, aes(x = Uvr, y = Valor, fill = Subgrupo)) +
  labs(title = "Thermus KEGG Annotation ") +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colores) +
  labs(x = "Cellular function", y = "", fill = "Genome\n type") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 15,face="bold"))+
  theme(axis.text.y = element_text(size = 14)) +
  coord_cartesian(ylim=c(1,65))

print(plot_thermus)


#DEINOCOCCUS

#SupplementaryFile5_2
deinococcus <- read.csv("/home/paablofdeez/TFM/matriz_general/Deinococcus/matrizkeggdeinococcus.csv")


deinococcus_cluster <- subset(deinococcus,!duplicated(Family))
colnames(deinococcus_cluster)[colnames(deinococcus_cluster)=="partition"] <- "PPanGGolin_Type"

deinococcus_persistent <- deinococcus_cluster %>%
  filter(PPanGGolin_Type == "persistent")

deinococcus_persistent_dna <- length(grep("DNA replication",deinococcus_persistent$Function_3))
deinococcus_persistent_BER <- length(grep("Base excision repair",deinococcus_persistent$Function_3))
deinococcus_persistent_NER <- length(grep("Nucleotide excision repair",deinococcus_persistent$Function_3))
deinococcus_persistent_mismatch <- length(grep("Mismatch repair",deinococcus_persistent$Function_3))
deinococcus_persistent_homologous <- length(grep("Homologous recombination",deinococcus_persistent$Function_3))
deinococcus_persistent_brite <- length(grep("DNA repair and recombination proteins",deinococcus_persistent$Function_3))


deinococcus_shell <- deinococcus_cluster %>%
  filter(PPanGGolin_Type == "shell")

deinococcus_shell_dna <- length(grep("DNA replication",deinococcus_shell$Function_3))
deinococcus_shell_BER <- length(grep("Base excision repair",deinococcus_shell$Function_3))
deinococcus_shell_NER <- length(grep("Nucleotide excision repair",deinococcus_shell$Function_3))
deinococcus_shell_mismatch <- length(grep("Mismatch repair",deinococcus_shell$Function_3))
deinococcus_shell_homologous <- length(grep("Homologous recombination",deinococcus_shell$Function_3))
deinococcus_shell_brite <- length(grep("DNA repair and recombination proteins",deinococcus_shell$Function_3))


deinococcus_cloud <- deinococcus_cluster %>%
  filter(PPanGGolin_Type == "cloud")

deinococcus_cloud_dna <- length(grep("DNA replication",deinococcus_cloud$Function_3))
deinococcus_cloud_BER <- length(grep("Base excision repair",deinococcus_cloud$Function_3))
deinococcus_cloud_NER <- length(grep("Nucleotide excision repair",deinococcus_cloud$Function_3))
deinococcus_cloud_mismatch <- length(grep("Mismatch repair",deinococcus_cloud$Function_3))
deinococcus_cloud_homologous <- length(grep("Homologous recombination",deinococcus_cloud$Function_3))
deinococcus_cloud_brite <- length(grep("DNA repair and recombination proteins",deinococcus_cloud$Function_3))



deinococcus_persistent_df <- data.frame(deinococcus_persistent_dna,deinococcus_persistent_BER,deinococcus_persistent_NER,deinococcus_persistent_mismatch,deinococcus_persistent_homologous,deinococcus_persistent_brite)
deinococcus_shell_df <- data.frame(deinococcus_shell_dna, deinococcus_shell_BER, deinococcus_shell_NER, deinococcus_shell_mismatch, deinococcus_shell_homologous, deinococcus_shell_brite)
deinococcus_cloud_df <- data.frame(deinococcus_cloud_dna, deinococcus_cloud_BER, deinococcus_cloud_NER, deinococcus_cloud_mismatch, deinococcus_cloud_homologous, deinococcus_cloud_brite)
names(deinococcus_persistent_df) <- c("uno","dos","tres","cuatro","cinco","seis")
names(deinococcus_shell_df) <- c("uno","dos","tres","cuatro","cinco","seis")
names(deinococcus_cloud_df) <- c("uno","dos","tres","cuatro","cinco","seis")
rownames(deinococcus_persistent_df) <- c("Persistent")
rownames(deinococcus_shell_df) <- c("Shell")
rownames(deinococcus_cloud_df)<- c("Cloud")
deinococcus_df <- rbind(deinococcus_cloud_df,deinococcus_shell_df,deinococcus_persistent_df)

total_cloud= sum(deinococcus_cloud_df)
total_persistent=sum(deinococcus_persistent_df)
total_shell=sum(deinococcus_shell_df)
total= total_cloud+total_shell+total_persistent

deinococcus_df$uno <- 100 * (deinococcus_df$uno /total)
deinococcus_df$dos <- 100 * (deinococcus_df$dos /total)
deinococcus_df$tres <- 100 * (deinococcus_df$tres /total)
deinococcus_df$cuatro<- 100 * (deinococcus_df$cuatro /total)
deinococcus_df$cinco<- 100 * (deinococcus_df$cinco/total)
deinococcus_df$seis <- 100 * (deinococcus_df$seis /total)


names(deinococcus_df) <- c("DNA \nreplication","Base excision \nrepair","Nucleotide excision \nrepair","Mismatch \nrepair","Homologous \nrecombination","DNA repair and \nrecombination")





deinococcus_family_transposed <- t(deinococcus_df)

# Convertir el dataframe transpuesto en un nuevo dataframe
df <- as.data.frame(deinococcus_family_transposed)

# Añadir una columna con los nombres de las columnas originales
df$Uvr <- rownames(df)

# Convertir el dataframe en formato largo
df_long <- reshape2::melt(df, id.vars = "Uvr", variable.name = "Subgrupo", value.name = "Valor")

# Definir los colores para cada subgrupo
colores <- c("Persistent" = "royalblue", "Shell" = "lightblue", "Cloud" = "green")

# Generar el gráfico de barras con apilamiento
plot_deinococcus <- ggplot(df_long, aes(x = Uvr, y = Valor, fill = Subgrupo)) +
  labs(title = "Deinococcus KEGG Annotation") +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colores) +
  labs(x = "Cellular function", y = "",fill="Genome\n type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15,face = "bold"))+
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 14)) +
  coord_cartesian(ylim=c(1,65))

print(plot_deinococcus)


grid.arrange(plot_thermus,plot_deinococcus,ncol=2)





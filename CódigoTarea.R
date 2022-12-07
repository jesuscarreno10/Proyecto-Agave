#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("phyloseq") # Install phyloseq

#install.packages(c("RColorBrewer", "patchwork")) 

library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")
library("vegan")
library("phyloseq")
library("ggplot2")
library("igraph")
library("vegan")

# Nos ubicamos en la carpeta adecuada e importamos los datos del .biom

setwd("~/ProyectoAgave/DatosProyectoAgave/Total/")
merged_metagenomes <- import_biom("Total.biom")

merged_metagenomes

View(merged_metagenomes@tax_table@.Data)


#Recortamos información de la tabla anterior y renombramos las columnas 

merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4) #cut the firts character of tax
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
unique(merged_metagenomes@tax_table@.Data[,"Phylum"])

View(merged_metagenomes@tax_table@.Data)


#Corremos la tabla de metadatos

met.sol <- read.csv2("MetaDatos.csv", header =  TRUE, row.names = 1, sep = ",")


#Creamos un archivo phyloseq con la tabla de metadatos y datos

met.sol <- sample_data(met.sol)
met.sol
merged_metagenomes<- merge_phyloseq(merged_metagenomes, met.sol)
merged_metagenomes #Archivo Phyloseq

View(merged_metagenomes@tax_table@.Data)
View(merged_metagenomes@otu_table@.Data)
View(merged_metagenomes@sam_data@.Data)


#Tomamos solo las muestras cuya secuenciación fueron hechos por Shotgun

merged_metagenomes<-subset_samples(merged_metagenomes, seq_type == "Shotgun") 


#Ahora realizamos anilisis gráfico de las medidas de Alpha-diversidad 

plot_richness(merged_metagenomes, x = "plant_part", color = "plant_part")

plot_richness(merged_metagenomes, measures = c("Chao1", "Shannon"), x="plant_part", color = "plant_part")

plot_richness(merged_metagenomes, x = "host_name", color = "host_name")


#Solo nos quedamos con las muestras de Agave

agave_salmiana<-subset_samples(merged_metagenomes, host_name == "Agave Salmiana") 

agave_tequilana<-subset_samples(merged_metagenomes, host_name == "Agave Tequilana")

agave<-merge_phyloseq(agave_salmiana, agave_tequilana)


#Hacemos algunos plots de riqueza de nuestro intereses

plot_richness(agave, x="host_name", color = "host_name")

plot_richness(agave, measures = c("Chao1", "Shannon"), x="plant_part", color = "plant_part")

plot_richness(agave, measures = c("Chao1", "Shannon"), x="host_name", color = "host_name")

plot_bar(agave, fill = "host_name")


#Mapa de calor de los taxones Top 5 de las muestras de agave

TopNGenusAgave<-names(sort(taxa_sums(agave), TRUE)[1:5])
Top5GenusAgave<-prune_taxa(TopNGenusAgave,agave) #Generamos un objeto phyloseq con solo con los 5 taxones que escogimos
plot_heatmap(Top5GenusAgave) #Gráficamos


#Mapa de calor de los taxones Top 5 de las muestras de todas las especies

TopNGenus_merged<-names(sort(taxa_sums(merged_metagenomes), TRUE)[1:5])
Top5Genus_merged<-prune_taxa(TopNGenus_merged,merged_metagenomes) #Generamos un objeto phyloseq con solo con los 5 taxones que escogimos
plot_heatmap(Top5Genus_merged) #Gráficamos


#Veamos con las medidas de Bray-Curtis con el método NMDS para las muestras de Agave
#Cambiando los rangos de colores

p<-plot_heatmap(Top5GenusAgave, "NMDS", "bray", low = "#000033", high="#CCFF66")
p


#Gráfico de redes con respecto a la municipalidad de las muestras

set.seed(123)

ig<-make_network(merged_metagenomes, max.dist = 0.8)
plot_network(ig, merged_metagenomes, color = "municipality_town", shape = "municipality_town")

#Abundancias relativas

top_5_relativo_merged <-transform_sample_counts(Top5Genus_merged,function(x) x/sum(x)*100)

plot_bar(top_5_relativo_merged, fill = "host_name", facet_grid = ~host_name)

top_5_relativo_merged_df <- psmelt(top_5_relativo_merged)

top_5_relativo_merged_df$OTU <-as.character(top_5_relativo_merged_df$OTU)

top_5_relativo_merged_df$Sample <- as.factor(top_5_relativo_merged_df$Sample)

OTU_colors_merged <- brewer.pal(length(levels(top_5_relativo_merged_df$OTU)),"Dark2")

#globla_rel_merged_plot <- ggplot(top_5_relativo_merged_df, aes(x=Sample, y=Abundance, fill=OTU))+geom_bar(aes(), stat = "identity", position = "stack")+scale_fill_manual(values = OTU_colors_merged)+labs(x="Sample", y="Relative abundance")+facet_grid(~host_name, scales = "free_x")+ggtitle("Relative abundance by host_name")#+theme(axis.title = element_text(size=14), axis.title.x = element_text(size = 10))

#globla_rel_merged_plot



#Mapa de calor de los 5 taxones top en los merged

plot_heatmap(top_5_relativo_merged)


#Agrupamiento (Clústers) Beta diversidad

#Agrupamiento sencillo

abund_table_norm <- decostand(merged_metagenomes@otu_table@.Data, "normalize")

abund_table_norm <-t(abund_table_norm)

bc_dist <- vegdist(abund_table_norm, method = "bray")

cluster_single <- hclust(bc_dist, method = 'single')

plot(cluster_single)


#Agrpamiento completo

#Bray-Curtis

cluster_complete <- hclust(bc_dist, method = 'complete')

plot(cluster_complete)


#Jaccard

jaccard <-vegdist(abund_table_norm, method = "jaccard")

jaccard

cluster_jaccard <- hclust(jaccard, method = 'complete')

plot(cluster_jaccard)


#Sorensen

sorensen <- vegdist(abund_table_norm, binary = TRUE)

sorensen

cluster_sorensen <-hclust(sorensen, method = 'complete')

plot(cluster_sorensen)

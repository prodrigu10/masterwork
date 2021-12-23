###########################################################################################################################
# # # # # # # # # # # # # # # # # # # # # # # #                         # # # # # # # # # # # # # # # # # # # # # # # # # # 
############################################### TRABAJO FINAL DE MÁSTER ###################################################
# # # # # # # # # # # # # # # # # # # # # # # #                         # # # # # # # # # # # # # # # # # # # # # # # # # # 
###########################################################################################################################
# #                                                                                                                     # # 
### Implementación de un nuevo método de genotipación para organismos no modelo y su comparación con un método de novo  ###
# #                                                                                                                     # # 
###########################################################################################################################







#A CONTINUACIÓN SE RECOGEN TODAS LAS LÍNEAS DE CÓDIGO USADAS PARA LLEVAR ACABO TODOS LOS ANÁLISIS Y OTRAS ACCIONES PRESENTES EN EL TFM PRESENTE.

#EN EL PRESENTE DOCUMENTO HAY CAMBIOS DEL DIRECTORIO DE TRABAJO YA QUE, PARA LLEVAR ACABO CADA UNO DE LOS APARTADOS DEL TRABAJO,
#SE USÓ UN DOCUMENTO DIFERENTE Y PARA PONER TODO EL CÓDIGO EN ESTE DOCUMENTO SE HA TENIDO QUE ESPECIFICAR EL DOCUMENTO PROCEDENTE EN CADA APARTADO. 





####################################
#################################################################################################################
#################################### IGUALAR FORMATOS: ARCHIVOS VCF

#############

# ESTABLECER EL MISMO FORMATO PARA LAS TRES PRIMERAS COLUMNAS

#############

#Establecemos el WORKING DIRECTORY
setwd("C:/Users/José Ramiro Guinea/Desktop/Master/TFM/DATOSSSSS")

#Cargamos los archivos necesarios: fil_genotyping y indices_POS
datos <- as.data.frame(read.table("fil_genotyping.vcf",comment.char = "",header = TRUE))

indices <- as.data.frame(read.table("indices_POS.bed",col.names = c("N.Locus","inicio","final")))


#Creamos un vector vacío donde poder albergar los loci de nuestros snps
x.chrom <- c()

#Creamos un vector donde recogeremos los datos de ID de cada SNPS
id <- c()

#Creamos un vector donde recogeremos los datos de nueva posicion
new_pos <- c()


#Creamos un bucle que establecerá el locus donde está cada snp presente
for (i in 1:length(datos$X.CHROM)) {
  #Determinamos la posición absoluta del snp dentro del "genoma"
  pos <- datos$POS[i]
  
  #Bucle para probar si está dentro del rango de cada loci
  for (e in 1:length(indices$N.Locus)) {
    
    #Determinamos el rango de cada loci
    rango <- c(indices$inicio[e]:indices$final[e])
    
    #¿Está la posición del snp dentro del rango?
    if(pos %in% rango){
      x.chrom[i] <- sub(".*_","",indices$N.Locus[e])
      id[i] <- paste(sub(".*_","",indices$N.Locus[e]),":",pos-indices$inicio[e], sep = "")
      new_pos[i] <- pos-indices$inicio[e]
      break
    }
    
  }
}

#Cambiamos las colmunas que X.Chrom, POS y ID
datos$X.CHROM <- x.chrom
datos$POS <- new_pos
datos$ID <- id

#Cambio el nombre X.CHROM por #CHROM
colnames(datos)[1] <- "#CHROM"




###########

#SELECCIONAMOS EL PRIMER SNP DE CADA LOCI

###########



#En este apartado se filtran todos aquellos SNPs detectados que no sean el primero de cada loci, de esta forma se acaban de establecer las mismas condiciones para ambos archivos VCF

#Seleccionamos el primer SNP de cada LOCI
loci_list <- c() #Guardamos el nombre del LOCI 
e <- 1
loci_id <- c() #Guardamos el id del LOCI

for (i in 1:length(datos$`#CHROM`)) {
  if(datos$`#CHROM`[i] %in% loci_list)
    next
  else {
    loci_list[e] <- datos$`#CHROM`[i]
    loci_id[e] <- i
    e <- e+1
  }
}


#Creamos un nuevo df con los primeros SNP por LOCI 
datos_finales <- datos[loci_id,]


#Creamos el nuevo archivo vcf
write.table(datos_finales, file = "pau_first_snp",sep = "\t",
            row.names = FALSE, quote = F)





####################################
#################################################################################################################
#################################### ANALISIS DESCRIPTIVO DEL ALINEAMIENTO

################

# EVALUACIÓN DE LAS CARACTERÍSTICAS DEL ALINEAMIENTO DE LA GENOTIPACIÓN CON GENOMA DE REFERENCIA

################





#WORKING DIRECTORY
setwd("C:/Users/José Ramiro Guinea/Desktop/Master/TFM")

#Cargamos paquetes
library(readr)
library(ggplot2)

#Cargamos los datos
datos_SAM_INFO <- read_csv("C:/Users/José Ramiro Guinea/Desktop/Master/TFM/datos_SAM_INFO.csv")



#Convertimos la variable porcentaje en numeric
e <- 1 
for (i in datos_SAM_INFO$PercentAlignment) {
  datos_SAM_INFO$PercentAlignment[e] <- as.numeric(substr(i,1,nchar(i)-1))
  e <- e+1
}
datos_SAM_INFO$PercentAlignment <- as.numeric(datos_SAM_INFO$PercentAlignment)





#Aplicamos un boxplot a los datos para ver outliers
ggplot(data = datos_SAM_INFO,
       aes(x = "",y = datos_SAM_INFO$PercentAlignment))+
  geom_boxplot(fill = 3,
               color = 1,
               outlier.colour = 6)+
  ylab("Porcentaje de alineamiento")+
  theme(axis.title = element_blank())+
  ggtitle("Porcentaje de alineamiento")



#Buscamos outliers con el rango intercuantil


riq <- IQR(datos_SAM_INFO$PercentAlignment) #Definimos el Rango Intercuartílico

cuantiles <- quantile(datos_SAM_INFO$PercentAlignment, c(0.25, 0.5, 0.75), type = 7) #Definimos los cuartiles

#Establecemos los límites para detectar outliers
out_min <- as.numeric(cuantiles[1])-1.5*riq
out_max <- as.numeric(cuantiles[3])+1.5*riq

#Creamos una variable de almacenamiento de outliers
muestras_out <- c()

#Contadores para el bucle
n <- 1
e <- 1

#Bucle de obtencion de outliers
for (i in datos_SAM_INFO$PercentAlignment) {
  
  #Out superiores
  if(i>out_max){
    muestras_out[n] <- datos_SAM_INFO$Muestra[e]
    n<-n+1
    
    cat("Outlier Superior:\t",datos_SAM_INFO$Muestra[e],"\t",i,"\n")
    
    e<-e+1
    
  }
  #Out inferiores
  if(i<out_min){
    muestras_out[n] <- datos_SAM_INFO$Muestra[e]
    n<-n+1
    
    cat("Outlier Inferior:\t",datos_SAM_INFO$Muestra[e],"\t",i,"\n")
    
    e<-e+1
    
  } else {
    e<-e+1
  }
  
}

#Visualizamos outliers
muestras_out



#ANALIZAMOS LOS READS QUE SE LOCALIZAN EN MÁS DE UN LOCUS

#Observamos un resumen de los datos
summary(datos_SAM_INFO$AlignedPlus)



#Visualizamos si existen outliers
ggplot(data = datos_SAM_INFO,
       aes(x = "",y = datos_SAM_INFO$AlignedPlus))+
  geom_boxplot(fill = 7,
               color = 1,
               outlier.colour = 6)+
  theme(axis.title = element_blank())+
  ggtitle("Alineamientos múltiples")


#Buscamos outliers con el rango intercuantil


riq <- IQR(datos_SAM_INFO$AlignedPlus) #Definimos el Rango Intercuartílico

cuantiles <- quantile(datos_SAM_INFO$AlignedPlus, c(0.25, 0.5, 0.75), type = 7) #Definimos los cuartiles

#Establecemos los límites para detectar outliers
out_min <- as.numeric(cuantiles[1])-1.5*riq
out_max <- as.numeric(cuantiles[3])+1.5*riq

#Creamos una variable de almacenamiento de outliers
muestras_mult_out <- c()

#Contadores para el bucle
n <- 1
e <- 1

#Bucle de obtencion de outliers
for (i in datos_SAM_INFO$AlignedPlus) {
  
  #Out superiores
  if(i>out_max){
    muestras_mult_out[n] <- datos_SAM_INFO$Muestra[e]
    n<-n+1
    
    cat("Outlier Superior:\t",datos_SAM_INFO$Muestra[e],"\t",i,"\n")
    
    e<-e+1
    
  }
  #Out inferiores
  if(i<out_min){
    muestras_mult_out[n] <- datos_SAM_INFO$Muestra[e]
    n<-n+1
    
    cat("Outlier Inferior:\t",datos_SAM_INFO$Muestra[e],"\t",i,"\n")
    
    e<-e+1
    
  } else {
    e<-e+1
  }
  
}

#Visualizamos outliers
muestras_mult_out


#Lo mismo para el porcentaje de alineamientos múltiples

#Valores en porcentaje
datos_SAM_INFO$`AlignedPlus%` <- datos_SAM_INFO$AlignedPlus*100/datos_SAM_INFO$TotalReads

#Visualizamos si existen outliers
ggplot(data = datos_SAM_INFO,
       aes(x = "",y = datos_SAM_INFO$`AlignedPlus%`))+
  geom_boxplot(fill = 5,
               color = 1,
               outlier.colour = 6)+
  theme(axis.title = element_blank())+
  ggtitle("Alineamientos múltiples (%)")


#Buscamos outliers con el rango intercuantil


riq <- IQR(datos_SAM_INFO$`AlignedPlus%`) #Definimos el Rango Intercuartílico

cuantiles <- quantile(datos_SAM_INFO$`AlignedPlus%`, c(0.25, 0.5, 0.75), type = 7) #Definimos los cuartiles

#Establecemos los límites para detectar outliers
out_min <- as.numeric(cuantiles[1])-1.5*riq
out_max <- as.numeric(cuantiles[3])+1.5*riq

#Creamos una variable de almacenamiento de outliers
muestras_multp_out <- c()

#Contadores para el bucle
n <- 1
e <- 1

#Bucle de obtencion de outliers
for (i in datos_SAM_INFO$`AlignedPlus%`) {
  
  #Out superiores
  if(i>out_max){
    muestras_multp_out[n] <- datos_SAM_INFO$Muestra[e]
    n<-n+1
    
    cat("Outlier Superior:\t",datos_SAM_INFO$Muestra[e],"\t",i,"\n")
    
    e<-e+1
    
  }
  #Out inferiores
  if(i<out_min){
    muestras_multp_out[n] <- datos_SAM_INFO$Muestra[e]
    n<-n+1
    
    cat("Outlier Inferior:\t",datos_SAM_INFO$Muestra[e],"\t",i,"\n")
    
    e<-e+1
    
  } else {
    e<-e+1
  }
  
}

#Visualizamos outliers
muestras_multp_out






####################################
#################################################################################################################
#################################### COMPARACIÓN ENTRE LA GENOTIPACION DE NOVO Y LA GENOTIPACION CON GENOMA DE REFERENCIA


#################################

# A NIVEL DE SNPS

#################################

#Declaramos el directorio de trabajo
setwd("C:/Users/José Ramiro Guinea/Desktop/Master/TFM/DATOSSSSS")

#Cargamos los archivos de datos

first_df <- as.data.frame(read.table("C:/Users/José Ramiro Guinea/Desktop/Master/TFM/DATOSSSSS/pau_first_snp",
                                     comment.char = "", sep = "\t", header = T))

stacks_df <- as.data.frame(read.table("C:/Users/José Ramiro Guinea/Desktop/Master/TFM/DATOSSSSS/stacks_first_snp_fil2403.vcf",
                                      comment.char = "", sep = "\t", header = T, skip = 14))



###################### SNPS Y LOCI COINCIDENTES




#Identificamos los ID de los snps que coinciden en ambos archivos VCF
ids <- c()
e <- 1

for (i in 1:length(first_df$ID)) {
  if(first_df$ID[i] %in% stacks_df$ID){
    ids[e] <- first_df$ID[i]
    e <- e+1
    
  }
}

#Identificamos los LOCI que coinciden en ambos archivos VCF
loci <- c()
e <- 1

for (i in 1:length(first_df$ID)) {
  if(first_df$X.CHROM[i] %in% stacks_df$X.CHROM){
    loci[e] <- first_df$X.CHROM[i]
    e <- e+1
  }
}

loci


#Creamos un df nuevo donde presentar los snps detectados por los dos vcf

ref1 <- c()
ref2 <- c()
alt1 <- c()
alt2 <- c()
e <- 1

for (i in ids) {
  ref1[e] <- first_df$REF[first_df$ID==i] 
  ref2[e] <- stacks_df$REF[stacks_df$ID==i]
  alt1[e] <- first_df$ALT[first_df$ID==i]
  alt2[e] <- stacks_df$ALT[stacks_df$ID==i]
  e <- e+1
}

#Documento que alberga los SNPs coincidentes entre Stacks y Freebayes
compare_vcf <- data.frame(cbind(ids,ref1,alt1,ref2,alt2))

#Exportamos el documento
write.table(compare_vcf, file = "COMMON_SITES",
            quote = F, row.names = F,
            sep = "\t")






#Buscamos los SNPS que solo aparecen en alguno de los dos VCF

#FREEBAYES
only_snps_freebayes <- c()
e <- 1

for (i in 1:length(first_df$ID)) {
  if(!first_df$ID[i] %in% ids){
    only_snps_freebayes[e] <- i
    e <- e+1
  }
  
}

#Conjunto de SNPs únicamente de FREEBAYES
only_snps_freebayes <- first_df[only_snps_freebayes,]

#CREAMOS EL DOCUMENTO
write.table(only_snps_freebayes, file = "ONLY_SNPS_FIRST_PAU",
            quote = F, row.names = F,
            sep = "\t")

#STACKS
only_snps_stacks <- c()
e <- 1

for (i in 1:length(stacks_df$ID)) {
  if(!stacks_df$ID[i] %in% ids){
    only_snps_stacks[e] <- i
    e <- e+1
  }
}

#SNPS solamente de STACKS
only_snps_stacks <- stacks_df[only_snps_stacks,]

#Creamos el documento 
write.table(only_snps_stacks, file = "ONLY_SNPS_STACKS",
            quote = F, row.names = F,
            sep = "\t")




################# %MISSING


#FREEBAYES
gt_table <- first_df[,-c(1:9)]
contador <- 0

for (i in 1:length(gt_table)) {
  vect <- gt_table[,i]
  for (e in 1:length(vect)) {
    gt <- substr(vect[e],1,3)
    if(gt=="./."){
      contador <- contador+1
    }
    
  }
}

total <- length(first_df[,10])*243

missing_freebayes <- paste(contador,round(contador*100/total,2),sep =" : ")



#STACKS
gt_table <- stacks_df[,-c(1:9)]
contador <- 0

for (i in 1:length(gt_table)) {
  vect <- gt_table[,i]
  for (e in 1:length(vect)) {
    gt <- substr(vect[e],1,3)
    if(gt=="./."){
      contador <- contador+1
    }
    
  }
}

total <- length(stacks_df[,10])*243

missing_stacks <- paste(contador,round(contador*100/total,2),sep =" : ")



########################### READS Y DEPTH


setwd("C:/Users/José Ramiro Guinea/Desktop/Master/TFM/DATOSSSSS")


#Cargamos los archivos procedentes del programa VCFtools

#FREBAYEES

depth_freebayes <- read.table("results_vcf_comparison/freebayes/freebayes.idepth", #COBERTURA MEDIA POR INDIVIDUO
                              header = T)
site_depth_freebayes <- read.table("results_vcf_comparison/freebayes/freebayes.ldepth", #COBERTURA TOTAL POR LOCUS
                                   header = T)
site_mean_depth_freebayes <- read.table("results_vcf_comparison/freebayes/freebayes.ldepth.mean", #COBERTURA MEDIA POR LOCUS
                                        header = T)


#STACKS

depth_stacks <- read.table("results_vcf_comparison/stacks/stacks.idepth", #COBERTURA MEDIA POR INDIVIDUO
                           header = T)
site_depth_stacks <- read.table("results_vcf_comparison/stacks/stacks.ldepth", #COBERTURA TOTAL POR LOCUS
                                header = T)
site_mean_depth_stacks <- read.table("results_vcf_comparison/stacks/stacks.ldepth.mean", #COBERTURA MEDIA POR LOCUS
                                     header = T)




#### Calculamos estadísticos

# Cobertura media por individuo
mean(depth_freebayes$MEAN_DEPTH)
mean(depth_stacks$MEAN_DEPTH)


range(depth_freebayes$MEAN_DEPTH)
range(depth_stacks$MEAN_DEPTH)


sd(depth_freebayes$MEAN_DEPTH)
sd(depth_stacks$MEAN_DEPTH)



# Cobertura total por locus
mean(site_depth_freebayes$SUM_DEPTH)
mean(site_depth_stacks$SUM_DEPTH)


range(site_depth_freebayes$SUM_DEPTH)
range(site_depth_stacks$SUM_DEPTH)

sd(site_depth_freebayes$SUM_DEPTH)
sd(site_depth_stacks$SUM_DEPTH)


# Cobertura media por locus
mean(site_mean_depth_freebayes$MEAN_DEPTH)
mean(site_mean_depth_stacks$MEAN_DEPTH)


range(site_mean_depth_freebayes$MEAN_DEPTH)
range(site_mean_depth_stacks$MEAN_DEPTH)

sd(site_mean_depth_freebayes$MEAN_DEPTH)
sd(site_mean_depth_stacks$MEAN_DEPTH)


#Gráficos cobertura

#FREEBAYES
ggplot(depth_freebayes, aes(x = MEAN_DEPTH))+
  geom_histogram(colour = "white", fill = 6,
                 bins = 150)+
  ggtitle("Histograma de Cobertura: Freebayes")+
  labs(x = "Cobertura media", y = "Frecuencia")+
  xlim(0,115)+ylim(0,18)

#STACKS
ggplot(depth_stacks, aes(x = MEAN_DEPTH))+
  geom_histogram(colour = "white", fill = 4,
                 bins = 150)+
  ggtitle("Histograma de Cobertura: Stacks")+
  labs(x = "Cobertura media", y = "Frecuencia")+
  xlim(0,115)+ylim(0,18)


# Valores IQR
IQR(depth_freebayes$MEAN_DEPTH)
IQR(depth_stacks$MEAN_DEPTH)

















####################################
#################################################################################################################
#################################### CAMBIO DE FORMATO DE LOS ARCHIVOS VCF PARA INTRODUCIRLOS EN PLINK2

##############

#PAU_FIRST_SNP

##############

#Cambiamos el nombre de los LOCI en la primera columna del dataframe

#Duplicamos el df con otro nombre para no modificar el anterior.
plink2_df <- datos_finales

#Cambiamos nombres
for (i in 1:length(plink2_df$`#CHROM`)) {
  plink2_df$`#CHROM`[i] <- paste("CLocus_",plink2_df$`#CHROM`[i],sep = "")
  
}


#Creamos el nuevo archivo vcf
write.table(plink2_df,file = "pau_first_plink2",sep = "\t",
            row.names = FALSE, quote = F)



##############

#STACKS_FIRST_SNP_FIL2403

##############


#Cargamos los datos de stacks
stacks_df <- as.data.frame(read.table("stacks_first_snp_fil2403.vcf",comment.char = "",header = TRUE,
                                      skip = 14))

#Cambiamos nombres de LOCI
for (i in 1:length(stacks_df$X.CHROM)) {
  stacks_df$X.CHROM[i] <- paste("CLocus_",i,sep = "")
}

#Cambiamos el nombre de la variable X.CHROM por #CHROM
colnames(stacks_df)[1] <- "#CHROM"

#Creamos el nuevo archivo .vcf
write.table(stacks_df,file = "stacks_plink2",sep = "\t",
            row.names = FALSE, quote = F)






























####################################
#################################################################################################################
#################################### VISUALIZAMOS LAS PCA

library(ggplot2)
library(tidyr)

setwd("C:/Users/José Ramiro Guinea/Desktop/Master/TFM/DATOSSSSS")








####################################

# FREEBAYES

####################################


################ FIRST SNP 

#cargamos archivos


first_pau_pca=read.table("plink2_first_pau/paau.eigenvec")

names(first_pau_pca) <- c('FI','ID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')

plot(first_pau_pca$PC1,first_pau_pca$PC2)

# plot pca
ggplot(first_pau_pca, aes(PC1, PC2, col = ID)) + geom_point(size = 3) + theme(legend.position="none")

pdf('plink2_first_pau/first_snps_pau.pdf')
ggplot(first_pau_pca, aes(PC1, PC2, col = ID)) + geom_point(size = 3) + theme(legend.position="none")
dev.off()




# BUscamos outliers
first_pau_pca$ID[first_pau_pca$PC1==sort(first_pau_pca$PC1,decreasing = F)[1]]

#Creamos un nuevo df sin outliers
pca_pau_first_no_outliers <- first_pau_pca[first_pau_pca$ID!=first_pau_pca$ID[first_pau_pca$PC1==sort(first_pau_pca$PC1,decreasing = F)[1]], ]

#Visualizamos la PCA SIN OUTLIERS
ggplot(prueba_first_pau, aes(PC1, PC2, col = ID)) + geom_point(size = 3) + theme(legend.position="none")








####################################################

# STACKS

####################################################


staacks_pca=read.table("plink2_stacks/staacks.eigenvec")

names(staacks_pca) <- c('FI','ID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')

plot(staacks_pca$PC1,staacks_pca$PC2)

# plot pca
ggplot(staacks_pca, aes(PC1, PC2, col = ID)) + geom_point(size = 3) + theme(legend.position="none")

pdf('plink2_stacks/snps_staacks.pdf')
ggplot(staacks_pca, aes(PC1, PC2, col = ID)) + geom_point(size = 3) + theme(legend.position="none")
dev.off()



#Buscamos los outliers
staacks_pca$ID[staacks_pca$PC2==sort(staacks_pca$PC2)[1]]
staacks_pca$ID[staacks_pca$PC2==sort(staacks_pca$PC2)[2]]

#Creamos un nuevo df sin los outliers
pca_staacks_no_outliers <- staacks_pca[!(staacks_pca$ID==staacks_pca$ID[staacks_pca$PC2==sort(staacks_pca$PC2)[1]] |
                         staacks_pca$ID==staacks_pca$ID[staacks_pca$PC2==sort(staacks_pca$PC2)[2]]),]

#Visualizamos la PCA SIN OUTLIERS
ggplot(prueba2, aes(PC1, PC2, col = ID)) + geom_point(size = 3) #+ theme(legend.position="none")










#################################################################

# COMPARAMOS GRUPOS, NO INDIVIDUOS

#################################################################

############################# FREEBAYES (HAY OUTLIERS)

#DEFINIMOS LOS GRUPOS
id_actual <- c()
ordinal <- c()

for (i in 1:length(first_pau_pca$ID)) {
  a <- gsub("0|1|2|3|4|5|6|7|8|9", "",first_pau_pca$ID[i])
  id_actual[i] <- a
  ordinal[i] <- i

}

#Creamos un data.frame con los datos de cada grupo
df1 <- data.frame(id_actual,ordinal)

#Miramos el nombre de cada grupo
tags <- unique(df1$id_actual)

#Creamos los grupos


A <- df1$ordinal[df1$id_actual=="A"]
length(A)
AK <- df1$ordinal[df1$id_actual=="AK"]
length(AK)
BEL <- df1$ordinal[df1$id_actual=="BEL"]
length(BEL)
DA <- df1$ordinal[df1$id_actual=="DA"]
length(DA)
IS <- df1$ordinal[df1$id_actual=="IS"]
length(IS)
KY <- df1$ordinal[df1$id_actual=="KY"]
length(KY)
LBN <- df1$ordinal[df1$id_actual=="LBN"]
length(LBN)
LI <- df1$ordinal[df1$id_actual=="LI"]
length(LI)
ME <- df1$ordinal[df1$id_actual=="ME"]
length(ME)
Mes <- df1$ordinal[df1$id_actual=="Mes"]
length(Mes)
RE <- df1$ordinal[df1$id_actual=="RE"]
length(RE)
Zak <- df1$ordinal[df1$id_actual=="Zak"]
length(Zak)

#Añadimos una columna con el nombre de cada grupo
first_pau_pca$GROUP <- c(rep("A",25), #A
                         rep("AK",25), #AK
                         rep("BEL",24), #BEL
                         rep("DA",24), #DA
                         rep("IS",21), #IS
                         rep("KY",25), #KY
                         rep("LBN",19), #LBN
                         rep("LI",23), #LI
                         rep("ME",18), #ME
                         rep("RE",23), #RE
                         rep("Zak",16)) #Zak


####################################### FREEBAYES (SIN OUTLIERS)
#DEFINIMOS LOS GRUPOS
id_actual <- c()
ordinal <- c()

for (i in 1:length(pca_pau_first_no_outliers$ID)) {
  a <- gsub("0|1|2|3|4|5|6|7|8|9", "",pca_pau_first_no_outliers$ID[i])
  id_actual[i] <- a
  ordinal[i] <- i
  
}

#Creamos un data.frame con los datos de cada grupo
df2<- data.frame(id_actual,ordinal)

#Miramos el nombre de cada grupo
tags <- unique(df2$id_actual)

#Creamos los grupos


A <- df2$ordinal[df2$id_actual=="A"]
length(A)
AK <- df2$ordinal[df2$id_actual=="AK"]
length(AK)
BEL <- df2$ordinal[df2$id_actual=="BEL"]
length(BEL)
DA <- df2$ordinal[df2$id_actual=="DA"]
length(DA)
IS <- df2$ordinal[df2$id_actual=="IS"]
length(IS)
KY <- df2$ordinal[df2$id_actual=="KY"]
length(KY)
LBN <- df2$ordinal[df2$id_actual=="LBN"]
length(LBN)
LI <- df2$ordinal[df2$id_actual=="LI"]
length(LI)
ME <- df2$ordinal[df2$id_actual=="ME"]
length(ME)
Mes <- df2$ordinal[df2$id_actual=="Mes"]
length(Mes)
RE <- df2$ordinal[df2$id_actual=="RE"]
length(RE)
Zak <- df2$ordinal[df2$id_actual=="Zak"]
length(Zak)

#Añadimos una columna con el nombre de cada grupo
pca_pau_first_no_outliers$GROUP <- c(rep("A",24), #A
                                   rep("AK",25), #AK
                                   rep("BEL",24), #BEL
                                   rep("DA",24), #DA
                                   rep("IS",21), #IS
                                   rep("KY",25), #KY
                                   rep("LBN",19), #LBN
                                   rep("LI",23), #LI
                                   rep("ME",18), #ME
                                   rep("RE",23), #RE
                                   rep("Zak",16)) #Zak








###################################### STACKS (CON OUTLIERS)

#DEFINIMOS LOS GRUPOS
id_actual <- c()
ordinal <- c()

for (i in 1:length(staacks_pca$ID)) {
  a <- gsub("0|1|2|3|4|5|6|7|8|9", "",staacks_pca$ID[i])
  id_actual[i] <- a
  ordinal[i] <- i
  
}

#Creamos un data.frame con los datos de cada grupo
df3<- data.frame(id_actual,ordinal)

#Miramos el nombre de cada grupo
tags <- unique(df3$id_actual)

#Creamos los grupos


A <- df3$ordinal[df3$id_actual=="A"]
length(A)
AK <- df3$ordinal[df3$id_actual=="AK"]
length(AK)
BEL <- df3$ordinal[df3$id_actual=="BEL"]
length(BEL)
DA <- df3$ordinal[df3$id_actual=="DA"]
length(DA)
IS <- df3$ordinal[df3$id_actual=="IS"]
length(IS)
KY <- df3$ordinal[df3$id_actual=="KY"]
length(KY)
LBN <- df3$ordinal[df3$id_actual=="LBN"]
length(LBN)
LI <- df3$ordinal[df3$id_actual=="LI"]
length(LI)
ME <- df3$ordinal[df3$id_actual=="ME"]
length(ME)
Mes <- df3$ordinal[df3$id_actual=="Mes"]
length(Mes)
RE <- df3$ordinal[df3$id_actual=="RE"]
length(RE)
ZA <- df3$ordinal[df3$id_actual=="ZA"]
length(ZA)

#Añadimos una columna con el nombre de cada grupo
staacks_pca$GROUP <- c(rep("AK",25), #A
                       rep("A",25), #AK
                       rep("BEL",24), #BEL
                       rep("DA",24), #DA
                       rep("IS",21), #IS
                       rep("KY",25), #KY
                       rep("LBN",19), #LBN
                       rep("LI",23), #LI
                       rep("ME",18), #ME
                       rep("RE",23), #RE
                       rep("ZA",16)) #ZA



###################################### DATOS STACKS (SIN OUTLIERS)

#DEFINIMOS LOS GRUPOS
id_actual <- c()
ordinal <- c()

for (i in 1:length(pca_staacks_no_outliers$ID)) {
  a <- gsub("0|1|2|3|4|5|6|7|8|9", "",pca_staacks_no_outliers$ID[i])
  id_actual[i] <- a
  ordinal[i] <- i
  
}

#Creamos un data.frame con los datos de cada grupo
df4<- data.frame(id_actual,ordinal)

#Miramos el nombre de cada grupo
tags <- unique(df4$id_actual)

#Creamos los grupos


A <- df4$ordinal[df4$id_actual=="A"]
length(A)
AK <- df4$ordinal[df4$id_actual=="AK"]
length(AK)
BEL <- df4$ordinal[df4$id_actual=="BEL"]
length(BEL)
DA <- df4$ordinal[df4$id_actual=="DA"]
length(DA)
IS <- df4$ordinal[df4$id_actual=="IS"]
length(IS)
KY <- df4$ordinal[df4$id_actual=="KY"]
length(KY)
LBN <- df4$ordinal[df4$id_actual=="LBN"]
length(LBN)
LI <- df4$ordinal[df4$id_actual=="LI"]
length(LI)
ME <- df4$ordinal[df4$id_actual=="ME"]
length(ME)
Mes <- df4$ordinal[df4$id_actual=="Mes"]
length(Mes)
RE <- df4$ordinal[df4$id_actual=="RE"]
length(RE)
ZA <- df4$ordinal[df4$id_actual=="ZA"]
length(ZA)

#Añadimos una columna con el nombre de cada grupo
pca_staacks_no_outliers$GROUP <- c(rep("AK",25), #A
                       rep("A",25), #AK
                       rep("BEL",24), #BEL
                       rep("DA",24), #DA
                       rep("IS",21), #IS
                       rep("KY",25), #KY
                       rep("LBN",19), #LBN
                       rep("LI",23), #LI
                       rep("ME",18), #ME
                       rep("RE",21), #RE
                       rep("ZA",16)) #ZA














################################

#PCAS POR GRUPOS

################################

#PCAS completas con colores por grupos (CON OUTLIERS)
ggplot(all_pau_pca, aes(PC1, PC2, col = GROUP)) + geom_point(size = 3) + ggtitle("ALL_SNPs")  #ALL_SNPS
ggplot(first_pau_pca, aes(PC1, PC2, col = GROUP)) + geom_point(size = 3) + ggtitle("FIRST_SNPs") #FIRST_SNPS
ggplot(staacks_pca, aes(PC1, PC2, col = GROUP)) + geom_point(size = 3) + ggtitle("STACKS") #STACKS


#PCAS completas con colores por grupos (SIN OUTLIERS)
ggplot(pca_pau_all_no_outliers, aes(PC1, PC2, col = GROUP)) + geom_point(size = 3) + ggtitle("ALL_SNPs") #ALL_SNPS
ggplot(pca_pau_first_no_outliers, aes(PC1, PC2, col = GROUP)) + geom_point(size = 3) + ggtitle("FIRST_SNPs") #FIRST_SNPS
ggplot(pca_staacks_no_outliers, aes(PC1, PC2, col = GROUP)) + geom_point(size = 3) + ggtitle("STACKS") #STACKS




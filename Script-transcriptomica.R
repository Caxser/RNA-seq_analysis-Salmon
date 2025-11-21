library(here)
library(tximeta)
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(vsn)

library(glmGamPoi)
library(IHW)


#######
library(DESeq2); library(pheatmap)

vsd <- vst(dds, blind=TRUE)
mat <- assay(vsd)

## Heatmap de distancias
d <- dist(t(mat))
pheatmap(as.matrix(d), annotation_col=as.data.frame(colData(dds)))

## PCA (marca outliers visuales)
plotPCA(vsd, intgroup="treatment") 
##
lib <- colSums(counts(dds))
genes_expres <- colSums(counts(dds) > 0)

boxplot(lib, main="Library size"); boxplot(genes_expres, main="# genes expresados")

suspect_lib <- names(lib[ lib < (median(lib) - 1.5*IQR(lib)) ])   # muy bajas
suspect_gen  <- names(genes_expres[ genes_expres < (median(genes_expres) - 1.5*IQR(genes_expres)) ])

# Matriz de Cooks (puede ser grande). Si no existe, vuelve a correr DESeq(dds)
cooks <- assays(dds)[["cooks"]]

# Proporción de genes por muestra con Cook’s > 5 (o el que prefieras)
prop_cooks_hi <- colMeans(cooks > 5, na.rm=TRUE)

# Muestras que concentran muchos genes con Cooks alto
suspect_cooks <- names(prop_cooks_hi[prop_cooks_hi > quantile(prop_cooks_hi, 0.95)])

sort(prop_cooks_hi, decreasing=TRUE)[1:5]
# Genes expresados por muestra
genes_expres <- colSums(counts(dds) > 0)
sort(genes_expres)

# Muestras con menos del Q1 - 1.5*IQR
Q1 <- quantile(genes_expres, 0.25)
IQRv <- IQR(genes_expres)
low_cutoff <- Q1 - 1.5*IQRv
names(genes_expres[genes_expres < low_cutoff])

plotPCA(vst(dds, blind=TRUE), intgroup=c("treatment", "names"))
suspects <- unique(c(suspect_cooks, suspect_lib, suspect_gen))
suspects
PCA <- pca(assay(rld), scale = T, metadata = colData(rld)) 
plotPCA(vst(dds, blind=TRUE), intgroup = "treatment")+
  geom_label_repel(aes(label = colnames(assay(vst(dds, blind=TRUE)))), 
                   segment.color = "grey50", 
                   box.padding = 0.35, 
                   point.padding = 0.5)
######
#Calculo de poder
#if (!require("BiocManager")) install.packages("BiocManager")
#BiocManager::install("RNASeqPower")
library(RNASeqPower)
gene <- "ENSG00000164867"  # es de NOS3 o "ESR1"
cts <- counts(dds)
gcounts <- cts[rownames(cts)==gene, ]
depth_gene <- median(gcounts)
# Gen Ensembl base (sin versión)
gene <- "ENSG00000164867"  # NOS3  (ESR1 sería ENSG00000091831)

# Extraer counts
cts <- counts(dds)

# Crear vector con IDs base (sin versión)
rownames_base <- sub("\\..*", "", rownames(cts))

# Buscar coincidencia
idx <- which(rownames_base == gene)

# Obtener cuentas del gen
if(length(idx) == 1) {
  gcounts <- cts[idx, ]
  depth_gene <- median(as.numeric(gcounts))
  depth_gene
} else {
  message("No se encontró una coincidencia exacta para el ID Ensembl base.")
}

rnapower(depth=depth_gene ,n=9, n2=10, cv=0.85, effect=2 , alpha= 0.05)

rnapower(depth=depth_gene , cv=0.54, effect=2 , alpha= 0.05, power= 0.8)


colSums(counts(dds)) / 1e6 


bm <- rowMeans(counts(dds, normalized = TRUE))
high <- bm >= quantile(bm, 0.5, na.rm=TRUE) 
disp <- dispersions(dds)           # dispersión NB por gen
bcv  <- sqrt(disp)                 # BCV ~ CV para rnapower
cv_global <- mean(bcv[bm], na.rm=TRUE)   # CV representativo (robusto)
cv_global

##########
#Intento para variables ortogonales
library(sva)

vsd <- vst(dds, blind=TRUE)
mat <- assay(vsd)        # genes x muestras
meta <- as.data.frame(colData(dds))

#generar modelos
# ejemplo: 'condition' es PE vs DMG
mod  <- model.matrix(~ base + treatment, data=meta)
mod0 <- model.matrix(~ 1, data=meta)

svobj <- svaseq(dat = as.matrix(mat), mod = mod, mod0 = mod0)
head(svobj$sv)
coldata$sv1 <- svobj$sv[,1]
coldata$sv2 <- svobj$sv[,2]
coldata$sv3 <- svobj$sv[,3]
coldata$sv4 <- svobj$sv[,4]
coldata$sv5 <- svobj$sv[,5]

#No reduce suficiente, se descarta

#___
#intento con IHW
#IHW es un método para ajustar p-valores en pruebas múltiples que incorpora 
# una covariable auxiliar además de los p-valores. La idea es la siguiente: 
# en estudios grandes de expresión diferencial, cada gen tiene distinta 
# “calidad” para detectar cambio (por ejemplo, diferente número de lecturas, 
# variabilidad, etc.). IHW aprovecha esa variabilidad para priorizar hipótesis 
# con más potencia o mayor probabilidad de cambio real.

# Necesario demostrar la indeéndencia (no correlacion de la variable con padj)
names(res)
#La variable a usar es baseMean
covariate <- res$baseMean
pval <- res$pvalue

keep <- !is.na(pval) & !is.na(covariate)
pval <- pval[keep]
covariate <- covariate[keep]
idx_null <- which(pval > 0.5)
cor.test(pval[idx_null], covariate[idx_null], method="spearman")

library(ggplot2)
df <- data.frame(pval=pval, baseMean=covariate)

ggplot(df, aes(x=log10(baseMean+1), y=pval)) +
  geom_point(alpha=0.3, size=0.7) +
  geom_smooth(method="loess", color="red", se=FALSE) +
  theme_minimal() +
  labs(x="log10(baseMean + 1)", y="p-value",
       title="Relación entre p-value y baseMean")

###ahora uso de IHW
library(DESeq2)
library(IHW)

# Supongamos que ya tienes tu objeto dds (DESeqDataSet) generado con tximeta
dds <- DESeq(dds)  
res <- results(dds)  # tabla de resultados: contiene baseMean, log2FoldChange, pvalue, padj, etc.

# Elige la covariable: por ejemplo baseMean (media normalizada de cuenta del gen)
# Verifica que no haya correlación entre baseMean y pvalue bajo hipótesis nula (puedes usar hist para pvalues>0.5, etc.)

# Aplica IHW
ihwRes <- ihw(pvalue ~ baseMean, data = as.data.frame(res), alpha = 0.05)

# Extrae los p-valores ajustados con IHW
res$ihw_padj <- adj_pvalues(ihwRes)

# Ver cuántos genes rechazamos con IHW
sum(res$ihw_padj <= 0.05, na.rm = TRUE)

# Comparar con método BH estándar
sum(res$padj <= 0.05, na.rm = TRUE)

## a partir de aquí, solo es filtrar por ihw_padj y guardar el csv
######
#Datos en wsl
# ruta: \\wsl.localhost\Ubuntu-24.04\home\cesar\trabajos-T\conteos-salmon\SRR12038075
#
setwd("d:/Proyecto-maestria/R_carpet/PE")
#setwd("d:/Proyecto-maestria/R_carpet/GDM")

#Necesitamos los metadatos para crear columnas names y dir
coldata <- as.data.frame(read.csv("2mod-fusion-placenta_Sraruntable.csv"))
#en este caso, la columna name fue creada de forma manual

##Posiblemente necesario saber la dirección del fragmento (para añadir covariable)
#despues de salmon, corre el script en R "stranded.r"
#manten la variable "tab" en el entorno para lo siguiente

#coldata <-left_join(coldata, tab, by = c("Run" = "sample"))

#Añadimos los nombres a las filas
#Importante, las columna names y treatment debe existir (son los nombres o ids de las muestras)
#recomendación, en names= añadir también qué tipo de tratamiento son
rownames(coldata) <- coldata$names

#ahora se buscan las direcciones de las cuantificaciones
#ubicamos la direccion de cuantificacion
dir <- file.path(here("d:/Proyecto-maestria/conteos-salmon"))

#revisamos si existen los archivos
list.files(dir)

#Añadimos una nueva columna con la direccion de los archivos
# se le agrega "_R1.clean.fastq.gz" por error en la creación de carpetas
####Error subsanado, para posteriores scripts solo colocar coldata$Run

coldata$files <- file.path(dir,paste0(coldata$Run), "quant.sf")

#verificamos que existan los archivos
data.frame(coldata$Run, file.exists(coldata$files))
# El resultado de la columna debe ser TRUE, si no, se debe revisar que paso

#-----------------------------
#Comenzamos con la importacion de los objetos
se <- tximeta(coldata)
#Hacemos que se promedie solo por gen, ya que no nos interesan transcritos
se <- summarizeToGene(se)
se

#Uso de sumarized experiments
#basicamente son "cajones" que contienen informacion
#colData, donde tienen los metadatos del experimento
#rowranges, hace referencia a las coordenadas de cada transcrito
#assay, almacena informacion de las cuentas para cada transcrito en tres niveles


#tabla de metadatos creada para importar las cuentas
colData(se)
#es un dataframe
class(colData(se))

#coordenadas de transcritos
rowRanges(se)

#cuentas de cada transcrito dividido en tres
assayNames(se)
#Counts= cuentas crudas
#Abundance= cuentas normalizadas por tpm
#Lenght= longitud

#Para acceder a counts (cuentas crudas)
#el head es solo para no saturar la vista (los 5 primerso genes)
head(assay(se),5)

#Para acceder a matriz cuentas normalizadas tpm
#se usa @ por ser una lista de listas
head(se@assays@data$abundance, 5)

#Es importante verificar que Rownames en coldata sea
# igual a los colnames en assay(se) - columnas de matriz de cuentas
#crudas. IMPORTANTE PARA IMPORTAR DATOS DESEQ2
row.names(colData(se))==colnames(assay(se))


#NORMALIZACION
#Criterios para normalizar
#*Tamaño de libreria (profundidad): cada muestra puede tener
#diferente tamaño de libreria
#*Tamaño del gen
#*Composicion de RNA
#Tipos:
#CPM, cuentas por millon
#TPM, transcritos por millon de lecturas
#Factor de normalizacion por DESeq2
#RPKM/FPKM:lecturas/fragmentos por kilobase 
#de exón por millón de lecturas/fragmentos mapeados
#EdgeR

#DESeq2
#Modelado lineal generalizado (binomial negativa)
#Necesitamos definir los grupos (control y tratamiento)
se$treatment <- factor(se$treatment)
se$base <- factor(se$base)
#se$call <- factor(se$call)
#asegurar que la referencia sea el control (untreated)
se$treatment <- relevel(se$treatment, ref = "Control")


#Checar la cantidad de muestras
table(se$treatment)
table(se$base)
#table(se$call)
View(as.data.frame(colData(se)))


#Generar objeto de DESeq
#En este punto se pueden añadir más covariables despues de
#treatment (recordar, es un modelo lineal)
dds <- DESeqDataSet(se, design= ~ treatment)

#######
#Ajuste de batch por comBAT
library(sva)

cts <- counts(dds)
batch <- se$base          # variable que indica el estudio o lote
group <- se$treatment     # variable biológica: Control / PE / GDM

cts_corr <- ComBat_seq(as.matrix(cts),
                       batch = batch,
                       group = group)

# Crear un nuevo DESeqDataSet con los counts ya corregidos
dds <- DESeqDataSetFromMatrix(countData = cts_corr,
                              colData = colData(se),
                              design = ~ treatment)

#######

#prefiltro para quitar genes con cuentas bajas
#primero se quitan aquellos que tengan valores en más de 3 genes
#
keep <- rowSums(counts(dds) >=3) >=6
dds <- dds[keep, ]

#Funcion para normalizar los datos y realiza analisis de
#expresión diferencial
#dds <- replaceOutliersWithTrimmedMean(dds)
dds <- estimateSizeFactors(dds, type= "poscount")
#Revision de sizefactor

sizeFactors(dds)



dds <- DESeq(dds)
vsd_before <- vst(DESeqDataSet(se, design=~ treatment), blind=TRUE)
vsd_after  <- vst(dds, blind=TRUE)   # este dds ya corregido


plotPCA(vsd_before, intgroup=c("base","treatment"))
plotPCA(vsd_after, intgroup=c("base","treatment"))

resultsNames(dds)
#Para ver los resultados se usa "results"
res_PE <- results(dds,contrast = c("treatment","Preeclampsia","Control"))
res_GDM <- results(dds,contrast = c("treatment","GDM","Control"))
#res

#Resumen de resultados
summary(res_PE)
summary(res_GDM)

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
#####
dds_param <- DESeq(dds, fitType="parametric")
dds_glm   <- DESeq(dds, fitType="glmGamPoi")

# Comparar log2FoldChange
plot(results(dds_param)$log2FoldChange,
     results(dds_glm)$log2FoldChange,
     pch=16, cex=0.6,
     xlab="parametric LFC", ylab="glmGamPoi LFC")
abline(0,1,col="red")
res_param <- results(dds_param)$log2FoldChange
res_glm   <- results(dds_glm)$log2FoldChange
cor(res_param, res_glm, use="complete.obs", method="pearson")

df_diff <- data.frame(
  baseMean = results(dds_param)$baseMean,
  LFC_param = res_param,
  LFC_glm = res_glm,
  diff = res_param - res_glm
)
df_diff <- df_diff[order(abs(df_diff$diff), decreasing=TRUE), ]
head(df_diff, 20)
######
#_:_________________________________________________________________________
#La funcion results acepta tambien filtros
#por ejemplo filtrar aquellos con lfc >1
#prueba <- results(dds,lfcThreshold = 1)
#filtro con pvalue adj menor a 0.1
#table(prueba$padj< 0.05)

#Ordenar por significancia
resOrdered_GDM <- res_GDM[order(res_GDM$padj), ]
resOrdered_PE <- res_PE[order(res_PE$padj), ]

#Genes con 0.05 en padj y con lfc mayor a 1 (absoluto)
resSig_GDM <- subset(resOrdered_GDM, padj < 0.05 & abs(log2FoldChange) > 1)
resSig_PE <- subset(resOrdered_PE, padj < 0.05 & abs(log2FoldChange) > 1)


#head(resSig)

#Pasarlo a un csv
write.csv(as.data.frame(resSig_GDM), file = "final-GDM_genes_diferencialmente_expresados.csv")
write.csv(as.data.frame(resSig_PE), file= "final-PE_genes_diferencialmente_expresados.csv")
write.csv(as.data.frame(resOrdered_GDM), file= "final-full-GDM_genes_diferencialmente_expresados.csv")
write.csv(as.data.frame(resOrdered_PE), file= "final-full-PE_genes_diferencialmente_expresados.csv")

#####
##que tanto se ajusta el modelo?
plotDispEsts(dds)

#cook
boxplot(log10(assays(dds)[["cooks"]]), las=2, range=0)

#Plot de 0
hist(res_GDM$pvalue, breaks=50, col="gray")
hist(res_PE$pvalue, breaks=50, col="gray")

#------------------------------------
#Los nombres estan en formato Ensembl Transcript IDs
#para pasarlos a simbolos de genes
library(biomaRt)
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

#los ids estan como rownames en resSig
enst_ids <- rownames(resSig_PE)
#Limpiarlos
# Extraer solo el ID base sin versión
enst_ids_clean <- sub("\\..*$", "", enst_ids)


#Se procede a su conversión
#ensmbl transcript id -> ensmbl gene id -> gen ID
conversion <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = enst_ids_clean,
  mart = ensembl
)



#Se procede a unir los resultados
resSig_df <- as.data.frame(resSig_PE)
# Agrega columna con ID original (su nombre como fila)
resSig_df$gene_id_raw <- rownames(resSig_df)

# Limpia el ID (sin versión)
resSig_df$gene_id_clean <- sub("\\..*$", "", resSig_df$gene_id_raw)

# Combinar por ID de transcrito
resSig_annotated <- merge(resSig_df, conversion,
                          by.x = "gene_id_clean",
                          by.y = "ensembl_gene_id",
                          all.x = TRUE)
#Guardar
write.csv(resSig_annotated, "final-PE_genes_diferencialmente_expresados_annotados.csv")

#-----------------------
#Transformaciones alternativas (normalizaciones)

#Grafica de cuentas crudas
meanSdPlot(counts(dds), ranks = F)
# Gráfica de cuentas en escala log2
meanSdPlot(log2(counts(dds) +1), ranks = F)

# variance stabilizing transformation (VST), (Anders and Huber 2010)
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

# regularized-logarithm transformation (rlog), (Love, Huber, and Anders 2014)
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

# Normalización con el factor de normalizacion 
dds <- estimateSizeFactors(dds, type="poscounts")

# Juntar los datos de las tres normalizaciones
df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog")
  )

# Renombrar columnas
colnames(df)[1:2] <- c("x", "y")  

# Nombre de las graficas
lvls <- c("log2(x + 1)", "vst", "rlog")

# Agrupar los tres tipos de normalizacion en grupos como factores
df$transformation <- factor(df$transformation, levels=lvls)

# Plotear los datos
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  


#--------------------------
#Seccion de graficos
#PCA, analisis de componentes principales
library(PCAtools)

#Se recomienda usar rlog (rld)
#para menos de 30 muestras
#para mayor a 30 se usa vst
vsd <- vst(dds, blind=FALSE)
#Comando proviene de DESEq2
#Este solo toma en cuenta los top 500 genes
plotPCA(vsd, intgroup= "base")

plotPCA(rld, intgroup = "treatment")

#esta es la grafica sencilla
# a continuacion se pone un ejemplo
# mas complejo
## Crear un objeto que contenga los datos del PCA
#comando interno de PCAtools
#Este toma todos los genes
PCA <- pca(assay(rld), scale = T, metadata = colData(rld)) 
## Graficar en 2D los resultados
biplot(PCA, colby = "treatment")
#Otro ejemplo
plotPCA(rld, intgroup = "treatment")+
  geom_label_repel(aes(label = colnames(assay(rld))), 
                   segment.color = "grey50", 
                   box.padding = 0.35, 
                   point.padding = 0.5)

## Generar un screeplot para visualizar la varianza asociada a cada componente
screeplot(PCA)

#MAplot
group <- coldata$treatment
dds$group <- group
##Es importante que recuerden que la hipótesis nula que se probó fue
##"El lfc del gen n es igual a 0" por lo tanto los genes coloreados son...
plotMA(res)

#Glima
#Genera un plot interactivo
library(Glimma)
glimmaMA(dds)


##############
#Pruebas para corregir el batch effect
library(limma)
library(ggfortify)
# Extrae la matriz de expresión normalizada
expr <- assay(vsd)

# Corrige por el efecto de "base"
expr_corrected <- removeBatchEffect(expr, batch = vsd$base, design = model.matrix(~ vsd$treatment))

# Ahora puedes volver a hacer el PCA
pca <- prcomp(t(expr_corrected))
autoplot(pca, data = as.data.frame(colData(vsd)), colour = 'base')
ggplot(pca_df, aes(PC1, PC2, color = base)) +
  geom_point(size = 3) + theme_bw()


DESeq2::plotMA(res_GDM)
DESeq2::plotMA(res_PE)
library(apeglm)
resLFC <- lfcShrink(dds, coef="treatment_GDM_vs_Control", type="apeglm")
DESeq2::plotMA(resLFC, ylim=c(-4,4))

######
#REVISIÓN PARA METODO DE NORMALIZACION
#Justificar uso postcounts
#Cuantos ceros por muestra
counts_mat <- counts(dds)
zero_prop <- sum(counts_mat == 0) / length(counts_mat)
zero_prop

#distribución de ceros por muestra
zero_by_sample <- colSums(counts_mat == 0) / nrow(counts_mat)
barplot(zero_by_sample, las=2,
        ylab = "Proporción de genes con cuenta cero",
        xlab = "Muestra")

#distribución de cuentas
hist(log10(rowSums(counts_mat) + 1), breaks=100,
     main="Distribución de sumas por gen", xlab="log10(suma de cuentas +1)")

#lecturas por muestra
barplot(colSums(counts(dds)),
        las = 2,
        ylab = "Total de lecturas asignadas",
        xlab = "Muestra")


#######
# Extraer la matriz de longitudes efectivas
len <- assays(se)$length

# Calcular la longitud promedio de cada transcrito entre todas las muestras
mean_len_per_tx <- rowMeans(len, na.rm = TRUE)

# Estadísticas generales
summary(mean_len_per_tx)

# Histograma de longitudes
hist(mean_len_per_tx,
     breaks = 80,
     col = "skyblue",
     main = "Distribución de longitudes promedio de transcritos",
     xlab = "Longitud efectiva promedio (nt)")


len <- assays(se)$length
mean_len_per_tx <- rowMeans(len, na.rm = TRUE)

summary(mean_len_per_tx)
sum(is.na(mean_len_per_tx))
sum(mean_len_per_tx < 200)
sum(mean_len_per_tx > 5000)

#######
###Profundidad
# dds creado a partir de tximeta/tximport
lib_depth_frag <- colSums(counts(dds))           # fragmentos asignados
lib_depth_M    <- lib_depth_frag / 1e6           # en millones
muestras_profundidad <-data.frame(sample=colnames(dds), Mfrags=round(lib_depth_M,2))
barplot(lib_depth_M, las=2, ylab="Millones de fragmentos asignados",
        main="Profundidad por muestra (DESeq2)")
sort(muestras_profundidad$Mfrags, decreasing = F)

muestras_profundidad <- muestras_profundidad[order(muestras_profundidad$Mfrags),]
muestras_profundidad

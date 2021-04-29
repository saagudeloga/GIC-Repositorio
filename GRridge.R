#=================================================================
# Codigo R - Articulo: Análisis estadístico para 
# la identificación de miRNAs de interés a partir
# de datos miRNAseq usando regresión logística ridge regularizada
# con inclusión de co-datos
#=================================================================

# Cargar librerias
require(edgeR)
require(limma)
require(DESeq2)
require(GRridge)
require(ggplot2)
require(factoextra)
require(glmnet)
require(ISLR)

# Lectura de datos
countdata_GIC <- read.table("miRNA_piRNA_HPV_UdeA.txt", row.names=1, header=T, check.names = F)
countdata_GIC <- as.data.frame(countdata_GIC)
colnames(countdata_GIC) <- NULL
calldata_GIC <- read.table("colnames.txt", header=T, row.names=1, sep="\t")
calldata_GIC$Grupo <- ifelse(calldata_GIC$Grupo == "Group1", "Normal", "Caso")
calldata_GIC$Grupo <- as.factor(calldata_GIC$Grupo)
table(calldata_GIC$Grupo)

# Profundidad genes
Profundidad_GIC <- apply(countdata_GIC, MARGIN = 1, sum)
max(Profundidad_GIC) # Maxima profundiad
min(Profundidad_GIC) # Minima profundiad
dim(countdata_GIC) # Numero de genes
table(calldata_GIC$Grupo) # Division muestras

# Densidad profundidad
#plot(density(rowSums(countdata_GIC > 1)), lwd = 3)

# Filtracion y eliminacion de genes de baja expresion
keep_GIC <- rowSums(cpm(countdata_GIC) > 1) >= 20
countdata_GIC <- countdata_GIC[keep_GIC,]
countdata_GIC <- countdata_GIC[rownames(countdata_GIC) != "hsa-miR-205-5p",]

# Control de calidad
total_lib <- apply(countdata_GIC, MARGIN = 2, sum)
epa <- data.frame(nombre = colnames(countdata_GIC, F), total = total_lib, grupo = calldata_GIC$Grupo)
ggplot(data = epa, aes(x = nombre, y = total, fill = grupo)) + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "") +
  xlab("Muestra") +
  ylab("Tamaño libreria") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Grafico de cajas
cpm_GIC <- cpm(countdata_GIC, log = T)
colnames(cpm_GIC) <- calldata_GIC$Grupo
az <- "dodgerblue2"
ro <- "firebrick3"
boxplot(cpm_GIC, las = 2, 
        col = ifelse(as.numeric(calldata_GIC$Grupo) == 1, ro, az), main = "", ylab = "Log2-CPM")
legend("topright", legend = c("Caso", "Normal"), fill = c(ro,az))
abline(h = mean(colMeans(cpm_GIC)), lwd = 2)

# Normalizacion
countdata_GIC_Norm <- DGEList(counts = countdata_GIC, group = calldata_GIC$Grupo)
countdata_GIC_Norm <- edgeR::calcNormFactors(countdata_GIC_Norm, method = "TMM")
max(countdata_GIC_Norm$samples$norm.factors)
min(countdata_GIC_Norm$samples$norm.factors)

# Grafico MDS
plotMDS(countdata_GIC_Norm, col = ifelse(as.numeric(countdata_GIC_Norm$samples$group) == 1, "firebrick2", "darkblue"), cex = 0.9,
        main = "")
grid()
abline(h = 0, v = 0, col = "gray60", lty = 1, lwd = 1)
legend("topright", legend = c("Normal", "Caso"), fill = c("darkblue","firebrick2"))

# Grafico MD
par(mfrow = c(1,2))
plotMD(DGEList(counts = countdata_GIC, group = calldata_GIC$Grupo), main = "Sin normalizar (A)")
abline(h = 0, col = "red")

plotMD(countdata_GIC_Norm, main = "Normalizado (B)")
abline(h = 0, col = "red")

#========================================
#     Analisis tendencia media-varianza
#========================================
design_GIC <- model.matrix(~ 0 + calldata_GIC$Grupo)
colnames(design_GIC) <- levels(countdata_GIC_Norm$samples$group)
rownames(design_GIC) <- rownames(calldata_GIC)
cont.matrix_GIC <- makeContrasts(CasovsNormal = Caso - Normal, levels = design_GIC)
v_GIC <- voom(counts = countdata_GIC_Norm, design = design_GIC, plot = T)

#=======================================================
#     Analisis coeficiente de variación biologica (BCV)
#=======================================================
x_GIC <- estimateDisp(countdata_GIC_Norm, design_GIC)
x_GIC$common.dispersion
BCV_GIC <- sqrt(x_GIC$common.dispersion)
plotBCV(x_GIC)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                       ANALISIS GRridge
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-------------------------------
# Preparacion final datos conteos
#-------------------------------
# Transformacion cuasi-gaussiana
funcion_1 <- function(x){
  sqrt(x + (3/8)) - sqrt(3/8)
}
conteoTransf_GIC <- apply(countdata_GIC, MARGIN = c(1,2), FUN = funcion_1)
conteoTransf_GIC <- as.data.frame(conteoTransf_GIC)

# Estandarizacion
conteoScale_GIC <- scale(t(conteoTransf_GIC))
conteoScale_GIC <- t(conteoScale_GIC)
conteoScale_GIC <- as.data.frame(conteoScale_GIC)

#======================================================================================
#                         INGRESO RESPUESTA Y CO-DATOS
#=======================================================================================

#------------------------------------------------
# Variable respuesta: NIC3: Caso, NORMAL: Control
#------------------------------------------------
respuesta <- calldata_GIC$Grupo

#------------
# Abundancia
#------------
abundancia <- rowSums(countdata_GIC)
parAbundancia <- CreatePartition(abundancia, mingr = 25, ngroup = 10, decreasing = TRUE, uniform = TRUE)

#--------------------
# Desviacion estandar
#--------------------
Destandar <- apply(conteoTransf_GIC, 1, sd)
parDestandar <- CreatePartition(Destandar, mingr = 25, ngroup = 10, decreasing = TRUE)

#--------------
# Conservacion
#---------------
miRNA_data <- read.delim("miR_Family_Info.txt", header = T)
miRNA_data <- subset(miRNA_data, Species.ID == 9606)
miRNA_data <- miRNA_data[,c("MiRBase.ID","Family.Conservation.")]
miRNA_data$MiRBase.ID <- gsub(pattern = "hsa-",replacement = "", x = miRNA_data$MiRBase.ID)

coincidir_con <- as.data.frame(rownames(countdata_GIC))
colnames(coincidir_con) <- "MiRBase.ID"
coincidir_con$MiRBase.ID <- gsub(pattern = "hsa-",replacement = "", x = coincidir_con$MiRBase.ID)


conservacion <- merge(x = coincidir_con, y = miRNA_data, by = "MiRBase.ID", all.x = T)
conservacion$Family.Conservation. <- ifelse(conservacion$Family.Conservation. == -1, 0,
                                            ifelse(conservacion$Family.Conservation. == 0, 1,
                                                   ifelse(conservacion$Family.Conservation. == 1,2,
                                                          ifelse(conservacion$Family.Conservation. == 2,3, "aja"))))

conservacion$Family.Conservation.[is.na(conservacion$Family.Conservation.)] <- 0
conservacion <- factor(conservacion$Family.Conservation.)
parConservacion <- CreatePartition(conservacion)

#-------------------------
# Combinacion de CO-DATOS
#--------------------------
ListaParticiones <- list(abundancia = parAbundancia, destandar = parDestandar,
                         conservacion = parConservacion)
ListaMonotona <- c(TRUE, FALSE, FALSE)

#=============================================
#             MODELACION GRride
#=============================================

#-------------------------
# Seleccion de particiones
#-------------------------
# Aplicar seleccion de particiones para optimizar el modelo GRridge
idTemp = 1:3
selecPar = PartitionsSelection(highdimdata = conteoScale_GIC, response = respuesta, 
                               partitions = ListaParticiones, monotoneFunctions = ListaMonotona,
                               optl = NULL, innfold = NULL)

# Lista de particiones que mejoran el rendimiento del modelo GRridge
particionActual <- ListaParticiones[2:3]

# Lista de funciones monotonas para la particion correspondiente
monotonaActual <- ListaMonotona[2:3]

#------------------------
#   Modelo GRridge
#------------------------
grM <- grridge(highdimdata = conteoScale_GIC, response = respuesta, partitions = particionActual,
               monotone = monotonaActual, optl = 326.9255, selectionEN = TRUE, maxsel = 15) 
grMCV <- grridgeCV(grr = grM, highdimdata = conteoScale_GIC, response = respuesta)
#-----------------------------------------------------------------------------------------------

#=================================
# RESULTADOS Y SELECCION POST-HOC
#=================================
# Ver resultados
grM$lambdamults # Parametros lambda

# Graficos curvas ROC y AUC
cutoffs <- rev(seq(0,1,by=0.01))
rocRidge <- roc(probs = grMCV[,2],
                true = grMCV[,1], cutoffs)
rocGRridge <- roc(probs = grMCV[,3],
                  true = grMCV[,1], cutoffs)
rocGRridgeEN <- roc(probs = grMCV[,4],
                    true = grMCV[,1], cutoffs)
plot(rocRidge[1,], rocRidge[2,], type = "l", lty = 2, 
     ann = T, col = "grey", ylab = "Sensibilidad", xlab = "Especificidad", lwd = 2)
points(rocGRridge[1,], rocGRridge[2,], type = "l",lty = 1,
       col = "black", lwd = 2)
points(rocGRridgeEN[1,], rocGRridgeEN[2,], type = "l",lty=1,
       col = "blue", lwd = 2)
#legend(0.0001,0.99, legend = c("ridge (AUC = 0.61)","GRridge (AUC = 0.61)","Panel 15 miRNA (AUC = 0.40)"),
#       col = c("grey","black","blue"), lty = c(1,1), horiz = F, box.lty = 1, cex = 0.72)
legend(0.43,0.19, legend = c("ridge (AUC = 0.61)","GRridge (AUC = 0.61)","Panel 15 miRNA (AUC = 0.40)"),
       col = c("grey","black","blue"), lty = c(1,1), horiz = F, box.lty = 1, cex = 1)

# Calculo Area bajo la curva
auc(rocRidge)
auc(rocGRridge)
auc(rocGRridgeEN)

# miRNAs seleccionados
betas <- as.data.frame(abs(grM$betas))
colnames(betas) <- "beta"
nombres <- rownames(countdata_GIC)
betas$miRNA <- nombres
betas1 <- betas[with(betas, order(betas$beta, decreasing = T)),][1:15,]

miRNAs_GRridge_GIC <- betas1$miRNA[1:15]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#             Graficos Boxplot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Grafico Boxplot - Seleccion GRridge
Filtro_boxplot <- countdata_GIC[miRNAs_GRridge_GIC, ]
Filtro_boxplot <- as.data.frame(t(Filtro_boxplot))
Filtro_boxplot$response <- respuesta
miRNAs_GRridge_GIC <- gsub(pattern = "hsa-",replacement = "", x = miRNAs_GRridge_GIC)


normales <- subset(Filtro_boxplot, response == "Normal")
normales <- normales[,-(dim(normales)[2])]

malos <- subset(Filtro_boxplot, response == "Caso")
malos <- malos[,-(dim(malos)[2])]

boxplot(sqrt(normales), boxwex = 0.9, col = "dodgerblue2",
        at = c(1,6,11,16,21,26,31,36,41,46,51,56,61,66,71), las = 2, ylim = c(0,9),
        xaxt = "n", main = "", ylab = "Raíz cuadrada conteos filtrados crudos")

grid()
boxplot(sqrt(normales), boxwex = 0.9, col = "dodgerblue2",
        at = c(1,6,11,16,21,26,31,36,41,46,51,56,61,66,71), las = 2, ylim = c(0,9),
        xaxt = "n", main = "Boxplot GIC", add = T)
boxplot(sqrt(malos),boxwex = 0.9, col = "firebrick3",
        at = c(2,7,12,17,22,27,32,37,42,47,52,57,62,67,72), las = 2, axes = F, add = TRUE)

axis(1, at = c(1.5,6.5,11.5,16.5,21.5,26.5,31.5,36.5,41.5,46.5,51.5,56.5,61.5,66.5,71.5), 
     labels = miRNAs_GRridge_GIC, cex.axis = 0.85, las = 2)        

legend("top", legend = c("Normal","Caso"), fill = c("dodgerblue2","firebrick3"),
       horiz = T, inset = 0.05, bty = "n")
#--------------------------------------------------------------------------------------------

##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################


### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
archivo <- file.choose()
datos <- read.delim(
  archivo,
  header = TRUE, 
  sep = "\t",
  check.names = FALSE
)
head(datos)
tail(datos)
dim(datos)

### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)

### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
control <- as.matrix(datos[, 2:4])
condicion <- as.matrix(datos[, 5:7])

datos$Media_Control <- rowMeans(control)
datos$DE_Control <- apply(control, 1, sd)

datos$Media_Condicion <- rowMeans(condicion)
datos$DE_Condicion <- apply(condicion, 1, sd)

head(datos)
dim(datos)

### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)

### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
color_control <- "#0072B2"    # Azul
color_condicion <- "#D55E00"  # Naranja

# Una fila con dos gráficos
par(mfrow = c(1, 2))

# Usar los mismos intervalos en ambos histogramas
cortes <- pretty(
  range(c(datos$Media_Control, datos$Media_Condicion)),
  n = 30
)
hist(
  datos$Media_Control,
  breaks = cortes,
  col = color_control,
  main = "Control",
  xlab = "Media por gen",
  ylab = "Número de genes"
)

hist(
  datos$Media_Condicion,
  breaks = cortes,
  col = color_condicion,
  main = "Condición",
  xlab = "Media por gen",
  ylab = "Número de genes"
)

# Volver a un único gráfico

par(mfrow = c(1, 1))

### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)

boxplot(
  list(
    Control = log10(datos$Media_Control + 1),
    Condición = log10(datos$Media_Condicion + 1)
  ),
  col = c(color_control, color_condicion),
  main = "Expresión media por gen",
  ylab = "log10(media + 1)"
)

legend(
  "topleft",
  legend = c("Control", "Condición"),
  fill = c(color_control, color_condicion),
  bty = "n"
)

### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.

shapiro.test(datos$Media_Control)
shapiro.test(datos$Media_Condicion)

### Se rechaza la normalidad de las distribuciones de las medias de Control y Condición

prueba_permutacion <- function(x, y) {
  
  valores <- c(x, y)
  
  # Diferencia observada entre los grupos
  diferencia_observada <- abs(mean(y) - mean(x))
  
  # Todos los repartos posibles de 6 valores en grupos de 3
  repartos <- combn(6, 3)
  
  # Diferencia de medias para cada reparto
  diferencias <- apply(repartos, 2, function(indices) {
    grupo1 <- valores[indices]
    grupo2 <- valores[-indices]
    
    abs(mean(grupo2) - mean(grupo1))
  })
  
  # Margen pequeño para evitar problemas de redondeo
  tolerancia <- 100 * .Machine$double.eps *
    max(1, abs(valores))
  
  # Proporción de diferencias iguales o mayores que la observada
  mean(diferencias >= diferencia_observada - tolerancia)
}

datos$p_value <- sapply(seq_len(nrow(datos)), function(i) {
  prueba_permutacion(
    control[i, ],
    condicion[i, ]
  )
})

head(datos)
range(datos$p_value)

# Interpretación del apartado 5:
# Los tests de Shapiro-Wilk de las medias de Control y Condición
# dan p < 0.05: se rechaza la normalidad de ambas distribuciones de medias entre genes.
# Esto no determina la normalidad de las réplicas de cada gen.
# Se utiliza una prueba exacta de permutación por gen,basada en la diferencia de medias.

### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

genes_DE <- datos[datos$p_value < 0.05, ]

# Nombres de los genes seleccionados
genes_DE[["Gene Name"]]

# Número total
nrow(genes_DE)

# Sobreexpresados: mayor media en la condición
sum(genes_DE$Media_Condicion > genes_DE$Media_Control)

# Reprimidos: menor media en la condición
sum(genes_DE$Media_Condicion < genes_DE$Media_Control)

# Respuesta de 6.¿Cuántos genes diferencialmente expresados hay?
# Con la prueba utilizada, ningún gen cumple p < 0.05.
# Por tanto, entre los genes seleccionados hay:
# 0 sobreexpresados y 0 reprimidos.
#
# Justificación: con tres valores por grupo existen 20 repartos.
# En esta prueba bilateral, el menor p-valor posible es 2/20 = 0.10.
# Por eso no se puede alcanzar p < 0.05.
# Esto no demuestra que no existan diferencias biológicas.

### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

### Con el umbral original p < 0.05 no se seleccionan genes.
### ponemos p a 0.10 para elaborar el heatmap.
### Usamos p <= 0.10 para incluir los valores iguales al corte.
### Se seleccionan 468 genes: 69 con mayor expresión media en la condición y 399 con menor expresión media.

genes_heatmap <- datos[datos$p_value <= 0.10, ]

# Número de genes seleccionados
nrow(genes_heatmap)  # 468

# Mayor y menor expresión media en la condición
sum(genes_heatmap$Media_Condicion >
      genes_heatmap$Media_Control)  # 69

sum(genes_heatmap$Media_Condicion <
      genes_heatmap$Media_Control)  # 399

# Matriz con las seis columnas de expresión
# Agrupación de genes y muestras, siguiendo el taller
matriz <- as.matrix(genes_heatmap[, 2:7])
rownames(matriz) <- genes_heatmap[["Gene Name"]]
rowv <- as.dendrogram(
  hclust(as.dist(1 - cor(t(matriz))))
)

colv <- as.dendrogram(
  hclust(as.dist(1 - cor(matriz)))
)

# Mostrar el gráfico en Plots
par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))

heatmap(
  matriz,
  Rowv = rowv,
  Colv = colv,
  scale = "row",
  labRow = NA,
  cexCol = 0.8,
  col = colorRampPalette(c("blue", "white", "red"))(100),
  ColSideColors = rep(
    c(color_control, color_condicion),
    each = 3
  ),
  margins = c(9, 4),
  main = "Genes seleccionados con p <= 0.10"
)

### El heatmap representa 468 genes seleccionados al ampliar
### el umbral. No son genes significativos con p < 0.05.
### Cada fila es un gen y cada columna una muestra.
### Los colores representan valores estandarizados por gen:
### rojo indica mayor expresión relativa y azul, menor.

### para guardar la tabla
write.table(
  datos,
  file = file.path(dirname(archivo), "resultados_trabajo.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

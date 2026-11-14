# Recuentos no ponderados de NHANES 2013-2018 (mujeres de 40 a 79 años).
library(foreign)
dir.create('nhanes_originales', showWarnings = FALSE)

leer <- function(year, suffix, module) {
  name <- paste0(module, '_', suffix, '.xpt')
  path <- file.path('nhanes_originales', name)
  url <- paste0('https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/', year, '/DataFiles/', name)
  if (!file.exists(path)) {
    tryCatch(download.file(url, path, mode = 'wb'),
             error = function(e) stop('No se pudo descargar ', name, ': ', conditionMessage(e)))
  }
  d <- read.xport(path)
  names(d) <- toupper(names(d))
  stopifnot(!anyDuplicated(d$SEQN))
  d
}

# Ciclos incluidos: 2013-14 (H), 2015-16 (I) y 2017-18 (J)
calendario <- data.frame(year = c(2013, 2015, 2017), sufijo = c('H', 'I', 'J'))
ciclos <- list()
for (i in seq_len(nrow(calendario))) {
  year <- calendario$year[i]
  suffix <- calendario$sufijo[i]
  demo <- leer(year, suffix, 'DEMO')
  rhq <- leer(year, suffix, 'RHQ')
  mcq <- leer(year, suffix, 'MCQ')
  d <- merge(demo[c('SEQN','RIAGENDR','RIDAGEYR','WTINT2YR','WTMEC2YR','SDMVSTRA','SDMVPSU','SDDSRVYR')],
             rhq[c('SEQN','RHQ420')], by = 'SEQN', all.x = TRUE)
  d <- merge(d, mcq[c('SEQN','MCQ220','MCQ230A','MCQ230B','MCQ230C','MCQ230D')],
             by = 'SEQN', all.x = TRUE)
  ciclos[[i]] <- d
}
d <- do.call(rbind, ciclos)
stopifnot(!anyDuplicated(d$SEQN))

# Tabla de flujo (registros que quedan en cada etapa, total y por ciclo)
flujo <- NULL
registrar <- function(etapa, d) {
  por_ciclo <- as.vector(table(factor(d$SDDSRVYR, levels = 8:10)))
  flujo <<- rbind(flujo, data.frame(etapa = etapa, n = nrow(d),
                                    H_2013 = por_ciclo[1],
                                    I_2015 = por_ciclo[2],
                                    J_2017 = por_ciclo[3]))
}
registrar('Registros totales', d)
d <- subset(d, RIAGENDR==2 & RIDAGEYR>=40 & RIDAGEYR<=79)
registrar('Mujeres de 40 a 79', d)
d <- subset(d, !is.na(WTMEC2YR) & WTMEC2YR>0)
registrar('Peso de examen positivo', d)
d <- subset(d, !is.na(SDMVSTRA) & !is.na(SDMVPSU))
registrar('Diseno muestral completo', d)
d <- subset(d, RHQ420 %in% c(1,2))
registrar('Exposicion valida', d)

# Diccionario de códigos MCQ230A-D (NHANES 2013-14): 10-39 = tipo de cáncer, 14 = mama,
# 77 = rechazo a responder, 99 = no sabe. 66 = "más de 3 tipos", registrado únicamente
# como cuarta respuesta, por lo que aparece solo en MCQ230D, nunca en A-C.
tipos <- as.matrix(d[c('MCQ230A','MCQ230B','MCQ230C')])
mama <- rowSums(tipos==14, na.rm = TRUE) > 0
tipo_valido <- rowSums(tipos>=10 & tipos<=39, na.rm = TRUE) > 0
mas_de_tres <- d$MCQ230D %in% 66

# Control de consistencia: quienes responden "no" a cáncer no deberían declarar ningún tipo.
cat('Mujeres sin cancer con un tipo declarado (deberia ser 0):',
    sum(d$MCQ220 %in% 2 & rowSums(!is.na(tipos)) > 0), '\n')

# Un diagnóstico de otro tumor no equivale a antecedente de mama.
# Los campos B/C vacíos pueden ser saltos previstos; no se excluyen por ello.
d$cancer_mama <- NA_integer_
d$cancer_mama[d$MCQ220 %in% 2] <- 0L
d$cancer_mama[d$MCQ220 %in% 1 & tipo_valido & !mama] <- 0L
d$cancer_mama[d$MCQ220 %in% 1 & mama] <- 1L

# Auditoría de respuestas ambiguas adicionales entre quienes se clasifican negativos.
ambiguos <- apply(tipos, 1, function(x) any(x %in% c(77,99)))
cat('Negativas con un tipo ambiguo o mas de tres tumores:',
    sum(d$cancer_mama==0 & (ambiguos | mas_de_tres), na.rm = TRUE), '\n')
d$respuesta_ambigua <- d$cancer_mama %in% 0 & (ambiguos | mas_de_tres)

# Resumen de las mujeres que se excluyen por no poder clasificar el antecedente
sin_clasificar <- is.na(d$cancer_mama)
cat('Mujeres excluidas por antecedente no clasificable:', sum(sin_clasificar), '\n')
print(table(MCQ220 = d$MCQ220[sin_clasificar],
            MCQ230A = d$MCQ230A[sin_clasificar], useNA = 'ifany'))

d <- subset(d, !is.na(cancer_mama))
registrar('Muestra final', d)

# Tabla principal: exposición (RHQ420) por antecedente de cáncer de mama
d$exposicion <- factor(d$RHQ420, levels = c(1,2), labels = c('Expuestas','No expuestas'))
tabla <- table(d$exposicion,
               factor(d$cancer_mama, levels = c(0,1),
                      labels = c('Sin antecedente','Con antecedente')))
print(flujo, row.names = FALSE)
print(addmargins(tabla))

# Análisis de sensibilidad: sin negativas con respuesta ambigua
d_sens <- subset(d, !respuesta_ambigua)
tabla_sens <- table(factor(d_sens$RHQ420, levels = c(1,2), labels = c('Expuestas','No expuestas')),
                    factor(d_sens$cancer_mama, levels = c(0,1),
                           labels = c('Sin antecedente','Con antecedente')))
cat('\nSensibilidad: sin negativas con respuesta ambigua\n')
print(addmargins(tabla_sens))

# Archivos de salida
write.csv(flujo, 'recuentos_verificados.csv', row.names = FALSE)
write.csv(as.data.frame.matrix(tabla), 'desglose_verificado.csv')
write.csv(as.data.frame.matrix(tabla_sens), 'desglose_sensibilidad.csv')
write.csv(d, 'muestra_analitica_verificada.csv', row.names = FALSE)
capture.output(sessionInfo(), file = 'sesion_R.txt')
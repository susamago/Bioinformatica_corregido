# Recuentos no ponderados de NHANES 2013-2018.
# Ejecutar desde la carpeta que contiene este script: Rscript comprobar_viabilidad.R
# Los originales se conservan en nhanes_originales; no se modifican.
# foreign forma parte de los paquetes recomendados de R.
library(foreign)
dir.create('nhanes_originales', showWarnings = FALSE)
leer <- function(year, suffix, module) {
  name <- paste0(module, '_', suffix, '.xpt')
  path <- file.path('nhanes_originales', name)
  url <- paste0('https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/', year, '/DataFiles/', name)
  if (!file.exists(path)) download.file(url, path, mode = 'wb')
  d <- read.xport(path)
  names(d) <- toupper(names(d))
  stopifnot(!anyDuplicated(d$SEQN))
  d
}
ciclos <- list()
for (i in seq_along(c(2013,2015,2017))) {
  year <- c(2013,2015,2017)[i]
  suffix <- c('H','I','J')[i]
  demo <- leer(year,suffix,'DEMO')
  rhq <- leer(year,suffix,'RHQ')
  mcq <- leer(year,suffix,'MCQ')
  d <- merge(demo[c('SEQN','RIAGENDR','RIDAGEYR','WTMEC2YR','SDMVSTRA','SDMVPSU','SDDSRVYR')],
             rhq[c('SEQN','RHQ420')], by='SEQN', all.x=TRUE)
  d <- merge(d, mcq[c('SEQN','MCQ220','MCQ230A','MCQ230B','MCQ230C','MCQ230D')], by='SEQN', all.x=TRUE)
  ciclos[[i]] <- d
}
d <- do.call(rbind,ciclos)
stopifnot(!anyDuplicated(d$SEQN))
flujo <- data.frame(etapa='Registros totales',n=nrow(d))
registrar <- function(etapa,d) {
  flujo <<- rbind(flujo,data.frame(etapa=etapa,n=nrow(d)))
}
d <- subset(d, RIAGENDR==2 & RIDAGEYR>=40 & RIDAGEYR<=79)
registrar('Mujeres de 40 a 79',d)
d <- subset(d,!is.na(WTMEC2YR) & WTMEC2YR>0)
registrar('Peso de examen positivo',d)
d <- subset(d,!is.na(SDMVSTRA) & !is.na(SDMVPSU))
registrar('Diseno muestral completo',d)
d <- subset(d,RHQ420 %in% c(1,2))
registrar('Exposicion valida',d)
tipos <- as.matrix(d[c('MCQ230A','MCQ230B','MCQ230C')])
mama <- rowSums(tipos==14,na.rm=TRUE)>0
tipo_valido <- rowSums(tipos>=10 & tipos<=39,na.rm=TRUE)>0
# Un diagnostico de otro tumor no equivale a antecedente de mama.
# Los campos B/C vacios pueden ser saltos previstos; no se excluyen por ello.
d$cancer_mama <- NA_integer_
d$cancer_mama[d$MCQ220 %in% 2] <- 0L
d$cancer_mama[d$MCQ220 %in% 1 & tipo_valido & !mama] <- 0L
d$cancer_mama[d$MCQ220 %in% 1 & mama] <- 1L
# Auditoria de respuestas ambiguas adicionales entre quienes se clasifican negativos.
ambiguos <- apply(tipos,1,function(x) any(x %in% c(66,77,99)))
cat('Negativas con un tipo ambiguo o mas de tres tumores:',
    sum(d$cancer_mama==0 & (ambiguos | d$MCQ230D %in% 66),na.rm=TRUE),'\n')
print(d[is.na(d$cancer_mama),c('MCQ220','MCQ230A','MCQ230B','MCQ230C')])
d <- subset(d,!is.na(cancer_mama))
registrar('Muestra final',d)
d$exposicion <- factor(d$RHQ420,levels=c(1,2),labels=c('Expuestas','No expuestas'))
tabla <- table(d$exposicion,factor(d$cancer_mama,levels=c(0,1),labels=c('Sin antecedente','Con antecedente')))
print(flujo,row.names=FALSE)
print(addmargins(tabla))
write.csv(flujo,'recuentos_verificados.csv',row.names=FALSE)
write.csv(as.data.frame.matrix(tabla),'desglose_verificado.csv')
write.csv(d,'muestra_analitica_verificada.csv',row.names=FALSE)
capture.output(sessionInfo(),file='sesion_R.txt')

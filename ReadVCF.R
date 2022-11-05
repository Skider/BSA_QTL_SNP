#-----------------------------------------------------------------------------
#Modulo 1: Analisis QTL-BSA. Programa "read_VCF_BSA"
#-----------------------------------------------------------------------------
#La funcion read_VCF_BSA permite leer archivos VCF creando un data frame
#con los datos necesario para el analisis BSA:
# - Region: Cromosoma [CHROM] y posicion [POS]
# - Lectura de referencia [REF] y alterada [ALT]
# - Valor de calidad de la mutacion [QUAL]
# - Tipo de Variante: SNP o INDEL [VARIANT]
# - Lecturas de referencia [AD_REF_PX] y alterada [AD_ALT_PX] para la muestra X
#-----------------------------------------------------------------------------
#La funcion devuelve un data frame con los datos necesarios para el analisis BSA
#-----------------------------------------------------------------------------
read_VCF_BSA <- function(file, f_QUAL = 0, f_DP = 0, output_file = ""){
  #Variables de la funcion:
  #file: Ruta y nombre del archivo VCF
  #f_QUAL: Valor de filtrado por QUAL (>=). Valor 0 por defecto
  #f_DP: Valor de filtrado por numero de lecturas (>=). Valor 0 por defecto
  #output_file: Ruta y nombre del archivo VCF filtrado. Por defecto, no se genera ningun archivo.
  #---------------------------------------------------------------------------
  #Librerias
  if (!require("dplyr", quietly = TRUE)){install.packages("dplyr")}
  if (!require("readr", quietly = TRUE)){install.packages("readr")}
  if (!require("tidyr", quietly = TRUE)){install.packages("tidyr")}
  if (!require("stringr", quietly = TRUE)){install.packages("stringr")}
  library(dplyr)
  library(readr, include.only = 'read_lines')
  library(tidyr, include.only = 'separate')
  library(stringr, include.only = c('str_split','str_detect'))
  defaultW <- getOption("warn")
  
  cat("Leyendo archivo VCF\n")
  
  #Saltar el encabezado del archivo VCF
  for (i in 0:10){
    vcf <- readr::read_lines(file, n_max = 5000, skip = i * 5000)
    for (n in 1:5000){
      if(substr(vcf[n],1,2) != "##"){
        break
      }
    }
    if (n != 5000){
      n <- n + (i * 5000)
      break
    }
  }
  rm(i, vcf)
  
  cat(sprintf("Lineas de metadatos: %i\n",n-1))
  #Nombre de las columnas del VCF
  namescol_df <- read.delim(file = file,
                            header = FALSE, sep = "\t",
                            skip = n-1,
                            nrows = 1)[1,] %>% unlist
  if(length(namescol_df) < 11){
    stop("No hay suficientes columnas en el archivo VCF.")
  }
  
  namescol_df[1] <- "CHROM"
  n_col_read <- length(namescol_df)
  n_muestras <- n_col_read - 9
  cat(sprintf("Columna de datos: %i\n",n_col_read))
  cat(sprintf("Muestras leidas: %i\n",n_muestras))
  for(i in 10:n_col_read){
    cat(sprintf(" - %s\n",namescol_df[i]))
  }


  cat("Archivo VCF leido. Generando tabla de resultados VCF...\n")
  pasos <- n_muestras + 3
  cat(sprintf("[1/%i pasos] Lectura de los datos...\n",pasos))
  
  #Lectura de los datos del archivo VCF
  vcf_df <- read.delim(file = file,
                       header = FALSE, sep = "\t",
                       col.names = namescol_df,
                       skip = n)
  rm(i, n, namescol_df, n_col_read)
  cat(sprintf("Variantes leidas: %i\n",nrow(vcf_df)))
  
  cat(sprintf("[2/%i pasos] Filtrado de datos por QUAL\n",pasos))
  #Creamos el data frame con los datos filtrados
  vcf_df <- vcf_df %>% filter(QUAL >= f_QUAL)
  cat(sprintf("Variantes leidas tras el filtrado: %i\n",nrow(vcf_df)))
  
  #Comprobacion de la profundidad alelica
  if(sum(stringr::str_detect(vcf_df$FORMAT,"AD")) == 0){
    stop("El archivo vcf no contiene datos de profundidad alelica para REF-ALT.")
  }
  
  #Obtencion del contenido de la columna formato
  formato <- vcf_df$FORMAT[1] %>% stringr::str_split(pattern = ":") %>% unlist
  
  vcf_df_output <- vcf_df %>% select(CHROM,POS,REF,ALT,QUAL) %>%
    mutate(VARIANT=ifelse(stringr::str_detect(vcf_df$INFO,"INDEL") == TRUE, "INDEL","SNP"))

  options(warn = -1)  
  for(i in 1:n_muestras){
    cat(sprintf("[%i/%i pasos] Lectura de la profundidad alelica de la Pool %i...\n",i+2,pasos,i))
    #Obtencion de los datos de AD (REF y ALT).
    #NOTA: Las profundidades de diferentes variantes para una misma posicion seran sumadas

    df_test <- data.frame(INFO = vcf_df[,(i+9)]) %>%
      tidyr::separate(col=INFO, into=formato, sep=":") %>%
      select(AD) %>%
      tidyr::separate(col = AD, into = c("REF","ALT1","ALT2"), sep = ",")
    df_test[] <- lapply(df_test, function(x) as.numeric(as.character(x)))
    df_test$ALT <- rowSums(df_test[,c(2:3)],na.rm=TRUE)
    
    vcf_df_output[,(7+(i-1)*2)] <- df_test$REF
    vcf_df_output[,(8+(i-1)*2)] <- df_test$ALT
    names(vcf_df_output)[(7+(i-1)*2)] <- sprintf("AD_REF_P%i",i)
    names(vcf_df_output)[(8+(i-1)*2)] <- sprintf("AD_ALT_P%i",i)
    if(sum(is.na(vcf_df_output[,(7+(i-1)*2)])) != 0 | sum(is.na(vcf_df_output[,(8+(i-1)*2)])) != 0){
      stop(sprintf("Alguna variante del Pool %i no tiene un valor de profundidad valido.",i-9))
    }
  }
  options(warn = defaultW)
  
  cat(sprintf("[%i/%i pasos] Filtrado de datos por numero de lecturas por posicion\n",pasos,pasos))
  if(f_DP != 0){
    vcf_df_output <- subset(vcf_df_output, rowSums(vcf_df_output[,c(7:ncol(vcf_df_output))]) >= f_DP)
  }
  cat(sprintf("Variantes leidas tras el filtrado por numero de lecturas: %i\n",nrow(vcf_df_output)))
  
  if(output_file != ""){
    cat(sprintf("Guardando los datos en el fichero %s\n",output_file))
    write.table(vcf_df_output, file = output_file, sep="\t", row.names = FALSE, quote = FALSE)
  }
  cat("Terminado!\n")
  
  return(vcf_df_output)
}

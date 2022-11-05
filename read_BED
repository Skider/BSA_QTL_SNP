read_BED <- function(file){
  #Librerias
  if (!require("dplyr", quietly = TRUE)){install.packages("dplyr")}
  if (!require("readr", quietly = TRUE)){install.packages("readr")}
  library(dplyr)
  library(readr, include.only = 'read_lines')
  
  cat(sprintf("Leyendo el archivo BED...\n"))
  
  #Saltar el encabezado del archivo BED
  for (i in 0:10){
    df_tmp <- readr::read_lines(file, n_max = 1000, skip = i * 1000)
    for (n in 1:1000){
      if(substr(df_tmp[n],1,7) != "browser" & substr(df_tmp[n],1,5) != "track"){
        break
      }
    }
    if (n != 1000){
      n <- n + (i * 1000)
      break
    }
  }
  
  #Comprobacion del numero de columnas del archivo BED
  df_tmp <- read.delim(file, nrows = 1, skip = n)
  if(ncol(df_tmp) < 6){
    stop("No hay suficientes columnas en el archivo BED.")
  }
  rm(df_tmp)
  
  #Lectura del archivo BED
  bed_df <- read.delim(file = file, sep="\t", header=FALSE, skip = n) %>%
    select(V1,V2,V3,V4,V6) %>%
    rename(Region = V1, Gene = V4, Dir = V6, Initial_POS = V2, Final_POS = V3) %>%
    #Algunos programas generan el archivo con "ID=Gene", para evitar problemas de lectura
    #de la secuencia del archivo FASTA de ARNm, se elimina "ID="
    mutate(Gene = gsub(".*=","",Gene))
  
  #Ordenado del data-frame por region y posicion
  bed_df <- bed_df[order(bed_df$Region, bed_df$Initial_POS),]
  
  return(bed_df)
}

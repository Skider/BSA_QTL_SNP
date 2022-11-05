#-----------------------------------------------------------------------------
#Modulo 2: Identificacion de SNP. Programa "BSA_SNP_Gene_2"
#-----------------------------------------------------------------------------
#La funcion BSA_SNP_Gene_2 realiza el proceso final de la funcion BSA_SNP_Gen
#comprobando que las mutaciones se encuentran dentro del ARNm del gen.
#-----------------------------------------------------------------------------
#La funcion devuelve un data-frame con 3 columnas adicionales para cada mutacion:
# - Gene: Gen asociado a la mutacion
# - Dir: Direccion del gen asociado con respecto al genoma de referencia
# - Check_Gene: Comprobacion de si la mutacion esta dentro del Gen (Valor 1), fuera
#del gen (Valor 0) o se encuentra en uno de los extremos o en direccion opuesta al gen (Valor ?)
#-----------------------------------------------------------------------------
BSA_SNP_Gene_2 <- function(vcf_df,
                         bed_df,
                         rna_fasta,
                         genome_fasta,
                         chr_list = "CHROM",
                         type = "Completo",
                         POS_Start = 1, POS_End = 1000000,
                         position = 1000000, range = 500000,
                         check_seq_pb = 20
){
  #Variables de entrada
  #vcf_df: Data Frame con los datos
  #bed_df: Data Frame o ruta del archivo BED de los genes del genoma de referencia
  #rna_fasta: Formato FASTA de seqinr o ruta del archivo FASTA del ARNm de los genes
  #genome_fasta: Formato FASTA de seqinr o ruta del archivo FASTA del genoma
  #chr_list: Lista de cromosomas a analizar:
  #   La opcion por defecto: "CHROM", utiliza solo los cromosomas
  #   "Todos", utiliza todos los valores en la columna del data frame de entrada
  #Eleccion de la region a analizar:
  #   type: Completo: Todo el cromosoma. Eleccion por defecto.
  #         Rango: Rango fijado entre POS_Start y POS_End. Por defecto entre 1-1.000.000
  #         Punto: Rango en un Intervalo a la derecha y a la izquierda de Posicion.
  #                 Por defecto en la position 1.000.000 con un range de mas-menos 500.000
  #---------------------------------------------------------------------------
  #Librerias
  if (!require("seqinr", quietly = TRUE)){install.packages("seqinr")}
  if (!require("dplyr", quietly = TRUE)){install.packages("dplyr")}
  library(seqinr, include.only = 'read.fasta')
  library(dplyr)
  
  #Comprobacion de si la mutacion esta dentro del ARNm del gen
  check_pc <- 1
  for (i in 1:nrow(vcf_df)){
    #Saltarse las mutaciones sin gen asignado
    if(vcf_df$Gene[i] == ""){
      next
    }
    if(i == as.integer(check_pc * nrow(vcf_df) * 0.1)){
      cat(sprintf("Analizado el %d%% de los SNP\n",check_pc*10))
      check_pc <- check_pc + 1
    }
    #Obtencion de la secuencia previo y posterior a la mutacion (20 pares de base)
    POS_SNP_tmp <- as.integer(vcf_df$POS[i])
    
    seq_Previa_SNP_tmp <- tolower(paste(genome_fasta["name" = vcf_df$CHROM[i]][[1]][(POS_SNP_tmp-check_seq_pb):(POS_SNP_tmp)],collapse=""))
    seq_Post_SNP_tmp <- tolower(paste(genome_fasta["name" = vcf_df$CHROM[i]][[1]][(POS_SNP_tmp):(POS_SNP_tmp+check_seq_pb)],collapse=""))
    ARNm_tmp <- tolower(paste(rna_fasta["name" = vcf_df$Gene[i]][[1]],collapse=""))
    
    #Busca la secuencia anterior y posterior a la mutacion en ambos sentidos del ARNm
    check_Gene <- c(0,0,0,0) #Previa_SNP_+;Post_SNP_+;Previa_SNP_-;Post_SNP_-
    if(unlist(gregexpr(seq_Previa_SNP_tmp,ARNm_tmp)) != -1){
      check_Gene[1] <- 1
    }
    if(unlist(gregexpr(seq_Post_SNP_tmp,ARNm_tmp)) != -1){
      check_Gene[2] <- 1
    }
    if(unlist(gregexpr(reverse_dna(seq_Post_SNP_tmp),ARNm_tmp)) != -1){
      check_Gene[3] <- 1
    }
    if(unlist(gregexpr(reverse_dna(seq_Previa_SNP_tmp),ARNm_tmp)) != -1){
      check_Gene[4] <- 1
    }
    
    #Si se ha obtenido coincidencia en la direccion positiva
    if(paste(check_Gene,collapse="") == "1100"){
      #y coincide con la direccion del gen
      if(vcf_df$Dir[i] == "+"){
        vcf_df$Check_Gene[i] <- 1
      } else {
        vcf_df$Check_Gene[i] <- "?"
      }
    #Si se ha obtenido coincidencia en la direccion negativa
    } else if(paste(check_Gene,collapse="") == "0011"){
      #y coincide con la direccion del gen
      if(vcf_df$Dir[i] == "-"){
        vcf_df$Check_Gene[i] <- 1
      } else {
        vcf_df$Check_Gene[i] <- "?"
      }
    #Si no se ha obtenido ninguna coincidencia
    } else if(paste(check_Gene,collapse="") == "0000"){
      vcf_df$Check_Gene[i] <- 0
    #Casos no completados anteriormente
    } else {
      vcf_df$Check_Gene[i] <- "?"
    }
  }
  
  return(vcf_df)
}

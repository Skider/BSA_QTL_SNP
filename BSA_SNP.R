#-----------------------------------------------------------------------------
#Modulo 2: Identificacion de SNP. Programa "BSA_SNP"
#-----------------------------------------------------------------------------
#La funcion BSA_SNP analiza el data-frame vcf_df buscando mutaciones que creen
#o rompan un codón de parada o inicio en uno de los 3 posibles marcos de lecturas
#a partir de la secuencia extraída del genome_fasta
#El data frame de entrada debe tener, como minimo, los siguientes datos:
# - Region: Cromosoma [CHROM] y posicion [POS]
#-----------------------------------------------------------------------------
#La funcion devuelve un data-frame con 8 columnas que indican si se ha producido
#alguno de los 4 posibles casos de relevancia en uno u otra direccion.
#-----------------------------------------------------------------------------
BSA_SNP <- function(vcf_df,
                     genome_fasta,
                     chr_list = "CHROM",
                     type = "Rango",
                     POS_Start = 1, POS_End = 1000000,
                     position = 1000000, range = 500000,
                     filter_type = "ge",
                     f_DeltaSNP_1 = 0.75, f_DeltaSNP_2 = -0.75,
                     filter = TRUE
                     ){
  #Variables de entrada
  #vcf_df: Data Frame con los datos
  #genome_fasta: Formato FASTA de seqinr o ruta del archivo FASTA del genoma
  #chr_list: Lista de cromosomas a analizar:
  #   La opcion por defecto: "CHROM", utiliza solo los cromosomas
  #   "Todos", utiliza todos los valores en la columna del data frame de entrada
  #Eleccion de la region a analizar:
  #   type: Completo: Todo el cromosoma. Eleccion por defecto.
  #         Rango: Rango fijado entre POS_Start y POS_End. Por defecto entre 1-1.000.000
  #         Punto: Rango en un Intervalo a la derecha y a la izquierda de Posicion.
  #                 Por defecto en la position 1.000.000 con un range de mas-menos 500.000
  #filter_type: Eleccion del filtro de Delta_SNPIndex a realizar.
  #   Opciones: gt (Mayor que), ge (Mayor o igual)
  #             lt (Menor que), le (Menor o igual)
  #             between (Entre f_DeltaSNP_1 y f_DeltaSNP_2, ambos inclusive),
  #             out (Fuera del intervalor entre f_DeltaSNP_1 y f_DeltaSNP_2)
  #   f_DeltaSNP_1 es el valor del filtro a realizar
  #   f_DeltaSNP_2 solo se utilizara si se escoge "between" o "out"
  #filter: En caso de ser TRUE, devuelve el data-frame solo con los SNP relevantes
  #---------------------------------------------------------------------------
  #Librerias
  if (!require("seqinr", quietly = TRUE)){install.packages("seqinr")}
  if (!require("dplyr", quietly = TRUE)){install.packages("dplyr")}
  library(seqinr, include.only = 'read.fasta')
  library(dplyr)
  
  #Comprobacion si estan los datos necesarios para el analisis.
  #Si no estan los valores SNPIndex calculado, se comprueba si estan los valores de AD y se
  #calculan los valores de SNPIndex.
  if(sum(colnames(vcf_df) == "SNPIndex_P1") != 1 | sum(colnames(vcf_df) == "SNPIndex_P2") != 1 |
     sum(colnames(vcf_df) == "Delta_SNPIndex") != 1){
    if(sum(colnames(vcf_df) == "AD_REF_P1") != 1 | sum(colnames(vcf_df) == "AD_ALT_P1") != 1 |
       sum(colnames(vcf_df) == "AD_REF_P2") != 1 | sum(colnames(vcf_df) == "AD_ALT_P2") != 1){
      stop("El data frame no contiene los datos de frecuencia alelica ni los SNP Index necesarios.")
    } else {
      cat(sprintf("Generando los SNP Index...\n"))
      vcf_df <- getSNPIndex(vcf_df)
    }
  }

  #Se comprueba si se ha dado el genoma en formata FASTA de seqinr, en caso contrario,
  #se lee un archivo fasta con el nombre y la ruta indicada
  if(class(genome_fasta) == "character"){
    cat(sprintf("Leyendo el genoma de referencia...\n"))
    genome_fasta <- seqinr::read.fasta(genome_fasta,seqtype = "DNA")
  } else if(class(genome_fasta) != "list") {
    stop("Introduzca un genoma valido o indica la ruta del archivo FASTA.")
  }
  
  
  #Obtencion del vector de los cromosomas
  if(chr_list[1] == "Todos"){
    chr_list <- unique(vcf_df$CHROM)
    cat(sprintf("Numero de entradas obtenidas: %i entradas\n",length(chr_list)))
  } else if(chr_list[1] == "CHROM"){
    chr_list <- unique(vcf_df$CHROM)
    chr_list <- chr_list[grepl("chr*",chr_list)]
    cat(sprintf("Numero de cromosomas obtenidos: %i cromosomas\n",length(chr_list)))
  }
  
  #Creacion del data frame solo con los cromosomas escogidos: vcf_df_chr y solo las mutaciones SNP
  cat(sprintf("Filtrando los SNP de las regiones de interes...\n"))
  vcf_df_chr <- subset(vcf_df, CHROM %in% chr_list & VARIANT == "SNP")
  #Creacion del data frame para regiones especificas
  if(type == "Punto"){
    POS_Start <- Posicion - Intervalo
    POS_End <- Posicion + Intervalo
    vcf_df_chr <- subset(vcf_df_chr, POS >= POS_Start & POS <= POS_End)
  } else if(type == "Rango"){
    vcf_df_chr <- subset(vcf_df_chr, POS >= POS_Start & POS <= POS_End)
  }
  
  #Comprobacion de errores
  if(typeof(vcf_df_chr) != "list" | class(vcf_df_chr) != "data.frame" | nrow(vcf_df_chr) == 0){
    stop("Error en la adquisicion de los datos para la grafica.")
  }
  
  #Filtrados de los SNP en relacion al valor de Delta_SNPIndex
  #Comprobacion de errores en los parametros de entrada
  if(class(f_DeltaSNP_1) != "numeric"){
    stop("Indique un valor numero para f_DeltaSNP_1.")
  }
  if(class(f_DeltaSNP_2) != "numeric" & (filter_type == "between" | filter_type == "out")){
    stop("Indique un valor numero para f_DeltaSNP_2.")
  }
  if(filter_type == "gt"){
    vcf_df_chr <- subset(vcf_df_chr, Delta_SNPIndex > f_DeltaSNP_1)
  } else if(filter_type == "ge"){
    vcf_df_chr <- subset(vcf_df_chr, Delta_SNPIndex >= f_DeltaSNP_1)
  } else if(filter_type == "lt"){
    vcf_df_chr <- subset(vcf_df_chr, Delta_SNPIndex < f_DeltaSNP_1)
  } else if(filter_type == "le"){
    vcf_df_chr <- subset(vcf_df_chr, Delta_SNPIndex <= f_DeltaSNP_1)
  } else if(filter_type == "between"){
    vcf_df_chr <- subset(vcf_df_chr, Delta_SNPIndex <= f_DeltaSNP_1 & Delta_SNPIndex >= f_DeltaSNP_2)
  } else if(filter_type == "out"){
    vcf_df_chr <- subset(vcf_df_chr, Delta_SNPIndex > f_DeltaSNP_1 | Delta_SNPIndex < f_DeltaSNP_2)
  } else{
    stop("Escoja un filtro de DeltaSNP valido.")
  }
  
  #Analisis de los SNP que generan o destruyen codones de parada o de inicio
  cat(sprintf("Analizando %d SNP en las regiones de interes...\n",nrow(vcf_df_chr)))
  new_seq_SNP_tmp <- ""
  seq_SNP_tmp <- ""
  check_pc <- 1
  for (i in 1:nrow(vcf_df_chr)){
    if(i == as.integer(check_pc * nrow(vcf_df_chr) * 0.1)){
      cat(sprintf("Analizado el %d%% de los SNP\n",check_pc*10))
      check_pc <- check_pc + 1
    }
    #Obtencion de los tres posibles marco de lecturas para el SNP
    POS_SNP_tmp <- as.integer(vcf_df_chr$POS[i])

    #Tres marcos de lecturas de la referencia
    seq_SNP_tmp[1] <- paste(genome_fasta["name" = vcf_df_chr$CHROM[i]][[1]][(POS_SNP_tmp):(POS_SNP_tmp+2)],collapse="")
    seq_SNP_tmp[2] <- paste(genome_fasta["name" = vcf_df_chr$CHROM[i]][[1]][(POS_SNP_tmp-1):(POS_SNP_tmp+1)],collapse="")
    seq_SNP_tmp[3] <- paste(genome_fasta["name" = vcf_df_chr$CHROM[i]][[1]][(POS_SNP_tmp-2):(POS_SNP_tmp)],collapse="")
    #Tres marcos de lecturas de la alternativa
    new_seq_SNP_tmp[1] <- paste(tolower(vcf_df_chr$ALT[i]),paste(genome_fasta["name" = vcf_df_chr$CHROM[i]][[1]][(POS_SNP_tmp+1):(POS_SNP_tmp+2)],collapse=""),sep="")
    new_seq_SNP_tmp[2] <- paste(paste(genome_fasta["name" = vcf_df_chr$CHROM[i]][[1]][(POS_SNP_tmp-1)],tolower(vcf_df_chr$ALT[i]),sep=""),genome_fasta["name" = vcf_df_chr$CHROM[i]][[1]][(POS_SNP_tmp+1)],sep="")
    new_seq_SNP_tmp[3] <- paste(paste(genome_fasta["name" = vcf_df_chr$CHROM[i]][[1]][(POS_SNP_tmp-2):(POS_SNP_tmp-1)],collapse=""),tolower(vcf_df_chr$ALT[i]),sep="")
    
    #Generacion de un codon de inicio ATG/CAT
    vcf_df_chr$Inicio_POS[i] <- 0
    vcf_df_chr$Inicio_NEG[i] <- 0
    for (j in 1:3){
      if(new_seq_SNP_tmp[j] == "atg"){
        vcf_df_chr$Inicio_POS[i] <- j
      }
      if(new_seq_SNP_tmp[j] == "cat"){
        vcf_df_chr$Inicio_NEG[i] <- j
      }
    }
    
    #Destruccion de un codon de inicio ATG/CAT
    vcf_df_chr$Ruptura_Inicio_POS[i] <- 0
    vcf_df_chr$Ruptura_Inicio_NEG[i] <- 0
    for (j in 1:3){
      if(seq_SNP_tmp[j] == "atg"){
        vcf_df_chr$Ruptura_Inicio_POS[i] <- j
      }
      if(seq_SNP_tmp[j] == "cat"){
        vcf_df_chr$Ruptura_Inicio_NEG[i] <- j
      }
    }
    
    #Generacion de un codon de parada TGA/TCA o TAG/CTA o TAA/TTA que anteriormente no lo fuese
    vcf_df_chr$Stop_POS[i] <- 0
    vcf_df_chr$Stop_NEG[i] <- 0
    for (j in 1:3){
      if(new_seq_SNP_tmp[j] == "tga" | new_seq_SNP_tmp[j] == "tag" | new_seq_SNP_tmp[j] == "taa"){
        if(seq_SNP_tmp[j] != "tga" & seq_SNP_tmp[j] != "tag" & seq_SNP_tmp[j] != "taa"){
          vcf_df_chr$Stop_POS[i] <- j
        }
      }
      if(new_seq_SNP_tmp[j] == "tca" | new_seq_SNP_tmp[j] == "cta" | new_seq_SNP_tmp[j] == "tta"){
        if(seq_SNP_tmp[j] != "tca" & seq_SNP_tmp[j] != "cta" & seq_SNP_tmp[j] != "tta"){
          vcf_df_chr$Stop_NEG[i] <- j
        }
      }
    }
    
    #Destruccion de un codon de parada TGA/TCA o TAG/CTA o TAA/TTA sin generar un nuevo codon de parada
    vcf_df_chr$Ruptura_Stop_POS[i] <- 0
    vcf_df_chr$Ruptura_Stop_NEG[i] <- 0
    for (j in 1:3){
      if(seq_SNP_tmp[j] == "tga" | seq_SNP_tmp[j] == "tag" | seq_SNP_tmp[j] == "taa"){
        if(new_seq_SNP_tmp[j] != "tga" & new_seq_SNP_tmp[j] != "tag" & new_seq_SNP_tmp[j] != "taa"){
          vcf_df_chr$Ruptura_Stop_POS[i] <- j
        }
      }
      if(seq_SNP_tmp[j] == "tca" | seq_SNP_tmp[j] == "cta" | seq_SNP_tmp[j] == "tta"){
        if(new_seq_SNP_tmp[j] != "tca" & new_seq_SNP_tmp[j] != "cta" & new_seq_SNP_tmp[j] != "tta"){
          vcf_df_chr$Ruptura_Stop_NEG[i] <- j
        }
      }
    }
  }
  if (filter == TRUE){
    vcf_df_chr <- subset(vcf_df_chr, (Inicio_POS + Inicio_NEG +
                           Ruptura_Inicio_POS + Ruptura_Inicio_NEG +
                           Stop_POS + Stop_NEG +
                           Ruptura_Stop_POS + Ruptura_Stop_NEG) != 0)
  }
  return(vcf_df_chr)
}

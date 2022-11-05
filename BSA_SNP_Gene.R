#-----------------------------------------------------------------------------
#Modulo 2: Identificacion de SNP. Programa "BSA_SNP_Gene"
#-----------------------------------------------------------------------------
#La funcion BSA_SNP_ORF analiza el data-frame vcf_df asignado a cada mutacion a un
#gen de los listados en un archivo BED y comprobando si se encuentra dentro de la
#secuencia del mismo utilizando un ficho FASTA con el ARNm de los genes.
#El data frame vcf_df de entrada debe tener, como minimo, los siguientes datos:
# - Region: Cromosoma [CHROM] y posicion [POS]
#El data frame bed de entrada debe tener, como minimo, los siguientes datos:
# - Region: Region del genoma [Columna 1]
# - Initial_POS: Posicion inicial del gen [Columna 2]
# - Final_POS: Posicion final del gen [Columna 3]
# - Gene: Nombre o identificador del gen [Columna 4]
# - Dir: Direccion del gen con respecto al genoma de referencia [Columna 6]
#-----------------------------------------------------------------------------
#La funcion devuelve un data-frame con 3 columnas adicionales para cada mutacion:
# - Gene: Gen asociado a la mutacion
# - Dir: Direccion del gen asociado con respecto al genoma de referencia
# - Check_Gene: Comprobacion de si la mutacion esta dentro del Gen (Valor 1), fuera
#del gen (Valor 0) o se encuentra en uno de los extremos o en direccion opuesta al gen (Valor ?)
#-----------------------------------------------------------------------------
BSA_SNP_Gene <- function(vcf_df,
                    bed_df,
                    rna_fasta,
                    genome_fasta,
                    chr_list = "CHROM",
                    type = "Completo",
                    POS_Start = 1, POS_End = 1000000,
                    position = 1000000, range = 500000,
                    check_seq_pb = 20,
                    filter = FALSE
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
  #check_seq_pb: Numero de pares de base de comprobacion de la secuencia antes y despues de la mutacion.
  #filter: En caso de ser TRUE, devuelve el data-frame solo con las mutaciones
  #que se encuentran dentro de un gen (Check_Gene = "1")
  #---------------------------------------------------------------------------
  #Librerias
  if (!require("seqinr", quietly = TRUE)){install.packages("seqinr")}
  if (!require("dplyr", quietly = TRUE)){install.packages("dplyr")}
  library(seqinr, include.only = 'read.fasta')
  library(dplyr)
  
  #Se comprueba si se ha dado el genoma en formato FASTA de seqinr, en caso contrario,
  #se lee un archivo fasta con el nombre y la ruta indicada
  if(class(genome_fasta) == "character"){
    cat(sprintf("Leyendo el archivo FASTA del genoma...\n"))
    genome_fasta <- seqinr::read.fasta(genome_fasta,seqtype = "DNA")
  } else if(class(genome_fasta) != "list") {
    stop("Introduzca un vector FASTA de seqinr del genoma valido o indica la ruta del mismo.")
  }
  
  #Se comprueba si se ha dado el ARNm en formato FASTA de seqinr, en caso contrario,
  #se lee un archivo fasta con el nombre y la ruta indicada
  if(class(rna_fasta) == "character"){
    cat(sprintf("Leyendo el archivo FASTA de ARNm...\n"))
    rna_fasta <- seqinr::read.fasta(rna_fasta,seqtype = "DNA")
  } else if(class(rna_fasta) != "list") {
    stop("Introduzca un vector FASTA de seqinr del ARNm valido o indica la ruta del mismo.")
  }
  
  #Comprobacion de los datos del data-frame BED o lectura del mismo
  if(class(bed_df) == "data.frame"){
    if(sum(colnames(bed_df) == "Region") != 1 | sum(colnames(bed_df) == "Initial_POS") != 1 |
       sum(colnames(bed_df) == "Final_POS") != 1 | sum(colnames(bed_df) == "Gene") != 1 |
       sum(colnames(bed_df) == "Dir") != 1){
      stop("El data frame no contiene los datos necesarios. Revise el nombre de las columnas.")
    }
  } else if (class(bed_df) == "character"){
    bed_df <- read_BED(bed_df)
  } else {
    stop("Por favor, introduzca una ruta valida para el archivo BED.")
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
  
  #Creacion del data frame vcf_df_chr y filtrado del data-frame BED para los
  #cromosomas escogidos
  cat(sprintf("Filtrando las mutaciones de las regiones de interes...\n"))
  vcf_df_chr <- subset(vcf_df, CHROM %in% chr_list)
  bed_df_chr <- subset(bed_df, Region %in% chr_list)
  #Creacion del data frame para regiones especificas
  if(type == "Punto"){
    POS_Start <- Posicion - Intervalo
    POS_End <- Posicion + Intervalo
    vcf_df_chr <- subset(vcf_df_chr, POS >= POS_Start & POS <= POS_End)
    bed_df_chr <- subset(bed_df_chr, Initial_POS >= POS_Start & Final_POS <= POS_End)
  }
  if(type == "Rango"){
    vcf_df_chr <- subset(vcf_df_chr, POS >= POS_Start & POS <= POS_End)
    bed_df_chr <- subset(bed_df_chr, Initial_POS >= POS_Start & Final_POS <= POS_End)
  }
  
  
  #Comprobacion de errores
  if(typeof(vcf_df_chr) != "list" | class(vcf_df_chr) != "data.frame" | nrow(vcf_df_chr) == 0){
    stop("Error en la adquisicion de los datos para la grafica.")
  }
  
  #Asignacion de la mutacion a un gen
  cat(sprintf("Asignando el gen a %d mutaciones...\n",nrow(vcf_df_chr)))
  vcf_df_chr$Gene <- ""
  vcf_df_chr$Dir <- ""
  check_pc <- 1
  start_count <- 1
  empty_region <- 0
  bed_df_tmp <- subset(bed_df_chr, Region == vcf_df_chr$CHROM[1])
  for (j in 1:nrow(bed_df_tmp)){
    if(vcf_df_chr$CHROM[1] == bed_df_tmp$Region[j] &
       as.integer(vcf_df_chr$POS[1]) >= as.integer(bed_df_tmp$Initial_POS[j]) &
       as.integer(vcf_df_chr$POS[1]) <= as.integer(bed_df_tmp$Final_POS[j])){
      vcf_df_chr$Gene[1] <- bed_df_tmp$Gene[j]
      vcf_df_chr$Dir[1] <- bed_df_tmp$Dir[j]
      start_count <- j
      break
    }
  }
  for (i in 2:nrow(vcf_df_chr)){
    if(i == as.integer(check_pc * nrow(vcf_df_chr) * 0.1)){
      cat(sprintf("Analizado el %d%% de los SNP\n",check_pc*10))
      check_pc <- check_pc + 1
    }
    #Si se cambia de region, se reinicia el contador "start_count" y se genera
    #un nuevo bed_df con la region
    if(vcf_df_chr$CHROM[i] != vcf_df_chr$CHROM[i-1] | empty_region == 1){
      bed_df_tmp <- subset(bed_df_chr, Region == vcf_df_chr$CHROM[1])
      #Comprobacion de que hay genes en la region analizada
      if(nrow(bed_df_tmp) == 0){
        empty_region <- 1
        next
      } else {
        empty_region <- 0
      }
      start_count <- 1
    }
    for (j in start_count:nrow(bed_df_tmp)){
      if(vcf_df_chr$CHROM[i] == bed_df_tmp$Region[j] & 
         as.integer(vcf_df_chr$POS[i]) >= as.integer(bed_df_tmp$Initial_POS[j]) &
         as.integer(vcf_df_chr$POS[i]) <= as.integer(bed_df_tmp$Final_POS[j])){
        vcf_df_chr$Gene[i] <- bed_df_tmp$Gene[j]
        vcf_df_chr$Dir[i] <- bed_df_tmp$Dir[j]
        start_count <- j
        break
      } else if (vcf_df_chr$POS[i] > bed_df_tmp$Final_POS[j]){
        next
      }
    }
  }
  
  #Filtrando el data-frame solo con mutaciones con gen asignado si filter=TRUE
  if (filter == TRUE){
    vcf_df_chr <- subset(vcf_df_chr, Gene != "")
  }
  
  #Comprobando si las mutaciones con gen asignado se encuentra dentro del ARNm del gen
  cat(sprintf("Comprobando si %d mutaciones con gen asignando se encuentran dentro del ARNm del gen...\n",nrow(subset(vcf_df_chr, Gene != ""))))
  vcf_df_chr$Check_Gene <- ""
  vcf_df_chr <- BSA_SNP_Gene_2(vcf_df = vcf_df_chr,
                               bed_df = bed_df_chr,
                               rna_fasta = rna_fasta,
                               genome_fasta = genome_fasta,
                               chr_list = chr_list,
                               type = type,
                               POS_Start = POS_Start, POS_End = POS_End,
                               position = position, range = range,
                               check_seq_pb = check_seq_pb)
  
  return(vcf_df_chr)
}

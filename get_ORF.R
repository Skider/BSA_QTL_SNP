get_ORF <- function(seq,
                     length_filter = 75
                     ){
  #Librerias
  if (!require("dplyr", quietly = TRUE)){install.packages("dplyr")}
  if (!require("seqinr", quietly = TRUE)){install.packages("seqinr")}
  library(dplyr)
  library(seqinr, include.only = "s2c")
  
  #Comprobacion de la secuencia
  length_seq <- length(seq)
  if(length_seq == 1){
    seq <- seqinr::s2c(seq)
    length_seq <- length(seq)
  }
  
  if (length_seq < (length_filter * 3)){
    stop("La secuencia no puede contenter un ORF del tamanho minimo requerido.")
  } else {
    seq <- tolower(seq)
  }
  
  #Creacion del data-frame de salida
  data_ORF <- data.frame(ORF = 1,
                         POS_Start = 0,
                         POS_End = 0,
                         Dir = "",
                         Length = 0,
                         CDS = "",
                         AA = "")
  
  N_ORF <- 1
  
  #Obtencion de los ORF en la direccion 5->3
  for (i in 1:(length_seq-2)){
    if (paste(seq[i:(i+2)],collapse="") == "atg"){
      for (j in seq(i+3,length_seq,3)){
        if (paste(seq[j:(j+2)],collapse="") == "tga" |
            paste(seq[j:(j+2)],collapse="") == "taa" |
            paste(seq[j:(j+2)],collapse="") == "tag"){
          data_ORF[nrow(data_ORF)+1,] <- c(as.integer(N_ORF),
                                           as.integer(i),
                                           as.integer(j) + 2,
                                           "+",
                                           (as.integer(j) + 2) - as.integer(i) + 1,
                                           paste(seq[i:(j+2)],collapse = ""),
                                           "")
          N_ORF <- N_ORF + 1
          break
        }
        if (j >= (length_seq - 3)){
          data_ORF[nrow(data_ORF)+1,] <- c(as.integer(N_ORF),
                                           as.integer(i),
                                           paste(">",length_seq,sep=""),
                                           "+",
                                           (as.integer(j) + 2) - as.integer(i) + 1,
                                           paste(seq[i:length_seq],collapse = ""),
                                           "")
          N_ORF <- N_ORF + 1
        }
      }
    }
  }
  
  #Obtencion de los ORF en la direccion 3->5
  for (i in 4:(length_seq-2)){
    if (paste(seq[i:(i+2)],collapse="") == "cat"){
      for (j in seq(i,1,-3)){
        if (paste(seq[j:(j+2)],collapse="") == "tca" |
            paste(seq[j:(j+2)],collapse="") == "tta" |
            paste(seq[j:(j+2)],collapse="") == "cta"){
          data_ORF[nrow(data_ORF)+1,] <- c(as.integer(N_ORF),
                                           as.integer(i) + 2,
                                           as.integer(j),
                                           "-",
                                           (as.integer(i) + 2) - as.integer(j) + 1,
                                           reverse_dna(paste(seq[j:(i+2)],collapse = "")),
                                           "")
          N_ORF <- N_ORF + 1
          break
        }
        if (j <= 3){
          data_ORF[nrow(data_ORF)+1,] <- c(as.integer(N_ORF),
                                           as.integer(i),
                                           "<2",
                                           "-",
                                           (as.integer(i) + 2) - as.integer(j) + 1,
                                           reverse_dna(paste(seq[1:(j+2)],collapse = "")),
                                           "")
          N_ORF <- N_ORF + 1
        }
      }
    }
  }

  #Filtrado por longitud
  data_ORF <- data_ORF %>% dplyr::filter(as.integer(Length) > length_filter)
  
  #Filtrado por ORF dentro de otros en el mismo marco de lectura
  data_ORF_POS <- data_ORF %>% dplyr::filter(Dir == "+") %>% dplyr::distinct(POS_End, .keep_all = TRUE)
  data_ORF_NEG <- data_ORF %>% dplyr::filter(Dir == "-")
  data_ORF_NEG <- data_ORF_NEG[order(as.integer(data_ORF_NEG$Length), decreasing = TRUE),]
  data_ORF_NEG <- data_ORF_NEG %>% dplyr::distinct(POS_End, .keep_all = TRUE)
  data_ORF <- rbind(data_ORF_POS,data_ORF_NEG)
  data_ORF <- data_ORF[order(as.integer(data_ORF$Length), decreasing = TRUE),]
  
  for(i in 1:nrow(data_ORF)){
    data_ORF$AA[i] <- translate_CDS(data_ORF$CDS[i])
  }
  return(data_ORF)
}

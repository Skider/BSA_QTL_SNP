#-----------------------------------------------------------------------------
#Modulo 1: Analisis QTL-BSA. Programa "plot_BSA_ctg"
#-----------------------------------------------------------------------------
#La funcion plot_BSA_ctg genera una grafica unica con todos los contigs/scaffolds
#representando los valores de SNP Index y Delta SNP, asi como una linea promedio
#de los valores de SNP.
#El archivo o data frame de entrada debe tener, como minimo, los siguientes datos:
# - Region: Cromosoma [CHROM] y posicion [POS]
# - Lectura de referencia [REF] y alternativa [ALT] para ambos pools o los valores Delta SNP
#-----------------------------------------------------------------------------
plot_BSA_ctg <- function(df,
                         windowStep = "AUTO",
                         plot_width = 1080, plot_height = 1080, text_size = 18,
                         P1_name = "Pool 1", P2_name = "Pool 2",
                         plot_type = "CHROM",
                         INDEL = FALSE,
                         output = "PNG",
                         max_region_graph = 25,
                         LOG = TRUE){
  #Variables de entrada
  #df: Data Frame con los datos
  #windowStep: Tamanho de cada intervalo para el calculo de la linea promedio
  #Informacion para la(s) grafica(s)
  #   plot_width: Anchura de cada grafica
  #   plot_height: Altura de cada grafica
  #   P1_name: Nombre asignado a la primera pool de datos
  #   P2_name: Nombre asignado a la segunda pool de datos
  #   text_size: Tamanho de las letras en la grafica
  #plot_type:
  #   "CHROM" Una grafica para cada cromosoma uniendo los tres SNP Index
  #   "Indiv" Una grafica individual para cada Pools y cromosoma
  #INDEL: Dibuja diferente las variantes que sean INDELs.
  #output:
  #   "PNG" para guardar las graficas en archivos PNG
  #   "plot" para devolver por pantalla el ultimo plot completo
  #   "lplot" para devolver los ultimos plot en una lista
  #---------------------------------------------------------------------------
  #Librerias
  if (!require("ggplot2", quietly = TRUE)){install.packages("ggplot2")}
  if (!require("plotly", quietly = TRUE)){install.packages("plotly")}
  if (!require("gridExtra", quietly = TRUE)){install.packages("gridExtra")}
  if (!require("scales", quietly = TRUE)){install.packages("scales")}
  library(ggplot2)
  library(plotly, include.only = 'subplot')
  library(gridExtra, include.only = 'grid.arrange')
  library(scales, include.only = 'label_number')
  
  vcf_df <- df
  #Comprobacion si estan los datos necesarios para la grafica: SNPIndex.
  #Si no estan los valores SNPIndex calculado, se comprueba si estan los valores de AD y se
  #calculan los valores de SNPIndex.
  if(sum(colnames(vcf_df) == "SNPIndex_P1") != 1 | sum(colnames(vcf_df) == "SNPIndex_P2") != 1 |
     sum(colnames(vcf_df) == "Delta_SNPIndex") != 1){
    if(sum(colnames(vcf_df) == "AD_REF_P1") != 1 | sum(colnames(vcf_df) == "AD_ALT_P1") != 1 |
       sum(colnames(vcf_df) == "AD_REF_P2") != 1 | sum(colnames(vcf_df) == "AD_ALT_P2") != 1){
      stop("El data frame no contiene los datos de frecuencia alelica ni los SNP Index necesarios.")
    } else {
      cat(sprintf("Generando los SNP Index\n"))
      vcf_df <- getSNPIndex(vcf_df)
    }
  }
  
  #Obtencion del vector de los cromosomas
  chr_list <- unique(vcf_df$CHROM)
  chr_list <- chr_list[grepl("chr*",chr_list)]
  
  #Creacion del data frame sin los cromosomas: vcf_df_ctg
  vcf_df_ctg <- subset(vcf_df, !CHROM %in% chr_list)
  ctg_list <- unique(vcf_df_ctg$CHROM)
  cat(sprintf("Numero de contigs/scaffold: %i regiones\n",length(ctg_list)))
  
  #Comprobacion de errores
  if(typeof(vcf_df_ctg) != "list" | class(vcf_df_ctg) != "data.frame" | nrow(vcf_df_ctg) == 0){
    stop("Error en la adquisicion de los datos para la grafica.")
  }
  
  #Creacion de la tabla con las ventanas para los valores promedios de SNP Index
  SNPindex_windows <- data.frame(CHROM = "",
                                 start = 0, end = 2, 
                                 mid = 1,
                                 SNPIndex_P1 = 0.0,
                                 SNPIndex_P2 = 0.0,
                                 Delta_SNPIndex = 0.0)
  
  
  #Generacion del archivo log
  if(LOG){write_lines("Archivo log de plot_BSA_ctg","plot_BSA_ctg.log", append = FALSE)}
  
  check_pc <- 1
  for (i in 1:length(ctg_list)){
    if(length(ctg_list) < 100){
      cat(sprintf("[%i/%i] Generando los datos para la grafica. Region %s\n",i,length(ctg_list),ctg_list[i]))
    } else {
      if(i == as.integer(check_pc * length(ctg_list) * 0.1)){
        cat(sprintf("Generado los datos para el %d%% de las regiones.\n",check_pc*10))
        check_pc <- check_pc + 1
      }
    }
    #Creacion del data frame con los datos de un cromosoma de la lista
    vcf_df_ctg_tmp <- subset(vcf_df_ctg, CHROM == ctg_list[i])
    
    #Si alguna region no tiene mutaciones se salta
    if(nrow(vcf_df_ctg_tmp) == 0){
      next
    }
    
    #Obtencion de la longitud de la region y parametros relacionados con la ventana para la linea promedio
    ctgLength <- max(vcf_df_ctg_tmp$POS)
    POS_Start <- 0
    POS_End <- ctgLength
    if(windowStep == "AUTO"){
      windowStep_tmp <- signif(as.integer(ctgLength)*0.1,1)
    }
    windowSize <- (ctgLength%/%windowStep_tmp)*windowStep_tmp
    windowStart <- seq(from = 1, to = ctgLength, by = windowStep_tmp)
    windowEnd <- seq(from = windowStep_tmp, to = ctgLength, by = windowStep_tmp)

    #Comprobacion de errores
    if(windowSize == 0){
      stop("El tamaño de ventana es nulo. Revise los parametros de entrada.")
    }
    if(length(windowStart) <= 3 | length(windowStart) <= 3){
      stop("Insuficientes puntos para el promedio. Aumenta el tamaño de la ventana o disminuye el intervalo de salto.")
    }
    if(length(windowStart) != length(windowEnd)){
      if(length(windowStart)+1 != length(windowEnd)){
        windowEnd <- c(windowEnd,POS_End)
      } else {
        stop("Error en la adquisicion de la ventana.")
      }
    }
    
    #Creacion de la tabla con las ventanas para los valores promedios de SNP Index para una region
    SNPindex_windows_tmp <- data.frame(CHROM = ctg_list[i],
                                       start = windowStart, end = windowEnd, 
                                       mid = windowStart + (windowEnd-windowStart)/2,
                                       SNPIndex_P1 = replicate(length(windowStart),NA),
                                       SNPIndex_P2 = replicate(length(windowStart),NA),
                                       Delta_SNPIndex = replicate(length(windowStart),NA))
    
    #Calculo de los SNP Index promedios
    for (n in 1:nrow(SNPindex_windows_tmp)) {
      #Data frame temporal con los datos a promediar para una ventana
      vcf_df_window <- vcf_df_ctg_tmp[which(vcf_df_ctg_tmp$POS >= SNPindex_windows_tmp$start[n] &
                                            vcf_df_ctg_tmp$POS <= SNPindex_windows_tmp$end[n]),]
      
      #Calculo del promedio de los SNPIndex para cada ventana
      mean_SNPIndex_P1 <- mean(vcf_df_window$SNPIndex_P1[vcf_df_window$SNPIndex_P1>0],na.rm = TRUE)
      mean_SNPIndex_P2 <- mean(vcf_df_window$SNPIndex_P2[vcf_df_window$SNPIndex_P2>0],na.rm = TRUE)
      mean_Delta_SNPIndex <- mean(vcf_df_window$Delta_SNPIndex,na.rm = TRUE)
      
      #Remplazo de los valores de NA por 0
      if (is.nan(mean_SNPIndex_P1)) {mean_SNPIndex_P1 <- 0}
      if (is.nan(mean_SNPIndex_P2)) {mean_SNPIndex_P2 <- 0}
      if (is.nan(mean_Delta_SNPIndex)) {mean_Delta_SNPIndex <- 0}
      
      #Guardado de los datos en el data frame
      SNPindex_windows_tmp$SNPIndex_P1[n] <- mean_SNPIndex_P1
      SNPindex_windows_tmp$SNPIndex_P2[n] <- mean_SNPIndex_P2
      SNPindex_windows_tmp$Delta_SNPIndex[n] <- mean_Delta_SNPIndex
    }
    #Comprobacion de posibles errores. Se guardan en un fichero si LOG es  TRUE
    if (LOG){
      if(nrow(subset(SNPindex_windows_tmp, SNPIndex_P1 == 0)) == nrow(SNPindex_windows_tmp)){
        write_lines(sprintf("Todos los valores promedios son 0 para el pool 1 en la region %s.\n",ctg_list[i]),
                    "plot_BSA_ctg.log", append = T)
      } else if(nrow(subset(SNPindex_windows_tmp, SNPIndex_P2 == 0)) == nrow(SNPindex_windows_tmp)){
        write_lines(sprintf("Todos los valores promedios son 0 para el pool 2 en la region %s.\n",ctg_list[i]),
                    "plot_BSA_ctg.log", append = T)
      }
    } else {
      if(nrow(subset(SNPindex_windows_tmp, SNPIndex_P1 == 0)) == nrow(SNPindex_windows_tmp)){
        cat(sprintf("Todos los valores promedios son 0 para el pool 1 en la region %s.\n",ctg_list[i]))
      } else if(nrow(subset(SNPindex_windows_tmp, SNPIndex_P2 == 0)) == nrow(SNPindex_windows_tmp)){
        cat(sprintf("Todos los valores promedios son 0 para el pool 2 en la region %s.\n",ctg_list[i]))
      }
    }
    SNPindex_windows <- rbind(SNPindex_windows,SNPindex_windows_tmp)
  }
  
  SNPindex_windows <- SNPindex_windows[-1,]
  
  #Probando con 25 regiones
  if(length(ctg_list) > max_region_graph){
    
  }
  ctg_list <- ctg_list[1:25]
  SNPindex_windows <- subset(SNPindex_windows, CHROM %in% ctg_list)
  vcf_df_ctg <- subset(vcf_df_ctg, CHROM %in% ctg_list)
  
  p_SNPIndex_P1 <- ggplot() + 
    geom_point(data=vcf_df_ctg,aes(x=POS,y=SNPIndex_P1),color="darkgreen") + 
    geom_line(data=SNPindex_windows,aes(x=mid,y=SNPIndex_P1), color="red", size=1.5) +
    scale_x_continuous(name="Chromosome Position [kb]",labels = scales::label_number(scale = 1e-3)) +
    scale_y_continuous(name="SNP-Index", limits = c(0,1)) +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size=as.numeric(text_size)), legend.position = "none")  + 
    ggtitle(P1_name) +
    facet_wrap(vars(CHROM), nrow = 1, scales = "free")

  p_SNPIndex_P2 <- ggplot() + 
    geom_point(data=vcf_df_ctg,aes(x=POS,y=SNPIndex_P2),color="orange") + 
    geom_line(data=SNPindex_windows,aes(x=mid,y=SNPIndex_P2), color="red", size=1.5) +
    scale_x_continuous(name="Chromosome Position [kb]",labels = scales::label_number(scale = 1e-3)) +
    scale_y_continuous(name="SNP-Index", limits = c(0,1)) +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size=as.numeric(text_size)), legend.position = "none")  + 
    ggtitle(P2_name) +
    facet_wrap(vars(CHROM), nrow = 1, scales = "free")
  
  p_Delta_SNPIndex <- ggplot() + 
    geom_point(data=vcf_df_ctg,aes(x=POS,y=Delta_SNPIndex),color="blue") + 
    geom_line(data=SNPindex_windows,aes(x=mid,y=Delta_SNPIndex), color="red", size=1.5) +
    scale_x_continuous(name="Chromosome Position [kb]",labels = scales::label_number(scale = 1e-3)) +
    scale_y_continuous(name="SNP-Index", limits = c(-1,1)) +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size=as.numeric(text_size)), legend.position = "none")  + 
    ggtitle("Delta") +
    facet_wrap(vars(CHROM), nrow = 1, scales = "free")
  

  png("plot_SNPIndex_ctg.png", width=as.numeric(plot_width*25), height=as.numeric(plot_height))
  gridExtra::grid.arrange(p_SNPIndex_P1, p_SNPIndex_P2, p_Delta_SNPIndex, ncol = 1)
  dev.off()

  #return()
}

#-----------------------------------------------------------------------------
#Modulo 1: Analisis QTL-BSA. Programa "plot_BSA"
#-----------------------------------------------------------------------------
#La funcion plot_BSA generea graficas de los valores de SNP Index entre dos pools
#y del Delta SNP entre ambos, asi como una linea promedio de los valores de SNP en cada
#grafica para cada una de las regiones del genoma: cromosomas y contigs/scaffold.
#Adicionalmente, se puede generar un archivo de texto con los valores de la linea promedio
#El archivo o data frame de entrada debe tener, como minimo, los siguientes datos:
# - Region: Cromosoma [CHROM] y posicion [POS]
# - Lectura de referencia [REF] y alternativa [ALT] para ambos pools o los valores Delta SNP
#-----------------------------------------------------------------------------
plot_BSA <- function(df,
                     chr_list = "CHROM",
                     type = "Completo",
                     POS_Start = 1, POS_End = 1000000,
                     position = 1000000, range = 500000,
                     windowStep = 100000,
                     plot_width = 1080, plot_height = 1080, text_size = 18,
                     P1_name = "Pool 1", P2_name = "Pool 2",
                     plot_type = "CHROM",
                     INDEL = FALSE,
                     output = "PNG"){
  #Variables de entrada
  #df: Data Frame con los datos
  #chr_list: Lista de cromosomas a analizar:
  #   La opcion por defecto: "CHROM", utiliza solo los cromosomas
  #   "Todos", utiliza todos los valores en la columna del data frame de entrada
  #Eleccion del tipo de grafica o rango utilizado:
  #   type: Completo: Todo el cromosoma. Eleccion por defecto.
  #         Rango: Rango fijado entre POS_Start y POS_End. Por defecto entre 1-1.000.000
  #         Punto: Rango en un 'range' a la derecha y a la izquierda de 'position'.
  #                 Por defecto en la posicion 1.000.000 con un rango de mas-menos 500.000
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
  if(chr_list[1] == "Todos"){
    chr_list <- unique(vcf_df$CHROM)
    cat(sprintf("Numero de entradas obtenidas: %i entradas\n",length(chr_list)))
  } else if(chr_list[1] == "CHROM"){
    chr_list <- unique(vcf_df$CHROM)
    chr_list <- chr_list[grepl("chr*",chr_list)]
    cat(sprintf("Numero de cromosomas obtenidos: %i cromosomas\n",length(chr_list)))
  }
  
  #Creacion del data frame solo con los cromosomas escogidos: vcf_df_chr
  vcf_df_chr <- subset(vcf_df, CHROM %in% chr_list)
  #Comprobacion de errores
  if(typeof(vcf_df_chr) != "list" | class(vcf_df_chr) != "data.frame" | nrow(vcf_df) == 0){
    stop("Error en la adquisicion de los datos para la grafica.")
  }
  
  for (i in 1:length(chr_list)){
    cat(sprintf("[%i/%i] Creando la grafica para la entrada %s\n",i,length(chr_list),chr_list[i]))
    #Creacion del data frame con los datos de un cromosoma de la lista
    vcf_df_chr_tmp <- subset(vcf_df_chr, CHROM == chr_list[i])
    
    #Obtencion de la longitud del cromosoma
    chrLength <- max(vcf_df_chr_tmp$POS)
    
    #Obtencion de los valores de posicion y el tamanho de la ventana de analisis
    if(type == "Completo"){
      #Para cromosomas completos
      windowSize <- (chrLength%/%windowStep)*windowStep
      windowStart <- seq(from = 1, to = chrLength, by = windowStep)
      windowEnd <- seq(from = windowStep, to = chrLength, by = windowStep)
      POS_Start <- 0
      POS_End <- chrLength
    } else {
      #Para region especifica
      if(type == "Punto"){
        POS_Start <- position - range
        POS_End <- position + range
        if(POS_Start <= 0){POS_Start = 1}
        if(POS_End > chrLength){POS_End = chrLength}
      }
      #Para rango
      if(type == "Rango"){
        windowSize <- POS_End - POS_Start
        windowStart <- seq(from = POS_Start, to = POS_End, by = windowStep)
        windowEnd <- seq(from = POS_Start+windowStep, to = POS_End, by = windowStep)
      }
    }
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
    
    #Creacion de la tabla con las ventanas para los valores promedios de SNP Index
    SNPindex_windows <- data.frame(start = windowStart, end = windowEnd, 
                                   mid = windowStart + (windowEnd-windowStart)/2,
                                   SNPIndex_P1 = replicate(length(windowStart),NA),
                                   SNPIndex_P2 = replicate(length(windowStart),NA),
                                   Delta_SNPIndex = replicate(length(windowStart),NA))
    
    #Calculo de los SNP Index promedios
    for (n in 1:nrow(SNPindex_windows)) {
      #Data frame temporal con los datos a promediar para una ventana
      vcf_df_window <- vcf_df_chr_tmp[which(vcf_df_chr_tmp$POS >= SNPindex_windows$start[n] & vcf_df_chr_tmp$POS <= SNPindex_windows$end[n]),]
      
      #Calculo del promedio de los SNPIndex para cada ventana
      mean_SNPIndex_P1 <- mean(vcf_df_window$SNPIndex_P1[vcf_df_window$SNPIndex_P1>0],na.rm = TRUE)
      mean_SNPIndex_P2 <- mean(vcf_df_window$SNPIndex_P2[vcf_df_window$SNPIndex_P2>0],na.rm = TRUE)
      mean_Delta_SNPIndex <- mean(vcf_df_window$Delta_SNPIndex,na.rm = TRUE)
      
      #Remplazo de los valores de NA por 0
      if (is.nan(mean_SNPIndex_P1)) {mean_SNPIndex_P1 <- 0}
      if (is.nan(mean_SNPIndex_P2)) {mean_SNPIndex_P2 <- 0}
      if (is.nan(mean_Delta_SNPIndex)) {mean_Delta_SNPIndex <- 0}
      
      #Guardado de los datos en el data frame
      SNPindex_windows$SNPIndex_P1[n] <- mean_SNPIndex_P1
      SNPindex_windows$SNPIndex_P2[n] <- mean_SNPIndex_P2
      SNPindex_windows$Delta_SNPIndex[n] <- mean_Delta_SNPIndex
    }
    #Comprobacion de errores
    if(nrow(subset(SNPindex_windows, SNPIndex_P1 == 0)) == nrow(SNPindex_windows)){
      stop("Todos los valores promedios son 0 para el pool 1.")
    } else if(nrow(subset(SNPindex_windows, SNPIndex_P2 == 0)) == nrow(SNPindex_windows)){
      stop("Todos los valores promedios son 0 para el pool 2.")
    }
    
    #Grafico de los resultados
    if(INDEL == FALSE){
      p_SNPIndex_P1 <- ggplot() +
        geom_point(data=subset(vcf_df_chr_tmp, POS >= POS_Start & POS <= POS_End),
                   aes(x=POS,y=SNPIndex_P1), color="darkgreen") + #Todos los SNP Index por posicion
        geom_line(data=SNPindex_windows,
                  aes(x=mid,y=SNPIndex_P1), color="red", size = 2) + #Linea promedio de SNP Index
        scale_x_continuous(name="Chromosome Position [Mb]",labels = scales::label_number(scale = 1e-6)) +
        scale_y_continuous(name="SNP-Index", limits = c(0,1)) +
        ggtitle(sprintf("%s. %s",P1_name,chr_list[i])) +
        theme(plot.title = element_text(hjust = 0.5), text = element_text(size=as.numeric(text_size)))
      p_SNPIndex_P2 <- ggplot() +
        geom_point(data=subset(vcf_df_chr_tmp, POS >= POS_Start & POS <= POS_End),
                   aes(x=POS,y=SNPIndex_P2), color="orange") + #Todos los SNP Index por posicion
        geom_line(data=SNPindex_windows,
                  aes(x=mid,y=SNPIndex_P2), color="red",size = 2) + #Linea promedio de SNP Index
        scale_x_continuous(name="Chromosome Position [Mb]",labels = scales::label_number(scale = 1e-6)) +
        scale_y_continuous(name="SNP-Index", limits = c(0,1)) +
        ggtitle(sprintf("%s. %s",P2_name,chr_list[i])) +
        theme(plot.title = element_text(hjust = 0.5), text = element_text(size=as.numeric(text_size)))
      p_Delta_SNPIndex <- ggplot() +
        geom_point(data=subset(vcf_df_chr_tmp, POS >= POS_Start & POS <= POS_End),
                   aes(x=POS,y=Delta_SNPIndex), color="blue") + #Todos los SNP Index por posicion
        geom_line(data=SNPindex_windows,
                  aes(x=mid,y=Delta_SNPIndex), color="red", size=2) + #Linea promedio de SNP Index
        scale_x_continuous(name="Chromosome Position [Mb]",labels = scales::label_number(scale = 1e-6)) +
        scale_y_continuous(name="SNP-Index", limits = c(-1,1)) +
        ggtitle(sprintf("Delta. %s",chr_list[i])) +
        theme(plot.title = element_text(hjust = 0.5), text = element_text(size=as.numeric(text_size)))
    } else {
      p_SNPIndex_P1 <- ggplot() +
        geom_point(data=subset(vcf_df_chr_tmp, POS >= POS_Start & POS <= POS_End),
                   aes(x=POS,y=SNPIndex_P1,color=VARIANT,shape=VARIANT)) + #Todos los SNP Index por posicion
        scale_color_manual(values = c("green","darkgreen")) +
        geom_line(data=SNPindex_windows,
                  aes(x=mid,y=SNPIndex_P1), color="red", size = 1) + #Linea promedio de SNP Index
        scale_x_continuous(name="Chromosome Position [Mb]",labels = scales::label_number(scale = 1e-6)) +
        scale_y_continuous(name="SNP-Index", limits = c(0,1)) +
        ggtitle(sprintf("%s. %s",P1_name,chr_list[i])) +
        theme(plot.title = element_text(hjust = 0.5), text = element_text(size=as.numeric(text_size)))
      p_SNPIndex_P2 <- ggplot() +
        geom_point(data=subset(vcf_df_chr_tmp, POS >= POS_Start & POS <= POS_End),
                   aes(x=POS,y=SNPIndex_P2,color=VARIANT,shape=VARIANT)) + #Todos los SNP Index por posicion
        scale_color_manual(values = c("yellow","orange")) +
        geom_line(data=SNPindex_windows,
                  aes(x=mid,y=SNPIndex_P2), color="red",size = 1) + #Linea promedio de SNP Index
        scale_x_continuous(name="Chromosome Position [Mb]",labels = scales::label_number(scale = 1e-6)) +
        scale_y_continuous(name="SNP-Index", limits = c(0,1)) +
        ggtitle(sprintf("%s. %s",P2_name,chr_list[i])) +
        theme(plot.title = element_text(hjust = 0.5), text = element_text(size=as.numeric(text_size)))
      p_Delta_SNPIndex <- ggplot() +
        geom_point(data=subset(vcf_df_chr_tmp, POS >= POS_Start & POS <= POS_End),
                   aes(x=POS,y=Delta_SNPIndex,color=VARIANT,shape=VARIANT)) + #Todos los SNP Index por posicion
        scale_color_manual(values = c("lightblue","blue")) +
        geom_line(data=SNPindex_windows,
                  aes(x=mid,y=Delta_SNPIndex), color="red", size=1) + #Linea promedio de SNP Index
        scale_x_continuous(name="Chromosome Position [Mb]",labels = scales::label_number(scale = 1e-6)) +
        scale_y_continuous(name="SNP-Index", limits = c(-1,1)) +
        ggtitle(sprintf("Delta. %s",chr_list[i])) +
        theme(plot.title = element_text(hjust = 0.5), text = element_text(size=as.numeric(text_size)))
    }
    #Guardado de las graficas en un png
    if(plot_type == "CHROM"){
      png(sprintf("plot_SNPIndex_%s_%i_%i.png",chr_list[i],POS_Start,POS_End), width=as.numeric(plot_width), height=as.numeric(plot_height))
      gridExtra::grid.arrange(p_SNPIndex_P1, p_SNPIndex_P2, p_Delta_SNPIndex, ncol = 1)
      dev.off()
    } else if(plot_type == "Indiv") {
      png(sprintf("plot_SNPIndex_%s_%s_%i_%i.png",P1_name,chr_list[i],POS_Start,POS_End), width=as.numeric(plot_width), height=as.numeric(plot_height))
      p_SNPIndex_P1
      png(sprintf("plot_SNPIndex_%s_%s_%i_%i.png",P2_name,chr_list[i],POS_Start,POS_End), width=as.numeric(plot_width), height=as.numeric(plot_height))
      p_SNPIndex_P2
      png(sprintf("plot_Delta_SNPIndex_%s_%i_%i.png",chr_list[i],POS_Start,POS_End), width=as.numeric(plot_width), height=as.numeric(plot_height))
      p_Delta_SNPIndex
      dev.off()
    }
  }
  if(output == "plot"){
    return(plotly::subplot(p_SNPIndex_P1,p_SNPIndex_P2,p_Delta_SNPIndex,nrows=3))
  } else if(output == "lplot"){
    return(list(p_SNPIndex_P1,p_SNPIndex_P2,p_Delta_SNPIndex))
  }
}

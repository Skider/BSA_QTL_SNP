#-----------------------------------------------------------------------------
#Modulo 1: Analisis QTL-BSA. Programa "getSNPIndex"
#-----------------------------------------------------------------------------
#La funcion getSNPIndex calcula los valores de SNP Index para los dos pools,
#asi como el valor del Delta SNP Index de estas Pools para un analisis BSA
#El archivo o data frame de entrada debe tener, como minimo, los siguientes datos
#para cada mutacion:
# - AD_REF_P1. Lectura de la base de referencia de la Pool 1
# - AD_ALT_P1. Lectura de la base de alternativa de la Pool 1
# - AD_REF_P2. Lectura de la base de referencia de la Pool 2
# - AD_ALT_P2. Lectura de la base de alternativa de la Pool 2
#NOTA: Tener en cuenta que el delta SNP Index sera la pool 1 menos la pool 2
#-----------------------------------------------------------------------------
#La funcion devuelve el data frame de entrada con los SNP Index en tres columnas nuevas
# - SNPIndex_P1. SNP Index para la Pool 1
# - SNPIndex_P2. SNP Index para la Pool 2
# - Delta_SNPIndex. Delta SNP Index de la Pool 1 menos la Pool 2
#-----------------------------------------------------------------------------
getSNPIndex <- function(df){
  #Variables de entrada
    #df: Data Frame con los datos
  #---------------------------------------------------------------------------
  vcf_df <- df
  if(sum(colnames(vcf_df) == "AD_REF_P1") != 1 | sum(colnames(vcf_df) == "AD_ALT_P1") != 1 |
     sum(colnames(vcf_df) == "AD_REF_P2") != 1 | sum(colnames(vcf_df) == "AD_ALT_P2") != 1){
    stop("El data frame no contiene los datos de frecuencia alelica necesarios.")
  }
  
  #Obtencion de los SNP Index para las dos Pools.
  vcf_df$SNPIndex_P1 <- vcf_df$AD_ALT_P1/(vcf_df$AD_REF_P1 + vcf_df$AD_ALT_P1)
  vcf_df$SNPIndex_P2 <- vcf_df$AD_ALT_P2/(vcf_df$AD_REF_P2 + vcf_df$AD_ALT_P2)
  vcf_df$SNPIndex_P1[is.na(vcf_df$SNPIndex_P1)] = 0
  vcf_df$SNPIndex_P2[is.na(vcf_df$SNPIndex_P2)] = 0
  vcf_df$Delta_SNPIndex <- vcf_df$SNPIndex_P1 - vcf_df$SNPIndex_P2
  
  #Chequeo de errores
  if(sum(is.na(vcf_df[,c("SNPIndex_P1","SNPIndex_P2","Delta_SNPIndex")])) != 0){
    stop("Error en el calculo de los SNP Index")
  }
  
  return(vcf_df)
}
  

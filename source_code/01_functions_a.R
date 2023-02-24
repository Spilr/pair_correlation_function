read_mpx <- function(sample_group, id){
  df <- read.table(file = sample_group[["sample_id"]][id],
                   sep = "\t",
                   check.names = F,
                   quote = "",
                   fill = T,
                   header = T)
  df
}

InterRDF <- function(g1, g2, range, dr=1){
  library(raster)
  library(plyr)
  #g1: 1st group of cells
  #g2: 2nd group of cells
  #range: c(xmin,xmax,ymin,ymax) coordinate range of the rectangle space
  #dr: increment for increasing radius
  
  rdf <- list()
  
  #iterate through list of center point cells
  for(i in 1:nrow(g1)){
    
    R <- min(c(abs(g1[i,1]-range[1:2]),abs(g1[i,2]-range[3:4])))
    if(R > dr){
      dis_n <- pointDistance(g1[i,,drop=F], g2, lonlat=F, allpairs=T)
      dis_n <- dis_n[dis_n <= R]
      
      breaks <- seq(0, R-dr, by=dr)
      obs <- table(cut(dis_n, c(breaks,breaks[length(breaks)]+dr)))
      exp <- sapply(breaks,function(x) ((x+dr)^2-x^2)/R^2*length(dis_n))
      
      rdf[[i]] <- obs %/% exp
    }
  }
  rdf <- ldply(rdf, rbind)
  
  rdf_summary <- as.data.frame(cbind(
    seq(0, ncol(rdf)-dr, by=dr),
    apply(rdf,2,function(x) sd(x,na.rm=T)),
    apply(rdf,2,function(x) mean(x,na.rm=T))
  ))
  names(rdf_summary) <- c("radius","mean","sd")
  
  rdf_summary
}

pcf_tam <- function(df1){
  library(dplyr)
  
  # Define B Cell Population
  B_CELLS <- df1[df1$`Phenotype CD20`=="CD20+",c("Cell X Position", "Cell Y Position")]
  
  # Define TAM Population
  TAM <- df1 %>%
    filter(
      df1$`Phenotype CD163`=="CD163+",
      df1$`Phenotype TMEM119` == "TMEM119-"
      ) %>%
    dplyr::select("Cell X Position", "Cell Y Position")
  
  
  range <- c(min(df1$`Cell X Position`),
             max(df1$`Cell X Position`),
             min(df1$`Cell Y Position`),
             max(df1$`Cell Y Position`))
  
  # Calculate the pair correlation function, set NA values to zero
  rdf <- InterRDF(B_CELLS, TAM, range = range, dr=1)
  rdf[is.na(rdf)]=0
  rdf
}

pcf_tumor <- function(df1){
  library(dplyr)
  # Define B Cell Population
  B_CELLS <- df1[df1$`Phenotype CD20`=="CD20+",c("Cell X Position","Cell Y Position")]
  
  # Define Tumor Population
  TUMOR <- df1 %>%
    filter(
      df1$`Phenotype SOX2`=="SOX2+"
      ) %>%
    dplyr::select("Cell X Position", "Cell Y Position")
  
  
  range <- c(min(df1$`Cell X Position`),
             max(df1$`Cell X Position`),
             min(df1$`Cell Y Position`),
             max(df1$`Cell Y Position`))
  
  # Calculate the pair correlation function, set NA values to zero
  rdf <- InterRDF(B_CELLS, TUMOR, range = range, dr=1)
  rdf[is.na(rdf)]=0
  rdf
}

pcf_tcell <- function(df1){
  # Define B Cell Population
  B_CELLS <- df1[df1$`Phenotype CD20`=="CD20+",c("Cell X Position","Cell Y Position")]
  
  # Define T Cell Population
  T_CELLS <- df1 %>%
    filter(
      df1$`Phenotype CD8`=="CD8+"
      ) %>%
    dplyr::select("Cell X Position", "Cell Y Position")
  
  
  range <- c(min(df1$`Cell X Position`),
             max(df1$`Cell X Position`),
             min(df1$`Cell Y Position`),
             max(df1$`Cell Y Position`))
  
  # Calculate the pair correlation function, set NA values to zero
  rdf <- InterRDF(B_CELLS, T_CELLS, range = range, dr=1)
  rdf[is.na(rdf)]=0
  rdf
}

pcf_mcg <- function(df1){
  library(dplyr)
  # Define B Cell Population
  B_CELLS <- df1[df1$`Phenotype CD20`=="CD20+", c("Cell X Position","Cell Y Position")]
  
  
  # Define Microglia Population
  MCG <- df1 %>%
    filter(
      df1$`Phenotype CD163`=="CD163-",
      df1$`Phenotype TMEM119` == "TMEM119+"
      ) %>%
    dplyr::select("Cell X Position", "Cell Y Position")
  
  
  range <- c(min(df1$`Cell X Position`),
             max(df1$`Cell X Position`),
             min(df1$`Cell Y Position`),
             max(df1$`Cell Y Position`))
  
  # Calculate the pair correlation function, set NA values to zero
  rdf <- InterRDF(B_CELLS, MCG, range = range, dr=1)
  rdf[is.na(rdf)]=0
  rdf
}

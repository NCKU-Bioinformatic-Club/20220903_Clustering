##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Current path and new folder setting* #####
  ProjectName = "CC"
  Sampletype = "DNAMeth"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)
  
  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }
  
##### Load Packages #####
  Package.set <- c("tidyverse","circlize","ComplexHeatmap","stringr")
  ## Check whether the installation of those packages is required from basic
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  # devtools::install_github("jokergoo/circlize")
  
##### load sample #####
  Input_1.df <- read.delim("cachexia_heatmap/cachexia_12_samples_285K_ANOVA.txt", sep = "\t")
  Input_2.df <- read.delim("cachexia_heatmap/cachexia 15 samples_285K _raw data (m-value).txt", sep = "\t")
  group.df <- read.delim("cachexia_heatmap/muscle group information.txt", sep = "\t")

##### Data preprocessing #####
  ## Matrix
  result_1.df <- Input_1.df %>% select(Probeset.ID,CHR,MAPINFO,p.value.Severe.vs..Mild.,
                                       Difference.Severe.vs..Mild.,Beta.Difference.Severe.vs..Mild.) %>% 
                                       filter(abs(Beta.Difference.Severe.vs..Mild.) > 0.1) %>%
                                       filter(p.value.Severe.vs..Mild. < 0.05) %>% 
                                       rename(ID=Probeset.ID) %>% 
                                       left_join(Input_2.df, by="ID")  %>% 
                                       distinct(ID,.keep_all=T)
  ## Extract matrix
  matrix.df <- result_1.df[,-2:-6]
  
  ## Convert "?" to NA
  matrix.df[matrix.df=="?"] <- NA
  
  ## Convert dataframe element to numeric
  matrix.df <- data.frame(apply(matrix.df, 2, function(x) as.numeric(as.character(x))))
  
  ## Annotation col
  anno_colum.df <- data.frame(Sentrix.Barcode=colnames(Input_2.df)[-1]) %>% 
                              left_join(group.df,by="Sentrix.Barcode")
  ## Annotation row
  anno_row.df <- result_1.df[,1:6]

  
##### Export data #####  
  write.table(matrix.df,"matrix.tsv",row.names = F,sep = "\t")
  write.table(anno_colum.df,"annotation.tsv",row.names = F,sep = "\t")


##### Heatmap plotting #####
  ## Set column annotation
  column_ha_T = HeatmapAnnotation(
    Condition = anno_colum.df$Body.weight,
    Condition2 = anno_colum.df$Area,
    col = list(Condition = c("Mild"="#9b6ab8", "Severe"="#6e6970"),
               Condition2= c("Mild"="#9b6ab8", "Severe"="#6e6970","Medium"="#b57545")),
    show_legend = T
  )
  
  ## generate color of top annotation 
  col_exp <-  colorRamp2(
    c(0, 0.025, 0.05),
    c("#248a5c", "#52bf8e","#bbedd7")
  ) 
  col_exp2 <-  colorRamp2(
    c(-17, 0, 17),
    c("#1a5691", "#96cbff", "#d1e8ff")
  ) 
  
  row_ha = rowAnnotation(
    p.value = anno_row.df$p.value.Severe.vs..Mild.,
    LogFC = anno_row.df$Difference.Severe.vs..Mild.,
    col = list(p.value = col_exp, LogFC = col_exp2),
    show_legend = T
  )
  
  
  Heatmap(
    matrix.df[-1],
    # column_title = target_gene,
    # column_title_side = "top",
    cluster_rows = T,
    cluster_columns = T,
    show_column_names = F,
    show_row_names = F,
    name = "M-Value",
    # set color
    col = colorRamp2(
      c(0, 0.5, 1),
      c("#1c77d9", "#1a2938", "#ffe182")
    ),
    show_heatmap_legend = T,
    use_raster = F,
    top_annotation = column_ha_T,
    right_annotation = row_ha
  ) -> P.Heatmap
  
  P.Heatmap %>% print
  
  
##### Export PDF #####
  
  pdf(
    file = paste0(getwd(), "/",Version,"/", Sys.Date(), "_MValue_Heatmap.pdf"),
    width = 7, height = 7
  )
  P.Heatmap
  
  graphics.off()
  
  
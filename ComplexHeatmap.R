##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Current path and new folder setting* #####
  ProjectName = "TCGA" #*
  Sampletype = "BRCA" #*
  ExportName <- "HeatmapTest" #*
 
  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype,"_",ExportName)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){dir.create(Save.Path)}


##### Import setting and Import* #####
  ## File setting*
  InFOLName_GE <- "Input_TCGA"  # Input Folder Name #*
  SampleName <- "Xena_TCGA_BRCA_GE" #*
  SamplePhenoName <- "TCGA.BRCA.sampleMap_BRCA_clinicalMatrix" #*
  
  ## Import genetic data file
  GeneExp.df <- read.table(paste0(InFOLName_GE,"/",SampleName), header=T, row.names = 1, sep="\t")
  # GeneExp.df <- read.table(paste0("D:/Dropbox/##_GitHub/#_NCKU_Bioinformatic_Club/20220903_Clustering/Input_TCGA/Xena_TCGA_BRCA_GE"), header=T, row.names = 1, sep="\t") #*
  
  ## Rename the colnames
  colnames(GeneExp.df) <-  gsub("\\.", "-", colnames(GeneExp.df))
  GeneExp_ORi.df <- GeneExp.df # Save Ori
  
  Anno.df <- read.table(paste0(InFOLName_GE,"/",SamplePhenoName), header=T, row.names = 1, sep="\t")
  # Anno.df <- read.table(paste0("D:/Dropbox/##_GitHub/#_NCKU_Bioinformatic_Club/20220903_Clustering/Input_TCGA/TCGA.BRCA.sampleMap_BRCA_clinicalMatrix"), header=T, row.names = 1, sep="\t") #*
  
  ## Rename the sample type
  # Anno.df$sample_type <-  gsub(" ", "", Anno.df$sample_type)
  
##### Conditions setting* #####
  ## Set grouping mode
  Group_Mode <- "GoupByPheno"   # c("GoupByPheno","GoupByGeneExp") #*
  
  # for GoupByPheno
  PhenoGroupType = "sample_type"
  GroupCompare_Pheno <- c("Primary Tumor","Solid Tissue Normal") #*
  # for GoupByGeneExp
  TarGene_name <- "TOP2A" #*
  GeneExpSet.lt <- list(GeneExpMode = "Mean", # c("Mean","Mean1SD","Mean2SD","Mean3SD","Median","Quartile","Customize"))
                        UpCutoff = 1, LowerCutoff = 1) #*
  
  # Annotation set by previous setting
  if(Group_Mode == "GoupByGeneExp"){
    ## Group by GeneExp
    AnnoSet.lt <- list(GroupType = TarGene_name, GroupCompare = c("High","Low") )   ## DEG by GeneExp group
  }else{
    ## Group by Pheno
    AnnoSet.lt <- list(GroupType = PhenoGroupType, GroupCompare = c(GroupCompare_Pheno) )
  }
  
  ## Set threshold for DEG
  Thr.lt <- list(LogFC = c("logFC",1), pVal = c("PValue",0.05) ) #*
  SampleID = "X_INTEGRATION"
  SampleNum = 100
  GeneNum = 2000
  FDRSet = 0.01
  LogFCSet = 1
  
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
 
##### Function setting #####
  ## Call function 
  source("FUN_Group_GE.R")
  source("FUN_DEG_Analysis.R")
  
##### Data preprocess #####
  Anno_Ori.df <- Anno.df # Save Ori
  colnames(Anno.df)
  
  #### Selection (Optional) ####
  ## Select Pheno column
  PhenoColKeep.set <- c("X_INTEGRATION","histological_type","sample_type","gender") #*
  Anno.df <- Anno.df[,c(PhenoColKeep.set)]
  colnames(Anno.df)
 
  # ## Select Pheno row
  # PhenoRowKeep.set <- list(col="sample_type",row=c("Primary Tumor","Recurrent Tumor"))
  # Anno.df <- Anno.df[Anno.df[,PhenoRowKeep.set[["col"]]] %in% PhenoRowKeep.set[["row"]], ]
  
  
  #### Random select sample ####
  # Extract Group1
  Anno_Grp1.df <- Anno.df[Anno.df[,AnnoSet.lt[["GroupType"]]] %in% AnnoSet.lt[["GroupCompare"]][1], ]
  # Extract Group2
  Anno_Grp2.df <- Anno.df[Anno.df[,AnnoSet.lt[["GroupType"]]] %in% AnnoSet.lt[["GroupCompare"]][2], ]
  
  # Anno.df <- rbind(Anno_Grp2.df[sample(1:nrow(Anno_Grp2.df),nrow(Anno_Grp1.df)),],Anno_Grp1.df)
  Anno.df <- rbind(Anno_Grp2.df[sample(1:nrow(Anno_Grp2.df),SampleNum),],Anno_Grp1.df[sample(1:nrow(Anno_Grp1.df),SampleNum),])
  
  GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% Anno.df[,SampleID]] 
  Anno.df <- Anno.df[Anno.df[,SampleID] %in% colnames(GeneExp.df),]
  
  rm(Anno_Grp1.df,Anno_Grp2.df)
  
  #### Reorder the Anno.df #### ## Important!!!
  # Anno.df <- left_join(data.frame("X_INTEGRATION"=colnames(GeneExp.df)), Anno.df)
  GeneExpCol.df <- colnames(GeneExp.df) %>% as.data.frame()
  colnames(GeneExpCol.df) <- SampleID
  Anno.df <- left_join(GeneExpCol.df, Anno.df)
  rm(GeneExpCol.df)
  
##### Grouping by GeneExp #####
  source("FUN_Group_GE.R")
  ##### Group by gene expression 1: CutOff by total  #####
  GeneExp_group.set <- FUN_Group_GE(GeneExp.df, Anno.df,
                                    TarGeneName = TarGene_name, GroupSet = GeneExpSet.lt,
                                    Save.Path = Save.Path, SampleName = ExportName)
  Anno.df <- GeneExp_group.set[["AnnoNew.df"]]
  rm(GeneExp_group.set)

##### Run DEG #####
  source("FUN_DEG_Analysis.R")
  DEG_ANAL.lt <- FUN_DEG_Analysis(GeneExp.df, Anno.df,
                                  GroupType = AnnoSet.lt[["GroupType"]], GroupCompare = AnnoSet.lt[["GroupCompare"]],
                                  ThrSet = Thr.lt,
                                  TarGeneName = TarGene_name, GroupMode = GeneExpSet.lt, SampleID =  SampleID,
                                  Save.Path = Save.Path, SampleName = ExportName, AnnoName = "")
  DE_Extract.df <- DEG_ANAL.lt[["DE_Extract.df"]]
  
  #### Filter genes ####
  ## Set selectedGenes
  selectedGenes <- DE_Extract.df
  selectedGenes <- selectedGenes[selectedGenes$FDR < FDRSet,]
  selectedGenes <- selectedGenes[abs(selectedGenes$logFC) > LogFCSet,]
  selectedGenes <- selectedGenes[rev(order(abs(selectedGenes$logFC)))[1:GeneNum],]
  # selectedGenes <- selectedGenes[rev(order(selectedGenes$logFC))[1:GeneNum],]

  ## Filter GeneExp matrix by selectedGenes
  matrix.df <- GeneExp.df[row.names(GeneExp.df) %in% selectedGenes$Gene,]
  
  ## Annotation col
  anno_colum.df <- Anno.df
  rm(Anno.df)
  
  ## Annotation row
  anno_row.df <- DE_Extract.df[DE_Extract.df$Gene %in% selectedGenes$Gene, ]


##### Heatmap plotting #####
  ## Set column annotation
  ha_column_T = HeatmapAnnotation(
    Sample = anno_colum.df$sample_type,
    Gender = anno_colum.df$gender,
    TarGene = anno_colum.df[,TarGene_name],
    col = list(Sample = c("Primary Tumor"="#9b6ab8", "Solid Tissue Normal"="#6e6970"),
               Gender = c("MALE"="#4382b5", "FEMALE"="#c25988"), #,"Medium"="#b57545"
               TarGene = c("High"="#db8051", "Low"="#c26334")), # #b6d4ca
    show_legend = T
  )
  
  
  ## Set row annotation
  ## Color setting
  col_exp <-  colorRamp2(
    c(min(anno_row.df$PValue), mean(anno_row.df$PValue), max(anno_row.df$PValue)),
    c("#3f705a", "#52bf8e","#b6d4ca")

  ) 
  col_exp2 <-  colorRamp2(
    c(min(anno_row.df$logFC), mean(anno_row.df$logFC), max(anno_row.df$logFC)),
    c("#488c67", "#333333","#edd493")
  ) 
  
  ha_row = rowAnnotation(
    p.value = anno_row.df$PValue,
    LogFC = anno_row.df$logFC,
    col = list(p.value = col_exp, LogFC = col_exp2),
    show_legend = T
  )
  
  
  Heatmap(
    matrix.df,
    # column_title = target_gene,
    # column_title_side = "top",
    cluster_rows = T,
    cluster_columns = T,
    show_column_names = F,
    show_row_names = F,
    name = "GeneExp",
    # set color
    col = colorRamp2(
      c(min(matrix.df), matrix.df %>% unlist() %>% mean() , max(matrix.df)),
      c("#416db0", "#1a2938", "#bf627e")
    ),
    show_heatmap_legend = T,
    use_raster = F,
    top_annotation = ha_column_T,
    right_annotation = ha_row
  ) -> P.Heatmap
  
  P.Heatmap %>% print
  
  
  Heatmap(
    matrix.df,
    cluster_rows = F,
    cluster_columns = F,
    show_column_names = F,
    show_row_names = F,
    name = "GeneExp",
    # set color
    col = colorRamp2(
      c(min(matrix.df), matrix.df %>% unlist() %>% mean() , max(matrix.df)),
      c("#416db0", "#1a2938", "#bf627e")
    ),
    show_heatmap_legend = T,
    use_raster = F,
    top_annotation = ha_column_T,
    right_annotation = ha_row
  ) -> P.Heatmap2
  
  P.Heatmap2 %>% print
  
  
  # Reorder
  # https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#row-and_column_orders
  Heatmap(
    matrix.df,
    cluster_rows = T,
    cluster_columns = F,
    column_order = order(anno_colum.df$sample_type),
    show_column_names = F,
    show_row_names = F,
    name = "GeneExp",
    # set color
    col = colorRamp2(
      c(min(matrix.df), matrix.df %>% unlist() %>% mean() , max(matrix.df)),
      c("#416db0", "#1a2938", "#bf627e")
    ),
    show_heatmap_legend = T,
    use_raster = F,
    top_annotation = ha_column_T,
    right_annotation = ha_row
  ) -> P.Heatmap3
  
  P.Heatmap3 %>% print
  
  
##### Export PDF #####
  pdf(
    file = paste0(getwd(), "/",Version,"/", Sys.Date(), "_GeneExp_Heatmap.pdf"),
    width = 12, height = 7
  )
    P.Heatmap
    P.Heatmap2
    P.Heatmap3
  graphics.off()
  
  
##### Export data #####  
  write.table(matrix.df,"matrix.tsv",row.names = F,sep = "\t")
  write.table(anno_colum.df,"annotation.tsv",row.names = F,sep = "\t")
  write.table(anno_row.df,"DEG.tsv",row.names = F,sep = "\t")  
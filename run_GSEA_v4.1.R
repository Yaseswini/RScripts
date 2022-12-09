## R Wrapper to run Gene set enrichment analysis 
## Input :
### gseaExec : executable .jar file for gsea ( can be downloaded from GSEA site )
### PreRankedDf : 2 column data frame with colnames : genes and scores ( scores are ranked )
### gmxFile : gmt/gmx file of genesets ( ex msigdb C2/Hallmark etc )
### analysisLabel : label for analysis 
### analysisDir : Output directory to write the gsea output files 
### nperm : number of permutations 
### setMinSize : minimum geneset size 
### setMaxSize : maximum geneset size ( genesets with more than the `setMaxSize` number of genes are removed ) 


GSEA_PreRanked_v4.1 <- function( gseaExec , PreRankedDf , gmxFile , analysisLabel , analysisDir  , nperm = 1000 , setMinSize , setMaxSize , numPlots = 50  )
{  
   #options( java.parameters = "-Xmx8g" )
   
   rankedFile = paste( analysisLabel , "_RankedList.rnk" , sep = "" ) 
   if( ! is.character(PreRankedDf[,1]) & !is.numeric(PreRankedDf[,2]) )
   {  
      stop( "PreRankedDf should be a 2 column dataframe with first column being genes and second logFC " )
   }else
   {  
      PreRankedDf =  PreRankedDf[ order(  PreRankedDf$scores , decreasing = T ) , ]
      write.table( PreRankedDf , file = rankedFile , sep = "\t" , quote = F , col.names = F , row.names = F )
   }
   
   # ./gsea-cli.sh GSEAPreRanked 
   cmd = paste( gseaExec , " GSEAPreRanked " , sep = ""  )
   
   cmd = paste( cmd , "-gmx " , gmxFile ,  sep = "" ) 
   cmd = paste( cmd , " -rnk ", rankedFile , sep = "" )
   cmd = paste( cmd , " -rpt_label " , analysisLabel , sep = "" )
   cmd = paste( cmd , " -out " , analysisDir , sep = "" )
   cmd = paste( cmd , " -set_max " , setMaxSize , sep = "" )
   cmd = paste( cmd , " -set_min " , setMinSize , sep = "" )
   
   system( cmd , wait = T )
   
   system( paste0( "rm " , rankedFile ) )

}

GSEA_LeadingEdge_v4.1 <- function( gseaExec , gseaAnalysisDIR , gsets = NULL , outDIR  )
{
   cmd = paste( gseaExec , " LeadingEdgeTool " , sep = "" )
   cmd = paste( cmd , "-dir " , gseaAnalysisDIR , sep = "" )
   cmd = paste( cmd , " -enrichment_zip " , gseaAnalysisDIR , sep = "" )
   #cmd = paste( cmd , " -output_file_name " , "xx" , sep = "" )
   #cmd = paste( cmd , " -zip_report " , "True" , sep = "" )
   cmd = paste( cmd , " -out " , outDIR , sep = "" )
   #cmd = paste( cmd , " -run_as_genepattern true " , sep = "" )
   if( ! is.null(gsets) )
   {
      cmd = paste( cmd , " -gsets " , gsets , sep = "" )
   }
   system( cmd )
}

## ------------------------------------------------------------ ## 
##                  Parse one GSEA analysis DIR             ## 
## ------------------------------------------------------------ ## 
## input : gsea analysis directory ; 
##       qvalue threshold , 
##         nes threshold
parse_GSEADir <- function( GSEA.analysisDir , qvalueThresh = 0.05 , nesThresh = 1.2, rename_column = T )
{
   if( is.null( qvalueThresh ) | is.null( nesThresh ) ) { stop( "Specify nes threshold and qvalue threshold" ) }
   #analysisLabel = gsub( "\\.GseaPreranked.*" , "" , basename( GSEA.analysisDir ) )
   analysisLabel = basename( GSEA.analysisDir )

   gseaReportFiles = list.files( GSEA.analysisDir , pattern = "gsea_report_for_na_.*tsv" , full.names =  T )
   gseaReportData = lapply( gseaReportFiles , function( reportFile ) { reportData = read.table( reportFile , head = T , sep = "\t" ) ;
                                                                       reportData = reportData[ , c("NAME","NES","FDR.q.val") ] ;
                                                                       return( reportData ) } )

   gseaReportData = do.call( "rbind" ,  gseaReportData )
   # Extracting significant genesets 
   gseaReportData.siggenesets = gseaReportData %>% dplyr::filter( ! is.na(NES) & abs(as.numeric(NES)) > nesThresh & FDR.q.val < qvalueThresh ) %>% dplyr::select( NAME , NES )
   # Genesets that arent significantly enriched 
   #not.significant.genesets = gseaReportData$NAME[ ! gseaReportData$NAME %in% gseaReportData.siggenesets$NAME ]  
   #
   gseaTab = gseaReportData.siggenesets
   #gseaTab = rbind( gseaReportData.siggenesets , data.frame( NAME = not.significant.genesets , NES = rep( NA , length(not.significant.genesets) ) ) ) 
   if( rename_column )
   {
   colnames(gseaTab)[ colnames(gseaTab) == "NES" ] = analysisLabel
   }
return( gseaTab )
}

## --------------------------------------------------------------------------- ## 
##        Parse GSEA directory with multiple analysis directories          ## 
## --------------------------------------------------------------------------- ## 
## input : gsea output directory with multiple analysis directory          
##          qvalue threshold 
##          nes threshold 
## FUN : Parse GSEA Results and filter for qvalue and nes 
## single dir or set of gsea outdir 
parse_GSEAResults <- function( gseaOutDir , qvalueThresh , nesThresh  )
{  
   gseaAnalysisDirs = list.dirs( gseaOutDir , recursive = F , full.names = T )
   
   gseaTab = data.frame( NAME = character() , NES = numeric() , Label = character() )
   analysisDir.NOsigsets = NULL
   for( analysisDir in gseaAnalysisDirs )
   {  
      if( grepl( "\\.GseaPreranked." , basename(analysisDir) ) )
      {  
         #analysisLabel = gsub( "\\.GseaPreranked.*" , "" , basename(analysisDir) )
         analysisLabel = basename(analysisDir)
         analysisData =  parse_GSEADir( GSEA.analysisDir = analysisDir , nesThresh = nesThresh , qvalueThresh = qvalueThresh )
         if( nrow(analysisData) > 0 )
         {      
                analysisData$Label = analysisLabel
                colnames(analysisData)[ colnames(analysisData) == analysisLabel ] = "NES"
                gseaTab = rbind( gseaTab , analysisData )
         }
         else
         {  
            analysisDir.NOsigsets = c( analysisDir.NOsigsets , analysisLabel )
         }
      }
    }
   gseaTab1 = dcast( gseaTab , NAME~Label , value.var = "NES" )
   gseaTab1[ , analysisDir.NOsigsets ] = NA
 
 return( gseaTab1 )
}


## ------------------------------------------------
## Volcanoplot of gsea results 
## ------------------------------------------------

plot_GSEAResults <- function( GSEA.analysisDir , qvalueThresh = 0.05 , nesThresh = 1.2, topn = 5 , outDIR , outLabel )
{

   if( is.null( qvalueThresh ) | is.null( nesThresh ) ) { stop( "Specify nes threshold and qvalue threshold" ) }
   #analysisLabel = gsub( "\\.GseaPreranked.*" , "" , basename( GSEA.analysisDir ) )
   analysisLabel = basename( GSEA.analysisDir )

   gseaReportFiles = list.files( GSEA.analysisDir , pattern = "gsea_report_for_na_.*tsv" , full.names =  T )
   gseaReportData = lapply( gseaReportFiles , function( reportFile ) { reportData = read.table( reportFile , head = T , sep = "\t" ) ;
                                                                       reportData = reportData[ , c("NAME","NES","FDR.q.val") ] ;
                                                                       return( reportData ) } )

   gseaReportData = do.call( "rbind" ,  gseaReportData )
   # Extracting significant genesets 
   gseaReportData.siggenesets = gseaReportData %>% 
                                       dplyr::filter( ! is.na(NES) & abs(as.numeric(NES)) > nesThresh & FDR.q.val < qvalueThresh ) %>% 
                                       dplyr::select( NAME , NES )

   sigGenesets_topn = gseaReportData.siggenesets %>% arrange(desc(NES)) %>% head( n = topn )
   if( min(sigGenesets_top5$NES) < 0 ){
      sigGenesets_topn = rbind( sigGenesets_topn , gseaReportData.siggenesets %>% arrange(NES) %>% head( n = topn ) )
   }

   ## Plot the genesets 
   gseaReportData = gseaReportData %>% mutate( negLogFDR = -log10(FDR.q.val+0.000001) )

   plotData = gseaReportData %>% 
                     mutate( point_color = if_else( NAME %in% sigGenesets_top5$NAME , "yes" , "no" ) ) %>% 
                     mutate( point_label = if_else( NAME %in% sigGenesets_top5$NAME , NAME , "" ) ) %>%
                     arrange( desc(point_color) )

   xmax = ceiling( max(abs(plotData$NES)) )

   pdf( file.path( outDIR , paste0( outLabel , "GSEAplots.pdf" ) ) )
   p = ggplot( plotData , aes( x = NES , y = negLogFDR , fill = point_color , label = point_label ) ) + geom_point( pch = 21 , color = "black" ) + 
            geom_text_repel(    size = 2 , fontface = 'bold', color = 'black',
                               box.padding = unit(0.35, "lines"),
                               point.padding = unit(0.5, "lines"),
                               segment.color = 'grey50' , max.overlaps = nrow(plotData) )
   p = p + scale_fill_manual( values = c("yes"="red","no"="grey") )
   p = p +  theme_bw() + xlim( -xmax , xmax )
   p = p + theme( axis.text = element_text( size = 14 , face = "bold" ) , 
                                         axis.title = element_text( size = 14 , face = "bold" ) ,
                                         panel.grid.major = element_blank() ,   
                                         panel.grid.minor = element_blank() , 
                                         legend.position = "none" )  
   print(p)
   dev.off()

}


## -------------------- 
## Generate bar plots of NES scores with top n enriched pathways 
## ------------------- 

plot_GSEAResults_as_barPlots <- function( GSEA.analysisDir , qvalueThresh = 0.05 , nesThresh = 1.2, topn = 5 , outDIR , outLabel , posNES_color = "red", negNES_color = "blue" , plotwidth = 7 , plotheight = 7 )
{

   if( is.null( qvalueThresh ) | is.null( nesThresh ) ) { stop( "Specify nes threshold and qvalue threshold" ) }
   #analysisLabel = gsub( "\\.GseaPreranked.*" , "" , basename( GSEA.analysisDir ) )
   analysisLabel = basename( GSEA.analysisDir )

   gseaReportFiles = list.files( GSEA.analysisDir , pattern = "gsea_report_for_na_.*tsv" , full.names =  T )
   gseaReportData = lapply( gseaReportFiles , function( reportFile ) { reportData = read.table( reportFile , head = T , sep = "\t" ) ;
                                                                       reportData = reportData[ , c("NAME","NES","FDR.q.val") ] ;
                                                                       return( reportData ) } )

   gseaReportData = do.call( "rbind" ,  gseaReportData )
   # Extracting significant genesets 
   gseaReportData.siggenesets = gseaReportData %>% 
                                       dplyr::filter( ! is.na(NES) & abs(as.numeric(NES)) > nesThresh & FDR.q.val < qvalueThresh ) %>% 
                                       dplyr::select( NAME , NES )

   gseaReportPosNES = gseaReportData.siggenesets %>% filter( NES > 0 ) %>% arrange( desc(NES) ) %>% head( n = topn )
   gseaReportNegNES = gseaReportData.siggenesets %>% filter( NES < 0 ) %>% arrange( NES ) %>% head( n = topn )



   ## Plot the genesets 
   plotData = rbind( gseaReportPosNES , gseaReportNegNES ) %>% mutate( NES_group = if_else( NES > 0 , "posNES" , "negNES" ) )

   xmax = ceiling( max(abs(plotData$NES)) )
   
   pdf( file.path( outDIR , paste0( outLabel , paste0( "GSEA_top_" , topn , "_significantPathways_barplot.pdf" ) ) ) , width = plotwidth , height = plotheight )
   
   p = ggplot( plotData , aes( x = NES , y = fct_reorder( NAME , NES ) , fill = NES_group ) ) + 
         geom_bar( stat = 'identity' ) +
         scale_fill_manual( values = c("posNES"= posNES_color,"negNES"=negNES_color ) )
   p = p + theme_bw() + xlim( -xmax , xmax )
   p = p + theme( axis.text = element_text( size = 14 , face = "bold" ), 
                  axis.title = element_text( size = 14 , face = "bold" ) , 
                  panel.grid.major = element_blank() , 
                  panel.grid.minor = element_blank() , 
                  legend.position = "none" ) 
   p = p + labs( x = "Negative Enrichment Score(NES)" , y = "" )
   print(p)

   dev.off()

}


## --- given a gsea analysis directory , plot the leading edge genes in top n pathways ( max value for n is 50 which is used to generate geneerate gsea reuslts )

get_leadingEdge_heatmap <- function( gsea_analysisDIR   , outDIR , outLabel ){

   outFilePath = file.path( outDIR , paste0( outLabel , "_LeadingEdge_OverlapMatrix.pdf" ) )

   leadingEdgeFiles = list.files( gsea_analysisDIR , full.names = T , pattern = "*.tsv" )
   names(leadingEdgeFiles) = basename(leadingEdgeFiles)
   leadingEdgeFiles = leadingEdgeFiles[ ! names(leadingEdgeFiles) %in% grep("^gsea_report|ranked|sizes",names(leadingEdgeFiles) , value = T ) ]

   leadingEdgeData = lapply( leadingEdgeFiles , function(fn) read_tsv( fn  ) %>% filter( `CORE ENRICHMENT` == "Yes" ) )
   leadingEdgeData = lapply( seq_along(leadingEdgeData) , function(idx,Lst) { dat = Lst[[idx]] ; 
                                                                                 gn = gsub("*.tsv","",names(Lst)[[idx]]) ; 
                                                                                 dat[ , gn ] = 1
                                                                                 return(dat %>% dplyr::select( SYMBOL , gn )) } , leadingEdgeData )

   LeadingEdge_geneToPathwayTab = data.frame( Reduce( function(...) merge(...,all=T) , leadingEdgeData  ) )
   LeadingEdge_geneToPathwayTab = data.frame( LeadingEdge_geneToPathwayTab ) %>%  column_to_rownames( var = "SYMBOL" )
   LeadingEdge_geneToPathwayTab[ is.na(LeadingEdge_geneToPathwayTab) ] = 0 

   LeadingEdge_geneToPathwayTab_long = LeadingEdge_geneToPathwayTab %>% rownames_to_column( var = "SYMBOL" ) %>% pivot_longer( cols = -SYMBOL ) %>% filter( value == 1 )

   LeadingEdgeOverlapTab = table( LeadingEdge_geneToPathwayTab_long[[1]] , LeadingEdge_geneToPathwayTab_long[[2]] )
   LeadingEdgeOverlapTab = t(LeadingEdgeOverlapTab) %*% LeadingEdgeOverlapTab
   LeadingEdgeOverlapTab = as.matrix( LeadingEdgeOverlapTab )

   if( nrow(LeadingEdgeOverlapTab) > 1 & ncol(LeadingEdgeOverlapTab) > 1 ){

   corrPlot_data = LeadingEdgeOverlapTab
   corrPlot_data[ upper.tri(corrPlot_data) ] = 0 
  
   pdf( outFilePath , width = 10 , height = 10 )
   corrplot( corrPlot_data , 
                       is.corr = FALSE , 
                       type = "lower" , 
                       method = "number" , 
                       col.lim = c( floor(min(corrPlot_data)) , ceiling(max(corrPlot_data)) ) , 
                       col = c("red","blue") , number.cex = 0.5 , tl.cex = 0.5 )

   pheatmap( LeadingEdgeOverlapTab )

   dev.off()
   }

   return( LeadingEdge_geneToPathwayTab )


}




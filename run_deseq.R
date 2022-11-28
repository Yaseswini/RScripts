## Yaseswini Neelamraju
## Performs deseq2 analysis given the following inputs -- 
## phenoTab : phenotype table ( rownames are samples and multiple columns can be your attributes 
## countTab : rownames are genes and column names are samples 
## formula : design formula to calculate diff exp ( ex : ~ condition ) 
## coef.name : if condition has two levels ( ex Tumor and Normal ), coef.name is condition_Tumor_vs_Normal . 
##	If coef.name == NULL , the script will calculate log2foldchange for the last variable in the design formula
## 

run_deseq2 <- function( phenoTab , countTab , formula , coef.name , outDIR = NULL , outLabel = NULL , volcanoPlot = TRUE , log2FC_thresh = 1 , qvalue_thresh = 0.05 )
{
	## Generating the output directory ## 
	outDIR1 = file.path( outDIR , paste0( "DifferentialExpression_" , outLabel  ) )
	if( ! dir.exists( outDIR1 ) ){
		dir.create( outDIR1 ) }

	#outDIR1 = file.path( outDIR1 , Sys.Date() )
	#if( ! dir.exists( outDIR1 ) ){ 
	#	dir.create( outDIR1 ) }
  
	if( ! all(match( rownames(phenoTab) , colnames(countTab) ) ) )
	{
		print( "Order of rownames & colnames is different\n")
		countTab = countTab[ , rownames(phenoTab) ]
	}
	dds = DESeqDataSetFromMatrix( countData = countTab , colData = phenoTab , design = formula )
	## Filter genes with normalized counts of greater than 3 
	#dds = estimateSizeFactors(dds)
	#idx = rowSums( counts(dds, normalized=TRUE) >= 3 ) >= ncol(countTab)
	#dds = dds[idx,]
	## Differential expression
	dds = DESeq(dds)
  if( is.null(coef.name) ){
    deseq2Res = results( dds , alpha = qvalue_thresh )
  }else{
    deseq2Res = results( dds, name = coef.name , alpha = qvalue_thresh ) 
  }
	
	resultsDF = data.frame(deseq2Res) %>% 
          					rownames_to_column(var="geneID") %>% 
          					mutate( geneID1 = geneID ) %>% 
          					left_join( countTab %>% rownames_to_column( var = "geneID" ) ) %>% 
          					separate( geneID1 , c("chr","gene_id","gene_name") , "!" ) %>% 
          					left_join( geneTypes ) %>% arrange( gene_name )

	resultRDS = list( "phenotypeTable" = phenoTab , "deseq2result" = resultsDF )

	deg_q0.05 = resultsDF %>% filter( as.numeric(padj) < 0.05 ) %>% arrange( gene_name )
  deg = resultsDF %>% filter( abs(log2FoldChange) > log2FC_thresh & as.numeric(padj) < qvalue_thresh  ) %>% arrange( gene_name )
  deg_simple = deg %>% dplyr::select( gene_name , log2FoldChange ) %>% arrange( desc(log2FoldChange) ) %>% arrange( gene_name )
	# Identify differentially expressed genes 
	#sigGenes =  resultsDF %>%
					#dplyr::filter( abs(log2FoldChange) > logFCThresh & padj < pvalueThresh ) %>%
					#mutate( status = if_else( log2FoldChange > 0 , "Up-Regulated" , "Down-Regulated") )

	#sigGenes.splitByStatus = split( sigGenes , sigGenes$status )
	if( ! is.null( outLabel ) )
	{
		saveRDS( resultRDS , file = file.path( outDIR1 , paste( basename(outLabel) , "_DESeq2_Result" , ".rds" , sep = "" ) ) )
		saveRDS( dds , file = file.path( outDIR1 , paste0( basename(outLabel) , "_DESeq2_dds_object.rds" ) ) )
    write.table( resultsDF , file = file.path( outDIR1 , paste( basename(outLabel) , "_DESeq2_Result"  , ".txt" , sep = "" ) ), quote = F , sep = "\t" , row.names = F , col.names = T )
		write.table( phenoTab , file = file.path( outDIR1 , paste0( basename(outLabel) , "_DESeq2_PhenoTable" , ".txt" ) ) , quote = F , sep = "\t" , row.names = T , col.names = NA )
		write.table( deg_q0.05 , file = file.path( outDIR1 , paste0( basename(outLabel) , "_DESeq2_DEG_q0.05" , ".txt" ) ) , quote = F , sep = "\t" , row.names = T , col.names = NA )
		write.table( deg , file = file.path( outDIR1 , paste0( basename(outLabel) , "_DESeq2_DEG_log2FC_" , log2FC_thresh , "_qvalue_" , qvalue_thresh , ".txt" ) ) , quote = F , sep = "\t" , row.names = T , col.names = NA )
    write.table( deg_simple , file = file.path( outDIR1 , paste0( basename(outLabel) , "_DESeq2_DEG_log2FC_" , log2FC_thresh , "_qvalue_" , qvalue_thresh , "_simple" ,".txt" ) ) , quote = F , sep = "\t" , row.names = T , col.names = NA )
    #saveRDS( sigGenes.splitByStatus , file = paste( outLabel , "_DESeq2_DEG_logFC_" , logFCThresh , "_pvalue_" , pvalueThresh , "_" , Sys.Date() , ".rds" , sep = "" ) )
		
		writeLines(capture.output(sessionInfo()) , file.path( outDIR1 , paste( basename(outLabel) , "_DESeq2_sessionInfo.txt" , sep = "" ) ) )
	
		# plot  mean of normalized counts on x-axis and log2foldchange on the y-axis
		pdf( file.path( outDIR1 , paste( basename(outLabel) , "_DESeq2_Result_MAplot" , ".pdf" , sep = "" ) ) )
		# p = ggplot( resultsDF , aes( x = baseMean , y = log2FoldChange ) ) + geom_point() 
		# p = p + theme_bw() 
		# print(p)
		DESeq2::plotMA( deseq2Res )
    dev.off()
    # ## plot examples of up and downregulated to ensure foldchange interprettation
    # up1 = deg_q0.05 %>% filter( log2FoldChange > 0 ) %>% arrange( desc(log2FoldChange) ) %>% head(n=1) 
    # down1 = deg_q0.05 %>% filter( log2FoldChange < 0 ) %>% arrange( log2FoldChange ) %>% head(n=1) 

    # up1_plots =  plotCounts( dds , gene = up1 %>% .$geneID  , intgroup = "PrognosticGroup" , returnData = TRUE , normalized = TRUE , transform = FALSE ) %>% 
    #                     rownames_to_column( var = "RNAseq.Identifier" ) %>% 
    #                     dplyr::rename( "Normalizedcounts" = "count" ) %>% 
    #                     mutate( "Log10_NormalizedCounts" = log10( Normalizedcounts + 1 ) ) %>% 
    #                     left_join( phenoTab  %>% rownames_to_column( var = "RNAseq.Identifier" ) ) %>% 
    #                     ggplot( aes( x = PrognosticGroup  , y =  Log10_NormalizedCounts ) ) + geom_boxplot() +
    #                     labs( title = paste0( "log2FC = " , up1 %>% .$log2FoldChange ) ) +
    #                     theme_bw()
    # print( up1_plots )
    # down1_plots = plotCounts( dds , gene = down1 %>% .$geneID , intgroup = "PrognosticGroup" , returnData = TRUE , normalized = TRUE , transform = FALSE ) %>% 
    #                     rownames_to_column( var = "RNAseq.Identifier" ) %>% 
    #                     dplyr::rename( "Normalizedcounts" = "count" ) %>% 
    #                     mutate( "Log10_NormalizedCounts" = log10( Normalizedcounts + 1 ) ) %>% 
    #                     left_join( phenoTab  %>% rownames_to_column( var = "RNAseq.Identifier" ) ) %>% 
    #                     ggplot( aes( x = PrognosticGroup , y =  Log10_NormalizedCounts ) ) + geom_boxplot() +
    #                     labs( title = paste0( "log2FC = " , down1 %>% .$log2FoldChange ) ) +
    #                     theme_bw()
    # print( down1_plots )

    if( volcanoPlot ){

      plotDat = resultsDF %>% 
                      dplyr::filter( ! is.na(log2FoldChange) ) %>% 
                       mutate( group = case_when( log2FoldChange > log2FC_thresh & padj < qvalue_thresh ~ "UpRegulated" , 
                                         log2FoldChange < -log2FC_thresh & padj < qvalue_thresh ~ "DownRegulated" , 
                                         padj > qvalue_thresh ~ "NonDEG") , 
                               negLog10padj = -log10(padj) ) 

      deg_count = plotDat %>% 
                filter( group %in% c("UpRegulated","DownRegulated") ) %>% 
                dplyr::select( group , geneID ) %>% 
                group_by( group ) %>% 
                tally() %>% 
                mutate( label = paste0( gsub("Regulated","",group) , "(n=" , n ,")") )

      TopUP = plotDat %>% filter( group == "UpRegulated" ) %>% arrange( desc(log2FoldChange) ) %>% head( n = 20 ) %>% .$gene_id
      TopDOWN = plotDat %>% filter( group == "DownRegulated" ) %>% arrange( log2FoldChange ) %>% head( n = 20 ) %>% .$gene_id

      plotDat = plotDat %>% mutate( pointLabel = if_else( is.element( gene_id , c(TopUP,TopDOWN) ) , gene_name , "" ) )
      

      ## Plot params 
      hline_intercept = -log10(qvalue_thresh)
      vline_intercepts = c( -log2FC_thresh , log2FC_thresh ) 
      maxFC_in_dat = ceiling( max( abs(plotDat$log2FoldChange) ) )
      maxPvalue = ceiling( max( abs(plotDat$negLog10padj) , na.rm = T ) )
      downRegulated_text_coord = -log2FC_thresh - 2
      upRegulated_text_coord = log2FC_thresh + 2 

      volcanoPlot = ggplot( plotDat , aes( x = log2FoldChange , y = negLog10padj  ) ) + 
                      geom_point( aes( fill = factor(group) , color=factor(group) , size = abs(log2FoldChange) ) , alpha = 0.5  )
      volcanoPlot = volcanoPlot + scale_color_manual( values = c("UpRegulated"="darkmagenta","DownRegulated"="darkgreen","NonDEG"="grey") ) 
      volcanoPlot = volcanoPlot + scale_size_continuous( range = c(0.2,3) )
      volcanoPlot = volcanoPlot + theme_bw() + xlim( -maxFC_in_dat , maxFC_in_dat )
      volcanoPlot = volcanoPlot + theme( axis.text = element_text( size = 14 , face = "bold" ) , 
                                         axis.title = element_text( size = 14 , face = "bold" ) ,
                                         panel.grid.major = element_blank() ,   
                                         panel.grid.minor = element_blank() , 
                                         legend.position = "bottom" )
      volcanoPlot = volcanoPlot + geom_hline( yintercept = 1.30103 , linetype = "dashed" , color = "darkgrey" )
      volcanoPlot = volcanoPlot + geom_vline( xintercept = vline_intercepts , linetype = "dashed" , color = "darkgrey" )
      volcanoPlot = volcanoPlot + labs( y = "-log10( Adjusted P-value )" ) 
      #volcanoPlot = volcanoPlot + geom_text_repel( aes(label = pointLabel) , size = 2  , max.overlaps = Inf , force=1, point.padding=unit(1,'lines'),
      #                               direction='y', nudge_x=0.1,segment.size=0.2)
      volcanoPlot = volcanoPlot + geom_text_repel( aes(label = pointLabel) , nudge_x = .15,
          box.padding = 0.5,
          nudge_y = 1,
          segment.curvature = -0.1,
          segment.ncp = 3,  , color = "red"  , segment.color = "black" , 
          segment.angle = 20 , size = 4 , max.overlaps = Inf )

      volcanoPlot = volcanoPlot + 
                        annotate("text", x = downRegulated_text_coord , y = maxPvalue , 
                                 label = deg_count %>% filter( group == "DownRegulated" ) %>% .$label , 
                                 label.padding=unit(4, "lines") ) + 
                        annotate("text", x = upRegulated_text_coord , y = maxPvalue , 
                                 label =  deg_count %>% filter( group == "UpRegulated" ) %>% .$label  )
      pdf( file.path( outDIR1 , paste( basename(outLabel) , "_DESeq2_Result_VolcanoPlot" , ".pdf" , sep = "" ) ) )
      print( volcanoPlot )
      dev.off()
	}

	return( resultRDS )
 }
}

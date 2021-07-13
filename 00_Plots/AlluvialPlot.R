# / Script params 
# Author : Yaseswini Neelamraju 
# Script name : AlluvialPlot.R
# Description : The script is used to generate alluvial plots
# / 
#
# input data :
# The input data should like this : 
# gene        deg_comaprison1 deg_comparison2 
# g1            value1     value2     value3   ... valuen
# g2 
# g3

library( ggalluvial )

inData = read_tsv( <Enter the input file name> )

colnames_as_alluvia = colnames(inData)[ colnames(inData) != "gene" ]

## The lines below will achieve these things :  
## Add up/down regulated group annotations 
## s
data_alluvial = inData %>% 
					pivot_longer( cols = colnames_as_alluvia , names_to = "group" , values_to = "log2FC" ) %>% 
					mutate( log2FC_group = case_when( log2FC > 0 ~ "UpRegulated" , 
													  log2FC < 0 ~ "DownRegulated" , 
													  is.na(log2FC) ~ "Undetected" ) ) %>% 
					pivot_wider( names_from = "group" , values_from = log2FC_group )


data_alluvial = inData %>% 
					pivot_longer( cols = colnames_as_alluvia ) %>% 
					mutate( value = if_else( is.na(value) , "Undetected" , value ) ) %>% 
					pivot_wider( names_from = name , values_from = value )
## Convert this format to long format to include group wise comparison : unite is used to achieve that 

data_alluvial1 = data_alluvial %>% 
					unite( "combo" , colnames_to_unite , remove = F  ) %>% 
					pivot_longer( cols = all_of( comparisons ) , names_to = "comparison" , values_to = "group" ) %>% 
					select( -geneID ,-gene_id ) %>% 
					group_by( comparison , group , combo ) %>%
					tally() %>% 
					mutate( group = factor( group , levels = rev(c("Undetected","Unchanged","UpRegulated","DownRegulated"))) , 
							comparison = factor( comparison , levels = comparisons ) )

# colors : 
colorPalette_stratum = data.frame( combo = unique( data_alluvial1$combo ) ) %>% 
								mutate( col = case_when( grepl("DownRegulated",combo) ~ "forestgreen" , 
														 grepl("UpRegulated",combo) ~ "darkmagenta" , 
														 combo == "Undetected_Undetected" ~ "grey" ) )
colorPalette_stratum = setNames( colorPalette_stratum$col , colorPalette_stratum$combo )
colorPalette_alluvival = c("Undetected"="grey","UpRegulated"="darkmagenta","DownRegulated"="forestgreen")

p = ggplot(  data_alluvial1 ,
           aes(x = comparison , 
               stratum = group ,
               alluvium = combo,
               y = n,
               label = group )) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_alluvium(aes(fill = combo)) +
  geom_stratum(aes( fill = group ) ) +
  scale_fill_manual( values = c(colorPalette_stratum,colorPalette_alluvival) ) + 
  geom_text( stat = "stratum", size = 4 ) +
  theme_minimal() +
  ylab(" ") +
  xlab(" ") +
  guides(fill= "none") + 
  theme(axis.text.x = element_text(face="bold", size=14))

print(p)

}
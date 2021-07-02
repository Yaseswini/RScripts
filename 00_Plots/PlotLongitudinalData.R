# ===== plot longitudinal data ===== # 
# The script is used to plot average expression/any values per gene across time points/other categorical groups 

# input data :
# The input data should like this : 
# gene        timePoint1 timePoint2 timePoint3 ... timePointn
# g1            value1     value2     value3   ... valuen
# g2 
# g3 

## Inputs from user : 
inData = read_tsv( "<path_to_your_file>" )
yaxisLabel = "choose your yaxis label"
plotTitle = "choose your plot title"
outFile = "enter the output filename"

p = inData %>% pivot_longer( cols = -contains("gene") , names_to = "xaxis_group" ) %>%
        ggplot(  aes( x = xaxis_group , y = value , color = gene ) ) + 
        geom_point() + 
        geom_line( aes(group=1) ) + 
        facet_wrap( gene ~ . ) + 
        theme_bw() + 
        theme( legend.position = "none" ) + 
        labs( x = "" , y = yaxisLabel , title = plotTitle )
ggsave( p , file = outFile )

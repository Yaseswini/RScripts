# @author : Yaseswini Neelamraju
# This function takes an XML file from http://software.broadinstitute.org/gsea/downloads.jsp , 
# extracts GO ids and associated MSigDB terms 

# Suggestions are welcome!

# args xmlFile is the path to the downloaded XML file from MSigDB
parse_MSigDB_for_GOTerms <- function( xmlFile )
{
     # Load neccessary packages 
     require( XML )
     require( dplyr )

     xmlData = xmlTreeParse( xmlFile )
     xmlData.root = xmlRoot( xmlData )

     xmlData.df = lapply( 1:xmlSize(xmlData.root) , function(i) return( c( xmlAttrs( xmlData.root[[i]] )[1] , xmlAttrs( xmlData.root[[i]] )[9] , xmlAttrs( xmlData.root[[i]] )[10] , xmlAttrs( xmlData.root[[i]] )[12] ) ) )
     xmlData.df = data.frame(do.call( "rbind" , xmlData.df ))

     # C5 is the category code for GO terms
     xmlData.df.C5 = dat1 %>% filter( CATEGORY_CODE == "C5" )
     xmlData.df.C5$EXTERNAL_DETAILS_URL = basename(xmlData.df.C5$EXTERNAL_DETAILS_URL)
                         
  return( xmlData.df.C5 )
}
#Danielle De La Pascua
#saving ugly gbif files
#6-17-21

#pull up gbif file
stgl <- read.delim("glandulosus_herb_spec.csv")
#save it again
write.csv(stgl, file = "glandulosus_herb_spec.csv")



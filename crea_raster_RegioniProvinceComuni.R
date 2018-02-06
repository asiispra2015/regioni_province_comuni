#6 Febbraio 2018
#
#Problema: l'informazione su regioni, province e comuni serve in diversi contesti dell'analisi.
#Ad esempio: i raster degli eventi di dust sono creati associando l'evento alle aree "nord","centro","sud","sicilia" e "sardegna"
#aree individuate dai nomi delle regioni. Il file GADM (shapefile) che contiene la suddivisione amministrativa dell'Italia
#non copre interamente la griglia (shapefile/raster) utilizzata per l'Italia. Lo shapefile GADM infatti segue la linea di 
#costa mentre la netfish e la sua versione raster per ogni centroide hanno un quadrato associato.
#
#Da questa differenza ne deriva che se rasterizziamo lo shapefile GADM ci perdiamo molti punti sulla costa. Ad esempio:
#i raster generati per l'evento dust non coprono tutta la griglia/netfish sulla linea di costa.

#Quindi necessario recuperare queste informazioni che vanno perse passando dallo shapefile GADM alla sua versione raster.

#Soluzione possibile: fare un buffer (1km) attorno al GADM (ad esempio a livello 1, regioni) e quindi rasterizzare questo
#shapefile più ampio in modo da associare a ciascuna cella della netfish/griglia l'info sulla regione di appartenenza.

#Problema: se faccop un buffer su tutta l'Italia nel suo insieme i confini regionali si vanno a sovrapporre. A me interessa
#invece solo dilatare le regioni verso la costa.

#Soluzione: prendo una singola regione, ne faccio un buffer a 1km. A questo punto la singola regione copre bene
#la mia netfish/griglia. Quindi posso rasterizzare la regione senza il timore di perdere l'informazione lungo la costa.

#Della regione rasterizzata voglio solo la linea di costa. Semplice: maschero la regione rasterizzata con un raster che descrive la
#la linea di costa dell'Italia (raster::boundaries).

#Alla fine del ciclo in purrr::map avrò una lista dove:
# - le regioni interne senza confini con il mare o con altri state saranno per costruzione NULL (sono regioni di cui non debbo recuperare alcuna linea di costa)
# - le regioni con confine lungo il mare o con altri stati saranno dei raster fatti di NA e con valore ID_1 (nel caso delle regioni) lungo la linea di confine

#Passo finale: faccio il mosaic delle linee di costa delle regioni per poi aggiungere il tutto alla versione rasterizzata
#del GADM. QUesto nuovo raster ora sarà uguale alla versione rasterizzata del GADM ma con in più la linea di costa.
#In questo modo posso costruire dei raster che per ogni cella della griglia/netfish mi danno l'info su regiioni, province comuni
#senza perdere informazioni a causa di una differenza tra lo shapefile del GADM e la netfish (o la sua versione rasterizzata)

#ATTENZIONE: IL FILE GADM ORIGINALE DELLE PROVINCE RIPORTAVA DUE VOLTE LA PROVINCIA DI BARLETTA. IN QGIS MEDIANTE "DISSOLVE" I DUE POLIGONI RELATIVI
#ALLO STESSO ID_2 sono stati unificati.
rm(list=objects())
library("raster")
library("rgdal")
library("rgeos")
library("purrr")
library("dplyr")
options(warn = 2,error=recover)


# Alcuni parametri --------------------------------------------------------

FEATURES<-c("regioni","province","comuni")[2] #scegliere 1,2,3
stopifnot(length(FEATURES)==1)
BUFFER_SIZE<-1000 #1km, in quanto GADM va convertito in epsg 32632 in quanto gBuffer richiede dati proiettati
"../../../copernicus/griglia_tif_nc/griglia.tif"->GRIGLIA_RASTER

#GADM DSN ovvero la directory dove si trovano gli shapefile di GADM
GADM_DSN<-"."

if(FEATURES=="regioni"){

  indice<-1
  LAYER<-"Reg_2016_WGS84"
  qualeCampo<-"COD_REG"
  
}else if(FEATURES=="province"){

  indice<-2
  LAYER<-"CMprov2016_WGS84"
  qualeCampo<-"COD_PRO"  
  
}else if(FEATURES=="comuni"){
  
  indice<-3
  LAYER<-"Com2016_WGS84"
  qualeCampo<-"PRO_COM"  
  
}else{
  stop("Quale feature?")
}

qualeLayer<-paste0("ITA_adm",indice) 


#griglia.tif è la versione raster della griglia shapefile (netfish) dell'ITALIA passata da Stafoggia
tryCatch({
  raster(GRIGLIA_RASTER)
},error=function(e){
  stop("Errore lettura file raster griglia.tif")  
})->griglia

#individuo le linee di costa del raster
raster::boundaries(griglia)->lineaCosta
lineaCosta[lineaCosta!=1]<-NA

#leggiamo i dati del GADM. Se il programma non viene fatto girare nella directory degli shapefile
#modificare il dsn 
tryCatch({
  readOGR(GADM_DSN,LAYER)
},error=function(e){
  stop(sprintf("Errore lettura %s",qualeLayer))
})->italia

#il GADM ha epsg 4326
spTransform(italia,CRS("+init=epsg:32632"))->italiaUTM

# #SE Si TRATTA DELLE PROVINVE DOBBIAMO AGGIUSTARE LA PROVINCIA DI Barletta-Andria-Trani che nello shapefile appare due volte creando un problema
# #a gBuffer
# if(indice==2){
#   italiaUTM[italiaUTM[[c(qualeCampo)]]==6,]->barletta
#   italiaUTM[italiaUTM[[c(qualeCampo)]]!=6,]->tmpItalia
# 
#   gUnaryUnion(barletta)->UnitedBarletta
#   SpatialPolygonsDataFrame(UnitedBarletta,data.frame(ID_2=6))->spBarletta
#   spBarletta@data<-barletta@data[1,]
#   gUnion(tmpItalia,spBarletta,byid =TRUE)->zz
#   rm(tmpItalia)
#   rm(barletta)
# }#fine if

#importante convertire gli ID in character e quindi in integer
as.integer(as.character(italiaUTM[[c(qualeCampo)]]))->italiaUTM[[c(qualeCampo)]]

#rasterizziamo l'italia in base a campo
rasterize(italiaUTM,griglia,field=qualeCampo)->rItaliaUTM

#ciclo su id regioni/province/comuni
purrr::map(italiaUTM[[qualeCampo]],.f=function(id){

  italiaUTM[italiaUTM[[c(qualeCampo)]]==id,c(qualeCampo) ]->singleFeature

  gBuffer(singleFeature,width=BUFFER_SIZE,id=qualeCampo,byid=TRUE)->bufferedFeature
  rasterize(bufferedFeature,lineaCosta,field=qualeCampo)->featureRaster
  #prendiamo solo la linea di confine con mare/stato estero
  mask(featureRaster,mask = lineaCosta)->maskedFeature
  #se maskedFeature tutta NA significa che si tratta di una regione interna quale l'Umbria ( o provincia o comune)
  #In questo caso non ci interessa il buffer
  if(cellStats(maskedFeature,stat = function(x,na.rm){
    
    ifelse(length(which(!is.na(x))),TRUE,FALSE)
    
  },na.rm=TRUE)) return(maskedFeature)

  NULL #trattasi di feature che non si affaccia al mare o all'esterno

}) %>% compact->listaCoste

#listaCoste è una lista che contiene per le regioni (oppure per le province o per i comuni) le linee di confine con il mare
#o con gli stati esteri. Sono le linee che si perdono facendo una semplice rasterizzazione del file GADM senza prima farne il buffer

#funzione per fare il mosaic delle linee di costa
sovrapponi<-function(x,na.rm=T){ 
  
  if(all(is.na(x))) return(NA) #ad esempio per l'Umbria o per una provincia o un comune che non confina con l'esterno
  
  if(all(!is.na(x))) return(min(x)) #può capitare che nel processo di rasterizzazione post-buffer una cella vada ad appartenere a due regioni
                                    #Si tratta comunque di linee di confine. Prendiamo l'ID minimo
  
  x[!is.na(x)]   #in questo caso sto facendo il mosaic tra due raster che non hanno alcuna cella in comune..in questo caso è semplice:
                 #per unire una cella NA con una non NA prendo il valore della cella non NA
}

#costaItalia: linea di costa dell'Italia ottenuta unendo solo lre regioni/province/comuni che si affacciano all'esterno (sul mare o su stati esteri)
reduce(listaCoste,mosaic,fun=sovrapponi)->costaItalia

#ultimo step: all'Italia rasterizzata associo la sua linea di costa..il pezzo che perdo facendo una semplice rasterizzazione del GADM
mosaic(costaItalia,rItaliaUTM,fun=sovrapponi)->rasterFinale

#funzione per riempire gli NA
fill.na <- function(x, i=5) {

  if(all(is.na(x))) return(NA)

  if(!is.na(x[i])) return(x[i])

  as.integer(names(table(x)[table(x)==max(table(x))])[1])

}  

#riempiamo gli NA, passiamo tre volte refill
focal(rasterFinale,w=matrix(1,3,3),fun=fill.na,pad=TRUE,NAonly=TRUE)->refilled
focal(refilled,w=matrix(1,3,3),fun=fill.na,pad=TRUE,NAonly=TRUE)->refilled
focal(refilled,w=matrix(1,3,3),fun=fill.na,pad=TRUE,NAonly=TRUE)->refilled

#ora eliminiamo gli NA che abbiamo riempito ma che non fanno parte della griglia
mask(refilled,mask=griglia)->daScrivere
writeRaster(daScrivere,paste0("istat_raster_",FEATURES,".tif"),overwrite=TRUE)

#Aquest paquet implementa la classe ExpressionSet, està preparat per contindre 
#informació de microarrays
library(Biobase)
#Generem una matriu amb valors d'expressió de cada gen del experiment 
expressionValues <- matrix (rnorm (300), nrow=30)
#Afegim noms a les columnes, que seran els identificadors de cada mostra
colnames(expressionValues) <- paste0("sample",1:10)
head(expressionValues)

#Es important que per crear un ExpressionSet els noms de les columnes del 
#objecte que contindrà els valors d'expressio i que es guardaran en assayData,
#concideixin amb els noms de les files del objecte que contindrà les covariancies
#que es guaradara en phenoData.

#Generem la taula de covariancies on hi ha informació de cada mostra del estudi
targets <- data.frame(sampleNames = paste0("sample",1:10),
                      group=c(paste0("CTL",1:5),paste0("TR",1:5)),
                      age = rpois(10, 30), 
                      sex=as.factor(sample(c("Male", "Female"),10,replace=TRUE)),
                      row.names=1)
head(targets, n=10)

#Creem un ojbecte que conte informació dels gens, en aqeust cas contindrà
#el nom simplement
myGenes <-  paste0("gene",1:30)

#Finalment creem un objecte que contindrà informació addicional
myInfo=list(myName="Guillem", 
            myLab="Casa",
            myContact="guillemcasa13@gmail.com", 
            myTitle="Exercici exemple")
show(myInfo)

#El format amb el qual tenim la informació es correcte per la inmensa majoria 
#de aplicacions, pero tot i això hi ha algunes que hem de fer alguna cosa extra
#Per exemple per fer un PCA
pcs<-prcomp(expressionValues)
names(pcs)
barplot(pcs$sdev)
plot(pcs$rotation[,1], pcs$rotation[,2], 
     main="Representation of first two principal components")
text(pcs$rotation[,1], pcs$rotation[,2], targets$group, cex=0.8, pos=3)

#Per a poder fer servir ExpressionSet cal generar un objecte d'aquesta classe
#Començarem per la part del AssayData
#Es l'element principal de la classe, es una matriu amb tantes files com gens
#estudiats i tantes columnes com mostres.
#Usant la funció ExpressionSet indiquem quin es l'objecte AssayData
myEset <- ExpressionSet(expressionValues)
class(myEset)

#Veiem com l'objecte es de la classe Biobase
show(myEset)
#Al fer show ens mostra quina informació hem introduit de la que accepta aquest
#objecte.

#Ara ja passem a les covariancies
#Per poder passar el nostre objecte "targets" al objecte "myEset" cal fer uns 
#passos.
#La classe "AnnotatedDataFrame" conte un data frame que volem que tingui 
#informació extra de les columnes, per això primer caldrà passar-ho a aquesta
#classe.
#En primer lloc creearem un data frame amb la informació extra que voldrem
columnDesc<- data.frame(labelDescription=c("Treatment/Control",
                                           "Age at disease onset",
                                           "Sex of patient(Male/Female"))
#Creem un objecte de classe "AnnotatedDataFrame" i afegim com a dades l'objecte
#targets i com a dades extra la informació extra
myAnnotDF<-new("AnnotatedDataFrame", data=targets, varMetadata=columnDesc)
show(myAnnotDF)

#No hem afegit els noms de les mostres com a labels ja que no es una columna del
#objecte phenoData
#Ara ja ho podem afegir al ExpressionSet
phenoData(myEset)<-myAnnotDF
show(myEset)
#Veim com ara ja tenim la phenoData afegida

#Passant als feature Data
#El nombre de files en featureData ha de coincidir amb el de assayData i els 
#nomes de featureData ha de coincidir amb els noms de la matriu assayData

#Una forma alternativa es si tenim els features en un vector ho podem afegir 
#diractmanet al objecte "featureNames"

myEset <- ExpressionSet(assayData=expressionValues,
                        phenoData=myAnnotDF,
                        featureNames =myGenes)
show(myEset)

#Ara passem a afegir informació del experiment, es fa amb la classe MIAME,
#Minimum Information About a Microrarray Experiment

myDesc <- new("MIAME", name= myInfo[["myName"]],
              lab= myInfo[["myLab"]],
              contact= myInfo[["myContact"]] ,
              title=myInfo[["myTitle"]])
print(myDesc)

#De nou el podem afegir al objecte clau

myEset <- ExpressionSet(assayData=expressionValues,
                        phenoData=myAnnotDF,
                        fetureNames =myGenes,
                        experimentData = myDesc)
show(myEset)

#Ara ja hem introduit tota la informació que voliem al objecte "ExperessionSet"
#pel que ara ja podrem començar a trebllar-hi

#Per accedir als valors introduits en ExpressionSet el que fem es usar "accessors"
#que son funcions especials.
#Per exemple, per accedir a la phenoData, usarem dos accessors, phenoData per 
#accedir a ExpressionSet'sphenoData i pData per accedir a la data.

dim(exprs(myEset))
class(phenoData(myEset))
class(pData(phenoData(myEset)))
head(pData(phenoData(myEset)))
head(pData(myEset))

#L'objecte ExpressionSet es va crear per fer que la manipulació de dades fos 
#consistent amb altres objectes de R. Per exmple crear un subset d'aquest 
#objecte vol dir fer un subset de la matriu d'expressió, informació de la mostra 
# i features simulatnament si ho fem correctament.
smallEset<-myEset[1:15, c(1:3,6:8)]
dim(exprs(smallEset))
dim(pData(smallEset))
head(pData(smallEset))

#El paquet GEOquery
#NCBI GEO es un repositori per dades experimentals high-throuhput. D'alla podem
#descarregar les dades de qualsevol experiment en diferent formats.
#Nosaltres ara parlarem de GEOquery que es un paquet basat en Bioconductor.
#El paquet permet descarregar i manipular dades de GEO per poder usar-les en
#Bioconductor. Ara mostrarem com descarregar datasets.

#En primer lloc instal·lem el paquet 

if (!require(GEOquery)) {
  BiocManager::install("GEOquery")
}
require(GEOquery)
#Amb la funció getGEO descarreguem el fitxer amb ID GSE27174 si esta en format
#GSE
gse <- getGEO("GSE27174", GSEMatrix=TRUE, AnnotGPL=TRUE)
class(gse)
names(gse)
length(gse)
gse[[1]]
#Ara capturem la informació que tenim en l'objecte ExpressionSet
essetFromGeo<-gse[[1]]
#Veiem els gens anotats, les mostres i els seus valors d'expressió.
head(exprs(essetFromGeo))
#Mirem les covariancies
colnames(pData(essetFromGeo))
#Mirem les que ens han cridat l'atenció
pData(essetFromGeo)[,39:40]
#Per a descarregar un fitxer en format GSD
gds <- getGEO("GDS4155")
class(gds)
slotNames(gds)
head(Meta(gds))
#Per a convertir-ho a objecte ExpressionSet
eset<-GDS2eSet(gds,do.log=FALSE)
eset

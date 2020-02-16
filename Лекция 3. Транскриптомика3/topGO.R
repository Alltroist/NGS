library(topGO)
options(stringsAsFactors = FALSE)
go = read.csv('~/skoltech/projects/evo.devo/input/GO/Homo_sapiens.GRCh37.74.GO.csv.gz')
go[1:2,]
# Ensembl.Gene.ID GO.Term.Accession
# 1 ENSG00000261657        GO:0006810
# 2 ENSG00000261657        GO:0016021
go = split(go$GO.Term.Accession,go$Ensembl.Gene.ID)
u = setNamesnames(go)
s = setNames(factor(as.integer(u %in% sample(u,1000))),u)
str(s)

# use topGO annotation
tgo1 <- new("topGOdata", ontology = "BP",
                    allGenes = s,
                    nodeSize = 10,
                    annotationFun = annFUN.org,mapping='org.Hs.eg.db',ID='Ensembl')

r1<- runTest(tgo1, algorithm = "classic", statistic = "fisher")
hist(score(r1))
table(p.adjust(score(r1),m='BH')<0.5)
GenTable(tgo1,r1,topNodes=10)
showSigOfNodes(tgo1, score(r1), firstSigNodes = 5, useInfo ='all')

# use external annotation
tgo2 <- new("topGOdata", ontology = "BP",
            allGenes = s,
            nodeSize = 10,
            annotationFun = annFUN.gene2GO,gene2GO=go)
r2<- runTest(tgo2, algorithm = "classic", statistic = "fisher")
showSigOfNodes(tgo2, score(r2), firstSigNodes = 5, useInfo ='all')

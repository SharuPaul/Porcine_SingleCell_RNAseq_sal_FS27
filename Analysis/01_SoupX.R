library(DropletUtils)
library(ggplot2)
library(SoupX)
set.seed(5)

## Perform ambient RNA estimation and removal on each individual sample 

### Sample S1
#### Provide file path to Cell Ranger outputs:
sc = load10X('CellRangerOutputs/S1/outs') 

#### Check the meta data to make sure we also have cluster information and t-SNE coordinates:
head(sc$metaData, n = 3)

#### Plot data according to Cell Ranger cluster assignments and t-SNE coordinates:
dd = sc$metaData # create an object with all the metadata
mids = aggregate(cbind(tSNE1,tSNE2) ~ clusters,data=dd,FUN=mean) # determine t-SNE coordinates for middle of each cluster
gg = ggplot(dd,aes(tSNE1,tSNE2)) + # make a t-SNE plot
  geom_point(aes(colour=factor(clusters)),size=0.2) +
  geom_label(data=mids,aes(label=clusters)) 
plot(gg) # show plot

#### Check expression patterns for some canonical genes:
dd$CD3E = sc$toc["CD3E", ] 
dd$JCHAIN = sc$toc["JCHAIN", ] 
dd$CD79B = sc$toc["CD79B", ] 
dd$FLT3 = sc$toc["FLT3", ] 
dd$AIF1 = sc$toc["AIF1", ] 
dd$HBB = sc$toc["HBB", ] 

a1 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = CD3E > 0)) + ggtitle('CD3E') # which cells express this gene?
a2 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = JCHAIN > 0)) + ggtitle('JCHAIN')
a3 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = CD79B > 0)) + ggtitle('CD79B')
a4 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = FLT3 > 0)) + ggtitle('FLT3')
a5 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = AIF1 > 0)) + ggtitle('AIF1')
a6 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = HBB > 0)) + ggtitle('HBB')
(a1+a2+a3) / (a4+a5+a6) # show 6 plots

## if we assumed all cells were nothing but soup, which cells still show higher than
## expected expression for the gene (TRUE = expression levels higher than expected if 
## cell was just soup, so likely real expression). This just gives us an idea of soup 
## expression, this is NOT a formal analysis used for removing the soup RNA.
a1 <- plotMarkerMap(sc, "CD3E") + ggtitle('CD3E')  
a2 <- plotMarkerMap(sc, "JCHAIN") + ggtitle('JCHAIN')
a3 <- plotMarkerMap(sc, "CD79B") + ggtitle('CD79B')
a4 <- plotMarkerMap(sc, "FLT3") + ggtitle('FLT3')
a5 <- plotMarkerMap(sc, "AIF1") + ggtitle('AIF1')
a6 <- plotMarkerMap(sc, "HBB") + ggtitle('HBB')
(a1+a2+a3) / (a4+a5+a6) # show 6 plots


#### Calculate the RNA soup fraction:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

#### See which genes were most affected by our soup correction:
cntSoggy = rowSums(sc$toc > 0)
cntStrained = rowSums(out > 0) 
tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) 
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) 

#### See how soup removal affects the genes we assessed expression patterns for earlier:
a1 <- plotChangeMap(sc, out, "CD3E") + ggtitle('CD3E')  
a2 <- plotChangeMap(sc, out, "JCHAIN") + ggtitle('JCHAIN')
a3 <- plotChangeMap(sc, out, "CD79B") + ggtitle('CD79B')
a4 <- plotChangeMap(sc, out, "FLT3") + ggtitle('FLT3')
a5 <- plotChangeMap(sc, out, "AIF1") + ggtitle('AIF1')
a6 <- plotChangeMap(sc, out, "HBB") + ggtitle('HBB')
(a1+a2+a3) / (a4+a5+a6) # show 6 plots

#### Save our strained count matrix to a new location:
write10xCounts("S1strainedCounts", out, version = "3", overwrite = TRUE)
rm(dd, gg, mids, out, sc, cntSoggy, cntStrained)

#Repeat with the rest of the samples

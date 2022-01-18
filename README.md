This repository contains metadata and scripts for the following manuscript. Manuscript has not been finalised to this repo may still be messy incomplete. 

Stuart KC†, Edwards RJ†, Cheng Y, Warren WC, Burt DW, Sherwin WB, Hofmeister NR, Werner SJ, Ball GF, Bateson M, Brandley MC, Buchanan KL, Cassey P, Clayton DF, De Meyer T, Meddle SL, Rollins LA (2021). Transcript- and annotation-guided genome assembly of the European starling. BioRXIV, https://doi.org/10.1101/2021.04.07.438753 † joint first author.

# VIGNETTE: Circularize genome plots

<p>It is a truth universally acknowledged, that a manuscript in possession of a genome, must be in want of a circular plot.</p>



<ul>
<li><a href="https://academic.oup.com/jhered/article/108/2/207/2726867">Talbot et al. 2016</a></li>
<li><a href="https://www.nature.com/articles/439803a">Phillips et al. 2006</a></li>
<li><a href="https://onlinelibrary.wiley.com/doi/full/10.1111/mec.13184">Rollins et al. 2015</a></li>
<li><a href="https://www.ncbi.nlm.nih.gov/pubmed/18467157">Prentis et al. 2008</a></li>
<li><a href="https://onlinelibrary.wiley.com/doi/full/10.1111/mec.13162">Colautti and Lau 2015</a></li>
</ul>

<h2>Installing Circlize</h2>

<pre class="r"><code>install.packages("circlize")
library(circlize)</code></pre>

<h2>Loading Data</h2>

You can find examples in the files located in X. I have also included screenshots of what the files should look like.

<pre class="r"><code>
cytoband.df = read.csv("cyto.csv", colClasses = c("character", "numeric", "numeric", "character", "character"), sep = ",", na.strings='NULL')
str(cytoband.df)
</code></pre>


<h2>Plotting tracks</h2>

Below is examples of the code I have used to plot the tracks on the first circlize example.

<pre class="r"><code>
circos.initializeWithIdeogram(cytoband.df, 
                              chromosome.index = paste0("chr", c(1,"1A",2,3,4,"4A",5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,"Z")),
                              labels.cex = 2, axis.labels.cex = 1.5)
                              
                              circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim},
  bg.col = c("#088A68", "#088A85", "#086A87","#084B8A", "#04B45F", "#04B486","#04B4AE", "#0489B1", "#045FB4","#01DF74", "#01DFA5", "#01DFD7","#01A9DB", "#0174DF", "#00FF80","#00FFBF", "#00FFFF", "#00BFFF","#0080FF", "#2EFE9A", "#2EFEC8","#2EFEF7", "#2ECCFA", "#2E9AFE","#58FAAC","#58FAD0", "#58FAF4", "#58D3F7","#58ACFA","#81F7BE", "#81F7D8", "#81F7F3", "#81DAF5"),
  track.height = 0.05, bg.border = NA)
  
  circos.genomicTrackPlotRegion(gffmrna.df, ylim = c(0, 250),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, ..., area = FALSE, col = c("#A30606", "#060BA3"))
                              }, track.height = 0.15)

circos.genomicTrackPlotRegion(gffgene.df, ylim = c(0, 80),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, ..., area = TRUE, col = "#E8DAEF")
                              }, track.height = 0.10)

circos.genomicTrackPlotRegion(vd.df, ylim = c(0, 16900),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, ..., area = TRUE, col = "#F2D7D5")
                              }, track.height = 0.10)

circos.genomicTrackPlotRegion(gc.df, ylim = c(0.30, 0.52),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, ..., type ='s', area = TRUE, col = "#FAE5D3")
                              }, track.height = 0.10)


circos.info()
circos.clear()
</code></pre>

![ScreenShot](/Sv3_vignette/Starling_genome.png)



<h2>Another Plot</h2>

<pre class="r"><code>

circos.initializeWithIdeogram(cytoband.df, 
                              chromosome.index = paste0("chr", c(1,"1A",2,3,4,"4A",5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,"Z",29)),
                              labels.cex = 2, axis.labels.cex = 1.5)
                              
col_fun = colorRamp2(c(0, 0.00003, 0.00008), c("blue", "green", "red"))
circos.genomicHeatmap(XtX.df.vals, col = col_fun, side = "inside", border = NA, heatmap_height = 0.1, connection_heigh=0.001, line_col = "white")

circos.genomicTrack(bed_list.A8, ylim = c(15, 30), track.height = 0.3, panel.fun = function(region, value, ...) {
  cex = (value[[19]])
  i = getI(...)
  circos.genomicPoints(region, value, cex = cex, pch = 20, col = i, ...)
})

</code></pre>

and as always, once you plot is complete, use the below:

<pre class="r"><code>
circos.clear()
</code></pre>


![ScreenShot](/Sv3_vignette/Starling_morphology.png)

<h2>Experiment with style</h2>

I would highly recomment that once you are feeling comfortable with how to plot tracks, that you look into the graphical options in Circlize. You can make very sleek and beautiful plots, as with this example below made by Stephanie Chen.

Link paper

![ScreenShot](/Sv3_vignette/waratah_genome.png)



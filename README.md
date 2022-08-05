# Sv3_StarlingGenome

![ScreenShot](/Sv3_vignette/Starling_genome.png)

This repository contains metadata and scripts for the following manuscript. 

Stuart KC†, Edwards RJ†, Cheng Y, Warren WC, Burt DW, Sherwin WB, Hofmeister NR, Werner SJ, Ball GF, Bateson M, Brandley MC, Buchanan KL, Cassey P, Clayton DF, De Meyer T, Meddle SL, Rollins LA (2021). Transcript- and annotation-guided genome assembly of the European starling. *Molecular Ecology Resources*. [https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13679](doi/full/10.1111/1755-0998.13679)


<h2><i>Circular genome plots</i></h2>

<p>It is a truth universally acknowledged, that a manuscript in possession of a genome, must be in want of a circular plot.</p>

Below contains code for some circular genome plots using the R package [https://jokergoo.github.io/circlize_book/book/](circlize). The cookbook for this is very complete and very informative. I provide some example plotting code below for those who may be struggling to get their initial plot up and running. If I use this program in the future to make more and diverse circular plots I'll add them below to provide further examples of the types of tracks than can be created.

## Plotting the circular plot

Install and load Circlize.

<pre class="r"><code>install.packages("circlize")
library(circlize)</code></pre>

## Example 1

Loading Data. 

<pre class="r"><code>cytoband.df = read.csv("cyto.csv", colClasses = c("character", "numeric", "numeric", "character", "character"), sep = ",", na.strings='NULL')
str(cytoband.df)
</code></pre>

You can find examples in the files located in X. I have also included screenshots of what the files should look like.

<pre class="r"><code>

</code></pre>


This is the first track, that sets the lengths of your chromosomes/scaffolds along the 360 degrees of the circular plot
<pre class="r"><code>
circos.initializeWithIdeogram(cytoband.df, 
                              chromosome.index = paste0("chr", c(1,"1A",2,3,4,"4A",5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,"Z")),
                              labels.cex = 2, axis.labels.cex = 1.5)
</code></pre>

Next, I like to lay down some simple colours to make my plot look a bit more visual. You can skip this track, but might be useful if you seek to specifically investigate specific chromosomes - you can chose to colour these differently and link the colour scheme throughoout your other manuscript figures.


<pre class="r"><code>circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim},
  bg.col = c("#088A68", "#088A85", "#086A87","#084B8A", "#04B45F", "#04B486","#04B4AE", "#0489B1", "#045FB4","#01DF74", "#01DFA5", "#01DFD7","#01A9DB", "#0174DF", "#00FF80","#00FFBF", "#00FFFF", "#00BFFF","#0080FF", "#2EFE9A", "#2EFEC8","#2EFEF7", "#2ECCFA", "#2E9AFE","#58FAAC","#58FAD0", "#58FAF4", "#58D3F7","#58ACFA","#81F7BE", "#81F7D8", "#81F7F3", "#81DAF5"),
  track.height = 0.05, bg.border = NA)
</code></pre>

No we start loading each track individually. This first track sets "area = FALSE" to produce two lines, using the two colours manually specified.

<pre class="r"><code>circos.genomicTrackPlotRegion(gffmrna.df, ylim = c(0, 250),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, ..., area = FALSE, col = c("#A30606", "#060BA3"))
                              }, track.height = 0.15)
</code></pre>

This next track uses area = TRUE to colour the area under the line.

<pre class="r"><code>circos.genomicTrackPlotRegion(gffgene.df, ylim = c(0, 80),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, ..., area = TRUE, col = "#E8DAEF")
                              }, track.height = 0.10)
</code></pre>

<pre class="r"><code>circos.genomicTrackPlotRegion(vd.df, ylim = c(0, 16900),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, ..., area = TRUE, col = "#F2D7D5")
                              }, track.height = 0.10)
</code></pre>

This final track has type ='s', which causes the plotted line to be less smoothed.

<pre class="r"><code>circos.genomicTrackPlotRegion(gc.df, ylim = c(0.30, 0.52),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, ..., type ='s', area = TRUE, col = "#FAE5D3")
                              }, track.height = 0.10)
</code></pre>

If at any stage you make a mistake when plotting your tracks, simply wipe the plot and start from the top.

<pre class="r"><code>circos.clear()
</code></pre>

And done!
![ScreenShot](/Sv3_vignette/Starling_genome.png)

Extra visuals like legends are best added in an image editor in my opinion, as this gives you better control of the style.

## Example 2

Same genome, different project.

<pre class="r"><code>circos.initializeWithIdeogram(cytoband.df, 
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

![ScreenShot](/Sv3_vignette/Starling_morphology.png)

(gene labels were added manually after during final figure editing).

## Experiment with style

I would highly recommend that once you are feeling comfortable with how to plot tracks, that you look into the graphical options in Circlize. You can make very sleek and beautiful plots, for example this very visually beautiful example published [https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13574](HERE).

![ScreenShot](/Sv3_vignette/waratah_genome.png)



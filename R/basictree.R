
library(dplyr)
library(phyloseq)
library(phylogeo)
library(leaflet)

#' interactivetree
#'
#' interactive phlogenetic tree that returns a phylogenetic tree as an htmlwidget
#'
#' @import leaflet
#' @import dplyr
#' @import phyloseq
#'
#' @param  physeq a phyloseg object
#' @param ladderize boolean. ladderize or not
#' @param xscale scale to correct APE's x/y discrepancy
#' @param yscale scale to correct APE's x/y discrepancy
#' @param linecolor color lines
#' @param lineweight
#' @param lineopacity
#' @param method
interactivetree <- function(physeq,
                            ladderize=TRUE,
                            xscale=NULL,
                            yscale=NULL,
                            #lineoptions
                            linecolor= "black",
                            lineweight = 2,
                            lineopacity = 1,
                            method="tree"){

    #get tree layout data
    treeSegs <- phyloseq::tree_layout(physeq, ladderize = ladderize)

    # set x and y scale values. Scale xmax and ymax to 1
    if (is.null(xscale)) {
        x_max <- max(treeSegs$edgeDT$xright)
        x_min <- min(treeSegs$edgeDT$xleft)
        y_max <- max(treeSegs$edgeDT$y)
        y_min <- min(treeSegs$edgeDT$y)
        xscale = 1/x_max
        yscale = 1/y_max

    }

    #' helper function to create a leaflet-compatable two-column matrix of lines
    #' from Ape's tree layout
    make_edgematrix <-function(){
        #get the edges
        edges <- treeSegs$edgeDT
        edges <- edges %>% add_rownames()
        start <- edges %>% mutate(lng = xleft*xscale, lat = y*yscale)
        end   <- edges %>% mutate(lng = xright*xscale, lat = y*yscale)
        space <- edges
        space$lng = NA
        space$lat = NA
        edges2 <- rbind(start, end, space) %>%
            arrange(rowname) %>%
            select(lng,lat) %>%
            as.matrix()
    }
    #' helper function to create a leaflet-compatable two-column matrix of lines
    #' from Ape's tree layout
    make_vertexmatrix <- function(){
        #and the vertices
        verts <- treeSegs$vertDT
        verts <- verts %>% add_rownames()
        start <- verts %>% mutate(lng=x*xscale, lat=vmin*yscale)
        end   <- verts %>% mutate(lng=x*xscale, lat=vmax*yscale)
        space <- verts
        space$lng = NA
        space$lat = NA
        verts2 <- rbind(start, end, space) %>%
            arrange(rowname) %>%
            select(lng,lat) %>%
            as.matrix()
    }

    #get edges, vertices, and points
    edges <- make_edgematrix()
    vertices <- make_vertexmatrix()
    points <- treeSegs$edgeDT %>%
        data.frame() %>%
        filter(!is.na(OTU)) %>%
        mutate(lng = xright*xscale, lat = y*yscale) %>%
        select(OTU,lng,lat)


    #set the basemap and add the tree
    m <- leaflet() %>%
        addPolylines(data=edges,    color=linecolor, weight=lineweight, opacity=lineopacity) %>%
        addPolylines(data=vertices, color=linecolor, weight=lineweight, opacity=lineopacity)

    #add the points at the end of the tree
    if(method=="tree"){
        taxdata <- tax_table(physeq) %>%
            data.frame() %>%
            mutate(OTU = rownames(.))
        #print(taxdata)
        pointdata <- merge(points, taxdata, by="OTU")
        #add points to the map
        m <- m %>% addCircles(data=pointdata, lng=~lng,lat=~lat, popup=~OTU)
    }


    return(m)
}

data("epoxamicin_KS")
interactivetree(epoxamicin_KS)
#p <- plot_tree(esophagus)
#p$data

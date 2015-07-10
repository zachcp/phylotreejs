library(ape)
library(phyloseq)
library(leaflet)

data(esophagus)


interactivetree <- function(physeq, xscale=100, yscale=1){

    #layout using ape/phyloseq
    treeSegs <- phyloseq::tree_layout(esophagus, ladderize = T)

    #get the edges
    edges <- treeSegs$edgeDT
    edges <- edges %>% add_rownames()
    start <- edges %>% mutate(lng = y, lat = xleft)
    end   <- edges %>% mutate(lng = y, lat = xright)
    #scale
    start <- start %>% mutate(lat = lat*xscale)
    end   <- end   %>% mutate(lat = lat*xscale)
    space <- edges
    space$lng = NA
    space$lat = NA
    edges2 <- rbind(start, end, space) %>%
        arrange(rowname) %>%
        select(lat,lng) %>%
        as.matrix()

    #and the vertices
    verts <- treeSegs$vertDT
    verts <- verts %>% add_rownames()
    start <- verts %>% mutate(lng = vmin, lat = x)
    end   <- verts %>% mutate(lng = vmax, lat = x)
    #scale
    start <- start %>% mutate(lng = lng*yscale, lat=lat*xscale)
    end   <- end   %>% mutate(lng = lng*yscale, lat=lat*xscale)
    space <- verts
    space$lng = NA
    space$lat = NA
    verts2 <- rbind(start, end, space) %>%
        arrange(rowname) %>%
        select(lat,lng) %>%
        as.matrix()

    #put them into  a leaflet map
    m <- leaflet() %>%
        addPolylines(data=edges2) %>%
        addPolylines(data=verts2)

    return(m)
}

interactivetree(esophagus)
#p <- plot_tree(esophagus)
#p$data

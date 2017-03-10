## netcontext.R - Using ggnetwork library to place gene models onto their network context
##
## by Artem Sokolov

library( optparse )

## Define accepted command-line arguments
opts <- list(
    make_option(
        c( "-n", "--network" ), type="character", default=NULL,
        help="A .sif file defining the interaction graph. Can be gzipped." ),
    make_option(
        c( "-w", "--weights" ), type="character", default=NULL,
        help="A tab-delimited two-column file that maps genes to their weights. Can be gzipped." ),
    make_option(
        c( "-s", "--seed" ), type="character", default=NULL,
        help="(Optional) A file that contains the seed genes, one per line. If not provided, the seed genes are taken to be the first column of the -w argument." ),
    make_option(
        c( "-o", "--output" ), type="character", default="out.pdf",
        help="Output filename. Default: out.pdf" ) )

## Parse the provided arguments and verify the presence of mandatory ones
op <- OptionParser( option_list = opts )
argv <- parse_args(op)
if( is.null( argv$network ) || is.null( argv$weights ) )
{
    print_help(op)
    stop( "Network and gene weights must be provided." )
}

## Load the remaining libraries
suppressMessages( library( dplyr ) )
suppressMessages( library( igraph ) )
suppressMessages( library( sna ) )
suppressMessages( library( network ) )
library( ggplot2 )
library( ggrepel )
library( ggnetwork )

## "Projects" a set of selected genes (seed) onto a gene-gene interaction network:
##   w - a named vector of gene weights (can be differential expression scores, model weights, etc.)
##   g - an igraph object that encapsulates the interaction network
##   seed - a character vector specifying genes to use as a plot seed
## Returns a ggplot object
network.context <- function( w, g, seed )
{
    ## Generate a connected subgraph over the selected genes
    sp <- shortest_paths( g, seed, seed ) %>% unlist() %>% unique()
    gs <- induced_subgraph( g, sp )

    ## Set up the associated weight vector and populate it with known values
    ww <- setNames( rep(0, length(V(gs))), names(V(gs)) )
    j <- intersect( names(ww), names(w) )
    if( length(j) > 0 )
        ww[j] <- w[j]

    ## Associate the new weight vector with the graph
    V(gs)$Weight <- ww

    ## Distinguish seed and linker nodes
    V(gs)$Type <- c( "Linker", "Seed" )[ names(V(gs)) %in% seed + 1 ]

    ## Plot all the things
    ggplot( gs, aes(x = x, y = y, xend = xend, yend = yend ) ) +
        geom_edges( size=1, aes(linetype=Interaction), color="gray" ) +
        geom_nodes( aes( fill = Weight, size=Type ), color="black", shape=21 ) +
        geom_nodetext_repel( aes( label = vertex.names ), size=4, fontface="bold" ) +
        scale_fill_gradient2( low = "steelblue", high = "tomato", mid="white", midpoint = 0 ) +
        scale_size_manual( values=c( "Linker" = 5, "Seed" = 10 ) ) +
        theme_blank()
}

main <- function( argv )
{
    ## Load the provided network
    cat( "Loading the provided .sif\n" )
    X <- read.delim( argv$network, header=FALSE, as.is=TRUE )

    ## Generate a graph from the list of edges
    cat( "Generating a network\n" )
    g <- X %>% select( V1, V3 ) %>% as.matrix() %>% graph_from_edgelist( directed = FALSE )
    E(g)$Interaction <- X$V2

    ## Handle duplicate edges by defining a new "composite-edge" interaction class
    g <- simplify( g, edge.attr.comb=function(v) {ifelse( length(v) > 1, "composite-edge", v )} )

    ## Isolate the largest connected component of the graph
    if( is_connected( g ) == FALSE )
    {
        cat( "Reducing the graph to the largest connected component\n" )
        lg <- decompose( g )
        g <- lg[[ lapply( lg, vcount ) %>% which.max() ]]
    }
    
    ## Load the weights
    w <- read.delim( argv$weights, header=FALSE ) %>%
        tibble::column_to_rownames( "V1" ) %>% as.matrix() %>% drop()

    ## Load the seed
    if( is.null( argv$seed ) )
        s <- names(w)
    else
        s <- scan( argv$seed, what=character() )

    ## Reduce the seed set to the available vertex set
    s <- V(g) %>% names() %>% intersect( s )
    if( length(s) < 1 )
        stop( "None of the requested genes are present in the graph." )
    if( length(s) > 50 )
        stop( 'The tool is not meant to plot more than 50 nodes, as this creates \"hairball\" networks that are difficult to interpret.' )

    ## Generate and save the plot
    gg <- network.context( w, g, s )
    ggsave( argv$output, gg )
}

main( argv )

#' @noRd
get_skeleton <- function(g) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop(
            "Package \"igraph\" must be installed to use this function.",
            call. = FALSE
        )
    }
    get_deg2 <- function(x) {
        dd <- igraph::degree(x)
        trim <- igraph::V(x)[names(dd[dd==2])]    
    }
    ng <- g
    trim <- get_deg2(ng)
    while(length(trim)) {
        tv <- trim[1]
        touch <- igraph::adjacent_vertices(ng, tv)[[1]]
        ng <- igraph::delete_edges(ng, igraph::E(ng)[tv %--% touch])
        ng <- igraph::add_edges(ng, touch$name)
        ng <- igraph::delete_vertices(ng, igraph::V(ng)[tv])
        trim <- get_deg2(ng)
    }
    ng
}

#' @noRd
assign_branch_to_nodes <- function(skeleton, principal_igraph) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop(
            "Package \"igraph\" must be installed to use this function.",
            call. = FALSE
        )
    }
    # a branch is defined as the shortest path between any two nodes of degree!=2 (i.e. a terminus or a branch point)
    
    branches <- igraph::as_edgelist(skeleton, names = TRUE)
    principal_nodes <- data.frame(vertex=igraph::V(principal_igraph)$name, branch=NA)
    colnames(principal_nodes) <- c("vertex","assigned_branch")
    num_branches <- nrow(branches)
    for (i in 1:num_branches) {
        branch_endpoints = branches[i,]
        # get the two skeletonnodes (P1 and P2) defining that edge (e.g. a branch point and a terminus)
        # on the mst draw shortest path between P1 and P2 (Segment1_2): shortest_path
        p1=branch_endpoints[1]
        p2=branch_endpoints[2]
        branch = igraph::shortest_paths(principal_igraph,
                                from = p1,
                                to = p2,
                                mode = 'all',
                                output ="both",
                                weights=NA)
        branch_nodes <- names(branch$vpath[[1]])
        branch_name <- paste(p1, p2, sep=".")
        indices <- match( branch_nodes, principal_nodes$vertex) 
        principal_nodes[indices,]$assigned_branch <- branch_name
    }
    return(principal_nodes)
}

#' @title Get cells between nodes.
#' @description cells that are on the shortest path on the minimal spanning tree between nodes.
#' @param endpoint_nodes Vector of length 2 denoting nodes to search between.
#' @param cds A `monocle3` cell_data_set object with `monocle3::order_cells()` run on it.
#' @return a vector of cell names.

#' @export
cells_between_vertices <- function(endpoint_nodes, cds) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop(
            "Package \"igraph\" must be installed to use this function.",
            call. = FALSE
        )
    }
    # a path is defined as the shortest path between any two nodes.
    cell_to_vertex_df <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex %>% as.data.frame()
    cell_to_vertex_df$vertex_name <- paste0("Y_",cell_to_vertex_df[,1])
    principal_igraph <- cds@principal_graph$UMAP
    skeleton <- get_skeleton(principal_igraph)
    
    principal_nodes <- data.frame(vertex=igraph::V(principal_igraph)$name, branch=NA)
    colnames(principal_nodes) <- c("vertex","assigned_branch")

    p1=endpoint_nodes[1]
    p2=endpoint_nodes[2]
    path = igraph::shortest_paths(principal_igraph,
                                    from = p1,
                                    to = p2,
                                    mode = 'all',
                                    output ="both",
                                  weights=NA)
    path_nodes <- names(path$vpath[[1]])
    cell_to_vertex_df %>%
        filter(vertex_name %in% path_nodes) %>%
        rownames() -> cell_names
    return(cell_names)
}

#' @noRd
assign_nodes_to_single_path <- function(branch_endpoints, skeleton, principal_igraph) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop(
            "Package \"igraph\" must be installed to use this function.",
            call. = FALSE
        )
    }
    # a branch is defined as the shortest path between any two nodes of degree!=2 (i.e. a terminus or a branch point)
    # get the two skeletonnodes (P1 and P2) defining that edge (e.g. a branch point and a terminus)
    # on the mst draw shortest path between P1 and P2 (Segment1_2): shortest_path
    p1 <- branch_endpoints[1]
    p2 <- branch_endpoints[2]
    branch  <- igraph::shortest_paths(principal_igraph,
                            from = p1,
                            to = p2,
                            mode = 'all',
                            output ="both")
    branch_nodes <- names(branch$vpath[[1]])
    branch_nodes <- sub('..', '', branch_nodes)     #strip first two characters"Y_" from the vector
    return(branch_nodes)
}

#' @noRd
get_branch_membership <- function(cds, reduction_method="UMAP") {
    principal_igraph <- cds@principal_graph[[reduction_method]]
    skeleton <- get_skeleton(principal_igraph) #create skeleton graph
    # each edge on skeleton is a branch.
    
    principal_nodes_assigned_branch <- assign_branch_to_nodes(skeleton, principal_igraph) 
    
    cell_assigned_vertices <- as.data.frame(cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex)
    cell_assigned_vertices <- as.data.frame(apply(cell_assigned_vertices, FUN=function(x) { paste("Y_",x,sep="")}, MARGIN=1))
    cell_assigned_vertices$cell <- rownames(cell_assigned_vertices)
    colnames(cell_assigned_vertices) <- c("vertex", "cell")
    
    #assign branch label to cells with a left join
    cells_assigned_branch = dplyr::left_join(cell_assigned_vertices, principal_nodes_assigned_branch)
    
    result <- as.character(cells_assigned_branch$assigned_branch)
    names(result) <- cells_assigned_branch$cell
    
    #update cds 
    cds@principal_graph_aux[[reduction_method]]$cell_branch <- result
    cds@colData$branch <- result
    cds@colData$pseudotime<-cds@principal_graph_aux$UMAP$pseudotime
    cds
}

#' @noRd
get_branch_genes_covar<- function(branch_expr, pseudotime, quantile_cutoff = 0.90, slot='data', binary=FALSE){
    covar <- cov(pseudotime,t(branch_expr))
    covar[is.na(covar)] <- 0
    covar <- abs(covar)
    if (binary==TRUE) {
        print("Obtaining genes covarying along the selected branch")
        result<-data.frame(as.logical(covar>quantile(covar,quantile_cutoff)), as.numeric(covar), row.names=rownames(branch_expr))
        colnames(result)<-c("Covar_classification", "Covariances")
    } else {
        result <- data.frame(as.numeric(covar), row.names=rownames(branch_expr))
        colnames(result) <- c("Covariances")
    }
    return(result)
}

#' @noRd
get_branch_correlated_genes_regress <- function(branch_expr, branch_cells){
    branch_pseudotime <- branch_cells$pseudotime
    
    corr <- cor(branch_pseudotime,t(branch_expr))
    corr[is.na(corr)] <- 0
    
    corr <- abs(corr)
    branch_corr_genes_regress <- as.data.frame(corr)
}

#' @noRd
filter_lr_genes <- function(cells, lr_network, filter) {
    if (filter=="l") {
        ligands <- lr_network$from %>% unique()
        common_genes <- intersect(ligands, row.names(cells))
        result <- cells[ligands,]
    } else if (filter=="r") {
        result_genes <- lr_network %>% pull(to) %>% unique()
    }
    result_genes
}

#' @noRd
get_expressed_receptors <- function(receivers_seuobj, expression_proportion_cutoff, lr_network) {
    expressed_genes_receivers <- get_expressed_genes(Idents(receivers_seuobj) %>% unique(),
                                                                receivers_seuobj,
                                                                pct = expression_proportion_cutoff)
    all_receptors <- lr_network$to %>% unique()
    expressed_receptors <- base::intersect(expressed_genes_receivers, all_receptors)
}

#' @noRd
get_active_ligand_receptor <- function(receivers_seuobj, expressed_ligands, expression_proportion_cutoff = 0.10, lr_network, ligand_target_matrix) {
    
    #subset to genes that are expressed in at least <expression_proportion_cutoff>% of the cells
    expressed_receptors <- get_expressed_receptors(receivers_seuobj, expression_proportion_cutoff, lr_network)
    lr_network_expressed <- lr_network %>% dplyr::filter(from %in% expressed_ligands & to %in% expressed_receptors) 
    assertthat::assert_that(nrow(lr_network_expressed) > 0, msg = paste("No active ligand-receptor pairs found. ",
                                                                      "To fix this, reduce the value of expression_proportion_cutoff",
                                                                      "or check that your sender and receiver cells are in fact communicating",
                                                                      "or check that your lr_network gene names correspond to the right species."))
    lr_network_expressed
}

#' @noRd
get_ligand_trajectory_scores_regression <- function(pseudotime_genes, active_ligand_potentials){
    assertthat::assert_that(ncol(active_ligand_potentials) > 0, msg = paste("Subsetted ligand_target_matrix has no ligands",
                                                                            "Check that your thresholds have captured active ligands"))
    
    intersection_genes <- base::intersect(rownames(active_ligand_potentials), row.names(pseudotime_genes))
    potentials_intersection <- active_ligand_potentials[intersection_genes,]
    pseudotime_genes_intersection <- pseudotime_genes[intersection_genes,]
    rf <- randomForest::randomForest(x=potentials_intersection, y=pseudotime_genes_intersection, importance=TRUE)
    
    list(genes=pseudotime_genes_intersection, model=rf)
}

#' @noRd
get_ligand_trajectory_scores <- function(pseudotime_associated_genes_bool, active_ligand_potentials, importance_measure = "MeanDecreaseGini"){
    intersection_genes <- base::intersect(rownames(active_ligand_potentials), names(pseudotime_associated_genes_bool))
    
    potentials_intersection<-active_ligand_potentials[intersection_genes,]
    pseudotime_associated_genes_intersection<-as.factor(pseudotime_associated_genes_bool[intersection_genes])
    rf <- randomForest::randomForest(x=potentials_intersection, y=pseudotime_associated_genes_intersection, importance=TRUE)
    
    importances<-sort(rf$importance[,importance_measure], decreasing=TRUE)
    list(importances, pseudotime_associated_genes_bool, rf)
}

#' @noRd
add_pseudotime_data_to_seuobj <- function(cds,
                                          obj,
                                          reduction_method = "UMAP") {
    obj@meta.data$pseudotime <- cds@principal_graph_aux$UMAP$pseudotime
    obj@meta.data$vertex <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex[,1]
    obj@misc$entrain$dp_mst <- cds@principal_graph_aux[[reduction_method]]$dp_mst
    return(obj)
}

#' @noRd
add_branch_data_to_seuobj <- function(cds, 
                                      obj,
                                      reduction_method = "UMAP") {
    
    
    if (is.null(obj@misc$entrain$dp_mst) == TRUE) {
        obj@misc$entrain$dp_mst <- list()
    }
    
    obj@meta.data$branch <- cds@colData$branch
    return(obj)
}

#' @noRd
get_background_genes_path <- function(path_expr, pct = 0.05) {
    all_genes <- row.names(path_expr)
    
    n_cells <- ncol(path_expr)
    n_cells_expressing<-Matrix::rowSums( path_expr > 0 )  
    expressed_genes <- names(n_cells_expressing[n_cells_expressing > pct*n_cells])
}

#Get UMAP xy coordinates of a principal node. Needed to know where to plot ligand labels.
#' @noRd
get_node_xy <- function(node_name, obj) {
    node_xy <- t(obj@misc$entrain$dp_mst) %>%
        as.data.frame() %>%
        dplyr::select(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2)
    x=node_xy[node_name,1]
    y=node_xy[node_name,2]
    return(c(x,y))
}

# Get UMAP xy coordinates of a cell. Needed to know where to plot ligand labels.
#' @noRd
get_cell_xy <- function(cell_name, obj, reduction_name) {
    xy <- obj[[reduction_name]]@cell.embeddings %>% as.data.frame()
    x = xy[cell_name, 1]
    y = xy[cell_name, 2]
    return(c(x,y))
}

#' @title Interactive selection of sender clusters. 
#' @description Choosing potential sender cell clusters in a GUI. Inspired from Monocle3's interactive workflow. 
#' @param sender_obj Seurat object containing niche cell data
#' @param sender_cluster_key column name in sender_obj@meta.data that denotes the clusters.
#' @param reduction_key Seurat reduction key for dimension reduction visualization.
#' @param lr_network NicheNet ligand-receptor pairs data file.
#' @return A vector of cluster names.
#' @import dplyr
#' @import shiny

#' @export
get_senders_interactive <- function(sender_obj,
                                    sender_cluster_key,
                                    reduction_key = NULL,
                                    lr_network
                                    ) {
    
    if (is.null(reduction_key)) {
        reduction_key <- sender_obj@reductions[[1]]@key
        message(paste("No Reduction key supplied, defaulting to", reduction_key))
    }
    
    cluster_df <- data.frame(matrix(NA, nrow = length(unique(sender_obj@meta.data[[sender_cluster_key]])), ncol=2))
    colnames(cluster_df) <- c("cluster", "colour")
    cluster_df$cluster <- sender_obj@meta.data[[sender_cluster_key]] %>% unique() %>% sort()
    cluster_df$colour <- scales::hue_pal()(nrow(cluster_df))
    
    ui <- shiny::fluidPage(
        shiny::titlePanel("Choose sender clusters interactively."),
        
        # Sidebar layout with input and output definitions ----
        shiny::sidebarLayout(
            
            # Sidebar panel for inputs ----
            shiny::sidebarPanel(
                # Visualization controls
                shiny::actionButton("zoom", "Zoom to selection"),
                shiny::actionButton("reset_zoom", "Reset Zoom"),
                shiny::br(),
                shiny::checkboxGroupInput("selected_clusters", "Selected Clusters:",
                                   choices = cluster_df$cluster),
                shiny::actionButton("preview", "Preview expressed ligands"),
                
                shiny::numericInput("proportion_cutoff", "Expression proportion cutoff:", 
                             0.01, min = 0, max = 1,
                             step = 0.01),
                
                shiny::actionButton("clear", "Clear all"),
                # done button
                shiny::actionButton("done", "Done"),
                shiny::br(),
                shiny::h3("Instructions:"),
                shiny::tags$ol(
                    shiny::tags$li("Select clusters using checkboxes."),
                    shiny::tags$li("Click Preview to view ligands expressed by your clusters"),
                    shiny::tags$li("Click 'Done' to return the selected ligands.")
                )
            ),
            
            # Main panel for displaying outputs ----
            shiny::mainPanel(
                shiny::plotOutput("plot1", height="auto",
                                  click = "plot1_click",
                                  brush = shiny::brushOpts(id = "plot1_brush", resetOnNew=TRUE)),
                #shiny::uiOutput("ligands", height="auto",
                #                resetOnNew=TRUE)
                DT::dataTableOutput("ligands", width = 500)
            )
        )
    )
    
    server <- function(input, output, session) {

        limits <- shiny::reactiveValues(
            x_range=NULL,
            y_range=NULL
        )
        
        selected_df <- shiny::reactive({
            cluster_selection <- input$selected_clusters
            cluster_df %>% dplyr::filter(cluster %in% cluster_selection)
        })
        
        expressed_ligands <- shiny::reactiveValues(
            expressed_ligands_df = NULL,
            expressed_ligands_df_visual = NULL
        )
        observeEvent(input$preview, {
            updated <- update_ligands_and_plot()
            expressed_ligands$expressed_ligands_df <- updated[[1]]
            expressed_ligands$expressed_ligands_df_visual <- updated[[2]]
            
            output$ligands <- DT::renderDataTable(
                expressed_ligands$expressed_ligands_df_visual,
                options = list(scrollY="300px", 
                               scrollX=TRUE,
                               pageLength = 30, info = FALSE,
                               lengthMenu = list(c(30, -1), c("30", "All")))
            )
            
        })
            
        output$plot1 <- shiny::renderPlot({
            
            cluster_df[cluster_df$cluster %in% selected_df()$cluster, "colour"] <- "#3B3B3B"
            colourscale <- cluster_df$colour
            names(colourscale) <- cluster_df$cluster
            seurat_plot <- Seurat::DimPlot(sender_obj,
                                           reduction = reduction_key,
                                           group.by = sender_cluster_key,
                                           cols = colourscale,
                                           combine=FALSE)[[1]]
            
            suppressMessages(seurat_plot + 
                ggplot2::coord_cartesian(xlim = limits$x_range, ylim = limits$y_range, expand = FALSE)
            )

        }, height = function() {
            session$clientData$output_plot1_width
        })
        
        if (is.null(expressed_ligands)) {
            output$ligands <- shiny::renderText( {
                "Click preview to view expressed ligands"}
            )
        }
       
    
        observeEvent(input$clear, {
            updateCheckboxGroupInput(session, "selected_clusters",
                                     choices=cluster_df$cluster,
                                     selected=NULL)
            expressed_ligands$expressed_ligands_df = NULL
            expressed_ligands$expressed_ligands_df_visual = NULL
        })
        
        shiny::observeEvent(input$done, {
            shiny::stopApp( update_ligands_and_plot()[[1]]) # value here is the returned data (expressed ligands)
        })

        shiny::observeEvent(input$zoom, {
            brush <- input$plot1_brush
            if (!is.null(brush)) {
                limits$x_range <- c(brush$xmin, brush$xmax)
                limits$y_range <- c(brush$ymin, brush$ymax)

            } else {
                limits$x_range <- NULL
                limits$y_range <- NULL
            }
        })

        shiny::observeEvent(input$reset_zoom, {
            limits$x_range <- NULL
            limits$y_range <- NULL
        })
        
        update_ligands_and_plot = function() {
            cutoff <- input$proportion_cutoff
            sender_cluster_names <- selected_df()$cluster
                
            if (length(sender_cluster_names) > 0) {
                expressed_ligands_df <- data.frame(matrix(0, nrow=0, ncol=0))
                expressed_ligands_df_visual <- data.frame(matrix(0, nrow=1, ncol=0))
                for (cluster in sender_cluster_names) {
                    cluster_ligands <- get_expressed_ligands(sender_obj, sender_cluster_key = sender_cluster_key,
                                                                   sender_cluster_names = cluster,
                                                                   proportion_cutoff = cutoff,
                                                                   lr_network = lr_network) %>%
                    data.frame() %>% setNames("ligands")
                    expressed_ligands_df_visual[[cluster]] <- toString(cluster_ligands[,1])
                    expressed_ligands_df <- rbind(expressed_ligands_df, cluster_ligands)
                }
                    #expressed_ligands_df <- t(expressed_ligands_df)
                rownames(expressed_ligands_df_visual) <- NULL
                return(list(expressed_ligands_df, expressed_ligands_df_visual))
            } else {
                return(NULL)
            }            
        }
    }
    sel <- shiny::runApp(shiny::shinyApp(ui, server))
}

#' @noRd
get_expressed_ligands <- function(obj, sender_cluster_key,
                             sender_cluster_names,
                             proportion_cutoff = 0.01,
                             lr_network) {
    ligands <- lr_network$from %>% unique()
    Seurat::Idents(obj) <- sender_cluster_key
    expressed_genes_senders <- get_expressed_genes(ident = sender_cluster_names,
                                                              seurat_obj = obj,
                                                              pct = proportion_cutoff)
    expressed_ligands = intersect(ligands,expressed_genes_senders)
    return(expressed_ligands)
}

#' @title Identify ligands responsible for trajectory dynamics using Monocle graph learning.
#' @description Runs a standard Monocle3 trajectory reconstruction pipeline, followed by Entrain to identify ligands responsible for trajectory dynamics.
#' @param obj A seurat object containing receiver cells undergoing a trajectory
#' @param sender_obj A seurat obj containing sender cells expressing ligands. Set to NULL if your sender cells are contained in 'obj'.
#' @param sender_cluster_key Column name in metadata that denotes where sender cells are specified. If NULL, this defaults to the default Seurat Idents given by `levels(obj)`
#' @param sender_cluster_names Unique values in column `sender_cluster_key` that denote sender cell clusters.
#' @param cds Optional. You can supply a `monocle3` `cds` object here if you have run `monocle3::order_cells()` on it already.
#' @param num_dim Number of dimensions for monocle reduce_dimension function.
#' @param root_cells Optional. Character vector of root cell names that comprise the beginning of the trajectory (must be in colnames(obj)). If not given or NULL, will open a window to interactively select a root cell (Monocle's default behaviour).
#' @param path_cell_names Optional. Character vector of cell names that comprise the differentiation trajectory to be analyzed by Entrain (must be in colnames(obj)). If not given or NULL, will open a window to interactively select the trajectory path to analyze (Monocle's default behaviour).
#' @param root_pr_nodes Argument for `monocle3::order_cells()`. Use if you want to specify root nodes programmatically.
#' @param path_nodes Use if you want to specify the path without user interaction. e.g. `path_nodes = c("Y_2", "Y_50")` tells Entrain to analyze cells between the nodes `Y_2` and `Y_50`
#' @param expressed_ligands Character vector of active ligands in dataset (must be in rownames(obj))
#' @param reduction_key Seurat reduction key for dimension reduction visualization.
#' @param ncenter Parameter for Monocle graph learning.
#' @param use_partition Parameter for Monocle graph learning
#' @param prune_graph Parameter for Monocle graph learning.
#' @param close_loop Parameter for Monocle graph learning.
#' @param ncenter Parameter for Monocle graph learning.
#' @param minimal_branch_len Parameter for Monocle graph learning.
#' @param overwrite If there is an existing Monocle trajectory in the obj (e.g. if you have run this function before), overwrite = TRUE will delete the old trajectory. Otherwise, the old trajectory will be used for Entrain.
#' @param export_cds If TRUE, saves the Monocle cell_data_set trajectory object in obj@misc$monocle_graph. Set TRUE if you want to plot Entrain results afterwards.
#' @param lr_network NicheNet ligand-receptor pairs data file.
#' @param ligand_target_matrix NicheNet ligand-target data file.
#' @param expression_proportion_cutoff Pct cutoff to threshold a ligand as active. A ligand is 'active' if it is expressed in more than `expression_proportion_cutoff` fraction of cells in the sender cluster.
#' @param covariance_cutoff Remove the bottom `covariance_cutoff` fraction of covariances (e.g. 0.10 = bottom 10 percent of genes, ranked by covariance, removed from later analysis).
#' @param precomputed_umap If you have a precomputed UMAP on your receiver cells that you wish to build a monocle trajectory on, then set this to the name of your precomputed embedding e.g. precomputed_umap = 'umap'
#' @param sender_reduction_key Only applicable if sender_obj is not NULL: The dimension reduction key of sender_obj that you wish to use,
#' @import dplyr
#' @return a Seurat object with ligand velocity results in obj$misc$entrain$velocity_result. 
#' @export
#' 
get_traj_ligands_monocle <- function(obj, sender_obj = NULL,
                                     cds = NULL,
                                     sender_cluster_key = NULL,
                                     sender_cluster_names = NULL,
                                     num_dim=10,
                                     root_cells = NULL,
                                     root_pr_nodes = NULL,
                                     expressed_ligands = NULL,
                                     path_cell_names = NULL,
                                     path_nodes = NULL,
                                     minimal_branch_len = 5,
                                     ncenter = 50,
                                     use_partition = TRUE,
                                     precomputed_umap=NULL,
                                     prune_graph=T, close_loop = F,
                                     overwrite=FALSE, 
                                     reduction_key = "MonocleUMAP_",
                                     sender_reduction_key = "umap",
                                     expression_proportion_cutoff = 0.01,
                                     covariance_cutoff = 0.05,
                                     lr_network, ligand_target_matrix, export_cds=TRUE
) {
    if (!requireNamespace("monocle3", quietly = TRUE) | requireNamespace("SeuratWrappers", quietly = TRUE) | requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("Please install monocle3, SeuratWrappers and SingleCellExperiment to use get_traj_ligands_monocle")  
    }
    
    if (is.null(obj@misc$monocle_graph) || overwrite==TRUE ) {
        if (!is.null(obj@misc$entrain$paths)) {
            message("You are generating a new trajectory but existing Entrain analysis found. The existing analysis will be overwritten.")
            obj@misc$entrain$paths <- NULL
            obj@misc$entrain$dp_mst <- NULL
        }
        # Receiver trajectory
        if (is.null(cds)) {
            cds <- SeuratWrappers::as.cell_data_set(obj,  assays = DefaultAssay(obj))
            #below line inserts Seurat object original umap embedding into the monocle workflow.
            cds <- monocle3::preprocess_cds(cds, num_dim=num_dim)
            
            if (!is.null(precomputed_umap)) {
                cds@int_colData@listData$reducedDims$UMAP <- obj@reductions[[precomputed_umap]]@cell.embeddings
            } else {
                cds <- monocle3::reduce_dimension(cds)
                message(paste("Adding Monocle reductions to Seurat object under reduction key: ", reduction_key,sep=""))
                cds_embedding <- cds@int_colData@listData$reducedDims$UMAP
                colnames(cds_embedding) <- c("Monocle_1", "Monocle_2")
                obj[[reduction_key]]<-Seurat::CreateDimReducObject(embeddings = cds_embedding , key = reduction_key, assay = DefaultAssay(obj))
                
            }

            cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP")

            graph_settings=list(minimal_branch_len=minimal_branch_len, ncenter=ncenter, prune_graph=prune_graph)
            cds_graph <- monocle3::learn_graph(cds,
                                               use_partition = use_partition,
                                               close_loop=close_loop,
                                               learn_graph_control=graph_settings
            )
            
            if (is.null(root_cells) == TRUE & is.null(root_pr_nodes) == TRUE ) {
                cds <- monocle3::order_cells(cds_graph)
            } else if (!is.null(root_cells) & is.null(root_pr_nodes)) {
                message("Ordering cells with argument root_cells")
                cds <- monocle3::order_cells(cds_graph, root_cells = root_cells)
            } else if (!is.null(root_pr_nodes) & is.null(root_cells) ) {
                message("Ordering cells with argument root_pr_nodes")
                cds <- monocle3::order_cells(cds_graph, root_pr_nodes = root_pr_nodes)
            } else {
                stop("Please supply one of root_cells or root_pr_nodes as arguments, not both.")
            }
        }
        else {
            message("Using pre-generated Monocle trajectory in argument cds.")
            if (is.null(root_cells) == TRUE & is.null(root_pr_nodes) == TRUE ) {
                cds <- monocle3::order_cells(cds_graph)
            } else if (!is.null(root_cells) & is.null(root_pr_nodes) ) {
                message("Ordering cells with argument root_cells. If you have already run monocle3::order_cells(), set root_cells = NULL and root_pr_nodes = NULL.")
                cds <- monocle3::order_cells(cds, root_cells = root_cells)
            } else if (!is.null(root_pr_nodes) & is.null(root_cells) ) {
                message("Ordering cells with argument root_pr_nodes. If you have already run monocle3::order_cells(), set root_cells = NULL and root_pr_nodes = NULL.")
                cds <- monocle3::order_cells(cds, root_pr_nodes = root_pr_nodes)
            } else {
                stop("Please supply one of root_cells or root_pr_nodes as arguments, not both.")
            }
            message(paste("Adding Monocle reductions to Seurat object under reduction key: ", reduction_key,sep=""))
            cds_embedding <- cds@int_colData@listData$reducedDims$UMAP
            colnames(cds_embedding) <- c("Monocle_1", "Monocle_2")
            obj[[reduction_key]]<-Seurat::CreateDimReducObject(embeddings = cds_embedding , key = reduction_key, assay = DefaultAssay(obj))
        }
        #cds <- get_branch_membership(cds)
        #obj <- add_branch_data_to_seuobj(cds, obj)
        obj <- add_pseudotime_data_to_seuobj(cds, obj)
        obj@misc$monocle_graph <- cds@principal_graph_aux$UMAP
    }
    else {
        message("Existing Monocle trajectory found, skipping graph learning. If you wish to build a new Monocle trajectory, set overwrite=TRUE")
        cds <- obj@misc$entrain$monocle_cds
        
    }
    
    if (identical(obj, sender_obj)) {
        stop("obj is the same as sender_obj. Set sender_obj = NULL if your sender and receiver cells are contained in the single seurat object.")
    } else if (is.null(sender_obj)) {
        sender_obj <- obj
        sender_reduction_key = reduction_key
    }
    
    if (is.null(expressed_ligands)) {
        message("expressed_ligands is NULL, selecting sender clusters")
        
        if (is.null(sender_cluster_key)) {
            sender_obj <- 
            stop("No sender_cluster_key supplied. Please supply a sender_cluster_key: The column name in metadata that denotes sender cell clusters.")
        }
        
        if (is.null(sender_cluster_names)) {
            message("sender_cluster_names is NULL, selecting sender clusters interactively.")
            assertthat::assert_that(!is.null(lr_network), msg = "If you are in interactive mode, please supply lr_network to get_traj_ligands_monocle()")
            expressed_ligands <- get_senders_interactive(sender_obj, 
                                                         sender_cluster_key,
                                                         reduction_key = sender_reduction_key,
                                                         lr_network = lr_network)
            expressed_ligands <- expressed_ligands 
        } else {

            expressed_ligands <- get_expressed_ligands(obj = sender_obj, 
                                                       sender_cluster_key = sender_cluster_key,
                                                       sender_cluster_names = sender_cluster_names,
                                                       proportion_cutoff = expression_proportion_cutoff,
                                                       lr_network = lr_network)
        }

    }
    if (is.null(path_cell_names) & is.null(path_nodes)) {
        selected_paths <- choose_path(cds, return_list=TRUE)
        for (path in selected_paths) {
            path_cell_names <- path$cells
            obj <- get_path_ligands(obj, expressed_ligands=expressed_ligands,  
                                  path_cell_names = path_cell_names,
                                  return_model=TRUE,
                                  reduction_name = reduction_key,
                                  lr_network=lr_network, ligand_target_matrix=ligand_target_matrix,
                                  covariance_cutoff = covariance_cutoff
                                  ) 
        }
    } else {
        if ( !is.null(path_nodes) & is.null(path_cell_names) ) {
            path_cell_names = cells_between_vertices(endpoint_nodes = path_nodes,
                                                     cds = cds)
            message("Non-interactive mode: Calculating for the branch specified by cells in argument path_nodes")
        } else if (is.null(path_nodes) & !is.null(path_cell_names) ) {
            message("Non-interactive mode: Calculating for the branch specified by cells in argument path_cell_names")
        } else {
            stop("Please supply only one of path_nodes and path_cell names")
        }
    
        obj<-get_path_ligands(obj, expressed_ligands=expressed_ligands,
                                  path_cell_names = path_cell_names,
                                  return_model=TRUE,
                                  reduction_name = reduction_key,
                                  lr_network=lr_network, ligand_target_matrix=ligand_target_matrix,
                                  covariance_cutoff = covariance_cutoff
            )
    }
  
    if (export_cds==TRUE) {
        message("Exporting cds...")
        obj@misc$entrain$monocle_cds <- cds
    }
    return(obj)
}

#' @title Core function for Entrain ligand assignment that is agnostic to the trajectory inference algorithm used.
#' @description Runs Entrain to identify ligands responsible for trajectory dynamics. Requires pseudotime labels for each cell in a column in metadata with name pseudotime_key
#' @param obj A seurat object  with cell pseudotimes contained in a metadata column with column name pseudotime_key.
#' @param pseudotime_key Column name in colnames(obj@meta.data) that contains cell pseudotime values.
#' @param cluster_key Column name in colnames(obj@meta.data) corresponding to 
#' @param path_cell_names a vector of cell names (must be in colnames(obj)), defining cells in a trajectory path for which associated ligands will be calculated.
#' @param return_model If TRUE, will additionally save the raw random forest object in obj@misc$entrain$paths$path_name$model
#' @param expressed_ligands Character vector of active ligands in dataset (must be in rownames(obj))
#' @param lr_network NicheNet ligand-receptor pairs data file.
#' @param ligand_target_matrix NicheNet ligand-target data file.
#' @param reduction_name Seurat reduction key for dimension reduction visualization.
#' @param covariance_cutoff Remove the bottom `covariance_cutoff` fraction of covariances (e.g. 0.10 = bottom 10 percent of genes, ranked by covariance, removed from later analysis).
#' @import dplyr
#' @return a Seurat object with ligand trajectory results in obj$misc$entrain$paths$path_name
#' @export
get_path_ligands <- function(obj, expressed_ligands, 
                             lr_network, ligand_target_matrix, 
                             path_cell_names, 
                             pseudotime_key="pseudotime",
                             cluster_key = NULL,
                             reduction_name = "MonocleUMAP",
                             covariance_cutoff = 0.05,
                             return_model = FALSE) {
    
    path_seuobj <- subset(obj, cells=path_cell_names)
    path_pseudotime <- path_seuobj@meta.data[[pseudotime_key]]

    if (!is.null(path_seuobj@meta.data$vertex)) {
        start_node <- path_seuobj@meta.data[path_seuobj@meta.data[,pseudotime_key]==min(path_seuobj@meta.data[,pseudotime_key]),"vertex"] %>% as.character()
        end_node <- path_seuobj@meta.data[path_seuobj@meta.data[,pseudotime_key]==max(path_seuobj@meta.data[,pseudotime_key]),"vertex"] %>% as.character()
    }
    
    if (!is.null(cluster_key)) {
        start_cluster <- path_seuobj@meta.data[path_seuobj@meta.data[,pseudotime_key]==min(path_seuobj@meta.data[,pseudotime_key]), cluster_key, drop=FALSE]
        end_cluster <- path_seuobj@meta.data[path_seuobj@meta.data[,pseudotime_key]==max(path_seuobj@meta.data[,pseudotime_key]), cluster_key, drop=FALSE] 
        path_name <- paste(start_cluster, end_cluster, sep="-")
        
        if (!is.null(path_seuobj@meta.data$vertex)) {
            path_name <- paste(path_name,"_", start_node[1],"-", end_node[1], sep="")
        }
        
    } else if (!is.null(path_seuobj@meta.data$vertex)) {
        path_name <- paste(start_node[1],"-", end_node[1], sep="")
    } else {
        stop("get_path_ligands requires either a pseudotime column passed to pseudotime_key or a cell cluster column passed to cluster_key")
    }
    

    
    midpoint_cell <- path_seuobj@meta.data[path_seuobj@meta.data$pseudotime==sort(path_pseudotime)[round(length(path_pseudotime)/2)],] %>% row.names()  %>% as.character()
    midpoint_xy <- get_cell_xy(midpoint_cell, path_seuobj, reduction_name)
    
    message(paste("Running Entrain on path:", path_name, sep=" "))
    
    n_inf_pseudotime_cells <-sum(is.infinite(path_pseudotime))
    if (n_inf_pseudotime_cells > 0 ) {
        message(paste(as.character(n_inf_pseudotime_cells), " cells in branch ", path_name[1], " possess Inf pseudotime. Removing cells from analysis. To fix this, reassign root nodes such that all cells are connected to a root node.", sep="") )
        non_inf_pseudotime_indices <- is.finite(path_pseudotime)
        path_seuobj <- path_seuobj[,non_inf_pseudotime_indices]
        path_pseudotime <- path_seuobj@meta.data[[pseudotime_key]]
    }
    
    path_expr <- as.matrix(GetAssayData(path_seuobj, slot = "data"))
    
    active_lr <- get_active_ligand_receptor(receivers_seuobj = path_seuobj, expressed_ligands = expressed_ligands, lr_network = lr_network, ligand_target_matrix = ligand_target_matrix)
    active_ligands<-active_lr$from %>% unique()
    active_receptors<-active_lr$to %>% unique()
    active_ligand_potentials <- ligand_target_matrix[,active_ligands]
    
    background_genes <- get_background_genes_path(path_expr, pct = 0.10)
    path_expr_bg<-path_expr[background_genes,]
    pseudotime_covars <- get_branch_genes_covar(path_expr_bg, path_seuobj@meta.data$pseudotime)
    
        # Bottom 5% of covariances (by absolute value) likely noise.
    cutoff<-quantile(pseudotime_covars[,"Covariances"], covariance_cutoff) 
    top_pseudotime_covars <- subset(pseudotime_covars, Covariances > cutoff) 
    ligand_scores_result<-get_ligand_trajectory_scores_regression(top_pseudotime_covars, active_ligand_potentials)
    names(ligand_scores_result) <- c(path_name[1], "model")
    importances<-sort(ligand_scores_result$model$importance[,"IncNodePurity"], decreasing=TRUE)
    
    ## add to seurat
    if (is.null(obj@misc$entrain$paths) == TRUE) {
        obj@misc$entrain$paths <- list()
    }
    
    path_cell_names<-colnames(path_seuobj)
    results <- list(cell_pseudotimes = data.frame(path_cell_names, 
                                                  path_pseudotime), 
                    pseudotime_associated_genes = pseudotime_covars,
                    top_pseudotime_genes = rownames(top_pseudotime_covars),
                    path_ligands_importances = importances,
                    path_midpoint = list(midpoint_cell, midpoint_xy),
                    active_network = active_lr)
    
    if (return_model == TRUE) {
        results[["model"]]=ligand_scores_result$model
    }
    
    if (is.null(obj@misc$entrain$paths[[path_name]])) {
        obj@misc$entrain$paths[[path_name]] <- results
    } else {
        i <- 1
        while(!is.null(obj@misc$entrain$paths[[path_name]])) {
            i <- i+1
            path_name <- paste(path_name,i,sep="_")
        }
        message(paste("Duplicate path name found. Appending current path as: ", path_name) )
        obj@misc$entrain$paths[[path_name]] <- results
    }
    
    return(obj)
}



#' @title Evaluates where ligands are most strongly driving trajectory dynamics.
#' @description Performs a sliding window analysis along the trajectory to identify windows where environmental influence is most apparent.
#' @param obj A seurat object that has been analyzed with get_path_ligands or get_traj_ligands_monocle.
#' @param ligand_target_matrix NicheNet ligand-target data file.
#' @param ligands_list Optional, a character vector denoting ligands. If supplied, `cellwise_influences` will calculate influences for the ligands in `ligands_list` instead of the top ranked ligands.
#' @param step_size Size of the step, or interval, as a percentage of the number of cells in a path.
#' @param window_pct Size of the window, as a percentage of the number of cells in a path.
#' @param n_top_ligands Number of top ligands to consider in the cellwise analysis.
#' @param n_top_genes Number of top dynamical genes to consider in the cellwise analysis. Because we have already evaluated the ligand importances for a branch, we restrict our more granular analysis to the top genes and ligands.
#' @param slot Slot of Seurat object from which to retrieve expression data. 

#' @return a Seurat object updated with per-cell ligand influence score in `obj@misc$entrain$path$path_name$cellwise_influences`, as well as a column in `obj@meta.data`
#' @export

cellwise_influences <- function(obj,
                                ligand_target_matrix,
                                ligands_list=NULL,
                                step_size=0.10,
                                window_pct=0.30,
                                n_top_ligands=5,
                                n_top_genes=500,
                                slot="data"
) {
    if (!requireNamespace("imputeTS", quietly = TRUE)) {
        stop(
            "Package \"imputeTS\" must be installed to use this function.",
            call. = FALSE
        )
    }
    if (is.null(ligand_target_matrix)) {
        stop("Please supply a ligand_target_matrix.")
    }
    
    paths<-names(obj@misc$entrain$paths)
    
    for (path in paths) {
        pathdata<-obj@misc$entrain$paths[[path]]
        #genes<-pathdata$pseudotime_associated_genes %>% arrange(. ,desc(.)) %>% rownames() %>% head(n_top_genes)
        #genes<-genes[genes %in% rownames(ligand_target_matrix)]
        pseudotimes<-pathdata$cell_pseudotimes
        pseudotimes <- pseudotimes %>% arrange(path_pseudotime)
        half_win <- pseudotimes %>% nrow(.)*(window_pct/2)
        half_win <- round(half_win)
        
        if (step_size < 1) {
            i <- round(nrow(pseudotimes)*step_size)
        } else {
            i <- step_size
        }
        
        # get expr matrix
        path_obj<-subset(obj, cells=pseudotimes$path_cell_names)
        path_expr <- Seurat::GetAssayData(path_obj, slot = slot)
        if (is.null(ligands_list) == TRUE) {
            importances <- pathdata$path_ligands_importances %>% head(n_top_ligands)
            meta_colname <- paste("Influences", gsub("-", ".", path), sep="_")
            cellwise_influences <- data.frame(matrix(nrow=nrow(pseudotimes), ncol=1+n_top_ligands))
            colnames(cellwise_influences) <- c(names(importances), meta_colname )
        } else {
            ligands <- ligands_list[[path]]
            importances <- pathdata$path_ligands_importances[ligands]
            message(paste("Calculating cellwise influences for ligands: ", paste(ligands,collapse=", "), sep=""))
            meta_colname <- paste("Influences", paste(ligands,collapse="_"), sep="_")
            cellwise_influences <- data.frame(matrix(nrow=nrow(pseudotimes), ncol=1+length(ligands)))
            colnames(cellwise_influences) <- c(names(importances), meta_colname )
        }
        
        rownames(cellwise_influences) <- pseudotimes$path_cell_names
        
        n=1
        while (n<nrow(pseudotimes)+1) {
            if (n<half_win) {
                window_pseudotimes<-pseudotimes[seq(1,n+half_win),]
            } else if (n > (nrow(pseudotimes) - half_win) ) {
                window_pseudotimes<-pseudotimes[seq(n-half_win, nrow(pseudotimes)),]
            } else {
                window_pseudotimes=pseudotimes[seq(n-half_win,n+half_win),]
            }
            
            window_cells<-window_pseudotimes$path_cell_names
            background_genes <- get_background_genes_path(path_expr, pct = 0.10)
            
            window_expr <- path_expr[background_genes,window_cells]
            
            #get correlation for each row. 
            window_expr <- as.matrix(window_expr)
            
            covars <- get_branch_genes_covar(window_expr, window_pseudotimes$path_pseudotime)
            
            if (sum(is.na(covars)) > 0) {
                covars <- na.omit(t(covars)) 
                covars<-covars[,1]
            }
            covars <- abs(covars)
            covars <- covars[order(covars[,1]),,drop=FALSE] %>% tail(n_top_genes)
            genes<-rownames(covars)
            genes<-base::intersect(genes, rownames(ligand_target_matrix))
            
            w<-ligand_target_matrix[genes, names(importances)]
            if (length(importances) == 1){
                w <- as.matrix(w)
            }

            covars<-covars[genes,]
            covars<-(covars-min(covars))/(max(covars)-min(covars))
            names(covars) <- genes
            w<-(w-min(w))/(max(w)-min(w))

            rf = randomForest::randomForest(y=covars, x=w, importance=TRUE)
            
            rsq<-rf$rsq %>% tail(1)
            imps<-rf$importance[,"IncNodePurity"]
            
            cellwise_influences[n,]<-c(imps, rsq)
            n<-n+i
        }
        
        if (i > 1) {
            cellwise_influences <- apply(cellwise_influences,
                                         MARGIN=2,
                                         imputeTS::na_interpolation, option ="linear")
        }

        # store data
        obj@misc$entrain$paths[[path]]$cellwise_influences <- cellwise_influences

        overwrite_cols <- colnames(cellwise_influences)[colnames(cellwise_influences) %in% colnames(obj@meta.data)]
        if ( length(overwrite_cols) > 0  ) {
            message("Existing cellwise-influences data found. Overwriting...")
            metadata <- obj@meta.data %>% select(-one_of(overwrite_cols))
        } else {
            metadata <- obj@meta.data
        }
        temp <- merge(metadata, cellwise_influences, by.x=0, by.y=0, all.x=TRUE)
        temp <- transform(temp, row.names=Row.names, Row.names=NULL)
        temp[is.na(temp)] <- 0
        
        # preserve row order with original object
        temp<-temp[match(rownames(temp), rownames(metadata)),]
        
        obj@meta.data <- temp
        
    }
    
    return(obj)
}

#' @title Intersect seurat object genes
#' @description Subsets the genes of obj to the genes that are the intersection of obj and the anndata h5ad
#' @param obj Seurat object to be subsetted
#' @param h5ad String denoting filename of h5ad from with to intersect genes with.
#' @return A seurat object with features subsetted to the intersection of genes in obj and h5ad.
intersect_genes = function(obj, h5ad) {
    
    anndata <- reticulate::import("anndata")
    ad <- anndata$read_h5ad(h5ad)
    ad_genes <- ad$var_names$values
    obj_genes <- rownames(obj)
    intersection <- intersect(ad_genes, obj_genes)
    
    result_obj <- subset(obj, features=intersection)
    return(result_obj)
}

#' @title Determine expressed genes of a cell type.
#' @description  NicheNetR function to determine expressed genes of a cell type from a Seurat object single-cell RNA seq dataset or Seurat spatial transcriptomics dataset. 
#' Function is modified from `nichenetr` v1.1.1 (GPL3.0): Browaeys, R., Saelens, W. & Saeys, Y. (2019). Nat Methods.
#' Modified to deal with non-sparse matrices, and avoids the (non-trivial) installation of NicheNet.
#' @param ident Name of cluster identity/identities of cells
#' @param seurat_obj Single-cell expression dataset as Seurat object https://satijalab.org/seurat/.
#' @param pct We consider genes expressed if they are expressed in at least a specific fraction of cells of a cluster. This number indicates this fraction. Default: 0.10. Choice of this parameter is important and depends largely on the used sequencing platform. We recommend to require a lower fraction (like the default 0.10) for 10X data than for e.g. Smart-seq2 data.
#' @param assay_oi If wanted: specify yourself which assay to look for. Default this value is NULL and as a consequence the 'most advanced' assay will be used to define expressed genes.
#'
#' @return A character vector with the gene symbols of the expressed genes
#'
#' @import Seurat
#' @import dplyr
#'
#' @export
#'
get_expressed_genes = function(ident, seurat_obj, pct = 0.1, assay_oi = NULL){
    requireNamespace("Seurat")
    requireNamespace("dplyr")
    
    # input check
    
    if (!"RNA" %in% names(seurat_obj@assays)) {
        if ("Spatial" %in% names(seurat_obj@assays)) {
            if (class(seurat_obj@assays$Spatial@data) != "matrix" &
                class(seurat_obj@assays$Spatial@data) != "dgCMatrix") {
                warning("Spatial Seurat object should contain a matrix of normalized expression data. Check 'seurat_obj@assays$Spatial@data' for default or 'seurat_obj@assays$SCT@data' for when the single-cell transform pipeline was applied")
            }
            if (sum(dim(seurat_obj@assays$Spatial@data)) == 0) {
                stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$Spatial@data'")
            }
        }
    }
    else {
        if (!("dgCMatrix" %in% class(seurat_obj@assays$RNA@data)))
            warning("RNA matrix is not sparse. Converting to sparse..")
            seurat_obj@assays$RNA@data <- Seurat::as.sparse(seurat_obj@assays$RNA@data)
        
        if (class(seurat_obj@assays$RNA@data) != "matrix" &
            class(seurat_obj@assays$RNA@data) != "dgCMatrix") {
            warning("Seurat object should contain a matrix of normalized expression data. Check 'seurat_obj@assays$RNA@data' for default or 'seurat_obj@assays$integrated@data' for integrated data or seurat_obj@assays$SCT@data for when the single-cell transform pipeline was applied")
        }
        if ("integrated" %in% names(seurat_obj@assays)) {
            if (sum(dim(seurat_obj@assays$RNA@data)) == 0 & sum(dim(seurat_obj@assays$integrated@data)) ==
                0)
                stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA@data' for default or 'seurat_obj@assays$integrated@data' for integrated data")
        }
        else if ("SCT" %in% names(seurat_obj@assays)) {
            if (sum(dim(seurat_obj@assays$RNA@data)) == 0 & sum(dim(seurat_obj@assays$SCT@data)) ==
                0) {
                stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA@data' for default or 'seurat_obj@assays$SCT@data' for data corrected via SCT")
            }
        }
        else {
            if (sum(dim(seurat_obj@assays$RNA@data)) == 0) {
                stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA@data'")
            }
        }
    }
    if (sum(ident %in% unique(Idents(seurat_obj))) != length(ident)) {
        stop("One or more provided cell clusters is not part of the 'Idents' of your Seurat object")
    }
    
    if(!is.null(assay_oi)){
        if(! assay_oi %in% Seurat::Assays(seurat_obj)){
            stop("assay_oi should be an assay of your Seurat object")
        }
    }
    
    # Get cell identities of cluster of interest
    
    
    cells_oi = Seurat::Idents(seurat_obj) %>% .[Seurat::Idents(seurat_obj) %in%
                                            ident] %>% names()
    
    # Get exprs matrix: from assay oi or from most advanced assay if assay oi not specifcied
    
    if(!is.null(assay_oi)){
        cells_oi_in_matrix = intersect(colnames(seurat_obj[[assay_oi]]@data), cells_oi)
        exprs_mat = seurat_obj[[assay_oi]]@data %>% .[, cells_oi_in_matrix]
    } else {
        if ("integrated" %in% names(seurat_obj@assays)) {
            warning("Seurat object is result from the Seurat integration workflow. The expressed genes are now defined based on the integrated slot. You can change this via the assay_oi parameter of the get_expressed_genes() functions. Recommended assays: RNA or SCT")
            cells_oi_in_matrix = intersect(colnames(seurat_obj@assays$integrated@data),
                                           cells_oi)
            if (length(cells_oi_in_matrix) != length(cells_oi))
                stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$integrated@data). Please check that the expression matrix contains cells in columns and genes in rows.")
            exprs_mat = seurat_obj@assays$integrated@data %>% .[,
                                                                cells_oi_in_matrix]
        }
        else if ("SCT" %in% names(seurat_obj@assays) & !"Spatial" %in%
                 names(seurat_obj@assays)) {
            warning("Seurat object is result from the Seurat single-cell transform workflow. The expressed genes are defined based on the SCT slot. You can change this via the assay_oi parameter of the get_expressed_genes() functions. Recommended assays: RNA or SCT")
            cells_oi_in_matrix = intersect(colnames(seurat_obj@assays$SCT@data),
                                           cells_oi)
            if (length(cells_oi_in_matrix) != length(cells_oi))
                stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$SCT@data). Please check that the expression matrix contains cells in columns and genes in rows.")
            exprs_mat = seurat_obj@assays$SCT@data %>% .[, cells_oi_in_matrix]
        }
        else if ("Spatial" %in% names(seurat_obj@assays) &
                 !"SCT" %in% names(seurat_obj@assays)) {
            warning("Seurat object is result from the Seurat spatial object. The expressed genes are defined based on the Spatial slot. If the spatial data is spot-based (mixture of cells) and not single-cell resolution, we recommend against directly using nichenetr on spot-based data (because you want to look at cell-cell interactions, and not at spot-spot interactions! ;-) )")
            cells_oi_in_matrix = intersect(colnames(seurat_obj@assays$Spatial@data),
                                           cells_oi)
            if (length(cells_oi_in_matrix) != length(cells_oi))
                stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$Spatial@data). Please check that the expression matrix contains cells in columns and genes in rows.")
            exprs_mat = seurat_obj@assays$Spatial@data %>% .[, cells_oi_in_matrix]
        }
        else if ("Spatial" %in% names(seurat_obj@assays) &
                 "SCT" %in% names(seurat_obj@assays)) {
            warning("Seurat object is result from the Seurat spatial object, followed by the SCT workflow. If the spatial data is spot-based (mixture of cells) and not single-cell resolution, we recommend against directly using nichenetr on spot-based data (because you want to look at cell-cell interactions, and not at spot-spot interactions! The expressed genes are defined based on the SCT slot, but this can be changed via the assay_oi parameter.")
            cells_oi_in_matrix = intersect(colnames(seurat_obj@assays$SCT@data),
                                           cells_oi)
            if (length(cells_oi_in_matrix) != length(cells_oi))
                stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$Spatial@data). Please check that the expression matrix contains cells in columns and genes in rows.")
            exprs_mat = seurat_obj@assays$SCT@data %>% .[, cells_oi_in_matrix]
        }
        else {
            if (sum(cells_oi %in% colnames(seurat_obj@assays$RNA@data)) ==
                0)
                stop("None of the cells are in colnames of 'seurat_obj@assays$RNA@data'. The expression matrix should contain cells in columns and genes in rows.")
            cells_oi_in_matrix = intersect(colnames(seurat_obj@assays$RNA@data),
                                           cells_oi)
            if (length(cells_oi_in_matrix) != length(cells_oi))
                stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$RNA@data). Please check that the expression matrix contains cells in columns and genes in rows.")
            exprs_mat = seurat_obj@assays$RNA@data %>% .[, cells_oi_in_matrix]
        }
        
    }
    
    # use defined cells and exprs matrix to get expressed genes
    
    n_cells_oi_in_matrix = length(cells_oi_in_matrix)
    if (n_cells_oi_in_matrix < 5000) {
        genes = exprs_mat %>% apply(1, function(x) {
            sum(x > 0)/n_cells_oi_in_matrix
        }) %>% .[. >= pct] %>% names()
    }
    else {
        splits = split(1:nrow(exprs_mat), ceiling(seq_along(1:nrow(exprs_mat))/100))
        genes = splits %>% lapply(function(genes_indices, exprs,
                                           pct, n_cells_oi_in_matrix) {
            begin_i = genes_indices[1] 
            end_i = genes_indices[length(genes_indices)]
            exprs = exprs[begin_i:end_i, , drop = FALSE]
            genes = exprs %>% apply(1, function(x) {
                sum(x > 0)/n_cells_oi_in_matrix
            }) %>% .[. >= pct] %>% names()
        }, exprs_mat, pct, n_cells_oi_in_matrix) %>% unlist() %>%
            unname()
    }
    return(genes)
}

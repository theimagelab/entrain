
# MIT License

#Copyright (c) 2019 Cole Trapnell

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
    
#    The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

#' Choose cells along the path of a principal graph. Modified from Monocle3 function `choose_graph_segments` to allow for selection of multiple branches and zooming in.
#'
#' @param cds CDS object to be subsetted.
#' @param reduction_method The reduction method to plot while choosing cells.
#'   Currently only "UMAP" is supported.
#' @param starting_pr_node NULL, or a string with the name of the starting
#'   principal node to be used. You can see the principal nodes in your dataset
#'   by using plot_cells with label_principal_points = TRUE.
#' @param ending_pr_nodes NULL, or one or more strings with the name(s) of the
#'   ending principal node(s) to be used. You can see the principal nodes in
#'   your dataset by using plot_cells with label_principal_points = TRUE.
#' @param return_list Logical, return a list of cells instead of a subsetted
#'   CDS object.
#' @param clear_cds Logical, clear CDS slots before returning.
#'   After clearing the cds, re-run processing from preprocess_cds(), ...
#'   Default is TRUE.
#'
#' @return A subset CDS object. If return_list = FALSE, a list of cell and
#'   graph node names.
#' @export
#'
choose_path <- function(cds,
                                  reduction_method = "UMAP",
                                  starting_pr_node = NULL,
                                  ending_pr_nodes = NULL,
                                  return_list = FALSE,
                                  clear_cds = TRUE) {
    
    assertthat::assert_that(methods::is(cds, "cell_data_set"))
    assertthat::assert_that(assertthat::are_equal("UMAP", reduction_method),
                            msg = paste("Currently only 'UMAP' is accepted as a",
                                        "reduction_method."))
    assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                            msg = paste0("No dimensionality reduction for ",
                                         reduction_method, " calculated. ",
                                         "Please run reduce_dimension with ",
                                         "reduction_method = ", reduction_method,
                                         ", cluster_cells and learn_graph ",
                                         "before running choose_graph_segments."))
    assertthat::assert_that(!is.null(cds@clusters[[reduction_method]]),
                            msg = paste("No cell clusters for",
                                        reduction_method, "calculated.",
                                        "Please run cluster_cells with",
                                        "reduction_method =", reduction_method,
                                        "and run learn_graph before running",
                                        "choose_graph_segments."))
    assertthat::assert_that(!is.null(monocle3::principal_graph(cds)[[reduction_method]]),
                            msg = paste("No principal graph for",
                                        reduction_method, "calculated.",
                                        "Please run learn_graph with",
                                        "reduction_method =", reduction_method,
                                        "before running choose_graph_segments."))
    assertthat::assert_that(is.logical(return_list))
    if (all(c(is.null(starting_pr_node),
              is.null(ending_pr_nodes)))) {
        interactive <- TRUE
        assertthat::assert_that(interactive(),
                                msg = paste("Interactive mode not working on",
                                            "this system (known issue with",
                                            "iPython notebooks and on remote",
                                            "servers). Please provide starting",
                                            "and ending principal nodes - see",
                                            "documentation."))
    } else {
        interactive <- FALSE
        assertthat::assert_that(!is.null(starting_pr_node),
                                msg = paste("If not using interactive mode, you",
                                            "must provide a starting_pr_node"))
        assertthat::assert_that(length(starting_pr_node) <= 1,
                                msg = paste("choose_graph_segments only supports",
                                            "1 starting_pr_node. You can pass",
                                            "multiple ending nodes."))
        
        assertthat::assert_that(!is.null(ending_pr_nodes),
                                msg = paste("If not using interactive mode, you",
                                            "must provide ending_pr_nodes."))
    }
    
    
    if (!interactive) {
        sel <- get_principal_path(cds, reduction_method,
                                  starting_cell = starting_pr_node,
                                  end_cells = ending_pr_nodes)
    } else {
        
        
        dp_mst <- cds@principal_graph[[reduction_method]]
        
        princ_points <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
            as.data.frame() %>%
            dplyr::select_(x = 1, y = 2) %>%
            dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))
        row.names(princ_points) <- princ_points$sample_name
        
        edge_df <- dp_mst %>%
            igraph::as_data_frame() %>%
            dplyr::select_(source = "from", target = "to") %>%
            dplyr::left_join(princ_points %>%
                                 dplyr::select_(source="sample_name",
                                                source_prin_graph_dim_1="x",
                                                source_prin_graph_dim_2="y"),
                             by = "source") %>%
            dplyr::left_join(princ_points %>%
                                 dplyr::select_(target="sample_name",
                                                target_prin_graph_dim_1="x",
                                                target_prin_graph_dim_2="y"),
                             by = "target")
        
        data_df <- data.frame(SingleCellExperiment::reducedDims(cds)[[reduction_method]])
        
        colnames(data_df) <- c("data_dim_1", "data_dim_2")
        data_df$sample_name <- row.names(data_df)
        
        data_df <- as.data.frame(cbind(data_df, SingleCellExperiment::colData(cds)))
        
        data_df$cell_color = tryCatch({
            partitions(cds, reduction_method = reduction_method)[data_df$sample_name]},
            error = function(e) {NULL})
        data_df$chosen_cells <- FALSE
        
        ui <- shiny::fluidPage(
            shiny::titlePanel("Choose trajectory branches to analyze."),
            
            # Sidebar layout with input and output definitions ----
            shiny::sidebarLayout(
                
                # Sidebar panel for inputs ----
                shiny::sidebarPanel(
                    shiny::actionButton("choose_start", "Choose starting node"),
                    shiny::actionButton("connect_end", "Connect ending node"),
                    shiny::br(),
                    shiny::actionButton("add_path_to_selection", "Add path"),
                    shiny::br(),
                    # Visualization controls
                    shiny::actionButton("zoom", "Zoom to selection"),
                    shiny::actionButton("reset_zoom", "Reset Zoom"),
                    shiny::br(),
                    # clear button
                    shiny::actionButton("reset", "Clear"),
                    shiny::br(),
                    # done button
                    shiny::actionButton("done", "Done"),
                    shiny::br(),

                    shiny::h3("Instructions:"),
                    shiny::tags$ol(
                        shiny::tags$li("Highlight starting principal_graph node."),
                        shiny::tags$li("Click 'Choose starting node' to highlight."),
                        shiny::tags$li("Highlight ending principal_graph node(s)."),
                        shiny::tags$li(paste("Click 'Connect ending node' to highlight connecting nodes",
                                       "and cells")),
                        shiny::tags$li("Click 'Add path' to add the highlighted cells to the selection"),
                        shiny::tags$li("To select >1 paths to analyze, repeat from beginning"),
                        shiny::tags$li("Click 'Done' to return the chosen selected path/s")
                    ),
                    shiny::h4("Details:"),
                    shiny::tags$ul(
                        shiny::tags$li("To start over, click 'Clear'"),
                        shiny::tags$li(paste("You can choose only 1 starting node.",
                                             "Multiple ending nodes can be selected by", 
                                             "clicking 'Connect ending node' multiple times"))
                    )
                ),
                
                # Main panel for displaying outputs ----
                shiny::mainPanel(
                    shiny::plotOutput("plot1", height="auto",
                                      click = "plot1_click",
                                      brush = shiny::brushOpts(id = "plot1_brush", resetOnNew=TRUE))
                )
            )
        )
        
        server <- function(input, output, session) {
            
            vals <- shiny::reactiveValues(
                start = rep(FALSE, nrow(princ_points)),
                end =  rep(FALSE, nrow(princ_points)),
                chosen = rep(FALSE, nrow(princ_points)),
                chosen_cells = rep(FALSE, nrow(data_df)),
                all_selected_cells = rep(FALSE, nrow(data_df)),
                selected_paths = list(),
                x_range = NULL, y_range = NULL
            )
            
            output$plot1 <- shiny::renderPlot({
                # Plot the selected, current, and unselected points as separate data sets
                princ_points$start <- vals$start
                princ_points$end <- vals$end
                princ_points$chosen <- "Unchosen"
                princ_points$chosen[vals$start] <- "Start"
                princ_points$chosen[vals$end] <- "End"
                princ_points$chosen[vals$chosen] <- "Chosen"
                
                data_df$chosen_cells <- "gray"
                data_df$chosen_cells[vals$all_selected_cells] <- "lightslateblue"
                data_df$chosen_cells[vals$chosen_cells] <- "purple"
                
                suppressMessages(plot_principal_graph(cds, data_df, princ_points,
                                                      label_branch_points = FALSE,
                                                      label_leaves = FALSE,
                                                      label_roots = FALSE) + 
                                     ggplot2::coord_cartesian(xlim = vals$x_range, ylim = vals$y_range, expand = FALSE))
            }, height = function() {
                session$clientData$output_plot1_width
            })
            
            # Toggle points that are brushed, when button is clicked
            shiny::observeEvent(input$choose_start, {
                vals$start = rep(FALSE, nrow(princ_points))
                res <- shiny::brushedPoints(princ_points, xvar = "x", yvar = "y",
                                            input$plot1_brush, allRows = TRUE)
                vals$start <- res$selected_
            })
            
            # Toggle points that are brushed, when button is clicked
            shiny::observeEvent(input$connect_end, {
                vals$end =  rep(FALSE, nrow(princ_points))
                res <- shiny::brushedPoints(princ_points, xvar = "x", yvar = "y",
                                            input$plot1_brush, allRows = TRUE)
                vals$end <- vals$end | res$selected_
                
                chosen <- tryCatch(
                    get_principal_path(cds, reduction_method,
                                       starting_cell = row.names(princ_points)[vals$start],
                                       end_cells = row.names(princ_points)[vals$end]),
                    error = function(e) print(e))
                vals$chosen <- vals$chosen | row.names(princ_points) %in% chosen$nodes
                vals$chosen_cells <- vals$chosen_cells | row.names(pData(cds)) %in%
                    chosen$cells
            })
            
            # Reset all points
            shiny::observeEvent(input$reset, {
                vals$start = rep(FALSE, nrow(princ_points))
                vals$end =  rep(FALSE, nrow(princ_points))
                vals$chosen =  rep(FALSE, nrow(princ_points))
                vals$chosen_cells = rep(FALSE, nrow(data_df))
                vals$all_selected_cells = rep(FALSE, nrow(data_df))
            })
            
            shiny::observeEvent(input$done, {
                shiny::stopApp(vals$selected_paths)
            })

            shiny::observeEvent(input$add_path_to_selection, {
                current_path <- list(nodes = row.names(princ_points)[vals$chosen],
                                    cells = row.names(data_df)[vals$chosen_cells])
                
                current_path_name <- paste(row.names(princ_points)[vals$start][1],
                                           row.names(princ_points)[vals$end][1],
                                           sep="_")
                vals$selected_paths[[ current_path_name ]] <- current_path
                vals$all_selected_cells[vals$chosen_cells] <- TRUE
                
                vals$start = rep(FALSE, nrow(princ_points))
                vals$end =  rep(FALSE, nrow(princ_points))
                vals$chosen =  rep(FALSE, nrow(princ_points))
                vals$chosen_cells = rep(FALSE, nrow(data_df))
                })
            
            shiny::observeEvent(input$zoom, {
                brush <- input$plot1_brush
                if (!is.null(brush)) {
                    vals$x_range <- c(brush$xmin, brush$xmax)
                    vals$y_range <- c(brush$ymin, brush$ymax)
                    
                } else {
                    vals$x_range <- NULL
                    vals$y_range <- NULL
                }
            })
            
            shiny::observeEvent(input$reset_zoom, {
                    vals$x_range <- NULL
                    vals$y_range <- NULL
            })
        }
        sel <- shiny::runApp(shiny::shinyApp(ui, server))
    }
    if(return_list) {
        return(sel)
    } else {
        cds<-cds[,sel$cells]
        monocle3::principal_graph(cds)[[reduction_method]] <- igraph::induced_subgraph(monocle3::principal_graph(cds)[[reduction_method]], sel[[1]])
        if( clear_cds )
            cds<-clear_cds_slots(cds)
        return(cds)
    }
}

plot_principal_graph <- function(cds,
                                 data_df,
                                 princ_points,
                                 reduction_method = "UMAP",
                                 trajectory_graph_color="black",
                                 trajectory_graph_segment_size=0.75,
                                 label_groups_by_cluster=TRUE,
                                 group_label_size=2,
                                 labels_per_group=1,
                                 label_branch_points=TRUE,
                                 label_roots=TRUE,
                                 label_leaves=TRUE,
                                 graph_label_size=2,
                                 cell_size=0.35,
                                 alpha = 1,
                                 min_expr=0.1,
                                 rasterize=FALSE) {
    
    gene_short_name <- NA
    sample_name <- NA
    #sample_state <- colData(cds)$State
    data_dim_1 <- NA
    data_dim_2 <- NA
    plotting_func <- ggplot2::geom_point
    
    
    ## Graph info
    
    ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
        as.data.frame() %>%
        dplyr::select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
        dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))
    
    dp_mst <- cds@principal_graph[[reduction_method]]
    
    edge_df <- dp_mst %>%
        igraph::as_data_frame() %>%
        dplyr::select_(source = "from", target = "to") %>%
        dplyr::left_join(ica_space_df %>%
                             dplyr::select_(
                                 source="sample_name",
                                 source_prin_graph_dim_1="prin_graph_dim_1",
                                 source_prin_graph_dim_2="prin_graph_dim_2"),
                         by = "source") %>%
        dplyr::left_join(ica_space_df %>%
                             dplyr::select_(
                                 target="sample_name",
                                 target_prin_graph_dim_1="prin_graph_dim_1",
                                 target_prin_graph_dim_2="prin_graph_dim_2"),
                         by = "target")
    
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))
    
    g <- g + geom_point(color=data_df$chosen_cells, size=I(cell_size),
                        na.rm = TRUE, alpha = I(alpha))
    
    message(paste("cluster_cells() has not been called yet, can't color cells",
                  "by cluster"))
    
    g <- g + geom_point(aes(x = x, y = y, color = chosen), data=princ_points) +
        scale_color_manual(values = c("Start" = "green", "End" = "blue",
                                      "Unchosen" = "black", "Chosen" = "purple"))
    g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1",
                                     y="source_prin_graph_dim_2",
                                     xend="target_prin_graph_dim_1",
                                     yend="target_prin_graph_dim_2"),
                          size=trajectory_graph_segment_size,
                          linetype="solid",
                          na.rm=TRUE,
                          data=edge_df)
    
    
    if (label_branch_points){
        mst_branch_nodes <- branch_nodes(cds)
        branch_point_df <- ica_space_df %>%
            dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
            dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
        
        g <- g +
            geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                       shape = 21, stroke=I(trajectory_graph_segment_size),
                       color="white",
                       fill="black",
                       size=I(graph_label_size * 1.5),
                       na.rm=TRUE, branch_point_df) +
            
            geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                                 label="branch_point_idx"),
                      size=I(graph_label_size), color="white", na.rm=TRUE,
                      branch_point_df)
    }
    
    if (label_leaves){
        mst_leaf_nodes <- leaf_nodes(cds)
        leaf_df <- ica_space_df %>%
            dplyr::slice(match(names(mst_leaf_nodes), sample_name)) %>%
            dplyr::mutate(leaf_idx = seq_len(dplyr::n()))
        
        g <- g +
            geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                       shape = 21, stroke=I(trajectory_graph_segment_size),
                       color="black",
                       fill="lightgray",
                       size=I(graph_label_size * 1.5),
                       na.rm=TRUE,
                       leaf_df) +
            geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                                 label="leaf_idx"),
                      size=I(graph_label_size), color="black", na.rm=TRUE, leaf_df)
    }
    
    if (label_roots){
        mst_root_nodes <- root_nodes(cds)
        root_df <- ica_space_df %>%
            dplyr::slice(match(names(mst_root_nodes), sample_name)) %>%
            dplyr::mutate(root_idx = seq_len(dplyr::n()))
        
        g <- g +
            geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                       shape = 21, stroke=I(trajectory_graph_segment_size),
                       color="black",
                       fill="white",
                       size=I(graph_label_size * 1.5),
                       na.rm=TRUE,
                       root_df) +
            geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                                 label="root_idx"),
                      size=I(graph_label_size), color="black", na.rm=TRUE, root_df)
    }
    
    g <- g +
        monocle3:::monocle_theme_opts() +
        xlab(paste(reduction_method, 1)) +
        ylab(paste(reduction_method, 2)) +
        theme(legend.key = element_blank()) +
        theme(panel.background = element_rect(fill='white'))
    g
}

get_principal_path <- function(cds, reduction_method,
                               starting_cell, end_cells) {
    
    subset_principal_nodes <- c()
    dp_mst <- principal_graph(cds)[[reduction_method]]
    
    for(end_cell in end_cells) {
        traverse_res <- traverse_graph(dp_mst, starting_cell, end_cell)
        path_cells <- names(traverse_res$shortest_path[[1]])
        if(length(path_cells) == 0) {
            stop(paste0("Starting and ending nodes are not connected"))
        }
        
        subset_principal_nodes <- c(subset_principal_nodes, path_cells)
    }
    subset_principal_nodes <- unique(subset_principal_nodes)
    corresponding_cells <- which(paste0("Y_", cds@principal_graph_aux[[
        reduction_method]]$pr_graph_cell_proj_closest_vertex) %in%
            subset_principal_nodes)
    subset_cells <- row.names(cds@principal_graph_aux[[
        reduction_method]]$pr_graph_cell_proj_closest_vertex)[corresponding_cells]
    
    return(list(nodes = subset_principal_nodes, cells = subset_cells))
}


traverse_graph <- function(g, starting_cell, end_cells){
    distance <- igraph::shortest.paths(g, v=starting_cell, to=end_cells)
    branchPoints <- which(igraph::degree(g) == 3)
    path <- igraph::shortest_paths(g, from = starting_cell, end_cells)
    
    return(list(shortest_path = path$vpath, distance = distance,
                branch_points = intersect(branchPoints, unlist(path$vpath))))
}
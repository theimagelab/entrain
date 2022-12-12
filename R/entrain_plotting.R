
#' @title Plots a heatmap of ligand-target regulatory potentials extracted from NicheNet.
#' @description Visualize the ligand-target regulatory relationships that are driving the analyzed trajectory paths.
#' @param obj A seurat object that has been analyzed with get_path_ligands or get_traj_ligands_monocle.
#' @param path A path name. Must be from in `names(obj@misc$entrain$paths)`
#' @param n_top_targets Number of top importance target genes to visualize. Default 10.
#' @param n_top_ligands Number of top ligands to visualize. Default 10.
#' @param ligand_target_matrix NicheNet ligand-target data file.
#' @param metric One of "Covariances" or "Correlations". Default "Covariances". Determines whether to use pseudotime covariance or pseudotime correlation in the calculation of TRAINing genes.


#' @return A plotly plot object.
#' @export
plot_ligand_targets <- function(obj,
                                path,
                                n_top_targets = 10,
                                n_top_ligands = 10,
                                ligand_target_matrix,
                                metric = "Covariances") {

    if (is.null(path) | is.null(obj@misc$entrain$paths[[path]])) {
        stop(paste( "Please select a path from available paths: \n", paste0(names(obj@misc$entrain$paths), collapse=", \n")))
    }

    pathdata<-obj@misc$entrain$paths[[path]]
    path_ligands<- pathdata$model$importance
    top_ligands <- sort(path_ligands[,"IncNodePurity"], decreasing=TRUE) %>% head(n_top_ligands)
    pct_pseudotime_explained <- pathdata$model$rsq %>% tail(1)
    pseudotime_genes <- pathdata$pseudotime_associated_genes
    pseudotime_genes_vector <- pseudotime_genes[[metric]]
    names(pseudotime_genes_vector) <- row.names(pseudotime_genes)
    top_pseudotime_genes <- pseudotime_genes_vector %>% sort(., decreasing=TRUE) %>% head(n_top_targets)
    ligand_target_matrix_mm %>% .[row.names(.) %in% names(top_pseudotime_genes),
                                  colnames(.) %in% names(top_ligands)] %>% as.matrix() %>% t() -> top_ligands_potentials


    df<-reshape2::melt(top_ligands_potentials, id.vars = c("Ligand", "TargetGene")) %>% dplyr::rename("Ligand" = "Var1", "TargetGene" = "Var2")
    hm<-ggplot(df, aes(x=TargetGene, y=Ligand, fill=value)) +
        geom_tile() +
        theme_minimal_big_font()
    return(hm)
}

#' @title Scatter plot of top trajectory covarying genes.
#' @description Recovers the expression of the top genes that covary with pseudotime and plots it as a scatter plot.
#' @param obj A seurat object that has been analyzed with get_path_ligands or get_traj_ligands_monocle.
#' @param path Name of list element in obj@misc$entrain$paths corresponding to the trajectory path that you wish to visualize.
#' @param n_top_genes Number of top importance target genes to visualize. Default 10.
#' @param color_key A column name in obj@meta.data corresponding to a label you wish to color by e.g., a cluster label.
#' @param nrow Number of rows for the facet plot grid.
#' @param ncol Number of columns for the facet plot grid.
#' @param metric One of "Covariances" or "Correlations". Default "Covariances". Determines whether to use pseudotime covariance or pseudotime correlation in the calculation of TRAINing genes.

#' @import ggplot2
#' @return A ggplot2 plot object.
#' @export
plot_covarying_genes_scatter <- function(obj,
                                         path = paths[1],
                                         n_top_genes = 10,
                                         color_key = NULL,
                                         nrow=5,
                                         ncol=2,
                                         metric = "Covariances") {
    if (is.null(path) | is.null(obj@misc$entrain$paths[[path]])) {
        stop(paste( "Please select a path from available paths: \n", paste0(names(obj@misc$entrain$paths), collapse=", \n")))
    }

    if (is.null(color_key)) {
        message("No color key supplied. Using default SeuratObject Idents,")
        color_key <- levels(x = obj)
    }

    pathdata<-obj@misc$entrain$paths[[path]]

    pseudotime_genes <- pathdata$pseudotime_associated_genes
    pseudotime_genes_vector <- pseudotime_genes[[metric]]
    names(pseudotime_genes_vector) <- row.names(pseudotime_genes)
    top_pseudotime_genes <- pseudotime_genes_vector %>% sort(., decreasing=TRUE)
    genes_to_plot <- top_pseudotime_genes %>% head(n_top_genes) %>% names()

    # Subset out expression data of these genes to a separate dataframe.
    expr <- Seurat::GetAssayData(subset(obj, cells = pathdata$cell_pseudotimes$path_cell_names))
    expr_pseu <- expr[names(top_pseudotime_genes),pathdata$cell_pseudotimes$path_cell_names] %>% as.matrix() %>% t() %>% as.data.frame()
    expr_pseu$pseudotime <- pathdata$cell_pseudotimes$path_pseudotime

    if (color_key != "pseudotime") {
        expr_pseu <- merge(expr_pseu, obj@meta.data[color_key], by=0)
        df<-expr_pseu[c("pseudotime", color_key, genes_to_plot)]
        df <- reshape2::melt(df ,  id.vars = c('pseudotime',color_key), variable.name = 'Target')

    } else {
        df<-expr_pseu[c("pseudotime", genes_to_plot)]
        df <- reshape2::melt(df , id.vars = c('pseudotime'), variable.name = 'Target')

    }


    # Plot
    covarying_genes_plot<-ggplot2::ggplot(df, aes_string('pseudotime','value', color=color_key)) +
        ggplot2::geom_point(size=0.5, position = position_jitter(w = 0.7, h = 0)) +
        ggplot2::facet_wrap(Target ~ ., nrow=nrow, ncol=ncol) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank(),  strip.text = element_text(size = 14))
}

#' @title Heatmap plot of top trajectory gene covariances.
#' @description Recovers the covariance of the top genes that covary with pseudotime and plots it as a 1d-heatmap.
#' @param obj A seurat object that has been analyzed with get_path_ligands or get_traj_ligands_monocle.
#' @param path Name of list element in obj@misc$entrain$paths corresponding to the trajectory path that you wish to visualize.
#' @param n_top_genes Number of top covarying genes to visualize. Default 10.
#' @param metric One of "Covariances" or "Correlations". Metric used during fitting. Default "Covariance"
#'
#' @return A ggplot2 plot object.
#' @export
plot_covarying_genes_heatmap <- function(obj,
                                         path = paths[1],
                                         n_top_genes = 10,
                                         metric = "Covariances") {
    if (is.null(path) | is.null(obj@misc$entrain$paths[[path]])) {
        stop(paste( "Please select a path from available paths: \n", paste0(names(obj@misc$entrain$paths), collapse=", \n")))
    }

    pathdata<-obj@misc$entrain$paths[[path]]

    pseudotime_genes <- pathdata$pseudotime_associated_genes
    pseudotime_genes_values <- pseudotime_genes[[metric]]
    names(pseudotime_genes_values) <- row.names(pseudotime_genes)
    pseudotime_genes_values <- pseudotime_genes_values %>% sort(., decreasing=TRUE)

    df <- pseudotime_genes_values %>% head(n_top_genes) %>% as.data.frame()
    colnames(df) <- metric
    df$gene <- rownames(df)

    plot <- ggplot(df) +
        geom_tile(aes(x = 1,
                      y = factor(gene, levels = rev(rownames(df))),
                      fill = Covariance)) +
        theme_minimal_big_font() + theme(axis.text.x = element_blank(),
                                         axis.ticks.x = element_blank(),
                                         axis.title.y = element_blank(),
                                         axis.title.x = element_blank())
    return(plot)
}

#' @title Draw circos plot of Entrain ligands.
#' @description 'Given a seurat object containing Entrain results, draw the ligands driving the observed velocities or trajectories..
#' @param obj A seurat object containing the Entrain result from entrain_velocity or get_traj_ligands()
#' @param path_name Optional. String indicating the name of the trajectory path. Use this if you want to plot ligands associated with pseudotime branches.
#' @param vcluster Optional. String indicating the name of the velocity cluster. Use this if you wish to plot ligands associated with velocity clusters.
#' @param n_top_ligands Number of top ligands to plot
#' @param lr_network NicheNet ligand-receptor pairs data file.
#' @param weighted_networks NicheNet weighted networks.
#' @param palette RColorBrewer color palette for plotting sender cell clors.

#' @import dplyr
#' @return NULL. Plots the circlize object to graphics output.

#' @export
draw_entrain_circos <- function(obj,
                                path_name = NULL,
                                vcluster = NULL,
                                n_top_ligands = 10,
                                lr_network,
                                weighted_networks,
                                palette = "Spectral") {
    if (!is.null(path_name) & is.null(vcluster)){
        ligand_weights <- obj@misc$entrain$paths[[path_name]]$path_ligands_importances
        top_ligands <- ligand_weights %>% names() %>% head(n_top_ligands)
        expressed_receptors <- obj@misc$entrain$paths[[path_name]]$active_network$to %>% unname() %>% unique()

        top_ligands_weights <- ligand_weights[top_ligands] %>% as.data.frame()
        colnames(top_ligands_weights) <- "entrain_weight"
        top_ligands_weights$ligand <- rownames(top_ligands_weights)

    } else if (!is.null(vcluster) & is.null(path_name)) {
        ligand_weights <- obj@misc$entrain$velocity_result$vclust_ligand_importances[c("ligand", vcluster)] %>%
            arrange(., desc(get(vcluster)))
        top_ligands<-ligand_weights$ligand %>% head(n_top_ligands)
        expressed_receptors <- obj@misc$entrain$velocity_model[[vcluster]]$active_network$to %>% unname() %>% unique()

        top_ligands_weights <- ligand_weights %>% filter(ligand %in% top_ligands)
        top_ligands_weights[is.na(top_ligands_weights)] <- 0
        colnames(top_ligands_weights) <- c("ligand", "entrain_weight")
    } else {
        stop("Please supply either a vcluster_name or a path_name, not both.")
    }
    # get the ligand-receptor network of the top-ranked ligands
    lr_network_top = lr_network %>% dplyr::filter(from %in% top_ligands & to %in% expressed_receptors) %>% dplyr::distinct(from,to)
    best_upstream_receptors = lr_network_top %>% dplyr::pull(to) %>% unique()

    weighted_networks_lr = weighted_networks$lr_sig %>% dplyr::inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

    # ligand to cell type associations, for coloring ligands in the circos plot.
    avg_expression_ligands = Seurat::AverageExpression(obj, features = top_ligands)

    sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
        ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
    }) %>% t()
    sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})

    all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
    unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
    general_ligands = top_ligands %>% setdiff(unique_ligands)

    ligand_celltype_assignments <- list()
    ligand_type_indication_df <- data.frame(matrix(nrow=0, ncol=2))
    colnames(ligand_type_indication_df) <- c("ligand_type", "ligand")

    # Dataframe associating ligands to cell types
    for (celltype in names(sender_ligand_assignment)) {
        ligand <- sender_ligand_assignment[[celltype]] %>% names %>% setdiff(general_ligands)
        if (length(ligand) > 0) {
            ligand_celltype_assignments[[celltype]] <- ligand
            df<-data.frame(ligand)
            df$ligand_type <- celltype
            ligand_type_indication_df <- rbind(ligand_type_indication_df, df)
        }
    }

    ligand_type_indication_df <- rbind(ligand_type_indication_df,
                                       data.frame(
                                           ligand_type = c(rep("General/other", times = general_ligands %>% length())),
                                           ligand = general_ligands
                                       )
    )

    ligand_type_indication_df <- tibble::tibble(ligand_type_indication_df)

    lr_network_top_df = weighted_networks_lr %>% dplyr::filter(from %in% top_ligands & to %in% best_upstream_receptors) %>% dplyr::rename(ligand = from, receptor = to)
    lr_network_top_df = lr_network_top_df %>% dplyr::mutate(receptor_type = "receptor") %>% dplyr::inner_join(ligand_type_indication_df, by = "ligand")

    # Instead of nichenet weights, we use entrain importance scores scaled by nichenet scores.
    colnames(lr_network_top_df) <- c("ligand", "receptor", "nichenet_weight", "receptor_type", "ligand_type")

    lr_network_top_df<-dplyr::left_join(lr_network_top_df, top_ligands_weights, by="ligand")
    lr_network_top_df <- dplyr::group_by(lr_network_top_df, ligand) %>%
        mutate(weight = entrain_weight*nichenet_weight/max(nichenet_weight)) %>%
        ungroup()

    col.dist <- function(inp, comp) sum( abs(inp - col2rgb(comp) ) )

    ligand_types <- unique(ligand_type_indication_df$ligand_type) %>% setdiff("General/other")
    if (length(palette) > 1 & all(ligand_types %in% names(palette))) {
        colorscale <- palette
        colorscale <- colorscale[names(colorscale) %in% ligand_types]
    } else if (length(palette) == 1) {
        colorscale <- RColorBrewer::brewer.pal(length(ligand_types), palette)
    } else {
        message("Please supply to argument 'palette' either 1: RColorBrewer palette name, or 2: a named vector consisting of hex code colours, with names corresponding to a value in the active identity class of the seurat object")
        message("Valid names for palette: ")
        message(paste(ligand_types, sep=", "))
        stop()
    }
    grid_col_ligand<- colors()[ apply(col2rgb(colorscale), 2,
                                      function(z) which.min( sapply(colors(),
                                                                    function(x) col.dist(inp=z, comp=x) ) ) ) ]

    names(grid_col_ligand) <- ligand_types

    grid_col_ligand <- c(grid_col_ligand, "General/other" = "gray50")
    grid_col_receptor =c(
        "receptor" = "darkred")

    grid_col_tbl_ligand = tibble::tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
    grid_col_tbl_receptor = tibble::tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

    circos_links = lr_network_top_df %>% dplyr::mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as receptor!
    circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand, by="ligand_type") %>% dplyr::inner_join(grid_col_tbl_receptor, by="receptor_type")
    links_circle = circos_links %>% dplyr::select(ligand,receptor, weight)

    ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
    grid_ligand_color = ligand_color$color_ligand_type %>% purrr::set_names(ligand_color$ligand)
    receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
    grid_receptor_color = receptor_color$color_receptor_type %>% purrr::set_names(receptor_color$receptor)

    grid_col = c(grid_ligand_color,grid_receptor_color)

    # give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
    transparency = circos_links %>%
        dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>%
        dplyr::mutate(transparency = 1-weight) %>%
        .$transparency

    receptor_order = circos_links$receptor %>% unique()
    ligand_order = c(top_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
    order = c(ligand_order,receptor_order)

    width_same_cell_same_ligand_type = 0.5
    width_different_cell = 6
    width_ligand_receptor = 15
    width_same_cell_same_receptor_type = 0.5

    gaps <- c()
    for (celltype in unique(circos_links$ligand_type)) {
        gaps<-c(
            gaps,
            rep(width_same_cell_same_ligand_type,
                times = (circos_links %>% dplyr::filter(ligand_type == celltype) %>% dplyr::distinct(ligand) %>% nrow() -1)),
            width_different_cell
        )
    }
    gaps<-c(
        gaps,
        rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receptor_type == "receptor") %>% dplyr::distinct(receptor) %>% nrow() -1)),
        width_ligand_receptor
    )

    cutoff_include_all_ligands = 0
    circlize::circos.clear()
    circlize::circos.par(gap.degree = gaps)

    circlize::chordDiagram(links_circle,
                           directional = 1,
                           order=order,link.sort = TRUE,
                           link.decreasing = FALSE, grid.col = grid_col,
                           transparency = transparency,
                           diffHeight = 0.005,
                           direction.type = c("diffHeight", "arrows"),
                           link.arr.type = "big.arrow",
                           link.visible = links_circle$weight >= cutoff_include_all_ligands,
                           annotationTrack = "grid",
                           preAllocateTracks = list(track.height = 0.075))
    circlize::circos.track(
        track.index = 1,
        panel.fun = function(x, y) {
            circlize::circos.text(
                circlize::CELL_META$xcenter,
                circlize::CELL_META$ylim[1],
                circlize::CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.55), cex = 0.8
            )
        },
        bg.border = NA
    )
}

#' @title Plots ligand influences and variance explained scores on a Monocle trajectory plot.
#' @description Requires get_path_ligands or get_traj_ligands_monocle to be run on dataset beforehand. Specifically for plotting Entrain results analyzed on a Monocle trajectory.
#' @param obj A seurat object that has been analyzed with get_path_ligands or get_traj_ligands_monocle.
#' @param n_top_ligands Number of top ligands to write on each plot label.
#' @param trajectory_graph_segment_size Parameters to be passed to Monocle plot_cells function.
#' @param trajectory_graph_color Parameters to be passed to Monocle plot_cells function.
#' @param cell_size Parameters to be passed to Monocle plot_cells function.
#' @param group_label_size Parameters to be passed to Monocle plot_cells function.
#' @param cds Optional cell_data_set object to use instead of the default cds at `obj@misc$entrain$monocle_cds`
#' @param color_cells_by Argument for `monocle3::plot_cells()`. A column in `colData(obj@misc$entrain$monocle_cds)` or `colData(cds)` denoting how to color the cells. Default `pseudotime`
#' @param label_groups_by_cluster Argument for `monocle3::plot_cells()`. Default FALSE.
#' @param label_font_size Argument for `monocle3::plot_cells()`. Default 4.

#' @return a ggplot object with cells, trajectories, and ligand labels plotted on a UMAP.
#' @export
plot_ligand_trajectories <- function(obj, color_cells_by = 'pseudotime', n_top_ligands = 5,
                                     trajectory_graph_color = "grey70",
                                     group_label_size = 5,
                                     trajectory_graph_segment_size = 1,
                                     cell_size=0.6,
                                     label_groups_by_cluster=FALSE,
                                     label_font_size=4,
                                     cds=NULL) {
    paths <- names(obj@misc$entrain$paths)

    if (is.null(cds)) {
        cds<-obj@misc$entrain$monocle_cds
    }

    labels_df <- get_ligand_plot_labels(obj, n_top_ligands, paths, reduction = NULL, label_position=NULL)

    g <- monocle3::plot_cells(cds, color_cells_by=color_cells_by,
                              show_trajectory_graph=TRUE,
                              group_label_size = group_label_size,
                              label_principal_points = FALSE,
                              cell_size=cell_size,
                              label_groups_by_cluster = label_groups_by_cluster,
                              trajectory_graph_segment_size = trajectory_graph_segment_size,
                              trajectory_graph_color = trajectory_graph_color,
                              label_branch_points = FALSE,
                              label_roots = FALSE,
                              label_leaves = FALSE) +
        ggrepel::geom_label_repel(data=labels_df,
                                  aes(x=as.numeric(labels_df$x), y=as.numeric(labels_df$y),
                                      label=labels_df$label,
                                  ),
                                  fill = alpha(c("white"),0.7),
                                  size = label_font_size,
                                  inherit.aes = FALSE)
    return(g)
}



#' @title Plots ligand cellwise inlfuences on a UMAP.
#' @description Requires cellwise_influences to be run on obj beforehand.
#' @param obj A seurat object that has been analyzed with get_path_ligands or get_traj_ligands_monocle.
#' @param reduction Name of reduction to be used for plotting. Defaults to "MonocleUMAP"; the UMAP coordinates generated by the Monocle workflow.
#' @param colormap String denoting a brewer.pal compatible color palette.
#' @param custom_colors Optional vector of colors to be passed to ggplot. Will override colormap parameter.
#' @param ligand_labels Boolean. Overlay the ligand label on the plot?
#' @param n_top_ligands Intheger denoting how many of the top ligands to overlay on the label. Default 5.
#' @param label_font_size Font size of label. Default 4.
#' @param label_position Where in the trajectory to place the label. Default `'endpoint'`, (i.e., later in pseudotime)

#' @return A ggplot object with cells plotted on a UMAP with cell color corresponding to the degree of ligand influence for a particular branch.
#' @export

plot_ligand_influences <- function(obj,
                                   colormap="Set1",
                                   reduction="MonocleUMAP_",
                                   custom_colors=NULL,
                                   ligand_labels=TRUE,
                                   n_top_ligands = 5,
                                   label_font_size = 4,
                                   label_position = "endpoint") {

    meta_cols <- colnames(obj@meta.data)
    influences_columns <- meta_cols[grep("Influences_", meta_cols)]

    influence_data<-obj@meta.data[c(influences_columns)]

    reduction_names<-names(obj@reductions[[reduction]])

    if (is.null(custom_colors) == FALSE) {
        discrete_colors<-custom_colors
    } else {
        discrete_colors <- RColorBrewer::brewer.pal(n=length(influences_columns), colormap) %>% shades::saturation(2)
    }
    g<- ggplot()
    color_index = 1
    for (column in influences_columns) {
        cells<-row.names(obj@meta.data[obj@meta.data[[column]] != 0,])
        cells<-obj@reductions[[reduction]]@cell.embeddings[cells,] %>% as.data.frame()
        cells<-merge(cells, influence_data, by.x=0, by.y=0, all.y=FALSE, all.x=TRUE)
        cells<-transform(cells, row.names=Row.names, Row.names=NULL)
        g <- g + geom_point(data=cells,
                            aes(
                                x=.data[[reduction_names[1]]],
                                y=.data[[reduction_names[2]]],
                                color=.data[[column]]),
                            size=0.3
        ) +
            scale_color_gradient(low="gray95", high=discrete_colors[color_index]) +
            ggnewscale::new_scale_color() +
            theme_bw() + theme(panel.border = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black"),
                               axis.text=element_text(size=12))

        color_index = color_index + 1
        }


    # plot non-analyzed cells
    other_cells<-obj@meta.data[c(influences_columns)] %>% dplyr::filter_all(., dplyr::all_vars(. == 0)) %>% row.names()
    other_cells<-obj@reductions[[reduction]]@cell.embeddings[other_cells,] %>% as.data.frame()
    other_cells<-merge(other_cells, influence_data, by.x=0, by.y=0, all.y=FALSE, all.x=TRUE)
    other_cells<-transform(other_cells, row.names=Row.names, Row.names=NULL)
    g <- g + geom_point(data=other_cells,
                        aes(
                            x=.data[[reduction_names[1]]],
                            y=.data[[reduction_names[2]]]),
                        color = "gray60",
                        size=0.3
    )

    if (ligand_labels==TRUE) {
        paths<-names(obj@misc$entrain$paths)

        labels_df <- get_ligand_plot_labels(obj,
                                            n_top_ligands,
                                            paths,
                                            reduction = reduction,
                                            label_position = label_position)
        g <- g +
            ggrepel::geom_label_repel(
                data=labels_df,
                aes(x=as.numeric(labels_df$x), y=as.numeric(labels_df$y),
                    label=labels_df$label,
                ),
                fill = alpha(c("white"),0.7),
                size = label_font_size,
                inherit.aes = FALSE
            )
    }
    return(g)
}

get_ligand_plot_labels <- function(obj,
                                   n_top_ligands,
                                   paths,
                                   reduction=NULL,
                                   label_position=NULL) {

    labels_df<-data.frame(matrix(ncol=5))
    colnames(labels_df) <- c("x", "y", "ligands", "var", 'label')
    for (path in paths) {
        x<-as.numeric(obj@misc$entrain$paths[[path]]$path_midpoint[[2]][1])
        y<-as.numeric(obj@misc$entrain$paths[[path]]$path_midpoint[[2]][2])

        if (!is.null(label_position)) {
            if (label_position == "endpoint") {
                pseu_df <- obj@misc$entrain$paths[[path]]$cell_pseudotimes
                pseu_max <- pseu_df[which.max(pseu_df$path_pseudotime),"path_cell_names"]
                xy <- obj@reductions[[reduction]]@cell.embeddings[pseu_max,]
                x <- xy[1]
                y <- xy[2]
            }
        }

        ligands<-obj@misc$entrain$paths[[path]]$path_ligands_importances %>% head(n_top_ligands) %>% names() %>% paste(., collapse = '\n')
        var<-obj@misc$entrain$paths[[path]]$model$rsq %>% tail(1)
        label<-paste(ligands, "\n", as.character(round(var, digits=3)*100), "% V.E.")
        labels_df[path,] <- c(x, y, ligands, var, label)
    }
    labels_df<-na.omit(labels_df)
    return(labels_df)
}

#' @title Plots a heatmap of ligand-target regulatory potentials extracted from NicheNet.
#' @description 'Visualize the ligand-target regulatory relationships that are driving the analyzed velocities.
#' @param obj A seurat object that has been analyzed with entrain_velocity
#' @param n_top_targets Number of top likelihood target genes to visualize. Default 10.
#' @param n_top_ligands Number of top ligands to visualize. Default 10.
#' @param ligand_target_matrix NicheNet ligand-target data file.
#' @param velocity_cluster Name of the velocity cluster whose regulatory relationships you wish to plot.
#' @param colorscale Name of an RColorBrewer palette. Default `Greens`
#' @param color_low Hex code of lower limit
#' @param color_high Hex code of upper limit.

#' @return A gglot2 plot object.
#' @export
plot_velocity_ligand_targets <- function(obj,
                                         n_top_targets = 10,
                                         velocity_cluster = NULL,
                                         n_top_ligands = 10,
                                         ligand_target_matrix,
                                         color_low = "#f6ffff",
                                         color_high = "#01bfc4",
                                         colorscale = "Greens") {
    res <- obj@misc$entrain$velocity_result
    imps <- res$vclust_ligand_importances
    likelihds <- res$vclust_gene_likelihoods

    if (is.null(velocity_cluster)){
        stop("Please supply the name of a velocity cluster")
    }

    top_imps <- imps[,c("ligand", velocity_cluster), drop=FALSE]
    top_imps <- top_imps[!is.na(top_imps[[velocity_cluster]]),] %>% arrange(desc(get(velocity_cluster))) %>% head(n_top_ligands)

    top_likelihds <- likelihds[,c("gene", velocity_cluster), drop=FALSE]
    top_likelihds <- top_likelihds[!is.na(top_likelihds[[velocity_cluster]]),] %>% arrange(desc(get(velocity_cluster))) %>% head(n_top_targets)

    ligand_target_matrix %>%
        .[row.names(.) %in% top_likelihds$gene,
          colnames(.) %in% top_imps$ligand] %>%
        as.matrix() %>%
        t() %>% reshape2::melt() %>%
        stats::setNames(c("ligand", "gene", "potential")) -> df

    plot <- ggplot2::ggplot(df) +
        geom_tile(aes(x = gene,
                      y = factor(ligand, levels = rev(top_imps$ligand)),
                      fill = potential)) +
        labs(y="Ligand", x = "Gene") + ggtitle(velocity_cluster) +
        theme_minimal_big_font()

    if (!is.null(colorscale)) {
        plot <- plot + scale_fill_distiller(palette = colorscale, direction=1)
    } else if (!is.null(color_high) & !is.null(color_low)) {
        plot <- plot + scale_fill_gradient(low = color_low, high = color_high)
    } else {
        stop("Please supply a colorscale or color_low and color_high")
    }
    return(plot)
}

#' @title Velocity Ligand Importance Heatmap
#' @description Visualize the ligands and their feature importances towards the velocity clusters.
#' @param obj A seurat object that has been analyzed with entrain_velocity
#' @param velocity_clusters Character vector of velocity clusters to include in the plot. Default `NULL` - include all velocity clusters,
#' @param n_top_ligands Number of top ligands to visualize. Default 5.
#' @param ligand_target_matrix NicheNet ligand-target data file.
#' @param color_low Hex code of lower limit
#' @param color_high Hex code of upper limit.
#' @param colorscale Name of a color palette to be passed to the palette argument of scale_fill_distiller. Overrides color_low and color_high.
#' @param rescale_columns Boolean. If TRUE, rescales all ligand importances between 0 and 1. This is recommended because random forest importances are relative.
#' @return A gglot2 plot object.
#' @export
velocity_ligand_importance_heatmap <- function(obj,
                                         velocity_clusters = NULL,
                                         n_top_ligands = 5,
                                         ligand_target_matrix,
                                         color_low = "#f6ffff",
                                         color_high = "#01bfc4",
                                         colorscale = "BuPu",
                                         rescale_columns = TRUE) {
    res <- obj@misc$entrain$velocity_result
    lig_importances <- res$vclust_ligand_importances

    if (!is.null(velocity_clusters)) {
        vcluster_imps <- lig_importances[,c("ligand", velocity_clusters), drop=FALSE]
    } else {
        vcluster_imps <- lig_importances
    }

    # get top n ligands in each vcluster
    ligands<-list()
    x <- vcluster_imps %>% dplyr::select(., -ligand)
    for (col in colnames(x)) {
        idx <- vcluster_imps[, col] %>% order(., decreasing = TRUE) %>% head(n_top_ligands)
        top_ligs <- vcluster_imps[idx, "ligand"]
        ligands <- append(ligands, top_ligs)
    }

    selected_ligands <- ligands %>% unlist() %>% unique()
    lig_heatmap <- vcluster_imps %>% dplyr::filter(., ligand %in% selected_ligands)
    lig_heatmap[is.na(lig_heatmap)] <- 0

    if (rescale_columns == TRUE) {
        lig_heatmap[,-1] <- lig_heatmap %>% dplyr::select(., -ligand) %>% sapply(., scales::rescale, to=c(0,1))
    }

    lig_heatmap <- lig_heatmap %>% reshape2::melt(id.var="ligand")
    x_order <- lig_importances[,2:ncol(lig_importances)] %>% colnames %>% sort()
    plot <- ggplot(lig_heatmap) +
        geom_tile(aes(x=factor(variable, levels=x_order), y=ligand, fill=value)) +
        theme_minimal_big_font() +
        ggtitle("Ligand Importances") +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())
    if (!is.null(colorscale)) {
        plot <- plot + scale_fill_distiller(palette=colorscale, direction=1)
    } else if (!is.null(color_high) & !is.null(color_low)) {
        plot <- plot + scale_fill_gradient(low = color_low, high = color_high)
    } else {
        stop("Please supply a colorscale or color_low and color_high")
    }
    return(plot)
}

#' @title Plot velocity likelihoods as a heatmap.
#' @description Visualize the top gene likelihoods for each velocity clusters.
#' @param obj A seurat object that has been analyzed with entrain_velocity
#' @param velocity_clusters Character vector of velocity clusters to include in the plot. Default `NULL` - include all velocity clusters,
#' @param n_top_genes Number of top likelihood  genes to visualize. Default 10.
#' @param color_low Hex code of lower limit
#' @param color_high Hex code of upper limit.
#' @param colorscale Name of a color palette to be passed to the palette argument of scale_fill_distiller. Overrides color_low and color_high.
#' @param rescale Boolean. If TRUE, rescales all likelihoods between 0 and 1.

#' @return A gglot2 plot object.
#' @export
plot_velocity_likelihoods_heatmap <- function(obj,
                                              velocity_clusters = NULL,
                                              n_top_genes = 5,
                                              color_low = "#f6ffff",
                                              color_high = "#01bfc4",
                                              colorscale = "BuPu",
                                              rescale = TRUE) {

    res <- obj@misc$entrain$velocity_result
    lik <- res$vclust_gene_likelihoods

    if (!is.null(velocity_clusters)) {
        vcluster_likelihds <- lik[,c("gene", velocity_clusters), drop=FALSE]
    } else {
        vcluster_likelihds <- lik
    }

    likelihood_genes<-list()
    x<-vcluster_likelihds %>% dplyr::select(., -gene)
    for (col in colnames(x)) {
        idx <- vcluster_likelihds[, col] %>% order(., decreasing = TRUE) %>% head(n_top_genes)
        top_genes <- vcluster_likelihds[idx, "gene"]
        likelihood_genes <- append(likelihood_genes, top_genes)
    }

    selected_genes <- likelihood_genes %>% unlist() %>% unique()
    lik_heatmap <- vcluster_likelihds %>% dplyr::filter(., gene %in% selected_genes) %>% reshape2::melt(id.var="gene")
    min_val <- lik_heatmap$value[!is.na(lik_heatmap$value)] %>% min()
    lik_heatmap[is.na(lik_heatmap)] <- min_val
    lik_heatmap$variable %>% unique() %>% sort()

    x_order <- vcluster_likelihds[,2:ncol(vcluster_likelihds)] %>% colnames %>% sort()

    if (rescale == TRUE) {
        lik_heatmap <- lik_heatmap %>% dplyr::mutate(., value = scales::rescale(value, to=c(0,1)))
    }

    plot <- ggplot(lik_heatmap) +
        geom_tile(aes(x=factor(variable, levels=x_order), y=gene, fill=value)) +
        scale_fill_distiller(palette="BuPu", direction=1) +
        theme_minimal_big_font() +
        ggtitle("Velocity Likelihoods") +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())

    return(plot)
}

#' @title Plot sender influences overlayed on a RNA velocity plot
#' @description Visualizes the contribution of each cell in the dataset towards the observed RNA velocities.
#' @param adata Path to a `.h5ad` `anndata` object. Must have had entrain_velocity run on it.
#' @param velocity_clusters Values in `adata.obs.velocity_cluster_key` denoting which velocity clusters to visualize. The sender influences will represent contribution towards this velocity cluster only, not all velocities in the dataset.
#' @param velocity_cluster_key Column name of `adata.obs` metadata denoting velocity clusters.
#' @param frameon scvelo.pl.velocity_embedding argument: Draws a frame around the plot.
#' @param size scvelo.pl.velocity_embedding argument: Size of points.
#' @param sender_palette A matplotlib colorscale.
#' @param top_n_ligands Number of top ligands whose influence to visualize. Default 3.
#' @param title Plot title.
#' @param filename Filename of plot to save.

#' @return NULL. Saves a plot in working directory with name `filename`
#' @export
plot_sender_influence_velocity <- function(adata,
                                           velocity_clusters,
                                           velocity_cluster_key = "vcluster",
                                           top_n_ligands = 3,
                                           frameon=FALSE,
                                           size=40,
                                           sender_palette = "Reds",
                                           title = "Influence",
                                           filename = "entrain_velocity_influence.png") {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop(
            paste0("Package \"reticulate\" and a python installation with scvelo must be installed to use Entrain-Velocity.\n",
                   "You should run the following:\n",
                   "1. `install.packages(\"reticulate\")`\n",
                   "2. `reticulate::install_miniconda()`\n",
                   "3. `reticulate::conda_create(\"r-entrain\")\n",
                   "4. `reticulate::conda_install(\"r-entrain\", \"scanpy\", \"leidenalg\", \"python-igraph\", \"adjusttext\", channel=\"conda-forge\")`\n",
                   "5. `reticulate::conda_install(\"r-entrain\", \"scvelo\", channel=\"bioconda\")`\n",
                   "6. `reticulate::use_condaenv(\"r-entrain\")`",
                   "to set up your reticulate and python environment properly.", sep=""),
            call. = FALSE
        )
    }
    reticulate::source_python(system.file("python/entrain_scvelo.py", package = "entrain"))

    plot_sender_influence_velocity_python(adata = adata,
                                   velocity_cluster_key = velocity_cluster_key,
                                   velocity_clusters = velocity_clusters,
                                   top_n_ligands = top_n_ligands,
                                   frameon = frameon,
                                   size = size,
                                   sender_palette = sender_palette,
                                   title = title,
                                   filename = filename)
}

#' @title Plot top ligands overlayed on a RNA velocity plot
#' @description Visualizes the ligands most likely to be driving the observed RNA velocities.
#' @param adata Path to a `.h5ad` `anndata` object. Must have had entrain_velocity run on it.
#' @param top_n_ligands Number of top ligands to include in each label.
#' @param velocity_clusters Values in `adata.obs.velocity_cluster_key` denoting which velocity clusters and their ligands to visualize. Defaults to NULL, in which it will plot the ligands for all positive variance explained velocity clusters
#' @param velocity_cluster_key Column name of `adata.obs` metadata denoting velocity clusters.
#' @param velocity_cluster_palette A matplotlib color palette denoting colorscale of velocity clusters e.g., `"Spectral"`. OR a single color in `matplotlib.colors.CSS4_COLORS.keys()` e.g. `"black"`;  OR: A vector of hex codes with the same length as the number of velocity clusters being plotted e.g., `c("#bfe5a0", "#9e0242", "#d8434e", "#f67b4a", "#5e4fa2")`.
#' @param color scvelo.pl argument: Key for annotations of observations/cells or variables/genes. Use if you want dot color to denote a cell type or other annotation. Default is `NULL` - all cells are colored a single color denoted by argument `cell_palette`
#' @param cell_palette A matplotlib color palette for annotations denoted in argument `color`.
#' @param alpha Opacity of points.
#' @param density Density of RNA velocity vectors.
#' @param arrow_size Size of RNA velocity vector arrowhead.
#' @param linewidth_scale Width of RNA velocity vector body.
#' @param vector_type Whether to plot a `scvelo.pl.velocity_embedding_stream()` plot or a `scvelo.pl.velocity_embedding_grid()` plot. Default `stream`
#' @param label_fontsize Font size of ligand labels
#' @param size scvelo.pl.velocity_embedding argument: Size of points.
#' @param figsize Figure size. Default is `c(10,10)`
#' @param plot_output_path Filename of plot to save.

#' @return NULL. Saves a plot in working directory with name `plot_output_path`
#' @export
plot_velocity_ligands <- function(adata,
                                   color='dimgrey',
                                   alpha=0.10,
                                   density=3,
                                   arrow_size = 1.0,
                                   linewidth_scale=2,
                                   vector_type="stream",
                                   velocity_cluster_key = "vcluster",
                                   velocity_clusters = NULL,
                                   cell_palette=NULL,
                                   size=100,
                                   velocity_cluster_palette = "gnuplot",
                                   label_fontsize = 15,
                                   top_n_ligands = 5L,
                                   figsize = c(10,10),
                                   plot_output_path = "entrain_velocity_ligands.png") {

    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop(
            paste0("Package \"reticulate\" and a python installation with scvelo must be installed to use Entrain-Velocity.\n",
                   "You should run the following:\n",
                   "1. `install.packages(\"reticulate\")`\n",
                   "2. `reticulate::install_miniconda()`\n",
                   "3. `reticulate::conda_create(\"r-entrain\")\n",
                   "4. `reticulate::conda_install(\"r-entrain\", \"scanpy\", \"leidenalg\", \"python-igraph\", \"adjusttext\", channel=\"conda-forge\")`\n",
                   "5. `reticulate::conda_install(\"r-entrain\", \"scvelo\", channel=\"bioconda\")`\n",
                   "6. `reticulate::use_condaenv(\"r-entrain\")`",
                   "to set up your reticulate and python environment properly.", sep=""),
            call. = FALSE
        )
    }
    reticulate::source_python(system.file("python/entrain_scvelo.py", package = "entrain"))

    plot_velocity_ligands_python(adata=adata,
                          color=color,
                          alpha=alpha,
                          density=density,
                          arrow_size = arrow_size,
                          linewidth_scale=linewidth_scale,
                          vector_type=vector_type,
                          velocity_cluster_key = velocity_cluster_key,
                          cell_palette=cell_palette,
                          size=size,
                          velocity_clusters=velocity_clusters,
                          velocity_cluster_palette = velocity_cluster_palette,
                          label_fontsize = label_fontsize,
                          top_n_ligands = top_n_ligands,
                          figsize = figsize,
                          plot_output_path = plot_output_path)
}


#' @title Plot top ligands overlayed on a RNA velocity plot
#' @description Visualizes the ligands most likely to be driving the observed RNA velocities.
#' @param adata_clustered Path to a `.h5ad` `anndata` object. Must have had `entrain_cluster_velocities` run on it.
#' @param plot_output_path Filename of plot to save.
#' @param velocity_cluster_key Column name of `adata.obs` metadata denoting velocity clusters.
#' @param vector_type Whether to plot a `scvelo.pl.velocity_embedding_stream()` plot or a `scvelo.pl.velocity_embedding_grid()` plot. Default `stream`

#' @param ... `**kwargs` to be passed to `scvelo.pl.velocity_embedding_*`
#' @return NULL. Saves a plot in working directory with name `plot_output_path`
#' @export
plot_velocity_clusters <- function(adata_clustered,
                                   plot_output_path,
                                   velocity_cluster_key,
                                   vector_type = "stream",
                                   ...
) {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop(
            paste0("Package \"reticulate\" and a python installation with scvelo must be installed to use Entrain-Velocity.\n",
                   "You should run the following:\n",
                   "1. `install.packages(\"reticulate\")`\n",
                   "2. `reticulate::install_miniconda()`\n",
                   "3. `reticulate::conda_create(\"r-entrain\")\n",
                   "4. `reticulate::conda_install(\"r-entrain\", \"scanpy\", \"leidenalg\", \"python-igraph\", \"adjusttext\", channel=\"conda-forge\")`\n",
                   "5. `reticulate::conda_install(\"r-entrain\", \"scvelo\", channel=\"bioconda\")`\n",
                   "6. `reticulate::use_condaenv(\"r-entrain\")`",
                   "to set up your reticulate and python environment properly.", sep=""),
            call. = FALSE
        )
    }
    reticulate::source_python(system.file("python/entrain_scvelo.py", package = "entrain"))

    plot_velocity_clusters_python(adata_clustered = adata_clustered,
                           plot_file = plot_output_path,
                           velocity_cluster_key = velocity_cluster_key,
                           vector_type = vector_type,
                           ...
                           )
}

theme_minimal_big_font <- function() {
    theme_minimal() %+replace%
        theme(
            axis.text.x=element_text(size=rel(2),vjust=1),
            axis.text.y=element_text(size=rel(2)),
            axis.title.x=element_text(size=rel(2)),
            axis.title.y=element_text(size=rel(2), angle=90),
            axis.line = element_line(color="black",size=rel(2)),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_line(colour = "black"),
            axis.ticks.length = unit(0.15, "cm")
        )
}


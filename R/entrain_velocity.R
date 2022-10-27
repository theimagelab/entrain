
#' @title Identify ligands responsible for observed velocity vectors.
#' @description Given an anndata object containing clustered velocity vectors, and NicheNet ligand database, identify ligands responsible for velocity vectors. 
#' This function first calls a wrapper for scvelo's recover_dynamics function, obtaining gene-wise velocity dynamics for each cluster. Folloewd by Entrain ligand analysis.
#' Note that velocity clusters are generated from velocity data (i.e. time derivative of expression), rather than expression data like a typical usual cell clustering workflow such as Seurat FindClusters.
#' Those wishing to run this on their own pre-generated clusters (E.g. if you have a bespoke set of marker genes) can do so by inputting an anndata with the custom cluster annotations 
#' in the anndata.obs metadata, and passing the column name to cluster_key But this is not quite recommended because cells that are similar in expression could exhibit very different velocities. In turn,
#' this will input a noisy dataset to Entrain's randomForest, and any ligands that come out of such an analysis will be spurious.
#' @param receiver_obj A seurat object.
#' @param adata An anndata file containing velocities, preclustered with cluster_velocities().
#' @param sender_cluster_key Column name in `adata.obs` that denotes the sender cluster annotation column.
#' @param sender_clusters Optional. Character vector denoting sender cluster cell types. If not supplied, will launch an interactive workflow selecting them manually.
#' @param sender_obj Optional. A seurat object. Use this if your sender cells are in a different dataset than your receiver cells. E.g. if the microenvironment data is in a different file than the differentiating cells.
#' @param expressed_ligands Optional. A list of ligands expressed by the sender cells. This is calculated by default but you can specify if you have specific ligands of interest.
#' @param save_adata Optional. Filename of anndata object to write results to. If not given or NULL, anndata will not be saved. Recommended if you want to re-run or continue analysis in scanpy/scvelo.
#' @param velocity_cluster_key Column name of the cluster labels in the metadata of adata_clustered. Can be calculated manually in scvelo with scvelo.tl.velocity_clusters(). 
#' @param reduction_key Seurat reduction key for dimension reduction visualization.
#' @param num_jobs Long integer. Number of parallel jobs to run for scvelo recover dynamics function. 
#' @param num_top_genes Long integer. Number of top velocity genes to calculate likelihoods for. Default 500L.
#' @param expression_proportion_cutoff Pct cutoff to threshold a ligand as active. A ligand is 'active' if it is expressed in more than `expression_proportion_cutoff` fraction of cells in the sender cluster.
#' @param ligand_target_matrix NicheNet ligand-target data file.
#' @param lr_network NicheNet ligand-receptor pairs data file.
#' @param resolution Optional argument defining resolution of velocity clustering. Default 0.05.
#' @param ... Arguments to pass to `entrain_cluster_velocities()` and `scvelo.pl.velocity_embedding_stream()`.

#' @return a Seurat object with ligand velocity results in obj$misc$entrain$velocity_result. 
#' @export

entrain_velocity <- function(receiver_obj, sender_obj = NULL,
                             adata,
                             sender_cluster_key,
                             expressed_ligands = NULL,
                             sender_clusters = NULL,
                             reduction_key = NULL,
                             save_adata = NULL,
                             velocity_cluster_key = NULL,
                             num_jobs = 10L,
                             num_top_genes = 500L,
                             expression_proportion_cutoff = 0.10,
                             resolution = 0.05,
                             lr_network = NULL, ligand_target_matrix = NULL,
                             ...) {
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
    
    if (class(num_jobs) != "integer") {
        num_jobs <- as.integer(num_jobs)
    }
    if (identical(receiver_obj, sender_obj)) {
        error("receiver_obj is the same as sender_obj. Set sender_obj = NULL if your sender and receiver cells are contained in the single seurat object.")
    } else if (is.null(sender_obj)) {
        sender_obj <- receiver_obj
    }
    
    if (is.null(expressed_ligands)) {
        if (is.null(sender_cluster_key)) {
            error("Please supply a sender_cluster_key: The column name in metadata that denotes sender cell clusters.")
        }
        
        if (is.null(sender_clusters)) {
            message("sender_clusters is NULL, selecting sender clusters interactively.")
            assertthat::assert_that(!is.null(lr_network), msg = "If you are in interactive mode, please supply lr_network.")
            
            if(is.null(reduction_key)) {
                message("reduction_key not provided. Defaulting to umap. Please ensure a dimension reduction exists")
                reduction_key <- "umap"
            }
            expressed_ligands <- get_senders_interactive(sender_obj, sender_cluster_key, reduction_key = reduction_key, lr_network = lr_network)
            expressed_ligands <- expressed_ligands[!duplicated(expressed_ligands),]
            
        } else {
            expressed_ligands <- get_expressed_ligands(sender_obj,
                                                       sender_cluster_key = sender_cluster_key,
                                                       sender_cluster_names = sender_clusters,
                                                       lr_network = lr_network)
        }
    }
    
    receiver_obj<-intersect_seurat_adata(receiver_obj, adata)
    
    if ("vcluster" %in% names(adata$obs)) {
        adata_clustered <- adata
        message("Existing velocity clusters found as column 'vcluster' in adata.obs") 
        velocity_cluster_key <- "vcluster"
    } else if (is.null(velocity_cluster_key)) {
        message("velocity_cluster_key is NULL, proceeding to cluster velocities by default.")
        message("If you want to tune the clustering parameters, run entrain_cluster_velocities()")
        velocity_cluster_key <- "vcluster"
        
        adata_clustered <- entrain_cluster_velocities(adata,
                                                      cluster_key = velocity_cluster_key,
                                                      resolution=resolution,
                                                      ...)
        message("Velocities have been clustered, with default cluster_key = \"vcluster\"")
    } else {
        message(paste0("Using precomputed velocity clusters with argument ", velocity_cluster_key))
        adata_clustered <- adata
    }
    
    likelihood_cell_data <- recover_dynamics_clusters(adata_clustered,
                                                      cluster_key = velocity_cluster_key,
                                                      n_jobs = num_jobs,
                                                      n_top_genes = num_top_genes)

    fit_likelihoods <- likelihood_cell_data[[1]]
    velocity_clusters <- likelihood_cell_data[[2]] %>% as.data.frame() %>% stats::setNames(velocity_cluster_key)
    assertthat::assert_that(!is.null(ligand_target_matrix), msg = "Please supply ligand_target_matrix")
    obj_v <- get_velocity_ligands(receiver_obj,
                                  fit_likelihoods = fit_likelihoods,
                                  velocity_clusters = velocity_clusters,
                                  expressed_ligands = expressed_ligands,
                                  expression_proportion_cutoff = expression_proportion_cutoff,
                                  lr_network = lr_network,
                                  ligand_target_matrix = ligand_target_matrix)
    
    
    velocity_cluster_cell_ids <- data.frame(adata_clustered$obs[[velocity_cluster_key]])
    rownames(velocity_cluster_cell_ids) <- rownames(adata_clustered$obs)
    colnames(velocity_cluster_cell_ids) <- velocity_cluster_key
    
    obj_v@meta.data <- merge(obj_v@meta.data, velocity_cluster_cell_ids, by = 0 )
    message(paste("Added velocity cluster assignments to Seurat metadata column key:" , velocity_cluster_key))
    
    if (!is.null(save_adata)) {

        res <- obj_v@misc$entrain$velocity_result
        
        if (is.character(adata_clustered)) {
            anndata <- reticulate::import("anndata")
            adata_clustered <- anndata$read_h5ad(adata_clustered)
        }
        
        adata_clustered$uns$data$entrain_velocity_ligands <- res
        adata_clustered$write(save_adata)
        message(paste("Writing adata at", save_adata))
    }
    return( obj_v )
}

#' @title Cluster the velocity matrix
#' @description Use Leiden algorithm to cluster the velocities into discrete groups representing differentiation activity.
#' @param adata An anndata file containing velocities.
#' @param save_adata Optional. Filename of anndata object to write results to. If not given or NULL, anndata will not be saved. Recommended if you want to re-run or continue analysis in scanpy/scvelo.
#' @param velocity_cluster_key Column name of the cluster labels in the metadata of adata_clustered. Can be calculated manually in scvelo with scvelo.tl.velocity_clusters(). 
#' @param resolution Optional argument defining resolution of velocity clustering. Default 0.05.
#' @param plot_file Plot file name. Default is "velocity_clusters.png".
#' @param vector_type Vector types to plot. Only used with plot_file is not NULL. Either "grid" or "stream". Default "stream".
#' @param reduction_key Dimension reduction to plot. Only use when plot_file is not NULL.
#' @param ... Arguments to pass to `scvelo.pl.velocity_embedding_stream()` or `scvelo.pl.velocity_embedding_grid()`, depending on `vector_type`

#' @return Anndata object with velocity clusters denoted in adata.obs.<velocity_cluster_key>
#' @export
#' 
entrain_cluster_velocities <- function(adata,
                                       velocity_cluster_key = "vcluster",
                                       resolution = 0.05,
                                       save_adata = NULL,
                                       plot_file = "velocity_clusters.png",
                                       vector_type = "stream",
                                       reduction_key = NULL,
                                       ...) {
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
    
    message(paste0("Clustering velocity matrix with resolution: ", resolution))
    adata_clustered <- cluster_velocities(adata,
                                          cluster_key = velocity_cluster_key,
                                          resolution = resolution
                                          )
    if (!is.null(plot_file)) {
        plot_velocity_clusters(adata_clustered = adata_clustered,
                               plot_file = plot_file,
                               velocity_cluster_key = velocity_cluster_key,
                               vector_type = vector_type,
                               ...)
        message(paste0("Clustered velocities plotted: ", plot_file))
    }
    
    return(adata_clustered)
}

get_velocity_ligands <- function(obj,
                                 fit_likelihoods,
                                 velocity_clusters,
                                 expressed_ligands,
                                 expression_proportion_cutoff,
                                 lr_network, ligand_target_matrix) {
    if (is.null(obj@misc$entrain$velocity) == TRUE) {
        obj@misc$entrain$velocity <- list()
    }
    
    reticulate::source_python(system.file("python/entrain_scvelo.py", package = "entrain"))
    
    vclust_names <- velocity_clusters[,1] %>% unique()
    
    vclust_df <- data.frame(matrix(NA, nrow=0, ncol=1))
    colnames(vclust_df)<-"ligand"
    vclust_df$ligand <- vclust_df$ligand %>% as.character()
    
    likelihood_df <- data.frame(matrix(NA, nrow=0, ncol=1))
    colnames(likelihood_df) <- "gene"
    likelihood_df$gene <- likelihood_df$gene %>% as.character()
    
    var_exp_df <- data.frame(matrix(NA, nrow=0, ncol=2))
    colnames(var_exp_df) <- c("cluster", "variance_explained")
    for (vclust in vclust_names) {
        likelihoods <- fit_likelihoods %>% dplyr::select(vclust)
        likelihoods <- likelihoods %>% tidyr::drop_na() 
        likelihoods <- likelihoods[ order(rownames(likelihoods)), , drop=FALSE]
        vclust_cell_names <- velocity_clusters[velocity_clusters[,1] == vclust, , drop=FALSE] %>% rownames()
        vclust_seuobj <- subset(obj, cells=vclust_cell_names)
        active_lr <- get_active_ligand_receptor(receivers_seuobj = vclust_seuobj,
                                                expressed_ligands = expressed_ligands,
                                                lr_network = lr_network,
                                                ligand_target_matrix = ligand_target_matrix,
                                                expression_proportion_cutoff = expression_proportion_cutoff)
        active_ligands <- active_lr$from %>% unique()
        active_receptors <- active_lr$to %>% unique()
        active_ligand_potentials <- ligand_target_matrix[,active_ligands]
        
        rf_model <- get_ligand_trajectory_scores_regression(likelihoods, active_ligand_potentials)
        names(rf_model) <- c(vclust, "model")
        importances <- sort(rf_model$model$importance[,"IncNodePurity"], decreasing=TRUE)
        var_exp <- rf_model$model$rsq %>% tail(1)
        
        obj@misc$entrain$velocity_model[[vclust]] <- list(model = rf_model, 
                                                          velocity_ligands_importances = importances,
                                                          top_likelihood_genes = likelihoods,
                                                          variance_explained = var_exp,
                                                          active_network = active_lr)
        
        importances <- importances %>% as.data.frame %>% setNames(vclust)
        importances$ligand <- rownames(importances)
        vclust_df <- dplyr::full_join(vclust_df, importances, by="ligand")
        
        likelihoods$gene <- rownames(likelihoods)
        likelihood_df <- dplyr::full_join(likelihood_df, likelihoods, by="gene")
        
        var_exp_df <- rbind(var_exp_df, c(vclust, var_exp))
    }
    
    obj@misc$entrain$velocity_result <- list(vclust_ligand_importances = vclust_df, vclust_gene_likelihoods = likelihood_df, variance_explained = var_exp_df)
    return(obj)
}


#' @title Subset cells in Seurat and Anndata objects such that they have identical cellIDs.
#' @description Given a Seurat object obj, and an anndata object adata, calculates the intersection of cell names 
#' between colnames(obj) and adata.obs_names, and subsets both objects to their intersection.
#' @param obj A seurat object
#' @param adata An anndata file.

#' @return a Seurat object subset to cells that it has in common with adata.

#' @export
intersect_seurat_adata <- function(obj,
                                 adata) {
    reticulate::source_python(system.file("python/entrain_scvelo.py", package = "entrain"))
    
    seurat_cells <- colnames(obj)
    adata_cells <- get_adata_obs_names(adata)
    
    intersection_cells <- intersect(seurat_cells, adata_cells)
    int_len <- length(intersection_cells)
    seu_len <- length(seurat_cells)
    adata_len <- length(adata_cells)
    
    if (seu_len != adata_len) {
        message(paste0("n cells in seurat object: ", seu_len, ". n cells in anndata object: ", adata_len))
        if (int_len != seu_len) {
            diff <- seu_len - int_len
            message("Seurat object and anndata object contain different cells")
            message("Subsetting Seurat object to cells that it has in common with adata..")
            obj<-subset(obj, cells = intersection_cells)
            message(paste0("Subsetting finished: ", diff, " cells removed."))
        }

        if (int_len != adata_len) {
            stop("There are cells in adata (adata.obs_names) that are not in the Seurat object (colnames(obj)).")
        }
    } else if (!all(seurat_cells == adata_cells)) {
        message("Warning: Seurat object and anndata object contain the same cells but in different order.")
    }
    
    return(obj)
}

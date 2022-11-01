
import anndata as ad
import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import os

def cluster_velocities(h5ad_path,
                       resolution=0.05,
                       prefix="vcluster",
                       cluster_key="vcluster"):
    adata = ad.read_h5ad(h5ad_path)

    v_dat = adata.layers["velocity"]
    v = sc.AnnData(v_dat)
    sc.pp.pca(v)
    sc.pp.neighbors(v)
    sc.tl.umap(v)
    sc.tl.leiden(v, resolution=resolution)


    v.obs["leiden"] = prefix + v.obs["leiden"].astype(str)
    v_clusts = v.obs["leiden"].to_frame()
    v_clusts = v_clusts.set_index(adata.obs_names.values)
    adata.obs[cluster_key] = v_clusts
    print("Velocity clusters found: ")
    value_counts = adata.obs[cluster_key].value_counts()
    print(value_counts)
    significant_vclusters = value_counts[value_counts > (0.01 * adata.n_obs)]
    print("Found " + str(len(significant_vclusters)) + " significant velocity clusters (clusters with more than 1% of total cells.)")


    return adata


def recover_dynamics_clusters(adata,
                              return_adata = False,
                              cluster_key="vcluster",
                              **kwargs):
    if type(adata) == str:
        adata = ad.read_h5ad(adata)

    vclust_likelihoods = []
    vclusts = adata.obs[cluster_key].unique()
    for vclust in vclusts:
        vclust_cell_idx = np.where(adata.obs[cluster_key] == vclust)
        ad_clust = adata[vclust_cell_idx].copy()
        scv.tl.recover_dynamics(ad_clust, **kwargs)

        top_likelihoods_clust = ad_clust.var['fit_likelihood']
        top_likelihoods_clust.rename(vclust, inplace=True)
        vclust_likelihoods.append(top_likelihoods_clust)


    likelihoods_df = pd.concat(vclust_likelihoods, axis=1)
    cluster_df = pd.DataFrame(adata.obs[cluster_key])

    if (return_adata == True):
        adata.uns["velocity_cluster_likelihoods"] = likelihoods_df
        return(adata)
        
    return([likelihoods_df, cluster_df])

def plot_sender_influence_velocity_python(adata, 
                                  sender_adata = None,
                                  velocity_clusters = None, 
                                  velocity_cluster_key = "velocity_clusters", 
                                  top_n_ligands=3, 
                                  vkey = "velocity",
                                  basis="umap",
                                  arrow_color="black",
                                  dpi=300,
                                  figsize = None,
                                  size=15,
                                  density=None,
                                  arrow_length=2,
                                  arrow_size=3,
                                  velocity_cluster_palette = "Spectral",
                                  sender_palette = "OrRd",
                                  filename = "entrain_velocity_influence.png",
                                  vector_type = "stream",
                                  signif_cluster_pct_cutoff = 0.01,
                                  frameon=None,
                                  **kwargs):
    import scvelo as scv
    import matplotlib.pyplot as plt
    import math
    if type(adata) == str: adata = ad.read_h5ad(adata)
    if type(sender_adata) == str: sender_adata = ad.read_h5ad(sender_adata) 
    
    if velocity_clusters:
        if np.all(np.isin(velocity_clusters, adata.obs[velocity_cluster_key].unique())):
            clusters = velocity_clusters
        else:
            raise ValueError(velocity_clusters + "not found in values of column: " + velocity_cluster_key)
    elif signif_cluster_pct_cutoff > 0.0:
        vcluster_counts = adata.obs[velocity_cluster_key].value_counts()
        clusters = vcluster_counts[vcluster_counts > (signif_cluster_pct_cutoff * adata.n_obs)].index.to_numpy()
    else:    
         clusters = adata.obs[velocity_cluster_key].unique()          
    
    velocity_basis=vkey + "_" + basis
    if velocity_basis not in adata.obsm.keys():
        scv.tl.velocity_embedding(adata, basis = basis, vkey = vkey)

    if vector_type not in ["stream", "grid"]:
        raise ValueError("vector_type must be one of \"stream\" or \"grid\".")

    adjust_for_stream = True if vector_type == "stream" else False

    fig = plt.figure()
    if sender_adata:
        
        nrows=2
        ncols=len(clusters)
        ax = []
        
        figsize = figsize if figsize is not None else (8*ncols, 16)
        gs = plt.GridSpec(nrows,
                          ncols,
                          plt.figure(figsize=figsize)
                          ) 
        for i, velocity_cluster in enumerate(clusters):
            velocity_cluster = str(velocity_cluster)
            print("Plotting velocities and ligands for cluster: " + velocity_cluster)
            cells = adata.obs.loc[adata.obs[velocity_cluster_key] == velocity_cluster].index
            cells_idx = np.where(adata.obs_names.isin(cells))[0].tolist()
        
            X_grid, V_grid = get_embedding_subset_grid(adata = adata,
                                                       cells_idx = cells_idx,
                                                       basis = basis,
                                                       vkey = vkey,
                                                       adjust_for_stream=adjust_for_stream
                                                       )
            
            if vector_type == "grid":
                ax.append(
                    scv.pl.velocity_embedding_grid(adata, 
                                                      palette=sender_palette,
                                                      color = velocity_cluster_key,
                                                      dpi=300, 
                                                      size=size,
                                                      arrow_length=arrow_length,
                                                      arrow_size=arrow_size,
                                                      density=density,
                                                      arrow_color=arrow_color, 
                                                      show=False,
                                                      ax=plt.subplot(gs[0,i]),
                                                      X_grid = X_grid,
                                                      V_grid = V_grid,
                                                      autoscale=True,
                                                      frameon=frameon,
                                                      title = velocity_cluster
                                                      )   
                    )
            elif vector_type == "stream":
                ax.append(
                    scv.pl.velocity_embedding_stream(adata, 
                                                      palette=sender_palette,
                                                      color = velocity_cluster_key,
                                                      dpi=300, 
                                                      size=size,
                                                      arrow_size=arrow_size,
                                                      density=density,
                                                      arrow_color=arrow_color, 
                                                      show=False,
                                                      ax=plt.subplot(gs[0,i]),
                                                      X_grid = X_grid,
                                                      V_grid = V_grid,
                                                      frameon=frameon,
                                                      title = velocity_cluster
                                                      )   
                    )       
            entrain_dat = adata.uns["entrain_velocity_ligands"]
            ligand_importances = entrain_dat["vclust_ligand_importances"][velocity_cluster].fillna(0)
            ligand_importances.index = entrain_dat["vclust_ligand_importances"]["ligand"]
            ligand_importances = ligand_importances.sort_values(ascending=False)
            ligands = ligand_importances.index.values
            
            sender_adata.obs["influence_"+velocity_cluster] = get_sender_influence(sender_adata, ligands)

            ax.append(
                sc.pl.umap(sender_adata,
                       color="sender_influence",
                       label="on data",
                       palette=sender_palette,
                       title=velocity_cluster,
                       show=False,
                       ax=plt.subplot(gs[1,i])
                       )
                )
            
        ax[0].figure.savefig(filename, dpi=dpi)

    else:
        nrows = 1 if len(clusters) < 5 else 1 + (math.floor(len(clusters) / 4) )
        ncols = len(clusters) if len(clusters) < 4 else 4
        ax = []
        
        figsize = figsize if figsize is not None else (8*ncols, 8*nrows)
        gs = plt.GridSpec(nrows,
                          ncols,
                          plt.figure(figsize=figsize)
                          ) 
        for i, velocity_cluster in enumerate(clusters):
            velocity_cluster = str(velocity_cluster)
            gs_row = math.floor(i/4)
            gs_col = i % 4
            
            print("Plotting velocities and ligands for cluster: " + velocity_cluster)
            if velocity_cluster not in adata.obs[velocity_cluster_key].unique():
                raise ValueError(velocity_cluster + " not found in column name: " + velocity_cluster_key)
            cells = adata.obs.loc[adata.obs[velocity_cluster_key] == velocity_cluster].index
            cells_idx = np.where(adata.obs_names.isin(cells))[0].tolist()
        
            X_grid, V_grid = get_embedding_subset_grid(adata = adata,
                                                       cells_idx = cells_idx,
                                                       basis = basis, vkey = vkey,
                                                       adjust_for_stream = adjust_for_stream)
            
            entrain_dat = adata.uns["entrain_velocity_ligands"]
            ligand_importances = entrain_dat["vclust_ligand_importances"][velocity_cluster].fillna(0)
            ligand_importances.index = entrain_dat["vclust_ligand_importances"]["ligand"]
            ligand_importances = ligand_importances.sort_values(ascending=False)
            ligands = ligand_importances.index.values
            
            influence_col = "influence" + velocity_cluster
            adata.obs[influence_col] = get_sender_influence(adata, ligands)

            if vector_type == "grid":
                ax.append(
                    scv.pl.velocity_embedding_grid(adata, 
                                                   color = influence_col,
                                                   dpi=dpi, 
                                                   size=size,
                                                   arrow_length=arrow_length,
                                                   arrow_size=arrow_size,
                                                   color_map=sender_palette,
                                                   density=density,
                                                   arrow_color=arrow_color, 
                                                   ax=plt.subplot(gs[gs_row,gs_col]),
                                                   X_grid = X_grid,
                                                   V_grid = V_grid,
                                                   title = velocity_cluster,
                                                   autoscale=True,     
                                                   frameon=frameon,
                                                   show=False)
                    )
            elif vector_type == "stream":
                ax.append(
                    scv.pl.velocity_embedding_stream(adata, 
                                                   color = influence_col,
                                                   dpi=dpi, 
                                                   size=size,
                                                   arrow_size=arrow_size,
                                                   color_map=sender_palette,
                                                   density=density,
                                                   arrow_color=arrow_color, 
                                                   ax=plt.subplot(gs[gs_row,gs_col]),
                                                   X_grid = X_grid,
                                                   V_grid = V_grid,
                                                   title = velocity_cluster,
                                                   frameon=frameon,
                                                   show=False)
                    ) 

            ax[0].figure.savefig(filename, dpi=dpi)

    print("Saved at "+ os.getcwd() + "/" + filename)

    return(None)


def plot_velocity_ligands_python(adata,
                          color = None,
                          top_n_ligands=3,
                          vector_type = "stream",
                          velocity_cluster_key="velocity_clusters",
                          velocity_clusters = None,
                          vkey="velocity",
                          basis="umap",
                          velocity_cluster_palette = "Spectral",
                          dpi=300,
                          arrow_length=1,
                          linewidth_scale = None,
                          density=0.5,
                          size = 15,
                          alpha=0.2,
                          arrow_size=1,
                          figsize=(16,9),
                          label_fontsize = 10,
                          cell_palette="tab20",
                          signif_cluster_pct_cutoff = 0.01,
                          plot_output_path="entrain_plot.png",
                          plot_negative_VE = False,
                          **kwargs):

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib import cm
    from adjustText import adjust_text
    if type(adata) == str:
        adata = ad.read_h5ad(adata)
    
    top_n_ligands = int(top_n_ligands)  

    if velocity_clusters is not None:
        if np.all(np.isin(velocity_clusters, adata.obs[velocity_cluster_key].unique())):
            velocity_clusters = velocity_clusters
        else:
            raise ValueError(velocity_clusters + "not found in values of column: " + velocity_cluster_key)
    elif signif_cluster_pct_cutoff > 0.0:
        vcluster_counts = adata.obs[velocity_cluster_key].value_counts()
        velocity_clusters = vcluster_counts[vcluster_counts > (signif_cluster_pct_cutoff * adata.n_obs)].index.to_numpy()
    else:    
        velocity_clusters = adata.obs[velocity_cluster_key].unique()          

    if plot_negative_VE == False:
        var_exp_df =  adata.uns["entrain_velocity_ligands"]["variance_explained"]
        var_exp_df.columns = [velocity_cluster_key, "var_exp"]
        positive_VE_idx = var_exp_df["var_exp"].astype("float") > 0
        positive_ve_clusters = var_exp_df.loc[positive_VE_idx, velocity_cluster_key].tolist()
        velocity_clusters = np.intersect1d(velocity_clusters, positive_ve_clusters)
        if len(velocity_clusters) == 0 :
            raise ValueError("None of the velocity clusters had a positive VE. If you wish to visualize negative VE clusters, set plot_negative_VE = True.")

    label_df = get_velocity_cluster_label_positions(adata = adata, 
                                                    velocity_clusters = velocity_clusters,
                                                    velocity_cluster_key = velocity_cluster_key,
                                                    top_n_ligands = top_n_ligands,
                                                    figsize = figsize,
                                                    dpi=dpi
                                                    )
    plt.clf()

    num_velo_clusters = len(velocity_clusters)  

    emb_x_range = adata.obsm[f"X_{basis}"][:,0].max() - adata.obsm[f"X_{basis}"][:,0].min()
    emb_y_range = adata.obsm[f"X_{basis}"][:,1].max() - adata.obsm[f"X_{basis}"][:,1].min()      

    ax = plt.axes()

    # arrow colors
    if type(velocity_cluster_palette) == str:
        arrow_colormap = cm.get_cmap(velocity_cluster_palette, num_velo_clusters)
        arrow_colors = [arrow_colormap(float(i) / num_velo_clusters) for i in range(num_velo_clusters)]
    elif type(velocity_cluster_palette) == list:
        arrow_colors = velocity_cluster_palette
        if len(velocity_cluster_palette) != num_velo_clusters:
            raise ValueError("Number of arrow colors: " + str(arrow_colors) + " does not match the number of velocity clusters being plotted: " + str(velocity_clusters))
        
    adjust_for_stream = True if vector_type=="stream" else False
    for i,velocity_cluster in enumerate(velocity_clusters):
        cells = adata.obs.loc[adata.obs[velocity_cluster_key] == velocity_cluster].index
        cells_idx = np.where(adata.obs_names.isin(cells))[0].tolist()


        if i>0:
            size = 0 # Plot scatter points (cells, not vectors) once.

        if vector_type=="stream": 
                    
            X_grid, V_grid = get_embedding_subset_grid(adata = adata,
                                                       cells_idx = cells_idx,
                                                       basis = basis, vkey = vkey,
                                                       adjust_for_stream = adjust_for_stream)
            
            lengths = np.sqrt((V_grid ** 2).sum(0))
            linewidth = 1 if linewidth_scale is None else linewidth_scale
            linewidth *= 2 * lengths / lengths[~np.isnan(lengths)].max()
            
            grid_x_range = X_grid[0,:].max() - X_grid[0,:].min()
            grid_y_range = X_grid[1,:].max() - X_grid[1,:].min()

            # maintains consistent arrow plot density across clusters.
            density_x = density * (grid_x_range/emb_x_range) 
            density_y = density * (grid_x_range/emb_y_range)
            plot = scv.pl.velocity_embedding_stream(adata, 
                                           color = color,
                                           alpha = alpha,
                                           X_grid = X_grid,
                                           V_grid = V_grid,
                                           dpi = dpi, 
                                           size = size,
                                           linewidth = linewidth,
                                           density=(density_x, density_y),
                                           arrow_size = arrow_size,
                                           palette = cell_palette,
                                           arrow_color = arrow_colors[i], 
                                           ax = ax,
                                           figsize = figsize,
                                           legend_loc = 'none',
                                           show = False, **kwargs)
        elif vector_type == "grid":
            
            
            X_grid, V_grid = get_embedding_subset_grid(adata = adata,
                                                       density = density,
                                                       cells_idx = cells_idx,
                                                       basis = basis, vkey = vkey,
                                                       adjust_for_stream = adjust_for_stream)
            
            grid_x_range = X_grid[0,:].max() - X_grid[0,:].min()
            grid_y_range = X_grid[1,:].max() - X_grid[1,:].min()
    
            # maintains consistent arrow plot density across clusters.
            density_x = density * (grid_x_range/emb_x_range) 

            plot = scv.pl.velocity_embedding_grid(adata, 
                                           color = color,
                                           alpha=alpha,
                                           X_grid = X_grid,
                                           V_grid = V_grid,
                                           dpi = dpi, 
                                           size = size,
                                           density=density_x,
                                           arrow_size = arrow_size,
                                           palette = cell_palette,
                                           arrow_color = arrow_colors[i], 
                                           ax = ax,
                                           figsize = figsize,
                                           legend_loc = 'none',
                                           show = False, **kwargs)
        else:
            raise ValueError("vector_type needs to be one of \"stream\" or \"grid\". ")

    for t in plot.texts:
        t.set_visible(False)

    if top_n_ligands > 0:
        text = []
        for i in range(label_df.shape[0]):
            t = plot.text(x=label_df.loc[i, "pos_x"],
                          y=label_df.loc[i, "pos_y"],
                          s=label_df.loc[i, "label"],
                          bbox=dict(facecolor='white', alpha=0.7),
                          fontsize=label_fontsize)

            text.append(t)
        adjust_text(text)
    fig = plot.get_figure()
    fig.savefig(plot_output_path)
    return(None)

def get_velocity_cluster_label_positions(adata,
                                         velocity_clusters,
                                         velocity_cluster_key,
                                         top_n_ligands,
                                         figsize,
                                         dpi):
    
    import matplotlib as mpl
    labels_plot = scv.pl.velocity_embedding_stream(adata, 
                                       color = velocity_cluster_key,
                                       dpi=dpi, 
                                       figsize=figsize,
                                       legend_loc = "on data",
                                       show=False)
    
    if  adata.uns["entrain_velocity_ligands"] is not None:
        try: 
            entrain_dat = adata.uns["entrain_velocity_ligands"]
            vclust_ligand_importances = entrain_dat["vclust_ligand_importances"]
            vclust_likelihoods = entrain_dat["vclust_gene_likelihoods"]
            vclust_var_exp = entrain_dat["variance_explained"]
            vclust_var_exp.columns = [velocity_cluster_key, "var_exp"]
        except: 
            raise ValueError("adata.uns[\"entrain_velocity_ligands\" has unrecognized names in it's structure. Did Entrain run correctly or have you changed list element names?")
    else:
        raise ValueError("adata.uns[\"vclust_ligand_importances\"] not found. Have you run entrain_velocity on this dataset yet?")

    if velocity_clusters is not None:
        velocity_clusters = velocity_clusters 
    else:
        velocity_clusters = adata.obs[velocity_cluster_key].unique()

    text_labels = [c for c in labels_plot.get_children() if isinstance(c, mpl.text.Text)]
    cluster_label_positions = [c.get_position()
                               for c in text_labels if c.get_text() in velocity_clusters]
    clusters_labels = [c.get_text()
                       for c in text_labels if c.get_text() in velocity_clusters]
    label_df = pd.DataFrame(
        clusters_labels, cluster_label_positions).reset_index()
    label_df[["pos_x", "pos_y"]] = label_df.iloc[:, 0].apply(
        lambda x: pd.Series([x[0], x[1]]))
    label_df["pos_y"] = label_df["pos_y"] - 1
    label_df = label_df.drop("index", 1)
    label_df = label_df.rename({0: velocity_cluster_key}, axis=1)
    label_df = label_df.merge(vclust_var_exp, on = velocity_cluster_key)
    
    temp = []
    for clust in clusters_labels:
        ligands = vclust_ligand_importances.loc[:, ["ligand", clust]].sort_values(by=clust, ascending=False)
        top_ligands = ligands.head(top_n_ligands)["ligand"].values.astype(str).tolist()
        top_ligands.append(clust)
        temp.append(top_ligands)
    top_ligand_df = pd.DataFrame(temp)
    vclust_idx = len(top_ligand_df.columns)-1
    top_ligand_df.rename(columns= {top_ligand_df.columns[vclust_idx] : velocity_cluster_key}, inplace=True)
    label_df = label_df.merge(top_ligand_df, on=velocity_cluster_key)
    label_df["var_exp"] = label_df["var_exp"].astype("float")

    label_df["var_exp"] = round(
        label_df.loc[:, "var_exp"].multiply(100.0), 1).astype(str) + "% V.E."

    label_df["label"] = label_df.iloc[:, 4:].apply(
        lambda x: ', \n'.join(x), axis=1) + " \n" + label_df["var_exp"]
    
    
    return(label_df)

def plot_velocity_clusters(adata_clustered, plot_file, velocity_cluster_key, vector_type="stream", **kwargs):

    if vector_type == "stream": 
                    
        plot = scv.pl.velocity_embedding_stream(adata_clustered, 
                                           color = velocity_cluster_key,
                                           save = plot_file, **kwargs)
    elif vector_type == "grid":
        plot = scv.pl.velocity_embedding_grid(adata_clustered, 
                                           color = velocity_cluster_key,
                                           save = plot_file, **kwargs)
    else:
        raise ValueError("vector_type needs to be one of \"stream\" or \"grid\". ")

    return None


def get_sender_influence(adata, ligands, top_n_ligands = 3):
    sender_ligand_X = adata[:,ligands[0:top_n_ligands]].to_df()
    sender_influence = np.array(sender_ligand_X.sum(axis=1))
    sender_influence_col= pd.Series(sender_influence)
    sender_influence_col.index = adata.obs_names
    
    return(sender_influence_col)
    

def get_embedding_subset_grid(adata, cells_idx, vkey, basis, adjust_for_stream):
    from scvelo.plotting.velocity_embedding_grid import compute_velocity_on_grid
    obsm = adata.obsm
    
    if len(cells_idx) < 51:
        raise ValueError("Too few (less than 50 cells) found in velocity cluster. Ensure your clusters have more than 50 cells or exclude them from analysis with velocity_clusters argument.")
    X_emb = np.array(obsm[f"X_{basis}"][cells_idx,:])
    V_emb = np.array(obsm[f"{vkey}_{basis}"][cells_idx, :])
    
    X_grid, V_grid = compute_velocity_on_grid(
                    X_emb=X_emb,
                    V_emb=V_emb,
                    adjust_for_stream=adjust_for_stream
                    )
    return(X_grid, V_grid)

def get_adata_obs_names(adata):

    if type(adata) == str:
        adata = ad.read_h5ad(adata)

    return(adata.obs_names.values)     

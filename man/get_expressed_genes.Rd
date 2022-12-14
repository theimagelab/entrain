% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrain_fitting.R
\name{get_expressed_genes}
\alias{get_expressed_genes}
\title{Determine expressed genes of a cell type.}
\usage{
get_expressed_genes(ident, seurat_obj, pct = 0.1, assay_oi = NULL)
}
\arguments{
\item{ident}{Name of cluster identity/identities of cells}

\item{seurat_obj}{Single-cell expression dataset as Seurat object https://satijalab.org/seurat/.}

\item{pct}{We consider genes expressed if they are expressed in at least a specific fraction of cells of a cluster. This number indicates this fraction. Default: 0.10. Choice of this parameter is important and depends largely on the used sequencing platform. We recommend to require a lower fraction (like the default 0.10) for 10X data than for e.g. Smart-seq2 data.}

\item{assay_oi}{If wanted: specify yourself which assay to look for. Default this value is NULL and as a consequence the 'most advanced' assay will be used to define expressed genes.}
}
\value{
A character vector with the gene symbols of the expressed genes
}
\description{
NicheNetR function to determine expressed genes of a cell type from a Seurat object single-cell RNA seq dataset or Seurat spatial transcriptomics dataset.
Function is modified from \code{nichenetr} v1.1.1 (GPL3.0): Browaeys, R., Saelens, W. & Saeys, Y. (2019). Nat Methods.
Modified to deal with non-sparse matrices, and avoids the (non-trivial) installation of NicheNet.
}

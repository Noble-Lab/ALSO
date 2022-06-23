########################################
# Utility functions
########################################
# Loads 10x dataset into a sparse matrix object
# Args:
#   matrix_fn (string): name of mtx or mtx.gz file coontaining MM formatted matrix
#   peaks_fn (string): name of peaks BED file (matches rows of matrix)
#   barcodes_fn (string): name of file containing cell barcodes (matches columns of matrix)
# Returns:
#   sparse matrix: matrix represented by files above with row and column names set (chr_start_stop format used for rownames)
load_tenx_atac = function(matrix_fn, peaks_fn, barcodes_fn) {
  atac_matrix = readMM(matrix_fn)
  colnames(atac_matrix) = read.delim(barcodes_fn, header=FALSE)$V1
  peaks = read.delim(peaks_fn, header=FALSE)
  peaks = paste(peaks$V1, peaks$V2, peaks$V3, sep = '_')
  rownames(atac_matrix) = peaks
  return(atac_matrix)
}

# Allows filtering of sites measured as non-zero in less than a given number of cells
# Args:
#   bmat (sparse matrix): sparse matrix (binarized)
#   cells (int): filter sites if they have less than this number of cells above zero
# Returns:
#   sparse matrix: filtered sparse matrix
filter_features = function(bmat, cells=10) {
  bmat = bmat[Matrix::rowSums(bmat) >= cells, ]
  return(bmat)
}

# Allows filtering of cells with below a given number of non-zero features
# Args:
#   bmat (sparse matrix): sparse matrix (binarized)
#   features_above_zero (int): filter cells if they have less than this number of features above zero
# Returns:
#   sparse matrix: filtered sparse matrix
filter_cells = function(bmat, features_above_zero=100) {
  bmat = bmat[, Matrix::colSums(bmat > 0) >= features_above_zero]
  return(bmat)
}

# Takes sparse matrix object and downsamples to a given fraction of entries remaining.
# Args:
#   bmat (sparse matrix): sparse matrix to downsample
#   fraction_remaining (float): float (0, 1) that indicates fraction of non-zero entries to retain
#   cells_per_site_min (int): min cells a site must be measured in to retain the site in matrix
#   sites_per_cell_min (int): min sites a cell must have non-zero entries in to retain the cell in matrix
# Returns:
#   sparse matrix: downsampled sparse matrix
downsample_matrix = function(bmat, fraction_remaining=0.5, cells_per_site_min=1, sites_per_cell_min=1) {
  set.seed(2019)
  non_zero_entries = which(bmat@x > 0)
  indices_to_zero = sample(non_zero_entries, size=ceiling(length(non_zero_entries) * (1 - fraction_remaining)))
  bmat@x[indices_to_zero] = 0
  
  # Make sure to get rid of stuff that has gone to ~ 0 after downsampling
  bmat = filter_features(bmat, cells=cells_per_site_min)
  bmat = filter_cells(bmat, features_above_zero=sites_per_cell_min)
  return(bmat)
}

########################################
# Functions for LSI
########################################
# Helper function to do fast version of row scaling of sparse TF matrix by IDF vector.
# Exploits the way that data is stored within sparseMatrix object. Seems to be much more memory efficient than tf * idf and faster than DelayedArray.
# Args:
#   tf (sparse matrix): term frequency matrix
#   idf (vector): IDF vector
# Returns:
#   sparse matrix: TF-IDF matrix
safe_tfidf_multiply = function(tf, idf) {
  tf = t(tf)
  tf@x <- tf@x * rep.int(idf, diff(tf@p))
  tf = t(tf)
  return(tf)
}

# Perform TF-IDF on binary matrix
# Args:
#   bmat (sparse matrix): sparse matrix to downsample
#   frequencies (bool): divide bmat by colSums (if FALSE simply use bmat for TF matrix)
#   log_scale_tf (bool): log scale TF matrix if TRUE
#   scale_factor (float): multiply terms in TF matrix by scale_factor prior to log1p. Equivalent to adding small pseudocount but doesn't cast to dense matrix at any point.
# Returns:
#   sparse matrix: TF-IDF matrix
tfidf = function(bmat, frequencies=TRUE, log_scale_tf=TRUE, scale_factor=100000) {
  # Use either raw counts or divide by total counts in each cell
  if (frequencies) {
    # "term frequency" method
    tf = t(t(bmat) / Matrix::colSums(bmat))
  } else {
    # "raw count" method
    tf = bmat
  }
  
  # Either TF method can optionally be log scaled
  if (log_scale_tf) {
    if (frequencies) {
      tf@x = log1p(tf@x * scale_factor)
    } else {
      tf@x = log1p(tf@x * 1)
    }
  }
  
  # IDF w/ "inverse document frequency smooth" method
  idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  
  # TF-IDF
  tf_idf_counts = safe_tfidf_multiply(tf, idf)
  rownames(tf_idf_counts) = rownames(bmat)
  colnames(tf_idf_counts) = colnames(bmat)
  return(tf_idf_counts)
}

# Perform current version of TF-IDF used by 10x on binary matrix
# Args:
#   bmat (sparse matrix): sparse matrix to downsample
# Returns:
#   sparse matrix: TF-IDF matrix
tenx_tfidf = function(bmat) {
  idf = log(ncol(bmat) + 1) - log(1 + Matrix::rowSums(bmat))
  tf_idf_counts = safe_tfidf_multiply(bmat, idf)
  
  rownames(tf_idf_counts) = rownames(bmat)
  colnames(tf_idf_counts) = colnames(bmat)
  tf_idf_counts = as(tf_idf_counts, "sparseMatrix")
  return(tf_idf_counts)
}

# Perform fast PCA (irlba) on matrix, retaining observation names
# Args:
#   mat (sparse matrix): matrix to use for PCA (no further scaling or centering done)
#   dims (int): number of PCs to calculate
# Returns:
#   sparse matrix: TF-IDF matrix
do_pca = function(mat, dims=50) {
  pca.results = irlba(t(mat), nv=dims)
  final_result = pca.results$u %*% diag(pca.results$d)
  rownames(final_result) = colnames(mat)
  colnames(final_result) = paste0('PC_', 1:dims)
  return(final_result)
}

########################################
# Helper functions for dim reduction
########################################
# Wrapper for performing further dim reduction (tSNE/UMAP) and clustering given PCA space via Seurat.
# Args:
#   atac_matrix (sparse matrix): matrix to store in Seurat object (not used in computations)
#   cell_embeddings (matrix): typically PCA coordinates of cells but could be any set of reduced dim coordinates
#   dims (vector of int): vector of dims to use from cell_embeddings in downstream analysis
#   metadata (dataframe): dataframe of metadata (rowonames are cell names) to add to Seurat object
#   reduction (string): reduction to use for downstream steps. Can be 'pca' (cell_embeddings) or 'pca.l2' (L2 normalized cell_embeddings)
# Returns:
#   Seurat object: seurat object
run_dim_reduction = function(atac_matrix, cell_embeddings, dims, metadata=NULL, reduction='pca.l2', n.neighbors, n.components, min.dist, k.param) {
  if (is.null(metadata)) {
    seurat_obj = Seurat::CreateSeuratObject(atac_matrix)
  } else {
    seurat_obj = Seurat::CreateSeuratObject(atac_matrix, meta.data = metadata)
  }
  
  seurat_obj[['pca']] = Seurat::CreateDimReducObject(embeddings=cell_embeddings, key='PC_', assay='RNA')
  seurat_obj = seurat_obj %>%
    Seurat::L2Dim(reduction='pca') %>%
    Seurat::RunUMAP(reduction = reduction, dims = dims, n.neighbors=n.neighbors, n.components=n.components, min.dist=min.dist) %>%
    Seurat::RunTSNE(reduction = reduction, dims = dims) %>%
    Seurat::FindNeighbors(reduction=reduction, nn.eps=0, dims=dims, k.param=k.param)
  return(seurat_obj)
}

# Helper function for plotting Spearman correlations of a given metadata column with all dimensions in a reduced space.
# Args:
#   seurat_obj (seurat object): Seurat object to use
#   reduction (string): name of reduction to use
#   column (string): name of column in metadata to use
# Returns:
#   ggplot object: plot object
plot_pc_correlation = function(seurat_obj, reduction, column='nCount_RNA') {
  coords = Seurat::Embeddings(seurat_obj, reduction=reduction)
  column_value = seurat_obj@meta.data[, column]
  correlations = apply(coords, 2, function(x) {cor(x, column_value, method='spearman')})
  correlations_df = data.frame(correlation=correlations, PC=1:ncol(coords))
  
  plot_obj = ggplot(correlations_df, aes(PC, correlation)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept = 0, linetype='dashed', color='red')
  
  return(plot_obj)
}

########################################
# Wrapper functions for workflows
########################################
# Wrapper for full LSI workflow (TF-IDF and PCA + clustering + further dim reduction)
# Args:
#   bmat (sparse matrix): sparse matrix (binarized)
#   dims (vector of int): vector of dims to use from cell_embeddings in downstream analysis
#   metadata: dataframe of metadata (rowonames are cell names) to add to Seurat object
#   log_scale_tf (bool): log scale TF matrix if TRUE
#   reduction (string): reduction to use for downstream steps. Can be 'pca' (cell_embeddings) or 'pca.l2' (L2 normalized cell_embeddings)
#   resolution (float): resolution parameter to Seurat Louvain clustering
# Returns:
#   Seurat object: Seurat object. clustering + tSNE + UMAP done on PCA results from TF-IDF matrix.
lsi_workflow = function(bmat, dims, metadata=NULL, log_scale_tf=TRUE, reduction='pca.l2', resolution=0.3, n.neighbors=50, n.components=3, min.dist=0.01, k.param=20) {
  tfidf_mat = tfidf(bmat, log_scale_tf=log_scale_tf)
  pca_mat = do_pca(tfidf_mat, dims=max(dims))
  
  seurat_obj = run_dim_reduction(bmat, pca_mat, dims, metadata, reduction=reduction, n.neighbors=n.neighbors, n.components=n.components, min.dist=min.dist, k.param=k.param) %>%
    Seurat::FindClusters(n.start=20, resolution=resolution)
  return(seurat_obj)
}

# Wrapper for 10x version of full LSI workflow (TF-IDF and PCA + clustering + further dim reduction). Only TF-IDF step is modified.
# Args:
#   bmat (sparse matrix): sparse matrix (binarized)
#   dims (vector of int): vector of dims to use from cell_embeddings in downstream analysis
#   metadata: dataframe of metadata (rownames are cell names) to add to Seurat object
#   reduction (string): reduction to use for downstream steps. Can be 'pca' (cell_embeddings) or 'pca.l2' (L2 normalized cell_embeddings)
#   resolution (float): resolution parameter to Seurat Louvain clustering
# Returns:
#   Seurat object: Seurat object. clustering + tSNE + UMAP done on PCA results from TF-IDF matrix.
tenx_lsi_workflow = function(bmat, dims, metadata=NULL, reduction='pca.l2', resolution=0.3) {
  tfidf_mat = tenx_tfidf(bmat)
  pca_mat = do_pca(tfidf_mat, dims=max(dims))
  
  seurat_obj = run_dim_reduction(bmat, pca_mat, dims, metadata) %>%
    Seurat::FindClusters(reduction=reduction, n.start=20, resolution=resolution)
  return(seurat_obj)
}

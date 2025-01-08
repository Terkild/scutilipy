def write_10x_h5(adata, file):
    """Writes adata to a 10X-formatted h5 file.
    
    Note that this function is not fully tested and may not work for all cases.
    It will not write the following keys to the h5 file compared to 10X:
    '_all_tag_keys', 'pattern', 'read', 'sequence'

    Args:
        adata (AnnData object): AnnData object to be written.
        file (str): File name to be written to. If no extension is given, '.h5' is appended.

    Raises:
        FileExistsError: If file already exists.

    Returns:
        None
    """
    
    import numpy as np
    import pandas as pd
    import h5py
    
    if "feature_type" not in adata.var.columns:
        adata.var["feature_type"] = "Gene Expression"
    if "genome" not in adata.var.columns:
        adata.var["genome"] = "GRCh38"
    if "gene_ids" not in adata.var.columns:
        adata.var["gene_ids"] = adata.var.index
    
    if '.h5' not in file: file = f'{file}.h5'
    if Path(file).exists():
        raise FileExistsError(f"There already is a file `{file}`.")
    def int_max(x):
        return int(max(np.floor(len(str(int(max(x)))) / 4), 1) * 4)
    def str_max(x):
        return max([len(i) for i in x])

    w = h5py.File(file, 'w')
    grp = w.create_group("matrix")
    grp.create_dataset("barcodes", data=np.array(adata.obs_names, dtype=f'|S{str_max(adata.obs_names)}'))
    grp.create_dataset("data", data=np.array(adata.X.data, dtype=f'<i{int_max(adata.X.data)}'))
    ftrs = grp.create_group("features")
    # this group will lack the following keys:
    # '_all_tag_keys', 'feature_type', 'genome', 'id', 'name', 'pattern', 'read', 'sequence'
    ftrs.create_dataset("feature_type", data=np.array(adata.var.feature_type, dtype=f'|S{str_max(adata.var.feature_type)}'))
    ftrs.create_dataset("genome", data=np.array(adata.var.genome, dtype=f'|S{str_max(adata.var.genome)}'))
    ftrs.create_dataset("id", data=np.array(adata.var.gene_ids, dtype=f'|S{str_max(adata.var.gene_ids)}'))
    ftrs.create_dataset("name", data=np.array(adata.var.index, dtype=f'|S{str_max(adata.var.index)}'))
    grp.create_dataset("indices", data=np.array(adata.X.indices, dtype=f'<i{int_max(adata.X.indices)}'))
    grp.create_dataset("indptr", data=np.array(adata.X.indptr, dtype=f'<i{int_max(adata.X.indptr)}'))
    grp.create_dataset("shape", data=np.array(list(adata.X.shape)[::-1], dtype=f'<i{int_max(adata.X.shape)}'))

# ------------------------------------------------------------------
# 2) Define a function to export an obsm to HDF5
# ------------------------------------------------------------------
def export_obsm_to_hdf5(adata, obsm_key, filename, chunk_max=1024):
    """
    Save the obsm matrix (adata.obsm[obsm_key]) to an HDF5 file with:
      - "embedding" dataset (float32, chunked+gzipped)
      - "barcodes" dataset (cell names)
    """
    import numpy as np
    import pandas as pd
    import h5py
    
    embedding = adata.obsm[obsm_key]  # shape: (n_cells, n_dims)

    with h5py.File(filename, "w") as f:
        shape = embedding.shape
        # Define chunk sizes: (rows, cols). Tweak for performance if needed.
        chunk_rows = min(shape[0], chunk_max)
        chunk_cols = min(shape[1], chunk_max)
        
        dset = f.create_dataset(
            name="embedding",
            shape=shape,
            dtype="float32",
            chunks=(chunk_rows, chunk_cols),
            compression="gzip",
            compression_opts=4,
        )
        # If it fits in memory, write at once:
        dset[...] = embedding
        
        # Store barcodes (cell names). 
        # Convert strings to fixed-length format for HDF5
        f.create_dataset(
            "barcodes",
            data=adata.obs_names.values.astype("S")  # S = fixed-length string
        )
    print(f"Saved obsm['{obsm_key}'] to '{filename}'")
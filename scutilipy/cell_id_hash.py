def cell_id_hash(barcode_column, run_column, length=16):
    """
    Generate a truncated SHA-256 hash for each row of the provided columns.

    Parameters:
        barcode_column (pd.Series): The column containing barcodes.
        run_column (pd.Series): The column containing run identifiers.
        length (int): Length of the truncated hash. Default is 16.

    Returns:
        pd.Series: A Series containing the hash for each row.
    """
    import pandas as pd
    from hashlib import sha256
    
    # Ensure barcode and run are strings and truncate barcode to 16 characters
    truncated_barcodes = barcode_column.astype(str).str[:16]
    run_strings = run_column.astype(str)
    
    # Generate hash for each row
    hashes = truncated_barcodes + run_strings
    
    return hashes.apply(lambda x: sha256(x.encode()).hexdigest()[:length])
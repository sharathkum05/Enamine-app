"""
SPEI and PPEI metrics calculation functions for HTS Screening App
"""

import numpy as np
import pandas as pd


def calculate_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate SPEI and PPEI metrics.
    
    Adds columns:
        - Per_one: Pct_Inhibition / 100
        - SPEI: Per_one / (MW * 0.001)
        - PPEI: Per_one / (TPSA * 0.01)
    
    Args:
        df: DataFrame with Pct_Inhibition, MW, and TPSA columns
    
    Returns:
        DataFrame with added metrics columns
    """
    df = df.copy()
    
    # Convert %Inhibition to per-one scale (0-1)
    df['Per_one'] = df['Pct_Inhibition'] / 100.0
    
    # Calculate SPEI (Size-normalized Per-one Efficiency Index)
    # Higher values = more potent per unit mass
    if 'MW' in df.columns:
        df['SPEI'] = df.apply(
            lambda row: row['Per_one'] / (row['MW'] * 0.001) 
            if pd.notna(row['MW']) and row['MW'] > 0 else np.nan,
            axis=1
        )
    else:
        df['SPEI'] = np.nan
    
    # Calculate PPEI (Polarity-normalized Per-one Efficiency Index)
    # Higher values = better membrane penetration
    if 'TPSA' in df.columns:
        df['PPEI'] = df.apply(
            lambda row: row['Per_one'] / (row['TPSA'] * 0.01)
            if pd.notna(row['TPSA']) and row['TPSA'] > 0 else np.nan,
            axis=1
        )
    else:
        df['PPEI'] = np.nan
    
    return df


def get_top_candidates(df: pd.DataFrame, percentile: float, metric: str = "combined") -> pd.DataFrame:
    """
    Get top candidates based on threshold.
    
    Args:
        df: DataFrame with SPEI and PPEI
        percentile: Top X% to select (e.g., 5 for top 5%)
        metric: "combined" (SPEI+PPEI), "spei", or "ppei"
    
    Returns:
        DataFrame with Is_Top column marking top candidates
    """
    df = df.copy()
    
    # Calculate ranking score based on metric
    if metric == "combined":
        df['Rank_Score'] = df['SPEI'].fillna(0) + df['PPEI'].fillna(0)
    elif metric == "spei":
        df['Rank_Score'] = df['SPEI'].fillna(0)
    elif metric == "ppei":
        df['Rank_Score'] = df['PPEI'].fillna(0)
    else:
        raise ValueError(f"Unknown metric: {metric}")
    
    # Calculate threshold for top percentile
    threshold = np.percentile(df['Rank_Score'].dropna(), 100 - percentile)
    
    # Mark top candidates
    df['Is_Top'] = df['Rank_Score'] >= threshold
    
    return df


def aggregate_replicates(df: pd.DataFrame, method: str) -> pd.DataFrame:
    """
    Aggregate replicate data based on user selection.
    
    Args:
        df: DataFrame with Rep column (1 or 2)
        method: One of "average", "rep1", "rep2", "both"
    
    Returns:
        Aggregated DataFrame (for "both", keeps both with Rep column)
    """
    if method == "both":
        return df.copy()
    
    if method == "rep1":
        return df[df['Rep'] == 1].reset_index(drop=True)
    
    if method == "rep2":
        return df[df['Rep'] == 2].reset_index(drop=True)
    
    if method == "average":
        # Group by plate and well, average numeric columns
        group_cols = ['Plate', 'Well']
        
        # Add library columns if present
        library_cols = ['Plate_ID', 'catalog_number', 'Smiles', 'MW', 'TPSA', 'Chemical_name']
        for col in library_cols:
            if col in df.columns:
                group_cols.append(col)
        
        # Numeric columns to average
        numeric_cols = ['Signal', 'Pct_Inhibition']
        
        # Group and aggregate
        agg_dict = {}
        for col in numeric_cols:
            if col in df.columns:
                agg_dict[col] = 'mean'
        
        # For non-numeric columns, take first value
        result = df.groupby(['Plate', 'Well'], as_index=False).agg(
            {**agg_dict, **{col: 'first' for col in group_cols if col not in ['Plate', 'Well'] and col in df.columns}}
        )
        
        # Add Rep column indicating averaged
        result['Rep'] = 'avg'
        
        return result
    
    raise ValueError(f"Unknown method: {method}")


def calculate_summary_stats(df: pd.DataFrame, column: str) -> dict:
    """
    Calculate summary statistics for a column.
    
    Args:
        df: DataFrame
        column: Column name to summarize
    
    Returns:
        Dictionary with count, mean, median, std, min, max
    """
    data = df[column].dropna()
    
    return {
        'count': len(data),
        'mean': data.mean(),
        'median': data.median(),
        'std': data.std(),
        'min': data.min(),
        'max': data.max()
    }


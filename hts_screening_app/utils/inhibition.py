"""
%Inhibition and Z' calculation functions for HTS Screening App
"""

from string import ascii_uppercase
from typing import List, Dict, Tuple

import numpy as np
import pandas as pd


def expand_control_span(span: str) -> List[str]:
    """
    Expand a span like 'B02-H02' to list of well IDs.
    
    Example: "B02-H02" → ["B02", "C02", "D02", "E02", "F02", "G02", "H02"]
    
    Args:
        span: Control span string like "B02-H02"
    
    Returns:
        List of well IDs
    
    Raises:
        ValueError: If span format is invalid
    """
    try:
        left, right = span.split("-")
        row_start, col_start = left[0], int(left[1:])
        row_end, col_end = right[0], int(right[1:])
    except Exception:
        raise ValueError(f"Bad span format: {span}. Expected like 'B02-H02'.")
    
    if col_start != col_end:
        raise ValueError(f"Span must be same column (got {span})")
    
    rows = ascii_uppercase[
        ascii_uppercase.index(row_start):ascii_uppercase.index(row_end) + 1
    ]
    
    return [f"{r}{col_start:02d}" for r in rows]


def calculate_control_stats(df: pd.DataFrame, neg_wells: List[str], 
                           pos_wells: List[str]) -> Dict:
    """
    Calculate statistics for control wells.
    
    Args:
        df: DataFrame with Well and Signal columns
        neg_wells: List of negative control well IDs
        pos_wells: List of positive control well IDs
    
    Returns:
        Dictionary with neg_mean, neg_sd, neg_n, pos_mean, pos_sd, pos_n
    
    Raises:
        ValueError: If control wells are not found
    """
    neg_signals = df.loc[df["Well"].isin(neg_wells), "Signal"].astype(float)
    pos_signals = df.loc[df["Well"].isin(pos_wells), "Signal"].astype(float)
    
    if neg_signals.empty:
        raise ValueError("Negative control wells not found in plate")
    if pos_signals.empty:
        raise ValueError("Positive control wells not found in plate")
    
    return {
        'neg_mean': neg_signals.mean(),
        'neg_sd': neg_signals.std(ddof=1),
        'neg_n': len(neg_signals),
        'pos_mean': pos_signals.mean(),
        'pos_sd': pos_signals.std(ddof=1),
        'pos_n': len(pos_signals)
    }


def calculate_percent_inhibition(signal: float, neg_mean: float, pos_mean: float) -> float:
    """
    Calculate %Inhibition for a single well.
    
    Formula: 100 * (neg_mean - signal) / (neg_mean - pos_mean)
    
    Args:
        signal: Luminescence reading for the test well
        neg_mean: Average signal of negative control wells
        pos_mean: Average signal of positive control wells
    
    Returns:
        %Inhibition value
    """
    denom = neg_mean - pos_mean
    if denom == 0 or np.isclose(denom, 0):
        return np.nan
    return 100.0 * (neg_mean - signal) / denom


def calculate_z_prime(neg_mean: float, neg_sd: float, 
                      pos_mean: float, pos_sd: float) -> float:
    """
    Calculate Z' factor for plate quality.
    
    Formula: 1 - 3 * (neg_sd + pos_sd) / |neg_mean - pos_mean|
    
    Args:
        neg_mean: Mean of negative control signals
        neg_sd: Standard deviation of negative control signals
        pos_mean: Mean of positive control signals
        pos_sd: Standard deviation of positive control signals
    
    Returns:
        Z' factor (good quality if >= 0.5)
    """
    denom = abs(neg_mean - pos_mean)
    if denom == 0 or np.isclose(denom, 0):
        return np.nan
    return 1.0 - 3.0 * (neg_sd + pos_sd) / denom


def process_plate_inhibition(df: pd.DataFrame, 
                             neg_control_spans: List[str],
                             pos_control_spans: List[str]) -> Tuple[pd.DataFrame, Dict]:
    """
    Calculate %Inhibition for all wells in a plate and compute QC stats.
    
    Args:
        df: Melted plate data with Plate, Well, Rep, Signal columns
        neg_control_spans: List of negative control span strings
        pos_control_spans: List of positive control span strings
    
    Returns:
        Tuple of (DataFrame with Pct_Inhibition column, QC stats dict)
    """
    # Expand control spans
    neg_wells = []
    for span in neg_control_spans:
        neg_wells.extend(expand_control_span(span))
    
    pos_wells = []
    for span in pos_control_spans:
        pos_wells.extend(expand_control_span(span))
    
    # Calculate control statistics
    stats = calculate_control_stats(df, neg_wells, pos_wells)
    
    # Calculate Z' factor
    z_prime = calculate_z_prime(
        stats['neg_mean'], stats['neg_sd'],
        stats['pos_mean'], stats['pos_sd']
    )
    
    # Calculate %Inhibition for all wells
    df = df.copy()
    df['Pct_Inhibition'] = df['Signal'].apply(
        lambda x: calculate_percent_inhibition(x, stats['neg_mean'], stats['pos_mean'])
    )
    
    # Mark control wells
    df['Is_Control'] = df['Well'].isin(neg_wells + pos_wells)
    
    # Create QC summary
    plate = df['Plate'].iloc[0]
    rep = df['Rep'].iloc[0]
    qc_stats = {
        'Plate': plate,
        'Rep': rep,
        'Neg_Mean': stats['neg_mean'],
        'Neg_SD': stats['neg_sd'],
        'Neg_N': stats['neg_n'],
        'Pos_Mean': stats['pos_mean'],
        'Pos_SD': stats['pos_sd'],
        'Pos_N': stats['pos_n'],
        'Z_Prime': z_prime,
        'QC_Pass': z_prime >= 0.5 if not np.isnan(z_prime) else False
    }
    
    return df, qc_stats





"""
Plate file parsing functions for HTS Screening App
"""

import os
import re
from string import ascii_uppercase
from typing import Tuple, Optional
from io import BytesIO

import numpy as np
import pandas as pd


def parse_plate_filename(filename: str) -> Tuple[int, int]:
    """
    Extract plate number and replicate from filename.
    
    Args:
        filename: e.g., "301-1.xlsx"
    
    Returns:
        Tuple of (plate_number, replicate): e.g., (301, 1)
    
    Raises:
        ValueError: If filename doesn't match expected format
    """
    base = os.path.splitext(os.path.basename(filename))[0]
    match = re.match(r"^(\d+)-(\d+)$", base)
    if not match:
        raise ValueError(f"Filename '{filename}' doesn't match pattern '<plate>-<rep>.xlsx'")
    return int(match.group(1)), int(match.group(2))


def find_plate_block(df: pd.DataFrame) -> pd.DataFrame:
    """
    Locate the 384-well plate block in the Excel file.
    
    The data block starts at the row where column 0 contains "A".
    Returns a 16×24 DataFrame with rows A-P and columns 1-24.
    
    Args:
        df: Raw DataFrame from Excel file
    
    Returns:
        16x24 DataFrame with row labels A-P and column labels 1-24
    
    Raises:
        ValueError: If plate block cannot be found or has wrong dimensions
    """
    start_idx = None
    
    # Find the row where column 0 contains "A"
    for i, v in df.iloc[:, 0].items():
        if isinstance(v, str) and v.strip().upper() == "A":
            start_idx = i
            break
    
    if start_idx is None:
        raise ValueError("Could not find row label 'A' in column 0.")
    
    # Extract the 16x24 block
    block = df.iloc[start_idx:start_idx + 16, 1:25].copy()
    
    if block.shape != (16, 24):
        raise ValueError(f"Detected block is {block.shape}, expected (16, 24).")
    
    # Set proper row and column labels
    block.index = list(ascii_uppercase[:16])  # A-P
    block.columns = list(range(1, 25))         # 1-24
    
    return block


def melt_plate_block(block: pd.DataFrame, plate: int, rep: int) -> pd.DataFrame:
    """
    Convert plate matrix to tidy long format.
    
    Args:
        block: 16x24 plate matrix DataFrame
        plate: Plate number
        rep: Replicate number
    
    Returns:
        DataFrame with columns: Plate, Well, Rep, Signal
    """
    records = []
    for row in block.index:
        for col in block.columns:
            well = f"{row}{col:02d}"
            value = block.loc[row, col]
            signal = float(value) if pd.notna(value) else np.nan
            records.append({
                'Plate': plate,
                'Well': well,
                'Rep': rep,
                'Signal': signal
            })
    
    return pd.DataFrame(records)


def filter_perimeter_wells(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove perimeter wells (row A, P and column 1, 24).
    Keep only rows B-O and columns 02-23.
    
    Args:
        df: DataFrame with Well column
    
    Returns:
        Filtered DataFrame with perimeter wells removed
    """
    valid_rows = set("BCDEFGHIJKLMNO")  # B-O
    valid_cols = {f"{i:02d}" for i in range(2, 24)}  # 02-23
    
    rows = df["Well"].str[0]
    cols = df["Well"].str[1:]
    
    keep = rows.isin(valid_rows) & cols.isin(valid_cols)
    
    return df.loc[keep].reset_index(drop=True)


def process_plate_file(file_content: BytesIO, filename: str) -> Tuple[pd.DataFrame, int, int]:
    """
    Process a single plate file and return the melted data.
    
    Args:
        file_content: BytesIO object containing the Excel file
        filename: Original filename for parsing plate/rep info
    
    Returns:
        Tuple of (melted_data, plate_number, replicate)
    """
    # Parse filename
    plate, rep = parse_plate_filename(filename)
    
    # Read Excel file
    df = pd.read_excel(file_content, header=None)
    
    # Find and extract plate block
    block = find_plate_block(df)
    
    # Melt to long format
    melted = melt_plate_block(block, plate, rep)
    
    # Remove perimeter wells
    filtered = filter_perimeter_wells(melted)
    
    return filtered, plate, rep


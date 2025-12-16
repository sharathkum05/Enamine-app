"""
Library merging functions for HTS Screening App
"""

import sys
from pathlib import Path

import pandas as pd
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import ENAMINE_PLATE_PREFIX, TEST_COMPOUND_COLS


def get_enamine_plate_id(plate_number: int) -> str:
    """
    Convert plate number to ENAMINE Plate_ID format.
    
    Args:
        plate_number: Plate number (e.g., 301)
    
    Returns:
        Plate_ID string (e.g., "2096462-Y10-301")
    """
    return f"{ENAMINE_PLATE_PREFIX}{plate_number:03d}"


def merge_with_library(plate_data: pd.DataFrame, library: pd.DataFrame) -> pd.DataFrame:
    """
    Merge plate data with ENAMINE library on Plate_ID and Well.
    
    Args:
        plate_data: DataFrame with Plate, Well, Rep, Signal, Pct_Inhibition columns
        library: ENAMINE library DataFrame with Plate_ID, Well, and compound info
    
    Returns:
        Merged DataFrame with compound information
    """
    # Create Plate_ID from plate number
    plate_data = plate_data.copy()
    plate_data['Plate_ID'] = plate_data['Plate'].apply(get_enamine_plate_id)
    
    # Standardize well format in library if needed
    library = library.copy()
    
    # Ensure Well column exists and is properly formatted
    if 'Well' not in library.columns:
        raise ValueError("Library must have 'Well' column")
    
    # Standardize well format (e.g., "B3" → "B03")
    def standardize_well(well):
        if pd.isna(well):
            return well
        well = str(well).strip()
        if len(well) >= 2:
            row = well[0].upper()
            col = well[1:]
            try:
                col_num = int(col)
                return f"{row}{col_num:02d}"
            except ValueError:
                return well
        return well
    
    library['Well'] = library['Well'].apply(standardize_well)
    
    # Select relevant columns from library
    library_cols = ['Plate_ID', 'Well']
    
    # Add available compound info columns
    optional_cols = ['Smiles', 'MW', 'TPSA', 'catalog number', 'Chemical name']
    for col in optional_cols:
        if col in library.columns:
            library_cols.append(col)
    
    library_subset = library[library_cols].drop_duplicates()
    
    # Merge
    merged = pd.merge(
        plate_data,
        library_subset,
        on=['Plate_ID', 'Well'],
        how='left'
    )
    
    # Rename columns for consistency
    column_map = {
        'catalog number': 'catalog_number',
        'Chemical name': 'Chemical_name'
    }
    merged = merged.rename(columns=column_map)
    
    return merged


def filter_test_compounds(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter to keep only test compound wells (exclude controls).
    Test compounds are in columns 03-22, rows B-O.
    
    Args:
        df: Merged DataFrame
    
    Returns:
        Filtered DataFrame with only test compounds
    """
    df = df.copy()
    
    # Filter out control wells (columns 02 and 23)
    cols = df['Well'].str[1:]
    is_test_compound = cols.isin(TEST_COMPOUND_COLS)
    
    # Also exclude if marked as control
    if 'Is_Control' in df.columns:
        is_test_compound = is_test_compound & (~df['Is_Control'])
    
    return df[is_test_compound].reset_index(drop=True)


def load_library(file_content) -> pd.DataFrame:
    """
    Load the ENAMINE library file.
    
    Args:
        file_content: File-like object or path to Excel file
    
    Returns:
        Library DataFrame
    """
    library = pd.read_excel(file_content)
    
    # Standardize column names
    library.columns = library.columns.str.strip()
    
    return library


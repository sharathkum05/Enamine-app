"""
Configuration constants for ITR ENAMINE LIBRARY Screening App
"""

from string import ascii_uppercase

# Control well configurations
NEG_CONTROL_SPANS = ["B02-H02", "B23-O23"]  # 0% inhibition reference (bacteria grow normally)
POS_CONTROL_SPANS = ["I02-O02"]             # 100% inhibition reference (bacteria killed)

# Plate dimensions
PLATE_ROWS = list(ascii_uppercase[:16])  # A-P
PLATE_COLS = list(range(1, 25))           # 1-24

# Valid wells (excluding perimeter)
VALID_ROWS = list("BCDEFGHIJKLMNO")       # B-O (exclude A and P)
VALID_COLS = [f"{i:02d}" for i in range(2, 24)]  # 02-23 (exclude 01 and 24)

# Test compound columns (excluding controls)
TEST_COMPOUND_COLS = [f"{i:02d}" for i in range(3, 23)]  # 03-22

# QC thresholds
Z_PRIME_THRESHOLD = 0.5  # Minimum Z' for good quality plate

# ENAMINE library plate ID format
ENAMINE_PLATE_PREFIX = "2096462-Y10-"

# Plot colors
COLOR_REP1 = "#1f77b4"  # Blue
COLOR_REP2 = "#ff7f0e"  # Orange
COLOR_HIGHLIGHT = "#e74c3c"  # Red for top candidates
COLOR_NORMAL = "#95a5a6"  # Gray for normal points

# App styling
APP_TITLE = "🧬 ITR ENAMINE LIBRARY Screening Data Analysis"
APP_ICON = "🧬"


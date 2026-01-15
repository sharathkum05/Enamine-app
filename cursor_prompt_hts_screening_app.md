# Cursor Prompt: HTS Screening Data Analysis Streamlit App

## Project Overview

Build a Streamlit application for analyzing High-Throughput Screening (HTS) data from ENAMINE compound library plates. The app will process 384-well plate luminescence data, calculate %Inhibition, merge with compound library information, and generate analysis plots for drug discovery research.

---

## Context: Drug Discovery Screening

This app supports antibiotic discovery research similar to the Stokes 2020 Halicin paper. Scientists screen thousands of compounds against bacteria, measuring luminescence (bacterial growth). Lower luminescence = more bacterial death = better antibiotic candidate.

**Key Metrics:**
- **%Inhibition**: How much a compound inhibits bacterial growth (0% = no effect, 100% = complete kill)
- **SPEI (Size-normalized Per-one Efficiency Index)**: Activity normalized by molecular weight - smaller molecules with same activity are better
- **PPEI (Polarity-normalized Per-one Efficiency Index)**: Activity normalized by polar surface area - less polar molecules penetrate membranes better

---

## Input Files

### 1. Plate Data Files (Multiple .xlsx files)
- **Naming convention**: `{PlateNumber}-{Replicate}.xlsx` (e.g., `301-1.xlsx`, `301-2.xlsx`)
- **Format**: CLARIOstar plate reader output, 384-well format
- **Structure**: Raw luminescence values in a 16×24 matrix (rows A-P, columns 1-24)
- **Data location**: Matrix starts at the row where column 0 contains "A"
- **Batch size**: ~116 plates × 2 replicates = ~232 files per batch

### 2. ENAMINE Library File (Enamine_library.xlsx)
- **Size**: ~100,160 compounds across 358 plates
- **Key columns**:
  - `Plate_ID`: Format "2096462-Y10-XXX" where XXX is the plate number
  - `Well`: Position like "B03", "C04"
  - `Smiles`: SMILES molecular structure string
  - `MW`: Molecular Weight (pre-calculated)
  - `TPSA`: Topological Polar Surface Area (pre-calculated)
  - `catalog number`: ENAMINE compound ID (e.g., "Z27358882")
  - `Chemical name`: Full compound name

### 3. Plate Layout Reference (Enamine_plate_layout1.xlsx)
- Shows control well positions (for reference only, logic is hardcoded)

---

## Plate Layout (384-well)

```
        1    2    3    4   ...  21   22   23   24
      ┌────┬────┬────┬────┬────┬────┬────┬────┬────┐
   A  │ ❌ │ ❌ │ ❌ │ ❌ │    │ ❌ │ ❌ │ ❌ │ ❌ │  ← Perimeter (excluded)
      ├────┼────┼────┼────┼────┼────┼────┼────┼────┤
   B  │ ❌ │ N  │ ✅ │ ✅ │    │ ✅ │ ✅ │ N  │ ❌ │  
      ├────┼────┼────┼────┼────┼────┼────┼────┼────┤
  ... │    │ N  │ ✅ │ ✅ │    │ ✅ │ ✅ │ N  │    │
      ├────┼────┼────┼────┼────┼────┼────┼────┼────┤
   H  │ ❌ │ N  │ ✅ │ ✅ │    │ ✅ │ ✅ │ N  │ ❌ │
      ├────┼────┼────┼────┼────┼────┼────┼────┼────┤
   I  │ ❌ │ P  │ ✅ │ ✅ │    │ ✅ │ ✅ │ N  │ ❌ │
      ├────┼────┼────┼────┼────┼────┼────┼────┼────┤
  ... │    │ P  │ ✅ │ ✅ │    │ ✅ │ ✅ │ N  │    │
      ├────┼────┼────┼────┼────┼────┼────┼────┼────┤
   O  │ ❌ │ P  │ ✅ │ ✅ │    │ ✅ │ ✅ │ N  │ ❌ │
      ├────┼────┼────┼────┼────┼────┼────┼────┼────┤
   P  │ ❌ │ ❌ │ ❌ │ ❌ │    │ ❌ │ ❌ │ ❌ │ ❌ │  ← Perimeter (excluded)
      └────┴────┴────┴────┴────┴────┴────┴────┴────┘

Legend:
  ❌ = Perimeter wells (excluded from analysis)
  N  = Negative control (0% inhibition) - B02-H02 and B23-O23
  P  = Positive control (100% inhibition) - I02-O02
  ✅ = Test compounds (columns 03-22, rows B-O)
```

---

## Core Formulas

### %Inhibition Calculation
```python
%Inhibition = 100 × (Neg_Mean - Signal) / (Neg_Mean - Pos_Mean)

# Where:
# - Neg_Mean = average signal of negative control wells (bacteria grow normally)
# - Pos_Mean = average signal of positive control wells (bacteria killed)
# - Signal = luminescence reading for the test well
```

### Z' Factor (Quality Control)
```python
Z_prime = 1 - 3 × (Neg_SD + Pos_SD) / |Neg_Mean - Pos_Mean|

# Z' ≥ 0.5 = Good quality plate
# Z' < 0.5 = Poor quality plate (may need to exclude)
```

### SPEI and PPEI
```python
Per_one = %Inhibition / 100  # Convert to 0-1 scale

SPEI = Per_one / (MW × 0.001)    # Size efficiency (MW in kDa)
PPEI = Per_one / (TPSA × 0.01)   # Polarity efficiency
```

---

## Data Model

### Backend Storage (Keeps Replicates Separate)
```python
# Primary DataFrame schema
main_df = pd.DataFrame({
    'Plate': int,              # 301, 302, etc. (extracted from filename)
    'Well': str,               # B03, B04, etc.
    'Rep': int,                # 1 or 2 (replicate number)
    'Signal': float,           # Raw luminescence value
    'Pct_Inhibition': float,   # Calculated %Inhibition
    
    # Merged from ENAMINE library
    'Plate_ID': str,           # 2096462-Y10-301
    'catalog_number': str,     # Z27358882
    'Smiles': str,             # SMILES string
    'MW': float,               # Molecular Weight
    'TPSA': float,             # Polar Surface Area
    'Chemical_name': str,      # Full compound name
})

# QC Summary DataFrame
qc_df = pd.DataFrame({
    'Plate': int,
    'Rep': int,
    'Neg_Mean': float,
    'Neg_SD': float,
    'Pos_Mean': float,
    'Pos_SD': float,
    'Z_Prime': float,
    'QC_Pass': bool,           # Z' >= 0.5
})

# Analysis DataFrame (after replicate selection)
analysis_df = pd.DataFrame({
    'Plate': int,
    'Well': str,
    'Pct_Inhibition': float,   # Based on user selection (avg/rep1/rep2)
    'Per_one': float,          # Pct_Inhibition / 100
    'MW': float,
    'TPSA': float,
    'SPEI': float,             # per_one / (MW * 0.001)
    'PPEI': float,             # per_one / (TPSA * 0.01)
    'catalog_number': str,
    'Smiles': str,
    'Chemical_name': str,
})
```

---

## Plate ID Mapping Logic

```python
def get_enamine_plate_id(plate_number: int) -> str:
    """Convert plate number to ENAMINE Plate_ID format."""
    return f"2096462-Y10-{plate_number:03d}"

# Example:
# Filename: 301-1.xlsx
# → Plate = 301, Rep = 1
# → Plate_ID = "2096462-Y10-301"
# → Merge with library on Plate_ID + Well
```

---

## App Structure (Step-by-Step Flow)

### Step 1: File Upload
```
┌─────────────────────────────────────────────────────────────────┐
│  UPLOAD FILES                                                   │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  📁 ENAMINE Library File (required, one-time):                 │
│     [Upload Enamine_library.xlsx]                              │
│                                                                 │
│  📁 Plate Data Files (multiple .xlsx):                         │
│     [Upload 301-1.xlsx, 301-2.xlsx, 302-1.xlsx, ...]          │
│                                                                 │
│  ─── OR ───                                                    │
│                                                                 │
│  📂 Load Previously Processed Data:                            │
│     [Load processed_data.pkl]                                  │
│                                                                 │
│  [Process Files]  ← Button with progress bar                   │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Step 2: QC Summary
```
┌─────────────────────────────────────────────────────────────────┐
│  QUALITY CONTROL SUMMARY                                        │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  ┌───────┬─────┬──────────┬─────────┬─────────┬───────┬──────┐ │
│  │ Plate │ Rep │ Neg_Mean │ Pos_Mean│ Z'      │ Status│ Use  │ │
│  ├───────┼─────┼──────────┼─────────┼─────────┼───────┼──────┤ │
│  │ 301   │ 1   │ 65432.1  │ 543.2   │ 0.72    │ ✅    │ ☑️   │ │
│  │ 301   │ 2   │ 64123.4  │ 521.1   │ 0.74    │ ✅    │ ☑️   │ │
│  │ 302   │ 1   │ 61234.5  │ 612.3   │ 0.38    │ ⚠️    │ ☐   │ │
│  │ ...   │ ... │ ...      │ ...     │ ...     │       │      │ │
│  └───────┴─────┴──────────┴─────────┴─────────┴───────┴──────┘ │
│                                                                 │
│  Summary: 230/232 plates passed QC (Z' ≥ 0.5)                  │
│                                                                 │
│  [Save Processed Data]  ← Saves to .pkl for future use         │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Step 3: Replicate Selection
```
┌─────────────────────────────────────────────────────────────────┐
│  REPLICATE HANDLING                                             │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  How to handle replicates for analysis?                        │
│                                                                 │
│  ◉ Average (mean of Rep1 & Rep2)     ← Default                 │
│  ○ Rep1 only                                                   │
│  ○ Rep2 only                                                   │
│  ○ Show both (overlay on plots)                                │
│                                                                 │
│  [Apply Selection]                                             │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Step 4: Task 1A - %Inhibition Histogram
```
┌─────────────────────────────────────────────────────────────────┐
│  TASK 1A: %INHIBITION DISTRIBUTION                              │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  [Interactive Histogram - Plotly]                              │
│                                                                 │
│  Statistics:                                                   │
│  • Count: 32,480 compounds                                     │
│  • Mean: 45.2%  |  Median: 42.1%  |  SD: 28.3%                │
│  • Min: -15.2%  |  Max: 112.5%                                 │
│                                                                 │
│  [Download Plot PNG] [Download Data CSV]                       │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Step 5: Task 1B - Normalized Histogram
```
┌─────────────────────────────────────────────────────────────────┐
│  TASK 1B: SIZE-NORMALIZED ACTIVITY                              │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  Formula: (% Inhibition / 100) / (MW × 0.001)                  │
│                                                                 │
│  [Interactive Histogram - Plotly]                              │
│                                                                 │
│  Interpretation: Higher values = more potent per unit mass     │
│                                                                 │
│  [Download Plot PNG] [Download Data CSV]                       │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Step 6: SPEI & PPEI Calculation
```
┌─────────────────────────────────────────────────────────────────┐
│  EFFICIENCY METRICS                                             │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  SPEI = per_one / (MW × 0.001)   [Size efficiency]             │
│  PPEI = per_one / (TPSA × 0.01)  [Polarity efficiency]         │
│                                                                 │
│  Data Preview:                                                 │
│  ┌────────────────┬─────────────┬─────────┬────────┬───────┐   │
│  │ catalog_number │ %Inhibition │ MW      │ SPEI   │ PPEI  │   │
│  ├────────────────┼─────────────┼─────────┼────────┼───────┤   │
│  │ Z27358882      │ 78.5        │ 407.51  │ 1.93   │ 0.72  │   │
│  │ Z133956054     │ 65.2        │ 342.46  │ 1.90   │ 1.30  │   │
│  │ ...            │ ...         │ ...     │ ...    │ ...   │   │
│  └────────────────┴─────────────┴─────────┴────────┴───────┘   │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Step 7: Cartesian Plot (PPEI vs SPEI)
```
┌─────────────────────────────────────────────────────────────────┐
│  CARTESIAN PLOT: BEST CANDIDATES                                │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  Highlight Top Candidates:                                     │
│  [────────●────────] 5%   (slider: 1% to 20%)                  │
│                                                                 │
│  Metric for ranking: ◉ SPEI + PPEI  ○ SPEI only  ○ PPEI only  │
│                                                                 │
│  [Interactive Scatter Plot - Plotly]                           │
│  • X-axis: PPEI (Polarity efficiency)                          │
│  • Y-axis: SPEI (Size efficiency)                              │
│  • Northeast corner = Best candidates                          │
│  • Hover: Shows compound details                               │
│  • Highlighted points: Top candidates based on threshold       │
│                                                                 │
│  [Download Plot PNG] [Download Top Candidates CSV]             │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Step 8: Export
```
┌─────────────────────────────────────────────────────────────────┐
│  EXPORT RESULTS                                                 │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  [Download Complete Dataset (Excel)]                           │
│  [Download Top Candidates Report (Excel)]                      │
│  [Download All Plots (ZIP)]                                    │
│  [Save Session Data (PKL)] ← For future loading                │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## Key Functions to Implement

### 1. File Processing Functions

```python
def parse_plate_filename(filename: str) -> Tuple[int, int]:
    """
    Extract plate number and replicate from filename.
    
    Args:
        filename: e.g., "301-1.xlsx"
    
    Returns:
        Tuple of (plate_number, replicate): e.g., (301, 1)
    """
    pass

def find_plate_block(df: pd.DataFrame) -> pd.DataFrame:
    """
    Locate the 384-well plate block in the Excel file.
    
    The data block starts at the row where column 0 contains "A".
    Returns a 16×24 DataFrame with rows A-P and columns 1-24.
    """
    pass

def melt_plate_block(block: pd.DataFrame, plate: int, rep: int) -> pd.DataFrame:
    """
    Convert plate matrix to tidy long format.
    
    Returns DataFrame with columns: Plate, Well, Rep, Signal
    """
    pass

def filter_perimeter_wells(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove perimeter wells (row A, P and column 1, 24).
    Keep only rows B-O and columns 02-23.
    """
    pass
```

### 2. %Inhibition Calculation Functions

```python
def expand_control_span(span: str) -> List[str]:
    """
    Expand a span like 'B02-H02' to list of well IDs.
    
    Example: "B02-H02" → ["B02", "C02", "D02", "E02", "F02", "G02", "H02"]
    """
    pass

def calculate_control_stats(df: pd.DataFrame, neg_wells: List[str], pos_wells: List[str]) -> Dict:
    """
    Calculate statistics for control wells.
    
    Returns: {
        'neg_mean': float, 'neg_sd': float, 'neg_n': int,
        'pos_mean': float, 'pos_sd': float, 'pos_n': int
    }
    """
    pass

def calculate_percent_inhibition(signal: float, neg_mean: float, pos_mean: float) -> float:
    """
    Calculate %Inhibition for a single well.
    
    Formula: 100 * (neg_mean - signal) / (neg_mean - pos_mean)
    """
    pass

def calculate_z_prime(neg_mean: float, neg_sd: float, pos_mean: float, pos_sd: float) -> float:
    """
    Calculate Z' factor for plate quality.
    
    Formula: 1 - 3 * (neg_sd + pos_sd) / |neg_mean - pos_mean|
    """
    pass
```

### 3. Data Merging Functions

```python
def get_enamine_plate_id(plate_number: int) -> str:
    """Convert plate number to ENAMINE Plate_ID format."""
    return f"2096462-Y10-{plate_number:03d}"

def merge_with_library(plate_data: pd.DataFrame, library: pd.DataFrame) -> pd.DataFrame:
    """
    Merge plate data with ENAMINE library on Plate_ID and Well.
    
    Adds columns: Smiles, MW, TPSA, catalog_number, Chemical_name
    """
    pass
```

### 4. Replicate Handling Functions

```python
def aggregate_replicates(df: pd.DataFrame, method: str) -> pd.DataFrame:
    """
    Aggregate replicate data based on user selection.
    
    Args:
        df: DataFrame with Rep column (1 or 2)
        method: One of "average", "rep1", "rep2", "both"
    
    Returns:
        Aggregated DataFrame (for "both", keeps both with Rep column)
    """
    pass
```

### 5. Metrics Calculation Functions

```python
def calculate_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate SPEI and PPEI metrics.
    
    Adds columns:
        - Per_one: Pct_Inhibition / 100
        - SPEI: Per_one / (MW * 0.001)
        - PPEI: Per_one / (TPSA * 0.01)
    """
    pass

def get_top_candidates(df: pd.DataFrame, percentile: float, metric: str) -> pd.DataFrame:
    """
    Get top candidates based on threshold.
    
    Args:
        df: DataFrame with SPEI and PPEI
        percentile: Top X% to select (e.g., 5 for top 5%)
        metric: "combined" (SPEI+PPEI), "spei", or "ppei"
    
    Returns:
        DataFrame with top candidates marked
    """
    pass
```

### 6. Plotting Functions

```python
def plot_inhibition_histogram(df: pd.DataFrame, replicate_mode: str) -> go.Figure:
    """
    Create histogram of %Inhibition.
    
    For "both" mode: overlay semi-transparent histograms for Rep1 and Rep2
    """
    pass

def plot_normalized_histogram(df: pd.DataFrame, replicate_mode: str) -> go.Figure:
    """
    Create histogram of Per_one / MW(kDa).
    
    This is equivalent to SPEI.
    """
    pass

def plot_cartesian(df: pd.DataFrame, top_percentile: float, 
                   replicate_mode: str, metric: str) -> go.Figure:
    """
    Create PPEI vs SPEI scatter plot.
    
    Features:
        - X-axis: PPEI, Y-axis: SPEI
        - Highlight top candidates
        - For "both" mode: different colors for Rep1 vs Rep2
        - Hover shows compound details
    """
    pass
```

---

## Control Well Configuration

```python
# Negative controls (0% inhibition - bacteria grow normally)
NEG_CONTROL_SPANS = ["B02-H02", "B23-O23"]
# Expands to: B02, C02, D02, E02, F02, G02, H02, B23, C23, D23, E23, F23, G23, H23, I23, J23, K23, L23, M23, N23, O23

# Positive controls (100% inhibition - bacteria killed)
POS_CONTROL_SPANS = ["I02-O02"]
# Expands to: I02, J02, K02, L02, M02, N02, O02
```

---

## Session State Management

```python
# Streamlit session state keys
st.session_state.library_df        # ENAMINE library DataFrame
st.session_state.main_df           # All processed plate data (replicates separate)
st.session_state.qc_df             # QC summary
st.session_state.analysis_df       # Current analysis DataFrame (after replicate selection)
st.session_state.replicate_mode    # "average", "rep1", "rep2", or "both"
st.session_state.processing_done   # Boolean flag
```

---

## File Structure

```
hts_screening_app/
├── app.py                    # Main Streamlit application
├── utils/
│   ├── __init__.py
│   ├── plate_parser.py       # Plate file parsing functions
│   ├── inhibition.py         # %Inhibition and Z' calculations
│   ├── data_merger.py        # Library merging functions
│   ├── metrics.py            # SPEI, PPEI calculations
│   └── plotting.py           # All plotting functions
├── config.py                 # Configuration constants
├── requirements.txt          # Dependencies
└── README.md                 # Documentation
```

---

## Requirements

```txt
streamlit>=1.28.0
pandas>=2.0.0
numpy>=1.24.0
plotly>=5.18.0
openpyxl>=3.1.0
xlrd>=2.0.0
```

---

## Important Implementation Notes

1. **Progress Bar**: Use `st.progress()` with status updates when processing 232 files

2. **Caching**: Use `@st.cache_data` for library loading and expensive computations

3. **Error Handling**: Gracefully handle:
   - Missing wells in library
   - Malformed plate files
   - Division by zero in calculations

4. **Memory Efficiency**: Process plates in batches if memory is a concern

5. **Replicate "Show Both" Mode**:
   - Histograms: Overlay with 50% opacity, different colors (blue for Rep1, orange for Rep2)
   - Cartesian: Different colored points, same legend

6. **Hover Information on Cartesian Plot**:
   ```
   Compound: Z27358882
   Name: 2-[4-(dimethylsulfamoyl)benzamido]...
   %Inhibition: 78.5%
   SPEI: 1.93
   PPEI: 0.72
   MW: 407.51
   TPSA: 109.57
   ```

7. **Top Candidates Highlighting**:
   - Use larger marker size
   - Different color (e.g., red/gold)
   - Add to legend as "Top X%"

---

## Reference: Existing %Inhibition Script

The existing `hts_percent_inhibition_zprime.py` script contains tested logic for:
- Finding plate blocks in Excel files
- Expanding control well spans
- Filtering perimeter wells
- Calculating %Inhibition and Z'
- Handling replicates

**Reuse this logic** but adapt for Streamlit integration and the new data model.

---

## Testing Checklist

- [ ] Single plate file upload works
- [ ] Multiple plate files with progress bar
- [ ] Library merging produces correct compound info
- [ ] %Inhibition calculation matches existing script
- [ ] Z' calculation is correct
- [ ] All 4 replicate modes work
- [ ] "Show both" overlays histograms correctly
- [ ] Cartesian plot interactive hover works
- [ ] Top candidates slider updates plot
- [ ] All download buttons work
- [ ] Save/Load session data works
- [ ] Error handling for edge cases

---

## Example Usage Flow

1. User uploads `Enamine_library.xlsx`
2. User uploads 232 plate files (or loads saved data)
3. App shows progress bar during processing
4. QC summary displayed - user can exclude bad plates
5. User selects replicate handling method
6. Task 1A histogram generated
7. Task 1B normalized histogram generated  
8. SPEI/PPEI calculated and displayed
9. Cartesian plot with slider for top candidates
10. User downloads results and saves session

---

## Additional Context

This app is part of a drug discovery pipeline similar to the Stokes 2020 Cell paper that discovered Halicin using deep learning. The screening data quality assessment (this app) happens BEFORE AI modeling. Gary Liu from J. Stokes' lab will use the output for AI predictions.

The "best candidates" in the northeast corner of the Cartesian plot are compounds that:
- Are small (high SPEI = good activity per unit mass)
- Are non-polar (high PPEI = good membrane penetration)
- Have high %Inhibition

These properties correlate with better drug-likeness and bioavailability.

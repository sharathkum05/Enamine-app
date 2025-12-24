# 🧬 ITR ENAMINE LIBRARY Screening Data Analysis App

A comprehensive Streamlit application for analyzing High-Throughput Screening data from ITR ENAMINE compound library plates. Built for antibiotic discovery research.

## Overview

This app processes 384-well plate luminescence data, calculates %Inhibition, merges with compound library information, and generates analysis plots for drug discovery research.

### Key Features

- **User-Uploaded Library**: Upload ENAMINE library file each session (not stored on server)
- **Batch Processing**: Process hundreds of plate files with progress tracking
- **Quality Control**: Z' factor calculation with visual QC summaries
- **Replicate Handling**: Average, single replicate, or overlay both
- **Plate Selection**: Analyze individual plates or all plates combined
- **%Inhibition Calculation**: Automatic control well identification and normalization
- **Drug Efficiency Metrics**: SPEI and PPEI calculation with histograms
- **Interactive Plots**: Plotly-powered histograms and scatter plots
- **Top Candidate Identification**: Highlight best drug candidates
- **Session Persistence**: Save and load processed data

## Setup Instructions

### Required Data Files

⚠️ **Important:** This app requires the ENAMINE library file which is **NOT included** in the repository due to size/licensing restrictions.

**Required files to upload each session:**
1. `Enamine_library.xlsx` - ENAMINE compound library with SMILES, MW, TPSA data (~100MB+)

**Why upload every time?**
- Files are excluded from GitHub due to large size (100MB+)
- Data is not stored on the server for security/privacy
- Each session requires fresh upload

### Required Columns in ENAMINE Library

The ENAMINE library Excel file must contain these columns:
- `Plate_ID`: Plate identifier (e.g., "2096462-Y10-001")
- `Well`: Well position (e.g., "B03", "C04")
- `Smiles`: SMILES molecular structure string
- `MW`: Molecular weight
- `TPSA`: Topological polar surface area
- `catalog number`: Compound identifier (e.g., "Z27358882")
- `Chemical name`: Full compound name (optional)

### Installation Steps

1. **Clone this repository**
   ```bash
   git clone https://github.com/sharathkum05/Enamine-app.git
   cd Enamine-app/hts_screening_app
   ```

2. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

3. **Run the app**
   ```bash
   streamlit run app.py
   ```

4. **Prepare your ENAMINE library file**
   - Ensure you have `Enamine_library.xlsx` ready to upload
   - File should be accessible on your local machine

## Usage

1. **Start the app:**
   ```bash
   streamlit run app.py
   ```
   The app will open in your browser at `http://localhost:8501`

2. **Upload ENAMINE Library:**
   - Click "Upload ENAMINE Library (xlsx)" in Step 1
   - Select your `Enamine_library.xlsx` file
   - Wait for validation (checks for required columns)
   - ⚠️ **Note:** You'll need to upload this file every time you start a new session

3. **Upload Plate Data:**
   - After library is loaded, upload plate reader output files in Step 2
   - Files like `301-1.xlsx`, `301-2.xlsx`, etc.
   - Multiple files can be uploaded at once

4. **Process and Analyze:**
   - Review QC summary (Z' factors, control statistics)
   - Select replicate handling method (average, rep1, rep2, both)
   - Choose individual plates or analyze all plates
   - View histograms (%Inhibition, SPEI, PPEI)
   - Explore efficiency metrics
   - Identify top candidates in Cartesian plot
   - Export results

## File Structure

```
hts_screening_app/
├── app.py                    # Main Streamlit application
├── config.py                 # Configuration constants
├── requirements.txt          # Dependencies
├── README.md                 # This file
└── utils/
    ├── __init__.py
    ├── plate_parser.py       # Plate file parsing
    ├── inhibition.py         # %Inhibition calculations
    ├── data_merger.py        # Library merging
    ├── metrics.py            # SPEI/PPEI calculations
    └── plotting.py           # Visualization functions
```

## Input Files

### Plate Data Files
- **Format**: CLARIOstar plate reader output (384-well)
- **Naming**: `{PlateNumber}-{Replicate}.xlsx` (e.g., `301-1.xlsx`)
- **Structure**: 16×24 matrix with rows A-P, columns 1-24

### ENAMINE Library File
- **File**: `Enamine_library.xlsx`
- **Key Columns**: `Plate_ID`, `Well`, `Smiles`, `MW`, `TPSA`, `catalog number`

## Plate Layout

```
        1    2    3    4   ...  21   22   23   24
      ┌────┬────┬────┬────┬────┬────┬────┬────┬────┐
   A  │ ❌ │ ❌ │ ❌ │ ❌ │    │ ❌ │ ❌ │ ❌ │ ❌ │  ← Perimeter
   B  │ ❌ │ N  │ ✅ │ ✅ │    │ ✅ │ ✅ │ N  │ ❌ │  
  ... │    │    │    │    │    │    │    │    │    │
   H  │ ❌ │ N  │ ✅ │ ✅ │    │ ✅ │ ✅ │ N  │ ❌ │
   I  │ ❌ │ P  │ ✅ │ ✅ │    │ ✅ │ ✅ │ N  │ ❌ │
  ... │    │    │    │    │    │    │    │    │    │
   O  │ ❌ │ P  │ ✅ │ ✅ │    │ ✅ │ ✅ │ N  │ ❌ │
   P  │ ❌ │ ❌ │ ❌ │ ❌ │    │ ❌ │ ❌ │ ❌ │ ❌ │  ← Perimeter

Legend:
  ❌ = Perimeter (excluded)
  N  = Negative control (0% inhibition)
  P  = Positive control (100% inhibition)
  ✅ = Test compounds
```

## Calculations

### %Inhibition
```
%Inhibition = 100 × (Neg_Mean - Signal) / (Neg_Mean - Pos_Mean)
```

### Z' Factor (Quality Control)
```
Z' = 1 - 3 × (Neg_SD + Pos_SD) / |Neg_Mean - Pos_Mean|
```
- Z' ≥ 0.5 = Good quality plate

### Efficiency Metrics
```
Per_one = %Inhibition / 100

SPEI = Per_one / (MW × 0.001)    # Size efficiency
PPEI = Per_one / (TPSA × 0.01)   # Polarity efficiency
```

## Scientific Context

This app supports antibiotic discovery research similar to the Stokes 2020 Cell paper that discovered Halicin using deep learning. The screening data quality assessment happens BEFORE AI modeling.

**Best candidates** (northeast corner of Cartesian plot) are compounds that:
- Are small (high SPEI = good activity per unit mass)
- Are non-polar (high PPEI = good membrane penetration)
- Have high %Inhibition

These properties correlate with better drug-likeness and bioavailability.

## License

Research use only.



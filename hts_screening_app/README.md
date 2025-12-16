# 🧬 HTS Screening Data Analysis App

A comprehensive Streamlit application for analyzing High-Throughput Screening (HTS) data from ENAMINE compound library plates. Built for antibiotic discovery research.

## Overview

This app processes 384-well plate luminescence data, calculates %Inhibition, merges with compound library information, and generates analysis plots for drug discovery research.

### Key Features

- **Batch Processing**: Process hundreds of plate files with progress tracking
- **Quality Control**: Z' factor calculation with visual QC summaries
- **Replicate Handling**: Average, single replicate, or overlay both
- **%Inhibition Calculation**: Automatic control well identification and normalization
- **Drug Efficiency Metrics**: SPEI and PPEI calculation
- **Interactive Plots**: Plotly-powered histograms and scatter plots
- **Top Candidate Identification**: Highlight best drug candidates
- **Session Persistence**: Save and load processed data

## Installation

```bash
cd hts_screening_app
pip install -r requirements.txt
```

## Usage

1. **Start the app:**
   ```bash
   streamlit run app.py
   ```

2. **Upload Files:**
   - Upload `Enamine_library.xlsx` (compound library)
   - Upload plate data files (`301-1.xlsx`, `301-2.xlsx`, etc.)

3. **Process and Analyze:**
   - Review QC summary
   - Select replicate handling method
   - View histograms and efficiency metrics
   - Identify top candidates
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


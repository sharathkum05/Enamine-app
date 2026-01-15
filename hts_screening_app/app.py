"""
ENAMINE LIBRARY Screening Data Analysis Streamlit Application

A comprehensive application for analyzing High-Throughput Screening data
from ENAMINE compound library plates.
"""

import io
import os
import pickle
from io import BytesIO
from zipfile import ZipFile

import numpy as np
import pandas as pd
import plotly.io as pio
import streamlit as st

from config import (
    APP_TITLE, APP_ICON, NEG_CONTROL_SPANS, POS_CONTROL_SPANS,
    Z_PRIME_THRESHOLD, COLOR_REP1, COLOR_REP2
)
from utils.plate_parser import process_plate_file, parse_plate_filename
from utils.inhibition import process_plate_inhibition, expand_control_span
from utils.data_merger import merge_with_library, filter_test_compounds, load_library
from utils.metrics import calculate_metrics, get_top_candidates, aggregate_replicates, calculate_summary_stats
from utils.plotting import (
    plot_inhibition_histogram, plot_normalized_histogram, 
    plot_cartesian, plot_qc_heatmap
)
from utils.structure_viewer import smiles_to_image, get_molecule_info


def create_png_download_button(fig, filename, button_key):
    """Create a PNG download button for a Plotly figure."""
    try:
        img_bytes = pio.to_image(fig, format='png', width=1200, height=600, scale=2)
        st.download_button(
            label="📥 Download as PNG",
            data=img_bytes,
            file_name=f"{filename}.png",
            mime="image/png",
            key=button_key
        )
    except Exception as e:
        st.error(f"PNG export failed: {e}")


# Page configuration
st.set_page_config(
    page_title="ENAMINE LIBRARY Analysis",
    page_icon=APP_ICON,
    layout="wide",
    initial_sidebar_state="expanded"
)

def initialize_session_state():
    """Initialize session state variables."""
    if 'library_df' not in st.session_state:
        st.session_state.library_df = None
    if 'main_df' not in st.session_state:
        st.session_state.main_df = None
    if 'qc_df' not in st.session_state:
        st.session_state.qc_df = None
    if 'analysis_df' not in st.session_state:
        st.session_state.analysis_df = None
    if 'replicate_mode' not in st.session_state:
        st.session_state.replicate_mode = 'average'
    if 'processing_done' not in st.session_state:
        st.session_state.processing_done = False
    if 'excluded_plates' not in st.session_state:
        st.session_state.excluded_plates = set()


def process_plate_files(plate_files, library_df, progress_bar, status_text):
    """Process all plate files and return combined data."""
    all_data = []
    qc_data = []
    
    total_files = len(plate_files)
    
    for i, plate_file in enumerate(plate_files):
        try:
            # Update progress
            progress = (i + 1) / total_files
            progress_bar.progress(progress)
            status_text.text(f"Processing {plate_file.name}... ({i+1}/{total_files})")
            
            # Read file content
            file_content = BytesIO(plate_file.read())
            plate_file.seek(0)  # Reset for potential re-read
            
            # Process plate file
            plate_data, plate_num, rep = process_plate_file(file_content, plate_file.name)
            
            # Calculate %Inhibition and QC stats
            plate_data, qc_stats = process_plate_inhibition(
                plate_data, NEG_CONTROL_SPANS, POS_CONTROL_SPANS
            )
            
            all_data.append(plate_data)
            qc_data.append(qc_stats)
            
        except Exception as e:
            st.warning(f"⚠️ Error processing {plate_file.name}: {str(e)}")
            continue
    
    if not all_data:
        return None, None
    
    # Combine all plate data
    combined_df = pd.concat(all_data, ignore_index=True)
    qc_df = pd.DataFrame(qc_data)
    
    # Merge with library
    status_text.text("Merging with ENAMINE library...")
    combined_df = merge_with_library(combined_df, library_df)
    
    # Filter to test compounds only
    combined_df = filter_test_compounds(combined_df)
    
    # Calculate metrics
    status_text.text("Calculating efficiency metrics...")
    combined_df = calculate_metrics(combined_df)
    
    progress_bar.progress(1.0)
    status_text.text("✅ Processing complete!")
    
    return combined_df, qc_df


def render_upload_section():
    """Render the file upload section."""
    st.header("📁 Data Upload")
    
    st.info("""
    📋 **Required Files:**
    1. **ENAMINE Library** (xlsx) - Compound library with SMILES, MW, TPSA data
    2. **Plate Data Files** (xlsx) - Plate reader output files
    
    ⚠️ **Note:** You must upload the ENAMINE library file each session (file not stored on server for security/size reasons)
    """)
    
    st.divider()
    
    # ENAMINE Library uploader - REQUIRED every session
    st.subheader("1️⃣ ENAMINE Library File")
    
    if st.session_state.get('library_loaded', False):
        num_compounds = len(st.session_state.library_df)
        try:
            num_plates = st.session_state.library_df['Plate_ID'].nunique()
            st.success(f"✅ Library loaded: {num_compounds:,} compounds across {num_plates} plates")
        except:
            st.success(f"✅ Library loaded: {num_compounds:,} compounds")
        
        # Option to re-upload if needed
        if st.button("🔄 Upload Different Library File"):
            st.session_state.library_df = None
            st.session_state.library_loaded = False
            st.rerun()
    else:
        st.warning("⚠️ Please upload the ENAMINE library file to begin")
        
        with st.expander("📖 Required Columns in Library File", expanded=False):
            st.markdown("""
            The ENAMINE library Excel file must contain these columns:
            - `Plate_ID` - Plate identifier (e.g., "2096462-Y10-001")
            - `Well` - Well position (e.g., "B03")
            - `Smiles` - SMILES molecular structure
            - `MW` - Molecular weight
            - `TPSA` - Topological polar surface area
            - `catalog number` - Compound identifier
            """)
        
        uploaded_library = st.file_uploader(
            "Upload ENAMINE Library (xlsx)",
            type=['xlsx'],
            key='library_upload',
            help="Upload the Enamine_library.xlsx file containing compound information"
        )
        
        if uploaded_library:
            with st.spinner("Loading library file..."):
                try:
                    library_df = pd.read_excel(uploaded_library)
                    
                    # Validate required columns
                    required_cols = ['Plate_ID', 'Well', 'Smiles', 'MW', 'TPSA', 'catalog number']
                    missing_cols = [col for col in required_cols if col not in library_df.columns]
                    
                    if missing_cols:
                        st.error(f"❌ Missing required columns: {', '.join(missing_cols)}")
                        st.info("Please ensure your file has all required columns listed above.")
                    else:
                        st.session_state.library_df = library_df
                        st.session_state.library_loaded = True
                        st.success("✅ Library loaded successfully!")
                        st.rerun()
                        
                except Exception as e:
                    st.error(f"❌ Error loading file: {str(e)}")
    
    st.divider()
    
    # Plate Data Files uploader
    st.subheader("2️⃣ Plate Data Files")
    st.write("Upload plate reader output files (*.xlsx)")
    
    plate_files = st.file_uploader(
        "Upload plate data files",
        type=["xlsx"],
        accept_multiple_files=True,
        key='plate_uploader',
        label_visibility="collapsed",
        help="Upload one or more plate data files (e.g., 301-1.xlsx, 301-2.xlsx)"
    )
    
    if plate_files:
        st.info(f"📊 {len(plate_files)} plate files selected")
    
    st.divider()
    
    # Alternative: Load saved data
    st.subheader("💾 Or Load Saved Session")
    saved_file = st.file_uploader(
        "Load previously processed data (.pkl)",
        type=['pkl'],
        key='saved_uploader'
    )
    
    if saved_file:
        try:
            saved_data = pickle.load(saved_file)
            st.session_state.main_df = saved_data.get('main_df')
            st.session_state.qc_df = saved_data.get('qc_df')
            st.session_state.library_df = saved_data.get('library_df')
            st.session_state.processing_done = True
            st.success("✅ Session data loaded successfully!")
            st.rerun()
        except Exception as e:
            st.error(f"Error loading saved data: {str(e)}")
    
    # Process button - check if library is loaded
    if not st.session_state.get('library_loaded', False):
        st.warning("⚠️ Cannot process files: ENAMINE Library not found. Please ensure Enamine_library.xlsx is in the app folder.")
    elif plate_files and st.session_state.library_df is not None:
        if st.button("🚀 Process Files", type="primary", use_container_width=True):
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            main_df, qc_df = process_plate_files(
                plate_files, 
                st.session_state.library_df,
                progress_bar,
                status_text
            )
            
            if main_df is not None:
                st.session_state.main_df = main_df
                st.session_state.qc_df = qc_df
                st.session_state.processing_done = True
                st.rerun()
            else:
                st.error("❌ No plates were successfully processed.")


def render_qc_section():
    """Render the QC summary section."""
    st.markdown("## 📋 Quality Control Summary")
    
    qc_df = st.session_state.qc_df.copy()
    
    # Summary metrics
    col1, col2, col3, col4 = st.columns(4)
    
    total_plates = len(qc_df)
    passed_plates = qc_df['QC_Pass'].sum()
    failed_plates = total_plates - passed_plates
    avg_zprime = qc_df['Z_Prime'].mean()
    
    with col1:
        st.metric("Total Plates", total_plates)
    
    with col2:
        st.metric("Passed QC", passed_plates)
    
    with col3:
        st.metric("Failed QC", failed_plates)
    
    with col4:
        st.metric("Avg Z' Factor", f"{avg_zprime:.3f}")
    
    st.markdown("")
    
    # QC Table with exclusion checkboxes
    st.markdown("### Plate Quality Details")
    
    # Add status column
    qc_df['Status'] = qc_df['QC_Pass'].apply(lambda x: '✅ Pass' if x else '⚠️ Fail')
    qc_df['Z_Prime'] = qc_df['Z_Prime'].round(3)
    qc_df['Neg_Mean'] = qc_df['Neg_Mean'].round(1)
    qc_df['Neg_SD'] = qc_df['Neg_SD'].round(1)
    qc_df['Pos_Mean'] = qc_df['Pos_Mean'].round(1)
    qc_df['Pos_SD'] = qc_df['Pos_SD'].round(1)
    
    # Display columns
    display_cols = ['Plate', 'Rep', 'Neg_Mean', 'Neg_SD', 'Pos_Mean', 'Pos_SD', 'Z_Prime', 'Status']
    
    st.dataframe(
        qc_df[display_cols],
        use_container_width=True,
        height=400
    )
    
    # Z' Heatmap
    with st.expander("📊 Z' Factor Heatmap", expanded=False):
        fig = plot_qc_heatmap(st.session_state.qc_df)
        st.plotly_chart(fig, use_container_width=True)
    
    # Save processed data
    st.markdown("")
    col1, col2 = st.columns(2)
    
    with col1:
        if st.button("💾 Save Processed Data", use_container_width=True):
            save_data = {
                'main_df': st.session_state.main_df,
                'qc_df': st.session_state.qc_df,
                'library_df': st.session_state.library_df
            }
            buffer = BytesIO()
            pickle.dump(save_data, buffer)
            buffer.seek(0)
            
            st.download_button(
                label="📥 Download Session File",
                data=buffer,
                file_name="hts_processed_data.pkl",
                mime="application/octet-stream"
            )


def render_replicate_section():
    """Render the replicate selection section."""
    st.markdown("## 🔄 Replicate Handling")
    
    st.info("Choose how to handle replicate measurements for downstream analysis. Different approaches may be appropriate depending on your data quality and analysis goals.")
    
    replicate_mode = st.radio(
        "Select replicate handling method:",
        options=['average', 'rep1', 'rep2', 'both'],
        format_func=lambda x: {
            'average': '📊 Average (mean of Rep1 & Rep2) - Recommended',
            'rep1': '1️⃣ Rep1 only',
            'rep2': '2️⃣ Rep2 only',
            'both': '🔀 Show both (overlay on plots)'
        }[x],
        index=['average', 'rep1', 'rep2', 'both'].index(st.session_state.replicate_mode),
        horizontal=True
    )
    
    if replicate_mode != st.session_state.replicate_mode:
        st.session_state.replicate_mode = replicate_mode
        # Update analysis dataframe
        st.session_state.analysis_df = aggregate_replicates(
            st.session_state.main_df.copy(),
            replicate_mode
        )


def render_plate_selector():
    """Render plate selector for individual or all plate analysis."""
    st.markdown("## 📊 Plate Selection")
    st.write("Analyze individual plates or all plates combined")
    
    # Check if data is available
    if st.session_state.analysis_df is None:
        st.warning("⚠️ Please upload and process data first before selecting plates.")
        return
    
    # Get unique plate numbers from the data
    available_plates = sorted(st.session_state.analysis_df['Plate'].unique())
    
    # Create selection options: "All Plates" + individual plate numbers
    plate_options = ["All Plates"] + [f"Plate {p}" for p in available_plates]
    
    selected_plate = st.selectbox(
        "Select plate to analyze:",
        options=plate_options,
        index=0  # Default to "All Plates"
    )
    
    # Filter data based on selection
    if selected_plate == "All Plates":
        filtered_df = st.session_state.analysis_df.copy()
        st.success(f"📊 Analyzing ALL {len(available_plates)} plates ({len(filtered_df):,} compounds)")
    else:
        plate_num = int(selected_plate.replace("Plate ", ""))
        filtered_df = st.session_state.analysis_df[st.session_state.analysis_df['Plate'] == plate_num].copy()
        st.success(f"📊 Analyzing {selected_plate} only ({len(filtered_df):,} compounds)")
    
    # Store filtered data for use in all plot sections
    st.session_state.current_analysis_df = filtered_df


def render_histogram_section():
    """Render Task 1A: %Inhibition Histogram."""
    st.markdown("## 📈 Task 1A: %Inhibition Distribution")
    
    if st.session_state.analysis_df is None:
        st.session_state.analysis_df = aggregate_replicates(
            st.session_state.main_df.copy(),
            st.session_state.replicate_mode
        )
    
    # Use filtered data if plate selector has been used
    df = st.session_state.get('current_analysis_df', st.session_state.analysis_df)
    
    # Create histogram
    fig = plot_inhibition_histogram(df, st.session_state.replicate_mode)
    st.plotly_chart(fig, use_container_width=True)
    create_png_download_button(fig, "percent_inhibition_histogram", "png_inhibition")
    
    # Statistics
    stats = calculate_summary_stats(df, 'Pct_Inhibition')
    
    col1, col2, col3, col4, col5, col6 = st.columns(6)
    
    with col1:
        st.metric("Count", f"{stats['count']:,}")
    with col2:
        st.metric("Mean", f"{stats['mean']:.1f}%")
    with col3:
        st.metric("Median", f"{stats['median']:.1f}%")
    with col4:
        st.metric("Std Dev", f"{stats['std']:.1f}%")
    with col5:
        st.metric("Min", f"{stats['min']:.1f}%")
    with col6:
        st.metric("Max", f"{stats['max']:.1f}%")
    
    # Download data
    csv = df[['Plate', 'Well', 'Pct_Inhibition']].to_csv(index=False)
    st.download_button(
        "📥 Download Data (CSV)",
        data=csv,
        file_name="inhibition_data.csv",
        mime="text/csv"
    )


def render_normalized_histogram_section():
    """Render Task 1B: Size-Normalized Activity Histograms colored by 10xPSA/MW bins."""
    st.markdown("## 📈 Task 1B: Size-Normalized Activity")
    
    # Use filtered data if plate selector has been used
    df = st.session_state.get('current_analysis_df', st.session_state.analysis_df)
    
    # Ensure metrics are calculated
    if 'SPEI' not in df.columns or 'PSA_MW_ratio' not in df.columns:
        df = calculate_metrics(df)
        if st.session_state.get('current_analysis_df') is not None:
            st.session_state.current_analysis_df = df
        else:
            st.session_state.analysis_df = df
    
    # SPEI Histogram - Stacked by 10xPSA/MW bins
    st.subheader("Size-Normalized Activity (SPEI)")
    st.info("**Formula:** SPEI = (% Inhibition / 100) / (MW × 0.001)\n\n*Histogram colored by 10×PSA/MW ratio bins to show relationship between size efficiency and polarity/size ratio.*")
    
    # Filter to only positive SPEI values
    spei_df = df[df['SPEI'] > 0].copy()
    
    st.info(f"📊 Showing {len(spei_df):,} compounds with SPEI > 0")
    
    # Create bins for 10xPSA/MW
    num_bins = 6
    spei_df['10PSAoMW_bin'] = pd.cut(
        spei_df['PSA_MW_ratio'], 
        bins=num_bins, 
        precision=2
    ).apply(lambda x: f"{x.left:.2f} - {x.right:.2f}" if pd.notna(x) else "N/A")
    
    # Define color palette (matching professor's image)
    color_palette = ['#87CEEB', '#1a237e', '#7CB342', '#2E7D32', '#FF7043', '#FFD54F']
    
    # Get unique bins in order
    bin_order = spei_df.groupby('10PSAoMW_bin')['PSA_MW_ratio'].min().sort_values().index.tolist()
    
    # Create stacked histogram
    import plotly.express as px
    fig_spei = px.histogram(
        spei_df,
        x='SPEI',
        color='10PSAoMW_bin',
        category_orders={'10PSAoMW_bin': bin_order},
        color_discrete_sequence=color_palette,
        nbins=30,
        title='Count vs SPEI',
        labels={'SPEI': 'SPEI', 'count': 'Count', '10PSAoMW_bin': '10PSAoMW'}
    )
    
    fig_spei.update_layout(
        xaxis_title="SPEI",
        yaxis_title="Count",
        legend_title="10PSAoMW",
        barmode='stack',
        template='plotly_white',
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.2,
            xanchor="center",
            x=0.5
        )
    )
    fig_spei.update_xaxes(range=[0, None])
    
    st.plotly_chart(fig_spei, use_container_width=True)
    create_png_download_button(fig_spei, "spei_histogram_colored", "png_spei_colored")
    
    # Statistics using filtered data
    stats = calculate_summary_stats(spei_df, 'SPEI')
    
    col1, col2, col3, col4, col5, col6 = st.columns(6)
    
    with col1:
        st.metric("Count", f"{stats['count']:,}")
    with col2:
        st.metric("Mean", f"{stats['mean']:.2f}")
    with col3:
        st.metric("Median", f"{stats['median']:.2f}")
    with col4:
        st.metric("Std Dev", f"{stats['std']:.2f}")
    with col5:
        st.metric("Min", f"{stats['min']:.2f}")
    with col6:
        st.metric("Max", f"{stats['max']:.2f}")
    
    # Download filtered data
    csv = spei_df[['Plate', 'Well', 'Pct_Inhibition', 'MW', 'SPEI', 'PSA_MW_ratio', '10PSAoMW_bin']].to_csv(index=False)
    st.download_button(
        "📥 Download Data (CSV)",
        data=csv,
        file_name="spei_data_with_bins.csv",
        mime="text/csv",
        key="download_spei_binned"
    )
    
    st.divider()
    
    # PPEI Histogram - Stacked by 10xPSA/MW bins
    st.subheader("Polarity-Normalized Activity (PPEI)")
    st.info("**Formula:** PPEI = (% Inhibition / 100) / (TPSA × 0.01)\n\n*Histogram colored by 10×PSA/MW ratio bins to show relationship between polarity efficiency and polarity/size ratio.*")
    
    # Filter to only positive PPEI values
    ppei_df = df[df['PPEI'] > 0].copy()
    
    st.info(f"📊 Showing {len(ppei_df):,} compounds with PPEI > 0")
    
    # Create bins for 10xPSA/MW (same bins as SPEI)
    ppei_df['10PSAoMW_bin'] = pd.cut(
        ppei_df['PSA_MW_ratio'], 
        bins=num_bins, 
        precision=2
    ).apply(lambda x: f"{x.left:.2f} - {x.right:.2f}" if pd.notna(x) else "N/A")
    
    # Get unique bins in order
    bin_order_ppei = ppei_df.groupby('10PSAoMW_bin')['PSA_MW_ratio'].min().sort_values().index.tolist()
    
    # Create stacked histogram
    fig_ppei = px.histogram(
        ppei_df,
        x='PPEI',
        color='10PSAoMW_bin',
        category_orders={'10PSAoMW_bin': bin_order_ppei},
        color_discrete_sequence=color_palette,
        nbins=30,
        title='Count vs PPEI',
        labels={'PPEI': 'PPEI', 'count': 'Count', '10PSAoMW_bin': '10PSAoMW'}
    )
    
    fig_ppei.update_layout(
        xaxis_title="PPEI",
        yaxis_title="Count",
        legend_title="10PSAoMW",
        barmode='stack',
        template='plotly_white',
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.2,
            xanchor="center",
            x=0.5
        )
    )
    fig_ppei.update_xaxes(range=[0, None])
    
    st.plotly_chart(fig_ppei, use_container_width=True)
    create_png_download_button(fig_ppei, "ppei_histogram_colored", "png_ppei_colored")
    
    # PPEI Statistics using filtered data
    ppei_stats = calculate_summary_stats(ppei_df, 'PPEI')
    
    col1, col2, col3, col4, col5, col6 = st.columns(6)
    
    with col1:
        st.metric("Count", f"{ppei_stats['count']:,}")
    with col2:
        st.metric("Mean", f"{ppei_stats['mean']:.2f}")
    with col3:
        st.metric("Median", f"{ppei_stats['median']:.2f}")
    with col4:
        st.metric("Std Dev", f"{ppei_stats['std']:.2f}")
    with col5:
        st.metric("Min", f"{ppei_stats['min']:.2f}")
    with col6:
        st.metric("Max", f"{ppei_stats['max']:.2f}")
    
    # Download filtered PPEI data
    csv = ppei_df[['Plate', 'Well', 'Pct_Inhibition', 'TPSA', 'PPEI', 'PSA_MW_ratio', '10PSAoMW_bin']].to_csv(index=False)
    st.download_button(
        "📥 Download Data (CSV)",
        data=csv,
        file_name="ppei_data_with_bins.csv",
        mime="text/csv",
        key="download_ppei_binned"
    )


def render_psa_mw_histogram_section():
    """Render Task 1C: 10xPSA/MW Distribution Histogram."""
    st.markdown("## 📊 Task 1C: 10xPSA/MW Distribution")
    
    st.info("""
**Formula:** 10 × TPSA / MW (equivalent to SPEI / PPEI)

*This ratio indicates the balance between polarity and size. Lower values suggest compounds that are relatively small for their polarity - potentially better drug-like properties.*
""")
    
    # Use filtered data if plate selector has been used
    df = st.session_state.get('current_analysis_df', st.session_state.analysis_df)
    
    # Ensure metrics are calculated
    if 'PSA_MW_ratio' not in df.columns:
        df = calculate_metrics(df)
        if 'current_analysis_df' in st.session_state:
            st.session_state.current_analysis_df = df
        else:
            st.session_state.analysis_df = df
    
    # Filter to only positive values
    total_psa_mw = len(df)
    psa_mw_df = df[df['PSA_MW_ratio'] > 0].copy()
    filtered_psa_mw = len(psa_mw_df)
    excluded_psa_mw = total_psa_mw - filtered_psa_mw
    
    st.info(f"📊 Showing {filtered_psa_mw:,} compounds with 10xPSA/MW > 0 ({excluded_psa_mw:,} negative values excluded)")
    
    # Create histogram with filtered data
    import plotly.express as px
    fig_psa_mw = px.histogram(
        psa_mw_df,
        x='PSA_MW_ratio',
        nbins=50,
        title='Distribution of 10×PSA/MW',
        labels={'PSA_MW_ratio': '10 × TPSA / MW'},
        color_discrete_sequence=['#9B59B6']  # Purple color
    )
    
    fig_psa_mw.update_layout(
        xaxis_title="10 × TPSA / MW",
        yaxis_title="Count",
        showlegend=False,
        template='plotly_white'
    )
    
    fig_psa_mw.update_xaxes(range=[0, None])  # Start x-axis from 0
    st.plotly_chart(fig_psa_mw, use_container_width=True)
    create_png_download_button(fig_psa_mw, "10xPSA_MW_histogram", "png_psa_mw")
    
    # Statistics using filtered data
    psa_mw_stats = calculate_summary_stats(psa_mw_df, 'PSA_MW_ratio')
    
    col1, col2, col3, col4, col5, col6 = st.columns(6)
    
    with col1:
        st.metric("Count", f"{psa_mw_stats['count']:,}")
    with col2:
        st.metric("Mean", f"{psa_mw_stats['mean']:.2f}")
    with col3:
        st.metric("Median", f"{psa_mw_stats['median']:.2f}")
    with col4:
        st.metric("Std Dev", f"{psa_mw_stats['std']:.2f}")
    with col5:
        st.metric("Min", f"{psa_mw_stats['min']:.2f}")
    with col6:
        st.metric("Max", f"{psa_mw_stats['max']:.2f}")
    
    # Download filtered data
    available_cols = ['Plate', 'Well', 'catalog_number', 'PSA_MW_ratio', 'SPEI', 'PPEI', 'MW', 'TPSA']
    export_cols = [c for c in available_cols if c in psa_mw_df.columns]
    csv = psa_mw_df[export_cols].to_csv(index=False)
    st.download_button(
        "📥 Download Data (CSV)",
        data=csv,
        file_name="10xPSA_MW_data.csv",
        mime="text/csv",
        key="download_psa_mw_csv"
    )


def render_metrics_section():
    """Render SPEI & PPEI metrics overview."""
    st.markdown("## 🧮 Efficiency Metrics")
    
    st.info("**SPEI** = per_one / (MW × 0.001) — Size efficiency: Higher = more potent per unit mass\n\n**PPEI** = per_one / (TPSA × 0.01) — Polarity efficiency: Higher = better membrane penetration")
    
    # Use filtered data if plate selector has been used
    df = st.session_state.get('current_analysis_df', st.session_state.analysis_df)
    
    # Data preview
    preview_cols = ['Plate', 'Well', 'Pct_Inhibition', 'MW', 'TPSA', 'SPEI', 'PPEI', 'catalog_number']
    available_cols = [c for c in preview_cols if c in df.columns]
    
    st.dataframe(
        df[available_cols].head(20).round(2),
        use_container_width=True
    )
    
    # Summary statistics for both metrics
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### SPEI Statistics")
        spei_stats = calculate_summary_stats(df, 'SPEI')
        for key, value in spei_stats.items():
            if key == 'count':
                st.write(f"**{key.title()}:** {value:,}")
            else:
                st.write(f"**{key.title()}:** {value:.3f}")
    
    with col2:
        st.markdown("### PPEI Statistics")
        ppei_stats = calculate_summary_stats(df, 'PPEI')
        for key, value in ppei_stats.items():
            if key == 'count':
                st.write(f"**{key.title()}:** {value:,}")
            else:
                st.write(f"**{key.title()}:** {value:.3f}")


def render_cartesian_section():
    """Render Cartesian Plot (PPEI vs SPEI)."""
    st.markdown("## 🎯 Cartesian Plot: Best Candidates")
    
    st.info("Compounds in the **northeast corner** (high SPEI + high PPEI) are the best drug candidates: they're small, non-polar, and highly active.")
    
    # Use filtered data if plate selector has been used
    df = st.session_state.get('current_analysis_df', st.session_state.analysis_df)
    
    # Check if data is available
    if df is None or df.empty:
        st.warning("⚠️ No analysis data available. Please upload and process data first.")
        return
    
    # Ensure metrics are calculated
    if 'SPEI' not in df.columns or 'PPEI' not in df.columns or 'PSA_MW_ratio' not in df.columns:
        df = calculate_metrics(df)
        if st.session_state.get('current_analysis_df') is not None:
            st.session_state.current_analysis_df = df
        else:
            st.session_state.analysis_df = df
    
    # Filter data to only show positive values (cut beyond origin 0,0)
    total_points = len(df)
    cartesian_df = df[(df['SPEI'] > 0) & (df['PPEI'] > 0)].copy()
    filtered_points = len(cartesian_df)
    excluded_points = total_points - filtered_points
    
    st.info(f"📊 Showing {filtered_points:,} compounds with SPEI > 0 and PPEI > 0. ({excluded_points:,} compounds with negative values excluded)")
    
    # CONTROLS ROW 1: Top candidates threshold and ranking metric
    col1, col2 = st.columns(2)
    
    with col1:
        top_percentile = st.slider(
            "Highlight Top Candidates (%)",
            min_value=1,
            max_value=20,
            value=5,
            step=1,
            help="Select the top X% of compounds to highlight",
            key="top_percentile_slider"
        )
    
    with col2:
        metric = st.radio(
            "Ranking Metric",
            options=['combined', 'spei', 'ppei'],
            format_func=lambda x: {
                'combined': 'SPEI + PPEI (Combined)',
                'spei': 'SPEI only',
                'ppei': 'PPEI only'
            }[x],
            horizontal=True,
            key="ranking_metric_radio"
        )
    
    # ============================================
    # PLOT 1: Original - Top Candidates with Red Stars
    # ============================================
    st.markdown("### 📍 View 1: Top Candidates Highlighted")
    st.caption("Gray dots = all compounds | Red stars = top candidates")
    
    # Calculate top candidates
    from utils.metrics import get_top_candidates
    df_ranked = get_top_candidates(cartesian_df.copy(), top_percentile, metric)
    top_df_viz = df_ranked[df_ranked['Is_Top']].copy()
    regular_df = df_ranked[~df_ranked['Is_Top']].copy()
    
    # CREATE ORIGINAL PLOT with red stars
    import plotly.graph_objects as go
    
    fig1 = go.Figure()
    
    # Plot regular compounds (gray)
    fig1.add_trace(go.Scatter(
        x=regular_df['PPEI'],
        y=regular_df['SPEI'],
        mode='markers',
        marker=dict(
            size=6,
            color='rgba(150, 150, 150, 0.5)',
            opacity=0.6
        ),
        text=regular_df.get('catalog_number', ''),
        hovertemplate=(
            "<b>%{text}</b><br>" +
            "SPEI: %{y:.3f}<br>" +
            "PPEI: %{x:.3f}<br>" +
            "<extra></extra>"
        ),
        name='Compounds',
        showlegend=True
    ))
    
    # Plot top candidates as red stars
    if not top_df_viz.empty:
        fig1.add_trace(go.Scatter(
            x=top_df_viz['PPEI'],
            y=top_df_viz['SPEI'],
            mode='markers',
            marker=dict(
                size=12,
                symbol='star',
                color='#FF0000',
                line=dict(color='#8B0000', width=1.5)
            ),
            text=top_df_viz.get('catalog_number', ''),
            hovertemplate=(
                "<b>⭐ %{text}</b><br>" +
                "SPEI: %{y:.3f}<br>" +
                "PPEI: %{x:.3f}<br>" +
                "<extra></extra>"
            ),
            name=f'Top {top_percentile}% (Best Candidates)',
            showlegend=True
        ))
    
    # Update layout for original plot
    fig1.update_layout(
        title='SPEI vs PPEI: Drug Efficiency Plot',
        xaxis_title='PPEI (Polarity Efficiency) → Higher = Better Membrane Penetration',
        yaxis_title='SPEI (Size Efficiency) → Higher = More Potent per Unit Mass',
        xaxis=dict(range=[0, None]),
        yaxis=dict(range=[0, None]),
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
        ),
        template='plotly_white',
        hovermode='closest'
    )
    
    st.plotly_chart(fig1, use_container_width=True)
    create_png_download_button(fig1, "spei_vs_ppei_top_candidates", "png_cartesian_stars")
    
    # ============================================
    # PLOT 2: New - Colored by 10xPSA/MW Bins with Wedge Lines
    # ============================================
    st.markdown("---")
    st.markdown("### 🎨 View 2: Colored by 10×PSA/MW Ratio")
    st.caption("Categorical coloring showing polarity/size distribution with wedge boundaries")
    
    # Calculate PSA_MW_ratio and create bins
    if 'PSA_MW_ratio' not in cartesian_df.columns:
        cartesian_df['PSA_MW_ratio'] = (10 * cartesian_df['TPSA']) / cartesian_df['MW']
    
    # Create bins with FIXED 0.45 Å²/Dalton intervals (not equal-width)
    min_ratio = 0  # Start from 0
    max_ratio = cartesian_df['PSA_MW_ratio'].max()
    wedge_spacing = 0.45  # Fixed spacing as per professor's specification
    
    # Create bin edges with 0.45 intervals
    bin_edges = np.arange(min_ratio, max_ratio + wedge_spacing, wedge_spacing)
    
    cartesian_df['10PSAoMW_bin'] = pd.cut(
        cartesian_df['PSA_MW_ratio'], 
        bins=bin_edges, 
        precision=2
    ).apply(lambda x: f"{x.left:.2f} - {x.right:.2f}" if pd.notna(x) else "N/A")
    
    # Define color palette - extend to match number of bins
    base_colors = ['#87CEEB', '#1a237e', '#7CB342', '#2E7D32', '#FF7043', '#FFD54F', '#9C27B0', '#FF4081']
    num_actual_bins = len(cartesian_df['10PSAoMW_bin'].unique())
    color_palette = (base_colors * ((num_actual_bins // len(base_colors)) + 1))[:num_actual_bins]
    
    # Get unique bins in order
    bin_order = cartesian_df.groupby('10PSAoMW_bin')['PSA_MW_ratio'].min().sort_values().index.tolist()
    
    # Add bins to df_ranked for the second plot
    df_ranked_with_bins = df_ranked.copy()
    if '10PSAoMW_bin' not in df_ranked_with_bins.columns:
        # Merge bins from cartesian_df
        df_ranked_with_bins = df_ranked_with_bins.merge(
            cartesian_df[['Plate', 'Well', '10PSAoMW_bin']], 
            on=['Plate', 'Well'], 
            how='left'
        )
    
    # ============================================
    # PER-WEDGE TOP CANDIDATE SELECTION
    # ============================================
    # Calculate distance from origin for each compound
    df_ranked_with_bins['Distance_from_Origin'] = np.sqrt(
        df_ranked_with_bins['SPEI']**2 + df_ranked_with_bins['PPEI']**2
    )
    
    # Select top N compounds per wedge (furthest from origin in each wedge)
    compounds_per_wedge = 10  # Number of top compounds to select from each wedge
    
    # Initialize column
    df_ranked_with_bins['Is_Top_PerWedge'] = False
    
    for bin_name in bin_order:
        wedge_compounds_mask = df_ranked_with_bins['10PSAoMW_bin'] == bin_name
        wedge_compounds = df_ranked_with_bins[wedge_compounds_mask]
        
        if not wedge_compounds.empty:
            # Get indices of top compounds in this wedge
            top_indices = wedge_compounds.nlargest(compounds_per_wedge, 'Distance_from_Origin').index
            # Mark them as top per wedge
            df_ranked_with_bins.loc[top_indices, 'Is_Top_PerWedge'] = True
    
    # Options for display
    col_a, col_b = st.columns(2)
    with col_a:
        show_top_only = st.checkbox("Show top candidates only", value=False, key="show_top_only_binned")
    with col_b:
        show_wedge_lines = st.checkbox("Show wedge boundary lines", value=True, key="show_wedge_lines")
    
    # Filter if checkbox selected
    plot_df = df_ranked_with_bins[df_ranked_with_bins['Is_Top']] if show_top_only else df_ranked_with_bins
    
    # CREATE SCATTER PLOT colored by 10PSAoMW bins
    import plotly.express as px
    
    fig2 = px.scatter(
        plot_df,
        x='PPEI',
        y='SPEI',
        color='10PSAoMW_bin',
        category_orders={'10PSAoMW_bin': bin_order},
        color_discrete_sequence=color_palette,
        title='SPEI vs PPEI (Colored by 10×PSA/MW) - Wedges at 0.45 Å²/Dalton intervals',
        labels={'PPEI': 'PPEI', 'SPEI': 'SPEI', '10PSAoMW_bin': '10PSAoMW'},
        hover_data=['catalog_number', 'Pct_Inhibition', 'MW', 'TPSA', 'PSA_MW_ratio']
    )
    
    fig2.update_traces(marker=dict(size=8, opacity=0.7))
    
    # ============================================
    # ADD RADIATING WEDGE LINES FROM ORIGIN
    # ============================================
    if show_wedge_lines:
        # Get max x and y values for extending lines
        max_ppei = plot_df['PPEI'].max()
        max_spei = plot_df['SPEI'].max()
        max_extent = max(max_ppei, max_spei) * 1.1  # Extend slightly beyond data
        
        # Draw wedge boundary lines (diagonal lines from origin)
        for i, edge in enumerate(bin_edges[1:]):  # Skip first edge (0)
            # Line equation: SPEI = (10PSA/MW) × PPEI
            # or y = slope × x, where slope = edge value
            slope = edge
            x_end = max_extent
            y_end = slope * x_end
            
            # Limit to plot bounds
            if y_end > max_extent:
                y_end = max_extent
                x_end = y_end / slope if slope > 0 else max_extent
            
            # Add line from origin (0,0) to calculated endpoint
            fig2.add_shape(
                type="line",
                x0=0, y0=0,
                x1=x_end, y1=y_end,
                line=dict(
                    color="rgba(100, 100, 100, 0.4)",
                    width=1.5,
                    dash="dash"
                ),
                layer="below"
            )
        
        # ============================================
        # ADD MIDPOINT ANNOTATIONS WITH ARROWS
        # ============================================
        # Add annotations at the midpoint of each wedge
        for i in range(len(bin_edges) - 1):
            lower_edge = bin_edges[i]
            upper_edge = bin_edges[i + 1]
            midpoint_slope = (lower_edge + upper_edge) / 2
            
            # Position annotation at a reasonable distance from origin
            annotation_distance = max_extent * 0.6
            x_pos = annotation_distance / np.sqrt(1 + midpoint_slope**2)
            y_pos = midpoint_slope * x_pos
            
            # Get bin label
            bin_label = f"{lower_edge:.2f} - {upper_edge:.2f}"
            
            fig2.add_annotation(
                x=x_pos,
                y=y_pos,
                text=f"<b>{bin_label}</b>",
                showarrow=True,
                arrowhead=2,
                arrowsize=1,
                arrowwidth=1.5,
                arrowcolor="rgba(50, 50, 50, 0.6)",
                ax=30,
                ay=-30,
                font=dict(size=10, color="rgba(50, 50, 50, 0.8)"),
                bgcolor="rgba(255, 255, 255, 0.8)",
                bordercolor="rgba(150, 150, 150, 0.5)",
                borderwidth=1
            )
    
    # ============================================
    # ADD PER-WEDGE TOP CANDIDATES AS STARS
    # ============================================
    # Overlay red stars for per-wedge top candidates
    top_per_wedge_viz = plot_df[plot_df['Is_Top_PerWedge']] if 'Is_Top_PerWedge' in plot_df.columns else pd.DataFrame()
    
    if not top_per_wedge_viz.empty:
        fig2.add_trace(go.Scatter(
            x=top_per_wedge_viz['PPEI'],
            y=top_per_wedge_viz['SPEI'],
            mode='markers',
            marker=dict(
                size=14,
                symbol='star',
                color='red',
                line=dict(color='darkred', width=1.5)
            ),
            text=top_per_wedge_viz['catalog_number'],
            customdata=top_per_wedge_viz[['PSA_MW_ratio']],
            hovertemplate=(
                "<b>⭐ Per-Wedge Top: %{text}</b><br>" +
                "SPEI: %{y:.3f}<br>" +
                "PPEI: %{x:.3f}<br>" +
                "10PSA/MW: %{customdata[0]:.3f}<br>" +
                "<extra></extra>"
            ),
            name='Top Per Wedge',
            showlegend=True
        ))
    
    # Update layout
    fig2.update_layout(
        xaxis_title='PPEI (Polarity Efficiency) → Higher = Better Membrane Penetration',
        yaxis_title='SPEI (Size Efficiency) → Higher = More Potent per Unit Mass',
        xaxis=dict(range=[0, None]),
        yaxis=dict(range=[0, None]),
        legend_title="10PSAoMW",
        template='plotly_white',
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.15,
            xanchor="center",
            x=0.5
        ),
        hovermode='closest'
    )
    
    st.plotly_chart(fig2, use_container_width=True)
    create_png_download_button(fig2, "spei_vs_ppei_colored_by_bins", "png_cartesian_colored")
    
    # Top candidates table
    st.markdown("---")
    st.markdown("### 🏆 Top Candidates")
    
    # Use df_ranked_with_bins to include the bin column
    top_df = df_ranked_with_bins[df_ranked_with_bins['Is_Top']].sort_values('Rank_Score', ascending=False)
    
    display_cols = ['Plate', 'Well', 'catalog_number', 'Smiles', 'Chemical_name', 
                    'Pct_Inhibition', 'MW', 'TPSA', 'SPEI', 'PPEI', 'PSA_MW_ratio', '10PSAoMW_bin']
    available_cols = [c for c in display_cols if c in top_df.columns]
    
    st.dataframe(
        top_df[available_cols].head(50).round(2),
        use_container_width=True,
        height=400
    )
    
    # Download top candidates data
    csv = top_df[available_cols].to_csv(index=False)
    st.download_button(
        "📥 Download Top Candidates (CSV with SMILES)",
        data=csv,
        file_name="top_candidates_with_smiles.csv",
        mime="text/csv",
        key="download_top_candidates"
    )
    
    # ============================================
    # CHEMICAL STRUCTURE VIEWER (Optional - requires RDKit)
    # ============================================
    st.markdown("---")
    st.subheader("🔬 View Compound Structure")
    
    # Check if RDKit is available
    from utils.structure_viewer import is_rdkit_available
    
    if not is_rdkit_available():
        st.info("""
        **Structure viewer is not available on this deployment.**
        
        The chemical structure viewer requires RDKit, which is not installed on Streamlit Cloud.
        
        **To view molecular structures:**
        - Run the app locally with RDKit installed: `pip install rdkit`
        - Or use the SMILES strings from the exported CSV in other tools (ChemDraw, PubChem, etc.)
        
        **Note:** All compound data and SMILES strings are available in the table above and CSV exports.
        """)
    else:
        st.write("Select a top candidate to view its chemical structure:")
    
    if not top_df.empty:
        # Create display labels for the dropdown
        top_df_display = top_df.copy()
        top_df_display['display_label'] = (
            top_df_display['catalog_number'].astype(str) + " | " +
            "SPEI: " + top_df_display['SPEI'].round(2).astype(str) + " | " +
            "PPEI: " + top_df_display['PPEI'].round(2).astype(str) + " | " +
            "%Inh: " + top_df_display['Pct_Inhibition'].round(1).astype(str) + "%"
        )
        
        # Dropdown to select compound
        selected_label = st.selectbox(
            "Select a compound:",
            options=["-- Select a compound --"] + top_df_display['display_label'].tolist(),
            key="compound_selector"
        )
        
        if selected_label != "-- Select a compound --":
            # Get the selected compound data
            selected_row = top_df_display[
                top_df_display['display_label'] == selected_label
            ].iloc[0]
            
            # Get SMILES - check for correct column name
            smiles = selected_row.get('Smiles') or selected_row.get('SMILES') or selected_row.get('smiles')
            
            if smiles and not pd.isna(smiles):
                col1, col2 = st.columns([1, 1])
                
                with col1:
                    st.markdown("**Molecular Structure:**")
                    # Generate and display structure image
                    img_buffer = smiles_to_image(smiles, size=(350, 350))
                    if img_buffer:
                        st.image(img_buffer, caption=f"{selected_row['catalog_number']}", width=350)
                    else:
                        st.warning("Could not generate structure image")
                
                with col2:
                    st.markdown("**Compound Details:**")
                    st.markdown(f"**Catalog Number:** `{selected_row['catalog_number']}`")
                    st.markdown(f"**SMILES:** `{str(smiles)[:50]}{'...' if len(str(smiles)) > 50 else ''}`")
                    st.code(smiles, language=None)
                    
                    st.markdown("---")
                    st.markdown("**Screening Metrics:**")
                    
                    metrics_col1, metrics_col2 = st.columns(2)
                    with metrics_col1:
                        st.metric("% Inhibition", f"{selected_row['Pct_Inhibition']:.2f}%")
                        st.metric("SPEI", f"{selected_row['SPEI']:.3f}")
                    with metrics_col2:
                        st.metric("PPEI", f"{selected_row['PPEI']:.3f}")
                        st.metric("MW", f"{selected_row['MW']:.2f}")
                    
                    st.markdown("---")
                    st.markdown("**Location:**")
                    st.write(f"Plate: {selected_row['Plate']} | Well: {selected_row['Well']}")
                    
                    # Get additional molecule info
                    mol_info = get_molecule_info(smiles)
                    if mol_info:
                        st.markdown("---")
                        st.markdown("**Molecular Info:**")
                        st.write(f"Formula: {mol_info['formula']}")
                        st.write(f"Atoms: {mol_info['num_atoms']} | Bonds: {mol_info['num_bonds']} | Rings: {mol_info['num_rings']}")
            else:
                st.warning("No SMILES data available for this compound")
    else:
        st.info("No top candidates to display. Adjust the threshold slider above.")


def render_export_section():
    """Render export options."""
    st.markdown("## 📤 Export Results")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### Complete Dataset")
        
        # Excel export
        buffer = BytesIO()
        with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
            st.session_state.analysis_df.to_excel(writer, sheet_name='Analysis', index=False)
            st.session_state.qc_df.to_excel(writer, sheet_name='QC Summary', index=False)
        buffer.seek(0)
        
        st.download_button(
            "📥 Download Complete Dataset (Excel)",
            data=buffer,
            file_name="hts_analysis_complete.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
    
    with col2:
        st.markdown("### Session Data")
        
        # Pickle export
        save_data = {
            'main_df': st.session_state.main_df,
            'qc_df': st.session_state.qc_df,
            'library_df': st.session_state.library_df,
            'analysis_df': st.session_state.analysis_df,
            'replicate_mode': st.session_state.replicate_mode
        }
        buffer = BytesIO()
        pickle.dump(save_data, buffer)
        buffer.seek(0)
        
        st.download_button(
            "💾 Save Session (PKL)",
            data=buffer,
            file_name="hts_session.pkl",
            mime="application/octet-stream"
        )


def render_sidebar():
    """Render sidebar navigation."""
    with st.sidebar:
        st.markdown(f"# {APP_TITLE}")
        st.markdown("---")
        
        st.markdown("### 📍 Navigation")
        
        sections = []
        
        if not st.session_state.processing_done:
            sections = ['📁 Upload Data']
        else:
            sections = [
                '📋 QC Summary',
                '🔄 Replicate Selection',
                '📊 Plate Selection',
                '📈 Task 1A: %Inhibition',
                '📈 Task 1B: Normalized',
                '📊 Task 1C: 10xPSA/MW',
                '🧮 Metrics',
                '🎯 Cartesian Plot',
                '📤 Export'
            ]
        
        selected = st.radio(
            "Go to section:",
            sections,
            label_visibility='collapsed'
        )
        
        st.markdown("---")
        
        if st.session_state.processing_done:
            st.markdown("### 📊 Data Summary")
            st.write(f"**Compounds:** {len(st.session_state.main_df):,}")
            st.write(f"**Plates:** {st.session_state.qc_df['Plate'].nunique()}")
            st.write(f"**Mode:** {st.session_state.replicate_mode}")
        
        st.markdown("---")
        st.caption("ENAMINE LIBRARY v1.0 | Drug Discovery Pipeline")
        
        return selected


def main():
    """Main application entry point."""
    initialize_session_state()
    
    # Render sidebar and get selection
    selected = render_sidebar()
    
    # Title
    st.markdown(f"# {APP_TITLE}")
    st.markdown("*Antibiotic Discovery • ENAMINE Library • 384-well Plate Analysis*")
    st.divider()
    
    # Route based on selection
    if not st.session_state.processing_done:
        render_upload_section()
    else:
        if selected == '📋 QC Summary':
            render_qc_section()
        elif selected == '🔄 Replicate Selection':
            render_replicate_section()
        elif selected == '📊 Plate Selection':
            render_plate_selector()
        elif selected == '📈 Task 1A: %Inhibition':
            render_histogram_section()
        elif selected == '📈 Task 1B: Normalized':
            render_normalized_histogram_section()
        elif selected == '📊 Task 1C: 10xPSA/MW':
            render_psa_mw_histogram_section()
        elif selected == '🧮 Metrics':
            render_metrics_section()
        elif selected == '🎯 Cartesian Plot':
            render_cartesian_section()
        elif selected == '📤 Export':
            render_export_section()
        else:
            render_qc_section()
    
    # Footer
    st.divider()
    st.caption("Built for Antibiotic Discovery Research • Inspired by Stokes 2020 Halicin Paper")


if __name__ == "__main__":
    main()



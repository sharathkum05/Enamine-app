"""
ITR ENAMINE LIBRARY Screening Data Analysis Streamlit App

A comprehensive application for analyzing High-Throughput Screening data
from ITR ENAMINE compound library plates.
"""

import io
import os
import pickle
from io import BytesIO
from zipfile import ZipFile

import numpy as np
import pandas as pd
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


# Page configuration
st.set_page_config(
    page_title="ITR ENAMINE LIBRARY Analysis",
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
    """Render Task 1B: Size-Normalized Activity Histogram."""
    st.markdown("## 📈 Task 1B: Size-Normalized Activity")
    
    # Use filtered data if plate selector has been used
    df = st.session_state.get('current_analysis_df', st.session_state.analysis_df)
    
    # Ensure metrics are calculated
    if 'SPEI' not in df.columns:
        df = calculate_metrics(df)
        if st.session_state.get('current_analysis_df') is not None:
            st.session_state.current_analysis_df = df
        else:
            st.session_state.analysis_df = df
    
    # SPEI Histogram
    st.subheader("Size-Normalized Activity (SPEI)")
    st.info("**Formula:** SPEI = (% Inhibition / 100) / (MW × 0.001)\n\n*Higher values indicate more potent compounds per unit mass - smaller molecules with same activity are preferred.*")
    
    # Filter to only positive SPEI values
    total_spei = len(df)
    spei_df = df[df['SPEI'] > 0].copy()
    filtered_spei = len(spei_df)
    excluded_spei = total_spei - filtered_spei
    
    st.info(f"📊 Showing {filtered_spei:,} compounds with SPEI > 0 ({excluded_spei:,} negative values excluded)")
    
    # Create histogram with filtered data
    fig = plot_normalized_histogram(spei_df, st.session_state.replicate_mode)
    fig.update_xaxes(range=[0, None])  # Start x-axis from 0
    st.plotly_chart(fig, use_container_width=True)
    
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
    csv = spei_df[['Plate', 'Well', 'Pct_Inhibition', 'MW', 'SPEI']].to_csv(index=False)
    st.download_button(
        "📥 Download Data (CSV)",
        data=csv,
        file_name="spei_data.csv",
        mime="text/csv"
    )
    
    st.divider()
    
    # PPEI Histogram
    st.subheader("Polarity-Normalized Activity (PPEI)")
    st.info("**Formula:** PPEI = (% Inhibition / 100) / (TPSA × 0.01)\n\n*Higher values indicate more potent compounds per unit polarity - less polar molecules with same activity are preferred for membrane penetration.*")
    
    # Filter to only positive PPEI values
    total_ppei = len(df)
    ppei_df = df[df['PPEI'] > 0].copy()
    filtered_ppei = len(ppei_df)
    excluded_ppei = total_ppei - filtered_ppei
    
    st.info(f"📊 Showing {filtered_ppei:,} compounds with PPEI > 0 ({excluded_ppei:,} negative values excluded)")
    
    # Create PPEI histogram using plotly express with filtered data
    import plotly.express as px
    fig_ppei = px.histogram(
        ppei_df,  # Use filtered data
        x='PPEI',
        nbins=50,
        title='Polarity-Normalized Activity (PPEI)',
        labels={'PPEI': 'PPEI = (% Inhibition / 100) / (TPSA × 0.01)'},
        color_discrete_sequence=['#FF6B6B']  # Red color
    )
    fig_ppei.update_layout(
        xaxis_title="PPEI = (% Inhibition / 100) / (TPSA × 0.01)",
        yaxis_title="Count",
        showlegend=False,
        template='plotly_white'
    )
    fig_ppei.update_xaxes(range=[0, None])  # Start x-axis from 0
    st.plotly_chart(fig_ppei, use_container_width=True)
    
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
    csv = ppei_df[['Plate', 'Well', 'Pct_Inhibition', 'TPSA', 'PPEI']].to_csv(index=False)
    st.download_button(
        "📥 Download Data (CSV)",
        data=csv,
        file_name="ppei_data.csv",
        mime="text/csv"
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
    
    # Filter data to only show positive values (cut beyond origin 0,0)
    total_points = len(df)
    cartesian_df = df[(df['SPEI'] > 0) & (df['PPEI'] > 0)].copy()
    filtered_points = len(cartesian_df)
    excluded_points = total_points - filtered_points
    
    st.info(f"📊 Showing {filtered_points:,} compounds with SPEI > 0 and PPEI > 0. ({excluded_points:,} compounds with negative values excluded)")
    
    # Controls
    col1, col2 = st.columns(2)
    
    with col1:
        top_percentile = st.slider(
            "Highlight Top Candidates (%)",
            min_value=1,
            max_value=20,
            value=5,
            step=1,
            help="Select the top X% of compounds to highlight"
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
            horizontal=True
        )
    
    # Create plot using filtered cartesian_df
    fig = plot_cartesian(cartesian_df, top_percentile, st.session_state.replicate_mode, metric)
    
    # Set axes to start from 0
    fig.update_xaxes(range=[0, None])
    fig.update_yaxes(range=[0, None])
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Top candidates table
    st.markdown("### 🏆 Top Candidates")
    
    df_ranked = get_top_candidates(cartesian_df.copy(), top_percentile, metric)
    top_df = df_ranked[df_ranked['Is_Top']].sort_values('Rank_Score', ascending=False)
    
    display_cols = ['Plate', 'Well', 'catalog_number', 'Chemical_name', 
                    'Pct_Inhibition', 'MW', 'TPSA', 'SPEI', 'PPEI']
    available_cols = [c for c in display_cols if c in top_df.columns]
    
    st.dataframe(
        top_df[available_cols].head(50).round(2),
        use_container_width=True,
        height=400
    )
    
    # Download top candidates data
    csv = top_df[available_cols].to_csv(index=False)
    st.download_button(
        "📥 Download Top Candidates (CSV)",
        data=csv,
        file_name="top_candidates.csv",
        mime="text/csv"
    )


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
        st.caption("ITR ENAMINE LIBRARY v1.0 | Drug Discovery Pipeline")
        
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



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
import plotly.io as pio

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


@st.cache_data
def load_enamine_library():
    """Auto-load ENAMINE library from backend."""
    library_path = os.path.join(os.path.dirname(__file__), "Enamine_library.xlsx")
    if os.path.exists(library_path):
        try:
            df = pd.read_excel(library_path)
            return df
        except Exception as e:
            st.error(f"Error loading ENAMINE library: {str(e)}")
            return None
    else:
        return None


@st.cache_data
def load_plate_layout():
    """Auto-load plate layout from backend."""
    layout_path = os.path.join(os.path.dirname(__file__), "Enamine_plate layout1.xlsx")
    if os.path.exists(layout_path):
        try:
            df = pd.read_excel(layout_path, header=None)
            return df
        except Exception as e:
            st.error(f"Error loading plate layout: {str(e)}")
            return None
    else:
        return None


# Page configuration
st.set_page_config(
    page_title="ITR ENAMINE LIBRARY Analysis",
    page_icon=APP_ICON,
    layout="wide",
    initial_sidebar_state="expanded"
)

# Auto-load ENAMINE Library on startup
if 'library_df' not in st.session_state or st.session_state.library_df is None:
    library_df = load_enamine_library()
    if library_df is not None:
        st.session_state.library_df = library_df
        st.session_state.library_loaded = True
    else:
        st.session_state.library_loaded = False

# Auto-load Plate Layout on startup
if 'plate_layout' not in st.session_state:
    layout_df = load_plate_layout()
    if layout_df is not None:
        st.session_state.plate_layout = layout_df
        st.session_state.layout_loaded = True
    else:
        st.session_state.layout_loaded = False

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
    
    # Show status of auto-loaded files
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("ENAMINE Library")
        if st.session_state.get('library_loaded', False):
            num_compounds = len(st.session_state.library_df)
            try:
                num_plates = st.session_state.library_df['Plate_ID'].nunique()
                st.success(f"✅ Loaded: {num_compounds:,} compounds across {num_plates} plates")
            except:
                st.success(f"✅ Loaded: {num_compounds:,} compounds")
        else:
            st.error("❌ Library file not found: Enamine_library.xlsx")
            st.caption("Please ensure the file is in the same folder as app.py")
    
    with col2:
        st.subheader("Plate Layout")
        if st.session_state.get('layout_loaded', False):
            st.success("✅ Loaded: Enamine_plate layout1.xlsx")
        else:
            st.error("❌ Layout file not found: Enamine_plate layout1.xlsx")
            st.caption("Please ensure the file is in the same folder as app.py")
    
    st.divider()
    
    # Only show Plate Data Files uploader
    st.subheader("Plate Data Files")
    st.write("Upload plate reader output files (*.xlsx)")
    plate_files = st.file_uploader(
        "Upload plate data files",
        type=["xlsx"],
        accept_multiple_files=True,
        key='plate_uploader',
        label_visibility="collapsed"
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


def render_histogram_section():
    """Render Task 1A: %Inhibition Histogram."""
    st.markdown("## 📈 Task 1A: %Inhibition Distribution")
    
    if st.session_state.analysis_df is None:
        st.session_state.analysis_df = aggregate_replicates(
            st.session_state.main_df.copy(),
            st.session_state.replicate_mode
        )
    
    df = st.session_state.analysis_df
    
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
    
    # Download options
    col1, col2 = st.columns(2)
    with col1:
        # Download plot
        img_bytes = pio.to_image(fig, format='png', width=1200, height=600)
        st.download_button(
            "📷 Download Plot (PNG)",
            data=img_bytes,
            file_name="inhibition_histogram.png",
            mime="image/png"
        )
    with col2:
        # Download data
        csv = df[['Plate', 'Well', 'Pct_Inhibition']].to_csv(index=False)
        st.download_button(
            "📊 Download Data (CSV)",
            data=csv,
            file_name="inhibition_data.csv",
            mime="text/csv"
        )


def render_normalized_histogram_section():
    """Render Task 1B: Size-Normalized Activity Histogram."""
    st.markdown("## 📈 Task 1B: Size-Normalized Activity")
    
    st.info("**Formula:** SPEI = (% Inhibition / 100) / (MW × 0.001)\n\n*Higher values indicate more potent compounds per unit mass - smaller molecules with same activity are preferred.*")
    
    df = st.session_state.analysis_df
    
    # Ensure metrics are calculated
    if 'SPEI' not in df.columns:
        df = calculate_metrics(df)
        st.session_state.analysis_df = df
    
    # Create histogram
    fig = plot_normalized_histogram(df, st.session_state.replicate_mode)
    st.plotly_chart(fig, use_container_width=True)
    
    # Statistics
    stats = calculate_summary_stats(df, 'SPEI')
    
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
    
    # Download options
    col1, col2 = st.columns(2)
    with col1:
        img_bytes = pio.to_image(fig, format='png', width=1200, height=600)
        st.download_button(
            "📷 Download Plot (PNG)",
            data=img_bytes,
            file_name="normalized_histogram.png",
            mime="image/png"
        )
    with col2:
        csv = df[['Plate', 'Well', 'Pct_Inhibition', 'MW', 'SPEI']].to_csv(index=False)
        st.download_button(
            "📊 Download Data (CSV)",
            data=csv,
            file_name="normalized_data.csv",
            mime="text/csv"
        )


def render_metrics_section():
    """Render SPEI & PPEI metrics overview."""
    st.markdown("## 🧮 Efficiency Metrics")
    
    st.info("**SPEI** = per_one / (MW × 0.001) — Size efficiency: Higher = more potent per unit mass\n\n**PPEI** = per_one / (TPSA × 0.01) — Polarity efficiency: Higher = better membrane penetration")
    
    df = st.session_state.analysis_df
    
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
    
    df = st.session_state.analysis_df
    
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
    
    # Create plot
    fig = plot_cartesian(df, top_percentile, st.session_state.replicate_mode, metric)
    st.plotly_chart(fig, use_container_width=True)
    
    # Top candidates table
    st.markdown("### 🏆 Top Candidates")
    
    df_ranked = get_top_candidates(df.copy(), top_percentile, metric)
    top_df = df_ranked[df_ranked['Is_Top']].sort_values('Rank_Score', ascending=False)
    
    display_cols = ['Plate', 'Well', 'catalog_number', 'Chemical_name', 
                    'Pct_Inhibition', 'MW', 'TPSA', 'SPEI', 'PPEI']
    available_cols = [c for c in display_cols if c in top_df.columns]
    
    st.dataframe(
        top_df[available_cols].head(50).round(2),
        use_container_width=True,
        height=400
    )
    
    # Download options
    col1, col2 = st.columns(2)
    
    with col1:
        img_bytes = pio.to_image(fig, format='png', width=1200, height=800)
        st.download_button(
            "📷 Download Plot (PNG)",
            data=img_bytes,
            file_name="cartesian_plot.png",
            mime="image/png"
        )
    
    with col2:
        csv = top_df[available_cols].to_csv(index=False)
        st.download_button(
            "📊 Download Top Candidates (CSV)",
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



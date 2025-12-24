"""
HTS Screening Data Analysis Streamlit App

A comprehensive application for analyzing High-Throughput Screening (HTS) data
from ENAMINE compound library plates.
"""

import io
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


# Page configuration
st.set_page_config(
    page_title="HTS Screening Analysis",
    page_icon=APP_ICON,
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for beautiful UI
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Crimson+Pro:wght@400;600;700&family=JetBrains+Mono:wght@400;500&display=swap');
    
    /* Main styling */
    .main .block-container {
        padding-top: 2rem;
        max-width: 1400px;
    }
    
    /* Headers */
    h1, h2, h3 {
        font-family: 'Crimson Pro', serif !important;
        font-weight: 700;
    }
    
    h1 {
        background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        font-size: 2.5rem !important;
    }
    
    /* Metric cards */
    .metric-card {
        background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
        border-radius: 12px;
        padding: 1.5rem;
        border-left: 4px solid #0f3460;
        box-shadow: 0 4px 6px rgba(0,0,0,0.05);
    }
    
    .metric-value {
        font-family: 'JetBrains Mono', monospace;
        font-size: 1.8rem;
        font-weight: 600;
        color: #0f3460;
    }
    
    .metric-label {
        font-family: 'Crimson Pro', serif;
        color: #6c757d;
        font-size: 0.9rem;
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }
    
    /* Status badges */
    .status-pass {
        background-color: #d4edda;
        color: #155724;
        padding: 0.25rem 0.75rem;
        border-radius: 20px;
        font-weight: 500;
    }
    
    .status-fail {
        background-color: #f8d7da;
        color: #721c24;
        padding: 0.25rem 0.75rem;
        border-radius: 20px;
        font-weight: 500;
    }
    
    /* Section dividers */
    .section-divider {
        height: 3px;
        background: linear-gradient(90deg, #0f3460 0%, #e94560 50%, #0f3460 100%);
        margin: 2rem 0;
        border-radius: 2px;
    }
    
    /* File uploader */
    .stFileUploader > div > div {
        border: 2px dashed #0f3460;
        border-radius: 12px;
    }
    
    /* Buttons */
    .stButton > button {
        background: linear-gradient(135deg, #0f3460 0%, #16213e 100%);
        color: white;
        border: none;
        border-radius: 8px;
        padding: 0.75rem 2rem;
        font-family: 'Crimson Pro', serif;
        font-weight: 600;
        transition: all 0.3s ease;
    }
    
    .stButton > button:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(15, 52, 96, 0.3);
    }
    
    /* Expanders */
    .streamlit-expanderHeader {
        font-family: 'Crimson Pro', serif;
        font-weight: 600;
        color: #0f3460;
    }
    
    /* DataFrames */
    .dataframe {
        font-family: 'JetBrains Mono', monospace !important;
        font-size: 0.85rem;
    }
    
    /* Sidebar */
    .css-1d391kg {
        background: linear-gradient(180deg, #1a1a2e 0%, #16213e 100%);
    }
    
    /* Info boxes */
    .info-box {
        background: linear-gradient(135deg, #e3f2fd 0%, #bbdefb 100%);
        border-radius: 12px;
        padding: 1rem 1.5rem;
        border-left: 4px solid #1976d2;
        margin: 1rem 0;
    }
    
    /* Footer */
    .footer {
        text-align: center;
        padding: 2rem;
        color: #6c757d;
        font-family: 'Crimson Pro', serif;
    }
</style>
""", unsafe_allow_html=True)


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
    st.markdown("## 📁 Data Upload")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### ENAMINE Library File")
        st.caption("Upload the compound library (Enamine_library.xlsx)")
        library_file = st.file_uploader(
            "Upload Library",
            type=['xlsx'],
            key='library_uploader',
            label_visibility='collapsed'
        )
        
        if library_file:
            if st.session_state.library_df is None:
                with st.spinner("Loading library..."):
                    st.session_state.library_df = load_library(library_file)
                st.success(f"✅ Loaded {len(st.session_state.library_df):,} compounds")
            else:
                st.info(f"📚 Library loaded: {len(st.session_state.library_df):,} compounds")
    
    with col2:
        st.markdown("### Plate Data Files")
        st.caption("Upload plate reader output files (*.xlsx)")
        plate_files = st.file_uploader(
            "Upload Plates",
            type=['xlsx'],
            accept_multiple_files=True,
            key='plate_uploader',
            label_visibility='collapsed'
        )
        
        if plate_files:
            st.info(f"📊 {len(plate_files)} plate files selected")
    
    st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)
    
    # Alternative: Load saved data
    st.markdown("### 💾 Or Load Saved Session")
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
    
    # Process button
    if plate_files and st.session_state.library_df is not None:
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
        st.markdown(f"""
        <div class="metric-card">
            <div class="metric-label">Total Plates</div>
            <div class="metric-value">{total_plates}</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown(f"""
        <div class="metric-card">
            <div class="metric-label">Passed QC</div>
            <div class="metric-value" style="color: #28a745;">{passed_plates}</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        st.markdown(f"""
        <div class="metric-card">
            <div class="metric-label">Failed QC</div>
            <div class="metric-value" style="color: #dc3545;">{failed_plates}</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col4:
        st.markdown(f"""
        <div class="metric-card">
            <div class="metric-label">Avg Z' Factor</div>
            <div class="metric-value">{avg_zprime:.3f}</div>
        </div>
        """, unsafe_allow_html=True)
    
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
    
    st.markdown("""
    <div class="info-box">
        Choose how to handle replicate measurements for downstream analysis.
        Different approaches may be appropriate depending on your data quality and analysis goals.
    </div>
    """, unsafe_allow_html=True)
    
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
    
    st.markdown("""
    <div class="info-box">
        <strong>Formula:</strong> SPEI = (% Inhibition / 100) / (MW × 0.001)<br>
        <em>Higher values indicate more potent compounds per unit mass - smaller molecules with same activity are preferred.</em>
    </div>
    """, unsafe_allow_html=True)
    
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
    
    st.markdown("""
    <div class="info-box">
        <strong>SPEI</strong> = per_one / (MW × 0.001) — Size efficiency: Higher = more potent per unit mass<br>
        <strong>PPEI</strong> = per_one / (TPSA × 0.01) — Polarity efficiency: Higher = better membrane penetration
    </div>
    """, unsafe_allow_html=True)
    
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
    
    st.markdown("""
    <div class="info-box">
        Compounds in the <strong>northeast corner</strong> (high SPEI + high PPEI) are the best drug candidates:
        they're small, non-polar, and highly active.
    </div>
    """, unsafe_allow_html=True)
    
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
        st.markdown("""
        <div style="text-align: center; font-size: 0.8rem; color: #6c757d;">
            HTS Screening Analysis v1.0<br>
            Drug Discovery Pipeline
        </div>
        """, unsafe_allow_html=True)
        
        return selected


def main():
    """Main application entry point."""
    initialize_session_state()
    
    # Render sidebar and get selection
    selected = render_sidebar()
    
    # Title
    st.markdown(f"# {APP_TITLE}")
    st.markdown("*Antibiotic Discovery • ENAMINE Library • 384-well Plate Analysis*")
    st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)
    
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
    st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)
    st.markdown("""
    <div class="footer">
        Built for Antibiotic Discovery Research • Inspired by Stokes 2020 Halicin Paper
    </div>
    """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()



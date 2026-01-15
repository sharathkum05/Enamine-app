"""
Plotting functions for HTS Screening App
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import COLOR_REP1, COLOR_REP2, COLOR_HIGHLIGHT, COLOR_NORMAL


def plot_inhibition_histogram(df: pd.DataFrame, replicate_mode: str) -> go.Figure:
    """
    Create histogram of %Inhibition.
    
    For "both" mode: overlay semi-transparent histograms for Rep1 and Rep2
    
    Args:
        df: DataFrame with Pct_Inhibition and Rep columns
        replicate_mode: One of "average", "rep1", "rep2", "both"
    
    Returns:
        Plotly Figure
    """
    fig = go.Figure()
    
    if replicate_mode == "both":
        # Overlay histograms for both replicates
        rep1_data = df[df['Rep'] == 1]['Pct_Inhibition'].dropna()
        rep2_data = df[df['Rep'] == 2]['Pct_Inhibition'].dropna()
        
        fig.add_trace(go.Histogram(
            x=rep1_data,
            name='Replicate 1',
            opacity=0.6,
            marker_color=COLOR_REP1,
            nbinsx=50
        ))
        
        fig.add_trace(go.Histogram(
            x=rep2_data,
            name='Replicate 2',
            opacity=0.6,
            marker_color=COLOR_REP2,
            nbinsx=50
        ))
        
        fig.update_layout(barmode='overlay')
    else:
        # Single histogram
        data = df['Pct_Inhibition'].dropna()
        
        fig.add_trace(go.Histogram(
            x=data,
            name='%Inhibition',
            marker_color=COLOR_REP1,
            nbinsx=50
        ))
    
    fig.update_layout(
        title={
            'text': 'Distribution of %Inhibition',
            'font': {'size': 24, 'family': 'Crimson Pro, serif'}
        },
        xaxis_title='%Inhibition',
        yaxis_title='Count',
        template='plotly_white',
        font=dict(family='Crimson Pro, serif', size=14),
        showlegend=True if replicate_mode == "both" else False,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99
        ),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(248,249,250,1)'
    )
    
    # Add gridlines
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='rgba(0,0,0,0.1)')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='rgba(0,0,0,0.1)')
    
    return fig


def plot_normalized_histogram(df: pd.DataFrame, replicate_mode: str) -> go.Figure:
    """
    Create histogram of Per_one / MW(kDa) (equivalent to SPEI).
    
    Args:
        df: DataFrame with SPEI and Rep columns
        replicate_mode: One of "average", "rep1", "rep2", "both"
    
    Returns:
        Plotly Figure
    """
    fig = go.Figure()
    
    if replicate_mode == "both":
        rep1_data = df[df['Rep'] == 1]['SPEI'].dropna()
        rep2_data = df[df['Rep'] == 2]['SPEI'].dropna()
        
        fig.add_trace(go.Histogram(
            x=rep1_data,
            name='Replicate 1',
            opacity=0.6,
            marker_color=COLOR_REP1,
            nbinsx=50
        ))
        
        fig.add_trace(go.Histogram(
            x=rep2_data,
            name='Replicate 2',
            opacity=0.6,
            marker_color=COLOR_REP2,
            nbinsx=50
        ))
        
        fig.update_layout(barmode='overlay')
    else:
        data = df['SPEI'].dropna()
        
        fig.add_trace(go.Histogram(
            x=data,
            name='SPEI',
            marker_color='#2ecc71',
            nbinsx=50
        ))
    
    fig.update_layout(
        title={
            'text': 'Size-Normalized Activity (SPEI)',
            'font': {'size': 24, 'family': 'Crimson Pro, serif'}
        },
        xaxis_title='SPEI = (% Inhibition / 100) / (MW × 0.001)',
        yaxis_title='Count',
        template='plotly_white',
        font=dict(family='Crimson Pro, serif', size=14),
        showlegend=True if replicate_mode == "both" else False,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99
        ),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(248,249,250,1)'
    )
    
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='rgba(0,0,0,0.1)')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='rgba(0,0,0,0.1)')
    
    return fig


def plot_cartesian(df: pd.DataFrame, top_percentile: float, 
                   replicate_mode: str, metric: str = "combined") -> go.Figure:
    """
    Create PPEI vs SPEI scatter plot.
    
    Features:
        - X-axis: PPEI, Y-axis: SPEI
        - Highlight top candidates
        - For "both" mode: different colors for Rep1 vs Rep2
        - Hover shows compound details
    
    Args:
        df: DataFrame with SPEI, PPEI, and compound info
        top_percentile: Top X% to highlight
        replicate_mode: One of "average", "rep1", "rep2", "both"
        metric: Metric for ranking ("combined", "spei", "ppei")
    
    Returns:
        Plotly Figure
    """
    fig = go.Figure()
    
    # Add rank score for highlighting
    from .metrics import get_top_candidates
    df_ranked = get_top_candidates(df.copy(), top_percentile, metric)
    
    # Create hover text
    def create_hover_text(row):
        catalog = row.get('catalog_number', 'N/A')
        name = row.get('Chemical_name', 'N/A')
        if pd.isna(name) or name == 'N/A':
            name = 'N/A'
        elif len(str(name)) > 40:
            name = str(name)[:40] + '...'
        
        return (
            f"<b>Compound:</b> {catalog}<br>"
            f"<b>Name:</b> {name}<br>"
            f"<b>%Inhibition:</b> {row['Pct_Inhibition']:.1f}%<br>"
            f"<b>SPEI:</b> {row['SPEI']:.2f}<br>"
            f"<b>PPEI:</b> {row['PPEI']:.2f}<br>"
            f"<b>MW:</b> {row.get('MW', 'N/A')}<br>"
            f"<b>TPSA:</b> {row.get('TPSA', 'N/A')}"
        )
    
    df_ranked['hover_text'] = df_ranked.apply(create_hover_text, axis=1)
    
    if replicate_mode == "both":
        # Plot both replicates with different colors
        for rep, color, name in [(1, COLOR_REP1, 'Rep 1'), (2, COLOR_REP2, 'Rep 2')]:
            rep_data = df_ranked[df_ranked['Rep'] == rep]
            
            # Non-top candidates
            non_top = rep_data[~rep_data['Is_Top']]
            fig.add_trace(go.Scatter(
                x=non_top['PPEI'],
                y=non_top['SPEI'],
                mode='markers',
                name=f'{name}',
                marker=dict(
                    color=color,
                    size=6,
                    opacity=0.5
                ),
                hovertext=non_top['hover_text'],
                hoverinfo='text'
            ))
            
            # Top candidates
            top = rep_data[rep_data['Is_Top']]
            if len(top) > 0:
                fig.add_trace(go.Scatter(
                    x=top['PPEI'],
                    y=top['SPEI'],
                    mode='markers',
                    name=f'{name} - Top {top_percentile}%',
                    marker=dict(
                        color=COLOR_HIGHLIGHT,
                        size=10,
                        symbol='star',
                        line=dict(width=1, color='white')
                    ),
                    hovertext=top['hover_text'],
                    hoverinfo='text'
                ))
    else:
        # Single mode
        # Non-top candidates
        non_top = df_ranked[~df_ranked['Is_Top']]
        fig.add_trace(go.Scatter(
            x=non_top['PPEI'],
            y=non_top['SPEI'],
            mode='markers',
            name='Compounds',
            marker=dict(
                color=COLOR_NORMAL,
                size=6,
                opacity=0.6
            ),
            hovertext=non_top['hover_text'],
            hoverinfo='text'
        ))
        
        # Top candidates
        top = df_ranked[df_ranked['Is_Top']]
        if len(top) > 0:
            fig.add_trace(go.Scatter(
                x=top['PPEI'],
                y=top['SPEI'],
                mode='markers',
                name=f'Top {top_percentile}%',
                marker=dict(
                    color=COLOR_HIGHLIGHT,
                    size=10,
                    symbol='star',
                    line=dict(width=1, color='white')
                ),
                hovertext=top['hover_text'],
                hoverinfo='text'
            ))
    
    # Add quadrant annotations
    fig.add_annotation(
        x=0.95, y=0.95,
        xref='paper', yref='paper',
        text='⭐ Best Candidates',
        showarrow=False,
        font=dict(size=14, color='#e74c3c'),
        bgcolor='rgba(255,255,255,0.8)',
        borderpad=4
    )
    
    fig.update_layout(
        title={
            'text': 'SPEI vs PPEI: Drug Efficiency Plot',
            'font': {'size': 24, 'family': 'Crimson Pro, serif'}
        },
        xaxis_title='PPEI (Polarity Efficiency) - Higher = Better Membrane Penetration',
        yaxis_title='SPEI (Size Efficiency) - Higher = More Potent per Unit Mass',
        template='plotly_white',
        font=dict(family='Crimson Pro, serif', size=14),
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01,
            bgcolor='rgba(255,255,255,0.9)'
        ),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(248,249,250,1)',
        hovermode='closest'
    )
    
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='rgba(0,0,0,0.1)')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='rgba(0,0,0,0.1)')
    
    return fig


def plot_qc_heatmap(qc_df: pd.DataFrame) -> go.Figure:
    """
    Create a heatmap of Z' values for quality control visualization.
    
    Args:
        qc_df: DataFrame with Plate, Rep, Z_Prime columns
    
    Returns:
        Plotly Figure
    """
    # Pivot for heatmap - use pivot_table to handle duplicates with mean aggregation
    pivot = qc_df.pivot_table(index='Plate', columns='Rep', values='Z_Prime', aggfunc='mean')
    
    fig = go.Figure(data=go.Heatmap(
        z=pivot.values,
        x=[f'Rep {c}' for c in pivot.columns],
        y=pivot.index,
        colorscale=[
            [0, '#e74c3c'],      # Red for low Z'
            [0.5, '#f39c12'],    # Orange for borderline
            [1, '#2ecc71']       # Green for good Z'
        ],
        zmin=0,
        zmax=1,
        colorbar=dict(title="Z' Factor")
    ))
    
    fig.update_layout(
        title={
            'text': "Z' Factor Quality Control Heatmap",
            'font': {'size': 20, 'family': 'Crimson Pro, serif'}
        },
        xaxis_title='Replicate',
        yaxis_title='Plate',
        template='plotly_white',
        font=dict(family='Crimson Pro, serif', size=14)
    )
    
    return fig


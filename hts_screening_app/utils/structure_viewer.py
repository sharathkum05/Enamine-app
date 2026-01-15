"""
Chemical structure viewer using RDKit.
RDKit is optional - structure viewer is disabled if RDKit is not installed.
"""

import io
import pandas as pd

# Try to import RDKit, but make it optional
RDKIT_AVAILABLE = False
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem import rdMolDescriptors
    RDKIT_AVAILABLE = True
except ImportError:
    Chem = None
    Draw = None
    rdMolDescriptors = None


def smiles_to_image(smiles, size=(350, 350)):
    """
    Convert a SMILES string to a PNG image.
    Returns None if RDKit is not available or conversion fails.
    """
    if not RDKIT_AVAILABLE:
        return None
    
    try:
        if not smiles or pd.isna(smiles):
            return None
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
            
        # Generate image
        img = Draw.MolToImage(mol, size=size)
        
        # Convert to bytes
        buffer = io.BytesIO()
        img.save(buffer, format='PNG')
        buffer.seek(0)
        return buffer
        
    except Exception as e:
        print(f"Error converting SMILES to image: {e}")
        return None


def get_molecule_info(smiles):
    """
    Get basic molecular information from SMILES.
    Returns None if RDKit is not available.
    """
    if not RDKIT_AVAILABLE:
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        return {
            'formula': rdMolDescriptors.CalcMolFormula(mol),
            'num_atoms': mol.GetNumAtoms(),
            'num_bonds': mol.GetNumBonds(),
            'num_rings': rdMolDescriptors.CalcNumRings(mol),
        }
    except Exception:
        return None


def is_rdkit_available():
    """Check if RDKit is available."""
    return RDKIT_AVAILABLE

"""
Chemical structure viewer using RDKit.
Converts SMILES strings to molecular structure images.
"""

import io
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw


def smiles_to_image(smiles, size=(350, 350)):
    """
    Convert a SMILES string to a PNG image.
    
    Args:
        smiles: SMILES string of the molecule
        size: Tuple of (width, height) for the image
        
    Returns:
        BytesIO buffer containing the PNG image, or None if conversion fails
    """
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
    
    Returns dict with molecular formula, etc.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
            
        from rdkit.Chem import Descriptors, rdMolDescriptors
        
        return {
            'formula': rdMolDescriptors.CalcMolFormula(mol),
            'num_atoms': mol.GetNumAtoms(),
            'num_bonds': mol.GetNumBonds(),
            'num_rings': rdMolDescriptors.CalcNumRings(mol),
        }
    except Exception:
        return None

"""
Molecular visualization and rendering.
"""
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

IPythonConsole.ipython_useSVG = True


class MoleculeRenderer:
    """Handles molecular visualization."""

    def __init__(self):
        self.render_styles = {}

    def render_2d(self, mol, title="", highlight_atoms=None):
        """
        Render 2D molecular structure.

        Args:
            mol: Molecule to render
            highlight_atoms: Atoms to highlight
            title: Name of molecule

        Returns:
            2D visualization
        """
        if mol is None:
            print(f"{title}: None")
            return None

        if title is "":
            title = f"{Chem.MolToSmiles(mol)}"

        mol_copy = Chem.Mol(mol)
        for atom in mol_copy.GetAtoms():
            atom.SetProp('atomNote', str(atom.GetIdx()))

        coef = mol.GetNumAtoms()
        img = Draw.MolToImage(mol_copy, legend=title, size=(30 * coef, 20 * coef),
                              highlightAtoms=highlight_atoms if highlight_atoms else [])
        return img

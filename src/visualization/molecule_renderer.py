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

    def render_2d(self, mol, title="", highlight_atoms=None, highlight_bonds =None, display_props=None):
        """
        Render 2D molecular structure.

        Args:
            mol: Molecule to render
            highlight_atoms: Atoms to highlight
            title: Name of molecule
            display_props: List of property names to display on atoms

        Returns:
            2D visualization
        """
        if mol is None:
            print(f"{title}: None")
            return None

        if title == "":
            title = f"{Chem.MolToSmiles(mol)}"

        mol_copy = Chem.Mol(mol)
        self._set_atom_display_notes(mol_copy, display_props)

        coef = mol.GetNumAtoms()

        img = Draw.MolToImage(mol_copy, legend=title, size=(max(30 * coef, 300), max(20 * coef, 200)),
                              highlightAtoms=highlight_atoms if highlight_atoms else [],
                              highlightBonds=highlight_bonds if highlight_bonds else [])
        return img

    def _set_atom_display_notes(self, mol, display_props=None):
        """
        Set atom display notes based on specified properties.
        """
        for atom in mol.GetAtoms():
            numbers = []
            booleans = []
            strings = []

            # Always include atom index
            numbers.append(str(atom.GetIdx()))

            if display_props:
                for prop in display_props:
                    if not atom.HasProp(prop):
                        continue

                    value = atom.GetProp(prop)

                    if value.lstrip('-').replace('.', '', 1).isdigit():
                        numbers.append(value)
                    elif value.lower() in ['true', 'false']:
                        if value.lower() == 'true':
                            booleans.append(prop)
                    else:
                        strings.append(value)


            display_parts = []
            if numbers:
                display_parts.append(":".join(numbers))
            if booleans:
                display_parts.append(":".join(booleans))
            if strings:
                display_parts.append(":".join(strings))


            note_text = " | ".join(display_parts)
            atom.SetProp('atomNote', note_text)


    def visualize_molecular_subgraph(self, mol, connection_atom_idx):
        """
        Визуализирует подграф молекулы вокруг точки соединения.
        """
        connection_atom = mol.GetAtomWithIdx(connection_atom_idx)
        replacement_group_smiles = connection_atom.GetProp('replacement_group')
        replacement_mol = Chem.MolFromSmiles(replacement_group_smiles)
        max_depth = replacement_mol.GetNumAtoms() if replacement_mol else 0

        self._print_visualization_info(mol, connection_atom_idx, replacement_group_smiles, max_depth)

        atoms_to_highlight, bonds_to_highlight = self._collect_subgraph_atoms_and_bonds(
            mol, connection_atom_idx, max_depth
        )

        self._create_subgraph_visualization(atoms_to_highlight, connection_atom_idx, mol, replacement_group_smiles)
        self._print_subgraph_info(atoms_to_highlight)

        return atoms_to_highlight

    def _print_visualization_info(self, mol, connection_atom_idx, replacement_group_smiles, max_depth):
        """Выводит информацию о визуализации."""
        print(f"Визуализация подграфа молекулы")
        print(f"Атом соединения: {mol.GetAtomWithIdx(connection_atom_idx).GetSymbol()}({connection_atom_idx})")
        print(f"Replacement group: {replacement_group_smiles}")
        print(f"Максимальная глубина обхода: {max_depth}")

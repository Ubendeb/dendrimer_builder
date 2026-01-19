"""
Format conversion utilities for molecular data.
"""
from rdkit import Chem


class FormatConverters:
    """Provides format conversion utilities between different molecular representations."""

    @staticmethod
    def json_to_metadata_mol(json_data):
        """
        Convert JSON representation to molecule with metadata.

        Args:
            json_data: Dictionary with molecular data in JSON format

        Returns:
            Molecule object with connection metadata
        """
        smiles = json_data["smiles"]
        mol = Chem.MolFromSmiles(smiles)

        # Freeze original atom indices
        for i, atom in enumerate(mol.GetAtoms()):
            atom.SetProp("original_index", str(i))

        # Set connection metadata
        connection_atoms = json_data["connection_atoms"]
        replacement_groups = json_data["replacement_groups"]

        for i, atom_idx in enumerate(connection_atoms):
            atom = mol.GetAtomWithIdx(atom_idx)
            atom.SetProp("is_connection", "true")
            atom.SetProp("replacement_group", replacement_groups[i])

            # Set connection type based on position
            if i == 0:
                atom.SetProp("connection_type", "to_core")
            else:
                atom.SetProp("connection_type", "to_branch")

        return mol

    @staticmethod
    def metadata_mol_to_json(mol):
        """
        Convert molecule with metadata to JSON representation.

        Args:
            mol: Molecule object with connection metadata

        Returns:
            Dictionary with molecular data in JSON format
        """
        smiles = Chem.MolToSmiles(mol)

        # Extract connection information from metadata
        connection_atoms = []
        replacement_groups = []

        for atom in mol.GetAtoms():
            if atom.HasProp("is_connection") and atom.GetProp("is_connection") == "true":
                connection_atoms.append(atom.GetIdx())
                replacement_groups.append(atom.GetProp("replacement_group"))

        return {
            "name": "",  # Name would need to be provided separately
            "smiles": smiles,
            "connection_atoms": connection_atoms,
            "replacement_groups": replacement_groups,
            "num_atoms": mol.GetNumAtoms()
        }
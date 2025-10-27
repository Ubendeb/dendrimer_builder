"""
Loading of molecular components and data.
"""

from rdkit import Chem


class DataLoader:
    """Handles loading of molecular data from various sources."""

    @staticmethod
    def load_smiles(smiles_string):
        """
        Load molecule from SMILES string.

        Args:
            smiles_string: SMILES representation

        Returns:
            RDKit Mol object or None if invalid
        """
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles_string}")
        return mol

    @staticmethod
    def load_file(file_path, file_format=None):
        """
        Load molecule from file.

        Args:
            file_path: Path to molecular file
            file_format: File format (auto-detected if None)

        Returns:
            RDKit Mol object or None if error
        """
        if file_format is None:
            file_format = file_path.split('.')[-1].lower()

        try:
            if file_format in ['sdf', 'sd']:
                supplier = Chem.SDMolSupplier(file_path)
                mols = [mol for mol in supplier if mol is not None]
                return mols[0] if mols else None
            elif file_format in ['mol', 'mol2']:
                return Chem.MolFromMolFile(file_path)
            elif file_format == 'pdb':
                return Chem.MolFromPDBFile(file_path)
            else:
                raise ValueError(f"Unsupported file format: {file_format}")
        except Exception as e:
            raise IOError(f"Error loading file {file_path}: {e}")

    @staticmethod
    def load_fragment_library(library_path):
        """
        Load library of molecular fragments from JSON file.

        Args:
            library_path: Path to fragment library JSON file

        Returns:
            Dictionary of fragments
        """
        try:
            with open(library_path, 'r', encoding='utf-8') as f:
                fragments_data = json.load(f)

            # Convert SMILES to mol objects
            for frag_id, frag_data in fragments_data.items():
                if 'smiles' in frag_data:
                    frag_data['mol'] = Chem.MolFromSmiles(frag_data['smiles'])

            return fragments_data
        except Exception as e:
            raise IOError(f"Error loading fragment library {library_path}: {e}")
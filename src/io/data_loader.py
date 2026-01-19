"""
Loading of molecular components and data.
"""

from rdkit import Chem
import json


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
    def load_fragment_data(filename):
        """
        Load fragments data from JSON file.

        Args:
            filename: Path to JSON file

        Returns:
            List of fragment data for builder
        """
        with open(filename, 'r', encoding='utf-8') as f:
            loaded_data = json.load(f)

        fragments_data = []
        for frag_id, frag_data in loaded_data.items():
            mol = Chem.MolFromSmiles(frag_data['smiles'])
            if mol:
                fragment_data = {
                    'name': frag_id,
                    'smiles': frag_data['smiles'],
                    'connection_atoms': frag_data.get('connection_points', []),
                    'replacement_groups': frag_data.get('replacement_groups', []),
                    'num_atoms': frag_data['num_atoms'],
                    'mol': mol
                }
                fragments_data.append(fragment_data)
        return fragments_data

    @staticmethod
    def load_fragment_library(file_path, fragment_manager):
        """
        Load fragment library into FragmentManager.

        Args:
            file_path: Path to fragment library JSON file
            fragment_manager: FragmentManager instance to populate

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                loaded_data = json.load(f)

            fragment_manager.fragments.clear()
            loaded_count = 0

            for frag_id, frag_data in loaded_data.items():
                mol = Chem.MolFromSmiles(frag_data['smiles'])
                if mol:
                    fragment_manager.fragments[frag_id] = {
                        'molecule': mol,
                        'type': frag_data['type'],
                        'connection_points': frag_data['connection_points'],
                        'replacement_groups': frag_data['replacement_groups'],
                        'smiles': frag_data['smiles'],
                        'num_atoms': frag_data['num_atoms']
                    }
                    loaded_count += 1
                else:
                    print(f"Warning: Could not parse SMILES for fragment {frag_id}")

            print(f"Fragment library loaded from: {file_path}")
            print(f"Fragments loaded: {loaded_count}")
            return True

        except FileNotFoundError:
            print(f"File not found: {file_path}")
            return False
        except json.JSONDecodeError as e:
            print(f"Invalid JSON format in file {file_path}: {e}")
            return False
        except Exception as e:
            print(f"Error loading file: {e}")
            return False
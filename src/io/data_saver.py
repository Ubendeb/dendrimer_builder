"""
Saving of results and molecular structures.
"""

from rdkit import Chem
import json


class DataSaver:
    """Handles saving of molecular data to various formats."""

    @staticmethod
    def save_molecule(molecule, file_path, file_format=None):
        """
        Save molecule to file.

        Args:
            molecule: RDKit Mol object to save
            file_path: Output file path
            file_format: Output format (auto-detected if None)
        """
        if file_format is None:
            file_format = file_path.split('.')[-1].lower()

        try:
            if file_format in ['sdf', 'sd']:
                writer = Chem.SDWriter(file_path)
                writer.write(molecule)
                writer.close()
            elif file_format in ['mol']:
                Chem.MolToMolFile(molecule, file_path)
            elif file_format == 'pdb':
                Chem.MolToPDBFile(molecule, file_path)
            elif file_format == 'smiles':
                smiles = Chem.MolToSmiles(molecule)
                with open(file_path, 'w') as f:
                    f.write(smiles)
            else:
                raise ValueError(f"Unsupported output format: {file_format}")

            print(f"Molecule saved to: {file_path}")
        except Exception as e:
            raise IOError(f"Error saving molecule to {file_path}: {e}")

    @staticmethod
    def serialize_fragments_data(fragments_data):
        """
        Serialize fragments data to dictionary format.

        Args:
            fragments_data: List of fragment data from builder

        Returns:
            Dictionary suitable for JSON serialization
        """
        serializable_data = {}
        for frag in fragments_data:
            frag_id = frag['name']
            serializable_data[frag_id] = {
                'smiles': frag['smiles'],
                'type': frag.get('type', 'custom'),
                'connection_points': frag['connection_atoms'],
                'replacement_groups': frag['replacement_groups'],
                'num_atoms': frag['num_atoms']
            }
        return serializable_data

    @staticmethod
    def save_fragment_data(fragments_data, filename):
        """
        Save fragments data to JSON file.

        Args:
            fragments_data: List of fragment data from builder
            filename: Output file path
        """
        try:
            serializable_data = DataSaver.serialize_fragments_data(fragments_data)

            with open(filename, 'w', encoding='utf-8') as f:
                json.dump(serializable_data, f, indent=2, ensure_ascii=False)

            print(f"Data saved to file: {filename}")
            print(f"Fragments saved: {len(serializable_data)}")
        except Exception as e:
            raise IOError(f"Error saving fragment data to {filename}: {e}")

    @staticmethod
    def save_fragment_library(fragment_manager, file_path='fragments_library.json'):
        """
        Save fragments library to JSON file.

        Args:
            fragment_manager: FragmentManager instance
            file_path: Output file path
        """
        try:
            serializable_data = {}
            for frag_id, frag_data in fragment_manager.fragments.items():
                serializable_data[frag_id] = {
                    'smiles': frag_data['smiles'],
                    'type': frag_data['type'],
                    'connection_points': frag_data['connection_points'],
                    'replacement_groups': frag_data['replacement_groups'],
                    'num_atoms': frag_data['num_atoms']
                }

            with open(file_path, 'w', encoding='utf-8') as f:
                json.dump(serializable_data, f, indent=2, ensure_ascii=False)

            print(f"Fragment library saved to: {file_path}")
            return True
        except Exception as e:
            print(f"Error saving fragment library to {file_path}: {e}")
            return False

    @staticmethod
    def save_dendrimer_report(dendrimer, file_path):
        """
        Save comprehensive dendrimer report.

        Args:
            dendrimer: Dendrimer structure
            file_path: Output file path
        """
        try:
            report = {
                'generation': dendrimer.generation,
                'total_atoms': dendrimer.get_num_atoms(),
                'total_bonds': dendrimer.get_num_bonds(),
                'molecular_weight': dendrimer.get_molecular_weight(),
                'formula': dendrimer.get_molecular_formula(),
            }

            with open(file_path, 'w', encoding='utf-8') as f:
                json.dump(report, f, indent=2, ensure_ascii=False)

            print(f"Dendrimer report saved to: {file_path}")
        except Exception as e:
            raise IOError(f"Error saving dendrimer report to {file_path}: {e}")

    @staticmethod
    def save_construction_log(builder, file_path):
        """
        Save construction log.

        Args:
            builder: Builder instance
            file_path: Output file path
        """
        try:
            log_data = {
                'construction_steps': builder.construction_history,
                'fragments_used': [
                    {
                        'name': frag['name'],
                        'connection_atoms': frag['connection_atoms'],
                        'replacement_groups': frag['replacement_groups']
                    }
                    for frag in builder.get_fragments_data()
                ]
            }

            with open(file_path, 'w', encoding='utf-8') as f:
                json.dump(log_data, f, indent=2, ensure_ascii=False)

            print(f"Construction log saved to: {file_path}")
        except Exception as e:
            raise IOError(f"Error saving construction log to {file_path}: {e}")
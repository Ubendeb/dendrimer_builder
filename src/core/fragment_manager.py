"""
Management of molecular fragments and their properties.
"""

import json
from rdkit import Chem


class FragmentManager:
    """Manages molecular fragments and their properties."""

    def __init__(self):
        self.fragments = {}
        self.fragment_counter = 0

    def add_fragment(self, fragment, fragment_type="core", name=None, connection_atoms=None, replacement_groups=None):
        """Add a new fragment to the manager."""
        if fragment is None:
            raise ValueError("Fragment molecule cannot be None")

        fragment_id = f"{fragment_type}_{self.fragment_counter}" if not name else name

        self.fragments[fragment_id] = {
            'molecule': fragment,
            'type': fragment_type,
            'connection_points': connection_atoms or [],
            'replacement_groups': replacement_groups or [],
            'smiles': Chem.MolToSmiles(fragment),
            'num_atoms': fragment.GetNumAtoms()
        }
        self.fragment_counter += 1
        return fragment_id

    def get_fragment(self, fragment_id):
        """Retrieve a fragment by ID."""
        return self.fragments.get(fragment_id)

    def list_fragments(self):
        """List all available fragments."""
        return list(self.fragments.keys())

    def get_fragment_by_type(self, fragment_type):
        """Get all fragments of specific type."""
        return {fid: frag for fid, frag in self.fragments.items()
                if frag['type'] == fragment_type}

    def remove_fragment(self, fragment_id):
        """Remove a fragment from manager."""
        if fragment_id in self.fragments:
            del self.fragments[fragment_id]
            return True
        return False

    def save_to_file(self, filename='fragments_library.json'):
        """Save fragments library to JSON file."""
        serializable_data = {}
        for frag_id, frag_data in self.fragments.items():
            serializable_data[frag_id] = {
                'smiles': frag_data['smiles'],
                'type': frag_data['type'],
                'connection_points': frag_data['connection_points'],
                'replacement_groups': frag_data['replacement_groups'],
                'num_atoms': frag_data['num_atoms']
            }

        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(serializable_data, f, indent=2, ensure_ascii=False)

        print(f"Fragment library saved to: {filename}")
        return True

    def load_from_file(self, filename='fragments_library.json'):
        """Load fragments library from JSON file."""
        try:
            with open(filename, 'r', encoding='utf-8') as f:
                loaded_data = json.load(f)

            self.fragments.clear()
            loaded_count = 0

            for frag_id, frag_data in loaded_data.items():
                mol = Chem.MolFromSmiles(frag_data['smiles'])
                if mol:
                    self.fragments[frag_id] = {
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

            print(f"Fragment library loaded from: {filename}")
            print(f"Fragments loaded: {loaded_count}")
            return True

        except FileNotFoundError:
            print(f"File not found: {filename}")
            return False
        except json.JSONDecodeError as e:
            print(f"Invalid JSON format in file {filename}: {e}")
            return False
        except Exception as e:
            print(f"Error loading file: {e}")
            return False

    def get_fragment_statistics(self):
        """Get statistics about fragments in the library."""
        stats = {
            'total_fragments': len(self.fragments),
            'fragments_by_type': {},
            'total_atoms': 0,
            'fragments_with_connections': 0
        }

        for frag_data in self.fragments.values():
            # Count by type
            frag_type = frag_data['type']
            stats['fragments_by_type'][frag_type] = stats['fragments_by_type'].get(frag_type, 0) + 1

            # Total atoms
            stats['total_atoms'] += frag_data['num_atoms']

            # Fragments with connection points
            if frag_data['connection_points']:
                stats['fragments_with_connections'] += 1

        return stats

    def clear_library(self):
        """Clear all fragments from the library."""
        fragment_count = len(self.fragments)
        self.fragments.clear()
        self.fragment_counter = 0
        print(f"Library cleared. Removed {fragment_count} fragments.")
        return fragment_count

    def update_fragment(self, fragment_id, **kwargs):
        """Update fragment properties."""
        if fragment_id not in self.fragments:
            print(f"Fragment {fragment_id} not found")
            return False

        allowed_fields = ['type', 'connection_points', 'replacement_groups']
        for field, value in kwargs.items():
            if field in allowed_fields:
                self.fragments[fragment_id][field] = value
            else:
                print(f"Field {field} cannot be updated")

        return True
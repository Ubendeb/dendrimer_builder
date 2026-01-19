"""
Management of molecular fragments and their properties.
"""

import json
from rdkit import Chem

from src import DataLoader, DataSaver


class FragmentManager:
    """Manages molecular fragments and their properties."""

    def __init__(self):
        self.fragments = {}
        self.fragment_counter = 0

    def add_fragment(self, fragment, fragment_type="core", name=None,
                     connection_atoms=None, replacement_groups=None):
        """Add a new fragment to the manager with detailed connection information."""
        if fragment is None:
            raise ValueError("Fragment molecule cannot be None")

        fragment_id = self._generate_fragment_id(fragment_type, name)
        connection_details = self._extract_connection_details(
            fragment, connection_atoms, replacement_groups
        )

        self.fragments[fragment_id] = self._create_fragment_data(
            fragment, fragment_type, connection_atoms,
            replacement_groups, connection_details
        )

        self.fragment_counter += 1
        return fragment_id

    def _create_fragment_data(self, fragment, fragment_type, connection_atoms,
                              replacement_groups, connection_details):
        """Create the complete fragment data dictionary."""
        return {
            'molecule': fragment,
            'type': fragment_type,
            'connection_points': connection_atoms or [],
            'connection_details': connection_details,
            'replacement_groups': replacement_groups or [],
            'smiles': Chem.MolToSmiles(fragment),
            'num_atoms': fragment.GetNumAtoms(),
            'properties': self._extract_atom_properties(fragment)
        }

    def _extract_atom_properties(self, mol):
        """Extract properties from all atoms in the molecule."""
        properties = []
        for atom in mol.GetAtoms():
            atom_props = self._get_basic_atom_properties(atom)
            self._add_custom_atom_properties(atom, atom_props)
            properties.append(atom_props)
        return properties

    def _get_basic_atom_properties(self, atom):
        """Get basic chemical properties of an atom."""
        return {
            'index': atom.GetIdx(),
            'symbol': atom.GetSymbol(),
            'charge': atom.GetFormalCharge(),
            'hybridization': str(atom.GetHybridization()),
            'is_aromatic': atom.GetIsAromatic(),
            'degree': atom.GetDegree(),
            'explicit_valence': atom.GetExplicitValence()
        }

    def _add_custom_atom_properties(self, atom, atom_props):
        """Add custom properties (like is_connection, connection_type) to atom properties."""
        custom_props = ["is_connection", "connection_type", "replacement_group"]
        for prop in custom_props:
            if atom.HasProp(prop):
                atom_props[prop] = atom.GetProp(prop)

    def _generate_fragment_id(self, fragment_type, name):
        """Generate a unique fragment ID."""
        if name:
            return name
        return f"{fragment_type}_{self.fragment_counter}"

    def _extract_connection_details(self, fragment, connection_atoms, replacement_groups):
        """Extract detailed information about connection points."""
        if not connection_atoms:
            return []

        connection_details = []

        for i, atom_info in enumerate(self._get_atom_info_pairs(connection_atoms, replacement_groups)):
            atom_idx, replacement = atom_info
            atom = fragment.GetAtomWithIdx(atom_idx)

            connection_details.append(self._create_connection_detail(
                atom, atom_idx, replacement, i
            ))

        return connection_details

    def _create_connection_detail(self, atom, atom_idx, replacement, position):
        """Create a dictionary with detailed connection information for a single atom."""
        connection_type = self._get_atom_property(atom, "connection_type")

        return {
            'atom_index': atom_idx,
            'replacement_group': replacement,
            'is_connection': self._get_atom_property(atom, "is_connection", default="false"),
            'connection_type': connection_type or self._get_default_connection_type(position),
            'atom_symbol': atom.GetSymbol(),
            'atom_charge': atom.GetFormalCharge()
        }

    def _get_default_connection_type(self, position):
        """Get default connection type based on position in connection list. TODO  может приводить к ошибкам"""
        return "to_core" if position == 0 else "to_branch"

    def _get_atom_property(self, atom, prop_name, default=None):
        """Get a property from atom if it exists."""
        if atom.HasProp(prop_name):
            return atom.GetProp(prop_name)
        return default

    def _get_atom_info_pairs(self, connection_atoms, replacement_groups):
        """Generate pairs of atom indices and their replacement groups."""
        if connection_atoms and replacement_groups:
            return zip(connection_atoms, replacement_groups)
        elif connection_atoms:
            return [(idx, None) for idx in connection_atoms]
        return []


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
        return DataSaver.save_fragment_library(self, filename)

    def load_from_file(self, filename='fragments_library.json'):
        """Load fragments library from JSON file."""
        return DataLoader.load_fragment_library(filename, self)

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
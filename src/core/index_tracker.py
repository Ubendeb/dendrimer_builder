"""
Tracking and management of atom indices during dendrimer construction.
"""


class IndexTracker:
    """Tracks atom indices during molecular assembly."""

    def __init__(self):
        self.atom_mapping = {}
        self.current_index = 1

    def add_atoms(self, fragment_id, atom_indices):
        """Add atom indices for a fragment."""
        self.atom_mapping[fragment_id] = atom_indices
        self.current_index = max(max(atom_indices) + 1, self.current_index)

    def get_next_indices(self, count):
        """Get the next available atom indices."""
        indices = list(range(self.current_index, self.current_index + count))
        self.current_index += count
        return indices

    def get_fragment_indices(self, fragment_id):
        """Get atom indices for a specific fragment."""
        return self.atom_mapping.get(fragment_id, [])
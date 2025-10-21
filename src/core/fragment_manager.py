"""
Management of molecular fragments and their properties.
"""


class FragmentManager:
    """Manages molecular fragments and their properties."""

    def __init__(self):
        self.fragments = {}
        self.fragment_counter = 0

    def add_fragment(self, fragment, fragment_type="core"):
        """Add a new fragment to the manager."""
        fragment_id = f"{fragment_type}_{self.fragment_counter}"
        self.fragments[fragment_id] = {
            'molecule': fragment,
            'type': fragment_type,
            'connection_points': []
        }
        self.fragment_counter += 1
        return fragment_id

    def get_fragment(self, fragment_id):
        """Retrieve a fragment by ID."""
        return self.fragments.get(fragment_id)

    def list_fragments(self):
        """List all available fragments."""
        return list(self.fragments.keys())
"""
Interface for fragment input and management.
"""


class FragmentInterface:
    """Provides interface for molecular fragment input and management."""

    def __init__(self):
        self.fragment_library = {}

    def input_fragment_interactive(self):
        """
        Interactive fragment input.

        Returns:
            Fragment data
        """
        # Implementation for interactive input
        pass

    def import_fragment_file(self, file_path):
        """
        Import fragment from file.

        Args:
            file_path: Path to fragment file

        Returns:
            Fragment data
        """
        # Implementation for file import
        pass

    def define_connection_points(self, fragment, connection_atoms):
        """
        Define connection points for fragment.

        Args:
            fragment: Molecular fragment
            connection_atoms: List of connection atom indices

        Returns:
            Fragment with defined connection points
        """
        # Implementation for connection point definition
        pass
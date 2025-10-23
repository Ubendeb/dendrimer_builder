"""
Interface for dendrimer construction control.
"""


class ConstructionInterface:
    """Provides interface for controlling dendrimer construction."""

    def __init__(self, builder):
        self.builder = builder
        self.construction_parameters = {}

    def set_parameters(self, parameters):
        """
        Set construction parameters.

        Args:
            parameters: Dictionary of parameters
        """
        self.construction_parameters.update(parameters)

    def run_construction(self):
        """Execute dendrimer construction."""
        # Implementation for construction execution
        pass

    def preview_generation(self, generation_num):
        """
        Preview specific generation.

        Args:
            generation_num: Generation number to preview

        Returns:
            Preview data
        """
        raise NotImplementedError("Implementation for generation preview")

    def optimize_structure(self, dendrimer):
        """
        Optimize dendrimer structure.

        Args:
            dendrimer: Dendrimer to optimize

        Returns:
            Optimized dendrimer
        """
        raise NotImplementedError("Implementation for structure optimization")
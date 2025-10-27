"""
Abstract base class for molecular builders.
"""

from abc import ABC, abstractmethod


class BaseBuilder(ABC):
    """Abstract base class for molecular construction."""

    def __init__(self):
        self.components = {}
        self.construction_steps = []

    @abstractmethod
    def add_component(self, component_type, component):
        """Add a construction component."""
        raise NotImplementedError

    @abstractmethod
    def build(self):
        """Execute the construction process."""
        raise NotImplementedError

    def get_construction_steps(self):
        """Get list of construction steps."""
        return self.construction_steps.copy()

    def reset(self):
        """Reset builder to initial state."""
        self.components.clear()
        self.construction_steps.clear()
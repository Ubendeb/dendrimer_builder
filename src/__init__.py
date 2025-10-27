"""
Dendrimer Builder - A modular system for molecular dendrimer construction.
"""

__version__ = "1.0.0"
__author__ = "Your Name"
__description__ = "Modular dendrimer construction system"

from .core import FragmentManager, IndexTracker, StructureValidator, ChemistryOperations
from .builders import DendrimerBuilder
from .io import DataLoader, DataSaver, FormatConverters
from .visualization import MoleculeRenderer, ProgressTracker
from .interfaces import FragmentInterface, ConstructionInterface, FragmentDataManager

__all__ = [
    'FragmentManager', 'IndexTracker', 'StructureValidator', 'ChemistryOperations',
    'DendrimerBuilder',
    'DataLoader', 'DataSaver', 'FormatConverters',
    'MoleculeRenderer', 'ProgressTracker',
    'FragmentInterface', 'ConstructionInterface', 'FragmentDataManager'
]
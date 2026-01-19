"""
Tracking and management of atom indices during dendrimer construction.
"""

from rdkit import Chem
from typing import Dict, List, Optional, Union, Iterable


class IndexTracker:
    """Tracks atom indices during molecular assembly with flexible filtering."""

    def __init__(self, mol: Chem.Mol):
        self.mol = mol

    def get_atoms(self,
                 generations: Optional[Iterable[int]] = None,
                 branches: Optional[Iterable[int]] = None,
                 steps: Optional[Iterable[int]] = None) -> List[int]:
        """
        Get atoms by any combination of criteria.

        Args:
            generations: List of generations to include (e.g., [0, 1]). None = all generations
            branches: List of branches to include (e.g., [0, 1]). None = all branches
            steps: List of steps to include (e.g., [0, 1]). None = all steps

        Returns:
            List of atom indices matching ALL specified criteria
        """
        # Convert to sets for fast lookup
        gen_set = set(generations) if generations is not None else None
        branch_set = set(branches) if branches is not None else None
        step_set = set(steps) if steps is not None else None

        indices = []
        for atom in self.mol.GetAtoms():
            atom_idx = atom.GetIdx()

            # Check generation filter
            if gen_set is not None:
                atom_gen = self._get_atom_prop(atom, "generation")
                if atom_gen == "" or int(atom_gen) not in gen_set:
                    continue

            # Check branch filter
            if branch_set is not None:
                atom_branch = self._get_atom_prop(atom, "branch")
                if atom_branch == "" or int(atom_branch) not in branch_set:
                    continue

            # Check step filter
            if step_set is not None:
                atom_step = self._get_atom_prop(atom, "step")
                if atom_step == "" or int(atom_step) not in step_set:
                    continue

            indices.append(atom_idx)

        return indices

    def get_original_indices(self,
                             generations: Optional[Iterable[int]] = None,
                             branches: Optional[Iterable[int]] = None,
                             steps: Optional[Iterable[int]] = None) -> List[int]:
        """
        Возвращает замороженные исходные индексы атомов вместо текущих RDKit индексов.
        """
        current_indices = self.get_atoms(generations, branches, steps)
        original_indices = []

        for atom_idx in current_indices:
            atom = self.mol.GetAtomWithIdx(atom_idx)
            if atom.HasProp("original_index"):
                original_indices.append(int(atom.GetProp("original_index")))
            else:
                original_indices.append(atom_idx)  # fallback

        return original_indices

    def get_atom_metadata(self, atom_idx: int) -> Dict[str, str]:
        """Get all metadata for a specific atom."""
        if atom_idx >= self.mol.GetNumAtoms():
            return {}

        atom = self.mol.GetAtomWithIdx(atom_idx)
        return {key: atom.GetProp(key) for key in atom.GetPropNames()}

    def _get_atom_prop(self, atom, prop_name: str, default: str = "") -> str:
        """Safely get atom property."""
        if atom.HasProp(prop_name):
            return atom.GetProp(prop_name)
        return default

    def print_all_metadata(self):
        """Print information about all atoms with their metadata."""
        print(f"\nTotal atoms: {self.mol.GetNumAtoms()}")
        for atom in self.mol.GetAtoms():
            metadata = self.get_atom_metadata(atom.GetIdx())
            formatted_metadata = " | ".join(f"{k}:{v}" for k, v in sorted(metadata.items()))
            print(f"Atom {atom.GetIdx():2d}: {formatted_metadata}")
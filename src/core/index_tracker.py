"""
Tracking and management of atom indices during dendrimer construction.
"""

from typing import Dict, List

from rdkit import Chem


class IndexTracker:
    """Tracks atom indices during molecular assembly."""

    def __init__(self, mol: Chem.Mol):
        self.mol = mol

    def get_core_indices(self) -> List[int]:
        """Get indices of core atoms (generation 0)."""
        core_indices = []
        for atom in self.mol.GetAtoms():
            if self._get_atom_prop(atom, "generation") == "0":
                core_indices.append(atom.GetIdx())
        return core_indices

    def get_generation_indices(self, generation: int) -> List[int]:
        """Get indices of atoms from specific generation."""
        gen_indices = []
        for atom in self.mol.GetAtoms():
            if self._get_atom_prop(atom, "generation") == str(generation):
                gen_indices.append(atom.GetIdx())
        return gen_indices

    def get_branch_indices(self, generation: int, branch: int) -> List[int]:
        """Get indices of atoms from specific branch."""
        branch_indices = []
        for atom in self.mol.GetAtoms():
            if (self._get_atom_prop(atom, "generation") == str(generation) and
                    self._get_atom_prop(atom, "branch") == str(branch)):
                branch_indices.append(atom.GetIdx())
        return branch_indices

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
            print(f"Atom {atom.GetIdx()}: {metadata}")

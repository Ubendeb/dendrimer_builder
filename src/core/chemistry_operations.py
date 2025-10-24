"""
Basic chemistry operations for molecular manipulation.
"""
from rdkit import Chem


class ChemistryOperations:
    """Provides basic chemistry operations."""

    @staticmethod
    def merge_molecules( mol1, mol2, atom_idx1, atom_idx2, bond_order=1):
        """
        Merge two molecules into one.

        Args:
            mol1: First molecule
            mol2: Second molecule

        Returns:
            Merged molecule
        """
        if mol1 is None or mol2 is None:
            return None

            # Создаем редактируемую молекулу из первой
        combined = Chem.RWMol(mol1)

        # Добавляем вторую молекулу
        mol2_frag = Chem.Mol(mol2)
        combined.InsertMol(mol2_frag)

        return combined.GetMol()

    @staticmethod
    def remove_atoms(molecule, atom_indices):
        """
        Remove specified atoms from molecule.

        Args:
            molecule: Input molecule
            atom_indices: List of atom indices to remove

        Returns:
            Modified molecule
        """
        raise NotImplementedError("Implementation for atom removal")

    @staticmethod
    def add_bond(molecule, atom1_idx, atom2_idx, bond_order=1):
        """
        Add bond between two atoms.

        Args:
            molecule: Input molecule
            atom1_idx: Index of first atom
            atom2_idx: Index of second atom
            bond_order: Bond order (default: 1)

        Returns:
            Modified molecule with new bond
        """
        raise NotImplementedError("Implementation for adding bonds")

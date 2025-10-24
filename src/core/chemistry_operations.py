"""
Basic chemistry operations for molecular manipulation.
"""
from rdkit import Chem


class ChemistryOperations:
    """Provides basic chemistry operations."""

    def merge_molecules(self, mol1, mol2):
        """
        Объединяет две молекулы, соединяя их через точки connection_type='to_branch' и 'to_core'
        с сохранением метаданных.
        """
        # Находим точки соединения
        mol1_branch_points = self._get_connection_points(mol1, 'to_branch')
        mol2_core_points = self._get_connection_points(mol2, 'to_core')

        if not mol1_branch_points or not mol2_core_points:
            raise ValueError("Не найдены подходящие точки соединения")

        # Берем первую доступную пару для соединения
        mol1_branch_atom = mol1_branch_points[0]
        mol2_core_atom = mol2_core_points[0]

        # Создаем копии молекул для модификации
        mol1_copy = Chem.RWMol(mol1)
        mol2_copy = Chem.RWMol(mol2)

        # Объединяем молекулы
        combined = Chem.CombineMols(mol1_copy, mol2_copy)
        combined_rw = Chem.RWMol(combined)

        # Создаем связь между атомами
        mol1_idx = mol1_branch_atom['atom_index']
        mol2_idx = mol2_core_atom['atom_index'] + mol1.GetNumAtoms()

        combined_rw.AddBond(mol1_idx, mol2_idx, Chem.BondType.SINGLE)

        # Обновляем метаданные соединенных атомов
        self._update_connection_metadata(combined_rw, mol1_idx, mol2_idx)

        return combined_rw

    def _get_connection_points(self, mol, connection_type):
        """Находит атомы с указанным типом соединения"""
        points = []
        for atom in mol.GetAtoms():
            if (atom.HasProp("is_connection") and
                    atom.GetProp("is_connection") == "true" and
                    atom.HasProp("connection_type") and
                    atom.GetProp("connection_type") == connection_type):
                points.append({
                    'atom_index': atom.GetIdx(),
                    'replacement_group': atom.GetProp("replacement_group")
                })
        return points

    def _update_connection_metadata(self, mol, atom1_idx, atom2_idx):
        """Обновляет метаданные соединенных атомов"""
        atom1 = mol.GetAtomWithIdx(atom1_idx)
        atom2 = mol.GetAtomWithIdx(atom2_idx)

        # Помечаем атомы как соединенные
        atom1.SetProp("connected_to", str(atom2_idx))
        atom2.SetProp("connected_to", str(atom1_idx))

        # Убираем флаги соединения, так как они уже использованы
        atom1.SetProp("is_connection", "false")
        atom2.SetProp("is_connection", "false")

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

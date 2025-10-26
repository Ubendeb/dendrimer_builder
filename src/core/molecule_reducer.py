from rdkit import Chem
from src.core.replacement_analyzer import ReplacementAnalyzer


class MoleculeReducer:
    """Выполняет операции уменьшения молекулы."""

    def reduce(self, mol, connection_atom_idx):
        """
        Уменьшает молекулу согласно replacement group спецификации.
        Возвращает новую молекулу с сохраненными метаданными.
        """
        analyzer = ReplacementAnalyzer()
        analysis_result = analyzer.analyze_replacement_structure(mol, connection_atom_idx)
        atoms_to_remove = analysis_result['atoms']
        bond_to_break = analysis_result['bond_to_break']

        print(f"\n=== ВЫПОЛНЕНИЕ REDUCE ===")
        print(f"Атомы для удаления: {atoms_to_remove}")
        print(f"Связь для разрыва: {bond_to_break}")

        # Создаем редактируемую молекулу
        editable_mol = Chem.RWMol(mol)

        if atoms_to_remove:
            # Случай 1: Удаляем указанные атомы (хвост или ветка)
            self._remove_atoms_with_metadata(editable_mol, atoms_to_remove)
            print(f"Удалены атомы: {atoms_to_remove}")

        elif bond_to_break is not None:
            # Случай 2: Разрываем указанную связь (уменьшаем кратность)
            self._reduce_bond_order(editable_mol, bond_to_break)
            print(f"Уменьшена кратность связи: {bond_to_break}")

        else:
            # Случай 3: Удаляем водород (просто удаляем connection atom)
            self._remove_hydrogen_replacement(editable_mol, connection_atom_idx)

        reduced_mol = editable_mol.GetMol()

        print("Reduce завершен успешно")
        return reduced_mol

    def _remove_atoms_with_metadata(self, editable_mol, atom_indices):
        """Удаляет атомы с сохранением метаданных оставшихся атомов."""
        # Сортируем индексы в обратном порядке чтобы избежать проблем со сдвигом индексов
        for atom_idx in sorted(atom_indices, reverse=True):
            if atom_idx < editable_mol.GetNumAtoms():
                editable_mol.RemoveAtom(atom_idx)

    def _reduce_bond_order(self, editable_mol, bond_idx):
        """Уменьшает кратность связи."""
        if bond_idx < editable_mol.GetNumBonds():
            bond = editable_mol.GetBondWithIdx(bond_idx)
            current_order = bond.GetBondType()

            # Уменьшаем кратность связи
            if current_order == Chem.BondType.TRIPLE:
                new_order = Chem.BondType.DOUBLE
            elif current_order == Chem.BondType.DOUBLE:
                new_order = Chem.BondType.SINGLE
            elif current_order == Chem.BondType.SINGLE:
                new_order = Chem.BondType.SINGLE  # Оставляем одинарной
            else:
                new_order = Chem.BondType.SINGLE

            # Получаем индексы атомов связи
            begin_atom = bond.GetBeginAtomIdx()
            end_atom = bond.GetEndAtomIdx()

            # Удаляем старую связь и создаем новую с уменьшенной кратностью
            editable_mol.RemoveBond(begin_atom, end_atom)
            editable_mol.AddBond(begin_atom, end_atom, new_order)

    def _remove_hydrogen_replacement(self, editable_mol, connection_atom_idx):
        """Добавляет явные водороды только к атому соединения."""
        if connection_atom_idx >= editable_mol.GetNumAtoms():
            return

        connection_atom = editable_mol.GetAtomWithIdx(connection_atom_idx)

        # Получаем текущее количество водородов у целевого атома
        current_h_count = connection_atom.GetTotalNumHs()

        if current_h_count == 0:
            print(f"Предупреждение: у атома {connection_atom_idx} нет водородов")
            return
        current_h_count = current_h_count-1
        connection_atom.SetNumExplicitHs(current_h_count)

        connection_atom.SetNoImplicit(True)

        editable_mol.UpdatePropertyCache()

        try:
            Chem.SanitizeMol(editable_mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_ADJUSTHS)
        except:
            # Если санитизация не проходит, просто обновляем свойства
            editable_mol.UpdatePropertyCache()

        print(f"Добавлено {current_h_count} явных водородов к атому {connection_atom_idx}")

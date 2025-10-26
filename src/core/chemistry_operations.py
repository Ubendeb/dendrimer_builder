"""
Basic chemistry operations for molecular manipulation.
"""
from rdkit import Chem

from src.core.molecule_reducer import MoleculeReducer
from src.core.replacement_analyzer import ReplacementAnalyzer


class ChemistryOperations:
    """Provides basic chemistry operations."""

    def __init__(self):
        self.analyzer = ReplacementAnalyzer()
        self.reducer = MoleculeReducer()

    def analyze_replacement_structure(self, mol, connection_atom_idx):
        """
        Анализирует структуру replacement group и находит соответствующие атомы в молекуле.
        Возвращает словарь с ключами:
        - 'atoms': список атомов для удаления/замещения
        - 'bond_to_break': индекс связи для разрыва (если есть)
        """
        return self.analyzer.analyze_replacement_structure(mol, connection_atom_idx)

    def reduce(self, mol, connection_atom_idx):
        """
        Уменьшает молекулу согласно replacement group спецификации.
        Возвращает новую молекулу с сохраненными метаданными.
        """
        return self.reducer.reduce(mol, connection_atom_idx)

    def connect(self, base_mol, base_mol_connection_atom_idx, additional_mol, additional_mol_connection_atom_idx,
                generation=None, branch=None, step=None):
        """
        Соединяет две молекулы через указанные атомы соединения.
        Возвращает новую молекулу с сохраненными метаданными.

        Args:
            base_mol: базовая молекула
            base_mol_connection_atom_idx: индекс атома соединения в базовой молекуле
            additional_mol: добавляемая молекула
            additional_mol_connection_atom_idx: индекс атома соединения в добавляемой молекуле
            generation: значение для свойства generation атомов добавляемой молекулы
            branch: значение для свойства branch атомов добавляемой молекулы
            step: значение для свойства step атомов добавляемой молекулы
        """
        # Уменьшаем обе молекулы согласно их replacement group спецификациям
        reduced_base = self.reducer.reduce(base_mol, base_mol_connection_atom_idx)
        reduced_additional = self.reducer.reduce(additional_mol, additional_mol_connection_atom_idx)

        # Устанавливаем свойства для атомов добавляемой молекулы, если параметры заданы
        if any(param is not None for param in [generation, branch, step]):
            self.set_atom_properties(reduced_additional, generation, branch, step)

        # Соединяем уменьшенные молекулы
        connected_mol = self._connect_reduced(reduced_base, base_mol_connection_atom_idx,
                                              reduced_additional, additional_mol_connection_atom_idx)

        return connected_mol

    def set_atom_properties(self, mol, generation, branch, step):
        """
        Устанавливает свойства для всех атомов молекулы.
        """
        for atom in mol.GetAtoms():
            if generation is not None:
                atom.SetProp("generation", str(generation))
            if branch is not None:
                atom.SetProp("branch", str(branch))
            if step is not None:
                atom.SetProp("step", str(step))


    def _connect_reduced(self, base_mol, base_conn_idx, additional_mol, additional_conn_idx):
        """
        Соединяет две уже уменьшенные молекулы через указанные атомы соединения.
        """
        print(f"\n=== ВЫПОЛНЕНИЕ CONNECT ===")
        print(f"Базовая молекула: атом {base_conn_idx}")
        print(f"Дополнительная молекула: атом {additional_conn_idx}")

        # Создаем комбинированную молекулу
        combined_mol = Chem.CombineMols(base_mol, additional_mol)
        editable_mol = Chem.RWMol(combined_mol)

        # Вычисляем новые индексы атомов после объединения
        base_atoms_count = base_mol.GetNumAtoms()
        new_base_conn_idx = base_conn_idx
        new_additional_conn_idx = additional_conn_idx + base_atoms_count

        print(f"После объединения: базовый атом {new_base_conn_idx}, дополнительный атом {new_additional_conn_idx}")

        # Добавляем связь между атомами соединения
        editable_mol.AddBond(new_base_conn_idx, new_additional_conn_idx, Chem.BondType.SINGLE)

        # Обновляем свойства молекулы
        connected_mol = editable_mol.GetMol()

        try:
            # Санитизируем молекулу для корректного расчета свойств
            Chem.SanitizeMol(connected_mol)
        except:
            print("Предупреждение: санитизация не прошла, но молекула создана")

        print("Connect завершен успешно")
        return connected_mol

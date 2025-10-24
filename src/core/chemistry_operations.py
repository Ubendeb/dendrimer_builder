"""
Basic chemistry operations for molecular manipulation.
"""
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from collections import deque


class ChemistryOperations:
    """Provides basic chemistry operations."""

    def visualize_molecular_subgraph(self, mol, connection_atom_idx, replacement_group_smiles):
        """
        Визуализирует подграф молекулы вокруг точки соединения.
        """
        # Определяем максимальную глубину обхода
        replacement_mol = Chem.MolFromSmiles(replacement_group_smiles)
        max_depth = replacement_mol.GetNumAtoms() if replacement_mol else 0

        print(f"Визуализация подграфа молекулы")
        print(f"Атом соединения: {mol.GetAtomWithIdx(connection_atom_idx).GetSymbol()}({connection_atom_idx})")
        print(f"Replacement group: {replacement_group_smiles}")
        print(f"Максимальная глубина обхода: {max_depth}")

        # Собираем все атомы для визуализации
        atoms_to_highlight = set()
        bonds_to_highlight = set()

        # Обходим граф в ширину от connection atom
        visited = set([connection_atom_idx])
        queue = deque([(connection_atom_idx, 0)])

        while queue:
            current_idx, depth = queue.popleft()
            current_atom = mol.GetAtomWithIdx(current_idx)

            atoms_to_highlight.add(current_idx)

            if depth < max_depth + 1:
                for neighbor in current_atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()

                    # Добавляем связь
                    bond = mol.GetBondBetweenAtoms(current_idx, neighbor_idx)
                    if bond:
                        bonds_to_highlight.add(bond.GetIdx())

                    if neighbor_idx not in visited:
                        visited.add(neighbor_idx)
                        queue.append((neighbor_idx, depth + 1))

        # Создаем визуализацию
        self.subplot_graph(atoms_to_highlight,   connection_atom_idx, mol, replacement_group_smiles)

        # Выводим информацию
        print(f"\nНайдено атомов в подграфе: {len(atoms_to_highlight)}")
        print("Атомы:", sorted(atoms_to_highlight))

        return atoms_to_highlight

    def subplot_graph(self, atoms_to_highlight,  connection_atom_idx, mol, replacement_group_smiles):
        # Создаем визуализацию дерева
        fig, ax = plt.subplots(figsize=(12, 8))

        # Строим дерево графа
        tree_layout = {}
        levels = {}

        # Размещаем атомы по уровням
        def build_tree_layout(current_idx, depth=0, pos=0):
            if current_idx in tree_layout:
                return pos

            atom = mol.GetAtomWithIdx(current_idx)
            symbol = atom.GetSymbol()

            # Сохраняем позицию
            tree_layout[current_idx] = (pos, -depth)
            if depth not in levels:
                levels[depth] = []
            levels[depth].append(current_idx)

            # Рекурсивно размещаем соседей (исключая уже размещенных)
            new_pos = pos
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in tree_layout and neighbor_idx in atoms_to_highlight:
                    new_pos = build_tree_layout(neighbor_idx, depth + 1, new_pos + 1)

            return new_pos

        # Начинаем построение с connection atom
        build_tree_layout(connection_atom_idx)

        # Рисуем дерево
        for atom_idx, (x, y) in tree_layout.items():
            atom = mol.GetAtomWithIdx(atom_idx)
            symbol = atom.GetSymbol()

            # Цвета: connection atom - синий, остальные - красный
            color = 'lightblue' if atom_idx == connection_atom_idx else 'lightcoral'

            # Рисуем узел
            circle = plt.Circle((x, y), 0.3, fill=True, color=color, ec='black', lw=2)
            ax.add_patch(circle)

            # Подписываем атом
            ax.text(x, y, f'{symbol}({atom_idx})', ha='center', va='center', fontweight='bold')

            # Рисуем связи к соседям
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in tree_layout:
                    nx, ny = tree_layout[neighbor_idx]
                    # Рисуем линию только если сосед уже размещен и это связь вперед
                    if ny < y:  # сосед выше по дереву (ближе к корню)
                        continue
                    ax.plot([x, nx], [y, ny], 'k-', lw=2)

        # Настраиваем отображение
        ax.set_xlim(min(x for x, y in tree_layout.values()) - 1,
                    max(x for x, y in tree_layout.values()) + 1)
        ax.set_ylim(min(y for x, y in tree_layout.values()) - 1,
                    max(y for x, y in tree_layout.values()) + 1)
        ax.set_aspect('equal')
        ax.set_title(f'Дерево графа вокруг атома {connection_atom_idx}\n'
                     f'Replacement: {replacement_group_smiles} | '
                     f'Атомы: {sorted(atoms_to_highlight)}',
                     fontsize=14, pad=20)
        ax.axis('off')

    def analyze_replacement_structure(self, mol, connection_atom_idx, replacement_group_smiles):
        """
        Анализирует структуру replacement group и находит соответствующие атомы в молекуле.
        """
        print(f"\n{'='*60}")
        print(f"ДЕТАЛЬНЫЙ АНАЛИЗ СТРУКТУРЫ")
        print(f"{'='*60}")

        connection_atom = mol.GetAtomWithIdx(connection_atom_idx)
        replacement_mol = Chem.MolFromSmiles(replacement_group_smiles)

        print(f"Атом соединения: {connection_atom.GetSymbol()}({connection_atom_idx})")
        print(f"Replacement group: {replacement_group_smiles}")

        # Анализируем структуру replacement group
        print(f"\nСТРУКТУРА REPLACEMENT GROUP:")
        for i, atom in enumerate(replacement_mol.GetAtoms()):
            neighbors = [f"{n.GetSymbol()}({n.GetIdx()})" for n in atom.GetNeighbors()]
            print(f"   Атом {i}: {atom.GetSymbol()} → соседи: {neighbors}")

        # Анализируем соседей connection atom
        print(f"\nСОСЕДИ АТОМА СОЕДИНЕНИЯ:")
        neighbors = list(connection_atom.GetNeighbors())
        for i, neighbor in enumerate(neighbors):
            neighbor_neighbors = [f"{n.GetSymbol()}({n.GetIdx()})" for n in neighbor.GetNeighbors()
                                if n.GetIdx() != connection_atom_idx]
            print(f"   Сосед {i}: {neighbor.GetSymbol()}({neighbor.GetIdx()}) → соседи: {neighbor_neighbors}")

        # Находим правильную ветку на основе структуры replacement group
        print(f"\nПОИСК СООТВЕТСТВУЮЩЕЙ ВЕТКИ:")
        first_replacement_atom = replacement_mol.GetAtomWithIdx(0)
        expected_symbol = first_replacement_atom.GetSymbol()
        expected_degree = first_replacement_atom.GetDegree()

        print(f"   Первый атом replacement: {expected_symbol}, степень: {expected_degree}")

        candidate_atoms = []

        for neighbor in neighbors:
            if neighbor.GetSymbol() == expected_symbol:
                print(f"\n   Найден кандидат: {neighbor.GetSymbol()}({neighbor.GetIdx()})")

                # Проверяем структуру ветки
                branch_atoms = self._validate_branch_structure(mol, neighbor, connection_atom_idx, replacement_mol)
                if branch_atoms:
                    candidate_atoms = branch_atoms
                    print(f"   Правильная ветка найдена: {branch_atoms}")
                    break
                else:
                    print(f"   Структура не соответствует replacement group")

        return candidate_atoms

    def _validate_branch_structure(self, mol, start_atom, connection_atom_idx, replacement_mol):
        """
        Проверяет, соответствует ли структура ветки replacement group.
        """
        visited = set([connection_atom_idx])
        branch_atoms = []
        replacement_atoms = list(replacement_mol.GetAtoms())

        def validate(atom, replacement_idx):
            atom_idx = atom.GetIdx()

            if atom_idx in visited:
                return True

            if replacement_idx >= len(replacement_atoms):
                return False

            # Проверяем соответствие символа
            expected_atom = replacement_atoms[replacement_idx]
            if atom.GetSymbol() != expected_atom.GetSymbol():
                return False

            visited.add(atom_idx)
            branch_atoms.append(atom_idx)

            # Получаем непосещенных соседей (исключая connection atom)
            unvisited_neighbors = []
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in visited:
                    unvisited_neighbors.append(neighbor)

            # Проверяем количество связей
            if len(unvisited_neighbors) != expected_atom.GetDegree() - (1 if replacement_idx > 0 else 0):
                return False

            # Рекурсивно проверяем соседей
            for i, neighbor in enumerate(unvisited_neighbors):
                if not validate(neighbor, replacement_idx + 1):
                    return False

            return True

        # Запускаем проверку
        if validate(start_atom, 0):
            return branch_atoms
        else:
            return []
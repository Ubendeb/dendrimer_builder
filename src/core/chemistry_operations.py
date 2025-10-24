"""
Basic chemistry operations for molecular manipulation.
"""
from collections import deque

import matplotlib.pyplot as plt
from rdkit import Chem


class ChemistryOperations:
    """Provides basic chemistry operations."""

    def analyze_replacement_structure(self, mol, connection_atom_idx):
        """
        Анализирует структуру replacement group и находит соответствующие атомы в молекуле.
        """
        self._print_analysis_header()
        connection_atom = mol.GetAtomWithIdx(connection_atom_idx)

        replacement_group_smiles = connection_atom.GetProp('replacement_group')
        replacement_mol = Chem.MolFromSmiles(replacement_group_smiles)

        self._print_basic_info(connection_atom, connection_atom_idx, replacement_group_smiles)
        self._analyze_replacement_structure(replacement_mol)
        self._analyze_connection_neighbors(mol, connection_atom, connection_atom_idx)

        candidate_atoms = self._find_matching_branch(mol, connection_atom, connection_atom_idx, replacement_mol)
        return candidate_atoms

    def _collect_subgraph_atoms_and_bonds(self, mol, connection_atom_idx, max_depth):
        """Собирает атомы и связи для визуализации подграфа."""
        atoms_to_highlight = set()
        bonds_to_highlight = set()
        visited = set([connection_atom_idx])
        queue = deque([(connection_atom_idx, 0)])

        while queue:
            current_idx, depth = queue.popleft()
            current_atom = mol.GetAtomWithIdx(current_idx)

            atoms_to_highlight.add(current_idx)

            if depth < max_depth + 1:
                self._process_neighbors(mol, current_idx, current_atom, visited, bonds_to_highlight, queue, depth)

        return atoms_to_highlight, bonds_to_highlight

    def _process_neighbors(self, mol, current_idx, current_atom, visited, bonds_to_highlight, queue, depth):
        """Обрабатывает соседей текущего атома."""
        for neighbor in current_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            self._add_bond_to_highlight(mol, current_idx, neighbor_idx, bonds_to_highlight)
            self._add_unvisited_neighbor_to_queue(neighbor_idx, visited, queue, depth)

    def _add_bond_to_highlight(self, mol, current_idx, neighbor_idx, bonds_to_highlight):
        """Добавляет связь в множество для подсветки."""
        bond = mol.GetBondBetweenAtoms(current_idx, neighbor_idx)
        if bond:
            bonds_to_highlight.add(bond.GetIdx())

    def _add_unvisited_neighbor_to_queue(self, neighbor_idx, visited, queue, depth):
        """Добавляет непосещенного соседа в очередь."""
        if neighbor_idx not in visited:
            visited.add(neighbor_idx)
            queue.append((neighbor_idx, depth + 1))

    def _create_subgraph_visualization(self, atoms_to_highlight, connection_atom_idx, mol, replacement_group_smiles):
        """Создает визуализацию подграфа."""
        fig, ax = plt.subplots(figsize=(12, 8))
        tree_layout, levels = self._build_tree_layout(mol, atoms_to_highlight, connection_atom_idx)
        self._draw_tree(ax, mol, tree_layout, connection_atom_idx)
        self._configure_tree_plot(ax, tree_layout, atoms_to_highlight, connection_atom_idx, replacement_group_smiles)

    def _build_tree_layout(self, mol, atoms_to_highlight, connection_atom_idx):
        """Строит древовидную структуру для визуализации."""
        tree_layout = {}
        levels = {}

        def build_tree(current_idx, depth=0, pos=0):
            if current_idx in tree_layout:
                return pos

            tree_layout[current_idx] = (pos, -depth)
            if depth not in levels:
                levels[depth] = []
            levels[depth].append(current_idx)

            new_pos = pos
            for neighbor in mol.GetAtomWithIdx(current_idx).GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in tree_layout and neighbor_idx in atoms_to_highlight:
                    new_pos = build_tree(neighbor_idx, depth + 1, new_pos + 1)

            return new_pos

        build_tree(connection_atom_idx)
        return tree_layout, levels

    def _draw_tree(self, ax, mol, tree_layout, connection_atom_idx):
        """Рисует дерево графа."""
        for atom_idx, (x, y) in tree_layout.items():
            atom = mol.GetAtomWithIdx(atom_idx)
            symbol = atom.GetSymbol()
            color = 'lightblue' if atom_idx == connection_atom_idx else 'lightcoral'

            circle = plt.Circle((x, y), 0.3, fill=True, color=color, ec='black', lw=2)
            ax.add_patch(circle)
            ax.text(x, y, f'{symbol}({atom_idx})', ha='center', va='center', fontweight='bold')

            self._draw_connections(ax, mol, atom, tree_layout, x, y)

    def _draw_connections(self, ax, mol, atom, tree_layout, x, y):
        """Рисует связи между атомами."""
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in tree_layout:
                nx, ny = tree_layout[neighbor_idx]
                if ny < y:  # сосед выше по дереву (ближе к корню)
                    continue
                ax.plot([x, nx], [y, ny], 'k-', lw=2)

    def _configure_tree_plot(self, ax, tree_layout, atoms_to_highlight, connection_atom_idx, replacement_group_smiles):
        """Настраивает внешний вид графика."""
        x_coords = [x for x, y in tree_layout.values()]
        y_coords = [y for x, y in tree_layout.values()]

        ax.set_xlim(min(x_coords) - 1, max(x_coords) + 1)
        ax.set_ylim(min(y_coords) - 1, max(y_coords) + 1)
        ax.set_aspect('equal')
        ax.set_title(
            f'Дерево графа вокруг атома {connection_atom_idx}\n'
            f'Replacement: {replacement_group_smiles} | '
            f'Атомы: {sorted(atoms_to_highlight)}',
            fontsize=14, pad=20
        )
        ax.axis('off')

    def _print_subgraph_info(self, atoms_to_highlight):
        """Выводит информацию о подграфе."""
        print(f"\nНайдено атомов в подграфе: {len(atoms_to_highlight)}")
        print("Атомы:", sorted(atoms_to_highlight))

    def _print_analysis_header(self):
        """Выводит заголовок анализа."""
        print(f"\n{'=' * 60}")
        print(f"Детальный анализ структуры")
        print(f"{'=' * 60}")

    def _print_basic_info(self, connection_atom, connection_atom_idx, replacement_group_smiles):
        """Выводит базовую информацию."""
        print(f"Атом соединения: {connection_atom.GetSymbol()}({connection_atom_idx})")
        print(f"Replacement group: {replacement_group_smiles}")

    def _analyze_replacement_structure(self, replacement_mol):
        """Анализирует структуру replacement group."""
        print(f"\nСтруктура replacement_group:")
        for i, atom in enumerate(replacement_mol.GetAtoms()):
            neighbors = [f"{n.GetSymbol()}({n.GetIdx()})" for n in atom.GetNeighbors()]
            print(f"   Атом {i}: {atom.GetSymbol()} → соседи: {neighbors}")

    def _analyze_connection_neighbors(self, mol, connection_atom, connection_atom_idx):
        """Анализирует соседей атома соединения."""
        print(f"\nСоседи атома соединения:")
        neighbors = list(connection_atom.GetNeighbors())
        for i, neighbor in enumerate(neighbors):
            neighbor_neighbors = [
                f"{n.GetSymbol()}({n.GetIdx()})" for n in neighbor.GetNeighbors()
                if n.GetIdx() != connection_atom_idx
            ]
            print(f"   Сосед {i}: {neighbor.GetSymbol()}({neighbor.GetIdx()}) → соседи: {neighbor_neighbors}")

    def _find_matching_branch(self, mol, connection_atom, connection_atom_idx, replacement_mol):
        """Находит соответствующую ветку в молекуле."""
        print(f"\nПоиск соответствующей ветки:")
        first_replacement_atom = replacement_mol.GetAtomWithIdx(0)
        expected_symbol = first_replacement_atom.GetSymbol()
        expected_degree = first_replacement_atom.GetDegree()

        print(f"   Первый атом replacement: {expected_symbol}, степень: {expected_degree}")

        candidate_atoms = []
        for neighbor in connection_atom.GetNeighbors():
            if neighbor.GetSymbol() == expected_symbol:
                print(f"\n   Найден кандидат: {neighbor.GetSymbol()}({neighbor.GetIdx()})")
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

            expected_atom = replacement_atoms[replacement_idx]
            if atom.GetSymbol() != expected_atom.GetSymbol():
                return False

            visited.add(atom_idx)
            branch_atoms.append(atom_idx)

            unvisited_neighbors = self._get_unvisited_neighbors(atom, visited)
            expected_degree = expected_atom.GetDegree() - (1 if replacement_idx > 0 else 0)

            if len(unvisited_neighbors) != expected_degree:
                return False

            for neighbor in unvisited_neighbors:
                if not validate(neighbor, replacement_idx + 1):
                    return False

            return True

        return branch_atoms if validate(start_atom, 0) else []

    def _get_unvisited_neighbors(self, atom, visited):
        """Возвращает непосещенных соседей атома."""
        return [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetIdx() not in visited]

"""
Molecular visualization and rendering.
"""
from collections import deque
from matplotlib import pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

IPythonConsole.ipython_useSVG = True


class MoleculeRenderer:
    """Handles molecular visualization."""

    def __init__(self):
        self.render_styles = {}

    def render_2d(self, mol, title="", highlight_atoms=None, highlight_bonds =None, display_props=None):
        """
        Render 2D molecular structure.

        Args:
            mol: Molecule to render
            highlight_atoms: Atoms to highlight
            highlight_bonds: Bonds to highlight
            title: Name of molecule
            display_props: List of property names to display on atoms

        Returns:
            2D visualization
        """
        if mol is None:
            print(f"{title}: None")
            return None

        if title == "":
            title = f"{Chem.MolToSmiles(mol)}"

        mol_copy = Chem.Mol(mol)
        self._set_atom_display_notes(mol_copy, display_props)

        coef = mol.GetNumAtoms()

        img = Draw.MolToImage(mol_copy, legend=title, size=(max(60 * coef, 300), max(40 * coef, 200)),
                              highlightAtoms=highlight_atoms if highlight_atoms else [],
                              highlightBonds=highlight_bonds if highlight_bonds else [])
        return img

    def _set_atom_display_notes(self, mol, display_props=None):
        """
        Set atom display notes based on specified properties.
        """
        for atom in mol.GetAtoms():
            numbers = []
            booleans = []
            strings = []

            # Always include atom index
            numbers.append(str(atom.GetIdx()))

            if display_props:
                for prop in display_props:
                    if not atom.HasProp(prop):
                        continue

                    value = atom.GetProp(prop)

                    if value.lstrip('-').replace('.', '', 1).isdigit():
                        numbers.append(value)
                    elif value.lower() in ['true', 'false']:
                        if value.lower() == 'true':
                            booleans.append(prop)
                    else:
                        strings.append(value)


            display_parts = []
            if numbers:
                display_parts.append(":".join(numbers))
            if booleans:
                display_parts.append(":".join(booleans))
            if strings:
                display_parts.append(":".join(strings))


            note_text = " | ".join(display_parts)
            atom.SetProp('atomNote', note_text)

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


    def visualize_molecular_subgraph(self, mol, connection_atom_idx):
        """
        Визуализирует подграф молекулы вокруг точки соединения.
        """
        connection_atom = mol.GetAtomWithIdx(connection_atom_idx)
        replacement_group_smiles = connection_atom.GetProp('replacement_group')
        replacement_mol = Chem.MolFromSmiles(replacement_group_smiles)
        max_depth = replacement_mol.GetNumAtoms() if replacement_mol else 0

        self._print_visualization_info(mol, connection_atom_idx, replacement_group_smiles, max_depth)

        atoms_to_highlight, bonds_to_highlight = self._collect_subgraph_atoms_and_bonds(
            mol, connection_atom_idx, max_depth
        )

        self._create_subgraph_visualization(atoms_to_highlight, connection_atom_idx, mol, replacement_group_smiles)
        self._print_subgraph_info(atoms_to_highlight)

        return atoms_to_highlight

    def _print_visualization_info(self, mol, connection_atom_idx, replacement_group_smiles, max_depth):
        """Выводит информацию о визуализации."""
        print(f"Визуализация подграфа молекулы")
        print(f"Атом соединения: {mol.GetAtomWithIdx(connection_atom_idx).GetSymbol()}({connection_atom_idx})")
        print(f"Replacement group: {replacement_group_smiles}")
        print(f"Максимальная глубина обхода: {max_depth}")

    def visualize_whole_molecule(self, mol, title="Molecular Graph"):
        """
        Визуализирует всю молекулу в виде графа.

        Args:
            mol: Молекула для визуализации
            title: Заголовок графика

        Returns:
            Множество индексов всех атомов в молекуле
        """
        if mol is None:
            print("Молекула не задана")
            return set()

        print(f"Визуализация всей молекулы: {title}")
        print(f"Общее количество атомов: {mol.GetNumAtoms()}")
        print(f"Общее количество связей: {mol.GetNumBonds()}")

        # Собираем все атомы и связи
        all_atoms = set(range(mol.GetNumAtoms()))
        all_bonds = set(range(mol.GetNumBonds()))

        # Создаем визуализацию
        fig, ax = plt.subplots(figsize=(14, 10))

        # Строим граф используя первый атом как корень
        if all_atoms:
            root_atom_idx = 0
            graph_layout, levels = self._build_molecule_graph_layout(mol, all_atoms, root_atom_idx)
            self._draw_molecule_graph(ax, mol, graph_layout, all_bonds)
            self._configure_molecule_plot(ax, graph_layout, title, all_atoms)

        self._print_molecule_info(all_atoms, all_bonds)

        return all_atoms

    def _build_molecule_graph_layout(self, mol, all_atoms, root_atom_idx):
        """
        Строит layout для всей молекулы.

        Args:
            mol: Молекула
            all_atoms: Множество всех атомов
            root_atom_idx: Индекс корневого атома

        Returns:
            graph_layout: Словарь с координатами атомов {idx: (x, y)}
            levels: Уровни в дереве
        """
        graph_layout = {}
        levels = {}
        visited = set()

        def bfs_layout(start_idx):
            queue = deque([(start_idx, 0, 0)])
            max_pos_at_level = {}

            while queue:
                current_idx, depth, horizontal_pos = queue.popleft()

                if current_idx in visited:
                    continue

                visited.add(current_idx)

                # Определяем позицию по горизонтали для текущего уровня
                if depth not in max_pos_at_level:
                    max_pos_at_level[depth] = 0
                else:
                    max_pos_at_level[depth] += 1

                x = max_pos_at_level[depth]
                y = -depth

                graph_layout[current_idx] = (x, y)

                if depth not in levels:
                    levels[depth] = []
                levels[depth].append(current_idx)

                # Добавляем соседей в очередь
                current_atom = mol.GetAtomWithIdx(current_idx)
                neighbors = list(current_atom.GetNeighbors())

                # Сортируем соседей для более предсказуемого layout
                neighbors.sort(key=lambda x: x.GetIdx())

                for i, neighbor in enumerate(neighbors):
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited and neighbor_idx in all_atoms:
                        queue.append((neighbor_idx, depth + 1, i))

        # Запускаем BFS с корневого атома
        bfs_layout(root_atom_idx)

        # Обрабатываем несвязные компоненты (если есть)
        disconnected_atoms = all_atoms - visited
        if disconnected_atoms:
            print(f"Обнаружены несвязные компоненты: {len(disconnected_atoms)} атомов")
            start_x = max([x for x, y in graph_layout.values()]) + 2 if graph_layout else 0
            for i, atom_idx in enumerate(disconnected_atoms):
                graph_layout[atom_idx] = (start_x + i, 0)
                if 0 not in levels:
                    levels[0] = []
                levels[0].append(atom_idx)

        return graph_layout, levels

    def _draw_molecule_graph(self, ax, mol, graph_layout, all_bonds):
        """
        Рисует граф молекулы.

        Args:
            ax: Ось matplotlib
            mol: Молекула
            graph_layout: Расположение атомов
            all_bonds: Все связи молекулы
        """
        # Рисуем связи
        for bond_idx in all_bonds:
            bond = mol.GetBondWithIdx(bond_idx)
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()

            if begin_idx in graph_layout and end_idx in graph_layout:
                x1, y1 = graph_layout[begin_idx]
                x2, y2 = graph_layout[end_idx]

                # Определяем стиль линии в зависимости от типа связи
                bond_type = bond.GetBondType()
                if bond_type == Chem.rdchem.BondType.SINGLE:
                    linestyle = '-'
                    linewidth = 2
                elif bond_type == Chem.rdchem.BondType.DOUBLE:
                    linestyle = '--'
                    linewidth = 3
                elif bond_type == Chem.rdchem.BondType.TRIPLE:
                    linestyle = ':'
                    linewidth = 4
                else:
                    linestyle = '-'
                    linewidth = 2

                ax.plot([x1, x2], [y1, y2], linestyle, color='black',
                        linewidth=linewidth, alpha=0.7)

        # Рисуем атомы
        for atom_idx, (x, y) in graph_layout.items():
            atom = mol.GetAtomWithIdx(atom_idx)
            symbol = atom.GetSymbol()

            # Выбираем цвет в зависимости от элемента
            color = self._get_atom_color(atom)

            circle = plt.Circle((x, y), 0.4, fill=True, color=color,
                                ec='black', lw=2, alpha=0.8)
            ax.add_patch(circle)

            # Подписываем атом
            label = f'{symbol}\n({atom_idx})'
            ax.text(x, y, label, ha='center', va='center',
                    fontweight='bold', fontsize=9)

    def _get_atom_color(self, atom):
        """
        Возвращает цвет для атома в зависимости от элемента.

        Args:
            atom: Атом

        Returns:
            Цвет для отображения
        """
        element = atom.GetSymbol()
        color_map = {
            'C': 'lightgray',
            'O': 'red',
            'N': 'blue',
            'H': 'white',
            'S': 'yellow',
            'P': 'orange',
            'F': 'green',
            'Cl': 'lime',
            'Br': 'darkred',
            'I': 'purple'
        }
        return color_map.get(element, 'lightblue')

    def _configure_molecule_plot(self, ax, graph_layout, title, all_atoms):
        """
        Настраивает внешний вид графика молекулы.

        Args:
            ax: Ось matplotlib
            graph_layout: Расположение атомов
            title: Заголовок
            all_atoms: Все атомы молекулы
        """
        if not graph_layout:
            return

        x_coords = [x for x, y in graph_layout.values()]
        y_coords = [y for x, y in graph_layout.values()]

        padding = 1.5
        ax.set_xlim(min(x_coords) - padding, max(x_coords) + padding)
        ax.set_ylim(min(y_coords) - padding, max(y_coords) + padding)
        ax.set_aspect('equal')

        ax.set_title(f'{title}\nВсего атомов: {len(all_atoms)}',
                     fontsize=14, pad=20, fontweight='bold')
        ax.axis('off')

    def _print_molecule_info(self, all_atoms, all_bonds):
        """
        Выводит информацию о молекуле.

        Args:
            all_atoms: Все атомы
            all_bonds: Все связи
        """
        print(f"\nИнформация о молекуле:")
        print(f"Всего атомов: {len(all_atoms)}")
        print(f"Всего связей: {len(all_bonds)}")
        print(f"Атомы: {sorted(all_atoms)}")
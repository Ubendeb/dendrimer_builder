from src.core.replacers.replacer import ReplacementParser


class TailReplacementParser(ReplacementParser):
    """Парсер для хвостовых спецификаций в формате 3_4."""

    def parse(self, spec, mol, connection_atom_idx):
        """Парсит хвостовую спецификацию и находит соответствующие атомы в молекуле."""
        tail_atoms = []
        current_atom_idx = connection_atom_idx

        # Разбираем цепочку original_index'ов
        atom_specs = spec.split('_')

        for i, atom_spec in enumerate(atom_specs):
            target_original_idx = int(atom_spec)
            current_atom = mol.GetAtomWithIdx(current_atom_idx)

            # Ищем соседа с нужным original_index
            found = False
            for neighbor in current_atom.GetNeighbors():
                neighbor_original_idx = self._get_original_index(neighbor)
                if neighbor_original_idx == target_original_idx:
                    tail_atoms.append(neighbor.GetIdx())
                    current_atom_idx = neighbor.GetIdx()
                    found = True
                    break

            if not found:
                print(f"Предупреждение: не найден атом с original_index {target_original_idx}")
                break

        return {
            'type': 'tail',
            'tail_atoms': tail_atoms,
            'mol': None  # Для хвостов не нужна молекула replacement
        }

    def _get_original_index(self, atom):
        """Возвращает original index атома."""
        if atom.HasProp('original_index'):
            return int(atom.GetProp('original_index'))
        return atom.GetIdx()


def _analyze_replacement_structure(self, replacement_mol):
    """Анализирует структуру replacement group."""
    print(f"\nСтруктура replacement_group:")
    if replacement_mol is None:
        print("   Хвостовая спецификация - молекула не требуется")
        return

    for i, atom in enumerate(replacement_mol.GetAtoms()):
        neighbors = [f"{n.GetSymbol()}({n.GetIdx()})" for n in atom.GetNeighbors()]
        print(f"   Атом {i}: {atom.GetSymbol()} → соседи: {neighbors}")
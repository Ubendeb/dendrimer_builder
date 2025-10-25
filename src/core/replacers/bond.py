from src.core.replacers.replacer import ReplacementParser


class BondBreakReplacementParser(ReplacementParser):
    """Парсер для спецификаций с указанием связи для разрыва."""

    def parse(self, spec, mol=None, connection_atom_idx=None):
        from rdkit import Chem
        # Форматы: "=3", "#1", "C=3"
        if '=' in spec:
            parts = spec.split('=')
            smiles_part = parts[0] if parts[0] else 'C'
            bond_type = 'DOUBLE'
            target_idx = int(parts[1])
        elif '#' in spec:
            parts = spec.split('#')
            smiles_part = parts[0] if parts[0] else 'C'
            bond_type = 'TRIPLE'
            target_idx = int(parts[1])
        else:
            smiles_part = 'C'
            bond_type = 'SINGLE'
            target_idx = int(spec)

        return {
            'mol': Chem.MolFromSmiles(smiles_part),
            'type': 'bond_break',
            'break_bond': (bond_type, target_idx)
        }
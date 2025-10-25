from src.core.replacers.replacer import ReplacementParser


class HydrogenReplacementParser(ReplacementParser):
    """Парсер для замены на водород."""

    def parse(self, spec, mol=None, connection_atom_idx=None):
        from rdkit import Chem
        return {
            'mol': None,
            'type': 'hydrogen'
        }

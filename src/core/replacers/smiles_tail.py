from src.core.replacers.replacer import ReplacementParser

class SmilesReplacementParser(ReplacementParser):
    """Парсер для стандартных SMILES строк."""

    def parse(self, spec, mol=None, connection_atom_idx=None):
        from rdkit import Chem
        return {
            'mol': Chem.MolFromSmiles(spec),
            'type': 'smiles'
        }
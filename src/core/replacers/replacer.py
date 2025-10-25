class ReplacementParser:
    """Базовый класс для парсеров replacement specification."""

    def parse(self, spec, mol=None, connection_atom_idx=None):
        raise NotImplementedError
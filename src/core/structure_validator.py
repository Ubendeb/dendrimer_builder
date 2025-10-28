"""
Validation of molecular structures and connections.
"""

from rdkit import Chem


class StructureValidator:
    """Validates molecular structures and connections."""

    def __init__(self):
        self.validation_rules = {}

    def validate_smiles(self, smiles_string):
        """
        Validate SMILES string and return molecule if valid.

        Args:
            smiles_string: SMILES representation

        Returns:
            tuple: (is_valid, molecule, error_message)
        """
        if not smiles_string or not smiles_string.strip():
            return False, None, "Empty SMILES string"

        mol = Chem.MolFromSmiles(smiles_string.strip())
        if mol is None:
            return False, None, "Invalid SMILES format"

        return True, mol, "Valid SMILES"

    def validate_atom_indices(self, molecule, atom_indices):
        """
        Validate if atom indices are within molecule bounds.

        Args:
            molecule: RDKit Mol object
            atom_indices: List of atom indices to validate

        Returns:
            tuple: (is_valid, error_message)
        """
        if not atom_indices:
            return False, "Empty atom indices list"

        num_atoms = molecule.GetNumAtoms()
        for atom_idx in atom_indices:
            if atom_idx < 0 or atom_idx >= num_atoms:
                return False, f"Atom index {atom_idx} out of range (0-{num_atoms-1})"

        return True, "Valid atom indices"

    def validate_connection_points(self, fragment, connection_points):
        """
        Validate connection points for a fragment.

        Args:
            fragment: Molecular fragment data
            connection_points: List of connection atom indices

        Returns:
            tuple: (is_valid, error_message)
        """
        if 'mol' not in fragment:
            return False, "Fragment has no molecule object"

        return self.validate_atom_indices(fragment['mol'], connection_points)

    def validate_connection(self, fragment1, fragment2, connection_points):
        """
        Validate if two fragments can be connected at specified points.

        Args:
            fragment1: First molecular fragment data
            fragment2: Second molecular fragment data
            connection_points: Tuple of (atom1_idx, atom2_idx)

        Returns:
            tuple: (is_valid, error_message)
        """
        if 'mol' not in fragment1 or 'mol' not in fragment2:
            return False, "One or both fragments missing molecule object"

        atom1_idx, atom2_idx = connection_points

        # Validate atom indices
        valid1, msg1 = self.validate_atom_indices(fragment1['mol'], [atom1_idx])
        if not valid1:
            return False, f"Fragment1: {msg1}"

        valid2, msg2 = self.validate_atom_indices(fragment2['mol'], [atom2_idx])
        if not valid2:
            return False, f"Fragment2: {msg2}"

        # Check if connection atoms are suitable for bonding
        atom1 = fragment1['mol'].GetAtomWithIdx(atom1_idx)
        atom2 = fragment2['mol'].GetAtomWithIdx(atom2_idx)

        # Basic valence check
        if atom1.GetExplicitValence() >= atom1.GetNoImplicit() or \
           atom2.GetExplicitValence() >= atom2.GetNoImplicit():
            return False, "Connection atoms may have valence issues"

        return True, "Connection is valid"

    def validate_structure(self, molecule):
        """
        Validate overall molecular structure.

        Args:
            molecule: Molecular structure to validate

        Returns:
            tuple: (is_valid, error_message)
        """
        if molecule is None:
            return False, "Molecule is None"

        try:
            # Basic sanity checks
            if molecule.GetNumAtoms() == 0:
                return False, "Molecule has no atoms"

            # Check if molecule can be sanitized
            mol_copy = Chem.Mol(molecule)
            try:
                Chem.SanitizeMol(mol_copy)
            except:
                return False, "Molecule failed sanitization"

            # Check for disconnected fragments
            if len(Chem.GetMolFrags(mol_copy)) > 1:
                return False, "Molecule contains disconnected fragments"

            return True, "Valid molecular structure"

        except Exception as e:
            return False, f"Validation error: {str(e)}"

    def validate_fragment_data(self, fragment_data):
        """
        Validate complete fragment data structure.

        Args:
            fragment_data: Fragment data dictionary

        Returns:
            tuple: (is_valid, error_message)
        """
        required_fields = ['name', 'smiles', 'connection_atoms', 'mol']
        for field in required_fields:
            if field not in fragment_data:
                return False, f"Missing required field: {field}"

        if not fragment_data['name'].strip():
            return False, "Fragment name cannot be empty"

        if not fragment_data['connection_atoms']:
            return False, "Connection atoms list cannot be empty"

        # Validate molecule structure
        valid_struct, struct_msg = self.validate_structure(fragment_data['mol'])
        if not valid_struct:
            return False, f"Invalid molecule structure: {struct_msg}"

        # Validate connection atoms
        valid_atoms, atoms_msg = self.validate_connection_points(fragment_data, fragment_data['connection_atoms'])
        if not valid_atoms:
            return False, f"Invalid connection atoms: {atoms_msg}"

        return True, "Fragment data is valid"

    def add_validation_rule(self, rule_name, rule_function):
        """Add custom validation rule."""
        self.validation_rules[rule_name] = rule_function

    def apply_custom_rules(self, molecule, rule_names=None):
        """
        Apply custom validation rules to molecule.

        Args:
            molecule: Molecule to validate
            rule_names: List of rule names to apply (None for all)

        Returns:
            tuple: (is_valid, error_message)
        """
        if not self.validation_rules:
            return True, "No custom rules defined"

        rules_to_apply = rule_names if rule_names else self.validation_rules.keys()

        for rule_name in rules_to_apply:
            if rule_name in self.validation_rules:
                is_valid, message = self.validation_rules[rule_name](molecule)
                if not is_valid:
                    return False, f"Rule '{rule_name}' failed: {message}"

        return True, "All custom rules passed"
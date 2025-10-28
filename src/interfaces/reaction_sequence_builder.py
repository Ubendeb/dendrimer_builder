class ReactionSequenceBuilder:
    """Handles building of reaction sequences."""

    def __init__(self, molecule_renderer, fragment_manager, chemistry_operations):
        self.molecule_renderer = molecule_renderer
        self.fragment_manager = fragment_manager
        self.chemistry_operations = chemistry_operations

    def build_reaction_sequence(self, core_mol, reaction_steps, progress_tracker):
        """Build complete reaction sequence."""
        reaction_sequence = []

        # Find branch connection points in core
        branch_points = self._find_connection_points(core_mol, "to_branch")

        current_molecules = [(core_mol, 0, -1)]  # (mol, generation, branch)

        for step_idx, step_fragment in enumerate(reaction_steps):
            # Update progress
            if progress_tracker:
                progress_tracker.update_progress(f"Step_{step_idx}_{step_fragment}")

            step_data = self.fragment_manager.fragments[step_fragment]
            step_mol = step_data['molecule']

            next_generation_molecules = []

            for mol, generation, branch in current_molecules:
                # Find branch connection points in current molecule
                mol_branch_points = self._find_connection_points(mol, "to_branch")

                for branch_point_idx in mol_branch_points:
                    # Find core connection point in step fragment
                    step_core_points = self._find_connection_points(step_mol, "to_core")
                    if not step_core_points:
                        continue

                    step_core_point = step_core_points[0]  # Use first core connection point

                    # Create reaction
                    new_mol = self.chemistry_operations.connect(
                        base_mol=mol,
                        base_mol_connection_atom_idx=branch_point_idx,
                        additional_mol=step_mol,
                        additional_mol_connection_atom_idx=step_core_point,
                        generation=generation + 1,
                        branch=branch_point_idx if branch == -1 else branch,
                        step=step_idx
                    )

                    # Store reaction info
                    reaction_info = {
                        'generation': generation,
                        'branch': branch,
                        'step': step_idx,
                        'base_mol': mol,
                        'additional_mol': step_mol,
                        'result_mol': new_mol,
                        'connection_points': {
                            'base': branch_point_idx,
                            'additional': step_core_point
                        }
                    }
                    reaction_sequence.append(reaction_info)

                    next_generation_molecules.append((new_mol, generation + 1, branch_point_idx))

            current_molecules = next_generation_molecules

        return reaction_sequence, branch_points

    def _find_connection_points(self, mol, connection_type):
        """Find atoms with specified connection type."""
        connection_points = []
        for atom in mol.GetAtoms():
            if (atom.HasProp("connection_type") and
                    atom.GetProp("connection_type") == connection_type):
                connection_points.append(atom.GetIdx())
        return connection_points

"""
Dendrimer Reaction Builder Interface
"""
import ipywidgets as widgets
from IPython.core.display_functions import clear_output, display
from rdkit import Chem
from rdkit.Chem import Draw


class DendrimerReactionBuilder:
    """Interface for building dendrimer reaction sequences from fragments."""

    def __init__(self, fragment_manager, chemistry_operations, molecule_renderer=None, progress_tracker=None):
        self.fragment_manager = fragment_manager
        self.chemistry_operations = chemistry_operations
        self.molecule_renderer = molecule_renderer
        self.progress_tracker = progress_tracker
        self.core_fragment = None
        self.reaction_steps = []
        self.reaction_sequence = []

        self._create_widgets()
        self._setup_event_handlers()

    def _create_widgets(self):
        """Create reaction builder widgets."""
        # Core selection
        self.core_select = widgets.Dropdown(
            options=[],
            description='Core:',
            layout=widgets.Layout(width='300px')
        )

        # Reaction steps selection
        self.steps_select = widgets.SelectMultiple(
            options=[],
            description='Steps:',
            layout=widgets.Layout(width='300px', height='150px')
        )

        # Refresh button
        self.refresh_btn = widgets.Button(description='Refresh Library', button_style='warning')

        # Control buttons
        self.add_step_btn = widgets.Button(description='Add Step')
        self.remove_step_btn = widgets.Button(description='Remove Step')
        self.clear_steps_btn = widgets.Button(description='Clear Steps')

        # Build controls
        self.preview_btn = widgets.Button(description='Preview Reactions', button_style='info')
        self.execute_btn = widgets.Button(description='Execute Reactions', button_style='success')
        self.visualize_btn = widgets.Button(description='Visualize Reactions', button_style='primary')

        # Display areas
        self.reaction_preview = widgets.Output(
            layout={'border': '1px solid gray', 'margin_top': '10px', 'max_height': '400px', 'overflow_y': 'auto'}
        )
        self.molecule_display = widgets.Output(
            layout={'border': '1px solid gray', 'margin_top': '10px'}
        )
        self.status_output = widgets.Output(
            layout={'border': '1px solid gray', 'margin_top': '10px'}
        )

    def _setup_event_handlers(self):
        """Setup button event handlers."""
        self.refresh_btn.on_click(self._on_refresh_clicked)
        self.add_step_btn.on_click(self._on_add_step_clicked)
        self.remove_step_btn.on_click(self._on_remove_step_clicked)
        self.clear_steps_btn.on_click(self._on_clear_steps_clicked)
        self.preview_btn.on_click(self._on_preview_clicked)
        self.execute_btn.on_click(self._on_execute_clicked)
        self.visualize_btn.on_click(self._on_visualize_clicked)

    def _on_refresh_clicked(self, btn):
        """Handle refresh button click."""
        with self.status_output:
            clear_output()
            self.update_fragment_lists()
            print("Library refreshed")

    def update_fragment_lists(self):
        """Update dropdown and select lists with available fragments from library."""
        fragments = list(self.fragment_manager.fragments.keys())
        self.core_select.options = fragments
        self.steps_select.options = fragments

        with self.status_output:
            print(f"Available fragments in library: {len(fragments)}")
            for frag in fragments:
                frag_data = self.fragment_manager.fragments[frag]
                print(f"  - {frag}: {frag_data['smiles']} (connection points: {frag_data['connection_points']})")

    def _on_add_step_clicked(self, btn):
        """Handle add step button click."""
        selected = self.steps_select.value
        with self.status_output:
            clear_output()
            for frag in selected:
                if frag not in self.reaction_steps:
                    self.reaction_steps.append(frag)
                    print(f"Added step: {frag}")
            self._update_steps_display()

    def _on_remove_step_clicked(self, btn):
        """Handle remove step button click."""
        selected = self.steps_select.value
        with self.status_output:
            clear_output()
            for frag in selected:
                if frag in self.reaction_steps:
                    self.reaction_steps.remove(frag)
                    print(f"Removed step: {frag}")
            self._update_steps_display()

    def _on_clear_steps_clicked(self, btn):
        """Handle clear steps button click."""
        with self.status_output:
            clear_output()
            self.reaction_steps.clear()
            print("All steps cleared")
            self._update_steps_display()

    def _update_steps_display(self):
        """Update steps display."""
        with self.status_output:
            print(f"Current steps: {self.reaction_steps}")

    def _find_connection_points(self, mol, connection_type):
        """Find atoms with specified connection type."""
        connection_points = []
        for atom in mol.GetAtoms():
            if (atom.HasProp("connection_type") and
                    atom.GetProp("connection_type") == connection_type):
                connection_points.append(atom.GetIdx())
        return connection_points

    def _on_preview_clicked(self, btn):
        """Preview reaction sequence."""
        with self.reaction_preview:
            clear_output()
            with self.status_output:
                clear_output()

            if not self.core_select.value:
                print("Please select a core fragment")
                return

            if not self.reaction_steps:
                print("Please add reaction steps")
                return

            # Initialize progress tracker if available
            if self.progress_tracker:
                self.progress_tracker.start_construction(len(self.reaction_steps))

            self.core_fragment = self.core_select.value
            print(f"Building reaction sequence:")
            print(f"Core: {self.core_fragment}")
            print(f"Steps: {self.reaction_steps}")

            # Get core molecule from library
            core_data = self.fragment_manager.fragments[self.core_fragment]
            core_mol = core_data['molecule']

            # Set initial properties for core
            for atom in core_mol.GetAtoms():
                atom.SetProp("generation", "0")
                atom.SetProp("branch", "-1")
                atom.SetProp("step", "-1")

            # Find branch connection points in core
            branch_points = self._find_connection_points(core_mol, "to_branch")
            print(f"Found {len(branch_points)} branch connection points in core")

            current_molecules = [(core_mol, 0, -1)]  # (mol, generation, branch)
            self.reaction_sequence = []

            for step_idx, step_fragment in enumerate(self.reaction_steps):
                print(f"\n--- Step {step_idx}: {step_fragment} ---")

                # Update progress
                if self.progress_tracker:
                    self.progress_tracker.update_progress(f"Step_{step_idx}_{step_fragment}")

                step_data = self.fragment_manager.fragments[step_fragment]
                step_mol = step_data['molecule']

                next_generation_molecules = []

                for mol, generation, branch in current_molecules:
                    # Find branch connection points in current molecule
                    mol_branch_points = self._find_connection_points(mol, "to_branch")
                    print(f"Generation {generation}, branch {branch}: {len(mol_branch_points)} branch points")

                    for branch_point_idx in mol_branch_points:
                        # Find core connection point in step fragment
                        step_core_points = self._find_connection_points(step_mol, "to_core")
                        if not step_core_points:
                            print(f"Warning: No core connection points in {step_fragment}")
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
                        self.reaction_sequence.append(reaction_info)

                        next_generation_molecules.append((new_mol, generation + 1, branch_point_idx))

                        print(
                            f"  Reaction: branch point {branch_point_idx} + {step_fragment} -> generation {generation + 1}")

                current_molecules = next_generation_molecules

            # Display preview using enhanced visualization
            self._display_reaction_preview()

    def _display_reaction_preview(self):
        """Display reaction sequence preview using enhanced visualization."""
        if not self.reaction_sequence:
            print("No reactions to preview")
            return

        print(f"\n=== REACTION SEQUENCE PREVIEW ===")
        print(f"Total reactions: {len(self.reaction_sequence)}")

        # Generate progress report if tracker available
        if self.progress_tracker:
            progress_report = self.progress_tracker.generate_progress_report()
            print(f"Construction Progress: {progress_report['progress_percentage']:.1f}% complete")
            print(f"Steps completed: {progress_report['current_step']}/{progress_report['total_steps']}")

        for i, reaction in enumerate(self.reaction_sequence):
            print(f"\nReaction {i}:")
            print(f"  Generation: {reaction['generation']} -> {reaction['generation'] + 1}")
            print(f"  Branch: {reaction['branch']}")
            print(f"  Step: {reaction['step']}")
            print(
                f"  Connection: base atom {reaction['connection_points']['base']} + additional atom {reaction['connection_points']['additional']}")

    def _on_visualize_clicked(self, btn):
        """Visualize reaction sequence using molecule renderer."""
        with self.molecule_display:
            clear_output()
            with self.status_output:
                clear_output()

            if not self.reaction_sequence:
                print("No reaction sequence to visualize. Please preview first.")
                return

            if not self.molecule_renderer:
                print("Molecule renderer not available")
                return

            print("Visualizing reaction sequence...")

            # Visualize key reactions
            for i, reaction in enumerate(self.reaction_sequence[:5]):  # Limit to first 5 for performance
                print(f"\n--- Visualizing Reaction {i} ---")

                # Render base molecule
                base_img = self.molecule_renderer.render_2d(
                    reaction['base_mol'],
                    title=f"Base Mol - Gen {reaction['generation']}, Branch {reaction['branch']}",
                    highlight_atoms=[reaction['connection_points']['base']],
                    display_props=['generation', 'branch', 'step']
                )

                # Render additional molecule
                additional_img = self.molecule_renderer.render_2d(
                    reaction['additional_mol'],
                    title=f"Additional Mol - {self.reaction_steps[reaction['step']]}",
                    highlight_atoms=[reaction['connection_points']['additional']],
                    display_props=['connection_type']
                )

                # Render result molecule
                result_img = self.molecule_renderer.render_2d(
                    reaction['result_mol'],
                    title=f"Result Mol - Gen {reaction['generation'] + 1}",
                    display_props=['generation', 'branch', 'step']
                )

                # Display all three
                display(widgets.HBox([
                    widgets.Image(value=base_img.data if hasattr(base_img, 'data') else self._pil_to_widget_image(base_img),
                                layout=widgets.Layout(width='300px')),
                    widgets.Image(value=additional_img.data if hasattr(additional_img, 'data') else self._pil_to_widget_image(additional_img),
                                layout=widgets.Layout(width='300px')),
                    widgets.Image(value=result_img.data if hasattr(result_img, 'data') else self._pil_to_widget_image(result_img),
                                layout=widgets.Layout(width='300px'))
                ]))

    def _pil_to_widget_image(self, pil_img):
        """Convert PIL image to widget-compatible format."""
        from io import BytesIO
        import base64

        buffer = BytesIO()
        pil_img.save(buffer, format='PNG')
        buffer.seek(0)

        image_data = base64.b64encode(buffer.read()).decode()
        return f"data:image/png;base64,{image_data}"

    def _on_execute_clicked(self, btn):
        """Execute the reaction sequence."""
        with self.molecule_display:
            clear_output()
            with self.status_output:
                clear_output()

            if not self.reaction_sequence:
                print("No reaction sequence to execute. Please preview first.")
                return

            print("Executing reaction sequence...")

            # Get the final molecules from reaction sequence
            final_molecules = []
            for reaction in self.reaction_sequence:
                if reaction['generation'] == len(self.reaction_steps) - 1:  # Last generation
                    final_molecules.append(reaction['result_mol'])

            if final_molecules:
                print(f"Generated {len(final_molecules)} final molecules")

                # Use molecule renderer if available, otherwise fall back to RDKit
                if self.molecule_renderer:
                    molecules_to_display = min(6, len(final_molecules))  # Limit for performance
                    for i, mol in enumerate(final_molecules[:molecules_to_display]):
                        img = self.molecule_renderer.render_2d(
                            mol,
                            title=f"Final_Mol_{i}",
                            display_props=['generation', 'branch', 'step']
                        )
                        display(img)
                else:
                    # Fall back to original RDKit display
                    img = Draw.MolsToGridImage(
                        final_molecules,
                        molsPerRow=3,
                        subImgSize=(300, 300),
                        legends=[f"Mol_{i}" for i in range(len(final_molecules))]
                    )
                    display(img)
            else:
                print("No final molecules generated")

    def display_interface(self):
        """Display the reaction builder interface."""
        display(widgets.HTML("<h3>Dendrimer Reaction Builder</h3>"))

        # Refresh button
        display(widgets.HBox([self.refresh_btn]))

        # Core selection
        display(widgets.HTML("<strong>Select Core Fragment:</strong>"))
        display(self.core_select)

        # Steps selection
        display(widgets.HTML("<strong>Select Reaction Steps:</strong>"))
        display(widgets.HBox([
            self.steps_select,
            widgets.VBox([
                self.add_step_btn,
                self.remove_step_btn,
                self.clear_steps_btn
            ])
        ]))

        # Control buttons
        display(widgets.HBox([
            self.preview_btn,
            self.visualize_btn,
            self.execute_btn
        ]))

        # Display areas
        display(widgets.HTML("<strong>Reaction Preview:</strong>"))
        display(self.reaction_preview)

        display(widgets.HTML("<strong>Molecule Display:</strong>"))
        display(self.molecule_display)

        display(widgets.HTML("<strong>Status:</strong>"))
        display(self.status_output)
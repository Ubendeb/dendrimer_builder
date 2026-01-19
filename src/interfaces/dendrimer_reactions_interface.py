"""
Dendrimer Reaction Builder Interface
"""
import ipywidgets as widgets
from IPython.core.display_functions import clear_output, display

from src.interfaces.reaction_sequence_builder import ReactionSequenceBuilder
from src.visualization.reaction_renderer import ReactionVisualizer
from src.visualization.progress_tracker import ProgressTracker
from src.core.chemistry_operations import ChemistryOperations


class DendrimerReactionBuilder:
    """Interface for building dendrimer reaction sequences from fragments."""

    def __init__(self, fragment_manager, chemistry_operations, molecule_renderer):
        """

        :type chemistry_operations: ChemistryOperations
        """
        self.fragment_manager = fragment_manager
        self.chemistry_operations = chemistry_operations
        self.molecule_renderer = molecule_renderer
        self.reaction_visualizer = ReactionVisualizer(molecule_renderer)
        self.progress_tracker = ProgressTracker()
        self.sequence_builder = ReactionSequenceBuilder(molecule_renderer, fragment_manager, chemistry_operations)

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
        self.visualize_btn = widgets.Button(
            description='Visualize Reactions',
            button_style='primary',
            disabled=not self.molecule_renderer  # Disable if no renderer
        )

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

            self.core_fragment = self.core_select.value

            # Initialize progress tracker
            self.progress_tracker.start_construction(len(self.reaction_steps))

            print(f"Building reaction sequence:")
            print(f"Core: {self.core_fragment}")
            print(f"Steps: {self.reaction_steps}")

            # Get core molecule from library
            core_data = self.fragment_manager.fragments[self.core_fragment]
            self.core_mol = core_data['molecule']
            self.chemistry_operations.set_atom_properties(self.core_mol,0,-1,-1)


            # Build reaction sequence
            self.reaction_sequence, branch_points = self.sequence_builder.build_reaction_sequence(
                self.core_mol,
                self.reaction_steps,
                self.progress_tracker
            )

            print(f"Found {len(branch_points)} branch connection points in core")

            # Display preview
            self._display_reaction_preview()

    def _display_reaction_preview(self):
        """Display reaction sequence preview."""
        if not self.reaction_sequence:
            print("No reactions to preview")
            return

        print(f"\n=== REACTION SEQUENCE PREVIEW ===")
        print(f"Total reactions: {len(self.reaction_sequence)}")

        # Generate progress report
        progress_report = self.progress_tracker.generate_progress_report()
        print(f"Construction Progress: {progress_report['progress_percentage']:.1f}% complete")
        print(f"Steps completed: {progress_report['current_step']}/{progress_report['total_steps']}")

        # Show reaction statistics by generation
        reactions_by_gen = {}
        for reaction in self.reaction_sequence:
            gen = reaction['generation']
            if gen not in reactions_by_gen:
                reactions_by_gen[gen] = 0
            reactions_by_gen[gen] += 1

        for gen in sorted(reactions_by_gen.keys()):
            print(f"Generation {gen}: {reactions_by_gen[gen]} reactions")

        # Show first few reactions in detail
        for i, reaction in enumerate(self.reaction_sequence[:3]):
            print(f"\nSample Reaction {i}:")
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

            if not self.reaction_visualizer:
                print("Molecule renderer not available for visualization")
                return

            self.reaction_visualizer.visualize_reaction_sequence(self.reaction_sequence, self.reaction_steps)

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


                molecules_to_display = min(6, len(final_molecules))  # Limit for performance
                for i, mol in enumerate(final_molecules[:molecules_to_display]):
                    img = self.molecule_renderer.render_2d(
                        mol,
                        title=f"Final_Mol_{i}",
                        display_props=['generation', 'branch', 'step']
                    )
                    if img:
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
        control_buttons = [self.preview_btn, self.execute_btn]
        if self.molecule_renderer:
            control_buttons.insert(1, self.visualize_btn)

        display(widgets.HBox(control_buttons))

        # Display areas
        display(widgets.HTML("<strong>Reaction Preview:</strong>"))
        display(self.reaction_preview)

        display(widgets.HTML("<strong>Molecule Display:</strong>"))
        display(self.molecule_display)

        display(widgets.HTML("<strong>Status:</strong>"))
        display(self.status_output)
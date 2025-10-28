import ipywidgets as widgets
from IPython.core.display_functions import display


class ReactionVisualizer:
    """Handles visualization of reaction sequences using MoleculeRenderer."""

    def __init__(self, molecule_renderer):
        self.molecule_renderer = molecule_renderer

    def visualize_reaction_sequence(self, reaction_sequence, reaction_steps):
        """Visualize complete reaction sequence."""
        if not reaction_sequence:
            print("No reaction sequence to visualize")
            return

        print("Visualizing reaction sequence...")

        # Group reactions by generation
        reactions_by_gen = {}
        for reaction in reaction_sequence:
            gen = reaction['generation']
            if gen not in reactions_by_gen:
                reactions_by_gen[gen] = []
            reactions_by_gen[gen].append(reaction)

        # Visualize each generation
        for gen in sorted(reactions_by_gen.keys()):
            print(f"\n--- Generation {gen} Reactions ---")
            reactions = reactions_by_gen[gen][:3]  # Limit for performance

            for i, reaction in enumerate(reactions):
                self._visualize_single_reaction(reaction, reaction_steps, f"Gen{gen}_React{i}")

    def _visualize_single_reaction(self, reaction, reaction_steps, title_prefix):
        """Visualize a single reaction."""
        # Render base molecule
        base_img = self.molecule_renderer.render_2d(
            reaction['base_mol'],
            title=f"{title_prefix}_Base",
            highlight_atoms=[reaction['connection_points']['base']],
            display_props=['generation', 'branch', 'step']
        )

        # Render additional molecule
        step_fragment = reaction_steps[reaction['step']]
        additional_img = self.molecule_renderer.render_2d(
            reaction['additional_mol'],
            title=f"{title_prefix}_Add_{step_fragment}",
            highlight_atoms=[reaction['connection_points']['additional']],
            display_props=['connection_type']
        )

        # Render result molecule
        result_img = self.molecule_renderer.render_2d(
            reaction['result_mol'],
            title=f"{title_prefix}_Result",
            highlight_atoms=[reaction['connection_points']['base']],
            display_props=['generation', 'branch', 'step']
        )

        # Display all three
        self._display_reaction_images([base_img, additional_img, result_img])

    def _display_reaction_images(self, images):
        """Display reaction images in a horizontal layout."""
        from io import BytesIO
        import base64

        image_widgets = []
        for img in images:
            if img is not None:
                buffer = BytesIO()
                img.save(buffer, format='PNG')
                buffer.seek(0)

                image_data = base64.b64encode(buffer.read()).decode()
                widget_image = widgets.Image(
                    value=f"data:image/png;base64,{image_data}",
                    layout=widgets.Layout(width='300px', height='200px')
                )
                image_widgets.append(widget_image)

        if image_widgets:
            display(widgets.HBox(image_widgets))

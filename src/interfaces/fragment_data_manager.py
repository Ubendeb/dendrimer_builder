import json

import ipywidgets as widgets
from IPython.core.display_functions import clear_output, display
from rdkit import Chem

from src import DataSaver, DataLoader


class FragmentDataManager:
    """Manages fragment data storage and retrieval operations."""

    def __init__(self, fragment_builder, fragment_manager):
        self.builder = fragment_builder
        self.manager = fragment_manager
        self._create_widgets()
        self._setup_event_handlers()

    def _create_widgets(self):
        """Create data management widgets."""
        self.filename_input = widgets.Text(
            value='dendrimer_fragments.json',
            placeholder='File name',
            layout=widgets.Layout(width='200px')
        )

        self.save_btn = widgets.Button(description='Save Fragments')
        self.load_btn = widgets.Button(description='Load Fragments')
        self.export_library_btn = widgets.Button(description='Export to Library')
        self.import_library_btn = widgets.Button(description='Import from Library')

        self.status_output = widgets.Output(
            layout={'border': '1px solid gray', 'margin_top': '10px', 'max_height': '150px', 'overflow_y': 'auto'}
        )

    def _setup_event_handlers(self):
        """Setup button event handlers."""
        self.save_btn.on_click(self._on_save_clicked)
        self.load_btn.on_click(self._on_load_clicked)
        self.export_library_btn.on_click(self._on_export_library_clicked)
        self.import_library_btn.on_click(self._on_import_library_clicked)

    def _on_save_clicked(self, btn):
        """Handle save button click."""
        with self.status_output:
            clear_output()
            self.save_fragments_to_file(self.filename_input.value)

    def _on_load_clicked(self, btn):
        """Handle load button click."""
        with self.status_output:
            clear_output()
            self.load_fragments_from_file(self.filename_input.value)

    def _on_export_library_clicked(self, btn):
        """Handle export to library button click."""
        with self.status_output:
            clear_output()
            self.export_to_fragment_manager()

    def _on_import_library_clicked(self, btn):
        """Handle import from library button click."""
        with self.status_output:
            clear_output()
            self.import_from_fragment_manager()

    def save_fragments_to_file(self, filename='dendrimer_fragments.json'):
        """Save fragments to JSON file."""
        fragments_data = self.builder.get_fragments_data()
        DataSaver.save_fragment_data(fragments_data, filename)

    def load_fragments_from_file(self, filename='dendrimer_fragments.json'):
        """Load fragments from JSON file."""
        try:
            fragments_data = DataLoader.load_fragment_data(filename)
            self.builder.fragments = fragments_data
            self.builder.update_fragments_display()
            print(f"Data loaded from file: {filename}")
            print(f"Fragments loaded: {len(fragments_data)}")

        except FileNotFoundError:
            print(f"File not found: {filename}")
        except Exception as e:
            print(f"Loading error: {e}")

    def export_to_fragment_manager(self):
        """Export fragments from builder to fragment manager."""
        fragments_data = self.builder.get_fragments_data()

        if not fragments_data:
            print("No fragments to export")
            return

        exported_count = 0
        for frag in fragments_data:
            fragment_id = self.manager.add_fragment(
                fragment=frag['mol'],
                name=frag['name'],
                connection_atoms=frag['connection_atoms'],
                replacement_groups=frag['replacement_groups']
            )
            exported_count += 1
            print(f"Exported fragment: {fragment_id}")

        print(f"Total fragments exported: {exported_count}")

    def import_from_fragment_manager(self):
        """Import fragments from fragment manager to builder."""
        if not self.manager.fragments:
            print("No fragments in manager library")
            return

        imported_count = 0
        for frag_id, frag_data in self.manager.fragments.items():
            fragment_data = {
                'name': frag_id,
                'smiles': frag_data['smiles'],
                'connection_atoms': frag_data['connection_points'],
                'replacement_groups': frag_data['replacement_groups'],
                'num_atoms': frag_data['num_atoms'],
                'mol': frag_data['molecule']
            }
            self.builder.fragments.append(fragment_data)
            imported_count += 1
            print(f"Imported fragment: {frag_id}")

        self.builder.update_fragments_display()
        print(f"Total fragments imported: {imported_count}")

    def display_interface(self):
        """Display the data management interface."""
        display(widgets.HTML("<h3>Data Management</h3>"))
        display(widgets.HBox([
            self.filename_input,
            self.save_btn,
            self.load_btn,
            self.export_library_btn,
            self.import_library_btn
        ]))
        display(widgets.HTML("<strong>Operation status:</strong>"))
        display(self.status_output)
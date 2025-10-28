"""
Interactive interface for molecular fragment input and management.
"""

import ipywidgets as widgets
from IPython.display import display, clear_output
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw

from src import StructureValidator


class DendrimerFragmentBuilder:
    """Interactive builder for dendrimer fragments with visualization."""

    def __init__(self):
        self.fragments = []
        self.current_fragment = {}
        self.validator = StructureValidator()
        self.editing_index = None  # Track which fragment is being edited
        self._create_widgets()
        self._setup_event_handlers()

    def _create_widgets(self):
        """Create all interface widgets."""
        # SMILES input section
        self.smiles_input = widgets.Text(
            placeholder='Enter SMILES (e.g., CCO for ethanol)',
            layout=widgets.Layout(width='400px')
        )
        self.validate_btn = widgets.Button(
            description='Validate SMILES',
            button_style='primary',
            layout=widgets.Layout(width='150px')
        )
        self.mol_display = widgets.Output()

        # Fragment properties
        self.fragment_name = widgets.Text(
            placeholder='Fragment name',
            layout=widgets.Layout(width='400px')
        )
        self.connection_atoms = widgets.Text(
            placeholder='Connection atoms (comma-separated, e.g.: 0,1,2)',
            layout=widgets.Layout(width='400px')
        )
        self.replacement_groups = widgets.Text(
            placeholder='Replacement groups (comma-separated, e.g.: H,OH,Cl)',
            layout=widgets.Layout(width='400px')
        )

        # Action buttons
        self.add_fragment_btn = widgets.Button(
            description='Add Fragment',
            button_style='success',
            layout=widgets.Layout(width='150px')
        )
        self.update_fragment_btn = widgets.Button(
            description='Update Fragment',
            button_style='warning',
            layout=widgets.Layout(width='150px')
        )
        self.cancel_edit_btn = widgets.Button(
            description='Cancel Edit',
            button_style='',
            layout=widgets.Layout(width='150px')
        )
        self.clear_btn = widgets.Button(
            description='Clear Form',
            layout=widgets.Layout(width='150px')
        )

        # Button container
        self.button_container = widgets.HBox([
            self.add_fragment_btn,
            self.update_fragment_btn,
            self.cancel_edit_btn,
            self.clear_btn
        ])

        # Display areas
        self.fragments_display = widgets.Output()
        self.status_output = widgets.Output()

    def _setup_event_handlers(self):
        """Setup widget event handlers."""
        self.validate_btn.on_click(self.validate_smiles)
        self.add_fragment_btn.on_click(self.add_fragment)
        self.update_fragment_btn.on_click(self.update_fragment)
        self.cancel_edit_btn.on_click(self.cancel_edit)
        self.clear_btn.on_click(self.clear_form)

        # Initially hide update and cancel buttons
        self.update_fragment_btn.layout.visibility = 'hidden'
        self.cancel_edit_btn.layout.visibility = 'hidden'

    def draw_molecule_with_atom_indices(self, mol, size=(400, 300)):
        """Draw molecule with atom index labels."""
        mol_for_drawing = Chem.Mol(mol)

        for atom in mol_for_drawing.GetAtoms():
            atom.SetProp('atomNote', str(atom.GetIdx()))

        return Draw.MolToImage(mol_for_drawing, size=size, kekulize=True)

    def validate_smiles(self, btn):
        """Validate SMILES string and display molecule."""
        with self.status_output:
            clear_output()
            smiles = self.smiles_input.value.strip()

            if not smiles:
                print("Please enter SMILES string")
                return

            # Использование StructureValidator
            is_valid, mol, message = self.validator.validate_smiles(smiles)
            if not is_valid:
                print(f"SMILES validation failed: {message}")
                return

            self.current_fragment.update({
                'smiles': smiles,
                'mol': mol,
                'num_atoms': mol.GetNumAtoms()
            })

            with self.mol_display:
                clear_output()
                img = self.draw_molecule_with_atom_indices(mol)
                display(img)

                print("\nAtom Information:")
                print("-" * 40)
                for atom in mol.GetAtoms():
                    atom_idx = atom.GetIdx()
                    atom_symbol = atom.GetSymbol()
                    neighbors = [neighbor.GetIdx() for neighbor in atom.GetNeighbors()]
                    print(f"Atom {atom_idx}: {atom_symbol} (neighbors: {neighbors})")

            print(f"Valid SMILES: {smiles}")
            print(f"Number of atoms: {mol.GetNumAtoms()}")
            print(f"Molecular formula: {Descriptors.CalcMolFormula(mol)}")

    def add_fragment(self, btn):
        """Add fragment to the collection."""
        with self.status_output:
            clear_output()

            if 'smiles' not in self.current_fragment:
                print("Please validate SMILES first")
                return

            name = self.fragment_name.value.strip()
            if not name:
                print("Please enter fragment name")
                return

            conn_atoms = self.connection_atoms.value.strip()
            if not conn_atoms:
                print("Please specify connection atoms")
                return

            try:
                connection_list = [int(x.strip()) for x in conn_atoms.split(',')]

                # Использование StructureValidator для проверки атомов
                valid_atoms, atoms_msg = self.validator.validate_atom_indices(
                    self.current_fragment['mol'], connection_list
                )
                if not valid_atoms:
                    print(f"Connection atoms validation failed: {atoms_msg}")
                    return

            except ValueError:
                print("Invalid connection atoms format")
                return

            replacement_list = [x.strip() for x in
                                self.replacement_groups.value.split(',')] if self.replacement_groups.value else []

            fragment_data = {
                'name': name,
                'smiles': self.current_fragment['smiles'],
                'connection_atoms': connection_list,
                'replacement_groups': replacement_list,
                'num_atoms': self.current_fragment['num_atoms'],
                'mol': self.current_fragment['mol']
            }

            # Финальная валидация всех данных фрагмента
            is_valid, validation_msg = self.validator.validate_fragment_data(fragment_data)
            if not is_valid:
                print(f"Fragment validation failed: {validation_msg}")
                return

            self.fragments.append(fragment_data)
            print(f"Fragment '{name}' added successfully")
            self.update_fragments_display()
            self.clear_form()

    def edit_fragment(self, index):
        """Load fragment data into form for editing."""
        if index < 0 or index >= len(self.fragments):
            return

        fragment = self.fragments[index]
        self.editing_index = index

        # Fill form with fragment data
        self.smiles_input.value = fragment['smiles']
        self.fragment_name.value = fragment['name']
        self.connection_atoms.value = ','.join(map(str, fragment['connection_atoms']))
        self.replacement_groups.value = ','.join(fragment['replacement_groups'])

        # Update current fragment and validate
        self.current_fragment.update({
            'smiles': fragment['smiles'],
            'mol': fragment['mol'],
            'num_atoms': fragment['num_atoms']
        })

        # Show molecule
        with self.mol_display:
            clear_output()
            img = self.draw_molecule_with_atom_indices(fragment['mol'])
            display(img)

        # Switch to edit mode
        self.add_fragment_btn.layout.visibility = 'hidden'
        self.update_fragment_btn.layout.visibility = 'visible'
        self.cancel_edit_btn.layout.visibility = 'visible'

        with self.status_output:
            clear_output()
            print(f"Editing fragment: {fragment['name']}")

    def update_fragment(self, btn):
        """Update existing fragment."""
        if self.editing_index is None:
            return

        with self.status_output:
            clear_output()

            if 'smiles' not in self.current_fragment:
                print("Please validate SMILES first")
                return

            name = self.fragment_name.value.strip()
            if not name:
                print("Please enter fragment name")
                return

            conn_atoms = self.connection_atoms.value.strip()
            if not conn_atoms:
                print("Please specify connection atoms")
                return

            try:
                connection_list = [int(x.strip()) for x in conn_atoms.split(',')]
                valid_atoms, atoms_msg = self.validator.validate_atom_indices(
                    self.current_fragment['mol'], connection_list
                )
                if not valid_atoms:
                    print(f"Connection atoms validation failed: {atoms_msg}")
                    return
            except ValueError:
                print("Invalid connection atoms format")
                return

            replacement_list = [x.strip() for x in
                                self.replacement_groups.value.split(',')] if self.replacement_groups.value else []

            fragment_data = {
                'name': name,
                'smiles': self.current_fragment['smiles'],
                'connection_atoms': connection_list,
                'replacement_groups': replacement_list,
                'num_atoms': self.current_fragment['num_atoms'],
                'mol': self.current_fragment['mol']
            }

            is_valid, validation_msg = self.validator.validate_fragment_data(fragment_data)
            if not is_valid:
                print(f"Fragment validation failed: {validation_msg}")
                return

            # Update the fragment
            self.fragments[self.editing_index] = fragment_data
            print(f"Fragment '{name}' updated successfully")
            self.update_fragments_display()
            self.cancel_edit()

    def cancel_edit(self, btn=None):
        """Cancel editing mode."""
        self.editing_index = None
        self.add_fragment_btn.layout.visibility = 'visible'
        self.update_fragment_btn.layout.visibility = 'hidden'
        self.cancel_edit_btn.layout.visibility = 'hidden'
        self.clear_form()

    def delete_fragment(self, index):
        """Delete fragment from collection."""
        if 0 <= index < len(self.fragments):
            fragment_name = self.fragments[index]['name']
            del self.fragments[index]

            with self.status_output:
                clear_output()
                print(f"Fragment '{fragment_name}' deleted successfully")

            self.update_fragments_display()

            # If we were editing the deleted fragment, cancel edit mode
            if self.editing_index == index:
                self.cancel_edit()

    def clear_form(self, btn=None):
        """Clear the input form."""
        self.smiles_input.value = ''
        self.fragment_name.value = ''
        self.connection_atoms.value = ''
        self.replacement_groups.value = ''
        self.current_fragment = {}

        with self.mol_display:
            clear_output()
        with self.status_output:
            if btn:  # Only print message if called by button
                clear_output()
                print("Form cleared. Ready for new fragment input.")

    def _create_fragment_card(self, fragment, index):
        """Create a styled card for fragment display with action buttons."""
        card = widgets.VBox([
            widgets.HTML(f"<b>{index + 1}. {fragment['name']}</b>"),
            widgets.HTML(f"SMILES: <code>{fragment['smiles']}</code>"),
            widgets.HTML(f"Connection atoms: {fragment['connection_atoms']}"),
            widgets.HTML(f"Replacement groups: {fragment['replacement_groups']}"),
            widgets.HTML(f"Total atoms: {fragment['num_atoms']}"),
            widgets.HBox([
                widgets.Button(
                    description='Edit',
                    button_style='warning',
                    layout=widgets.Layout(width='80px', height='30px')
                ),
                widgets.Button(
                    description='Delete',
                    button_style='danger',
                    layout=widgets.Layout(width='80px', height='30px')
                )
            ]),
            widgets.HTML("<hr style='margin: 10px 0;'>")
        ])

        # Set button handlers
        card.children[-2].children[0].on_click(lambda btn, idx=index: self.edit_fragment(idx))
        card.children[-2].children[1].on_click(lambda btn, idx=index: self.delete_fragment(idx))

        return card

    def update_fragments_display(self):
        """Update the fragments list display with interactive cards."""
        with self.fragments_display:
            clear_output()
            if not self.fragments:
                print("No fragments added yet")
                return

            print("FRAGMENTS COLLECTION:")
            print("=" * 60)

            # Create cards for all fragments
            cards = []
            for i, frag in enumerate(self.fragments):
                cards.append(self._create_fragment_card(frag, i))

            # Display all cards in a scrollable area
            display(widgets.VBox(cards, layout=widgets.Layout(
                min_height='400px',
                overflow_y='auto',
                border='1px solid #ccc',
                padding='10px'
            )))

    def get_fragments_data(self):
        """Return collected fragments data."""
        return self.fragments.copy()

    def display_interface(self):
        """Display the main interface."""
        form = widgets.VBox([
            widgets.HTML("<h3>Molecular Fragment Input for Dendrimers</h3>"),
            widgets.HTML("<b>1. Molecule SMILES:</b>"),
            widgets.HBox([self.smiles_input, self.validate_btn]),
            self.mol_display,
            widgets.HTML("<b>2. Fragment Name:</b>"),
            self.fragment_name,
            widgets.HTML("<b>3. Connection Atoms (indices):</b>"),
            widgets.HTML("<i>Atom indices start from 0. Indices shown in image above</i>"),
            self.connection_atoms,
            widgets.HTML("<b>4. Replacement Groups (optional):</b>"),
            self.replacement_groups,
            self.button_container,
            self.status_output
        ])

        display(form)
        display(widgets.HTML("<h3>Fragment List:</h3>"))
        display(self.fragments_display)
        self.update_fragments_display()

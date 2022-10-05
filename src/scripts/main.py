import ipywidgets as wg
import sys
sys.path.append('./scripts/calculators/siesta')
sys.path.append('./scripts/calculators/vasp')
sys.path.append('./images/')
sys.path.append('./src/scripts/')
sys.path.append('./src/scripts/calculators/siesta')
sys.path.append('./src/scripts/calculators/vasp')
sys.path.append('./src/images/')
sys.path.append('./scripts/')
sys.path.append('./scripts/calculators/siesta')
sys.path.append('./scripts/calculators/vasp')
sys.path.append('./images/')

import numpy as np
import os
from ipywidgets import interact, interactive, fixed, interact_manual, Button, Layout, jslink, IntText, GridBox, Label, Box
from IPython.display import display, HTML, IFrame, clear_output
from ase.io import read, write
from utils import *
from create_widgets import *
import py3Dmol
from create_fdf_file import *
from siesta_check_outputs import *
from viewer import *
from siesta import *
import pandas as pd
import ipysheet
from ipysheet import column, row
from ipysheet import sheet, cell
from ase import Atoms
from ase.build import bulk, molecule, graphene_nanoribbon
from ase.build import nanotube as ase_nanotube
from ase.collections import g2
from openbabel import pybel
from io import StringIO
from html2image import Html2Image

#Dependencias externas - siesta conda install -c conda-forge siesta=4.1.5=mpi_openmpi_hfab99a0_1


class Spin(Siesta):

    """
    Graphical User interface using Ipywg for Atomic simulation Environments:
    (SPIN: [S]imple [P]ython [I]pywidgets Interface to obtain the optoelectronic properties of [N]anostructures)
    """
    
    def __init__(self):
        
        """
        Initial window, create the main tabs of app
        """
        
        self.create_structure_tabs()

        
    def create_parameter_tabs(self):
        
        """
        Create parameter tab for each calculator.
        """
    
        self.select_SIESTA_parameters()
        
        self.supported_calculators = wg.Accordion(children = [self.siesta_box])
        self.supported_calculators.set_title(0, 'SIESTA')
        
    def structure(self):
        """
        display structure blocks
        """
        display(self.structure_viewer_tabs)
        
    def preprocessing(self):
        """
        display preprocessing blocks
        """
        
        try:
            self.structure_label = self.structure_label_input.value 
        except:
            self.structure_label = 'Structure'

        self.create_parameter_tabs()
        logo = get_calculator_logo("./src/images/spin_logo.png", "png")
        display(logo)
        display(self.supported_calculators)

        
    def create_structure_tabs(self):
        
        """
        Create initial schema for structure viewer
        """
        
        self.title_current_structure = wg.HTML(
                 value=f"<b> Structure</b>")
        self.structure_label_input = create_text_input('Structure','Label: ',  width='450px')
        self.structure_label_input.observe(self.update_structure_label_input, 'value')
        #define box layout / box style
        box_layout = Layout(display='flex',
                    flex_flow='column',
                    align_items='stretch',
                    border='doted',
                    width='100%')
        
        #create button for the manual creation of structures
        create_structure_button = create_expanded_button('Create structure','success')     
        create_structure_button.on_click(self.create_structure_manually)
        
        #create widgets for upload structures (upload cif file)

        #predefined molecules
        self.molecules_list = create_dropdown_input(g2.names, 'H','')
        self.IS_MOLECULE = False

        generate_button = create_expanded_button('Generate structure','success')
        generate_button.on_click(self.generate_molecules)
        self.molecule_box = wg.VBox([HBox([Label('Select molecule: '), self.molecules_list]),
                                         generate_button],
                                        layout = {'width':'max-content'})
                                
        self.make_bulk()
        self.make_nanotube_with_buckling()
        self.make_nanotube_without_buckling()
        self.make_graphene_nanoribons()
        self.predefined_structures = wg.Accordion(children = [self.molecule_box, self.bulk_box,
                                            self.nanotube_buckling_Box, self.nanotube_Box, self.g_n_box])
        self.predefined_structures.set_title(0, 'Molecules')
        self.predefined_structures.set_title(1, 'Creating bulk systems')
        self.predefined_structures.set_title(2, 'Creating Nanotube with buckling')
        self.predefined_structures.set_title(3, 'Creating Nanotube without buckling')
        self.predefined_structures.set_title(4, 'Creating Graphene nanoribons')

        self.load_saved_structures()
        self.created_structures = wg.Accordion(children = [self.created_structures_box])
       

        self.structure_from_file = True
        self.geometric_structure()

        self.job_folder = FileChooser('./calculations/')#,layout=Layout(height='200px', width='500px'))
        self.job_folder.default_filename = self.structure_label_input.value
        self.geometric_structure_box = wg.Accordion(children = [VBox([HBox([self.structure_input,
                                                Label('Please upload a structure file')]),
                                        self.valid_structure_input]),
                                          VBox([HBox([Label('Total number of atoms: '), self.number_of_atoms_input]),
                                            HBox([self.create_structure_manually_button])]),
                                            self.predefined_structures,
                                            self.created_structures_box])
        
        self.geometric_structure_box.set_title(0, 'Upload file')
        self.geometric_structure_box.set_title(1, 'Create structure manually')
        self.geometric_structure_box.set_title(2, 'Predefined structures')
        self.geometric_structure_box.set_title(3, 'Load created structures')
        self.valid_structure_input.layout.display = 'none'
        
        #Divide window in two parts, left part: coordinates and unit cell,  right part: Structure visualization
        self.structure_grid = {}
        self.structure_grid['grid'] = wg.GridspecLayout(1, 4)
        
        box_layout = Layout(display='flex',
                    flex_flow='column',
                    align_items='stretch',
                    align_content='stretch',
                    flex='flex-grow',
                    border='solid 1px')
        
        
        self.structure_parameters = wg.Output(layout = box_layout)
        
        with self.structure_parameters:
            clear_output(True)
            logo = get_calculator_logo("./src/images/spin_logo.png", "png")
            display(logo)
            display(self.geometric_structure_box)
            
        #Display right part elements
        self.out = wg.Output(layout = box_layout)

        #Display left part elements
        show_labels_button = create_checkbox_input('')      
        show_labels_button.observe(self.show_labels, 'value')
        self.show_structure_labels = show_labels_button.value

        select_style = create_dropdown_input(['line','cross','stick','sphere','cartoon','clicksphere'], 'sphere','')
        select_style.observe(self.change_viewer_style, 'value')
        self.viewer_style = select_style.value

        select_visualizer = create_dropdown_input(['Py3Dmol','NGLviewer'], 'Py3Dmol','')
        select_visualizer.observe(self.change_visualizer, 'value')
        self.visualizer= select_visualizer.value

        self.get_atom_by_index = create_dropdown_input([0], 0, '')
        self.get_atom_by_index.observe(self.search_atom_by_index, 'value')
        self.selected_atom = wg.Output( descripcion = 'Edit atom',layout = box_layout)

        self.add_atom_button = create_checkbox_input('') 
        self.add_atom_button.observe(self.add_atom, 'value')

        self.delete_atom_button = create_checkbox_input('') 
        self.delete_atom_button.observe(self.delete_atom, 'value')
        self.selected_atom_to_drop= wg.Output( descripcion = 'Edit atom',layout = box_layout)

        self.visualizer_box = wg.HBox([Label('Select Visualizer'), select_visualizer])
        self.viewer_options = wg.VBox([
                                       wg.HBox([Label('Show labels'), show_labels_button]),
                                       wg.HBox([Label('Style'), select_style]),
                                       wg.HBox([Label('Edit Atom'),self.get_atom_by_index]),
                                       self.selected_atom,
                                       wg.HBox([Label('Add Atom'),self.add_atom_button]),
                                       wg.HBox([Label('Delete Atom'),self.delete_atom_button]),
                                       self.selected_atom_to_drop
                                       ])
        
        self.structure_grid['grid'][0, 0] = wg.VBox(children=[self.visualizer_box,self.viewer_options])
        self.structure_grid['grid'][0, 1:] = wg.VBox(children=[self.out])


        with self.out:
            display(wg.HTML(value=f"Any structure has been loaded yet."))
        

        
        self.structure_viewer_tabs = wg.Tab(children = [self.structure_parameters,
                                                        self.structure_grid['grid']], 
                    layout = box_layout)
        self.structure_viewer_tabs.set_title(0, 'Structure options')
        self.structure_viewer_tabs.set_title(1, 'Viewer')

    def load_saved_structures(self):
        
        list_calculations = os.listdir('./calculations/')
    
        if len(list_calculations) == 0:
            self.created_structures_list = create_dropdown_input(['None'],'None','')
        else:
            self.created_structures_list  = widgets.Select(
                                options=list_calculations,
                                value=list_calculations[0],
                                description='',
                                disabled=False
                            )
        
        load_button = create_expanded_button('Load structure','success')
        self.job_load_folder = FileChooser('./calculations/')#,layout=Layout(height='200px', width='500px'))
        self.job_load_folder.default_filename = ' '
        self.update_job_load_folder_button = create_expanded_button('Update','primary')
        self.update_job_load_folder_button.on_click(self.list_available_structures)
        #generate structure
        load_button.on_click(self.load_select_saved_structure) 
        self.created_structures_box = widgets.VBox([HBox([Label('Select saved structures folder (Default "Calculations"): '),self.job_load_folder,
         self.update_job_load_folder_button]),
                                                self.created_structures_list, load_button],
                                                layout = {'width':'max-content'})
    def update_structure_label_input(self,change):
        self.job_folder.default_filename = change.new
        try:
            self.label_input.value = change.new
        except:
            pass
        
    def list_available_structures(self, change):
        try:
            folder = self.job_load_folder.value.split('/ ')[0]
            list_calculations = os.listdir(folder)

            if len(list_calculations) == 0:
                self.created_structures_list.options = ['None']
                self.created_structures_list.value = 'None'
            else:
                self.created_structures_list.options = list_calculations
                self.created_structures_list.value = [list_calculations[0]]
        except:pass

    def load_select_saved_structure(self,aux):
        try:
            folder = self.job_load_folder.value.split('/ ')[0]
        except:
            folder = './calculations'
        label = self.created_structures_list.value
        structure_directory = f'{folder}/{label}/structure/'

        
        content = ''
        try:
            with open(f'{structure_directory}{label}-R.fdf', 'r') as file:
                content = file.read()
        except:
             with open(f'{structure_directory}{label}.fdf', 'r') as file:
                content = file.read()
        content = self.convert_any_format_to_cif([{f'{label}.fdf':{'name':f'{label}.fdf',
                                    'type':'type_',
                                    'size':'size',
                                    'lastModified':'lastModified',
                                    'content':bytes(content, 'utf-8')}}, f'{label}.fdf','fdf'])

        self.structure_input = Simulated_structure(content, filename = f'{label}.cif')
        uploaded_file = self.structure_input.value
        filename = next(iter(uploaded_file))
        self.structure_input.value[filename]['content'] = str.encode(content)
        self.structure_from_file = False
        self.check_uploaded_structure_file()


    def update_content_with_predefined_structure(self):
        """
        storage predefined structure content
        """

        write('temp.cif', self.predefined_structure, wrap=False)
        content = ''
        
        with open('temp.cif', 'r') as file:
            content = file.read()
        self.structure_input = Simulated_structure(content)
        uploaded_file = self.structure_input.value
        filename = next(iter(uploaded_file))
        self.structure_input.value[filename]['content'] = str.encode(content)
        self.structure_from_file = False
        self.check_uploaded_structure_file()



    def search_atom_by_index(self, change):
        """
        Method to search and show atom by its index
        """

        atom_info = self.df_coordinates.iloc[int(change.new)]
        self.x = create_float_text(atom_info['Atomic coordinate x (Å)'], '')
        self.y =  create_float_text(atom_info['Atomic coordinate y (Å)'], '')
        self.z = create_float_text(atom_info['Atomic coordinate z (Å)'], '')
        self.atom = create_text_input(atom_info['Atom'], '', width='70px')
        self.atom_attributes = wg.VBox([wg.HBox([Label('Atomic coordinate x (Å)'), self.x]),
        wg.HBox([Label('Atomic coordinate y (Å)'), self.y]),
        wg.HBox([Label('Atomic coordinate z (Å)'), self.z]),
        wg.HBox([Label('Atom'), self.atom])
        ])
        self.index = int(change.new)
        self.x.observe(self.edit_atom_by_index, 'value')
        self.y.observe(self.edit_atom_by_index, 'value')
        self.z.observe(self.edit_atom_by_index, 'value')
        self.atom.observe(self.edit_atom_by_index, 'value')
        self.new_atom = False

        with self.selected_atom:
            clear_output(True)
            display(self.atom_attributes)

    def edit_atom_by_index(self, change):
        """
        Method to edit atom by its index
        """
        try:
            if self.new_atom == False:
                self.df_coordinates['Atomic coordinate x (Å)'].iloc[self.index] = self.x.value
                self.df_coordinates['Atomic coordinate y (Å)'].iloc[self.index] = self.y.value
                self.df_coordinates['Atomic coordinate z (Å)'].iloc[self.index] = self.z.value
                self.df_coordinates['Atom'].iloc[self.index] = self.atom.value

            df_coordinates = self.df_coordinates
                
            unit_cell = self.structure_viewer.unit_cell
            cell_split = [x.split(' ') for x in unit_cell.split('\n')]
            df_cell = []
            for coo in cell_split:
                coo = np.array(coo, dtype=object)
                
                df_cell += [coo[coo!='']]
                
            df_cell = pd.DataFrame(df_cell, columns=['Lattice vector x (Å)','Lattice vector y (Å)','Lattice vector z (Å)']).dropna()
            df_cell = df_cell[['Lattice vector x (Å)','Lattice vector y (Å)','Lattice vector z (Å)']]
            df_cell.index = ['a','b','c']
            species = self.structure_viewer.chemical_species
            
            for specie in species:
                id_ = species[specie]['id']
                symbol = species[specie]['symbol']
                df_coordinates['Atom'] = df_coordinates['Atom'].replace(str(id_),symbol)
                
            self.sheet1 = ipysheet.from_dataframe(df_coordinates.dropna())
            self.sheet2 = ipysheet.from_dataframe(df_cell.dropna())
            create_structure_button = create_expanded_button('Update structure','success')
            create_structure_button.on_click(self.create_structure_manually)
            with self.structure_parameters:
                clear_output(True)
                logo = get_calculator_logo("./src/images/spin_logo.png", "png")
                display(logo)
                display(self.geometric_structure_box)
                display(self.structure_label_input)
                display(HBox([Label('Output job folder (Default "Calculations"): '), self.job_folder]))
                display(VBox([self.sheet1,self.sheet2, create_structure_button]))

            self.df_coordinates = df_coordinates
            
            self.get_atom_by_index.options = df_coordinates.dropna().index
            self.create_structure_manually('aux')

        except:
            pass

    def add_atom(self, change):
        """
        Method for adding new atom, default: H in 0, 0, 0, positions
        """
        if change.new == True:
            new_atom = {'Atomic coordinate x (Å)':0,
                        'Atomic coordinate y (Å)':0,
                        'Atomic coordinate z (Å)':0,
                        'Atom':'H'}

            self.df_coordinates = self.df_coordinates.append(new_atom,ignore_index = True).dropna()
            self.df_coordinates.index = np.arange(0,len(self.df_coordinates),1)
            self.new_atom = True
            self.edit_atom_by_index('value')

    def delete_atom(self, change):
        """
        Method for display and select delete atom widgets
        """
        if change.new == True:
            self.delete_atom_by_index = create_dropdown_input(self.df_coordinates.index, 0, '')
            self.delete_atom_button = create_expanded_button('delete','danger')   
            with self.selected_atom_to_drop:
                clear_output(True)
                display(self.delete_atom_by_index)   
            self.delete_atom_button.on_click(self.delete_atom_by_index_action)
            self.delete_atom_by_index.observe(self.display_delete_atom_properties, 'value')

        else:
            with self.selected_atom_to_drop:
                clear_output(True)

    def display_delete_atom_properties(self, change):
        """
        Method to display selected atom properties for delete
        """

        atom_info = self.df_coordinates.iloc[int(change.new)]
        self.x = create_float_text(atom_info['Atomic coordinate x (Å)'], '')
        self.y =  create_float_text(atom_info['Atomic coordinate y (Å)'], '')
        self.z = create_float_text(atom_info['Atomic coordinate z (Å)'], '')
        self.atom = create_text_input(atom_info['Atom'], '', width='70px')
        self.atom_attributes = wg.VBox([wg.HBox([Label('Atomic coordinate x (Å)'), self.x]),
        wg.HBox([Label('Atomic coordinate y (Å)'), self.y]),
        wg.HBox([Label('Atomic coordinate z (Å)'), self.z]),
        wg.HBox([Label('Atom'), self.atom])
        ])
        with self.selected_atom_to_drop:
            clear_output(True)
            display(self.atom_attributes)
            display(self.delete_atom_button)


    def delete_atom_by_index_action(self, change):
        """
        Method to delete atom after confirmation
        """
        index = self.delete_atom_by_index.value
        self.df_coordinates = self.df_coordinates.drop(index).reset_index()
        self.new_atom = True
        self.edit_atom_by_index('value')
        self.delete_atom_button.value = False 
        #self.create_structure_manually('aux')
        with self.selected_atom_to_drop:
            clear_output(True)


    def show_labels(self, change):
        """
        Method for show or hide structure labels
        """
        try:
            self.show_structure_labels = change.new 
        except:
            pass
        self.check_uploaded_structure_file()

    def change_viewer_style(self, change):
        """
        Method for change viewer style
        """
        self.viewer_style = change.new
        self.check_uploaded_structure_file()

    def change_visualizer(self,change):
        """
        Method for change visualizer
        """
        self.visualizer = change.new
        if self.visualizer == 'Py3Dmol':
            self.viewer_options.layout.display = 'block'
        else: 
            self.viewer_options.layout.display = 'none'
        self.check_uploaded_structure_file()
            
    def generate_molecules(self, aux):
        """
        Generate structure and display it on screen
        """
        self.predefined_structure = molecule(self.molecules_list.value)
        self.IS_MOLECULE = True
        self.update_content_with_predefined_structure()

    def make_bulk(self):
        """
        Method for create widgets for building crystal bulks
        """
        
        description_bulk = widgets.HTML(value="Crystal structure and lattice constant(s) will be guessed if not provided")
        
        #create bulk input widgets
        self.name = create_text_input(None, '')
        self.crystalstructure = create_dropdown_input(['sc', 'fcc', 'bcc', 'tetragonal',
                                                       'bct', 'hcp', 'rhombohedral', 
                                                       'orthorhombic', 'mlc', 'diamond', 
                                                       'zincblende', 'rocksalt', 'cesiumchloride', 
                                                       'fluorite','wurtzite','None'], 'None', '')
     
        #params
        self.a = create_text_input('None', '')
        self.b = create_text_input('None', '')
        self.c = create_text_input('None', '')
        self.alpha  = create_text_input('None', '')
        self.covera = create_text_input('None', '')
        self.u  = create_text_input('None', '')
        self.orthorhombic = create_boolean_input('Orthorhombic ', 'Construct orthorhombic unit cell instead of primitive cell which is the default')
        self.cubic = create_boolean_input('Cubic ','Construct cubic unit cell if possible.')
        generate_button = create_expanded_button('Generate structure','success')
        
        #display inputs in Vbox and Hbox for better visualization
        self.bulk_box = widgets.VBox([description_bulk,
                                     HBox([Label('Chemical symbol or symbols as in "MgO" or "NaCl": '),self.name]),
                                      HBox([Label('Crystal Structure: '),self.crystalstructure]),
                                      HBox([Label('Lattice constant a: '),self.a]), 
                                      HBox([Label('Lattice constant b. If only a and b is given, b will be interpreted as c instead: '),self.b]), 
                                      HBox([Label('Lattice constant c: '),self.c]), 
                                      HBox([Label('Angle in degrees for rhombohedral lattice: '),self.alpha]), 
                                      HBox([Label('c/a ratio used for hcp. Default is ideal ratio: sqrt(8/3): '),self.covera]),
                                      HBox([Label('Internal coordinate for Wurtzite structure: '),self.u]), 
                                      self.orthorhombic,
                                      self.cubic,
                                      generate_button],layout = {'width':'max-content'})

        #generate structure
        generate_button.on_click(self.generate_bulk) 

    def check_bulk_parameters(self):
        """
        Function to check inputs and convert them to corresponding type 
        """
        
        if self.a.value == 'None':
            self.avalue = None
        else:
            self.avalue = float(self.a.value)
            
        if self.b.value == 'None':
            self.bvalue = None
        else:
            self.bvalue = float(self.b.value)
            
        if self.c.value == 'None':
            self.cvalue = None
        else:
            self.cvalue = float(self.c.value)
            
        if self.alpha.value == 'None':
            self.alphavalue = None
        else:
            self.aplhavalue = float(self.alpha.value)
            
        if self.covera.value == 'None':
            self.coveravalue = None
        else:
            self.coveravalue = float(self.covera.value)
        
        if self.u.value == 'None':
            self.uvalue = None
        else:
            self.uvalue = float(self.u.value)
        
    def generate_bulk(self,aux):
        """
        check bulk params, generate structure and display it on screen
        """

        #Check bulk params
        self.check_bulk_parameters()
        self.IS_MOLECULE = False

        #Create Structure
        self.predefined_structure = bulk(self.name.value, crystalstructure=self.crystalstructure.value, a=self.avalue, 
                             b=self.bvalue, c=self.cvalue,alpha=self.alphavalue,covera= self.coveravalue, 
                             u=self.uvalue, orthorhombic=self.orthorhombic.value, cubic=self.cubic.value)

        self.update_content_with_predefined_structure()

    def make_nanotube_with_buckling(self):
        """
        Method for create widgets for building nanotubes with blucking
        """
        
        self.buckling_flag = True
        #Create nanotubes input widgets
        #n value
        self.n_w = create_int_text(1,"")
        #m value
        self.m_w = create_int_text(1,"")
        
        #bond
        self.bond_w = create_float_text(1,'')
        
        #buckling
        self.buckling_w = create_float_text(1,'')
        
        #Symbol
        self.symbol_w = create_text_input('P','')
        
        #genere button
        generate_button = create_expanded_button('Generate structure','success')
        
        self.nanotube_buckling_Box = widgets.VBox([
                                         HBox([Label('n value: '),self.n_w]),
                                         HBox([Label('m value: '),self.m_w]),
                                         HBox([Label('bond value: '),self.bond_w]),
                                         HBox([Label('Symbol: '),self.symbol_w]),
                                         HBox([Label('Buckling value: '),self.buckling_w]),
                                         generate_button],
                                         layout = {'width':'max-content'})

        generate_button.on_click(self.generate_nanotube_buckling)
    
    def make_nanotube_without_buckling(self):
        """
        Method for create widgets for building nanotubes without blucking
        """

        self.buckling_flag = False
        #Create nanotubes input widgets
        #n value
        self.n_wo = create_int_text(1,"")
        #m value
        self.m_wo = create_int_text(1,"")
        #length
        self.length_wo = create_int_text(1,'')
        #bond
        self.bond_wo = create_float_text(1,'')
        #Symbol
        self.symbol_wo = create_text_input('C','')
        #Vacuum
        self.vacuum_wo = create_text_input('None','')
        #genere button
        generate_button = create_expanded_button('Generate structure','success')
        
        self.nanotube_Box = widgets.VBox([
                                         HBox([Label('n value: '),self.n_wo]),
                                         HBox([Label('m value: '),self.m_wo]),
                                         HBox([Label('length value: '),self.length_wo]),
                                         HBox([Label('bond value: '),self.bond_wo]),
                                         HBox([Label('Symbol: '),self.symbol_wo]),
                                         HBox([Label('Vacumm value: '),self.vacuum_wo]),
                                         generate_button],
                                         layout = {'width':'max-content'})

        generate_button.on_click(self.generate_nanotube)
         
    def generate_nanotube_buckling(self, aux):
        """
        Method for generate nanotube with bucling with given attributes and display it on screen
        """
        self.IS_MOLECULE = False
        bond = self.bond_w.value
        buckling = self.buckling_w.value
        symbol = self.symbol_w.value
        a2 = ((3/2.)*bond, -np.sqrt(3)/2. * bond )
        a1 = ((3/2.)*bond, np.sqrt(3)/2. * bond )

        #atomos iniciales 

        atom1 = (0.0,0.0,5.0)
        atom2 = (bond,0.0,5 + buckling)
        atomsR, unitCell = Nanotube(bond, a1, a2, atom1, atom2,self.n_w.value,self.m_w.value,False)
        self.predefined_structure = Atoms(symbol + str(np.shape(atomsR)[0]),
                                            positions = atomsR)
        self.predefined_structure.cell = [[30,0,0],[0,30,0],[0,0,20]]
        self.predefined_structure.center()
        self.update_content_with_predefined_structure()
        
            
    def generate_nanotube(self,aux):
        """
        Method for generate nanotube withot bucling with given attributes and display it on screen
        """
        self.IS_MOLECULE = False
        buckling = self.buckling_flag

        if self.vacuum_wo.value == 'None':
            self.vacuumvalue = None
        else:
            self.vacuumvalue = float(self.vacuum.value)

        self.predefined_structure = ase_nanotube(n = self.n_wo.value, m = self.m_wo.value,
                                        length=self.length_wo.value, bond=self.bond_wo.value,
                                        symbol=self.symbol_wo.value, vacuum = self.vacuumvalue)
        self.predefined_structure.cell = [[30,0,0],[0,30,0],[0,0,20]]
        self.predefined_structure.center()
        self.update_content_with_predefined_structure()

            
    def make_graphene_nanoribons(self):
        """
        Method for create widgets for building graphene nanoribons
        """
        description_g_n = widgets.HTML(value="Creates a graphene nanoribbon in the x-z plane, with the nanoribbon running along the z axis.")
        
        self.n_g_n = create_int_text(1,'')
        self.m_g_n = create_int_text(1,'')
        self.type_g_n = create_dropdown_input(['zigzag','armchair'],'zigzag','')
        self.saturated_g_n = create_boolean_input('Saturated','If true, hydrogen atoms are placed along the edge.')
        self.C_H_g_n = create_float_text(1.09,'')
        self.C_C_g_n = create_float_text(1.42,'')
        self.vacuum_g_n = create_text_input('None','')
        self.magnetic_g_n = create_boolean_input('magnetic','Make the edges magnetic.')
        self.initial_mag_g_n = create_float_text(1.12, '')
        self.sheet_g_n = create_boolean_input('sheet','If true, make an infinite sheet instead of a ribbon ')
        #genere button
        generate_button = create_expanded_button('Generate structure','success')
        
        self.g_n_box = widgets.VBox([description_g_n,
                                    HBox([Label('n value: The width of the nanoribbon. For armchair nanoribbons, this n may be half-integer to repeat by half a cell.'),self.n_g_n]),
                                    HBox([Label('m value: The length of the nanoribbon.'),self.m_g_n]),
                                    HBox([Label('type value: The orientation of the ribbon '),self.type_g_n]),
                                    HBox([self.saturated_g_n]),
                                    HBox([Label('Carbon-hydrogen bond length (Angstrom): '),self.C_H_g_n]),
                                    HBox([Label('Carbon-carbon bond length (Angstrom): '),self.C_C_g_n]),
                                    HBox([Label('Vacuum: Amount of vacuum added to non-periodic directions '),self.vacuum_g_n]),
                                    HBox([self.magnetic_g_n]),
                                    HBox([Label('Magnitude of magnetic moment if magnetic.'),self.initial_mag_g_n]),
                                    HBox([self.sheet_g_n]), generate_button],
                                   layout = {'width':'max-content'})

        generate_button.on_click(self.generate_graphene_nanoribbon)
        
    def generate_graphene_nanoribbon(self,aux):
        """
        Method for generate graphene nanoribbon with given attributes and display it on screen
        """
        self.IS_MOLECULE = False
        if self.vacuum_g_n.value == 'None':
            self.vacuumvalue = None
        else:
            self.vacuumvalue = float(self.vacuum_g_n.value)
            
        self.predefined_structure = graphene_nanoribbon(n=self.n_g_n.value, m=self.m_g_n.value, type=self.type_g_n.value,
                                             saturated=self.saturated_g_n.value, C_H=self.C_H_g_n.value,
                                             C_C=self.C_C_g_n.value, vacuum=self.vacuumvalue,
                                             magnetic=self.magnetic_g_n.value, initial_mag=self.initial_mag_g_n.value,
                                         sheet=self.sheet_g_n.value)
        self.predefined_structure.center()
        self.predefined_structure.cell = [[15,0,0],[0,15,0],[0,0,15]]
        self.predefined_structure.center()
        self.update_content_with_predefined_structure()
        
          

    def make_structure_manually(self, aux):
        """
        Get the initial schema for the manual creation of structures.
        According to number of atoms sheet with atomic coordinates (sheet1)
        and sheet with unit cell (sheet2)
        """
        
        num_rows = self.number_of_atoms_input.value
        self.sheet1 = sheet(rows=num_rows,columns=4,
                     column_headers = ['Atomic coordinate x (Å)','Atomic coordinate y (Å)','Atomic coordinate z (Å)', 'Atom'])


        column(0,[0.0]*num_rows)
        column(1,[0.0]*num_rows)
        column(2,[0.0]*num_rows)
        column(3,['C']*num_rows)
        
        self.sheet2 = sheet(rows=3,columns=3,row_headers = ['a','b','c'], column_headers = ['Lattice vector x (Å)','Lattice vector y (Å)','Lattice vector z (Å)'])
        column(0,[0.0]*3)
        column(1,[0.0]*3)
        column(2,[0.0]*3)
        
        create_structure_button = create_expanded_button('Get structure','success')
        create_structure_button.on_click(self.create_structure_manually)
        
        with self.structure_parameters:
            clear_output(True)
            logo = get_calculator_logo("./src/images/spin_logo.png", "png")
            display(logo)
            display(self.structure_label_input)
            display(HBox([Label('Output job folder (Default "Calculations"): '), self.job_folder]))

            display(self.geometric_structure_box)
            display(VBox([self.sheet1,self.sheet2, create_structure_button]))
    
    def create_structure_manually(self, aux):
        """
        Method for make the structure that user enter manually when "Get Structure" button is clicked
        """
        update_content = False
        try:
            uploaded_file = self.structure_input.value
            filename = next(iter(uploaded_file))
            update_content = True
        except:
            pass
            
        datos = ipysheet.to_dataframe(self.sheet1)
        unit_cell = ipysheet.to_dataframe(self.sheet2)
        atoms = Atoms(symbols=datos.Atom.to_list(), positions=datos[['Atomic coordinate x (Å)','Atomic coordinate y (Å)','Atomic coordinate z (Å)']].values, pbc = False)
        atoms.set_positions(datos[['Atomic coordinate x (Å)','Atomic coordinate y (Å)','Atomic coordinate z (Å)']].values, apply_constraint = False)
        atoms.set_cell(cell = unit_cell[['Lattice vector x (Å)','Lattice vector y (Å)','Lattice vector z (Å)']].values, scale_atoms=False, apply_constraint=False)    

        if update_content == False:

            write('temp.xyz', atoms)
            content = ''
            with open('./temp.xyz', 'r') as file:
                content = file.read()
            
            content = self.convert_any_format_to_cif([{'temp.xyz':{'name':'Prueba.xyz',
                                    'type':'type_',
                                    'size':'size',
                                    'lastModified':'lastModified',
                                    'content':bytes(content, 'utf-8')}}, 'temp.xyz','extxyz'])
            filename = 'temp.extxyz'

        else:
            write('temp.cif', atoms)
            content = ''
            
            with open('temp.cif', 'r') as file:
                content = file.read()
            filename = 'structure.cif'
    
        self.structure_input = Simulated_structure(content, filename = filename)
        uploaded_file = self.structure_input.value
        filename = next(iter(uploaded_file))
        self.structure_input.value[filename]['content'] = str.encode(content)
        self.structure_from_file = False
        self.show_labels({'new':self.show_structure_labels})
        with self.out:
            clear_output(True)
            self.structure_viewer = Viewer('structure',content,
            [self.show_structure_labels,self.viewer_style,self.visualizer ],format_= filename.split(".")[-1])
            if self.visualizer == 'Py3Dmol':
                x = HTML(self.structure_viewer.build_viewer())
            else:
                x = self.structure_viewer.build_viewer()
            display(x)




    def convert_any_format_to_cif(self, args):
        """
        Method to convert any uploaded file to cif file.
        """
        
        uploaded_file = args[0]
        filename = args[1]
        format_ = args[2]

        if format_ != 'cif' and format_ != 'extxyz':
            
            with open(filename, 'w') as f:
                f.write("".join(map(chr, uploaded_file[filename]['content'])))
            f.close()

            if format_ == 'fdf':
                
                fdf = read_fdf(filename)
                dict_species_label = {}
                for i in fdf['chemicalspecieslabel']:
                    if len(i[2].split('-')) > 1:
                        dict_species_label[i[0]] = i[2].split('-')[0]
                    else:
                        dict_species_label[i[0]] = i[2].split('.')[0]

                atoms = []
                symbols = []
  
                for i in fdf['atomiccoordinatesandatomicspecies']:

                    atoms += [i[:3]]
                    symbols += [dict_species_label[i[3]]]

                try:
                    atoms = Atoms(symbols=symbols, positions=atoms,
                                    cell = fdf['latticevectors'])
                except:
                   atoms = Atoms(symbols=symbols, positions=atoms)


                write('temp.cif', atoms, wrap=False)
                content = ''
                
                with open('temp.cif', 'r') as file:
                    content = file.read()

                self.structure_input = Simulated_structure(content, f"{filename}.cif")
                uploaded_file = self.structure_input.value
                filename = next(iter(uploaded_file)).split(".cif")[0]
                self.structure_input.value[f"{filename}.cif"]['content'] = str.encode(content)
                os.system(f"chmod 777 temp.cif")
                os.remove('temp.cif')

            else:
                for molecule in pybel.readfile(format_,filename):
                    molecule
                content = molecule.write('cif')

            os.system(f"chmod 777 {filename}")
            os.remove(filename)
        else:
            content = uploaded_file[filename]['content'].decode("utf-8") 
        return content
            
        
    def check_uploaded_structure_file(self,*args):
        """
        Method to checking that structure file is correctly uploaded.
        If correct, displays the structure atomic coordinates and unit cell. 
        """
        self.structure_input._counter = 1
        if len(args) > 0:
            self.structure_from_file = True
   
            uploaded_file = args[0].new
            filename = next(iter(uploaded_file))
            
            try: format_ = filename.split('.')[1]
            except: format_ = filename.split('.')[0]
        else:
            uploaded_file = self.structure_input.value
            filename = next(iter(uploaded_file))
            format_ = 'cif'

        try:
            if filename.split('.')[-1] == 'extxyz':
                format_= 'extxyz'
            if self.structure_from_file == True:
                text='Structure file has been loaded successfully! You can see it in "Structure" tab.'
                self.valid_structure_input.description = text
                self.valid_structure_input.value = True
                self.valid_structure_input.layout.display = 'block'
            content = self.convert_any_format_to_cif([uploaded_file, filename,format_])
            
            #Display the uploaded structure visualization
            with self.out:
                clear_output(True)
                
                self.structure_viewer = Viewer(filename,content, [self.show_structure_labels, self.viewer_style,
                self.visualizer ], format_=format_)
                if self.visualizer == 'Py3Dmol':
                    x = HTML(self.structure_viewer.build_viewer())
                else:
                    x = self.structure_viewer.build_viewer()
                display(x)
               
               
            #Get and displays the atomic coordinates and unit cell
            create_structure_button = create_expanded_button('Update structure','success')
            create_structure_button.on_click(self.create_structure_manually)
            coordinates = self.structure_viewer.coordinates
            coordinates_split = [x.split(' ') for x in coordinates.split('\n')]
            df_coordinates = []
            for coo in coordinates_split:
                coo = np.array(coo, dtype=object)
                df_coordinates += [coo[coo!='']]
            
            df_coordinates = pd.DataFrame(df_coordinates, columns=['Atomic coordinate x (Å)','Atomic coordinate y (Å)','Atomic coordinate z (Å)', 'Atom'])
            
            unit_cell = self.structure_viewer.unit_cell
            cell_split = [x.split(' ') for x in unit_cell.split('\n')]
            df_cell = []
            for coo in cell_split:
                coo = np.array(coo, dtype=object)
                
                df_cell += [coo[coo!='']]
                
            df_cell = pd.DataFrame(df_cell, columns=['Lattice vector x (Å)','Lattice vector y (Å)','Lattice vector z (Å)']).dropna()
            df_cell = df_cell[['Lattice vector x (Å)','Lattice vector y (Å)','Lattice vector z (Å)']]
            df_cell.index = ['a','b','c']
            species = self.structure_viewer.chemical_species
            
            for specie in species:
                id_ = species[specie]['id']
                symbol = species[specie]['symbol']
                df_coordinates['Atom'] = df_coordinates['Atom'].replace(str(id_),symbol)
                
            self.sheet1 = ipysheet.from_dataframe(df_coordinates.dropna())
            self.sheet2 = ipysheet.from_dataframe(df_cell.dropna())
            
            with self.structure_parameters:
                clear_output(True)
                logo = get_calculator_logo("./src/images/spin_logo.png", "png")
                display(logo)
                display(self.geometric_structure_box)
                display(self.structure_label_input)
                display(HBox([Label('Output job folder (Default "Calculations"): '), self.job_folder]))

                display(VBox([self.sheet1,self.sheet2, create_structure_button]))

            self.df_coordinates = df_coordinates
            
            self.get_atom_by_index.options = df_coordinates.dropna().index
            if self.structure_from_file == True:
                self.create_structure_manually('aux')
            #self.structure_from_file = True
        except:
        #else:
            text=f'Something wrong! {format_} is not a recognised format, please try again!'
            self.valid_structure_input.description = text
            self.valid_structure_input.value = False
            self.valid_structure_input.layout.display = 'block'

class Calculations(Spin):
    def __init__(self):
        """
        display calculations blocks
        """
        self.create_run_tabs()
        logo = get_calculator_logo("./src/images/spin_logo.png", "png")
        display(logo)
        display(self.supported_runs)

    def create_run_tabs(self):
        """
        Create run tab for each calculator.
        """

        self.siesta_run()
        
        self.supported_runs = wg.Accordion(children = [self.run_siesta_box])
        self.supported_runs.set_title(0, 'SIESTA run')

class Postprocessing(Spin):

    def __init__(self):
        """
        display postprocessing blocks
        """
        self.create_postprocessing_tabs()
        logo = get_calculator_logo("./src/images/spin_logo.png", "png")
        display(logo)
        display(self.supported_postprocessing)

    def create_postprocessing_tabs(self):
        """
        Create postprocessing tab for each calculator.
        """

        self.siesta_postprocessing()
        
        self.supported_postprocessing = wg.Accordion(children = [self.postprocessing_siesta_box])
        self.supported_postprocessing.set_title(0, 'SIESTA post-processing') 

class Simulated_structure(object):
    """
    Simulate an uploaded file content when building predefined structure and manual structure
    """

    def __init__(self, content, filename='structure.cif'):
        self.content = content
        self.filename = filename
        self.get_simulated_structure_input()
    
    def get_simulated_structure_input(self):

        filename = self.filename
        type_ = ''
        size = ''
        lastModified = ''
        content = self.content

        self.value = {filename:{'metadata':{'name':filename,
                                'type':type_,
                                'size':size,
                                'lastModified':lastModified,
                                'content':content}}}


from multiprocessing.managers import ValueProxy
from operator import le
import ipywidgets as wg
import sys
from matplotlib.pyplot import savefig
sys.path.append('../../')
sys.path.append('../../../')
sys.path.append('../../../images/')
import numpy as np
import os
from ipywidgets import  Layout, Label
from IPython.display import display, HTML, IFrame, clear_output,Javascript
from utils import *
from create_widgets import *
import py3Dmol
from create_fdf_file import *
from siesta_check_outputs import *
from viewer import *
from siesta import *
import pandas as pd
from utils import *
import base64
import time
from io import StringIO
import random
from ipyfilechooser import FileChooser
import json

class Siesta(object):

    """
    Provide Siesta elements for main app
    """

    def select_SIESTA_parameters(self):
        
        """
        main method for build siesta parameters inputs
        """
        self.structure_input = wg.FileUpload(accept='*.cif', multiple=False)
        self.structure_input.observe(self.read_siesta_file, 'value')
        self.structure_input.observe(self.check_uploaded_structure_file,'value') 
        self.valid_structure_input = wg.Valid(value=True,description='Valid!')
        
        siesta_logo = get_calculator_logo("./src/images/siesta.png", "png")

        title_SIESTA = wg.HTML(
            value="SIESTA parameters")
        
        #define box layout / box style
        box_layout = Layout(display='flex',
                    flex_flow='column',
                    align_items='stretch',
                    align_content='stretch',
                    flex='flex-grow',
                    border='doted',
                    width='100%')
        
        #Get widgets for siesta calcule control
        self.calcule_control()
        

        #get siesta examples
        self.bandLinesExample, self.bandPointsExample, self.exampleProjectedDensityOfStates,self.examplePDOSkgrid_Monkhorst_Pack, self.exampleOpticalMesh,self.exampleOpticalVector =  get_siesta_examples()

        self.relax_structure = create_checkbox_input('')
        #Get widgets for siesta relax control
        self.system_relax_control()

        #bands parameters
        self.bands = create_checkbox_input('')
        self.create_siesta_bands_parameters()
        
        #dos parameters
        self.dos = create_checkbox_input('')
        self.create_siesta_dos_parameters()
        
        #optical parameters
        self.opt = create_checkbox_input('')
        self.create_siesta_opt_parameters()

        #aditional params
        self.create_aditional_parameters_box()
        

                #put siesta calcule control widgets in a box
        dft_calculator_control_box = wg.VBox([HBox([Label('Basis set: '),self.basis_set_input]),
                                       HBox([Label('Mesh cutoff: '),
                                             self.mesh_cutoff_input,
                                             self.mesh_cutoff_units]),
                                       HBox([Label('XC Functional and Author: '),self.xc_input, self.author]),
                                       HBox([Label('Energy shift: '),self.energy_shift_input,self.energy_shift_units]),
                                       HBox([Label('Max. SCF Iterations: '),self.max_SCF_iterations_input]),
                                       HBox([Label('Mixing weight: '), self.mixing_weigth_input]),
                                       HBox([Label('Number pulay: '), self.number_pulay_input]),
                                       HBox([Label('Spin polarization: '), self.spin_input]),
                                       HBox([Label('kgrid_Monkhorst_Pack: '),self.kpts_options_input,self.kpts_input]),
                                             HBox([self.relax_structure, Label('Relax structure')]),
                                            self.system_relax_control_box,
                                         HBox([self.bands, Label('Calculate Band structure')]),
                                       self.siesta_bands_parameters,
                                         HBox([self.dos, Label('Calculate Total and Partial (projected) Density of States')]),
                                       self.siesta_dos_parameters,
                                         HBox([self.opt, Label('Calculate Optical')]),
                                       self.siesta_opt_parameters],
                                              layout=box_layout)

        self.siesta_blocks = wg.Accordion(children = [dft_calculator_control_box, self.aditional_parameters_box])
        
        
        self.siesta_blocks.set_title(0, 'DFT calcule control')
        self.siesta_blocks.set_title(1, 'Aditional parameters')

        #create button for the creation of fdf file for run
        create_fdf_button = create_expanded_button('Make fdf','primary')
        create_fdf_button.on_click(self.make_fdf)
        self.create_fdf_message = widgets.Output()
   
        #Create button for the download od the fdf input file
        self.download_file = create_expanded_button('¿Do you want to download input file?','primary')
        self.download_file.on_click(self.download_fdf)
        self.download_file_ = wg.HTML(value='None')
        self.download_file_.layout.display = 'none'

        self.create_siesta_box = wg.VBox([self.label_input,
                                        HBox([Label('Output job folder (Default "Calculations"): '),self.job_folder]),
                                       self.siesta_blocks,
                                       create_fdf_button,self.create_fdf_message,
                                       HBox([self.download_file, self.download_file_]),
                                       ],
                                       layout = {'width':'max-content'})
        
        self.siesta_blocks_options = wg.Accordion(children = [self.structure_input, self.create_siesta_box])
        self.siesta_blocks_options.set_title(0, 'Upload input file (.fdf)')
        self.siesta_blocks_options.set_title(1, 'Create input file')

        siesta_manual = wg.Button(description='   Open SIESTA manual',
                        disabled=False,
                        button_style='', # 'success', 'info', 'warning', 'danger' or ''
                        tooltip='Open siesta manual in new browser tab',
                        icon='info', # (FontAwesome names without the `fa-` prefix)
                        layout=Layout(height='50px', width='200px'))
        siesta_manual.on_click(self.window_open)

        self.siesta_box = wg.VBox([wg.HBox([siesta_logo, siesta_manual]), self.siesta_blocks_options],layout = {'width':'max-content'})

    def window_open(self, aux):
        url = "https://departments.icmab.es/leem/SIESTA_MATERIAL/Docs/Manuals/siesta-4.1-b4.pdf"
        display(Javascript('window.open("{url}");'.format(url=url)))

    def read_siesta_file(self, *args):
            
        if len(args) > 0:
            uploaded_file = args[0].new
            filename = next(iter(uploaded_file))
            try: format_ = filename.split('.')[1]
            except: format_ = filename.split('.')[0]
        else:
            uploaded_file = self.structure_input.value
            filename = next(iter(uploaded_file))
            format_ = 'fdf'
            

        if format_ == 'fdf':
            with open(filename, 'w') as f:
                f.write("".join(map(chr, uploaded_file[filename]['content'])))
            f.close()

            fdf = read_fdf(filename) 
            if 'paobasissize' in fdf.keys():
                self.basis_set_input.value = fdf['paobasissize'][0]
            if  'paoenergyshift' in fdf.keys():
                self.energy_shift_input.value = fdf['paoenergyshift'][0]
                self.energy_shift_units.value = fdf['paoenergyshift'][1]
            if 'meshcutoff' in fdf.keys():
                self.mesh_cutoff_input.value = fdf['meshcutoff'][0]
                self.mesh_cutoff_units.value = fdf['meshcutoff'][1]
            if 'xcfunctional' in fdf.keys():
                self.xc_input.value = fdf['xcfunctional'][0]
            if 'xcauthors' in fdf.keys():
                self.author.value = fdf['xcauthors'][0]
            if 'maxscfiterations' in fdf.keys():
                self.max_SCF_iterations_input.value = fdf['maxscfiterations'][0]
            if 'dmmixingweight' in fdf.keys():
                self.mixing_weigth_input.value = fdf['dmmixingweight'][0]
            if 'dmnumberpulay' in fdf.keys():
                self.number_pulay_input.value = fdf['dmnumberpulay'][0]
            if 'spin' in fdf.keys():
                self.spin_input.value = fdf['spin'][0]
            if 'kgridmonkhorstpack' in fdf.keys():
                self.kpts_options_input.value = 'Enter kgrid_Monkhorst_Pack'
                self.check_siesta_kpoints_input()
                self.k_1.value = fdf['kgridmonkhorstpack'][0][0]
                self.k_2.value = fdf['kgridmonkhorstpack'][0][1]
                self.k_3.value = fdf['kgridmonkhorstpack'][0][2]
                self.k_4.value = fdf['kgridmonkhorstpack'][1][0]
                self.k_5.value = fdf['kgridmonkhorstpack'][1][1]
                self.k_6.value = fdf['kgridmonkhorstpack'][1][2]
                self.k_7.value = fdf['kgridmonkhorstpack'][2][0]
                self.k_8.value = fdf['kgridmonkhorstpack'][2][1]
                self.k_9.value = fdf['kgridmonkhorstpack'][2][2]
            


        os.system(f"chmod 777 {filename}")
        os.remove(filename)

    def download_fdf(self, aux):
        """
        Method for download the fdf file
        """
        html_fdf = self.make_fdf(None,download=True)
        self.download_file_.value=html_fdf
        self.download_file_.layout.display = 'block'
        time.sleep(20)
        self.download_file_.layout.display = 'none'
      
        
        
    def calcule_control(self):
        """
        Create input for siesta calcule control box for fdf creation
        """
        
        self.label_input = create_text_input(self.structure_label,'Label: ', width='350px')
        self.label_input.observe(self.update_label_input, 'value')
        self.job_folder = FileChooser('./calculations/')#,layout=Layout(height='200px', width='500px'))
        self.job_folder.default_filename = self.label_input.value
        self.mesh_cutoff_input = create_int_text(200,'')
        self.mesh_cutoff_units = create_dropdown_input(['eV','Ang','nm','Bohr','Hartree','kJ','kcal','mol','Ry','meV'],
                                                  'Ry','')
        self.converge_mesh_cutoff_min = create_int_text(50, '')
        self.converge_mesh_cutoff_max = create_int_text(350, '')
        self.converge_mesh_cutoff_interval = create_int_text(50, '')
        self.converge_mesh_cutoff_input = VBox([HBox([Label('min: '), self.converge_mesh_cutoff_min]),
                                          HBox([Label('max: '), self.converge_mesh_cutoff_max]),
                                          HBox([Label('interval: '), self.converge_mesh_cutoff_interval])])

        self.mesh_cutoff_options_input = wg.RadioButtons(
            options=['Enter Mesh cutoff'],
            value='Enter Mesh cutoff',
            description='',
            disabled=False
        )
        self.converge_mesh_cutoff_input.layout.display = "none"
        
        self.xc_input = create_dropdown_input(['GGA','LDA','VDW'],'GGA','')
        self.author = create_dropdown_input(['PBE','PW91','revPBE','RPBE','WC','AM05','PBEsol',
                                   'PBEJsJrLO','PBEGcGxLO','PBEGcGxHEG','BLYP'],'PBE','')
        self.energy_shift_input = create_float_text(100, '')
        self.energy_shift_units = create_dropdown_input(['eV','Ang','nm','Bohr','Hartree','kJ','kcal','mol','Ry','meV'],
                                                 'meV','')
        self.max_SCF_iterations_input = create_int_text(150, '', 50)
        self.mixing_weigth_input = create_float_text(0.1, '')
        self.number_pulay_input = create_int_text(5, '')
        self.basis_set_input = create_dropdown_input(['SZ','DZ','SZP','DZP','TZ','TZP'],'DZP','')
        
        
        self.spin_input = create_dropdown_input(['non-polarized','collinear','non-collinear','spin-orbit'],
                                                'non-polarized','','120px')

        self.k_1 = create_int_text(1, '')
        self.k_2 = create_int_text(0, '') 
        self.k_3 = create_int_text(0, '')
        self.k_4 = create_int_text(0, '')
        self.k_5 = create_int_text(1, '')
        self.k_6 = create_int_text(0, '')
        self.k_7 = create_int_text(0, '')
        self.k_8 = create_int_text(0, '')
        self.k_9 = create_int_text(2, '')

        uno = wg.HBox([self.k_1, self.k_2, self.k_3], layout = {'width':'100%'})
        dos = wg.HBox([self.k_4, self.k_5, self.k_6], layout = {'width':'100%'})
        tres = wg.HBox([self.k_7, self.k_8, self.k_9,], layout = {'width':'100%'})

        self.kpts_input = VBox([uno, dos, tres])
        self.converge_kpts_min = create_int_text(1, '')
        self.converge_kpts_max = create_int_text(12, '')
        self.converge_kpts_interval = create_int_text(1, '')
        self.converge_kpts_input = VBox([HBox([Label('min: '), self.converge_kpts_min]),
                                          HBox([Label('max: '), self.converge_kpts_max]),
                                          HBox([Label('interval: '), self.converge_kpts_interval])])
        self.kpts_input.layout.display = "none"
        self.converge_kpts_input.layout.display = "none"
        
        self.kpts_options_input = wg.RadioButtons(
            options=['Gamma', 'Enter kgrid_Monkhorst_Pack'],
            value='Gamma',
            description='',
            disabled=False)
        
        self.xc_input.observe(self.update_siesta_author,'value')
        self.kpts_options_input.observe(self.check_siesta_kpoints_input,'value')
        self.mesh_cutoff_options_input.observe(self.check_siesta_mesh_cutoff_input, 'value')

    def update_label_input(self,change):
        self.job_folder.default_filename = change.new

    def system_relax_control(self):
        """
        Create input for siesta calcule relax box for fdf creation
        """
        
        self.cg_steps_input =  create_int_text(200, '', 50)

        self.variable_cell_input = create_dropdown_input(['.False.','.True.'],
                                                '.True.','')

        self.typeRun_input = create_dropdown_input(['CG','Broyden','FIRE','Verlet',
                                                'Nose','ParrinelloRahman','NoseParrinelloRahman',
                                                'Anneal','FC','Master','Lua'],
                                                'FIRE','')

        self.max_force_tol_input = create_float_text(0.04, '')

        
        #put siesta relax control widgets in a box
        self.system_relax_control_box = wg.VBox([HBox([Label('Num CG steps: '), self.cg_steps_input ]),
                                HBox([Label('Variable cell: '),self.variable_cell_input]),
                                HBox([Label('Type of Run: '),self.typeRun_input]),
                                HBox([Label('Max force tol [eV/Ang]: '),self.max_force_tol_input])
                                ],
                                        layout = {'width':'max-content'}) 

        self.system_relax_control_box.layout.display = "none"
        self.relax_structure.observe(self.get_siesta_relax_parameters,'value')

    def geometric_structure(self):
        """
        Create input for geometric structure box for structure creation
        """
        
        self.structure_input = wg.FileUpload(accept='*.cif', multiple=False)
        self.structure_input.observe(self.check_uploaded_structure_file,'value')
        self.valid_structure_input = wg.Valid(value=True,description='Valid!')
        
        self.create_structure_manually_button = create_expanded_button('Create structure manually','primary')
        self.create_structure_manually_button.on_click(self.make_structure_manually)
        self.number_of_atoms_input = create_int_text(1,'')
        
    def make_fdf(self, aux, download = False):
        """
        Create fdf input according to user parameters input
        """   
        siesta_label = self.label_input.value
        try:
            number_of_atoms = self.structure_viewer.number_of_atoms
            number_of_species = self.structure_viewer.number_of_species
            chemical_species = self.structure_viewer.chemical_species
            chemical_species_label = ''
        except:
            return 'Please build or load a structure.'
        
    
        basis_set = self.basis_set_input.value
        energy_shift = self.energy_shift_input.value
        energy_shift_unit = self.energy_shift_units.value
        xc_functional = self.xc_input.value
        xc_author = self.author.value
        
        self.pseudopotentials = []
        for specie in chemical_species:
       
            chemical_species_label += str(chemical_species[specie]['id']) + '  '   
            chemical_species_label += str(chemical_species[specie]['atomic_number']) + '  '
            chemical_species_label += str(chemical_species[specie]['symbol'])+f'-{xc_functional.lower()}\n'
            self.pseudopotentials += [str(chemical_species[specie]['symbol'])+f'-{xc_functional.lower()}.psf']
        
            
        max_scf_iterations = self.max_SCF_iterations_input.value
        mesh_cutoff = self.mesh_cutoff_input.value
        mesh_cutoff_units = self.mesh_cutoff_units.value
        mixing_weight =  self.mixing_weigth_input.value
        number_pulay = self.number_pulay_input.value
        spin = self.spin_input.value 
        k1 = self.k_1.value
        k2 = self.k_2.value
        k3 = self.k_3.value
        k4 = self.k_4.value
        k5 = self.k_5.value
        k6 = self.k_6.value
        k7 = self.k_7.value
        k8 = self.k_8.value
        k9 = self.k_9.value

        cg_steps = self.cg_steps_input.value
        variable_cell = self.variable_cell_input.value
        type_run = self.typeRun_input.value
        max_force_tol = self.max_force_tol_input.value

        

        unit_cell = self.structure_viewer.unit_cell
        check_unit_cell = unit_cell.replace('  ',' ').split()
        check_unit_cell = [float(i) for i in check_unit_cell]
        if all(np.array(check_unit_cell) == 0):
            unit_cell = False
        coordinates = self.structure_viewer.coordinates
        
        args = [siesta_label,
                number_of_atoms,
                number_of_species,
                basis_set,
                energy_shift,
                energy_shift_unit,
                xc_functional,
                xc_author,
                max_scf_iterations,
                mesh_cutoff,
                mesh_cutoff_units,
                mixing_weight,
                number_pulay,
                k1,
                k2,
                k3,
                k4,
                k5,
                k6,
                k7,
                k8,
                k9,
                cg_steps,
                variable_cell,
                unit_cell,
                coordinates,
                spin,chemical_species_label, type_run,
                max_force_tol]

        
        if download is False:
            # if download is False, then we create a directory and save file on it
            try:
                self.job_folder.value.split('/')
                self.root_directory = self.job_folder.value
                self.structure_directory = f'{self.job_folder.value}/structure/'
            except:
                self.root_directory = f'./calculations/{siesta_label}'
                self.structure_directory = f'./calculations/{siesta_label}/structure/'

            if not os.path.exists(self.root_directory):
                os.makedirs(self.root_directory)
    
            if not os.path.exists(self.structure_directory):
                os.makedirs(self.structure_directory)
                
            get_fdf_file(args, self.structure_directory + f'{siesta_label}.fdf',IS_MOLECULE = self.IS_MOLECULE)
           
            if self.aditional_parameters.value != '':
                os.system(f"echo '{self.aditional_parameters.value}' >> {self.structure_directory}{siesta_label}.fdf")

            #Save configurations params
            if self.relax_structure.value == True:
                self.process_siesta_relax_input()
                relax_inputs = self.siesta_relax_inputs
            else:
                relax_inputs = None
              

            if self.bands.value == True: 
                self.process_siesta_bands_input()
                bands_inputs = self.siesta_bands_inputs
            else:
                bands_inputs = None

            
            if self.dos.value == True:
                self.process_siesta_dos_input()
                dos_inputs = self.siesta_dos_inputs
            else:
                dos_inputs = None
            
            if self.opt.value == True:
                self.process_siesta_opt_input()
                opt_inputs = self.siesta_opt_inputs
            else:
                opt_inputs = None
            

            self.config = {'siesta_label':siesta_label,
                            'relax_structure':self.relax_structure.value,
                            'bands':self.bands.value,
                            'dos':self.dos.value,
                            'optical':self.opt.value,
                            'relax_structure_inputs':relax_inputs,
                            'bands_inputs':bands_inputs,
                            'dos_inputs':dos_inputs,
                            'opt_inputs':opt_inputs,
                            'pseudopotentials':self.pseudopotentials}

            with open(f'{self.root_directory}/config.json', 'w') as configfile:
                json.dump(self.config, configfile)

            text = f"{siesta_label} job was successfully created in {self.root_directory}"
            self.IS_MOLECULE = False

            with self.create_fdf_message:
                clear_output(wait=True)
                display(wg.Valid(value=True,description=''),wg.HTML(value=text))
            
        else:
            # if download is True, then the user has the possiblity to download the fdf in election path
            filename = f'{siesta_label}.fdf'
            fdf_text = get_fdf_file(args,filename,download=True, IS_MOLECULE = self.IS_MOLECULE)
            if self.aditional_parameters.value != '':
                fdf_text += self.aditional_parameters.value
            b64 = base64.b64encode(fdf_text.encode())
            payload = b64.decode()

            #Download button
            html_buttons = f'''<html>
            <head>
            <meta name="viewport" content="width=device-width, initial-scale=1">
            </head>
            <body>
            <a download="{filename}" href="data:text/csv;base64,{payload}" download>
            <button class="p-Widget jupyter-widgets jupyter-button widget-button mod-warning">Download File</button>
            </a>
            </body>
            </html>
            '''
            html_button = html_buttons.format(payload=payload,filename=filename)
            return html_button

    def list_available_calculations(self, change):
        try:
            folder = self.job_run_folder.value.split('/ ')[0]
            if not os.path.exists(folder):
                os.makedirs(folder)
            list_calculations = os.listdir(folder)

            if len(list_calculations) == 0:
                self.available_calculations.options = ['None']
                self.available_calculations.value = 'None'
            else:
                self.available_calculations.options = list_calculations
                self.available_calculations.value = [list_calculations[0]]
        except:pass
 

    def siesta_run(self):
        """
        Create siesta run calculations block
        """
        
        siesta_logo = get_calculator_logo("./src/images/siesta.png", "png")
    
        #job parameters
        self.n_processors = create_int_text(1,'')
        self.siesta_command = create_text_input('siesta','')
        self.run_siesta_button = create_expanded_button('Run','Danger')



        #list already create works in calculations folder
        if not os.path.exists('./calculations/'):
            os.makedirs('./calculations/')

        self.job_run_folder = FileChooser('./calculations/')#,layout=Layout(height='200px', width='500px'))
        self.job_run_folder.default_filename = ' '

        list_calculations = os.listdir('./calculations/')
        if len(list_calculations) == 0:
            self.available_calculations = create_dropdown_input(['None'],'None','')
        else:
            self.available_calculations = widgets.SelectMultiple(
                                options=list_calculations,
                                value=[list_calculations[0]],
                                description='',
                                disabled=False
                            )
            #self.available_calculations = create_dropdown_input(list_calculations,list_calculations[0],'')
        self.update_job_run_folder_button = create_expanded_button('Update','primary')
        self.update_job_run_folder_button.on_click(self.list_available_calculations)

        
        self.available_calculations.observe(self.select_available_calculation, 'value')
        self.select_available_calculation('aux')

        self.out_workflow = widgets.Output()
  
        self.run_siesta_box = wg.VBox([siesta_logo,HBox([Label('Select calculations folder (Default "Calculations"): '),self.job_run_folder,
         self.update_job_run_folder_button]),
                                        Label('Select calculations (Multiple values can be selected with shift + click):'),
                                        HBox([self.available_calculations, self.available_actions]),
                                       HBox([Label('N° processors: '), self.n_processors]),
                                       HBox([Label('Siesta command: '), self.siesta_command]),
                                      self.run_siesta_button,self.out_workflow],
                                       layout = {'width':'max-content'})  
        
        self.run_siesta_button.on_click(self.launch_siesta_run)
    
    def create_aditional_parameters_box(self):

        description = """<p>Write additional parameters that are not available in the previous options.</p>
<p>Please follow the same format as Siesta, i.e. {parameter} {value} or</p>
<p>{parameter} {value} {unit} (if required).For example:</p>
<p>&nbsp;</p>
<p style="padding-left: 40px;">MD.TargetPressure 0 GPa</p>
<p style="padding-left: 40px;">PAO.SoftDefault false</p>
<p style="padding-left: 40px;">&nbsp;</p>
<p>Blocks are also supported, for example:</p>
<p>&nbsp;</p>
<p style="padding-left: 40px;">%block GeometryConstraints<br />stress 1 2 4 5 6<br />%endblock GeometryConstraints</p>
<p>or</p>
<p style="padding-left: 40px;" >%block Geometry.Constraints<br />Z 6 # constrain Carbon<br />Z 1 1. 0. 0. # constrain Hydrogen along x Cartesian vector<br />%endblock</p> """
        description = wg.HTML(
        value=description)

        placeholder = """MD.TargetPressure 0 GPa
PAO.SoftDefault false"""
        self.aditional_parameters = widgets.Textarea(
                                            value='',
                                            placeholder=placeholder,
                                            description='New params:',
                                            disabled=False
                                        )
        self.aditional_parameters_box = wg.VBox([description, self.aditional_parameters])


    def create_siesta_bands_parameters(self):
        """
        Method for create widgets for input bands parameters
        """
        
        self.bandLineScale = create_dropdown_input(['pi/a','ReciprocalLatticeVectors'],'ReciprocalLatticeVectors','', '150px')

        self.bands_blocks_options = wg.RadioButtons(
                options=['bandLines', 'bandPoints'],
                value = 'bandLines',
                description='',
                disabled=False
            )
        
        self.blockBandLines = create_text_area_input(self.bandLinesExample,'enter bandLines block','')
        self.blockBandPoints = create_text_area_input(self.bandPointsExample,'enter bandpoints block','')
        self.writeKbands = create_dropdown_input(['false','true'],'false','')
        self.writeBands = create_dropdown_input(['false','true'],'false','')
        self.siesta_bands_parameters = wg.VBox([wg.Box([Label('bandLineScale'),self.bandLineScale]),
                                               wg.Box([Label('Do you want enter bandLines or bandPoints?'),self.bands_blocks_options]),
                                               wg.Box([Label('%block BandLines  '),self.blockBandLines]),
                                               wg.Box([Label('%block BandPoints'), self.blockBandPoints]),
                                               wg.Box([Label('writeKbands'), self.writeKbands]),
                                               wg.Box([Label('writeBands'), self.writeBands])])
        
        self.siesta_bands_parameters.layout.display = "none"
        self.bands_blocks_options.observe(self.display_siesta_band_blocks, 'value')
        self.bands.observe(self.get_siesta_bands_parameters,'value')
        
        
    def create_siesta_dos_parameters(self):
        """
        Method for create widgets for input DOS parameters
        """
        self.blockProjectedDensityOfStates = create_text_area_input(self.exampleProjectedDensityOfStates,
                                                                    'enter ProjectedDensityOfStatesBlock block','')
        

        self.blockPDOSkgrid_Monkhorst_Pack = create_text_area_input(self.examplePDOSkgrid_Monkhorst_Pack,
                                                                    'enter PDOS.kgrid_Monkhorst_Pack block','')
        
        self.pdos_kgrid_options = wg.RadioButtons(
                options=['Yes', 'No'],
                value = 'No',
                description='',
                disabled=False
            )
        
        self.siesta_dos_parameters = wg.VBox([wg.Box([Label('Enter PDOS.kgrid_Monkhorst_Pack?'),self.pdos_kgrid_options]),
                                               wg.Box([Label('%block ProjectedDensityOfStates  '),
                                                       self.blockProjectedDensityOfStates]),
                                               wg.Box([Label('%block PDOS.kgrid_Monkhorst_Pack'),
                                                       self.blockPDOSkgrid_Monkhorst_Pack ])])
        
        self.siesta_dos_parameters.layout.display = "none"
        self.pdos_kgrid_options.observe(self.display_siesta_pdos_kgrid, 'value')
        self.dos.observe(self.get_siesta_dos_parameters,'value')
        
    def create_siesta_opt_parameters(self):
        """
        Method for create widgets for input optical parameters
        """   
        self.OpticalCalculation = create_dropdown_input(['false','true'],'true','')
        self.OpticalEnergyMinimum = create_float_text(0, '')
        self.OpticalEnergyMaximum = create_float_text(10, '')
        self.OpticalBroaden =  create_float_text(0, '')
        self.OpticalScissor = create_float_text(0, '')
        
        self.options_OpticalNumberOfBands = create_dropdown_input(['Yes','No'],'No','')
        self.OpticalNumberOfBands = create_int_text(0, '')
 
        self.options_opticalMesh = create_dropdown_input(['Yes','No'],'No','')
        self.blockOpticalMesh =  create_text_area_input(self.exampleOpticalMesh,
                                                                    'enter Optical.Mesh block','')
        self.OpticalOffsetMesh = create_dropdown_input(['false','true'],'false','')
        self.OpticalPolarizationType = create_dropdown_input(['polarized','unpolarized','polycrystal'],'polycrystal','')
        

        self.blockOpticalVector = create_text_area_input(self.exampleOpticalVector,
                                                                    'enter Optical.Vector block','')
        
        self.optical_units_min = create_dropdown_input(['eV','Ang','nm','Bohr','Hartree','kJ','kcal','mol','Ry','meV'],
                                                 'Ry','')
        self.optical_units_max = create_dropdown_input(['eV','Ang','nm','Bohr','Hartree','kJ','kcal','mol','Ry','meV'],
                                                 'Ry','')
        self.optical_units_broaden = create_dropdown_input(['eV','Ang','nm','Bohr','Hartree','kJ','kcal','mol','Ry','meV'],
                                                 'Ry','')
        self.optical_units_scissor = create_dropdown_input(['eV','Ang','nm','Bohr','Hartree','kJ','kcal','mol','Ry','meV'],
                                                 'Ry','')

        self.siesta_opt_parameters = wg.VBox([wg.Box([Label('OpticalCalculation'),self.OpticalCalculation]),
                                               wg.Box([Label('Optical.Energy.Minimum'),self.OpticalEnergyMinimum, self.optical_units_min]),
                                               wg.Box([Label('Optical.Energy.Maximum'),self.OpticalEnergyMaximum, self.optical_units_max]),
                                               wg.Box([Label('Optical.Broaden'),self.OpticalBroaden, self.optical_units_broaden]),
                                               wg.Box([Label('Optical.Scissor '),self.OpticalScissor, self.optical_units_scissor]),
                                               wg.Box([Label('Enter different number of bands (default all bands)'),
                                                       self.options_OpticalNumberOfBands]),
                                               wg.Box([Label('Optical.NumberOfBands'),self.OpticalNumberOfBands]),
                                               wg.Box([Label('Optical.OffsetMesh'),self.OpticalOffsetMesh]),
                                               wg.Box([Label('Enter Optical Mesh block?'),self.options_opticalMesh]),
                                               wg.Box([Label('%block Optical.Mesh'),self.blockOpticalMesh]),
                                               wg.Box([Label('Polarization type'),self.OpticalPolarizationType]),
                                               wg.Box([Label('%block Optical.Vector'),self.blockOpticalVector]),])
        
        self.siesta_opt_parameters.layout.display = "none"
        self.options_opticalMesh.observe(self.display_siesta_optical_mesh, 'value')
        self.options_OpticalNumberOfBands.observe(self.display_optical_number_of_bands, 'value')
        self.OpticalPolarizationType.observe(self.display_optical_vector, 'value')
        
        self.opt.observe(self.get_siesta_opt_parameters,'value')

    def select_available_calculation(self, aux):
        option_list = []
        if aux == 'aux':
            calcs = self.available_calculations.value
        else:
            calcs = aux.new
        try:
            folder = self.job_run_folder.value.split('/ ')[0]
        except:
            folder = './calculations/'
        for calc in calcs:

            self.root_directory = f'{folder}/{calc}'
            self.structure_directory = f'{folder}/{calc}/structure/'

            # load config params
            with open(f'{self.root_directory}/config.json', 'r') as fp:
                config = json.load(fp)

            siesta_label = config['siesta_label']
            relax_structure = config['relax_structure']
            bands = config['bands']
            dos = config['dos']
            opt = config['optical']

            if relax_structure:
                option_list += [f'{calc} relax_structure']
            if bands:
                option_list += [f'{calc} bands']
            if dos:
                option_list += [f'{calc} dos']
            if opt:
                option_list += [f'{calc} opt']

        try:
            value = option_list[0]
        except: 
            option_list += [None]
            value = None
        if aux == 'aux':
            description = wg.HTML(value="Select actions:")
            self.select_actions = widgets.SelectMultiple(
                    options=option_list,
                    value=[value],
                    #rows=10,
                    description='',
                    disabled=False
                )
            self.available_actions = widgets.VBox([description, self.select_actions])
        else:
            self.select_actions.options = option_list
            self.select_actions.value = [value]
      
              
    def launch_siesta_run(self, aux):
        """
        Method for siesta job execution according to user's requirements
        """
        
        os.environ['OMP_NUM_THREADS'] = "1"
        try:
            folder = self.job_run_folder.value.split('/ ')[0]
        except:
            folder = './calculations'
        for calc in self.available_calculations.value:
            self.siesta_label = calc
            self.root_directory = f'{folder}/{self.siesta_label}'
            self.structure_directory = f'{folder}/{self.siesta_label}/structure/'

            # load config params
            with open(f'{self.root_directory}/config.json', 'r') as fp:
                config = json.load(fp)

            siesta_label = config['siesta_label']
            relax_structure = config['relax_structure']
            bands = config['bands']
            dos = config['dos']
            opt = config['optical']
            self.siesta_relax_inputs = config['relax_structure_inputs']
            self.siesta_bands_inputs = config['bands_inputs']
            self.siesta_dos_inputs = config['dos_inputs']
            self.siesta_opt_inputs = config['opt_inputs']
            self.pseudopotentials = config['pseudopotentials']
            
            if f"{calc} relax_structure" in list(self.select_actions.value):
                relax_structure = True
            else: 
                relax_structure = False
            if f"{calc} bands" in list(self.select_actions.value):
                bands = True
            else:
                bands = False
            if f"{calc} dos" in list(self.select_actions.value):
                dos = True
            else:
                dos = False
            if f"{calc} opt" in list(self.select_actions.value):
                opt = True
            else:
                opt = False

            input_ = f'{siesta_label}.fdf'
            output = f'{siesta_label}.out'
            
            # Get total energy
            if (relax_structure != True) and (bands!= True) and (dos != True) and (opt != True):
                with self.out_workflow:
                    self.launch_siesta_total_energy(input_, output, siesta_label)

            # Relax structure
            if relax_structure == True:
                with self.out_workflow :
                    self.launch_relax_run(input_, output, siesta_label)

            # Get band structure
            if bands == True: 
                with self.out_workflow :
                    self.launch_bands_run(input_, output, siesta_label)

            # Get density of states 
            if dos == True:                
                with self.out_workflow :
                    self.launch_dos_run(input_, output, siesta_label)
            
            #Get optical results
            if opt == True:         
                with self.out_workflow :
                    self.launch_opt_run(input_, output, siesta_label)

    def launch_siesta_total_energy(self,input_, output, siesta_label):
        """
        execute if any action is selected.
        """
        
        display(wg.HTML(value=f"{siesta_label}: Getting total energy..."))
        change_system_relax_control(f"{self.structure_directory}{siesta_label}.fdf",
                           f"{self.structure_directory}", siesta_label)

        for pseudopotential in self.pseudopotentials:
            os.system(f"cp src/pseudopotentials/{pseudopotential} {self.structure_directory}{pseudopotential}")
        os.system(f"cd {self.structure_directory} && mpirun -np {self.n_processors.value} {self.siesta_command.value} < {input_} > {output}")
        clear_output(wait=True)
        display(wg.HTML(value="Done."))
                
    def launch_relax_run(self,input_, output, siesta_label):
        """
        execute if relax is selected
        """
        
        display(wg.HTML(value=f"{siesta_label}: Relaxing Structure..."))

        for i in self.siesta_relax_inputs.split('\n'):
            try:
                variable = i.split(' ')[-2]
                value = i.split(' ')[-1]
                if variable == 'MD.NumCGsteps':
                    self.cg_steps_input_value = value
                elif variable == 'MD.VariableCell':
                    self.variable_cell_input_value = value
                elif variable == 'MD.TypeOfRun':
                    self.typerun_input_value = value
                elif variable == 'MD.MaxForceTol':
                    self.max_force_tol_input_value = value

            except: pass

        change_system_relax_control(f"{self.structure_directory}{siesta_label}.fdf",
                           f"{self.structure_directory}", siesta_label, str(self.cg_steps_input_value),
                            self.variable_cell_input_value, self.typerun_input_value, self.max_force_tol_input_value)

        for pseudopotential in self.pseudopotentials:
            os.system(f"cp src/pseudopotentials/{pseudopotential} {self.structure_directory}{pseudopotential}")
        os.system(f"cd {self.structure_directory} && mpirun -np {self.n_processors.value} {self.siesta_command.value} < {input_} > {output}")
        clear_output(wait=True)
        

        
        #display(wg.HTML(value="Getting total energy..."))

        check(self.structure_directory)
        df = pd.read_csv(self.structure_directory + 'Status.csv')[['System','TotalEnergy','Status']]
        if df['Status'][0] == 'False':
            check_image = get_check_image("./src/images/check-failed.png", "png")
            status = 'Unrelaxed'
        else:
            check_image = get_check_image("./src/images/check-success.png", "png")
            status = 'Relaxed'
        display(wg.HTML(value="<p><em><strong>Relax structure</strong></em>&nbsp;</p>"))
        #display(wg.HTML(value="Done..."))
        #display(wg.HTML(value="Checking output..."))
        output = wg.HTML(value=f"<p><em><strong>System:</strong></em> {siesta_label}&nbsp; &nbsp; &nbsp;<strong>Status:&nbsp;</strong></p>")
        display(wg.HBox([output, check_image]))  
        total_energy =  wg.HTML(value=f"<p><strong>{status} structure total energy: </strong>{df['TotalEnergy'][0]} eV</p>")
        display(total_energy)
        
    def launch_bands_run(self,input_, output, siesta_label):
        """
        execute if bands job is selected
        """
        display(wg.HTML(value=f"{siesta_label}: Calculating bands Structure..."))
        if os.path.isdir(f'{self.root_directory}/bands'):
            pass
        else:
            os.system(f"mkdir {self.root_directory}/bands")

        for pseudopotential in self.pseudopotentials:
            os.system(f"cp src/pseudopotentials/{pseudopotential} {self.root_directory}/bands/{pseudopotential}")
        try:
            os.system(f"cp {self.structure_directory}{siesta_label}-R.fdf {self.root_directory}/bands/{siesta_label}.fdf")
            change_system_relax_control(f"{self.structure_directory}{siesta_label}-R.fdf",
                           f"{self.root_directory}/bands", siesta_label)
        except:
            os.system(f"cp {self.structure_directory}{siesta_label}.fdf {self.root_directory}/bands/{siesta_label}.fdf")
            change_system_relax_control(f"{self.structure_directory}{siesta_label}.fdf",
                           f"{self.root_directory}/bands", siesta_label)

        os.system(f"cp {self.structure_directory}*.psf {self.root_directory}/bands/")

        os.system(f"echo '{self.siesta_bands_inputs}' >> {self.root_directory}/bands/{siesta_label}.fdf")
        os.system(f"cd {self.root_directory}/bands/ && mpirun -np {self.n_processors.value} {self.siesta_command.value} < {input_} > {output}")

        display(wg.HTML(value="<p><em><strong>Band structure</strong></em>&nbsp;</p>"))

        try: 
            fig, gap = self.get_bands_figure(siesta_label)
            gap =  wg.HTML(value=f"<p><strong>structure Band Gap value: </strong>{gap}</p>")
            check_image = get_check_image("./src/images/check-success.png", "png")
            status = 'Ok'
        except:
            check_image = get_check_image("./src/images/check-failed.png", "png")
            status = 'False'

        output = wg.HTML(value=f"<p><em><strong>System:</strong></em> {siesta_label}&nbsp; &nbsp; &nbsp;<strong>Status:&nbsp;</strong></p>")
        display(wg.HBox([output, check_image]))  
        if status == 'Ok':
            display(gap)
            display(fig)

    def get_bands_figure(self, filename):


        f = open(f'{self.root_directory}/bands/{filename}.bands', "r")
        text = f.read()

        info = []
        y_data_up = []
        y_data_down = []
        for i,line in enumerate(text.split('\n')):    
            if i == 0:
                e_f = float(line)
            if i == 1:
                k_min, k_max = float(line.split()[0]),float(line.split()[1])
            if i == 2:
                e_min, e_max = float(line.split()[0]),float(line.split()[1])
            if i == 3:
                nbands, nspin, nk = int(line.split()[0]),int(line.split()[1]),int(line.split()[2])
                x_data = np.linspace(k_min, k_max, nk)
            if i > 3:
                for value in line.split():
                    info += [float(value)]
                if (i-3)%round(nbands/10) == 0:
                    if nspin == 1:
                        y_data_up += [np.array(info[1:])-e_f]
                    else:
                        y_data_up += [np.array(info[1:nbands])-e_f]
                        y_data_down += [np.array(info[-nbands:])-e_f]
                    info = []
                if len([i[0] for i in y_data_up]) == nk:
                    break


        fig = plt.figure(figsize=(5, 7))
        ax = fig.add_subplot(1,1,1)
      
        positiveY = []
        negativeY = []
        verbose = False	
        for i in range(nbands):
            ax.plot([band[i] for band in y_data_up],color='blue',linewidth=0.5)
            y = np.array([band[i] for band in y_data_up])
            try:
                if all(i >= 0  for i in y) or all(i <= 0  for i in y):    
                    pass
                else:               
                    verbose = True
            except: pass

            try:
                positiveY += [np.min(y[y>0])]
            except:
                pass
            try:
                negativeY += [np.max(y[y<0])]
            except:
                pass

        if verbose:
            gap = 0

        else:
            if np.min(positiveY) < 1e-04 or np.abs(np.max(negativeY)) < 1e-04:
                gap = 0
            elif np.min(positiveY) < 1e-03 and np.abs(np.max(negativeY)) < 1e-03:
                gap = 0
            else:
                gap = np.min(positiveY) + np.abs(np.max(negativeY))
                if gap < 0.02:
                    gap = 0
        if nspin == 2:
            ax.plot([band[i] for band in y_data_down],color='red',linewidth=0.5)
        
        ax.set_title(f"{filename} - Band structure",loc='center')

        ax.set_ylim(-8,8)
        ax.set_xlabel(r'$K$ points')
        ax.set_ylabel(r'$E-E_{F} (eV)$')
        ax.set_xticklabels([])
        fig.savefig(f'{self.root_directory}/bands/bands.jpg')
        plt.close()
        return fig, gap
        
    def launch_dos_run(self,input_, output, siesta_label):
        """
        execute if DOS job is selected
        """
        display(wg.HTML(value=f"{siesta_label}: Calculating Density of Sates..."))
        if os.path.isdir(f'{self.root_directory}/dos'):
            pass
        else:
            os.system(f"mkdir {self.root_directory}/dos")

        for pseudopotential in self.pseudopotentials:
            os.system(f"cp src/pseudopotentials/{pseudopotential} {self.root_directory}/dos/{pseudopotential}")
        try:
            os.system(f"cp {self.structure_directory}{siesta_label}-R.fdf {self.root_directory}/dos/{siesta_label}.fdf")
            change_system_relax_control(f"{self.structure_directory}{siesta_label}-R.fdf",
                           f"{self.root_directory}/dos", siesta_label)
        except:
            os.system(f"cp {self.structure_directory}{siesta_label}.fdf {self.root_directory}/dos/{siesta_label}-R.fdf")
            change_system_relax_control(f"{self.structure_directory}{siesta_label}.fdf",
                           f"{self.root_directory}/dos", siesta_label)

        os.system(f"cp {self.structure_directory}*.psf {self.root_directory}/dos/")
        

        os.system(f"echo '{self.siesta_dos_inputs}' >> {self.root_directory}/dos/{siesta_label}.fdf")
        os.system(f"cd {self.root_directory}/dos/ && mpirun -np {self.n_processors.value} {self.siesta_command.value} < {input_} > {output}")
        #clear_output(wait=True)
        #display(wg.HTML(value="Done..."))
        #display(wg.HTML(value="Checking output..."))
        display(wg.HTML(value="<p><em><strong>Density of States</strong></em>&nbsp;</p>"))
        check(f"{self.root_directory}/dos/")

        try: 
            fig, eig = self.get_dos_figure(siesta_label)
            fermi_energy =  wg.HTML(value=f"<p><strong>structure fermi energy: </strong>{eig} eV</p>")
            check_image = get_check_image("./src/images/check-success.png", "png")
            status = 'Ok'
        except:
            check_image = get_check_image("./src/images/check-failed.png", "png")
            status = 'False'

        output = wg.HTML(value=f"<p><em><strong>System:</strong></em> {siesta_label}&nbsp; &nbsp; &nbsp;<strong>Status:&nbsp;</strong></p>")
        display(wg.HBox([output, check_image]))  
        if status == 'Ok':
            display(fermi_energy)
            display(fig)
        
    
    def get_dos_figure(self, filename):

        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(1,1,1)
        val = np.loadtxt(f'{self.root_directory}/dos/{filename}.DOS')

        f = open(f'{self.root_directory}/dos/{filename}.EIG','r')
        eig= float(f.readline())

        val[:,0] = val[:,0]-eig

        val[:,1] = val[:,1]
        try:
            val[:,2] = -1*val[:,2]
        except:pass
        xmax = np.ceil(np.ndarray.max(val[:,0]))
        xmin = np.floor(np.ndarray.min(val[:,0]))
        ymax = np.ceil(np.ndarray.max(val[:,1]))
        ymin = np.floor(np.ndarray.min(val[:,1]))

        ax.set_title(f"{filename} - Density of States",loc='center')
        ax.plot(val[:,0],val[:,1],color='blue',linewidth=0.5)
        try:  
            ax.plot(val[:,0],val[:,2],color='red',linewidth=0.5)
        except:pass

        ax.set_xlabel(r'$E-E_{F} (eV)$')
        ax.set_ylabel(r'DOS $(a.u)$')
        ax.axvline(0, -800, 2600, color = 'green')
        fig.savefig(f'{self.root_directory}/dos/dos.jpg')
        plt.close()
        return fig, eig

    def launch_opt_run(self,input_, output, siesta_label):
        """
        execute if optical job is selected
        """
        display(wg.HTML(value=f"{siesta_label}: Calculating Optical"))
        if os.path.isdir(f'{self.root_directory}/optical'):
            pass
        else:
            os.system(f"mkdir {self.root_directory}/optical")

        for pseudopotential in self.pseudopotentials:
            os.system(f"cp src/pseudopotentials/{pseudopotential} {self.root_directory}/optical/{pseudopotential}")
        try:
            os.system(f"cp {self.structure_directory}{siesta_label}-R.fdf {self.root_directory}/optical/{siesta_label}.fdf")
            change_system_relax_control(f"{self.structure_directory}{siesta_label}-R.fdf",
                           f"{self.root_directory}/optical", siesta_label)
        except:
            os.system(f"cp {self.structure_directory}{siesta_label}.fdf {self.root_directory}/optical/{siesta_label}-R.fdf")
            change_system_relax_control(f"{self.structure_directory}{siesta_label}.fdf",
                           f"{self.root_directory}/optical", siesta_label)

        os.system(f"cp {self.structure_directory}*.psf {self.root_directory}/optical/")
        

        os.system(f"echo '{self.siesta_opt_inputs}' >> {self.root_directory}/optical/{siesta_label}.fdf")
        os.system(f"cd {self.root_directory}/optical/ && mpirun -np {self.n_processors.value} {self.siesta_command.value} < {input_} > {output}")

        #clear_output(wait=True)
        #display(wg.HTML(value="Done..."))
        #display(wg.HTML(value="Checking output..."))
        display(wg.HTML(value="<p><em><strong>Optical calculation</strong></em>&nbsp;</p>"))
        
        try: 
            fig = self.get_opt_figure(siesta_label)
            check_image = get_check_image("./src/images/check-success.png", "png")
            status = 'Ok'
        except:
             check_image = get_check_image("./src/images/check-failed.png", "png")
             status = 'False'

        output = wg.HTML(value=f"<p><em><strong>System:</strong></em> {siesta_label}&nbsp; &nbsp; &nbsp;<strong>Status:&nbsp;</strong></p>")
        display(wg.HBox([output, check_image]))  
        if status == 'Ok':
            display(fig)


    def get_opt_figure(self, filename):

        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(1,1,1)
        val = np.loadtxt(f'{self.root_directory}/optical/{filename}.EPSIMG')

        try:
            val[:,1] = (val[:,1] + val[:,2])
        except:pass
        xmax = np.ceil(np.ndarray.max(val[:,0]))
        xmin = np.floor(np.ndarray.min(val[:,0]))
        ymax = np.ceil(np.ndarray.max(val[:,1]))
        ymin = np.floor(np.ndarray.min(val[:,1]))

        ax.set_title(f"{filename} - Optical: imaginary part of the dielectric function",loc='center')
        ax.plot(val[:,0],val[:,1],color='blue',linewidth=0.5)


        ax.set_xlabel(r'Photon energy $eV$')
        ax.set_ylabel(r'$\epsilon_2$ $(a.u)$')
        ax.set_xlim(0,10)
        fig.savefig(f'{self.root_directory}/optical/imaginary-part.jpg')
        plt.close()
        return fig

    def get_geometric_constraints_parameters(self,aux):
        if self.constraints.value == True:
            self.geometric_constraints_parameters.layout.display = "block"
        else:
            self.geometric_constraints_parameters.layout.display  = "none"

    def get_siesta_relax_parameters(self, aux):

        if self.relax_structure.value == True:
            self.system_relax_control_box.layout.display = "block"
        else:
            self.system_relax_control_box.layout.display = "none"

    def get_siesta_opt_parameters(self, aux):
        """
        Method for display dinamically optical inputs according to user's interactions.
        'block' displays the widget, otherwise 'none' hidde it
        """
        if self.opt.value == True:
            self.siesta_opt_parameters.layout.display = "block"
            self.blockOpticalVector.layout.display = "none"
            self.OpticalNumberOfBands.layout.display = "none"
            self.blockOpticalMesh.layout.display = "none"
            
        else:
            self.siesta_opt_parameters.layout.display = "none"
            self.blockOpticalVector.layout.display = "none" 
            self.OpticalNumberOfBands.layout.display = "none"
            self.blockOpticalMesh.layout.display = "none"
        
    def get_siesta_dos_parameters(self, aux):
        """
        Method for display dinamically DOS inputs according to user's interactions.
        'block' displays the widget, otherwise 'none' hidde it
        """
        if self.dos.value == True:
            self.siesta_dos_parameters.layout.display = "block"
            self.blockPDOSkgrid_Monkhorst_Pack.layout.display = "none"
            
        else:
            self.siesta_dos_parameters.layout.display = "none"
            self.blockPDOSkgrid_Monkhorst_Pack.layout.display = "none"
                        
            
    def get_siesta_bands_parameters(self, aux):
        """
        Method for display dinamically bands inputs according to user's interactions.
        'block' displays the widget, otherwise 'none' hidde it
        """
        if self.bands.value == True:
            self.siesta_bands_parameters.layout.display = "block"
            self.blockBandLines.layout.display = "block"
            self.blockBandPoints.layout.display = "none"
        else:
            self.siesta_bands_parameters.layout.display = "none"
            self.blockBandLines.layout.display = "none"
            self.blockBandPoints.layout.display = "none"
            
            
    def process_siesta_bands_input(self):
        """
        Method for process the user's band inputs
        """
        
        bandLineScale = f"BandLinesScale {self.bandLineScale.value}"
        bands_blocks_options = self.bands_blocks_options.value
        
        if bands_blocks_options ==  'bandLines':
            blockBandx = self.blockBandLines.value
        elif bands_blocks_options == 'bandPoints':
            blockBandx = self.blockBandPoints.value
            
        writekbands = f"WriteKbands {self.writeKbands.value}"
        writebands = f"WriteBands {self.writeBands.value}"     
        self.siesta_bands_inputs = f"\n{bandLineScale}\n {blockBandx}\n {writekbands}\n {writebands}"
    
    def process_siesta_relax_input(self):
        num_cg_steps = f"MD.NumCGsteps {self.cg_steps_input.value}"
        variable_cell = f"MD.VariableCell {self.variable_cell_input.value}"
        type_run = f"MD.TypeOfRun {self.typeRun_input.value}"
        max_force_tol = f"MD.MaxForceTol {self.max_force_tol_input.value}"


        self.siesta_relax_inputs = f"\n{num_cg_steps}\n {variable_cell}\n {type_run}\n {max_force_tol}\n"
    
    def process_siesta_dos_input(self):
        """
        Method for process the user's dos inputs
        """  
        pdos_kgrid_options = self.pdos_kgrid_options.value
        
        if pdos_kgrid_options == 'Yes':
            blockPDOSkgrid_Monkhorst_Pack = self.blockPDOSkgrid_Monkhorst_Pack.value
        else:
            blockPDOSkgrid_Monkhorst_Pack = ''
            
        blockProjectedDensityOfStates = self.blockProjectedDensityOfStates.value
        
        self.siesta_dos_inputs = f"\n{blockProjectedDensityOfStates}\n {blockPDOSkgrid_Monkhorst_Pack}"
    
    def process_siesta_opt_input(self):
        """
        Method for process the user's optical inputs
        """
        opticalCalculation = f'OpticalCalculation {self.OpticalCalculation.value}'
        opticalEnergyMinimum = f"Optical.Energy.Minimum {self.OpticalEnergyMinimum.value} {self.optical_units_min.value}"
        opticalEnergyMaximum = f"Optical.Energy.Maximum {self.OpticalEnergyMaximum.value} {self.optical_units_max.value}"
        opticalBroaden = f"Optical.Broaden {self.OpticalBroaden.value} {self.optical_units_broaden.value}"
        opticalScissor = f"Optical.Scissor {self.OpticalScissor.value} {self.optical_units_scissor.value}"
        
        if self.options_OpticalNumberOfBands.value == 'Yes':
            opticalNumberOfBands = f"Optical.NumberOfBands {self.OpticalNumberOfBands.value}"
        else:
            opticalNumberOfBands = ''
        
        opticalOffsetMesh = f"Optical.OffsetMesh {self.OpticalOffsetMesh.value}"
        
        if self.options_opticalMesh.value == 'Yes':
            blockOpticalMesh = self.blockOpticalMesh.value
        else:
            blockOpticalMesh = ''
        
        polarizationType = self.OpticalPolarizationType.value
        
        if polarizationType in ['polarized', 'unpolarized']:
            blockOpticalVector = self.blockOpticalVector.value
        else:
            blockOpticalVector = ''
            
        self.siesta_opt_inputs = f"\n{opticalCalculation}\n {opticalEnergyMinimum}\n {opticalEnergyMaximum}\n {opticalBroaden}\n {opticalScissor}\n {opticalNumberOfBands}\n {opticalOffsetMesh}\n {blockOpticalMesh}\n Optical.PolarizationType {polarizationType}\n {blockOpticalVector}\n"                                                     
      
        
    def update_siesta_author(self,*args):
        """
        Method for update dinamically the author parameter according to xc input selection
        """
        if self.xc_input.value=='GGA':
            self.author.options = ['PBE','PW91','revPBE','RPBE','WC','AM05','PBEsol',
                                   'PBEJsJrLO','PBEGcGxLO','PBEGcGxHEG','BLYP']
        if self.xc_input.value=='LDA':
            self.author.options = ['CA','PW92']
        if self.xc_input.value=='VDW':
            self.author.options = ['DRSLL','LMKLL','KBM','C09','BH','VV']
            

    def check_siesta_kpoints_input(self,*args):
        """
        Method for check and display kpoints inputs according to user's selection
        """    
        if self.kpts_options_input.value == 'Gamma':
            self.kpts_input.layout.display = "none"
            self.k_9.value = 1
            self.converge_kpts_input.layout.display = "none"
            
        elif self.kpts_options_input.value == 'Enter kgrid_Monkhorst_Pack':
            self.kpts_input.layout.display = "block"
            self.converge_kpts_input.layout.display = "none"
            
        else:
            self.kpts_input.layout.display = "none"
            self.converge_kpts_input.layout.display = "none"
            
    def check_siesta_mesh_cutoff_input(self,*args):
        """
        Method for check and display Mesh cutoff inputs according to user's selection
        """       
        if self.mesh_cutoff_options_input.value == 'Enter Mesh cutoff':
            self.converge_mesh_cutoff_input.layout.display = "none"
            self.mesh_cutoff_input.layout.display = "block"
            
        else:
            
            self.mesh_cutoff_input.layout.display = "none"
            self.converge_mesh_cutoff_input.layout.display = "none"
      
    def display_siesta_pdos_kgrid(self, aux):   
        """
        Method for display (or not) the DOS kgrid input parameter according to user's señection
        """
        if self.pdos_kgrid_options.value == 'Yes':
            self.blockPDOSkgrid_Monkhorst_Pack.layout.display = "block"
            
        elif self.pdos_kgrid_options.value == 'No':
            self.blockPDOSkgrid_Monkhorst_Pack.layout.display = "none"
            
    def display_siesta_band_blocks(self, aux):   
        """
        Method for display (or not) bands inputs parameters according to user's selection
        """ 
        if self.bands_blocks_options.value == 'bandLines':
            
            self.blockBandLines.layout.display = "block"
            self.blockBandPoints.layout.display = "none"
            
        elif self.bands_blocks_options.value == 'bandPoints':
            
            self.blockBandLines.layout.display = "none"
            self.blockBandPoints.layout.display = "block"
            
    def display_siesta_optical_mesh(self, aux):
        """
        Method for display (or not) the optical Mesh block input parameter according to user's selection
        """ 
        if self.options_opticalMesh.value == 'Yes':
            self.blockOpticalMesh.layout.display = 'block'
            
        elif self.options_opticalMesh.value == 'No':
            self.blockOpticalMesh.layout.display = 'None'
            
    def display_optical_number_of_bands(self, aux):
        """
        Method for display (or not) the optical number of bands input parameter according to user's selection
        """
        if self.options_OpticalNumberOfBands.value == 'Yes':
            self.OpticalNumberOfBands.layout.display = 'block'
            
        elif self.options_OpticalNumberOfBands.value == 'No':
            self.OpticalNumberOfBands.layout.display = 'None'  
            
    def display_optical_vector(self, aux):
        """
        Method for display (or not) the optical vector input parameter according to user's selection
        """ 
        if self.OpticalPolarizationType.value in ['polarized','unpolarized']:
            self.blockOpticalVector.layout.display = 'block'
        else:
            self.blockOpticalVector.layout.display = 'none'
    
    def siesta_postprocessing(self):
        """
        Display post-processing options: Optical results, band structure, density of states 
        """
        
        siesta_logo = get_calculator_logo("./src/images/siesta.png", "png")
        
        self.siesta_optical_results()
        self.siesta_bands_results()
        self.siesta_dos_results()
        
        self.supported_postprocessing = wg.Accordion(children = [self.dos_result_box,
                                                              self.bands_result_box,
                                                              self.optical_result_box
                                                              ])

        self.supported_postprocessing.set_title(2, 'Optical results')
        self.supported_postprocessing.set_title(1, 'Band structure results')
        self.supported_postprocessing.set_title(0, 'Density of states results')
        self.postprocessing_siesta_box = wg.VBox([siesta_logo, self.supported_postprocessing])  

    def siesta_optical_results(self):
        """
        Display post-processing optical result block
        """
        self.graph_type = 'Optical'
        self.optical_result_box = wg.VBox([Postprocessing_graphs(self.graph_type).main_box]) 
        
    def siesta_bands_results(self):
        """
        Display post-processing band structure result block
        """
        self.graph_type = 'Bands'
        siesta_logo = get_calculator_logo("./src/images/siesta.png", "png")
        self.bands_result_box = wg.VBox([Postprocessing_graphs(self.graph_type).main_box]) 
        
    def siesta_dos_results(self):
        """
        Display post-processing density of states result block
        """
        self.graph_type = 'Dos'
        siesta_logo = get_calculator_logo("./src/images/siesta.png", "png")
        self.dos_result_box = wg.VBox([Postprocessing_graphs(self.graph_type).main_box]) 


def make_box_layout():
    """
    Function to create a predefined box layout
    """
    return widgets.Layout(
        border='solid 1px black',
        margin='0px 10px 10px 0px',
        padding='5px 5px 5px 5px')


def get_random_hex_colors(number_colors):
    """
    Function to choose random hexadecimal color 
    """
    hexadecimals = []
    for i in range(number_colors):
        hexadecimals += ["#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])]
    return hexadecimals

class Postprocessing_graphs(wg.HBox):
    """
    Creates / Edits postprocessing graphs
    """
     
    def __init__(self, graph_type):
        """
        graph_type = type of graph: optical, bands or dos
        """
        super().__init__()
        #create output box to display graph
        self.output = widgets.Output() 

        #init figure
        self.fig, self.ax = plt.subplots(constrained_layout=True, figsize=(5, 3.5))
        self.graph_type = graph_type
        self.ax.grid(False)
        plt.close(self.fig)
        self.fig.canvas.toolbar_position = 'bottom'   

        with self.output:
            clear_output(wait=True)
            display(self.fig)

        # get accepted_format according to graph type

        self.separation_lines = create_float_text(0, '')
        if self.graph_type == 'Optical':
            accepted_format = '*.EPSIMG'

        if self.graph_type == 'Dos':
            accepted_format = '*.DOS'

        if self.graph_type == 'Bands':
            accepted_format = '*.dat'

        # Create upload file widget
        self.upload_file = widgets.FileUpload(
                accept=accepted_format,  # Accepted file extension e.g. '.txt', '.pdf', 'image/*', 'image/*,.pdf'
                multiple=True  # True to accept multiple files upload else False
            )

        # Create options widgets
        text_xlabel = widgets.Text(
            value='', 
            description='xlabel', 
            continuous_update=False
        )
        
        text_ylabel = widgets.Text(
            value='', 
            description='ylabel', 
            continuous_update=False
        )
        
        text_title = widgets.Text(
            value='', 
            description='Title', 
            continuous_update=False
        )

        self.input_xticks_positions = widgets.Textarea(description = 'xticks positions:')
        self.input_xticks_labels = widgets.Textarea(description = 'xticks labels:')
        self.input_yticks_positions = widgets.Textarea(description = 'yticks positions:')
        self.input_yticks_labels = widgets.Textarea(description = 'yticks labels:')
        self.input_xlim = widgets.FloatRangeSlider(
            value=[5, 7.5],
            min=-50,
            max=50.0,
            step=0.1,
            description='xlim:',
            disabled=False,
            continuous_update=False,
            orientation='horizontal',
            readout=True,
            readout_format='.1f',)
        
        self.input_ylim = widgets.FloatRangeSlider(
            value=[5, 7.5],
            min=-50,
            max=50.0,
            step=0.1,
            description='ylim:',
            disabled=False,
            continuous_update=False,
            orientation='horizontal',
            readout=True,
            readout_format='.1f',)

        self.line_options_list = []
        self.line_options = widgets.VBox(self.line_options_list, layout=Layout(width='320px'))

        self.update_graph = create_expanded_button('Update graph','primary')
        self.update_graph.layout.display = 'none'

        self.display_legend = create_checkbox_input('Legend')
        display_grid = create_checkbox_input('Grid')
        
        #file chooser and save button for save plot as image
        self.fileChooser = FileChooser('./')#,layout=Layout(height='200px', width='500px'))
        self.fileChooser.default_filename = 'output_fig.eps'
        self.savebutton = widgets.Button(
                        description='',
                        disabled=False,
                        button_style='', # 'success', 'info', 'warning', 'danger' or ''
                        tooltip='Savefig',
                        icon='save', # (FontAwesome names without the `fa-` prefix)
                        layout=Layout(height='25px', width='50px')
                    )
        self.savebutton.on_click(self.savefig)

        #Display graphs according to its type
        if self.graph_type == 'Optical':
            
            controls = widgets.VBox([
                widgets.HBox([Label('Separation between lines:'),self.separation_lines]),
            widgets.HBox([Label('Select optical file result:'),self.upload_file]),
                self.line_options,self.update_graph, 
                text_xlabel, 
                text_ylabel,
                text_title,
                self.display_legend,
                display_grid,   
                self.input_xlim,
                self.input_ylim,
                self.input_xticks_positions,
                self.input_xticks_labels,
                self.input_yticks_positions,
                self.input_yticks_labels])

        if self.graph_type == 'Bands':
            self.is_other_file = False
            self.upload_spin_down_file = widgets.FileUpload(
                accept='.dat',  # Accepted file extension e.g. '.txt', '.pdf', 'image/*', 'image/*,.pdf'
                multiple=False  # True to accept multiple files upload else False
            )
            self.spin = create_checkbox_input('Spin polarized')
            self.spin_down = widgets.HBox([Label('Opcional if Spin: Select Spin down bands file result:'),
                                            self.upload_spin_down_file ])
            controls = widgets.VBox([
                self.spin,
            widgets.HBox([Label('Select bands or Spin Up (if Spin) file result:'),self.upload_file]),
            self.spin_down,
                self.line_options,self.update_graph, 
                text_xlabel, 
                text_ylabel,
                text_title,
                self.display_legend,
                display_grid,   
                self.input_xlim,
                self.input_ylim,
                self.input_xticks_positions,
                self.input_xticks_labels,
                self.input_yticks_positions,
                self.input_yticks_labels])
            self.spin_down.layout.display = 'None'
            self.spin.observe(self.bands_with_spin, 'value')
            self.upload_spin_down_file.observe(self.update_spin_down_data, 'value')

        if self.graph_type == 'Dos' :
            self.available_eigs = {}
            self.upload_eig_file = widgets.FileUpload(
                accept='.EIG',  # Accepted file extension e.g. '.txt', '.pdf', 'image/*', 'image/*,.pdf'
                multiple=True  # True to accept multiple files upload else False
            )
            self.eig_float = create_float_text(0,'')
            controls = widgets.VBox([
                widgets.HBox([Label('Separation between lines:'),self.separation_lines]),
             widgets.HBox([Label('Select EIG file:'),self.upload_eig_file]),
               widgets.HBox([Label('Select DOS file:'),self.upload_file]),
                self.line_options,self.update_graph, 
                text_xlabel, 
                text_ylabel,
                text_title,
                self.display_legend,
                display_grid,   
                self.input_xlim,
                self.input_ylim,
                self.input_xticks_positions,
                self.input_xticks_labels,
                self.input_yticks_positions,
                self.input_yticks_labels])
                
            self.upload_eig_file.observe(self.update_eig_value_from_file, 'value')

    
        controls.layout = make_box_layout()
         
        self.out_box = widgets.Box([self.output])
        self.output.layout = make_box_layout()
 
        # observe stuff

        self.upload_file.observe(self.update_data, 'value')

        text_xlabel.observe(self.update_xlabel, 'value')
        text_ylabel.observe(self.update_ylabel, 'value')
        text_title.observe(self.update_title, 'value')
        self.display_legend.observe(self.update_legend, 'value')
        display_grid.observe(self.update_grid, 'value')
        self.input_xlim.observe(self.update_xlim, 'value')
        self.input_ylim.observe(self.update_ylim, 'value')
        self.input_xticks_positions.observe(self.update_xticks_positions, 'value')
        self.input_xticks_labels.observe(self.update_xticks_labels, 'value')
        self.input_yticks_positions.observe(self.update_yticks_positions, 'value')
        self.input_yticks_labels.observe(self.update_yticks_labels, 'value')
        self.separation_lines.observe(self.update_separation_line_data, 'value')
        self.update_graph.on_click(self.manual_update_fig)
        
   
        text_xlabel.value = 'x'
        text_ylabel.value = 'y' 
        
        # add to children
        self.main_box = widgets.HBox([controls,
                                 self.out_box])

    def bands_with_spin(self, change):
        self.is_other_file = change.new
        self.spin_down.layout.display = 'block'

    def update_spin_down_data(self, change):
        self.update_data_bands(change)

    def update_eig_value_from_file(self,change):
        """
        Method for update fermin value from eig file if DOS graph is selected
        """

        uploaded_file = self.upload_eig_file.value
        filename = next(iter(uploaded_file))
        
        try: format_ = filename.split('.')[1]
        except: format_ = filename.split('.')[0]

        self.number_files = len(change.new.keys())
        self.upload_eig_file._counter = self.number_files
        for new_line in range(self.number_files):
            key = list(change.new.keys())[new_line]
            content = uploaded_file[key]['content'].decode("utf-8") 
            eig = float(content.split('\n')[0].replace(' ','')) 
            self.available_eigs[filename.split('.EIG')[0]] = eig
            

    def update_separation_line_data(self, change):
        """
        Method for update separation value between lines
        """
        self.update_data(self.upload_file.value)

    def update_data(self, change):
        """
        Method for update data according to selected graph type
        """

        if self.graph_type == 'Optical':
            self.update_data_optical(change)

        if self.graph_type == 'Bands': 
            self.update_data_bands(change)

        if self.graph_type == 'Dos' :
            uploaded_file = self.upload_file.value
            filename = next(iter(uploaded_file))
            display(self.available_eigs)
            self.available_eigs[filename.split('.DOS')[0]]
            self.update_data_dos(change)

    def update_data_bands(self, change):
        """
        Method for update data if bands graph is selected
        """
        
        try:
            change_new = change.new
        except:
            change_new = change

        
        self.number_files = len(change_new.keys())
        self.upload_file._counter = self.number_files
        
        
        
        if not self.is_other_file:
            self.line_options_list = []
            self.ax.clear()

        #Select initial color for each graph
        colors = get_random_hex_colors(self.number_files)
        #Separation value between lines or plots
        separation = self.separation_lines.value
        # loop over the band file or files uploaded
        for new_line in range(self.number_files):
            key = list(change_new.keys())[new_line]
            #get content of each file
            content ="".join(map(chr, change_new[key]['content']))
            content = StringIO(content)
            val = np.loadtxt(content)
            color = colors[new_line]
            label = key
            
            valx = val[:,0] 
            valy =val[:,1] 
            data = [valx,valy ]
            
            #Create optical line object with label, color picker and linestyle options widgets
            band_line = Line_graph(data, self.ax, label,self.graph_type , color = color)
            line, color_picker, label, linestyle_option = band_line.add_new_line()
            #Update values
            self.ax.set_xlim(np.min(valx), np.max(valx))
            self.input_xticks_positions.value = str(self.ax.get_xticks()).replace('[ ','').replace(']','').replace('[','').replace(' ]','')
            self.input_xticks_labels.value = str(self.ax.get_xticks()).replace('[ ','').replace(']','').replace('[','').replace(' ]','')
            self.input_yticks_positions.value = str(self.ax.get_yticks()).replace('[ ','').replace(']','').replace('[','').replace(' ]','')
            self.input_yticks_labels.value = str(self.ax.get_yticks()).replace('[ ','').replace(']','').replace('[','').replace(' ]','') 
            self.input_xlim.value = [np.min(valx), np.max(valx)]
            self.ax.set_ylim(np.min(valy), np.max(valy))
            self.input_ylim.value = [np.min(valy), np.max(valy)]
            self.line_options_list += [HBox([label,color_picker, linestyle_option])]
 
        #Update figure
        self.update_line_options('value')
        self.fig.canvas.draw()
        self.update_fig()

    def update_data_optical(self, change):
        """
        Method for update data if optical graph is selected
        """
        
        try:
            change_new = change.new
        except:
            change_new = change

        
        self.number_files = len(change_new.keys())
        self.upload_file._counter = self.number_files
        self.ax.clear()
        self.line_options_list = []
        #Select initial color for each graph
        colors = get_random_hex_colors(self.number_files)
        #Separation value between lines or plots
        separation = self.separation_lines.value
        # loop over the optical file or files uploaded
        for new_line in range(self.number_files):
            key = list(change_new.keys())[new_line]
            #get content of each file
            content ="".join(map(chr, change_new[key]['content']))
            content = StringIO(content)
            val = np.loadtxt(content)
            color = colors[new_line]
            label = key
            
            try:
                #Calculation with Spin polarization
                valx = val[:,0] 
                valy =val[:,1]+val[:,2] + separation*new_line
                data = [valx,valy ]
            except:
                #Calculation without spin polarization
                valx = val[:,0] 
                valy =val[:,1] + separation*new_line
                data = [valx,valy ]
            
            #Create optical line object with label, color picker and linestyle options widgets
            optical_line = Line_graph(data, self.ax, label,self.graph_type , color = color)
            line, color_picker, label, linestyle_option = optical_line.add_new_line()
            #Update values
            self.ax.set_xlim(np.min(valx), np.max(valx))
            self.input_xticks_positions.value = str(self.ax.get_xticks()).replace('[ ','').replace(']','').replace('[','').replace(' ]','')
            self.input_xticks_labels.value = str(self.ax.get_xticks()).replace('[ ','').replace(']','').replace('[','').replace(' ]','')
            self.input_yticks_positions.value = str(self.ax.get_yticks()).replace('[ ','').replace(']','').replace('[','').replace(' ]','')
            self.input_yticks_labels.value = str(self.ax.get_yticks()).replace('[ ','').replace(']','').replace('[','').replace(' ]','') 
            self.input_xlim.value = [np.min(valx), np.max(valx)]
            self.ax.set_ylim(np.min(valy), np.max(valy))
            self.input_ylim.value = [np.min(valy), np.max(valy)]
            self.line_options_list += [HBox([label,color_picker, linestyle_option])]
 
        #Update figure
        self.update_line_options('value')
        self.fig.canvas.draw()
        self.update_fig()

    def update_data_dos(self,change):
        """
        Method for update data if DOS graph is selected
        """
        try:
            change_new = change.new
        except:
            change_new = change

        
        self.number_files = len(change_new.keys())
        self.upload_file._counter = self.number_files
        self.ax.clear()
        self.line_options_list = []
        #Select initial color for each graph
        colors = get_random_hex_colors(self.number_files)
        #Separation value between lines or plots
        separation = self.separation_lines.value
        # loop over the optical file or files uploaded
        for new_line in range(self.number_files):
            key = list(change_new.keys())[new_line]
            #get content of each file
            content ="".join(map(chr, change_new[key]['content']))
            content = StringIO(content)
            val = np.loadtxt(content)
            color = colors[new_line]
            label = key
            try:
                #Calculation with Spin polarization
                valx = val[:,0] 
                valy_1 = val[:,1]+ separation*new_line
                valy_2 = -1*val[:,2] + separation*new_line
                valsy = [valy_1, valy_2]
            except:
                #Calculation without spin polarization
                valx = val[:,0] 
                valsy =val[:,1] + separation*new_line
            if len(valsy) ==2:
                labels = [f"{label}-Spin up", f"{label}-Spin down"]
            else:
                labels = [label]
            count = 0
            for valy in valsy:
                data = [valx, valy]
                label = labels[count]
                 #Create optical line object with label, color picker and linestyle options widgets
                dos_line = Line_graph(data, self.ax, label,self.graph_type, color = color)
                count += 1
                line, color_picker, label, linestyle_option = dos_line.add_new_line()
                #Update values
                self.ax.set_xlim(np.min(valx), np.max(valx))
                self.input_xticks_positions.value = str(self.ax.get_xticks()).replace('[ ','').replace(']','').replace('[','').replace(' ]','')
                self.input_xticks_labels.value = str(self.ax.get_xticks()).replace('[ ','').replace(']','').replace('[','').replace(' ]','')
                self.input_yticks_positions.value = str(self.ax.get_yticks()).replace('[ ','').replace(']','').replace('[','').replace(' ]','')
                self.input_yticks_labels.value = str(self.ax.get_yticks()).replace('[ ','').replace(']','').replace('[','').replace(' ]','') 
                self.input_xlim.value = [np.min(valx), np.max(valx)]
                self.ax.set_ylim(np.min(valy), np.max(valy))
                self.input_ylim.value = [np.min(valy), np.max(valy)]
                self.line_options_list += [HBox([label,color_picker, linestyle_option])]
 
        #Update figure
        self.update_line_options('value')
        self.fig.canvas.draw()
        self.update_fig()

    def update_line_options(self, change):
        """
        Method for update line options like label, color_picker and linestyle option
        """
        self.line_options.children = tuple(self.line_options_list)
        self.update_graph.layout.display = 'block'
        self.update_fig()

    def update_xlabel(self, change):
        """
        Method for update x label
        """
        self.ax.set_xlabel(change.new)
        self.update_fig()
  
    def update_ylabel(self, change):
        """
        Method for update y label
        """

        self.ax.set_ylabel(change.new)
        self.update_fig()
    
    def update_xticks_positions(self, change):
        """
        Method for update xticks positions
        """

        xticks = change.new.split()
        try:
            xticks = [float(i) for i in xticks]
        
            self.ax.set_xticks(xticks)
            self.update_fig()
        except:
             pass

    def update_xticks_labels(self, change):
        """
        Method for update xticks labels
        """
        xticks = change.new.split()
        try:
            self.ax.set_xticklabels(xticks)
            self.update_fig()
        except:
             pass
    def update_yticks_positions(self, change):
        """
        Method for update yticks positions
        """

        yticks = change.new.split()
        try:
            yticks = [float(i) for i in yticks]
        
            self.ax.set_yticks(yticks)
            self.update_fig()
        except:
             pass
    def update_yticks_labels(self, change):
        """
        Method for update yticks labels
        """

        yticks = change.new.split()
        try:
            self.ax.set_yticklabels(yticks)
            self.update_fig()
        except:
             pass
    def update_title(self, change):
        """
        Method for update title graph
        """
        
        self.ax.set_title(change.new)
        self.update_fig()
            
    def update_legend(self, change):
        """
        Method for update legend
        """
        if change.new == True:
            self.ax.legend()
        else:
            self.ax.get_legend().remove()
        self.update_fig()

    def update_grid(self, change):
        """
        Method for update grid (hide or display)
        """
        if change.new == True:
            self.ax.grid(True)
        else:
            self.ax.grid(False)
        self.update_fig()
            
    def update_xlim(self, change):
        """
        Method for update x limits
        """
        
        x_min = change.new[0]
        x_max = change.new[1]
        self.ax.set_xlim(x_min, x_max)
        
        self.update_fig()
            
    def update_ylim(self, change):
        """
        Method for update y limits
        """
        
        y_min = change.new[0]
        y_max = change.new[1]
        self.ax.set_ylim(y_min, y_max)
        self.update_fig()

    def update_fig(self): 
        """
        Method for update fig (main plot)
        """
    
        with self.output:
            clear_output(wait=True)
            display(wg.HBox([self.savebutton,self.fileChooser]),self.fig)

    def manual_update_fig(self,value):   
        """
        Method for update fig manually (in a force way)
        """

        if self.display_legend.value == True:
            self.ax.legend()

        with self.output:
            clear_output(wait=True)
            display(wg.HBox([self.savebutton,self.fileChooser]),self.fig)

    def savefig(self, value):
        """
        Method for save figure in selected fileChooser path
        """
        try:
            self.fig.savefig(self.fileChooser.value)
        except:pass


class Line_graph(Postprocessing_graphs):
    """
    Line Graph object
    """

    def __init__(self,data,ax,label, graph_type,color ='#FF00DD' ):
        super().__init__(graph_type)
        #initial data
        self.x = data[0]
        self.y = data[1]

        self.ax = ax
        self.initial_color = color
        self.label = label
        self.graph_type = graph_type
        
    def add_new_line(self):
        """
        Method for adding a new line to main plot (fig)
        """
        self.line, = self.ax.plot(self.x, self.y, self.initial_color, label = self.label)
        #Color picker widget 
        self.color_picker = widgets.ColorPicker(
            value=self.initial_color, 
            description='', layout=Layout(height='15px', width='100px')
        )

        #Label widget
        self.label_option = widgets.Text(value=self.label,
                         placeholder='Label',description='',
                        disabled=False,layout=Layout(height='15px', width='85px'))
        #Line style widget
        self.linestyle_option =  create_dropdown_input(['solid line','dashed line','dash-dotted line',
        'dotted line'],'solid line','', width = '85px')

        #Observe widgets to update when change
        self.linestyle_option.observe(self.linestyle, 'value')                                                  
        self.color_picker.observe(self.line_color, 'value')
        self.label_option.observe(self.label_text, 'value')
        return self.line, self.color_picker, self.label_option, self.linestyle_option
    
    def line_color(self, change):
        """
        Method for update line color
        """
        self.line.set_color(change.new) 

    def label_text(self, change):
        """
        Method for update label
        """
        self.line.set_label(change.new) 

    def linestyle(self, change):
        """
        Method for update line style
        """
        linestyles={'solid line':'solid','dashed line':'dashed','dash-dotted line':'dashdot',
        'dotted line':'dotted'}
        self.line.set_linestyle(linestyles[change.new]) 


  

import ipywidgets as wg
import sys
import nglview
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
import pandas as pd


class Viewer(object):

    """
    Viewer for GUI app
    """

    
    def __init__(self,filename,cif_content, args, format_ = 'cif'):
        """
        inputs:
            filename: name of uploaded file
            cif_content: cif content of uploaded file
            args: list of additional arguments
                args[0]: Flag, show labels (True or False)
                args[1]: String, selected viewer style
                args[2]: string, selected visualizer
        """

        self.format_ = format_
        if self.format_ != 'extxyz':
            self.filename = "temp.cif"
        else:
             self.filename = "temp.xyz"
        self.cif = cif_content
        
 
        self.periodic_table = pd.read_excel('./complements/periodic_table.xlsx')        
        self.show_structure_labels = args[0]    
        self.viewer_style = args[1]   
        self.visualizer = args[2]                                                                                      

    def build_viewer(self):
        """
        Method for create viewer structure
        """

        if self.visualizer == 'Py3Dmol':
            viewer = py3Dmol.view()
            if self.format_ == 'cif':
                viewer.addModel(self.cif,'cif',{'doAssembly':True,'duplicateAssemblyAtoms':True})
            else:
                viewer.addModel(self.cif,'xyz',{'doAssembly':True,'duplicateAssemblyAtoms':True})

            if self.viewer_style == 'sphere':
                viewer.setStyle({self.viewer_style:{'colorscheme':'Jmol','scale':.5}, 'stick':{'colorscheme':'Jmol'}})
            else:
                viewer.setStyle({self.viewer_style:{'colorscheme':'Jmol','scale':.5}})
            viewer.addUnitCell()
            viewer.zoomTo()

            if self.show_structure_labels:
                viewer.addPropertyLabels("index",{'not':{'elem':'J'}}, {'fontColor':'black','font': 'sans-serif', 'fontSize': 15, 'showBackground':'false','alignment':'center'})
            
            html = viewer._make_html()
        else:
    
            self._get_temporal_file()
            cif = read(self.filename)

            html = nglview.show_ase(cif, gui = True)
            html.center()
            
        self._get_number_of_atoms()
        self._get_number_of_species()
        self._get_chemical_species_label()
        self._get_unit_cell()
        self._get_atomic_coordinates()
        
        return html
    
    def _remove_file(self):
        """
        Method for remove temporal cif file created
        """
        
        if self.format_ != 'extxyz':
            os.system("chmod 777 temp.cif")
        else:
            os.system("chmod 777 temp.xyz")
       
        os.remove(self.filename)
        
    def _get_temporal_file(self):
        """
        Method for create temporal cif file
        """
        with open(self.filename, "w") as f1:
            for line in self.cif:
                f1.write(line)
        
    def _get_number_of_atoms(self):
        """
        Method for get number of atoms
        """
        
        self._get_temporal_file()
        cif = read(self.filename)
        self.number_of_atoms = cif.get_global_number_of_atoms()
        self._remove_file()
        
        
    def _get_number_of_species(self):
        """
        Method for get number of species
        """
        
        self._get_temporal_file()    
        cif = read(self.filename)
        self.number_of_species = len(set(cif.get_chemical_symbols()))

    
    def _get_unit_cell(self):
        """
        Method for get unit cell
        """
        
        self._get_temporal_file()
        cif = read(self.filename)
        self.unit_cell = ''
        for i in cif.cell.array:
            for j in i:
                self.unit_cell+= '    '+ str("%.12f"%j) 
            self.unit_cell+='\n'
        self._remove_file()

    def _get_chemical_species_label(self):
        """
        Method for get chemical species label
        """
        
        self._get_temporal_file()
        cif = read(self.filename)
        numbers = cif.arrays['numbers']
        self.chemical_species = {}
        self.coordinates_ids = {}
        count = 0
        for item,number in enumerate(numbers):
            if number not in self.chemical_species.keys():
                count += 1
                self.coordinates_ids[number] = count
                self.chemical_species[number] = {}
                self.chemical_species[number]['id'] = count
                self.chemical_species[number]['atomic_number'] = number
                self.chemical_species[number]['symbol'] = self.periodic_table[self.periodic_table['atomic_number']==number].symbol.values[0]
  
        self._remove_file()
        
    def _get_atomic_coordinates(self):
        """
        Method for get atomic coordinates
        """
        self._get_temporal_file()
        cif = read(self.filename)
        self.coordinates = ''
        for pos,i in enumerate(cif.arrays['positions']):
            for j in i:
                self.coordinates+= '    '+ str("%.12f"%j) 
            self.coordinates+='    '+str(self.coordinates_ids[cif.arrays['numbers'][pos]])+'\n'
        self._remove_file()
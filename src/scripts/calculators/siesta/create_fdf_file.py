# Defines templates
parallel_options = """#####################################
#       Parallel options            #
#    BlockSize   32                 #
#    ProcessorY  16                 #
#   DiagMemory  1                   #
#   TryMemoryIncrease  true         #
#   DiagScale   1                   #
#   ParallelOverK      false        #
#######################################\n"""

general_options = """#####################################
#       General Options             #
  UseSaveData	     False	        #	
  WriteCoorInital    false          #
  WriteCoorStep      false          #
  WriteForces	     false          #
  WriteKpoints	     false          #
  WriteEigenvalues   false          #
  WriteKbands	     false          #
  WriteBands	     false          #
  WriteMullikenPop   0              #
  WriteDM	     true               #
  WriteCoorXmol	    false           #
  WriteCoorCerius    false          #
  WriteMDXmol	     false          #
  WriteMDhistory     false          #
  WriteDenchar       true           #
  SaveRHO            true           #
#  WriteDenchar       true          #
#  COOP.Write         true          #
#  SaveElectrostaticPotential true  #
#####################################
"""


    
system_configuration = """SystemLabel      {}
SystemName       {}
NumberOfAtoms    {}
NumberOfSpecies  {}
%block ChemicalSpeciesLabel
{}%endblock ChemicalSpeciesLabel
"""


dft_calcule_control = """#---------------------------------------------------------------------------------------------
#    DFT   Calcule Control
#---------------------------------------------------------------------------------------------
PAO.BasisSize     {}
PAO.EnergyShift   {} {}
XC.functional     {}
XC.authors        {}
MaxSCFIterations  {}
MeshCutoff        {}  {}
DM.MixingWeight   {}
DM.Tolerance      1.000E-4
DM.NumberPulay    {}
SolutionMethod    diagon
ElectronicTemperature  25.86  meV
Spin    {}
{}
"""

system_control = """#----------------------------------------------------------------------------------------------
#     System  Relax Control
#----------------------------------------------------------------------------------------------
MD.TypeOfRun     {}
MD.NumCGsteps     {}
MD.VariableCell   {}
MD.MaxForceTol    {} eV/Ang
"""
geometric_structure = """#-----------------------------------------------------------------------------------------------
# Geometric  structure
#-----------------------------------------------------------------------------------------------
{}
AtomicCoordinatesFormat Ang
AtomCoorFormatOut  Ang
{}

AtomicCoordinatesFormat NotScaledCartesianAng
%block AtomicCoordinatesAndAtomicSpecies
{}%endblock AtomicCoordinatesAndAtomicSpecies
"""

def get_fdf_file(args,path, download = False, IS_MOLECULE = False):
    """
    create fdf file given params
    """

    siesta_label = args[0]
    number_of_atoms = args[1]
    number_of_species = args[2]

    basis_set = args[3]
    energy_shift = args[4]
    energy_shift_unit = args[5]
    xc_functional = args[6]
    xc_author = args[7]
    max_scf_iterations = args[8]
    mesh_cutoff = args[9]
    mesh_cuttoff_units = args[10]
    mixing_weight =  args[11]
    number_pulay = args[12]
    k1= args[13]
    k2= args[14]
    k3= args[15]
    k4= args[16]
    k5= args[17]
    k6= args[18]
    k7= args[19]
    k8= args[20]
    k9= args[21]

    cg_steps= args[22]
    variable_cell= args[23]

    unit_cell= args[24]

    if unit_cell == False or IS_MOLECULE == True:
      unit_cell = ''
    else:
      unit_cell = f"""%block LatticeVectors
{args[24]}%endblock LatticeVectors"""
    coordinates= args[25]
    spin = args[26]
    chemical_species_label= args[27]
    type_run = args[28]
    max_force_tol = args[29]

    if IS_MOLECULE == False:
      k_block = f"""%block kgrid_Monkhorst_Pack
        {k1}    {k2}     {k3}     0.0
        {k4}    {k5}     {k6}     0.0
        {k7}    {k8}     {k9}     0.0
%endblock kgrid_Monkhorst_Pack"""
      lattice_constant = "LatticeConstant 1.0 Ang"
    else:
      k_block = ""
      lattice_constant = ""

    
    fdf_file = parallel_options + general_options + system_configuration.format(siesta_label,
     siesta_label,number_of_atoms, number_of_species,chemical_species_label) + dft_calcule_control.format(basis_set,
     energy_shift, energy_shift_unit, xc_functional, xc_author,max_scf_iterations, mesh_cutoff,
     mesh_cuttoff_units,mixing_weight,number_pulay,spin,k_block) + system_control.format(type_run,cg_steps,
     variable_cell, max_force_tol)+geometric_structure.format(lattice_constant,unit_cell,coordinates)
    if download is False:
      with open(path, "w") as text_file:
          text_file.write(fdf_file)
    else:
      return fdf_file
    
def change_system_relax_control(input_path, output_path, name, numCGsteps='0', variableCell='.False.', typeRun = None,  max_force_tol= None):
  """
  Function to change MD.NumCGsteps and MD.VariableCell values in fdf file
  """
  f = open(input_path, "r")
  string = ''
  for i in f:
      if 'MD.NumCGsteps' in i:
        string+= f'MD.NumCGsteps  {numCGsteps}\n'
      elif 'MD.VariableCell' in i:
        string += f'MD.VariableCell  {variableCell}\n'
      elif typeRun and 'MD.TypeOfRun' in i:
        string += f"MD.TypeOfRun {typeRun}\n"
      elif max_force_tol and 'MD.MaxForceTol' in i:
        string += f"MD.MaxForceTol    {max_force_tol} eV/Ang\n"
      else:
        string += i
  f.close()

  #open text file
  text_file = open(f"{output_path}/{name}.fdf", "w")
  
  #write string to file
  text_file.write(string)
  return
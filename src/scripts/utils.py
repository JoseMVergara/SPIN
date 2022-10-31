import numpy as np
import matplotlib.pylab as plt
from math import atan2, cos, sin, acos, gcd
import ipywidgets as wg
from pathlib import Path
from re import compile
from ase import Atoms
from ase.utils import reader
from ase.units import Bohr
from threading import Thread

def get_calculator_logo(path, format_):
    """
    get logo of each calculator available
    """
    file = open(path, "rb")
    image = file.read()
    image_wg=wg.Image(
        value=image,
        format=format_,
        width=100,
        height=200,
    )
    return image_wg

def get_check_image(path, format_):
    """
    get logo of each calculator available
    """
    file = open(path, "rb")
    image = file.read()
    image_wg=wg.Image(
        value=image,
        format=format_,
        width=30,
        height=30,
    )
    return image_wg


class PropagatingThread(Thread):
    """
    Catch thread exceptions when "stop" button is pressed
    """
    def run(self):
        self.exc = None
        try:
            if hasattr(self, '_Thread__target'):
                # Thread uses name mangling prior to Python 3.
                self.ret = self._Thread__target(*self._Thread__args, **self._Thread__kwargs)
            else:
                self.ret = self._target(*self._args, **self._kwargs)
        except BaseException as e:
            self.exc = e
            raise self.exc

    def join(self, timeout=None):
        super(PropagatingThread, self).join(timeout)
        if self.exc:
            raise self.exc
        return self.ret

def get_siesta_examples():
    """
    initialize siesta examples
    """
    
    bandLinesExample = """#example:
%block BandLines
1 1.000 1.000 1.000 L # Begin at L
20 0.000 0.000 0.000 \Gamma # 20 points from L to gamma
25 2.000 0.000 0.000 X # 25 points from gamma to X
30 2.000 2.000 2.000 \Gamma # 30 points from X to gamma
%endblock BandLines"""

    
    bandPointsExample = """#example:
%block BandPoints
0.000 0.000 0.000 # This is a comment. eg this is gamma
1.000 0.000 0.000
0.500 0.500 0.500
%endblock BandPoints"""

    
    exampleProjectedDensityOfStates = """#example:
%block ProjectedDensityOfStates
-20.00 10.00 0.200 500 eV
%endblock ProjectedDensityOfStates"""

    
    examplePDOSkgrid_Monkhorst_Pack = """#example:
%block PDOS.kgrid_Monkhorst_Pack
1       0       0       0.0
0       1       0       0.0
0       0      101      0.0
%endblock PDOS.kgrid_Monkhorst_Pack"""

    
    exampleOpticalMesh = """#example:
%block Optical.Mesh
5 5 5
%endblock Optical.Mesh"""

    exampleOpticalVector = """#example:
%block Optical.Vector
1.0 0.0 0.5
%endblock Optical.Vector"""   

    return bandLinesExample, bandPointsExample, exampleProjectedDensityOfStates,examplePDOSkgrid_Monkhorst_Pack, exampleOpticalMesh,exampleOpticalVector


def mag(x,y):
    """
    Return the magnitude of a vector
    input:
        x: x-coordinates 
        y: y-coordinates
    output:
        rst: Magnitude of a vector 
    """
    rst = 0
    for i in range(len(x)):
        rst += (y[i] - x[i])**2
    return rst**0.5

def Nanotube(app,a1,a2,atom1,atom2,n,m,verbose,maxiter=300):
    """
    Create blue-phosphorene nanotube (with buckling)
    inputs:
        app: bond lenght between phosphorene atoms (Ångströms).
        a1 and a2: unit cell vectors.
        atom1 and atom2: initial atoms to be replicated.
        n and  m: integers to defined the nanotube quirality.
        verbose: if True, plots and prints are displayed.
    outputs:
        AtomsR: nanotube coordinates
        unit_cell: nanotube unit cell
    """
    
    #convert inital values to arrays
    a1 = np.array(a1)
    a2 = np.array(a2)
    atom1 = np.array(atom1)
    atom2 = np.array(atom2)

    #we define the chiral vector along which the nanotube is coiled
    Ch = n*a1 + m*a2

    #defined the gcd between (2m+n, 2n+m). 
    num1 = 2*m+n
    num2 = 2*n+m
    dR = gcd(num1,num2)

    #defined the traslation vector of nanotube
    T = num1/dR * a1 -num2/dR * a2

    #The parallelepiped formed by the chiral vector and the 
    #translational vector is the unit cell of the nanotube.
    unitCell = np.linalg.norm(T)

    #replicated atoms to get monolayer
    h = np.linalg.norm(a1)
    Ca = np.cos(30*np.pi/180.)*h
    Co = np.sin(30*np.pi/180.)*h
    atomsx = [atom1[0],atom2[0]]
    atomsy = [atom1[1],atom2[1]]
    atomsz = [atom1[2],atom2[2]]

    for i in range(-maxiter,maxiter):
        mult = np.arange(-i,i+2,2)
        mult = mult[::-1]
        for j in mult:
            x1,y1 = a1[0]+Ca*i,a1[1]+Co*j
            x2,y2 = a2[0]+Ca*i,a2[1]+Co*j

            atomsx += [x1,x2+app]
            atomsy += [y1,y2]
            atomsz += [atom1[2],atom2[2]]

    #get the chiral angle and reciprocal network vectors
    ChT = np.array([maxiter*2.,0])
    TT = np.array([maxiter*2.,0])
    ChI = ChT-Ch
    TI = TT-T
    Chc = np.array([ChI,ChT])
    Tc = np.array([TI,TT])
    Chc2 = np.array([ChI-T,ChT-T])
    Tc2 = np.array([TI-Ch,TT-Ch])

    f = 0.1

    for i in range(len(Chc)):
        if i == 0:
            Chc[i][0],Chc[i][1] = Chc[i][0]-f,Chc[i][1]+f
            Chc2[i][0],Chc2[i][1] = Chc2[i][0]+f,Chc2[i][1]-f
            Tc[i][0],Tc[i][1] = Tc[i][0]+f,Tc[i][1]-f
            Tc2[i][0],Tc2[i][1] = Tc2[i][0]+f,Tc2[i][1]-f
        else:
            Chc[i][0],Chc[i][1] = Chc[i][0]-f,Chc[i][1]+f
            Chc2[i][0],Chc2[i][1] = Chc2[i][0]+f,Chc2[i][1]-f
            Tc[i][0],Tc[i][1] = Tc[i][0]-f,Tc[i][1]+f
            Tc2[i][0],Tc2[i][1] = Tc2[i][0]-f,Tc2[i][1]+f
        
    points = [Chc,Tc,Chc2,Tc2]
    Atoms = []
    angles = []

    vCh = Chc[1]-Chc[0]
    mCh = mag(Chc[0],Chc[1])
    R = mCh/(2*np.pi)

    xCh = vCh[0]
    yCh = vCh[1]
    #anCh = atan2(yCh,xCh)

    for i in range(len(atomsx)):
        point = [atomsx[i],atomsy[i]]
        for j in points:
            l1 = mag(point,j[0])
            l2 = mag(point,j[1]) 
            l3 = mag(j[0],j[1])
            try:
                angle = acos((l3**2 - l1**2 - l2**2)/(-2.*l1*l2))
                angles += [angle]
            except:pass
  
        if abs(np.sum(angles)-2.*np.pi)<1e-3: 
            Atoms += [[atomsx[i],atomsy[i],atomsz[i]]]
        angles = [] 
    
    #cut atoms to be rolled up
    Atoms = np.array(Atoms)


    for i in range(len(Atoms-1)):
        Atoms[i][0] = np.array(Atoms[i][0]) - Chc2[0][0]
        Atoms[i][1] = np.array(Atoms[i][1]) - Chc2[0][1]
        
    if verbose == True:
        xs = []
        ys = []
        for i in range(len(Chc)):
            xs += [Chc[i][0],Tc[i][0],Chc2[i][0],Tc2[i][0]]
            ys += [Chc[i][1],Tc[i][1],Chc2[i][1],Tc2[i][1]]
            
        
        for i in range(len(xs)):
            xs[i] += -Chc2[0][0]
            ys[i] += -Chc2[0][1]
            
        x = []
        y = []
        
        for i in range(len(Atoms)):
            x += [Atoms[i][0]]
            y += [Atoms[i][1]]
            
        plt.plot(xs,ys,'ob')
        plt.plot(x,y, '.r')
        plt.grid()
        plt.show()
    
    #roll up atoms to get nanotube
    AtomsR = []
    thetaCh = None
    thetap = None

    vCh = Chc[1]-Chc[0]
    mCh = mag(Chc[0],Chc[1])
    #RCh = mCh/(2*np.pi)

    xCh = vCh[0]
    yCh = vCh[1]
    thetaCh = atan2(yCh,xCh)
    
    zmax = max(atom1[2],atom2[2])
    #zmin = min(atom1[2],atom2[2])
    #zdif = zmax-zmin
    for i in range(len(Atoms)):
        xp = Atoms[i][0]
        yp = Atoms[i][1]
        zp = Atoms[i][2]
    
            
        thetap = atan2(yp,xp)
        theta = thetap - thetaCh

        Vpoint = [xp,yp]
        Vmag = np.linalg.norm(Vpoint)
        alpha = Vmag*cos(theta)/R
        
        #we are taking into consideration the buckling value
        if abs(zp-zmax) < 1e-2:
            
            z = Vmag*sin(theta) 
            x = R*cos(alpha)
            y = R*sin(alpha)
        else:
            
            z = Vmag*sin(theta) 
            x = 0.9*R*cos(alpha)
            y = 0.9*R*sin(alpha)
        
        if z>np.linalg.norm(T):
            pass
        else: 
            AtomsR += [(x-25,y-25,z)]

    #define unit cell 
    unit_Cell = [[ 0.0000000000,    0.0000000000,    0.0000000000],
                [0.0000000000,   0.0000000000,    0.0000000000],
                [0.0000000000,    0.0000000000, unitCell]]

    return AtomsR,unit_Cell




"""

template of run job file in pbs cluster with ssh connection
"""

def get_ssh_pbs_run(args):

    number_processor = args[0]
    node_name = args[1]
    omp_num_threads = args[2]
    siesta_path= args[3]
    pseudopotentials_path= args[4]
    python_path= args[5]
    #folder = args[6]

    template_file = """#!/bin/bash
# Example for {} processors
#$ -N PREFIX
#$ -pe smp {}
#$ -l h={}
#$ -cwd
#
export OMP_NUM_THREADS={}
export ASE_SIESTA_COMMAND="mpirun -np {} {} < PREFIX.fdf > PREFIX.out"
export SIESTA_PP_PATH={}

{} run.py"""

    ssh_pbs_run = template_file.format(number_processor,number_processor, node_name, 
    omp_num_threads, number_processor, siesta_path, pseudopotentials_path, python_path)

    return ssh_pbs_run


"""
Template of python file to run calculations in cluster 
"""

def get_ssh_python_run(args):

    positions = args[0]
    ase_command = args[1]
    calculator = args[2]


    template_ssh_python = """from ase import Atoms
from ase.calculators.siesta import Siesta
from ase.build import bulk, molecule, graphene_nanoribbon
from ase.build import nanotube as ase_nanotube
from numpy import array
structure = {}
structure.set_positions({})
        """

    ssh_python = template_ssh_python.format(ase_command, positions)
    ssh_python = ssh_python + """
calculator = """ + calculator + """ 
structure.calc = calculator
structure.get_total_energy()"""

   
    return ssh_python.replace("'",'"')


"""Helper functions for read_fdf."""

_label_strip_re = compile(r'[\s._-]')


def _labelize(raw_label):
    # Labels are case insensitive and -_. should be ignored, lower and strip it
    return _label_strip_re.sub('', raw_label).lower()


def _is_block(val):
    # Tell whether value is a block-value or an ordinary value.
    # A block is represented as a list of lists of strings,
    # and a ordinary value is represented as a list of strings
    if isinstance(val, list) and \
       len(val) > 0 and \
       isinstance(val[0], list):
        return True
    return False


def _get_stripped_lines(fd):
    # Remove comments, leading blanks, and empty lines
    return [_f for _f in [L.split('#')[0].strip() for L in fd] if _f]


@reader
def _read_fdf_lines(file):
    # Read lines and resolve includes
    lbz = _labelize

    lines = []
    for L in _get_stripped_lines(file):
        w0 = lbz(L.split(None, 1)[0])

        if w0 == '%include':
            # Include the contents of fname
            fname = L.split(None, 1)[1].strip()
            parent_fname = getattr(file, 'name', None)
            if isinstance(parent_fname, str):
                fname = Path(parent_fname).parent / fname
            lines += _read_fdf_lines(fname)

        elif '<' in L:
            L, fname = L.split('<', 1)
            w = L.split()
            fname = fname.strip()

            if w0 == '%block':
                # "%block label < filename" means that the block contents
                # should be read from filename
                if len(w) != 2:
                    raise IOError('Bad %%block-statement "%s < %s"' %
                                  (L, fname))
                label = lbz(w[1])
                lines.append('%%block %s' % label)
                lines += _get_stripped_lines(open(fname))
                lines.append('%%endblock %s' % label)
            else:
                # "label < filename.fdf" means that the label
                # (_only_ that label) is to be resolved from filename.fdf
                label = lbz(w[0])
                fdf = read_fdf(fname)
                if label in fdf:
                    if _is_block(fdf[label]):
                        lines.append('%%block %s' % label)
                        lines += [' '.join(x) for x in fdf[label]]
                        lines.append('%%endblock %s' % label)
                    else:
                        lines.append('%s %s' % (label, ' '.join(fdf[label])))
                # else:
                #    label unresolved!
                #    One should possibly issue a warning about this!
        else:
            # Simple include line L
            lines.append(L)
    return lines


def read_fdf(fname):
    """Read a siesta style fdf-file.
    The data is returned as a dictionary
    ( label:value ).
    All labels are converted to lower case characters and
    are stripped of any '-', '_', or '.'.
    Ordinary values are stored as a list of strings (splitted on WS),
    and block values are stored as list of lists of strings
    (splitted per line, and on WS).
    If a label occurres more than once, the first occurrence
    takes precedence.
    The implementation applies no intelligence, and does not
    "understand" the data or the concept of units etc.
    Values are never parsed in any way, just stored as
    split strings.
    The implementation tries to comply with the fdf-format
    specification as presented in the siesta 2.0.2 manual.
    An fdf-dictionary could e.g. look like this::
        {'atomiccoordinatesandatomicspecies': [
              ['4.9999998', '5.7632392', '5.6095972', '1'],
              ['5.0000000', '6.5518100', '4.9929091', '2'],
              ['5.0000000', '4.9746683', '4.9929095', '2']],
         'atomiccoordinatesformat': ['Ang'],
         'chemicalspecieslabel': [['1', '8', 'O'],
                                  ['2', '1', 'H']],
         'dmmixingweight': ['0.1'],
         'dmnumberpulay': ['5'],
         'dmusesavedm': ['True'],
         'latticeconstant': ['1.000000', 'Ang'],
         'latticevectors': [
              ['10.00000000', '0.00000000', '0.00000000'],
              ['0.00000000', '11.52647800', '0.00000000'],
              ['0.00000000', '0.00000000', '10.59630900']],
         'maxscfiterations': ['120'],
         'meshcutoff': ['2721.139566', 'eV'],
         'numberofatoms': ['3'],
         'numberofspecies': ['2'],
         'paobasissize': ['dz'],
         'solutionmethod': ['diagon'],
         'systemlabel': ['H2O'],
         'wavefunckpoints': [['0.0', '0.0', '0.0']],
         'writedenchar': ['T'],
         'xcauthors': ['PBE'],
         'xcfunctional': ['GGA']}
    """
    fdf = {}
    lbz = _labelize
    lines = _read_fdf_lines(fname)
    while lines:
        w = lines.pop(0).split(None, 1)
        if lbz(w[0]) == '%block':
            # Block value
            if len(w) == 2:
                label = lbz(w[1])
                content = []
                while True:
                    if len(lines) == 0:
                        raise IOError('Unexpected EOF reached in %s, '
                                      'un-ended block %s' % (fname, label))
                    w = lines.pop(0).split()
                    if lbz(w[0]) == '%endblock':
                        break
                    content.append(w)

                if label not in fdf:
                    # Only first appearance of label is to be used
                    fdf[label] = content
            else:
                raise IOError('%%block statement without label')
        else:
            # Ordinary value
            label = lbz(w[0])
            if len(w) == 1:
                # Siesta interpret blanks as True for logical variables
                fdf[label] = []
            else:
                fdf[label] = w[1].split()
    return fdf
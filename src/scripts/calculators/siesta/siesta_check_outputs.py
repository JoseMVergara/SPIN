#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Author: @JMVergara

'''
Program that extract necessary information and automatize processes from SIESTA software outputs.

The program do the process for all files with ".out" extension that are in the same folder of the script

How to use:
	$python Check.py

outputs:
Note: here 'File' is the name of the ".out" file. example: if BPNT1414.out File = BPNT1414

	Status.csv file: Show iformation about the correct termination of job in SIESTA: OK, FAILED
	File-R.fdf file: replace relaxed coordinates in fdf input file for future calculations  
'''

import os
import pandas as pd
from itertools import dropwhile,takewhile   
 
def NotBlockStartRelaxedCoor(line):
    return not line.startswith('outcoor: Relaxed atomic coordinates (Ang):')

def NotBlockStartUnRelaxedCoor(line):
    return not line.startswith('outcoor: Final (unrelaxed) atomic coordinates (Ang):')

def NotBlockStartUnitCell(line):
	return not line.startswith('outcell: Unit cell vectors (Ang):')

def NotBlockEnd(line):
	return not line.startswith('\n')

def NotBlockStartInputCoordinates(line):
	return not line.startswith('%block AtomicCoordinatesAndAtomicSpecies')

def NotBlockInputEnd(line):
	return not line.startswith('%endblock AtomicCoordinatesAndAtomicSpecies')

def NotBlockStartInputUnit(line):
	return not line.startswith('%block LatticeVectors')

def NotBlockUnitCellEnd(line):
	return not line.startswith('%endblock LatticeVectors')

def WriteLines(lines,outputFile):
	'''Write lines in outputfile
	input: 
		lines = lines of text or numbers that you want to write
		outputFile = file where the lines will be written'''

	for line in lines:
		outputFile.write(line)

def WriteFindCoorFile(data):
	"""
	Function to write coordinates in temporal file
	"""
	blockInputCoor = dropwhile(NotBlockStartInputCoordinates, data)
	next(blockInputCoor)
	block=list(takewhile(NotBlockInputEnd, blockInputCoor))
	findFile = open('findCoor.txt',"w")
	WriteLines(block,findFile)

def WriteFindUnitFile(data):   
	"""
	Function to write unit cell in temporal file
	"""
	blockUnitCell = dropwhile(NotBlockStartInputUnit, data)
	next(blockUnitCell)
	block=list(takewhile(NotBlockUnitCellEnd, blockUnitCell))
	findFile = open('findUnit.txt',"w")
	WriteLines(block,findFile)
	

def GetInputInfo(inputFile,valueToObtain):
	"""
	funciton to search given value in file and obtain corresponding line
	"""
	with open(inputFile) as inputFile:
		inputLines = iter(inputFile)
		while True:
			inputLine = next(inputLines)
			if inputLine.startswith(valueToObtain,0,50):
				result = inputLine
				break
	return inputLine

def GetRelaxedCoor(file):
	"""
	Function to obtain the relaxed coordinates
	"""
	blockRelaxedCoor = dropwhile(NotBlockStartRelaxedCoor, file)
	next(blockRelaxedCoor)
	block=list(takewhile(NotBlockEnd, blockRelaxedCoor))
	replaceFile = open('replaceCoor.txt',"w")
	WriteLines(block,replaceFile)

def GetUnRelaxedCoor(file):
	"""
	Function to obtain the unrelaxed coordinates
	"""
	blockUnRelaxedCoor = dropwhile(NotBlockStartUnRelaxedCoor, file)
	next(blockUnRelaxedCoor)
	block=list(takewhile(NotBlockEnd, blockUnRelaxedCoor))
	replaceFile = open('replaceCoor.txt',"w")
	WriteLines(block,replaceFile)

def GetUnitCell(file):
	"""
	Function to obtain unit cell
	"""

	blockUnitCell = dropwhile(NotBlockStartUnitCell, file)
	next(blockUnitCell)
	block=list(takewhile(NotBlockEnd, blockUnitCell))
	try:
		GetUnitCell(file)
	except: 
		replaceFile = open('replaceUnit.txt',"w")
		WriteLines(block,replaceFile)

def CheckScf(inputFile,file):
	"""
	Function to check if calculation ends successfully acording to scf iterations
	"""
	lines = iter(file)
	while True:
		try:
			line = next(lines)
		except:
			break
		if line.startswith( 'SCF cycle converged after', 0, 30 ):
			lastScf = line
	splitLastScf = lastScf.split(' ')
	#print(splitLastScf)
	try:
		numIter = int(splitLastScf[7])
	except:
		numIter = int(splitLastScf[-2])
	maxSCF = GetInputInfo(inputFile,'MaxSCFIterations').split('\n')
	maxSCF = int(maxSCF[0].split(' ')[-1])

	if numIter < maxSCF:
		status = True
	else:
		status = False

	return status

def CheckMaxForceTol(inputFile,file):
	"""
	Function to check if calculations ends successfully according to max force tolerance value
	"""
	lines = iter(file)
	while True:
		try:
			line = next(lines)
		except:
			break
		if line.startswith( 'Max', 3, 30 ):
			lastForceTol = line
	splitLastForceTol = lastForceTol.split(' ')
	numLastForceTol = float(splitLastForceTol[7])

	maxForceTol = GetInputInfo(inputFile,'MD.MaxForceTol').split('\n')

	maxForceTol = float(maxForceTol[0].split(' ')[4].replace(',','.'))

	if numLastForceTol < maxForceTol:
		status =  True 
	else:
		status = False

	return status

def CheckCGsteps(inputFile,file):
	"""
	Function to check if calculations ends successfully according to num CG steps value
	"""
	lines = iter(file)
	while True:
		try:
			line = next(lines)
		except:
			break
		if line.startswith( 'Begin FIRE opt. move', 24, 50 ):
			lastCGsteps = line
	splitCGsteps = lastCGsteps.split('\n')
	splitCGsteps = splitCGsteps[0].split(' ')
	numCGsteps = int(splitCGsteps[-1])

	maxCGsteps = GetInputInfo(inputFile,'MD.NumCGsteps').split('\n')
	maxCGsteps = int(maxCGsteps[0].split(' ')[-1])

	if numCGsteps < maxCGsteps:
		status = True
	else:
		status = True

	return status

def GetTotalEnergy(file):
	"""
	Functions to get the total energy value
	"""
	lines = iter(file)
	
	while True:
		try:
			line = next(lines)
		except:
			break
		if line.startswith( 'Total', 16, 50 ):
			totalEnergy = line
	try:
		split = totalEnergy.split('\n')
		totalEnergy = float(split[0].split(' ')[-1])
	except:
		totalEnergy = 0
	return totalEnergy

def FindFilenames(pathToDir, suffix=".out"): 
	"""Search all the .out files found in a directory.

		input: Folder where the files are
		output: list of all files found out"""
	filenames = os.listdir(pathToDir) 
	filenamesWithSuffix = [ filename for filename in filenames if filename.endswith(suffix) ] 	
	filenamesWithoutSuffix = []
	for filename in filenamesWithSuffix:
		filenamesWithoutSuffix += [filename.split('.out')[0]]

	return filenamesWithoutSuffix

def ReplaceCoordinatesInputFile(inputFile,RelaxFile):
	"""
	Function to replace initial coordinates in relaxed coordinates
	"""
	RelaxFiletemp = RelaxFile+"tmp.fdf"
	RelaxFile = RelaxFile+"-R.fdf"
	with open(inputFile) as data:
		WriteFindCoorFile(data)
	try:
		with open(inputFile) as data:
			WriteFindUnitFile(data)
		UNIT_CELL = True
	except:
		UNIT_CELL = False

	findLines = open('findCoor.txt').read().split('\n')
	replaceLines = open('replaceCoor.txt').read().split('\n')
	findReplace = dict(zip(findLines, replaceLines))
	with open(inputFile) as data:
		with open(RelaxFiletemp,"w") as newData:
			for line in data:
				for key in findReplace:
					if key in line:
						line = line.replace(key, findReplace[key])
				newData.write(line)
	if UNIT_CELL == True:
		findLines = open('findUnit.txt').read().split('\n')
		replaceLines = open('replaceUnit.txt').read().split('\n')
		findReplace = dict(zip(findLines, replaceLines))
		with open(RelaxFiletemp) as data:
			with open(RelaxFile,"w") as newData:
				for line in data:
					for key in findReplace:
						if key in line:
							line = line.replace(key, findReplace[key])
					newData.write(line)
	return UNIT_CELL
	
def GetTypeFile(outputfile):
	"""
	Function to get type of calculation, relaxed or result (different to relaxed)
	"""
	maxCGsteps = GetInputInfo(outputfile,'MD.NumCGsteps').split('\n')
	maxCGsteps = int(maxCGsteps[0].split(' ')[-1])
	if maxCGsteps > 0:
		typeFile = 'Relaxed'
	else:
		typeFile = 'Result'
	return typeFile

def CheckIfFileIsEmpty(path):
	"""
	Function to check if file is empty
	"""
	try:
		return os.stat(path).st_size==0
	except:
		return True

def GetFerminEnergy(file):
	"""
	Function to get fermin energy value if exist
	"""
	lines = iter(file)
	while True:
		try:
			line = next(lines)
		except:
			break
		if line.startswith( '   scf:', 0, 7 ):
			lastScf = line
	splitLastScf = lastScf.split(' ')[-1]
	ferminEnergy = splitLastScf.split('\n')[0]
	return ferminEnergy
	
def CheckFile(filename):
	"""
	Function to check file, get total energy and fermin energy
	input: 
		filename: path of file to be checked
	output: 
		info: list of attributes, status, total energy, fermin energy etc.
	"""

	#print('=======================================')
	#print('Checking '+filename+'..................')
	outputFile = filename + '.out'
	inputFile = filename + '.fdf' 
	typeFile = GetTypeFile(outputFile)
	if typeFile == 'Relaxed':		
		with open(outputFile) as file:
			totalEnergy = GetTotalEnergy(file)
		with open(outputFile) as file:
			statusCG = CheckCGsteps(inputFile,file)
		with open(outputFile) as file:
			statusForceTol = CheckMaxForceTol(inputFile,file)
		with open(outputFile) as file:
			statusScf = CheckScf(inputFile,file)
		with open(outputFile) as file:
			GetUnitCell(file)
		flag = False
		with open(outputFile) as file:
			try:
				GetRelaxedCoor(file)
			except StopIteration:
				flag = True
		if flag == True:
			with open(outputFile) as file:
				GetUnRelaxedCoor(file)
		statusList = [statusCG,statusForceTol,statusScf]
		if False in statusList:
			fileStatus = 'FAILED'
		else:
			fileStatus = 'OK'
	
		info = [filename,totalEnergy,fileStatus,'NA']
		UNIT_CELL = ReplaceCoordinatesInputFile(inputFile,filename)

		os.chmod('findCoor.txt',0o777)
		os.chmod('replaceCoor.txt',0o777)
		if UNIT_CELL == True:	
			os.chmod('findUnit.txt',0o777)
			os.chmod('replaceUnit.txt',0o777)
		os.chmod(filename+"tmp.fdf",0o777)

		os.remove('findCoor.txt')
		os.remove('replaceCoor.txt')
		if UNIT_CELL == True:	
			os.remove('findUnit.txt')
			os.remove('replaceUnit.txt')
		os.remove(filename+"tmp.fdf")
	else:
		filenameAux = filename.split('-R')[0]
		bandsFile = filenameAux + '.dat'
		DOSFile = filenameAux + '.DOS'
		with open(outputFile) as file:
			totalEnergy = GetTotalEnergy(file)

		with open(outputFile) as file:
			ferminEnergy = GetFerminEnergy(file)
		if CheckIfFileIsEmpty(bandsFile) == False and CheckIfFileIsEmpty(DOSFile) == False:
			fileStatus='OK'
		else:
			fileStatus='False'
			
		info = [filename,totalEnergy,fileStatus,ferminEnergy]

	return info



def check(path):
	"""
	Check all *.out files in path
	output:
		database in csv and xlsx format containig System name, total energy, fermin energy and Status (OK if run ends successfully else FAILED)
	"""
	filenames = FindFilenames(path)
	Data = []
	for filename in filenames:
		#print(filename)
		info = CheckFile(path+filename)
		Data += [info]

	infoFile = pd.DataFrame(Data,columns=['System','TotalEnergy','Status','FerminEnergy'])
	infoFile.to_csv(path+"Status.csv")
	infoFile.to_excel(path+"Status.xlsx")

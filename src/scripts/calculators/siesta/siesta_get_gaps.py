#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Author: @JoseMVergara
import os
import numpy as np
import pandas as pd
from itertools import dropwhile,takewhile 

def NotBlockStart(line):
	return not line.startswith('\n')

def NotBlockEnd(line):
	return not line.startswith('\n')

def GetBands(file):
	try:	
		blockBand = dropwhile(NotBlockStart,file)
		next(blockBand)
		block=list(takewhile(NotBlockEnd,blockBand))
		return block
	except StopIteration:
		return False


def FindFilenames(pathToDir, suffix): 
	"""Search all the ."suffix" files found in a directory.

		input: Folder where the files are
		output: list of all files found """
	filenames = os.listdir(pathToDir) 
	filenamesWithSuffix = [ filename for filename in filenames if filename.endswith(suffix) ] 	
	filenamesWithoutSuffix = []
	for filename in filenamesWithSuffix:
		filenamesWithoutSuffix += [filename.split('.')[0]]

	return filenamesWithSuffix

def GetGap(datFile):
	with open(datFile) as file:
		bands = []
		while True:
			block = GetBands(file)
			if block is not False:
				bands += [block]
			else:
				break
	positiveY = []
	negativeY = []
	for band in bands:
		x = []
		y = []
		for line in band:
			splitBand = np.array(line.split(),dtype=float)
			x+=[splitBand[0]]
			y+=[splitBand[1]]
		y = np.array(y)
		try:
			positiveY += [np.min(y[y>0])]
		except:
			pass
		try:
			negativeY += [np.max(y[y<0])]
		except:
			pass

	gap = np.min(positiveY) + np.abs(np.max(negativeY))

	return gap


files=FindFilenames('./','.dat')
results = []
for filename in files:
	gap = GetGap(filename)
	results+=[[filename,gap]]

infoFile = pd.DataFrame(results,columns=['System','Gap'])
infoFile.to_csv("Gap.csv")
infoFile.to_excel("Gap.xlsx")
	

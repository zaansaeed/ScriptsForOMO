import pandas as pd
import math
from numpy import *
import os
import glob
input = '/Users/zaansaeed/Google Drive/My Drive/OMO Lab - Peptide Cyclization - Zaan Saeed/Data/Peptide Library/R6C7_XYZ' #input to where files .xyz files are stored.
atoms = [72, 44git] # select atom IDs to calculate distances
#22or28? c 56N for R1C1
#75N 59C R8C1
#72N  44C? R6C7
csv_of_qikprop= '/Users/zaansaeed/Downloads/OMO Lab - Peptide Cyclization - Peptide Descriptors.xlsx - PeptideR6C7_descriptors.csv' #csv of qikprop, only names then propeties

#csv file should be the same length as .xyz's in "input"

df = pd.read_csv(csv_of_qikprop)
df = df.iloc[:, 1:38]
###############################################################################
def calculatedistances():
	outputs = []
	#FileName = argv[1]
	#tempAtom = "\t".join(sys.argv[2:])
	#atoms = tempAtom.split()
	file_list = glob.glob(os.path.join(input, '**', '*.xyz'), recursive=True)

	file_list.sort()
	for file_path in file_list:

		FileName = file_path

		# function to get x,y,z coordinates of the atomic positions
		def getAtomPositions(offset):
			file = open(FileName)
			i = 0
			while (i < (offset+1)):
				file.readline()
				i = i + 1
			atomLine = file.readline()

			atomInfo = atomLine.split()
			#print atomInfo[1]
			file.close()
			# parsing string and returning float representation of x, y, z values
			return array([float(atomInfo[1]), float(atomInfo[2]), float(atomInfo[3])])

		### Distances ###
		if (len(atoms) == 2):
			a1 = getAtomPositions(int(atoms[0]))
			a2 = getAtomPositions(int(atoms[1]))
			v1 = a2-a1
			dist = linalg.norm(v1)
			#print(FileName + "\t" + '%.1f' % dist)
			outputs.append(dist)

		### Angles ###
		if (len(atoms) == 3):
			a1 = getAtomPositions(int(atoms[0]))
			a2 = getAtomPositions(int(atoms[1]))
			a3 = getAtomPositions(int(atoms[2]))

			# making appropriate vectors and normals
			v1 = a1-a2
			v2 = a3-a2
			da = arccos(dot(v1,v2)/(linalg.norm(v1)*linalg.norm(v2)))
			print(FileName + "\t" + '%.1f' % degrees(da))
			outputs.append(da)

		### Dihedrals ###
		if (len(atoms) == 4):
			# put positions in array
			a1 = getAtomPositions(int(atoms[0]))
			a2 = getAtomPositions(int(atoms[1]))
			a3 = getAtomPositions(int(atoms[2]))
			a4 = getAtomPositions(int(atoms[3]))

			# making appropriate vectors and normals
			v1 = a2-a1
			v2 = a3-a2
			n1 = cross(v1,v2)
			v3 = a2-a3
			v4 = a4-a3
			n2 = cross(v3,v4)

			# finding dihedral angle
			da = arccos(-dot(n1,n2)/(linalg.norm(n1)*linalg.norm(n2)))

			# checking the sign of the dihedral then displaying result
			if dot(n1,v4) < 0:
				print(FileName + "\t" + '%.1f' % degrees(da))
			else:
				print(FileName + "\t" + '%.1f' % -degrees(da))
			outputs.append(da)

	return outputs



def boltzmann(values, name = None,):
    numerator = 0
    denominator = 0
    for i in range(df.shape[0]):
        e_term = math.exp(-(df.loc[i, "Relative Potential Energy-OPLS-2005"])/(298*8.314 * 10**-3))
        denominator += e_term
        if name:
            numerator += e_term * df.loc[i,name]
        else:
            numerator += e_term * values[i]

    return numerator/denominator




if csv_of_qikprop != '':
	for column_name, column_data in df.items():
		print(f'Boltzmann Average of : {column_name} is ' + str(boltzmann(column_data,column_name)))
print(f'Boltzmann Average of N-C termini distance is : ' + str(boltzmann(calculatedistances())))







#!/usr/bin/env python


#This script is modified from the script read_growth.py by Matthew Jenior 
# read.cellindex.py takes data exported from the Xcelligence RTCA system and writes restructured data to a new file
# USAGE: read_cellindex.py raw_data output_table


import sys


# Function to read in cell index values and create a dictioaries of time series
def read_data(raw_data, well_dict):

	for line in raw_data:

		line = line.strip().split()

		if len(line) < 3:
			continue
			
		elif line[0] == 'A':
			well_dict['A1'].append(str(line[1]))
			well_dict['A2'].append(str(line[2]))
			well_dict['A3'].append(str(line[3]))
			well_dict['A4'].append(str(line[4]))
			well_dict['A5'].append(str(line[5]))
			well_dict['A6'].append(str(line[6]))
			well_dict['A7'].append(str(line[7]))
			well_dict['A8'].append(str(line[8]))
			well_dict['A9'].append(str(line[9]))
			well_dict['A10'].append(str(line[10]))
			well_dict['A11'].append(str(line[11]))
			well_dict['A12'].append(str(line[12]))

		elif line[0] == 'B':
			well_dict['B1'].append(str(line[1]))
			well_dict['B2'].append(str(line[2]))
			well_dict['B3'].append(str(line[3]))
			well_dict['B4'].append(str(line[4]))
			well_dict['B5'].append(str(line[5]))
			well_dict['B6'].append(str(line[6]))
			well_dict['B7'].append(str(line[7]))
			well_dict['B8'].append(str(line[8]))
			well_dict['B9'].append(str(line[9]))
			well_dict['B10'].append(str(line[10]))
			well_dict['B11'].append(str(line[11]))
			well_dict['B12'].append(str(line[12]))

		elif line[0] == 'C':
			well_dict['C1'].append(str(line[1]))
			well_dict['C2'].append(str(line[2]))
			well_dict['C3'].append(str(line[3]))
			well_dict['C4'].append(str(line[4]))
			well_dict['C5'].append(str(line[5]))
			well_dict['C6'].append(str(line[6]))
			well_dict['C7'].append(str(line[7]))
			well_dict['C8'].append(str(line[8]))
			well_dict['C9'].append(str(line[9]))
			well_dict['C10'].append(str(line[10]))
			well_dict['C11'].append(str(line[11]))
			well_dict['C12'].append(str(line[12]))

		elif line[0] == 'D':
			well_dict['D1'].append(str(line[1]))
			well_dict['D2'].append(str(line[2]))
			well_dict['D3'].append(str(line[3]))
			well_dict['D4'].append(str(line[4]))
			well_dict['D5'].append(str(line[5]))
			well_dict['D6'].append(str(line[6]))
			well_dict['D7'].append(str(line[7]))
			well_dict['D8'].append(str(line[8]))
			well_dict['D9'].append(str(line[9]))
			well_dict['D10'].append(str(line[10]))
			well_dict['D11'].append(str(line[11]))
			well_dict['D12'].append(str(line[12]))

		elif line[0] == 'E':
			well_dict['E1'].append(str(line[1]))
			well_dict['E2'].append(str(line[2]))
			well_dict['E3'].append(str(line[3]))
			well_dict['E4'].append(str(line[4]))
			well_dict['E5'].append(str(line[5]))
			well_dict['E6'].append(str(line[6]))
			well_dict['E7'].append(str(line[7]))
			well_dict['E8'].append(str(line[8]))
			well_dict['E9'].append(str(line[9]))
			well_dict['E10'].append(str(line[10]))
			well_dict['E11'].append(str(line[11]))
			well_dict['E12'].append(str(line[12]))

		elif line[0] == 'F':
			well_dict['F1'].append(str(line[1]))
			well_dict['F2'].append(str(line[2]))
			well_dict['F3'].append(str(line[3]))
			well_dict['F4'].append(str(line[4]))
			well_dict['F5'].append(str(line[5]))
			well_dict['F6'].append(str(line[6]))
			well_dict['F7'].append(str(line[7]))
			well_dict['F8'].append(str(line[8]))
			well_dict['F9'].append(str(line[9]))
			well_dict['F10'].append(str(line[10]))
			well_dict['F11'].append(str(line[11]))
			well_dict['F12'].append(str(line[12]))

		elif line[0] == 'G':
			well_dict['G1'].append(str(line[1]))
			well_dict['G2'].append(str(line[2]))
			well_dict['G3'].append(str(line[3]))
			well_dict['G4'].append(str(line[4]))
			well_dict['G5'].append(str(line[5]))
			well_dict['G6'].append(str(line[6]))
			well_dict['G7'].append(str(line[7]))
			well_dict['G8'].append(str(line[8]))
			well_dict['G9'].append(str(line[9]))
			well_dict['G10'].append(str(line[10]))
			well_dict['G11'].append(str(line[11]))
			well_dict['G12'].append(str(line[12]))

		elif line[0] == 'H':
			well_dict['H1'].append(str(line[1]))
			well_dict['H2'].append(str(line[2]))
			well_dict['H3'].append(str(line[3]))
			well_dict['H4'].append(str(line[4]))
			well_dict['H5'].append(str(line[5]))
			well_dict['H6'].append(str(line[6]))
			well_dict['H7'].append(str(line[7]))
			well_dict['H8'].append(str(line[8]))
			well_dict['H9'].append(str(line[9]))
			well_dict['H10'].append(str(line[10]))
			well_dict['H11'].append(str(line[11]))
			well_dict['H12'].append(str(line[12]))

		elif line[2] == 'at:': 
			time = str(line[3]).split(':')
			hours=int(time[0])*3600
			#pulling out the hours and converting it to seconds 
			minutes = int(time[1])*60
			#pulling out minutes and converting it to seconds 
			seconds = int(time[2])
			#converting seconds from a string to an integer 
			time = (hours + seconds + minutes)/3600
			#converting seconds to fraction of an hour 
			 
			well_dict['time_hrs'].append(str(time))

	return well_dict
	
	


# Function to interpret well dictionary and write it to a file
def print_wells(well_dict, output_file):

	wells = ['time_hrs','A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','D1','D2','D3','D4','D5','D6','D7','D8','D9','D10','D11','D12','E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12','F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12','G1','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12','H1','H2','H3','H4','H5','H6','H7','H8','H9','H10','H11','H12']
	output_file.write('\t'.join(wells + ['\n']))
	
	for index1 in range(0, len(well_dict['A1'])):
		temp_list = []
		for index2 in wells:
			temp_list.append(well_dict[index2][index1])
		
		output_file.write('\t'.join(temp_list + ['\n']))



# Initialize well dictionary
empty_well_dict = {}
empty_well_dict['A1'] = []
empty_well_dict['A2'] = []
empty_well_dict['A3'] = []
empty_well_dict['A4'] = []
empty_well_dict['A5'] = []
empty_well_dict['A6'] = [] 
empty_well_dict['A7'] = [] 
empty_well_dict['A8'] = [] 
empty_well_dict['A9'] = [] 
empty_well_dict['A10'] = []
empty_well_dict['A11'] = []
empty_well_dict['A12'] = []
empty_well_dict['B1'] = [] 
empty_well_dict['B2'] = [] 
empty_well_dict['B3'] = [] 
empty_well_dict['B4'] = []
empty_well_dict['B5'] = [] 
empty_well_dict['B6'] = [] 
empty_well_dict['B7'] = [] 
empty_well_dict['B8'] = [] 
empty_well_dict['B9'] = [] 
empty_well_dict['B10'] = []
empty_well_dict['B11'] = []
empty_well_dict['B12'] = []
empty_well_dict['C1'] = [] 
empty_well_dict['C2'] = [] 
empty_well_dict['C3'] = [] 
empty_well_dict['C4'] = []
empty_well_dict['C5'] = [] 
empty_well_dict['C6'] = [] 
empty_well_dict['C7'] = [] 
empty_well_dict['C8'] = [] 
empty_well_dict['C9'] = [] 
empty_well_dict['C10'] = []
empty_well_dict['C11'] = []
empty_well_dict['C12'] = []
empty_well_dict['D1'] = [] 
empty_well_dict['D2'] = [] 
empty_well_dict['D3'] = [] 
empty_well_dict['D4'] = []
empty_well_dict['D5'] = []
empty_well_dict['D6'] = []
empty_well_dict['D7'] = []
empty_well_dict['D8'] = []
empty_well_dict['D9'] = []
empty_well_dict['D10'] = []
empty_well_dict['D11'] = []
empty_well_dict['D12'] = []
empty_well_dict['E1'] = []
empty_well_dict['E2'] = []
empty_well_dict['E3'] = []
empty_well_dict['E4'] = []
empty_well_dict['E5'] = []
empty_well_dict['E6'] = []
empty_well_dict['E7'] = []
empty_well_dict['E8'] = []
empty_well_dict['E9'] = []
empty_well_dict['E10'] = []
empty_well_dict['E11'] = []
empty_well_dict['E12'] = []
empty_well_dict['F1'] = []
empty_well_dict['F2'] = []
empty_well_dict['F3'] = []
empty_well_dict['F4'] = []
empty_well_dict['F5'] = []
empty_well_dict['F6'] = []
empty_well_dict['F7'] = []
empty_well_dict['F8'] = []
empty_well_dict['F9'] = []
empty_well_dict['F10'] = []
empty_well_dict['F11'] = []
empty_well_dict['F12'] = []
empty_well_dict['G1'] = []
empty_well_dict['G2'] = []
empty_well_dict['G3'] = []
empty_well_dict['G4'] = []
empty_well_dict['G5'] = []
empty_well_dict['G6'] = []
empty_well_dict['G7'] = []
empty_well_dict['G8'] = []
empty_well_dict['G9'] = []
empty_well_dict['G10'] = []
empty_well_dict['G11'] = []
empty_well_dict['G12'] = []
empty_well_dict['H1'] = []
empty_well_dict['H2'] = []
empty_well_dict['H3'] = []
empty_well_dict['H4'] = []
empty_well_dict['H5'] = []
empty_well_dict['H6'] = []
empty_well_dict['H7'] = []
empty_well_dict['H8'] = []
empty_well_dict['H9'] = []
empty_well_dict['H10'] = []
empty_well_dict['H11'] = []
empty_well_dict['H12'] = []
empty_well_dict['time_hrs']= []

# Construct well time course dictionary
filled_well_dict = read_data(open(sys.argv[1], 'r'), empty_well_dict)

# Create output file
print_wells(filled_well_dict, open(sys.argv[2], 'w'))


# Further analysis should be done in R
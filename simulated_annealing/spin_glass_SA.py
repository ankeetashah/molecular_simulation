#spin glass
#version that goes through each point in lattice individual and picks a spin one by one, taking into account what the nearest neighbors spins and energies were 

import sys
import re
import numpy 
import random
import math

N = 0 				      #NxN lattice (filled in after reading from file)
spins = []
energy_list = []  		  #array of spin orientations and a list of most recent energy added to total energy per point in lattice 
energy = 0.0			  #total energy of spin glass 
#energy0 = 0.0 #initial energy
h = 0					  #Hamiltonian's magnetic field is given as 0 




spin_list = [] 




#caclualtes energy (that can be appended to total energy later)
def energy_calc(Jij, si, sj): 
	return Jij*si*sj
	

def init(): #may rename to init 
	global energy
	energy0 = 0.0
	global N
	line_number = 1							#line ocunter
	point_horizontal = 1 					#point counter (horizontal interactions)
	point_vertical = 1						#point counter (vertical intereactions)
	a = []									#interactions stored here as we parse through input file line by line 
	N_count_upper = 1						#dealing with boundaries of lattice as we parse through vertical interactions 
	calc = 0 


	#open input file 
	with open(sys.argv[1], 'r') as input:
		for line in input: 					#parse through N*2 lines of input 
			a = line.strip().split(' ')

			if line_number == 1:			#amend N as per first line of input file 			
				N = int(line)
				energy_list_matrix = [[0 for x in range (N*N+1)] for y in range(N*N+1)]  #intialize to 0 
				for i in range((N*N)+1):
					spin_list.append(0)

			#lines with horizontal interactions (N-1 entries)
			elif line_number % 2 == 0: 
				for i in range (0, N-1, 1):
					#append energy to matrix
					energy_list_matrix [point_horizontal] [point_horizontal+1] = float(a[i])
					energy_list_matrix [point_horizontal+1] [point_horizontal] = float(a[i])
					#assign a random spin 

					if (spin_list[point_horizontal] == 0):
						rand = numpy.random.random()*2 - 1.0
						if rand < 0: 
							spin_list[point_horizontal] = -1
						else:
							spin_list[point_horizontal] = 1

					if (spin_list[point_horizontal+1] == 0):
						rand = numpy.random.random()*2 - 1.0
						if rand < 0: 
							spin_list[point_horizontal+1] = -1
						else:
							spin_list[point_horizontal+1] = 1

					energy0 += energy_calc(float(a[i]), spin_list[point_horizontal], spin_list[point_horizontal+1] )
					calc += 1
					point_horizontal += 1

			#lines with vertical interactions (N entries)
			else:
				for i in range (0, N, 1): 
				#print "hi"
					if (point_vertical + N > N*N): #moved beyond point boundary of lattice, move to next column  
				 		N_count_upper += 1 
				 		point_vertical = N_count_upper
 					energy_list_matrix [point_vertical] [point_vertical+N] = float(a[i])
 					energy_list_matrix [point_vertical+N] [point_vertical] = float(a[i])

 					if (spin_list[point_vertical] == 0):
						rand = numpy.random.random()*2 - 1.0
						if rand < 0: 
							spin_list[point_vertical] = -1
						else:
							spin_list[point_vertical] = 1

					if (spin_list[point_vertical+N] == 0):
						rand = numpy.random.random()*2 - 1.0
						if rand < 0: 
							spin_list[point_vertical+N] = -1
						else:
							spin_list[point_vertical+N] = 1

					energy0 += energy_calc(float(a[i]), spin_list[point_vertical], spin_list[point_vertical+N])
					calc += 1
 					point_vertical +=N

			line_number += 1
	energy = -h-energy0
	#print "calc init"
	#print calc
	input.close() #close file 
	return energy_list_matrix, line_number, energy


def spin_glass(energy_list_matrix, num_of_lines):
	#perform MC 
	#pick a random spin to flip and see if it changes anything 
	
	line = 2
	point_horizontal = 1
	point_vertical = 1 
	N_count_upper = 1 
	global energy
	energyn = 0.0
	spin_num = random.randint(1, N*N)
	#print spin_num
	
	spin_list[spin_num] *= -1 #effectively flips spin
	#calculate energy of system now

	for i in range (2, num_of_lines):
		if line % 2 == 0: #horizontal interactions
			for i in range (0, N-1, 1):
				energyn += energy_calc(energy_list_matrix[point_horizontal][point_horizontal+1], spin_list[point_horizontal], spin_list[point_horizontal+1] )
				point_horizontal += 1

			
		else: #vertical interactions
			for i in range (0, N, 1): 
				if (point_vertical + N > N*N): #moved beyond point boundary of lattice, move to next column  
				 	N_count_upper += 1 
				 	point_vertical = N_count_upper
				energyn += energy_calc(energy_list_matrix[point_vertical][point_vertical+N], spin_list[point_vertical], spin_list[point_vertical+N])
 				point_vertical +=N

		line += 1

	
	diff = -h-energyn - energy 	
	if diff > 0:
		W = math.exp(-diff)
		#rand = random.randint(0,1)
		#if W < rand:
			#revert back to old spin and retain the value in energy already
		spin_list[spin_num] *= -1
		#else: #accept new spin
			#energy = -h-energyn
		
	else: #if this new energy is smaller than old:
		#keep this new spin and assign this new energyn to energy
		
		energy = -h-energyn
	return energy






def main():
	global energy
	global spin_list
	spin_list_old = []

	# initialize with a random configuration 
	for i in range((N*N)+1):
		spin_list_old.append(0)
	list = init()
	spin_list_old = spin_list
	prev_energy = 0.0

	
	#SA 
	T = 1.0
	while (T != 0):
		for i in range (1, 600, 1):
			if i == 1:
				prev_energy = list[2]
			else:
				temp = init()
				if temp[2] >= prev_energy: #old energy/config was fine - so must use acceptance probability
				#acceptance probability is: 
					alpha = math.exp((prev_energy - temp[2])/T) #acceptance probability
					#generate random number and compare to alpha
					R = random.uniform(0, 1)
					if alpha > R: #switch to new configuration
						prev_energy = temp[2]
						spin_list_old = spin_list
					else: #if alpha < R, keep old configuration 
						energy = prev_energy
						spin_list = spin_list_old
			
				else: #temp[2] < list[2]: #great - energy is lower in new one (switch)
					prev_energy = temp[2]
					spin_list_old = spin_list

		final_energy = spin_glass(list[0], list[1])
		T *= 0.5
	print final_energy



main()















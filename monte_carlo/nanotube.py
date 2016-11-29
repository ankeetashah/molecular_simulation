import random
import math
import numpy

#J = 2.0
J = 1.0
#eps = 1.0
eps = 0.5
kb = 1.0
total_energy = 0.0

N = 10
liquid = 1
vapor = -1
s = []
e = []

moves = 10000


def calculate_energy():
        energy = 0.0
        wall_interaction = 0.0
        neighbor_interaction = 0.0

        for i in range (0, N-1):
                wall_interaction += s[i]

        for i in range (0, N-2):
                neighbor_interaction += s[i]*s[i+1]
        energy = (4*J*wall_interaction)-(eps*neighbor_interaction)

        return energy


def try_move():
        global total_energy
        energy0 = calculate_energy()
        site = int(numpy.random.random()*(N-1))
        if (s[site] == liquid):
                s[site] = vapor
        else:
                s[site] = liquid
        energyn = calculate_energy()
        if (numpy.random.random() > math.exp((-1/kb)*(energyn-energy0))):
                total_energy += energy0 
                #print energy0
                #switch sites back at that point
                if (s[site] == liquid):
                        s[site] = vapor
                else:
                        s[site] = liquid
        else:
                total_energy += energyn
                #print energyn


def initialize():
        for i in range (0, N-1):
                rand = random.randint(0,1)
                if (rand == 0):
                        s.append(liquid)
                else:
                        s.append(vapor)


def main():
        initialize()
        
        for i in range (1, moves):
                try_move()
                #sample averages every so often
                if i % 10 == 0: #this was for testing
                        average_energy = (total_energy/i)
                        #print i, average_energy
        print("Average Energy: %f" % (total_energy/moves))
                        

main()

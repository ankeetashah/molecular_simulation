#Lennard-Jones fluid simulation

import random 
import numpy as np 
import math

box_total = 200
N = 40
eps = 1.0
sigma = 1.0
tf = 0.1                            #10 au, where 1 au is 2pi 
dt = 0.01
m = 1.0
kb = 0.8                            #determined by 100K/125K 
r = eps
V = (box_total/N)*N*(4.0/3)*(3.14)*((sigma)**3)
l = (V)**(1.0/3)
T = 1.0                             #temperature = 1 

PE = 0.0
KE = 0.0

xx=[]
xy=[]
xz=[]

vx=[]
vy=[]
vz=[]

rx=[]
ry=[]
rz=[]

fx=[]
fy=[]
fz=[]

for i in range (0, N):
        rx.append(0)
        ry.append(0)
        rz.append(0)
        fx.append(0)
        fy.append(0)
        fz.append(0)




def forces():
        global PE
                
        for i in range (0, N):
                fx[i] = 0;
                fy[i] = 0;
                fz[i] = 0;
        for i in range (0, N):
                for j in (0, i):
                        if j != i:
                                r2 = periodic_closest_r(i, j)
                                r2i = 1/r2
                                ri = r2i**0.5
                                r6i = r2i**3
                                ff = -48*ri*r6i*(r6i-0.5)
                                fx[i] += ff*rx[i]*ri 
                                fy[i] += ff*ry[i]*ri
                                fz[i] += ff*rz[i]*ri

                                fx[j] -= ff*rx[i]*ri 
                                fy[j] -= ff*ry[i]*ri
                                fz[j] -= ff*rz[i]*ri
                                
                                PE += -4.0*eps*(r6i-(r6i**2))


def square():
        total = 0.0
        for i in range (0, N):
                total += xx[i]**2 + xy[i]**2 + xz[i]**2
        return total

def integrate(dt):
        global KE, PE
        KE = 0.0
        PE = 0.0
        for i in range (0, N):
                vx[i] += fx[i]/m*dt/2.0
                vy[i] += fy[i]/m*dt/2.0
                vz[i] += fz[i]/m*dt/2.0

                xx[i] += vx[i]*dt
                xy[i] += vy[i]*dt
                xz[i] += vz[i]*dt

        forces()
        for i in range (0, N):
                vx[i] += fx[i]/m*dt/2.0
                vy[i] += fy[i]/m*dt/2.0
                vz[i] += fz[i]/m*dt/2.0
                KE += 0.5*m*(vx[i]**2 + vy[i]**2 + vz[i]**2)

def periodic_closest_r(i, j):
        rx[i] = (xx[j] - xx[i])
        ry[i] = (xy[j] - xy[i])
        rz[i] = (xz[j] - xz[i])
        
        
        while (rx[i] > (l/2)):
                rx[i] -= l
        while (ry[i] > (l/2)):
                ry[i] -= l
        while (rz[i] > (l/2)):
                rz[i] -= l

        while (rx[i] < -(l/2)):
                rx[i] += l
        while (ry[i] < -(l/2)):
                ry[i] += l
        while (rz[i] < -(l/2)):
                rz[i] += l

        return (rx[i]**2 + ry[i]**2 + rz[i]**2)

def generate_pos():
        return np.random.random(3)*l #list of length three with numbers between 0 and 1

def initialize_forces_positions_velocities():
        global KE
        KE = 0.0
        for i in range (0, N):
                new_pos = generate_pos()
                xx.append(new_pos[0])
                xy.append(new_pos[1])
                xz.append(new_pos[2])
                j = 0
                while j < i: #check
                        if periodic_closest_r(i, j) < sigma:
                                new_pos = generate_pos()
                                xx[i] = new_pos[0]
                                xy[i] = new_pos[1]
                                xz[i] = new_pos[2]
                                j = -1
                        j += 1
                        

                
                

                vx.append(np.random.uniform(-0.5, 0.5))      
                vy.append(np.random.uniform(-0.5, 0.5)) 
                vz.append(np.random.uniform(-0.5, 0.5)) 

                KE += 0.5*m*(vx[i]**2 + vy[i]**2 + vz[i]**2)

                #scale 
        target_KE = 1.5*N*kb*T
        scale_factor = (target_KE/KE)**0.5
        
        KE = 0
        #loop over particles and multiply particles by that scale factor
        for i in range(0, N):
                vx[i] *= scale_factor
                vy[i] *= scale_factor
                vz[i] *= scale_factor
                KE += 0.5*m*(vx[i]**2 + vy[i]**2 + vz[i]**2)

        
def main():
        msd = []
        energy_output = open("energy_data.txt","w")
        msd_output = open("msd_data.txt","w")
        max_error = 0.0
        sd = 0.0
        t = 0.0
        curr_E = start_E = KE+PE
        
        initialize_forces_positions_velocities()
        forces()
        start_E = KE+PE
        #energy_output.write("%f %f %f %f\n" % (t, KE, PE, (KE+PE)))
        print("%f %f %f %f\n" % (t, KE, PE, (KE+PE)))
        while t < tf:
                integrate(dt)
                sd += square()
                msd_output.write("%f %f\n" % (t, sd/N))
                t += dt 
                print("%f %f %f %f\n" % (t, KE, PE, (KE+PE)))
                #sample averages
                curr_E = KE + PE
                if abs(curr_E - start_E)/start_E*100.0 > max_error:
                        max_error = abs(curr_E - start_E)/start_E*100.0
                print ("%f %f %f %f\n" % (t, KE, PE, (KE+PE)))
                #energy_output.write("%f %f %f %f\n" % (t, KE, PE, (KE+PE)))

                
        msd = sd/N
        #print("Diffusion Constant = %f" % (msd/2*3*dt))
        print("Energy error = %f%%" % (max_error))


if __name__ == "__main__":
            main()

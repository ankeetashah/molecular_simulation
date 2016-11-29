#MD of a two particle system using Leap-Frog integration 

import math

#initialization
m = 1.0                             #mass
x0 = 1.0                            #bond length
k = 1.0                             #spring constant
x1 = 0.0                            #position1
x2 = 1.25                           #position 2
v1 = v2 = 0.0                       #velocities 
f1 = f2 = 0.0                       #forces
a1 = a2 = 0.0                       #acceleration

T = 2*3.14*math.sqrt((m*m)/(m+m))   #period of oscillation


def force(x1, x2):
        f1 = k*(x2-x1-x0)/m
        return f1                   #force on left particle from right particle

def integrate(t, dt):
        global x1, x2, a1, a2, v1, v2
        f1 = force(x1, x2)          #calculate forces
        f2 = -f1

        a1 = f1/m
        a2 = f2/m

                                    #leap-frog integrator
        v1 = v1 + (0.5*dt*a1)       #advance velocity by half step
        v2 = v2 + (0.5*dt*a2)

        x1 = x1 + v1*dt             #advance position by full step
        x2 = x2 + v2*dt

        v1 = v1 + (0.5*dt*a1)       #complete velocity step (so now it is a full step)
        v2 = v2 + (0.5*dt*a2)


def energy():
        return (0.5*m*((v1**2) + (v2**2)) + (0.5*k*(x2-x1-x0)**2))

#need to get maximum energy over teh entire trajectory
#plot energy error as a function of time
def main():
        #getting timestep
        dt = float(input("Timestep in atomic units: "))
        initialE = (0.5*m*((v1**2) + (v2**2)) + (0.5*k*(x2-x1-x0)**2))
        t = 0.0
        E = 0.0
        Eerror = 0.0
        while (t < (3*T)):
                integrate(t, dt)
                E = energy()
                if ((abs(E - initialE)/initialE*100.0) > Eerror):
                        Eerror = (abs(E - initialE)/initialE*100.0)
                t = t + dt
        print("Energy error: %f%%\n" % Eerror)


if __name__ == "__main__":
            main()

# ----------------------------
# Quantum well with an internal barrier
# ----------------------------
# Finite differences method as developed by Truhlar JCP 10 (1972) 123-132
#
# code by Jordi Faraudo
#  
# edited by Daniel Sarcanean
#
import numpy as np
import matplotlib.pyplot as plt

#Potential as a function of position
def getV(x):
    if np.any(x <= rb) & np.any(x >= lb):
        potvalue = V
    else:
        potvalue = 0
    return potvalue

#Discretized Schrodinger equation in n points (FROM 0 to n-1)
def Eq(n,h,x):
    F = np.zeros([n,n])
    for i in range(0,n):
        F[i,i] = -2*((h**2)*getV(x[i]) + 1)
        if i > 1:
            F[i,i-1] = 1
            if i < n-2:
                F[i,i+1] = 1
    return F

#-------------------------
# Main program
#-------------------------

# Choosing the barrier's conditions:
V = float(input("Select the barrier height (in hartrees, recommended 100 h): ")) 
len_barr = float(input("Select the barrier width (from 0 to 1): " ))
rb = 0.5 + len_barr/2
lb = 0.5 - len_barr/2

# Interval for calculating the wave function [-L/2,L/2]
xlower = 0
xupper = 1

#Discretization options
npoints = 1000 #discretization in space
h = 1 / npoints

#Create coordinates at which the solution will be calculated
x = np.linspace(xlower, xupper, npoints)
print(getV(x))

#Create an equal-size graph for the potential:
x0 = [xlower, lb-0.00001, lb, rb, rb+0.000001, xupper]
y0 = [0, 0, V, V, 0, 0]

#grid size (how many discrete points to use in the range [-xlower,xupper])   
print("Using",npoints, "grid points.")

#Calculation of discrete form of Schrodinger Equation
print("Calculating matrix...")
F=Eq(npoints,h,x)

#diagonalize the matrix F
print("Diagonalizing...")
eigenValues, eigenVectors = np.linalg.eig(F)

#Order results by eigenvalue
# w ordered eigenvalues and vs ordered eigenvectors
idx = eigenValues.argsort()[::-1]   
w = eigenValues[idx]
vs = eigenVectors[:,idx]

#Energy Level
E = - w/(2.0*h**2)
print(E)

#Print Energy Values
print("RESULTS:")
for k in range(0,9):
	print("State ",k," Energy = %.2f" %E[k])

#Init Wavefunction (empty list with npoints elements)
psi = [None]*npoints

#Calculation of normalised Wave Functions
for k in range(0,len(w)):
	psi[k] = vs[:,k]
	integral = h*np.dot(psi[k],psi[k])
	psi[k] = psi[k]/integral**0.5
 
#Y values for potential

#Plot Wave functions
print("Plotting")

#v = int(input("\n Quantum Number (enter 0 for ground state):\n>"))
for v in range(0,2):
    plt.fill_between(x0,0,y0, color="green")
    plt.plot(x,10*psi[v]+E[v],label=r'$\psi_v(x)$, k = ' + str(v))
    plt.plot(x, 15*(psi[v])**2 + E[v], color='red', linestyle = '--', label=r'$|\psi_v (x)|^2$, k =' + str(v))
    plt.axhline(E[v], color = 'black', linestyle = ':')
    plt.title(r'$n=$'+ str(v) + r', $E_n$=' + '{:.2f}'.format(E[v]))
    plt.xlim(xlower,xupper)
    plt.legend()
    plt.xlabel(r'$x$(dimensionless)')
    plt.ylabel(r'$E(hartrees)$')
    plt.show()
print("Bye")

#
# Example parabolic motion
#

# Here we import the mathematical library and the plots library
import numpy as np
import matplotlib.pyplot as plt
 
#
# FUNCTION DRAW A TRAJECTORY FOR PARABOLIC MOTION
# Input: velocity and angle 
#
def draw_trajectory(u, theta, reb):
    
    #convert angle in degrees to rad
    theta = np.radians(theta)
    
    #gravity acceleration in m/s2
    g = 9.8
    
    # Set the velocity components at initial conditions
    vx = u*np.cos(theta)
    vy = u*np.sin(theta)
    
    # Time of flight initially
    tf = 2 * vy / g
    
    # Find time intervals
    intervals = np.arange(0, tf, 0.001)
    
    # Create an empty list of x and y coordinates
    x = []
    y = []
    
    # Do a loop over time calculating the coordinates
    for t in intervals:
        x.append(vx * t)
        y.append(vy * t - 0.5 * g * t ** 2)

    # I define tprev just in case:
    tprev = 0

    # Loop for rebounds:
    for i in range(reb+1):
        if i == 0:
            continue
        
        vy = coef*(vy)
        tprev = tprev + tf
        print(vy * tf)
        tf = 2 * vy / g
        intervals = np.arange(tprev, tprev+tf, 0.001)
        for t in intervals:
            x.append( vx * t)
            y.append(vy * (t-tprev) - 0.5 * g * (t-tprev) ** 2)
           
    #Plot the results
    plt.plot(x, y, label=vel)
    plt.xlabel('Distance (m)')
    plt.ylabel('Height (m)')
    plt.title('Projectile motion')

#--------------------------------------------------------------------------------
# Main Program: give specific values and call to the function draw_trajectory
#--------------------------------------------------------------------------------

print("Parabolic motion of a projectile\n")

# Ask the user for angle
print("Enter desired launch angle in degrees (recommended 45 degrees):")
theta=float(input())

# Insert a velocity in m/s
vel = float(input("Insert your desired velocity: "))

# How many rebounds the user requires
reb = int(input("Insert a rebound number (Less than 10 is recommended): "))

# Rebound coefficient selection
coef = float(input("Insert a number from 0 to 1 for your rebound coefficient: "))
if coef > 1:
    coef = float(input("Insert a number from 0 to 1 for your rebound coefficient: "))

draw_trajectory(vel, theta, reb)

# Add a legend and show the graph
plt.legend()
plt.ylim(bottom=0)
plt.show()

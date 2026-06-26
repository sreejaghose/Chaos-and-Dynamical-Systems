import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def chen_ode(y, t, a, b, c):
    dydt = [a*(y[1]-y[0]),
            (c-a)*y[0] - y[0]*y[2] + c*y[1],
            y[0]*y[1] - b*y[2]]
    return dydt

def chen():
    # Initial conditions
    x0 = -0.1
    y0 = 0.5
    z0 = -0.6
    
    # Parameters
    a = 40
    b = 3
    c = 28
    
    # Time points
    t = np.linspace(0, 100, 10000)
    
    # Numerically solve
    sol = odeint(chen_ode, [x0, y0, z0], t, args=(a, b, c))
    
    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(sol[:,0], sol[:,1], sol[:,2])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('Chen double scroll attractor')
    plt.show()

# Call the function
chen()


#! /usr/bin/env python3

'''
Analytical solution to the Time-Dependent Schrodinger equation
for a particle in an infinite square well

author: Timothy Holmes
email: tpholmes7@gmail.com
website: http://timothypholmes.github.io


Still working to finish a few minor user friendly features. Please update how
would like.
'''

######################################################################
#Import libraries
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import style
import seaborn
import numpy as np
import pandas as pd
import warnings
import pandas as pd
#import scipy.special
import numpy.polynomial
import csv
from mpmath import mp



class quantum_harmonic_oscillator:
    '''
    Class that take an approach at solving the Time-Dependent Schrodinger
    equation for an unbounded particle in an infinite square well and
    generates a gaussian wave function for the psi initial.
    '''
    def __init__(self, x, state, E0, N, hbar, omega, m):
        self.x = x
        self.state = state
        self.E0_energy = E0
        self.N_steps = N
        self.hbar = hbar
        self.omega = omega
        self.mass = m


        '''
        Parameters
        ----------
        x : array, float
            length-N array of evenly spaced spatial coordinates
        psi_x0 : array, complex
            Initial wave function at time t0 = 0
        psi_xt : array, complex
            Time-dependent Schrodinger equation
        hbar : scalar, float
            value of planck's constant
        m : scalar, float
            particle mass
        quantum_number : scalar, integer
            values of conserved quantities in a quantum system
        total_time : float
        Time_step : scalar, integer
        '''

    def xi(self):
        self.xi = (np.sqrt(m * self.omega / self.hbar) * self.x)#.reshape(len(x),1)

    '''
    def hermite(self, xi, n):
        if (n == 0).any():
            return 1
        elif n == 1:
            return 2*xi
        else:
            return 2*xi*hermite(xi, n-1)-2*(n-1)*hermite(xi,n-2)

        for n in range(1, N + 1):
            print(hermite(xi, n))
    '''

    def hermite(n):

        herm = np.zeros((len(x), N))
        print(str(herm.shape))
        herm[:,0] = np.ones((len(x)))
        herm[:,1] = 2*self.xi

        if (n == 0):
            f = herm[:,0]
        if (n == 1):
            f = herm[:,1]
        elif (n > 1):

            a = 1
            b = 2

            for n in range(1, n):

                newHerm = 2 * self.xi * herm[:,b] - 2 * (n - 1) * herm[:,a]
                tem = a
                a = b
                b = tem
                herm[:,b] = newHerm;

            f = newHerm


    def state_normalized(self):
        A = np.trapz((np.conj(state)*state)[:,0], x[:,0])
        stateNorm = state/(np.sqrt(A))
        #print('scalar A: ' + str(A))
        #print(str(state_normalized.shape))


    def eigen_values(self):

        A = ((self.mass * self.omega / (np.pi * self.hbar)) ** (.25))

        phi = np.zeros((len(x),1),dtype=float).reshape(len(x),1)
        En = np.zeros((N,1),dtype=float).reshape(N,1)
        Cn = np.zeros((N,1),dtype=float).reshape(N,1)

        for n in range(0,N):


            phi[:,n] = A * (1 / np.sqrt(2 ** (n) * np.math.factorial(n))) * self.hermite(self.xi, n) * np.exp((-(self.xi[:,0]) ** 2) / 2)
            #En[:,n] = hbar*omega*(n-1+.5)
            #Cn[n] = np.trapz((np.conj(phi)*stateNorm)[:,n], x[:,0])


    def psi_xt(self):
        dt = 1
        timeTotal = (1*10**3)/2

        for j in range(0, 10000, dt):

            time = j*1**-18
            psi_xt = np.zeros((len(x),1),dtype=complex).reshape(len(x),1)

            for k in range(0, quantum_number):

                psi_xt[:,0] = psi_xt[:,0] + (Cn[0,k]*phi[:,k]*(np.exp((-1j*En[0,k]*time)/hbar))) #(2001,1)

                count += 1

            center = max(np.abs(psi_xt))
            plt.plot(x, np.real(psi_xt),'r', label=r'$\mathbb{R} \psi(x,t)$', linewidth = 0.75)
            plt.plot(x, np.imag(psi_xt),'b', label=r'$\mathbb{C} \psi(x,t)$', linewidth = 0.75)
            plt.plot(x, np.abs(psi_xt),'y', label=r'$|\psi(x,t)|$', linewidth = 0.75)

            x_min = min(x[:,0]-1)
            x_max = max(x[:,0]+1)
            psi_min = -A
            psi_max = A
            plt.xlim((x_min, x_max))
            #plt.ylim((psi_min, psi_max))

            plt.legend(prop=dict(size=6))
            center_line = plt.axvline(center, c='k', ls=':')
            left_wall_line = plt.axvline(0, c='k', linewidth=2)
            right_well_line = plt.axvline(x[-1], c='k', linewidth=2)


            #anim = FuncAnimation(fig, animate, interval=100, frames=len(t)-1)

            plt.pause(0.01)
            plt.draw()
            plt.clf()
            plt.cla()
            #plt.close()



df = pd.read_csv("/Users/timholmes/Desktop/Github/harmonic_oscillator/data.csv", header=None)
print(df[3])
x = df[0].apply(lambda x: np.complex(x), '%.100f')
state1 = df[1].apply(lambda x: np.complex(x), '%.100f')
state2 = df[2].str.replace(' ','').str.replace('i', 'j').apply(lambda x: np.complex(x), '%.100f')
state3 = df[3].apply(lambda x: np.complex(x), '%.100f')
state4 = df[4].str.replace(' ','').str.replace('i', 'j').apply(lambda x: np.complex(x), '%.100f')

x = np.array(x).reshape(len(x),1)
state1 = np.array(state1).reshape(len(x),1)
state2 = np.array(state2).reshape(len(x),1)
state3 = np.array(state3).reshape(len(x),1)
state4 = np.array(state4).reshape(len(x),1)

E0 = 10
mass = 511000
N = 50
hbar = 6.582*10**-16
omega = 2*E0/hbar
m =  int(mass/(3e8**2))
state = state1
n = np.arange(1,N+1).reshape(1,N)
xi = (np.sqrt(m * omega / hbar) * x)

harmonic = quantum_harmonic_oscillator(x, state, E0, N, hbar, omega, m)

harmonic.xi()
harmonic.hermite(n, xi)
harmonic.state_normalized()
harmonic.eigen_values()

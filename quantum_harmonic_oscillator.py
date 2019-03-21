#! /usr/bin/env python3

'''
Analytical solution to the Time-Dependent Schrodinger equation
for a particle in an infinite square well

author: Timothy Holmes
email: tpholmes7@gmail.com
website: http://timothypholmes.github.io

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

class quantum_harmonic_oscillator:
    '''
    Class that take an approach at solving the Time-Dependent Schrodinger
    equation for an unbounded particle in an infinite square well and
    generates a gaussian wave function for the psi initial.
    '''
    def __init__(self,hbar,m,quantum_number,total_time,dt,
        L,x,n,a,l):
        self.hbar = hbar
        self.mass = m
        self.quantum_number = quantum_number
        self.total_time = total_time
        self.time_step = dt
        self.length = L
        self.x = x
        self.n = n
        self.a = a

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

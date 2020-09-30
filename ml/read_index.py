#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 12:31:40 2020

@author: agus
"""

import numpy as np

'''

For example we have btp-index, we want to read it and give some format.
This code generates one list, each element is a list which contains tuples with (parameter, value) and in the last position the simulation file name.

[ 
[(parameter, value), (parameter, value), ..., "simulation_file_name"],
[(parameter, value), (parameter, value), ..., "simulation_file_name"],
....
[(parameter, value), (parameter, value), ..., "simulation_file_name"],
]

'''

index_file = "./simulations/btp-index.txt"


all_results=[]
simulation=[]
with open(index_file) as file:
    for line in file:
        line = line.strip() #preprocess line
        line = line.split(",") #split the tuples
        
        for i in range(len(line)-1):  # the last element is the simulation name
            tmp = line[i].split(" ")
            
            simulation.append((str(tmp[0]), float(tmp[1]))) 
        simulation.append(str(line[-1]))  #the simulation name
        all_results.append(simulation)        


# [ [(param, value), (param,value), ..., "simulation name"], [...]...]
# now we can look for centain parameter and value (or more than one)
# example, we want to look if we have some simulation with the 'P' value greater than 8

parameter='P'
value=8

for simu in all_results:
    for i in range(len(simu)-1):
        simu_name=simu[-1]
        tuple = simu[i]
        if tuple[0] == parameter and tuple[1]>value:
            print(tuple, "in Simulation ==> " + simu_name)
        

# When we have the simulation file name (i.e. ./simulations/btp-simu-4.csv), we can open it using numpy
      

simulation = np.genfromtxt('./simulations/btp-simu-4.csv', delimiter=',')

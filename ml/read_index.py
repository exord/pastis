#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 12:31:40 2020

@author: agus
"""

import numpy as np

## formating pla-index

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

# now we can look for centain parameter and value
# example, we want to look if we have some simulation with the 'P' value greater than 8

parameter='P'
value=8

for simu in all_results:
    for i in range(len(simu)-1):
        simu_name=simu[-1]
        tuple = simu[i]
        if tuple[0] == parameter and tuple[1]>value:
            print(tuple, "in Simulation ==> " + simu_name)
        
'''
 to read the simulation file
      
./simulations/btp-simu-4.csv

'''

simulation = np.genfromtxt('./simulations/btp-simu-4.csv', delimiter=',')

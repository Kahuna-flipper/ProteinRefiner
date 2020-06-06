import numpy as np
import math 
import random
import CostFunction

ngenerations = 20
npopulations = 20
Natoms = 1388
X_temp = np.zeros((Natoms,1))
Y_temp = np.zeros((Natoms,1))
Z_temp = np.zeros((Natoms,1))

X_trial = np.zeros((Natoms,1))
Y_trial = np.zeros((Natoms,1))
Z_trial = np.zeros((Natoms,1))

X_target = np.zeros((Natoms,1))
Y_target = np.zeros((Natoms,1))
Z_target = np.zeros((Natoms,1))




F = 0.5
recomb = 0.6


def Bound_Check(vector,bounds):
    for i in range(0,3):
        for j in range(0,Natoms):
            temp_magnitude = vector[j][i]
            if (temp_magnitude < bounds[i][0]):
                vector[j][i] = bounds[i][0]
            if (temp_magnitude > bounds[i][1]):
                vector[j][i] = bounds[i][1]
            

def Main(X,Y,Z,data):
    xmin = np.min(X)
    ymin = np.min(Y)
    zmin = np.min(Z)
    
    xmax = np.max(X)
    ymax = np.max(Y)
    zmax = np.max(Z)
    
    bounds = [(xmin,xmax),(ymin,ymax),(zmin,zmax)]
    population = []
    for i in range(0,npopulations):
        for j in range(0,Natoms):
            X_temp[j,0] = xmin + random.random()*(xmax-xmin)
            Y_temp[j,0] = ymin + random.random()*(ymax-ymin)
            Z_temp[j,0] = zmin + random.random()*(zmax-zmin)
            
        individual = np.concatenate((X_temp,Y_temp,Z_temp),axis=1)
        population.append(individual)
    
        
    # Mutation :
    for i in range(0,ngenerations):
        print ('GENERATION:'),i
        gen_scores = []

        for j in range(0,npopulations):
            temp_range = list(range(0,npopulations))
            temp_range.remove(j)
            indexes = random.sample(temp_range,3) # Choosing 3 random vectors for mutation
            rand_index1 = indexes[0]
            rand_index2 = indexes[1]
            rand_index3 = indexes[2]
            
            target_vector = population[j]
            donor_vector = population[rand_index1] + F*(population[rand_index2]-population[rand_index3])
            
            Bound_Check(donor_vector,bounds)
            
            # Crossover Operation : 
            trial_vector = []
            
            for i in range(0,len(donor_vector)):
                crossover_var = random.random()
                
                if(crossover_var <= recomb) :
                    trial_vector.append(donor_vector[i])
                
                
                else :
                    trial_vector.append(target_vector[i])
            
            for i in range(0,Natoms):
                X_trial[i] = trial_vector[i][0]
                Y_trial[i] = trial_vector[i][1]
                Z_trial[i] = trial_vector[i][2]
                X_target[i] = target_vector[i][0]
                Y_target[i] = target_vector[i][1]
                Z_target[i] = target_vector[i][2]
                
            temp_trial_vector = np.concatenate((X_trial,Y_trial,Z_trial),axis=1)
                
            #Selection :    
            trial_cost = CostFunction.CostFunction(X_trial,Y_trial,Z_trial,data)
            target_cost = CostFunction.CostFunction(X_target,Y_target,Z_target,data)
            
            if (trial_cost < target_cost) :
                population[j] = temp_trial_vector
                gen_scores.append(trial_cost)
                print ('   >'),trial_cost, trial_vector
            else:
                print ('   >'),target_cost, target_vector
                gen_scores.append(target_cost)
                
    gen_sol = population[gen_scores.index(min(gen_scores))]
    return gen_sol
        

            
            
            
        
                    
            
            
            
            
        
            
            
            
            
    
    
    
    
    

 

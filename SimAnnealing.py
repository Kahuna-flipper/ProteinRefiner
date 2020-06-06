import numpy as np
import matplotlib.pyplot as plt 
import CostFunction
import random

# Generate random start point 
Natoms = 1388
bounds = [(100,250),(-25,25),(-30,30)]
def SimAnneal(X,Y,Z,T,datalist) : 
    scale = np.sqrt(T)
    #generate random start coordinates :
    
    X0 = np.zeros((Natoms,1))
    Y0 = np.zeros((Natoms,1))
    Z0 = np.zeros((Natoms,1))
    
    X_n = np.zeros((Natoms,1))
    Y_n = np.zeros((Natoms,1))
    Z_n = np.zeros((Natoms,1))
    
    temp1 = np.random.normal()
    while(temp1<-1 or temp1<1):
        temp1 = np.random.normal()
    k1 = temp1
    
    temp1 = np.random.normal()
    while(temp1<-1 or temp1<1):
        temp1 = np.random.normal()
    k2 = temp1
    
    temp1 = np.random.normal()
    while(temp1<-1 or temp1<1):
        temp1 = np.random.normal()
    k3 = temp1
    
    X0 = X + 1.5*k1*np.ones((Natoms,1))
    Y0 = Y + 1.5*k2*np.ones((Natoms,1))
    Z0 = Z + 1.5*k3*np.ones((Natoms,1))
    
    coordinates = np.concatenate((X0,Y0,Z0),axis=1)
    temp_coordinates = coordinates 
    
    currenty = CostFunction.CostFunction(X0,Y0,Z0,datalist)
    
    history = [coordinates]
    
    neighbours = np.concatenate((X_n,Y_n,Z_n),axis=1)
    
    for i in range(0,5):
        print("Epoch : ",(i+1))
        k = np.random.normal()
        
        #Xn,Yn,Zn = np.zeros((Natoms,1)),np.zeros((Natoms,1)),np.zeros((Natoms,1))
        
        #neighbours = np.concatenate((Xn,Yn,Zn),axis=1)
        neighbours = coordinates + k*np.ones((Natoms,3))/50
        
        for i in range(0,Natoms):
            X_n[i]=neighbours[i,0]
            Y_n[i]=neighbours[i,1]
            Z_n[i]=neighbours[i,2]
        
        for j in range(0,3):
            while(np.max(neighbours[:,j]) > bounds[j][1] or np.min(neighbours[:,j]) < bounds[j][0]):
                k = np.random.normal()
                neighbours = coordinates + k*np.ones((Natoms,3))
                
        E_neighbour = CostFunction.CostFunction(X_n,Y_n,Z_n,datalist)
        
        delE = E_neighbour - currenty
        
        if(delE>0):
            prob = 1/(1+np.exp(delE/T))
            
            if(random.random()<prob):
                accept = True 
            else :
                accept = False 
                
        else :
            accept = True 
        
        if(accept == True):
            currenty = E_neighbour
            coordinates = neighbours
            
        T = T*random.random()
        
    return coordinates
            
        
            
        
       
        
        
        
        
        
        
        
                

                
                
                
        
        
        
    
    
    

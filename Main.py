import numpy as np

import pandas as pd
import math
import CostFunction
import DifferentialEvolution
import random
from sklearn.cluster import KMeans

data = pd.read_csv('myDataFile.csv')
pdb = pd.read_csv('pdbfile.csv')


# Loading necessary sizes for the parameters : 
Natoms = 1388
Ntypes = 14
Nbonh = 691
Mbona = 710
Ntheth = 1582
Mtheta = 951
Nphih = 3158
Mphia = 2917
Nhparm = 0 
Nparm = 0
Nnb = 7656
Nres = 88
Nbona = 710
Ntheta = 951
Nphia = 2917
Numbnd = 65
Numang = 152
Nptra = 191
Natyp = 33
Nphb = 0
Ifpert =0 
Nbper = 0
Ngper = 0
Ndper =0
Mbper = 0
Mgper = 0
Mdper = 0
Ifbox = 0
Nmxrs = 0
Ifcap = 24
Numextra =0 
Ncopy =   0



indexvar=pd.Index(data['DATE'])

p=indexvar.get_loc('ATOM_TYPE_INDEX')
temp_data=data[p+2:p+2+Natoms].astype(float)
Atom_type_index=temp_data['DATE'].values[0:Natoms]



p=indexvar.get_loc('CHARGE')
temp_data=data[p+2:p+2+Natoms].astype(float)
Charge=temp_data['DATE'].values[0:Natoms]


p=indexvar.get_loc('NONBONDED_PARM_INDEX')
temp_data=data[p+2:p+2+Ntypes*Ntypes].astype(float)
NBparmindex=temp_data['DATE'].values[0:Ntypes*Ntypes]



p=indexvar.get_loc('BOND_FORCE_CONSTANT')
temp_data=data[p+2:p+2+Numbnd].astype(float)
bond_force_constant=temp_data['DATE'].values[0:Numbnd]


p=indexvar.get_loc('BOND_EQUIL_VALUE')
temp_data=data[p+2:p+2+Numbnd].astype(float)
bond_equ_distance=temp_data['DATE'].values[0:Numbnd]


p=indexvar.get_loc('ANGLE_FORCE_CONSTANT')
temp_data=data[p+2:p+2+Numang].astype(float)
angle_force_constant=temp_data['DATE'].values[0:Numang]


p=indexvar.get_loc('ANGLE_EQUIL_VALUE')
temp_data=data[p+2:p+2+Numang].astype(float)
equ_angle=temp_data['DATE'].values[0:Numang]  # Changed this from numbnd to numang


p=indexvar.get_loc('DIHEDRAL_FORCE_CONSTANT')
temp_data=data[p+2:p+2+Nptra].astype(float)
dihedral_force_constant=temp_data['DATE'].values[0:Nptra]


p=indexvar.get_loc('DIHEDRAL_PERIODICITY')
temp_data=data[p+2:p+2+Nptra].astype(float)
dihedral_periodicity=temp_data['DATE'].values[0:Nptra]


p=indexvar.get_loc('DIHEDRAL_PHASE')
temp_data=data[p+2:p+2+Nptra].astype(float)
dihedral_phase=temp_data['DATE'].values[0:Nptra]


p=indexvar.get_loc('LENNARD_JONES_ACOEF')
t=int(Ntypes*(Ntypes+1)/2)
temp_data=data[p+2:p+2+t].astype(float)
lennard_a=temp_data['DATE'].values[0:t]


p=indexvar.get_loc('LENNARD_JONES_BCOEF')
temp_data=data[p+2:p+2+t].astype(float)
lennard_b=temp_data['DATE'].values[0:t]


p=indexvar.get_loc('BONDS_INC_HYDROGEN')
temp_data=data[p+2:p+2+(3*Nbonh)].astype(float)
bonds_inc_hyd=temp_data['DATE'].values[0:3*Nbonh]


p=indexvar.get_loc('BONDS_WITHOUT_HYDROGEN')
temp_data=data[p+2:p+2+(3*Nbona)].astype(float)
bonds_without_hyd=temp_data['DATE'].values[0:3*Nbona]


p=indexvar.get_loc('ANGLES_INC_HYDROGEN')
temp_data=data[p+2:p+2+(4*Ntheth)].astype(float)
angles_inc_hyd=temp_data['DATE'].values[0:4*Ntheth]


p=indexvar.get_loc('ANGLES_WITHOUT_HYDROGEN')
temp_data=data[p+2:p+2+(4*Ntheta)].astype(float)
angles_without_hyd=temp_data['DATE'].values[0:(4*Ntheta)]


p=indexvar.get_loc('DIHEDRALS_INC_HYDROGEN')
temp_data=data[p+2:p+2+(5*Nphih)].astype(float)
dihedrals_inc_hyd=temp_data['DATE'].values[0:(5*Nphih)]


p=indexvar.get_loc('DIHEDRALS_WITHOUT_HYDROGEN')
temp_data=data[p+2:p+2+(5*Nphia)].astype(float)
dihedrals_without_hyd=temp_data['DATE'].values[0:(5*Nphia)]

# Extracting PDB Data :
X = np.zeros((Natoms,1))
Y = np.zeros((Natoms,1))
Z = np.zeros((Natoms,1))

i=0
while i<Natoms :
    X[i,0] = pdb.values[i][0]
    Y[i,0] = pdb.values[i][1]
    Z[i,0] = pdb.values[i][2]
    i=i+1
    
    
    
data_list = [Atom_type_index,Charge,NBparmindex,bond_force_constant,bond_equ_distance,angle_force_constant,equ_angle,dihedral_force_constant,dihedral_periodicity,dihedral_phase,lennard_a,lennard_b,bonds_inc_hyd,bonds_without_hyd,angles_inc_hyd,angles_without_hyd,dihedrals_inc_hyd,dihedrals_without_hyd]   
#Energy = CostFunction.CostFunction(X,Y,Z,data_list)
#trial_vector = DifferentialEvolution.Main(X,Y,Z,data_list)\





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

final_trial_vector = []




F = 0.9
recomb = 0.5
K = 0.5
decider_list = [1,-1]
#bounds = [(np.min(X),np.max(X)),(np.min(Y),np.max(Y)),(np.min(Z),np.max(Z))]


def Bound_Check(vector,target):
    for i in range(0,3):
        for j in range(0,Natoms):
            temp_magnitude = vector[j][i]
            bound = [(target[j,0]-10,target[j,0]+10),(target[j,1]-5,target[j,1]+5),(target[j,2]-2.5,target[j,2]+2.5)]
            while (temp_magnitude < bound[i][0] or temp_magnitude > bound[i][1]):
                 temp_magnitude = np.median(target[0:Natoms,i])+random.choice(decider_list)*random.random()*(np.max(target[0:Natoms,i])-np.min(target[0:Natoms,i]))
                  #temp_magnitude = np.min(target[0:Natoms,i])+random.random()*(np.max(target[0:Natoms,i])-np.min(target[0:Natoms,i]))
            vector[j][i] = temp_magnitude
      
            

def Main(X,Y,Z,data,ngenerations,npopulations):
    xmin = np.min(X)
    ymin = np.min(Y)
    zmin = np.min(Z)
    
    xmax = np.max(X)
    ymax = np.max(Y)
    zmax = np.max(Z)
    
    coordinates = np.concatenate((X,Y,Z),axis=1)
    
    
    
    population = []
    for i in range(0,npopulations):
        for j in range(0,Natoms):
#            X_temp[j] = X[j] + random.choice(decider_list)*(random.random())*(X[j])/100
#            Y_temp[j] = Y[j] + random.choice(decider_list)*(random.random())*(Y[j])/100
#            Z_temp[j] = Z[j] + random.choice(decider_list)*random.random()*(Z[j])/100
            
          
            
#            X_temp[j] = xmin + (random.random())*(xmax-xmin)
#            Y_temp[j] = ymin + (random.random())*(ymax-ymin)
#            Z_temp[j] = zmin + (random.random())*(zmax-zmin)
            
            X_temp[j] = np.median(X) + random.choice(decider_list)*(random.random())*(xmax-xmin)
            Y_temp[j] = np.median(Y) + random.choice(decider_list)*(random.random())*(ymax-ymin)
            Z_temp[j] = np.median(Z) + random.choice(decider_list)*(random.random())*(zmax-zmin)
            
        individual = np.concatenate((X_temp,Y_temp,Z_temp),axis=1)
        population.append(individual)
    
        
    # Mutation :
    for i in range(0,ngenerations):
        print ('GENERATION:',i+1)
        gen_scores = []
        trial_vector_list = []

        for j in range(0,npopulations):
            temp_range = list(range(0,npopulations))
            temp_range.remove(j)
            indexes = random.sample(temp_range,3) # Choosing 3 random vectors for mutation
            rand_index1 = indexes[0]
            rand_index2 = indexes[1]
            rand_index3 = indexes[2]
            
            
            target_vector = population[j]
            donor_vector = population[rand_index1] + F*(population[rand_index2]-population[rand_index3])
            
            Bound_Check(donor_vector,coordinates)
            Bound_Check(target_vector,coordinates)
            #Bound_Check(donor_vector,bounds)
            
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
            final_trial_vector.append(temp_trial_vector)
            trial_vector_list.append(temp_trial_vector)
                
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

        

solution13 = Main(X,Y,Z,data_list,15,30)
solution14 = Main(X,Y,Z,data_list,15,40)
solution15 = Main(X,Y,Z,data_list,15,50)



        
        
        
    
    


    
    
    
    
    

    

    
    
    
    









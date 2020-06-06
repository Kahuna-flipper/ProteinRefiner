import numpy as np
import math 


def CostFunction(X,Y,Z,data):
    Atom_type_index = data[0]
    Charge = data[1]
    NBparmindex = data[2]
    bond_force_constant = data[3]
    bond_equ_distance = data[4]
    angle_force_constant = data[5]
    equ_angle = data[6]
    dihedral_force_constant = data[7]
    dihedral_periodicity = data[8]
    dihedral_phase = data[9]
    lennard_a = data[10]
    lennard_b = data[11]
    bonds_inc_hyd = data[12]
    bonds_without_hyd = data[13]
    angles_inc_hyd = data[14]
    angles_without_hyd = data[15]
    dihedrals_inc_hyd = data[16]
    dihedrals_without_hyd = data[17]
    Natoms = 1388
    Ntypes = 14
    
    
    energy_bond_inc = 0
    temp_length = np.size(bonds_inc_hyd)/3
    r_inc_hyd = np.zeros((int(temp_length),1),dtype=float)

    for i in range(0,int(temp_length)):
      index1 = (bonds_inc_hyd[(i-1)*3])
      index2 = (bonds_inc_hyd[1+(i-1)*3])
      index3 = int(bonds_inc_hyd[2+(i-1)*3])
      corrected_index1 = int(index1/3 + 1)
      corrected_index2 = int(index2/3 + 1)
      x1 = X[corrected_index1-1]
      y1 = Y[corrected_index1-1]
      z1 = Z[corrected_index1-1]
    
      x2 = X[corrected_index2-1]
      y2 = Y[corrected_index2-1]
      z2 = Z[corrected_index2-1]
    
      distance = np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
      r_inc_hyd[i]=distance
      bond_constant = bond_force_constant[index3-1]
      equ_distance = bond_equ_distance[index3-1]
      energy_bond_inc = energy_bond_inc + bond_constant*(distance-equ_distance)**2
    

#Bonds excluding Hydrogen : 
    energy_bond_exc = 0
    temp_length = np.size(bonds_without_hyd)/3
    r_exc_hyd = np.zeros((int(temp_length),1),dtype=float)

    for i in range(0,int(temp_length)):
      index1 = (bonds_without_hyd[(i-1)*3])
      index2 = (bonds_without_hyd[1+(i-1)*3])
      index3 = int(bonds_without_hyd[2+(i-1)*3])
      corrected_index1 = int(index1/3 + 1)
      corrected_index2 = int(index2/3 + 1)
      x1 = X[corrected_index1-1]
      y1 = Y[corrected_index1-1]
      z1 = Z[corrected_index1-1]
    
      x2 = X[corrected_index2-1]
      y2 = Y[corrected_index2-1]
      z2 = Z[corrected_index2-1]
    
      distance = np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
      r_exc_hyd[i]=distance
      bond_constant = bond_force_constant[index3-1]
      equ_distance = bond_equ_distance[index3-1]
      energy_bond_exc = energy_bond_exc + bond_constant*(distance-equ_distance)**2
    

# Angular Energy :     
# Angles including Hydrogen :

    energy_angle_inc = 0
    temp_length = np.size(angles_inc_hyd)/4
    theta_inc_hyd = np.zeros((int(temp_length),1),dtype=float)
 
    for i in range(0,int(temp_length)):
      index1 = angles_inc_hyd[(i-1)*4]
      index2 = angles_inc_hyd[1+(i-1)*4]
      index3 = angles_inc_hyd[2+(i-1)*4]
      index4 = int(angles_inc_hyd[3+(i-1)*4])
      corrected_index1 = int(index1/3 + 1)
      corrected_index2 = int(index2/3 + 1)
      corrected_index3 = int(index3/3 + 1)
    
      x1 = X[corrected_index1-1]
      y1 = Y[corrected_index1-1]
      z1 = Z[corrected_index1-1]
    
      x2 = X[corrected_index2-1]
      y2 = Y[corrected_index2-1]
      z2 = Z[corrected_index2-1]
    
      x3 = X[corrected_index3-1]
      y3 = Y[corrected_index3-1]
      z3 = Z[corrected_index3-1]
    
      vector1 = np.array([(x1-x2),(y1-y2),(z1-z2)])
      vector2 = np.array([(x3-x2),(y3-y2),(z3-z2)])
      vector1 = vector1.transpose()
    
      mod1 = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
      mod2 = np.sqrt((x3-x2)**2 + (y3-y2)**2 + (z3-z2)**2)
    
      theta_inc_hyd[i] = np.dot(vector1,vector2)/(mod1*mod2)
      theta_inc_hyd[i] = math.acos(theta_inc_hyd[i])
      angle_constant = angle_force_constant[index4-1]
      equivalent_angle = equ_angle[index4-1]
    
      energy_angle_inc = energy_angle_inc + angle_constant*(theta_inc_hyd[i]-equivalent_angle)**2
    
    
# Angles excluding Hydrogen :  
    
    energy_angle_exc = 0
    temp_length = np.size(angles_without_hyd)/4
    theta_exc_hyd = np.zeros((int(temp_length),1),dtype=float)

    for i in range(0,int(temp_length)):
      index1 = angles_without_hyd[(i-1)*4]
      index2 = angles_without_hyd[1+(i-1)*4]
      index3 = angles_without_hyd[2+(i-1)*4]
      index4 = int(angles_without_hyd[3+(i-1)*4])
      corrected_index1 = int(index1/3 + 1)
      corrected_index2 = int(index2/3 + 1)
      corrected_index3 = int(index3/3 + 1)
    
      x1 = X[corrected_index1-1]
      y1 = Y[corrected_index1-1]
      z1 = Z[corrected_index1-1]
    
      x2 = X[corrected_index2-1]
      y2 = Y[corrected_index2-1]
      z2 = Z[corrected_index2-1]
    
      x3 = X[corrected_index3-1]
      y3 = Y[corrected_index3-1]
      z3 = Z[corrected_index3-1]
    
      vector1 = np.array([(x1-x2),(y1-y2),(z1-z2)])
      vector2 = np.array([(x3-x2),(y3-y2),(z3-z2)])
      vector1 = vector1.transpose()
    
      mod1 = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
      mod2 = np.sqrt((x3-x2)**2 + (y3-y2)**2 + (z3-z2)**2)
    
      theta_exc_hyd[i] = np.dot(vector1,vector2)/(mod1*mod2)
      theta_exc_hyd[i] = math.acos(theta_exc_hyd[i])
      angle_constant = angle_force_constant[index4-1]
      equivalent_angle = equ_angle[index4-1]
    
      energy_angle_exc = energy_angle_exc + angle_constant*(theta_exc_hyd[i]-equivalent_angle)**2
    
# Energy summed over Dihedrals : 

    temp_length = np.size(dihedrals_inc_hyd)/5    
    energy_dihedral_inc = 0
    dihedral_angle_inc = np.zeros((int(temp_length),1),dtype=float)

    for i in range(0,int(temp_length)):
      index1 = dihedrals_inc_hyd[(i-1)*5]
      index2 = dihedrals_inc_hyd[1+(i-1)*5]
      index3 = dihedrals_inc_hyd[2+(i-1)*5]
      index4 = dihedrals_inc_hyd[3+(i-1)*5]
      index5 = int(dihedrals_inc_hyd[4+(i-1)*5])
    
      if(index3<0 or index4<0):
          energy_dihedral_inc = energy_dihedral_inc 
        
      else :
          corrected_index1 = int(index1/3 + 1) 
          corrected_index2 = int(index2/3 + 1)
          corrected_index3 = int(index3/3 + 1)
          corrected_index4 = int(index4/3 + 1)
          x1 = X[corrected_index1-1]
          y1 = Y[corrected_index1-1]
          z1 = Z[corrected_index1-1]
    
          x2 = X[corrected_index2-1]
          y2 = Y[corrected_index2-1]
          z2 = Z[corrected_index2-1]
    
          x3 = X[corrected_index3-1]
          y3 = Y[corrected_index3-1]
          z3 = Z[corrected_index3-1]
        
          x4 = X[corrected_index4-1]
          y4 = Y[corrected_index4-1]
          z4 = Z[corrected_index4-1]
        
          vector1 = np.array([(x2-x1),(y2-y1),(z2-z1)])
          vector2 = np.array([(x3-x2),(y3-y2),(z3-z2)])
          vector3 = np.array([(x4-x3),(y4-y3),(z4-z3)])
          vector1 = vector1.transpose()
          vector2 = vector2.transpose()
          vector3 = vector3.transpose()
        
        
          plane_vector1 = np.cross(vector1,vector2)
          plane_vector2 = np.cross(vector2,vector3)
        
        
          mod1 = np.sqrt(plane_vector1[0,0]**2 + plane_vector1[0,1]**2 + plane_vector1[0,2]**2)
          mod2 = np.sqrt(plane_vector2[0,0]**2 + plane_vector2[0,1]**2 + plane_vector2[0,2]**2)
          plane_vector1 = plane_vector1.transpose()
        
          dihedral_angle_inc[i] = np.dot(plane_vector2, plane_vector1)/(mod1*mod2)
          dihedral_angle_inc[i] = math.acos(dihedral_angle_inc[i])
        
          dihedral_constant = dihedral_force_constant[index5-1]
          periodicity = dihedral_periodicity[index5-1]
          phase = dihedral_phase[index5-1]
        
          energy_dihedral_inc = energy_dihedral_inc + dihedral_constant*(1+math.cos(periodicity*dihedral_angle_inc[i] - phase))
        
        

# Dihedral Energy without hydrogen : 
        
    temp_length = np.size(dihedrals_without_hyd)/5    
    energy_dihedral_exc = 0
    dihedral_angle_exc = np.zeros((int(temp_length),1),dtype=float)

    for i in range(0,int(temp_length)):
      index1 = dihedrals_without_hyd[(i-1)*5]
      index2 = dihedrals_without_hyd[1+(i-1)*5]
      index3 = dihedrals_without_hyd[2+(i-1)*5]
      index4 = dihedrals_without_hyd[3+(i-1)*5]
      index5 = int(dihedrals_without_hyd[4+(i-1)*5])
     
      if(index3<0 or index4<0):
          energy_dihedral_exc = energy_dihedral_exc 
        
      else :
          corrected_index1 = int(index1/3 + 1) 
          corrected_index2 = int(index2/3 + 1)
          corrected_index3 = int(index3/3 + 1)
          corrected_index4 = int(index4/3 + 1)
          x1 = X[corrected_index1-1]
          y1 = Y[corrected_index1-1]
          z1 = Z[corrected_index1-1]
    
          x2 = X[corrected_index2-1]
          y2 = Y[corrected_index2-1]
          z2 = Z[corrected_index2-1]
    
          x3 = X[corrected_index3-1]
          y3 = Y[corrected_index3-1]
          z3 = Z[corrected_index3-1]
         
          x4 = X[corrected_index4-1]
          y4 = Y[corrected_index4-1]
          z4 = Z[corrected_index4-1]
          
          vector1 = np.array([(x2-x1),(y2-y1),(z2-z1)])
          vector2 = np.array([(x3-x2),(y3-y2),(z3-z2)])
          vector3 = np.array([(x4-x3),(y4-y3),(z4-z3)])
          vector1 = vector1.transpose()
          vector2 = vector2.transpose()
          vector3 = vector3.transpose()
         
        
          plane_vector1 = np.cross(vector1,vector2)
          plane_vector2 = np.cross(vector2,vector3)
        
        
          mod1 = np.sqrt(plane_vector1[0,0]**2 + plane_vector1[0,1]**2 + plane_vector1[0,2]**2)
          mod2 = np.sqrt(plane_vector2[0,0]**2 + plane_vector2[0,1]**2 + plane_vector2[0,2]**2)
          plane_vector1 = plane_vector1.transpose()
        
          dihedral_angle_exc[i] = np.dot(plane_vector2, plane_vector1)/(mod1*mod2)
          dihedral_angle_exc[i] = math.acos(dihedral_angle_exc[i])
        
          dihedral_constant = dihedral_force_constant[index5-1]
          periodicity = dihedral_periodicity[index5-1]
          phase = dihedral_phase[index5-1]
        
          energy_dihedral_exc = energy_dihedral_exc + dihedral_constant*(1+math.cos(periodicity*dihedral_angle_exc[i] - phase))     
        
        
# Lennard Jones Potentials : 

    energy_lennard_potential = 0
    energy_charge_interactions = 0

    for i in range(0,Natoms):
      for j in range(i+1,Natoms):
          atom_index = int(Ntypes*(Atom_type_index[i]-1) + Atom_type_index[j])
          index1 = int(NBparmindex[atom_index-1])
          Acoeff = lennard_a[index1-1]
          Bcoeff = lennard_b[index1-1]
          radius = np.sqrt((X[j]-X[i])**2 + (Y[j]-Y[i])**2 + (Z[j]-Z[i])**2 )
          
          energy_lennard_potential = energy_lennard_potential + Acoeff/(radius)**12 - Bcoeff/(radius)**6
          energy_charge_interactions = energy_charge_interactions + (Charge[i]*Charge[j])/radius
        
        
# Total Energy :

    Bond_energy = energy_angle_inc + energy_angle_exc
    Angle_energy = energy_angle_inc + energy_angle_exc
    Dihedral_energy = energy_dihedral_inc + energy_dihedral_exc

    Energy = Bond_energy + Angle_energy + Dihedral_energy + energy_lennard_potential + energy_charge_interactions
    print(Energy)
    return Energy                                                  
    
  
    
    

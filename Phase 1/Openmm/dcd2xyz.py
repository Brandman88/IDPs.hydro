import numpy as np
import MDAnalysis as mda

u = mda.Universe('Running_Config.pdb')
u.add_TopologyAttr('masses') # Mass attribute is not present: Add the masses attribute to the universe
Tag3d = u.select_atoms('all')
mass=3.0
Tag3d.masses = mass # Assign mass 

print ("Details of the Input XYZ Trajectory:")
print ("======================================")        
print (u.atoms)
print (u.trajectory)

with mda.Writer("Running_Config.xyz", Tag3d.n_atoms) as W:
    for ts in u.trajectory:
        W.write(Tag3d) 

print ("======================================")        
print("Conversion Complete !!!")

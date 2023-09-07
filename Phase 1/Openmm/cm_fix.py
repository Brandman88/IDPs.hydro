import MDAnalysis as mda
import pandas as pd

def read_multi_run(parameters='multi_run.dat'):
    with open(parameters,'r') as para:
        lines = []
        for line in para:
            lines.append(line)
    para.close()
    num_run_raw=lines[0]
    start_raw=lines[1]
    marker_raw=lines[2]
    
    if '#Number of proteins you care to run' in num_run_raw:
        num_run=int(num_run_raw.replace('#Number of proteins you care to run','').strip())
    else:
        num_run=int(num_run_raw)

    if '#Start value (if you are just starting the process should be on 1 by default)' in start_raw:
        start=int(start_raw.replace('#Start value (if you are just starting the process should be on 1 by default)','').strip())
    else:
        start=int(start_raw)

    if '#Put S or E based on starting or end process (don\'t really need to touch), if its I then set the other 2 numbers and run the file again.' in marker_raw:
        marker=marker_raw.replace('#Put S or E based on starting or end process (don\'t really need to touch), if its I then set the other 2 numbers and run the file again.','')
    else:
        marker=marker_raw
    marker=marker.strip()
    return num_run,start,marker

u = mda.Universe('Running_Config.pdb')
Poly3d = u.select_atoms('all')
ParticleN = len(u.atoms)
ntime = len(u.trajectory)

print ("Number of Particles:", ParticleN, "Trajectory Length:", ntime)

count = 0
fi = open ("Running_Config.xyz", "w")
for ts in u.trajectory:
    center_of_mass = Poly3d.center_of_mass()
    print(ParticleN, file = fi)
    print ("frame", count, file = fi)
    for i in range (ParticleN):
        cm_fixed_pos = Poly3d[i].position - center_of_mass
        print (Poly3d[i].name, cm_fixed_pos[0], cm_fixed_pos[1], cm_fixed_pos[2], file = fi)
    count = count + 1

fi.close()


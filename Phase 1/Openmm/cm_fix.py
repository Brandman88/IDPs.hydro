import MDAnalysis as mda

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


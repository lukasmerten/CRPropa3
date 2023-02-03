from crpropa import *
import numpy as np

Id = np.array([nucleusId(1, 1)]*3)
E = np.array([1, 2, 3])
lon = np.ones(3)
lat = np.zeros(3)
weight = np.array([0.5]*3)


M = ParticleMapsContainer()
#M.addParticles(Id, E, lon, lat, weight)
for i in range(3):
    #M.addParticle(Id[i], E[i], lon[i], lat[i], weight[i])
    M.addParticle(1, 1, 1, 1, 1)

print(M.getParticleIds())
print(M.getSumOfWeights())
#print(M.getMap(nucleusId(1,1), 1))
print(M.getEnergies(1))
import numpy as np

#dz_min = 20.0
#dz_max = 280.0
dz_min = 50.0
dz_max = 950.0
dz_strectch_zmin = 0.0
dz_strectch_zmax = 12000.0

# model top
zmax = 20000.0 



# Lowest level
z = dz_min * 0.5
k = 1
while(z <= zmax):
  print(z, k)

  if z < dz_strectch_zmin:
     dz = dz_min
  elif z >= dz_strectch_zmax:
     dz = dz_max
  else:
     print("")
  k += 1
  z += dz

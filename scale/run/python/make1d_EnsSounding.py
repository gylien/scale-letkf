import numpy as np
import sys
import os

def ens_pert(mems, SOUNDING):
   kmax = len(SOUNDING["z"])   

   amp = 2.0 # m/s Coffer et al. 2017MWR
   rand2d = np.zeros((kmax,2, mems)) # U & V

   for k in range(kmax):
     for i in range(2):
       rand = np.random.uniform(amp*(-1),amp,mems)
       rand -= np.average(rand) # ensure the mean is 0
       #print(rand2d.shape, rand.shape)
       #sys.exit()
       rand2d[k,i,:] = rand[:]
   return(rand2d)

def write_txt(outfile, SOUNDING, rand2d=None, m=None):
   print("Write ", outfile)
  
   with open(outfile, 'w') as f:
      f.write(SOUNDING["psfc"] + " " + SOUNDING["ptsfc"] + " "+ SOUNDING["qvsfc"] + "\n")

      for k in range(len(SOUNDING["z"])):
 
         if rand2d is not None:
           urand = rand2d[k,0,m]
           vrand = rand2d[k,1,m]
         else:
           urand = 0.0
           vrand = 0.0

         line = str(SOUNDING["z"][k]) + " " + str(SOUNDING["pt"][k]) + " " + str(SOUNDING["qv"][k]) + " " + str(SOUNDING["u"][k] + urand) + " " + str(SOUNDING["v"][k] + vrand) + "\n"
         f.write(line)

def read_txt(infile):
   print("Read ",infile)
  
   with open(infile, 'r') as f:
      data = f.read()

   data = np.array(data.split("\n"))

   z = np.zeros(len(data)-2)
   pt = np.zeros(len(data)-2)
   qv = np.zeros(len(data)-2)
   u = np.zeros(len(data)-2)
   v = np.zeros(len(data)-2)

   for l, line in enumerate(data):

      ldata = line.split()
      if l >= ( len(data) - 1): 
         break
      elif l == 0:
         psfc = ldata[0] # surface pressure (hPa)
         ptsfc = ldata[1] # surface potential temperature (K) 
         qvsfc = ldata[2] # surface water vapor mixing ratio (g/kg) 
      else:
         z[l-1] = ldata[0]  # height (m)
         pt[l-1] = ldata[1] # potential temperature (K)
         qv[l-1] = ldata[2] # water vapor mixing ratio (g/kg)
         u[l-1] = ldata[3]  # zonal wind (m/s)
         v[l-1] = ldata[4]  # meridional wind (m/s)


   SOUNDING = {"z":z, "pt":pt, "qv":qv, "u":u, "v":v, "psfc":psfc, "ptsfc":ptsfc, "qvsfc":qvsfc}

   return(SOUNDING)

def hodograph(top, mems):

  import matplotlib.pyplot as plt

  fig, (ax1) = plt.subplots(1, 1, figsize=(5.0, 5.0))
  fig.subplots_adjust(left=0.15, bottom=0.1, right=0.94, top=0.92)#, wspace=0.15, hspace=0.1)

  for m in range(mems+1):
     if m == 0:
        mem = "mean"
        lw = 4.0
        lc = 'k'
     else:
        mem = str(m).zfill(4)
        lw = 0.1
        lc = 'b'

     print("draw ",mem)
     input = os.path.join(top,mem + "_sounding.txt")
     SOUNDING = read_txt(input)
     plt.plot(SOUNDING['u'],SOUNDING['v'], color=lc, lw=lw)



  xmin = -22
  xmax = xmin * (-1)
  ymin = -22
  ymax = ymin * (-1)

  ax1.set_xlim(xmin,xmax)
  ax1.set_ylim(ymin,ymax)

  # draw circles
  theta = np.linspace(-np.pi, np.pi, 200)
  dR = 5
  for R in range(dR,dR*10,dR):
    if R % 10 == 0:
      ls = 'solid'
      lw = 2.0
    else:
      ls = 'dashed'
      lw = 1.0
    plt.plot(R*np.sin(theta), R*np.cos(theta), color='k', ls='dashed', 
             lw = lw, alpha = 0.2)

  ax1.axhline(y=0.0, xmin=xmin, xmax=xmax, color='k', ls= 'dotted')
  ax1.axvline(x=0.0, ymin=ymin, ymax=ymax, color='k', ls= 'dotted')

  plt.xlabel("Zonal wind", fontsize=14)
  plt.ylabel("Meridional wind", fontsize=14)

  fig.suptitle("Hodograph",fontsize=18)

  plt.show()

def main(top, mems):

  input = os.path.join(top,"mean_sounding.txt")

  SOUNDING = read_txt(input)
  rand2d = ens_pert(mems, SOUNDING)

  for m in range(mems):
    mem = str(m+1).zfill(4)

    outfile = os.path.join(top, mem + "_sounding.txt")
    write_txt(outfile, SOUNDING, rand2d=rand2d, m = m)


top = "/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/input/WK"
mems = 80

#main(top, mems)

hodograph(top, mems)


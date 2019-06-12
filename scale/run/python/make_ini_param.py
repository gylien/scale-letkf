import numpy as np
from datetime import datetime
import sys

def phys2func(phys,pmin,pmax):
    pmean = 0.5 * (pmax + pmin)

    tmp_phys = (phys - pmean) / (pmax - pmean)
    return(0.5 * (np.log((1.0 + tmp_phys) / (1.0 - tmp_phys))))

def func2phys(func,pmin,pmax):
    pmean = 0.5 * (pmax + pmin)

    return(np.tanh(func) * (pmax - pmean) + pmean)


#def func2phys(func,pmin,pmax):
#  return (0.5 * ((pmax + pmin) + np.tanh(func) * (pmax - pmin) ))
#
#
#
#def phys2func(phys,pmin,pmax):
#   var = (2.0 * phys - (pmax + pmin) ) / (pmax + pmin)
#   return(0.5 * (np.log(1.0 + var) - np.log(1.0 - var)))
#

# ---------

if __name__ == "__main__":

  PNAME = "Cs"
  PNAME = "drag_g"
  #PNAME = "Cr"

  ens = 100

  if PNAME == "Cs":
    pmax = 8.0
    #pmax = 1.7
    pmin = 0.1

    mean = 0.9 # default
    mean = 2.5 # default
    #mean = 4.0

    vari_func = 0.7 # function space
    #vari_func = 1.2 # function space

    vari_phys = 0.7
    #vari_phys = 0.6

  elif PNAME == "Cr":
    pmax = 180.0
    pmin = 10.0

    mean = 130.0
    vari_func = 0.8 # function space
    vari_phys = 15.0

    mean = 58.0
    vari_phys = 30.0
    vari_phys = 25.0


  elif PNAME == "drag_g":
    pmax = 5.0
    pmin = 0.1
    #mean = 0.6
    mean = 2.5
    #vari_phys = 1.0
    vari_phys = 0.7
    vari_func = 0.8 # function space


  PNAME_HEAD = "PARAM_ATMOS_PHY_MP_TOMITA08_"
  time = datetime(2013, 7, 13, 6, 0, 0)

  ctime = time.strftime('%Y%m%d%H%M%S')

  INIT_PHYS = True
  #INIT_PHYS = False


  mean_func = phys2func(mean,pmin,pmax)

  if INIT_PHYS:
    pens_phys = np.random.normal(mean,vari_phys,ens)

    #print(pens_phys)
    #print(np.mean(pens_phys))
    #print(mean)
    #sys.exit()

    if np.max(pens_phys) >= pmax or np.min(pens_phys) <= pmin:
      print(pens_phys)
      print("Failed to construct initial ensemble in phys")
      sys.exit()
  else:
    pens_func = np.random.normal(mean_func,vari_func,ens)


  if INIT_PHYS:
    ref_ens = pens_phys
  else:
    ref_ens = pens_func

  #for idx, func in enumerate(pens_func):
  for idx, func in enumerate(ref_ens):
    if INIT_PHYS:
      print(PNAME + ", ",str(idx+1) + ",",str(func) + "D0")
    else:
      print(PNAME + ", ",str(idx+1) + ",",str(func2phys(func,pmin,pmax)) + "D0")

  if INIT_PHYS:
    print(PNAME + ",mean,",str(mean) + "D0")
    print(PNAME + ",mdet,",str(mean) + "D0")
  else:
    print(PNAME + ",mean,",str(mean) + "D0")
    print(PNAME + ",mdet,",str(mean) + "D0")


  print("!#, Name, mean, pmin, pmax, vari_func, vari_phys, ens, INIT_PHYS")
  print("!#, " + PNAME + ",",str(mean),pmin,pmax,vari_func,vari_phys,ens,INIT_PHYS)

#########
#main()


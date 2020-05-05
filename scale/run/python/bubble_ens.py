import numpy as np




def main(INFO, of):

    rand_x = np.random.normal(0.0, INFO["SIGMA_LOC"], INFO["ENS"])
    rand_x -= np.mean(rand_x) 
    rand_x += INFO["CX"]

    rand_y = np.random.normal(0.0, INFO["SIGMA_LOC"], INFO["ENS"])
    rand_y -= np.mean(rand_y) 
    rand_y += INFO["CY"]

    rand_amp = np.random.normal(0.0, INFO["SIGMA_PT"], INFO["ENS"])
    rand_amp -= np.mean(rand_amp) 
    rand_amp += INFO["MEAN_PT"]

    f = open(of, "w")
    f.write("mean, " + str(INFO["CX"]) + "d0, " + str(INFO["CY"]) + "d0, " +
            str(INFO["MEAN_PT"]) + "d0 \n")

    for i in range(len(rand_x)):
       mem = str(i+1).zfill(4)
       if rand_amp[i] > INFO["MAX_PT"]:
          rand_amp[i] = INFO["MAX_PT"]
          print("MAX")
       elif rand_amp[i] < INFO["MIN_PT"]:
          rand_amp[i] = INFO["MIN_PT"]
          print("MIN")

       f.write(mem + ", " + str(rand_x[i]) + "d0, " + str(rand_y[i]) + "d0, " +
               str(rand_amp[i]) + "d0 \n")

       #print(mem, rand_x[i], rand_y[i], rand_amp[i])
       #f.write(mem,rand_x[i], rand_y[i], randa_amp[i])

    f.write("#" + str(INFO) + "\n") 
    f.close()

##########

#IMAX = 24
#JMAX = 24

IMAX = 20
JMAX = 20

PX = 8
PY = 8

ENS = 80
ENS = 320
CX = 2.e3 * IMAX * PX * 0.5
CY = 2.e3 * IMAX * PX * 0.5
SIGMA_LOC = 20.e3
SIGMA_PT = 1.0
# Similar to Aksoy et al. (2009MWR)
#MEAN_PT = 4.0
MEAN_PT = 3.0
MAX_PT = 6.5
MIN_PT = 1.5

MEAN_PT = 3.0
SIGMA_PT = 0.5
MAX_PT = 4.0
MIN_PT = 2.0

INFO = {"ENS":ENS, "CX":CX, "CY":CY, "SIGMA_LOC":SIGMA_LOC, "SIGMA_PT":SIGMA_PT,
        "MEAN_PT":MEAN_PT, "MAX_PT":MAX_PT, "MIN_PT":MIN_PT}

of = "init_loc_ens.txt"
of = "init_loc_ens320.txt"
of = "init_loc_ens320_0221.txt"
of = "init_loc_ens320_3K.txt"
of = "init_loc_ens320_3K_0414.txt"
main(INFO, of)



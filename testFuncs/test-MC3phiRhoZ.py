# -*- coding: utf-8 -*-

"""
test-MC3phiRhoZ.py
 
  Date      Who  Comments
==========  ===  =========================================================================
2018-01-23  JRM  Initial for ZnO

Elapse: 0:09:01.8 jrmFastMac

"""

import os
import sys
import glob
import math
import shutil
import time
import java.io as jio
import gov.nist.microanalysis.EPQLibrary as epq
import dtsa2.dt2Conv as dt2c
import dtsa2.mcSimulate3 as mc3

print(dt2c.__revision__)

start = time.time()

e0        =    15     # kV
lt        =   100     # sec
pc        =     5.0   # nA
nTraj     = 20000     # trajectories
rhoMat    =     5.606 # density g/cm3
naMat     = "ZnO"     # name of material
przDepUm  =     1.35  # phirhoz depth in microns
det       = findDetector("Oxford p4 05eV 2K")



homDir = os.environ['HOME']
wrkDir = homDir + "/Documents/git/dt2Conv/testFuncs"
csvDir = homDir + "/Documents/git/dt2Conv/data/csv"
dt2c.ensureDir(csvDir)
simDir = homDir + "/Documents/git/dt2Conv/data/sim/"
dt2c.ensureDir(simDir)

# print(wrkDir)
mat    = material(naMat, density=rhoMat)
pyrDir = wrkDir + "/test-MC3phiRhoZ Results"

DataManager.clearSpectrumList()

xrts = mc3.suggestTransitions(mat, e0)
print(xrts)
dose = pc * lt  # nA-sec"


xtraParams={}
xtraParams.update(mc3.configurePhiRhoZ(przDepUm*1.0e-6))
xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
xtraParams.update(mc3.configureOutput(simDir))

simStr = "MC simulation of bulk %s at %.1f keV + CSF + BSF" % (naMat, e0)

przFi = simDir + simStr + "/PhiRhoZ.csv"
# print(simStr)
# print(przFi)

sim = mc3.simulate(mat, det, e0, dose, True, nTraj, True, True, xtraParams)
fmtS = "%s-%g-kV-MC3"
sName = fmtS % (mat.name, e0)
sim.rename(sName)
sim.setAsStandard(mat)
sim.display()


csv = csvDir + "/%s-%.1f-keV-MC3-%d-traj-prz.tsv" % (naMat, e0, nTraj)
print(csv)
shutil.copyfile(przFi, csv)

shutil.rmtree(simDir+simStr)
# clean up cruft
shutil.rmtree(pyrDir)
print "Done!"

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg



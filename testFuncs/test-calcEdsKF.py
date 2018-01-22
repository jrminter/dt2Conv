# -*- coding: utf-8 -*-

"""
test-calcEdsKF.py
 
  Date      Who  Comments
==========  ===  =======================================================
2018-01-22  JRM  Test the calcEdsKF function

"""

import os, shutil, time
import dtsa2.dt2Conv as dt2c

print(dt2c.__revision__)

start = time.time()

homDir = os.environ['HOME']
wrkDir = homDir + "/Documents/git/dt2Conv/testFuncs"

os.chdir(wrkDir)
print(wrkDir)
pyrDir="./test-calcEdsKF Results"


det  = findDetector("Oxford p4 05eV 2K")
e0   = 15.     # keV
sp   = epq.SpectrumProperties()
ver  = False

sio2 = material("SiO2", density=2.65)
si   = material("Si", density=2.648)


print("PAP1991-Chantler2005 Si-SiO2")
print("\n")

resPAP = dt2c.calcEdsKF(si, det, e0, epq.PAP1991(),
                        epq.MassAbsorptionCoefficient.Chantler2005,
                        xtra=sp, stds={"Si": sio2}, verbose=ver)
print(resPAP)


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



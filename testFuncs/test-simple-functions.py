# -*- coding: utf-8 -*-

"""
test-simple-functions.py
 
  Date      Who  Comments
==========  ===  =========================================================================
2018-01-01  JRM  Print dt2Conv revision

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

start = time.time()

print(dt2c.__revision__)

homDir = os.environ['HOME']
wrkDir = homDir + "/Documents/git/dt2Conv/testFuncs"
spcDir = homDir + "/Documents/git/dt2Conv/data/Sirion"
csvDir = homDir + "/Documents/git/dt2Conv/data/csv"

os.chdir(wrkDir)
print(wrkDir)
pyrDir="./test-simple-functions Results"

DataManager.clearSpectrumList()

print("Test printTimeNow()")
dt2c.printTimeNow()

print("\nTest calcMAN()")
dt2c.calcMAN("0.1933*MgO+0.0927*Al2O3+0.4535*SiO2+0.1525*CaO+0.0996*FeO", "K412")

print("\nTest writeCompo()")
mat = dt2c.writeCompo("0.1933*MgO+0.0927*Al2O3+0.4535*SiO2+0.1525*CaO+0.0996*FeO", 2.600, "K412", True)

print("\nTest addMatToDatabase()")
kapton = epq.Material(epq.Composition([ epq.Element.C,
                                        epq.Element.O,
                                        epq.Element.N,
                                        epq.Element.H],
                                        [ 0.69113,
                                          0.20924,
                                          0.07327,
                                          0.02636 ]
                                          ),
                                          epq.ToSI.gPerCC(1.420))

dt2c.addMatToDatabase(kapton, "Kapton", 1.4200)

print("\nTest summarizeMaterial()")
out = dt2c.summarizeMaterial(kapton, 5)
print(out)

print("\nTest multiFilmBSE()")

al  = material("Al",  density=2.70)
c   = material("C",   density=2.267)
zno = material("ZnO", density=5.61)
si  = material("Si",  density=2.3296)

layers = [ [c,   20*1.0e-9],
           [zno, 20*1.0e-9],
           [al,  20*1.0e-9],
           [si,  50.0e-6]
         ]
sLayers = """
           layers = [
               [c,   20*1.0e-9],
               [zno, 20*1.0e-9],
               [al,  20*1.0e-9],
               [si,  50.0e-6]
           ]"""

print(sLayers)

a = dt2c.multiFilmBSE(layers, 7.0, nTraj=10000, outDir=csvDir)
print(a)

print("\nTest addCompositionsToDatabase()")



aupd = epq.Material(epq.Composition([epq.Element.Au, epq.Element.Pd],
                                    [0.60,0.40]),
                                    epq.ToSI.gPerCC(0.6*19.30+0.4*11.9))
aupd.setName("Au-Pd")
ta = material("Ta", density=16.4)

compoList = [aupd, ta]
dt2c.addCompositionsToDatabase(compoList)

print("added: [AuPd, Ta]")

print("\nTest getMassFractions()")
print("For: Al2Zn98O98")
elements = [epq.Element.Al, epq.Element.Zn, epq.Element.O]
massFra = dt2c.getMassFractions("Al2Zn98O98", elements, 5)
print(massFra)



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

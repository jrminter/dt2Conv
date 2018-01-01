# -*- coding: utf-8 -*-

"""
test-comp-KRs.py
 
  Date      Who  Comments
==========  ===  =========================================================================
2017-12-31  JRM  Test the compKRs function

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

homDir = os.environ['HOME']
wrkDir = homDir + "/Documents/git/dt2Conv/testFuncs"
spcDir = homDir + "/Documents/git/dt2Conv/data/7kV-calc"
os.chdir(wrkDir)
print(wrkDir)
pyrDir="./test-comp-KRs Results"

DataManager.clearSpectrumList()

e0  = 7 # kV
det = findDetector("Oxford p4 05eV 2K")
trs = [ epq.XRayTransitionSet(epq.Element.C, epq.XRayTransitionSet.K_FAMILY),
        epq.XRayTransitionSet(epq.Element.Si, epq.XRayTransitionSet.K_FAMILY)
      ]

cFile = spcDir + "/C-Std-7kV.msa"
cSpc = wrap(readSpectrum(cFile, i=0, det=det))
cSpc.display()

siFile = spcDir + "/Si-Std-7kV.msa"
siSpc = wrap(readSpectrum(siFile, i=0, det=det))
siSpc.display()


unFile = spcDir + "/10-nm-C-on-Si-7kV.msa"
unSpc = wrap(readSpectrum(unFile, i=0, det=det))
unSpc.display()

cStd  = {"El":element("C"),  "Spc":cSpc}
siStd = {"El":element("Si"), "Spc":siSpc}

stds = [cStd, siStd]

res = dt2c.compKRs(unSpc, stds, trs, det, e0)
print(res)


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



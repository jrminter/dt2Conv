# -*- coding: utf-8 -*-

"""
test-phiRhoZ.py
 
  Date      Who  Comments
==========  ===  =========================================================================
2018-01-01  JRM  Initial

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

print(dt2c.__revision__)

start = time.time()

homDir = os.environ['HOME']
wrkDir = homDir + "/Documents/git/dt2Conv/testFuncs"
csvDir = homDir + "/Documents/git/dt2Conv/data/csv"

print(wrkDir)

pyrDir = wrkDir + "/test-phiRhoZ Results"

DataManager.clearSpectrumList()

e0  = 15 # kV
det = findDetector("Oxford p4 05eV 2K")
zno = material("ZnO", density=5.606)


csv = csvDir + '/ZnO-prz-%g-kV.csv' % (e0)
ret = dt2c.computePhiRhoZ(zno, det, e0, csv, nSteps=200, usePAP=True)

csv = csvDir + '/ZnO-prz-Z-%g-kV.csv' % (e0)
ret = dt2c.computePhiRhoReportZ(zno, det, e0, csv,
                                nSteps=200, usePAP=True)

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



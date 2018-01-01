# -*- coding: utf-8 -*-

"""
test-strip-background.py
 
  Date      Who  Comments
==========  ===  =========================================================================
2017-12-31  JRM  Test the StripBackground function

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
spcDir = homDir + "/Documents/git/dt2Conv/data/Sirion"
os.chdir(wrkDir)
print(wrkDir)
pyrDir="./test-strip-background Results"

DataManager.clearSpectrumList()
det = findDetector("Oxford p4 05eV 2K")


cuFile = spcDir + "/Cu-2017-12-14-5eV-ch-2K.msa"
cuSpc = wrap(readSpectrum(cuFile, i=0, det=det))
cuSpc.display()

res = dt2c.StripBackground(cuSpc)
res.display()

print("Note the pulse-pileup peaks. Oxford exports the **raw spectral data**, not the pulse pileup corrected data...")


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



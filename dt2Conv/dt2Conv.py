# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

dt2Conv.py

A series of wrapper scripts to make DTSA-II automation easy

Place this file in DTSA_ROOT/lib/dtsa2/
call with:

import dtsa2.dt2Conv as dt2c

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2017-12-31  JRM  0.0.1   First test for DTSA-II Jupiter


"""

__revision__ = "$Id: dt2Conv.py John R. Minter 2018-01-01 $"
__version__ = "0.0.11"

import sys
import os
import glob
import shutil
import fnmatch
import time
import math
import csv
import codecs

sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.Utility as epu
import gov.nist.microanalysis.EPQTools as ept
import java.io as jio
import java.lang as jl
import java.util as ju

from java.lang import Double


import dtsa2 as dt2
import gov.nist.microanalysis.dtsa2 as gdtsa2
import dtsa2.mcSimulate3 as mc3

def printTimeNow():
    """printTimeNow()

    Print the current date and time

    Input
    -----
    None

    Return
    ------
    None

    Example
    -------
    import dtsa2.dt2Conv as dt2c
    dt2c.printTimeNow()
    """
    from datetime import datetime
    a = datetime.now()
    print(a.strftime("%A, %d. %B %Y %I:%M%p"))

def calcMAN(cf, name):
    """
    calcMAN(cf, name)

    Create a material object from a chemical formula and output
    the mean atomic number.

    Input
    -----
    cf: string
        The chemical formula (can also be a sum of oxides...)
    name: string
        name for the material

    Return
    ------
    none: prints the MAN to the command line

    Example
    -------
    import dtsa2.dt2Conv as dt2c

    dt2c.calcMAN("0.1933*MgO+0.0927*Al2O3+0.4535*SiO2+0.1525*CaO+0.0996*FeO", "K412")

    """
    mat = dt2.parseChemicalFormula(cf, 0, name)
    man = mat.meanAtomicNumber()
    line = "%s mean atomic number: %.3f" % (name, man)
    print(line)

def writeCompo(cf, density, name, norm=False):
    """
    writeCompo(cf, density, name, norm=False)

    Create a material from a string and print out the composition as
    mass fraction and atom fraction.

    Input
    -----
    cf: string
        The chemical formula (can also be a sum of oxides...)
    density: float
        The density [g/cm3]
    name: string
        name for the material
    norm: Boolean (False)
        Force normalization of mass fraction

    Return
    ------
    mat: DTSA material
        The material object

    Example
    -------
    import dtsa2.dt2Conv as dt2c

    mat = dt2c.writeCompo("0.1933*MgO+0.0927*Al2O3+0.4535*SiO2+0.1525*CaO+0.0996*FeO", 2.600, "K412", True)

    """

    mat = dt2.parseChemicalFormula(cf, density, name)
    if norm:
        mat.forceNormalization()
    se = mat.getSortedElements() # an array highest mf to lowest
    print("")
    print(name)
    print("element, mass-fraction, at-fraction")
    sumMF = 0.
    sumAF = 0.
    for el in se:
        line = "%s, " % (el)
        d0 = mat.weightFractionU(el, True)
        line += "%.6f, " % (round(d0.floatValue(), 6))
        sumMF += round(d0.floatValue(), 6)
        d1 = mat.atomicPercentU(el, True)
        line += "%.6f" % (round(d1.floatValue(), 6))
        sumAF += round(d1.floatValue(), 6)
        print(line)
    print("")
    line = "sum mass fractions = %.6f" % round(sumMF, 6)
    print(line)
    line = "sum atom fractions = %.6f" % round(sumAF, 6)
    print(line)
    print("")
    return mat

def addMatToDatabase(mat, name, density):
    """addMatToDatabase(mat, name, density)

    Add (or correct) an entry in the database.

    Parameters
    ----------

    mat - a DTSA material
        The material to add
    name - string
        The name to associate
    density - float
        the density in g/cm3

    Returns
    -------
    none. It just prints a message.

    Example
    -------
    import gov.nist.microanalysis.EPQLibrary as epq
    import dtsa2.dt2Conv as dt2c
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

    """
    mat.setName(name)
    mat.setDensity(epq.ToSI.gPerCC(density))
    dt2.Database.addStandard(mat)
    print("Added %s" % mat)

def summarizeMaterial(mat, iDigits=5):
    """"
    summarizeMaterial(mat, iDigits)

    A utility function to summarize a material

    Parameters
    ----------

    mat - a DTSA material
        The material to list
    iDigits - integer
        decimal places to round

    Returns
    -------
    A tuple: ( name, dictionary{Symbol : {mf : mass-fraction, af: atom-fraction}}, density)

    Example
    -------
    import gov.nist.microanalysis.EPQLibrary as epq
    import dtsa2.dt2Conv as dt2c
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

    kapton.setName("Kapton")

    out = dt2c.summarizeMaterial(kapton, 5)
    print(out)
    """
    density = round(mat.density/1000.0,iDigits) # in g/cc
    name = mat.name
    mf = {}
    elemList = mat.getElementSet()
    for el in elemList:
        wf = round(mat.weightFractionU(el, True).doubleValue(), iDigits)
        af = round(mat.atomicPercentU(el,True).doubleValue(), iDigits)
        mf[el.toAbbrev().encode('ascii','ignore')] = {"wf": wf, "af": af}

    rv = (name, mf, density)
    
    return(rv)

def multiFilmBSE(layers, e0=20.0, nTraj=1000, bOutFull=False,
                 outDir='C:/Temp/',fullOutFile='BSE.csv'):
    """multiFilmBSE(layers, det, e0=20.0, withPoisson=True,
                    nTraj=1000, dose=60.0, bOutFull=False,
                    outDir='C:/Temp/',fullOutFile='BSE.csv')

    Monte Carlo simulate the backscattered fraction from a multilayer
    thin film.

    Parameters
    ----------
    layers - an iterable list of [material,thickness].
        Note the materials must have associated densities.
    e0 - float (20)
        The accelerating voltage in kV
    nTraj - integer
        The number of trajectories to run
    bOutFull - boolean (False)
        A flag to write the full data to a csv file
    outDir - string ('C:/Temp/')
        Path for full output
    fullOutFile - string ('BSE.csv')
        name for full output file


    Returns
    -------
    tuple - (layers, e0, nTraj, backscatterFraction, 
             forwardscatterFraction)

    Example
    -------
    import dtsa2.dt2cen2 as dt2c
    al      = material("Al",  density=2.70)
    c       = material("C",   density=2.267)
    zno     = material("ZnO", density=5.61)
    si      = material("Si",  density=2.3296)

    layers = [ [c,   20*1.0e-9],
               [zno, 20*1.0e-9],
               [al,  20*1.0e-9],
               [si,  50.0e-6]
             ]

    a = dt2c.multiFilmBSE(layers, 7.0, nTraj=10000)
    print(a)

    """
    # myRP = dt2.reportPath()
    # myRP.replace("\\", "/")
    # print(myRP)
    monte = nm.MonteCarloSS()
    monte.setBeamEnergy(epq.ToSI.keV(e0))
    p = {}
    p["Layers"] = layers
    sr = monte.getChamber()
    origin = (0.0, 0.0, 0.0)
    pos = origin
    for (mat, thickness,) in layers:
        if thickness <= 0.0:
            raise "The layer thickness must be larger than zero."
        monte.addSubRegion(sr, mat, nm.MultiPlaneShape.createFilm([0.0, 0.0, -1.0], pos, thickness))
        pos = epu.Math2.plus(pos, [0.0, 0.0, thickness + 1.0e-12])
    bs0 = nm.BackscatterStats(monte, 100)
    monte.addActionListener(bs0)
    ann = nm.AnnularDetector(1.0e-3, 10, (0.0, 0.0, -1.0e-3), (0.0, 0.0, 1.0)) 
    monte.addActionListener(ann)
    monte.runMultipleTrajectories(nTraj)
    if bOutFull:
        fi = outDir + fullOutFile
        print(u"Results -> %s" % fi)
        fos = jio.FileOutputStream(fi)
        try:
            osw = jio.OutputStreamWriter(fos)
            osw.append("Parameters:\n")
            osw.append("E0\t%g keV\n" % e0)
            for k, v in p.iteritems():
                osw.append("%s\t%s\n" % (k, v))
            ann.dump(osw)
            osw.flush()
            bs0.dump(fos)
        finally:
            fos.close()
    return (layers, e0, nTraj, bs0.backscatterFraction(), bs0.forwardscatterFraction())




def addCompositionsToDatabase(compoList):
    """
    addCompositionsToDatabase(compoList)

    A utility function to add compositions to the DTSA-II database
    This is based on an exemplar from N. Ritchie.

    Parameters
    ----------
    compoList - List of type epq.Composition
        The list of compositions to be added to the database

    Returns
    -------
    None

    Example
    -------
    # Add Au/Pd and Ta coating materials to the database
    import dtsa2.dt2Conv as dt2c
    aupd = epq.Material(epq.Composition([epq.Element.Au, epq.Element.Pd],
                                        [0.60,0.40]),
                                        epq.ToSI.gPerCC(0.6*19.30+0.4*11.9))
    aupd.setName("Au-Pd")
    ta = material("Ta", density=16.4)

    compoList = [aupd, ta]
    dt2c.addCompositionsToDatabase(compoList)

    """
    for comp in compoList:
        ses = gdtsa2.DTSA2.getSession()
        if not ses.findStandard(comp.getName()):
            ses.addStandard(comp)

def getMassFractions(compound, elemList, iDigits):
    """
    getMassFractions(compound, elemList, iDigits)

    A utility function to compute the mass fractions for a compound

    Parameters
    ----------
    compound - string
        The stoichiometry as a molecular formula
    elemList = a list of type epq.element.symbol
        the elements in the compound
    iDigits - integer
        decimal places to round

    Returns
    -------
    A dictionary of {Symbol : mass-fraction}

    Example
    -------
    import dtsa2.dt2Conv as dt2c
    elements = [epq.Element.Al, epq.Element.Zn, epq.Element.O]
    massFra = dt2c.getMassFractions("Al2Zn98O98", elements, 5)
    """
    mat = dt2.material(compound)
    mf = {}
    for el in elemList:
        wf = round(mat.weightFractionU(el, True).doubleValue(), iDigits)
        mf[el.toAbbrev().encode('ascii','ignore')] = wf
    
    return(mf)



def hasProbeCurrent(spc):
    """hasProbeCurrent(spc)

    Check if the spectrum has a probe current

    Parameters
    ----------
    spc: ScriptableSpectrum
        he spectrum to check

    Returns
    -------
    res: Boolean
        True if the spectrum has a recorded probe current, False otherwise
    """
    sp = spc.getProperties()
    spStr = sp.toString()
    res = spStr.find("Probe current")
    if (res >= 0):
        ret = True
    else:
        ret = False
    return(ret)

def StripBackground(spc):
    """
    StripBackground(spc)

    Strip the background from a spectrum

    Parameters
    ----------
    spc: ScriptableSpectrum
        The spectrum to process
    det: DTSA detector instance (None)
        If a detector is passed, the detector is set

    Returns
    -------
    res: ScriptabelSpectrum
        The stripped spectrum
    """

    sp = spc.getProperties()
    lt = sp.getNumericProperty(epq.SpectrumProperties.LiveTime)
    hp = hasProbeCurrent(spc)
    if (hp == True):
        pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)

    nam = spc.getProperties().getTextProperty(epq.SpectrumProperties.SpectrumDisplayName)
    new =  nam + '-bks'
    spc = spc.applyLLT()
    spc = epq.SpectrumUtils.applyZeroPeakDiscriminator(spc)
    spc = epq.PeakStripping.Clayton1987.getStrippedSpectrum(spc)
    spc = dt2.wrap(spc)
    sp = spc.getProperties()
    if (hp == True):
        sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
        sp.setNumericProperty(epq.SpectrumProperties.FaradayEnd, pc)
    sp.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
    spc.rename(new)
    return spc

def compPeakIntegral(spc, ePeak, wid, digits=1, display=False):
    """compPeakIntegral( spc, ePeak, wid, digits=1, display=False)

    Compute the background-corrected peak interval for an peak

    Integrates around the center of a transition (or family)

    Parameters
    ----------
    spc: A DTSA-II scriptable spectrum
        The spectrum to process
    ePeak: float
        The centroid of the peak 
    wid: float
        The width (in eV) to integrate
    digits: integer (1)
        The number of digits to round the integral
    display: Boolean (False)
        Display the BKS spectrum if True

    Returns
    -------
    ret: A 2 element list [mean-integral-cps, std-dev-cps]

    Example
    -------
    import dtsa2.dt2Conv as dt2c
    spc = readSpectrum("path/to/CuSpc.msa")
    spc.display()
    pI = dt2c.compPeakIntegral(spc, digits=3)
    print(pi)


    """
    spc = StripBackground(spc)
    if(display==True):
        spc.display()
    ev = 1000*ePeak
    hw = 0.5*wid
    start = ev - hw
    end = ev + hw
    res = epq.SpectrumUtils.backgroundCorrectedIntegral(spc, start, end)
    r0 = round(res[0],digits)
    r1 = round(res[1],digits)
    # note: under the hood, this returns the estimate, res[0] 
    # and the uncertainty, res[1]. Typically it is one digit.
    # r0    = '%.0f' % res[0]
    # r1    = '%.0f' % res[1]
    ret = [r0, r1]
    return ret

def intCKandSiK(spc, digits=3, display=True):
    """intCkandSiK(spc, digits=3, display=True)

    Integrate the background-subtracted spectrum with C-K and Si-K peaks.

    Provides a useful way to estimate the C thickness.

    The function first zeros all counts below the Zero Peak Discriminator
    and does a background-corrected integrals:

     C: centered on 0.282  keV with a 180 eV width
    Si: centered on 1.7397 keV with a 240 eV width

    Parameters
    ----------
    spc: A DTSA-II scriptable spectrum
        The spectrum to process
    det: The detector associated with the spectrum
    digits: integer (3)
        The number of digits to round the integral
    display: Boolean (True)
        A flag to display the background corrected spectrum

    Returns
    -------
    ret: A dictionary where the key is the X-ray transition
        and the mean integral an d standard deviation in
        units of cps-per-nA if spectrum contains a probe
        current or units of cps if the probe current is missing.

    Example
    -------
    import dtsa2.dt2Conv as dt2c
    spc = readSpectrum("path/to/C-Si.msa")
    spc.display()
    pI = dt2c.intSiK(spc, digits=3, display=True)
    print(pi)
    """
    sp = spc.getProperties()
    lt = sp.getNumericProperty(epq.SpectrumProperties.LiveTime)
    e0 = sp.getNumericProperty(epq.SpectrumProperties.BeamEnergy)
    hp = hasProbeCurrent(spc)
    if (hp == True):
        pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)
        pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)

    cInt  = compPeakIntegral(spc, 0.282, 180, 1, display)
    siInt = compPeakIntegral(spc, 1.746, 240, 1, display)

    if (hp == True):
        out = {  "C": (round(cInt[0]/(lt*pc), digits),
                       round(cInt[1]/(lt*pc), digits)),
                "Si": (round(siInt[0]/(lt*pc), digits),
                       round(siInt[1]/(lt*pc), digits),
                       "cps-per-nA", e0, "kV") }
    else:
        out = {  "C": (round(cInt[0]/lt, digits),
                       round(cInt[1]/lt, digits)),
                "Si": (round(siInt[0]/lt, digits),
                       round(siInt[1]/lt, digits),
                       "cps", e0, "kV") }

    return out

def intCK(spc, digits=3, display=True):
    """intCK(spc, digits=3, display=True)

    Integrate the background-subtracted CK peak.

    Provides a useful way to estimate the probe current.

    The function first zeros all counts below the Zero Peak Discriminator
    and does a background-corrected integral centered on 0.282 keV with
    a 180 eV width

    Parameters
    ----------
    spc: A DTSA-II scriptable spectrum
        The spectrum to process
    det: The detector associated with the spectrum
    digits: integer (3)
        The number of digits to round the integral
    display: Boolean (True)
        A flag to display the background corrected spectrum

    Returns
    -------
    ret: A dictionary where the key is the X-ray transition
        and the mean integral an d standard deviation in
        units of cps-per-nA if spectrum contains a probe
        current or units of cps if the probe current is missing.

    Example
    -------
    import dtsa2.dt2Conv as dt2c
    spc = readSpectrum("path/to/Si.msa")
    spc.display()
    pI = dt2c.intCK(spc, digits=3, display=True)
    print(pi)
    """
    sp = spc.getProperties()
    lt = sp.getNumericProperty(epq.SpectrumProperties.LiveTime)
    e0 = sp.getNumericProperty(epq.SpectrumProperties.BeamEnergy)
    hp = hasProbeCurrent(spc)
    if (hp == True):
        pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)
        pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)

    cInt = compPeakIntegral(spc, 0.282, 180, 1, display)

    if (hp == True):
        out = { "C": (round(cInt[0]/(lt*pc), digits),
                      round(cInt[1]/(lt*pc), digits), "cps-per-nA",
                       e0, "kV") }
    else:
        out = { "C": (round(siInt[0]/lt, digits),
                      round(siInt[1]/lt, digits),
                       "cps", e0, "kV") }

    return out

def intSiK(spc, digits=3, display=True):
    """intSiK(spc, digits=3, display=True)

    Integrate the background-subtracted SiK peak.

    Provides a useful way to estimate the probe current.

    The function first zeros all counts below the Zero Peak Discriminator
    and does a background-corrected integral centered on 1.7397 keV with
    a 240 eV width

    Parameters
    ----------
    spc: A DTSA-II scriptable spectrum
        The spectrum to process
    det: The detector associated with the spectrum
    digits: integer (3)
        The number of digits to round the integral
    display: Boolean (True)
        A flag to display the background corrected spectrum

    Returns
    -------
    ret: A dictionary where the key is the X-ray transition
        and the mean integral an d standard deviation in
        units of cps-per-nA if spectrum contains a probe
        current or units of cps if the probe current is missing.

    Example
    -------
    import dtsa2.dt2Conv as dt2c
    spc = readSpectrum("path/to/Si.msa")
    spc.display()
    pI = dt2c.intSiK(spc, digits=3, display=True)
    print(pi)
    """
    sp = spc.getProperties()
    lt = sp.getNumericProperty(epq.SpectrumProperties.LiveTime)
    e0 = sp.getNumericProperty(epq.SpectrumProperties.BeamEnergy)
    hp = hasProbeCurrent(spc)
    if (hp == True):
        pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)
        pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)

    siInt = compPeakIntegral(spc, 1.746, 240, 1, display)

    if (hp == True):
        out = { "Si": (round(siInt[0]/(lt*pc), digits),
                       round(siInt[1]/(lt*pc), digits), "cps-per-nA",
                       e0, "kV") }
    else:
        out = { "Si": (round(siInt[0]/lt, digits), round(siInt[1]/lt, digits),
                       "cps", e0, "kV") }

    return out


def intCuK(spc, digits=3, display=True):
    """intCuK(spc, digits=3, display=True)

    Integrate the background-subtracted CuK peak.

    Provides a useful way to estimate the probe current.

    The function first zeros all counts below the Zero Peak Discriminator
    and does a background-corrected integral centered on 8.0478 keV with
    a 440 eV width

    Parameters
    ----------
    spc: A DTSA-II scriptable spectrum
        The spectrum to process
    digits: integer (3)
        The number of digits to round the integral
    display: Boolean (True)
        A flag to display the background corrected spectrum

    Returns
    -------
    ret: A dictionary where the key is the X-ray transition
        and the mean integral an d standard deviation in
        units of cps-per-nA if spectrum contains a probe
        current or units of cps if the probe current is missing.

    Example
    -------
    import dtsa2.dt2Conv as dt2c
    spc = readSpectrum("path/to/Cu.msa")
    spc.display()
    det = findDetector("Oxford p4 05eV 2K")
    pI = dt2c.intCuK(spc, det, digits=3, display=True)
    print(pi)
    """
    sp = spc.getProperties()
    lt = sp.getNumericProperty(epq.SpectrumProperties.LiveTime)
    e0 = sp.getNumericProperty(epq.SpectrumProperties.BeamEnergy)
    spStr = sp.toString()
    res = spStr.find("Probe current")
    if (res >= 0):
        pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)

    cuInt = compPeakIntegral(spc, 8.0478, 440, 1, display)

    if (res >= 0):
        out = { "Cu-K-L3": (round(cuInt[0]/(lt*pc), digits),
                       round(cuInt[1]/(lt*pc), digits), "cps-per-nA",
                       e0, "kV") }
    else:
        out = { "Cu-K-L3": (round(cuInt[0]/lt, digits), round(cuInt[1]/lt, digits),
                       "cps", e0, "kV") }

    return out

def intCuL(spc, digits=3, display=True):
    """intCuL(spc, digits=3, display=True)

    Integrate the background-subtracted CuL peak.

    Provides a useful way to estimate the probe current.

    The function first zeros all counts below the Zero Peak Discriminator
    and does a background-corrected integral centered on 0.894 keV with
    a 478 eV width

    Parameters
    ----------
    spc: A DTSA-II scriptable spectrum
        The spectrum to process
    digits: integer (3)
        The number of digits to round the integral
    display: Boolean (True)
        A flag to display the background corrected spectrum

    Returns
    -------
    ret: A dictionary where the key is the X-ray transition
        and the mean integral an d standard deviation in
        units of cps-per-nA if spectrum contains a probe
        current or units of cps if the probe current is missing.

    Example
    -------
    import dtsa2.dt2Conv as dt2c
    spc = readSpectrum("path/to/Cu.msa")
    spc.display()
    det = findDetector("Oxford p4 05eV 2K")
    pI = dt2c.intCul(spc, det, digits=3, display=True)
    """
    sp = spc.getProperties()
    lt = sp.getNumericProperty(epq.SpectrumProperties.LiveTime)
    e0 = sp.getNumericProperty(epq.SpectrumProperties.BeamEnergy)
    spStr = sp.toString()
    res = spStr.find("Probe current")
    if (res >= 0):
        pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)

    cuInt = compPeakIntegral(spc, 0.894, 478, 1, display)

    if (res >= 0):
        out = { "Cu-L": (round(cuInt[0]/(lt*pc), digits),
                       round(cuInt[1]/(lt*pc), digits), "cps-per-nA",
                       e0, "kV") }
    else:
        out = { "Cu-L": (round(cuInt[0]/lt, digits), round(cuInt[1]/lt, digits),
                       "cps", e0, "kV") }

    return out

def computePhiRhoZ(mat, det, e0, xrts, csv, nSteps=200, usePAP=True):
    """
    computePhiRhoZ(mat, det, e0, xrts, csv, nSteps=200, usePAP=True):
    Compute a phi-rho-z curve

    Parameters
    ----------
    mat: a material
        A composition to compute
    det: a detector
        The detector to use
    e0: float
        The beam energy in kV
    xrts: XRayTransitionSet
        the transitions to compute
    csv: string
        Path to .csv file
    nSteps: int (200)
        The number of steps to compute
    usePAP: boolean (True)
        PAP1991 if true otherwise XPP1991

    Returns
    -------
    None - but writes output to file csv

    """
    rhoZmax = epq.ElectronRange.KanayaAndOkayama1972.compute(mat, epq.ToSI.keV(e0))
    if usePAP:
        alg = epq.PAP1991()
    else:
        alg = epq.XPP1991()
    sp = epq.SpectrumProperties(det.getProperties())
    sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    res = "Idx,rhoz(mg/cm^2)"
    for xrt in xrts:
        res = "%s,G(%s),E(%s)" % (res, xrt, xrt)
    fi = open(csv, 'w')
    line = res + '\n'
    fi.write(line)

    for step in range(0, nSteps):
        rz = step * rhoZmax / nSteps
        res = "%d,%g" % (step, 100.0 * rz) # in mg/cm^2
        for xrt in xrts:
            alg.initialize(mat, xrt.getDestination(), sp)
            res = "%s,%g,%g" % (res, alg.computeCurve(rz), alg.computeAbsorbedCurve(xrt, rz))
        line = res+"\n"
        fi.write(line)
        # print(res)
    fi.close()

def computePhiRhoReportZ(mat, det, e0, xrts, csv, nSteps=200, usePAP=True):
    """
    computePhiRhoReportZ(mat, det, e0, xrts, csv, nSteps=200, usePAP=True)

    Compute a phi-rho-z curve and report a depth profile in Z
    The Z values are in microns.

    Parameters
    ----------
    mat: a material
        A composition to compute
    det: a detector
        The detector to use
    e0: float
        The beam energy in kV
    xrts: XRayTransitionSet
        the transitions to compute
    csv: string
        Path to .csv file
    nSteps: int (200)
        The number of steps to compute
    usePAP: boolean (True)
        PAP1991 if true otherwise XPP1991

    Returns
    -------
    None - but writes output to file csv

    """
    rhoZmax = epq.ElectronRange.KanayaAndOkayama1972.compute(mat, epq.ToSI.keV(e0))
    zMax = rhoZmax / mat.getDensity()
    if usePAP:
        alg = epq.PAP1991()
    else:
        alg = epq.XPP1991()
    
    sp = epq.SpectrumProperties(det.getProperties())
    sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    res = "Z(um)"
    for xrt in xrts:
        res = "%s,G(%s),E(%s)" % (res, xrt, xrt)
    fi = open(csv, 'w')
    line = res + '\n'
    fi.write(line)

    for step in range(0, nSteps):
        rz = step * rhoZmax / nSteps
        dz = step * zMax / nSteps
        res = "%g" % (1.0e6 * dz) # in um
        for xrt in xrts:
            alg.initialize(mat, xrt.getDestination(), sp)
            res = "%s,%g,%g" % (res, alg.computeCurve(rz), alg.computeAbsorbedCurve(xrt, rz))
        line = res+"\n"
        fi.write(line)
        # print(res)
    fi.close()

def getSpectrumFromDataCube(hsVecCube, det, spcFile, x, y, mapTime, pc=1.0, bDebug=False):
    """getSpectrumFromDataCube(hsVecCube, det, spcFile, x, y, mapTime, pc=1.0, bDebug=False)

    Get a calibrated spectrum from an open hyperspectral vector datacube
    (hsVecCube) which was typically opened from a .rpl file. Calibrate the
    energy using DTSA detector (det) and the sum spectrum file (spcFile).
    Live time is computed from the mapTime (sec)
    Extract the spectrum from the coordinates x (width) and y (height).
    The probe current (pc) defaults to 1.0. A bDebug flag (default False)
    presents diagnostic messages.

    Return a DTSA-II scriptable spectrum."""
    ss = dt2.readSpectrum(spcFile)
    props=ss.getProperties()
    props.setDetector(det)
    cw = ss.getChannelWidth()
    zo = ss.getZeroOffset()
    nc = ss.getChannelCount()
    if(bDebug):
        strDetProps = "Detector properties: cw %.5f, zo %.3f" % (cw, zo)
        print(strDetProps)
    
    ss.getProperties().setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
    ss.getProperties().setNumericProperty(epq.SpectrumProperties.FaradayEnd, pc)

    if(bDebug):
        display(ss)
    newSP = epq.SpectrumProperties()
    newSP.addAll(props);
    newSP.setTextProperty(epq.SpectrumProperties.SourceFile, "paint-sum-dtsa.msa")
    lt = newSP.getNumericWithDefault(epq.SpectrumProperties.LiveTime, Double.NaN)
    ### y,x in oxford map
    w = hsVecCube.getWidth()
    h = hsVecCube.getHeight()
    d = hsVecCube.getDepth()
    if ((x < w) and (y < h)):
        hsVecCube.seek(y,x) # Ba
        dat = hsVecCube.readDouble(d)
        sp = epq.SpectrumUtils.toSpectrum(cw, -zo, nc, dat)
        sp = dt2.wrap(sp)
        sp.setEnergyScale(newSP.getNumericProperty(epq.SpectrumProperties.EnergyOffset), newSP.getNumericProperty(epq.SpectrumProperties.EnergyScale))
        # props.setTimestampProperty(epq.SpectrumProperties.AcquisitionTime, ts)
        sp.setProbeCurrent(pc)
    
        lt = mapTime / (w*h)
        sp.setLiveTime(lt)
        newName = "location (%d,%d)" % (x,y)
        sp.rename(newName)
        dt2.display(sp)
    else:
        print("check x and y")
        sp = None
    return sp

def getSirionSiCpsPerNa(e0):
    """getSirionSiCpsPerNa(e0)

    Get the Si cps/nA for a desired kV.
    These values were determined for PT 6, 5 eV/ch, 2K channels
    in 2014-09-12-Cu-Si.

    Parameters
    ----------
    e0: double
        The beam energy in kV

    Returns
    -------
    val: double
        The cps/nA

    Example:
    import dtsa2.dt2Conv as dt2c
    a = dt2c.getSirionSiCpsPerNa(7.0)
    """
    val = 0.0
    if(e0 == 5.0):
        val = 1832.4
    if(e0 == 7.0):
        val = 3884.2
    if(e0 == 10.0):
        val = 7573.1
    if(e0 == 15.0):
        val = 13921.0
    if(e0 == 20.0):
        val = 20059.1
    if(e0 == 30.0):
        val = 29362.6
    return val
    
def getSirionCuCpsPerNa(e0):
    """getSirionCuCpsPerNa(e0)

    Get the Cu cps/nA for a desired kV.
    These values were determined for PT 6, 5 eV/ch, 2K channels
    in 2014-09-12-QC.

    Parameters
    ----------
    e0: double
        The beam energy in kV

    Returns
    -------
    val: double
        The cps/nA

    Example:
    import dtsa2.dt2Conv as dt2c
    a = dt2c.getSirionCuCpsPerNa(7.0)
    """
    val = 0.0
    if(e0 == 5.0):
        val = 2525.7
    if(e0 == 7.0):
        val = 3890.7
    if(e0 == 10.0):
        val = 5391.1
    if(e0 == 15.0):
        val = 6786.4
    if(e0 == 20.0):
        val = 7145.9
    if(e0 == 30.0):
        val = 6702.2
    return val

def reportComps(comps,path,fString="%s,%6.5f"):
    """reportComps(comps,path,,fString="%s,%6.5f")

    Outputs a collection of compositions (as weight fraction) to a .csv file.

    Parameters
    ----------
    comps: list
        The list of compositions to process
    path: string
        The path to a csv file to write
    fString: string ("%s,%6.5f")
        The format string for output

    Example:
    --------
    import dtsa2.dt2Conv as dt2c
    dt2c.reportComps(comps,"C:/Temp/fooWf.csv")
    """
    f = open(path, 'w')
    all = set()
    for comp in comps:
        all = all | set(comp.getElementSet())
    tmp = "Material"
    for elm in all:
        tmp = "%s,%s" % (tmp, elm.toAbbrev())
    f.write(tmp+"\n")
    for comp in comps:
        tmp = str(comp)
        for elm in all:
            tmp = fString % (tmp, comp.weightFraction(elm, False))
        f.write(tmp+"\n")
    # close the file
    f.close()
    
def reportAtmPct(comps,path):
    """reportAtmPct(comps,path)
    Outputs a collection of compositions (as atm pct) to a .csv file.

    Parameters
    ----------
    comps: list
        The list of compositions to process
    path: string
        The path to a csv file to write

    Returns
    _______
    None:   Does write a csv file


    Example:
    --------
    import dtsa2.dt2Conv as dt2c
    dt2c.reportAtmPct(comps,"C:/Temp/fooWf.csv")
    """
    f = open(path, 'w')
    all = set()
    for comp in comps:
        all = all | set(comp.getElementSet())
    tmp = "Material"
    for elm in all:
        tmp = "%s,%s" % (tmp, elm.toAbbrev())
    f.write(tmp+"\n")
    for comp in comps:
        tmp = str(comp)
        for elm in all:
            tmp = "%s,%6.5f" % (tmp, comp.atomicPercent(elm))
        f.write(tmp+"\n")
    # close the file
    f.close()

def isNaN(num):
    """isNaN(num)
    Check if a number is NaN, returning True of False"""
    return num != num

def checkNaN(x):
    """checkNaN(x)
    This checks if a value (e.g. K-ratio) is NaN and sets the value to
    zero if it is. This really helps when writing data frames to be
    read by R."""
    if isNaN(x):
        x = 0.0
    return x
    
def ensureDir(d):
    """ensureDir(d)
    Check if the directory, d, exists, and if not create it."""
    if not os.path.exists(d):
        os.makedirs(d)
        
def clearAllSpectra():
    """clearAllSpectra()
    Clear all spectra from the data manager."""
    DataManager = dt2.DataManager.getInstance()
    DataManager.clearSpectrumList()

def compKRs(unSpc, stds, trs, det, e0, digits=5):
    """compKRs(unSpc, stds, trs, det, e0, digits=5)

    Performs a MLLSQ filter-fit of unknown spectrum (unSpc),
    to the standard spectra in a list of dictionaries of standards.

    Parameters
    ----------
    unSpc: DTSA scriptable spectrum
        The unknown spectrum to analyze
    stds: a list of dictionaries of standards.
        See the example below
    trs: a list of transition sets
        See the example below
    det: a DTSA-II detector
        See the example below
    e0: float
        The beam energy (in kV)
    digits: int (5)
        The number of digits to round the k-ratios

    Example:
    --------
    import dtsa2.dt2Conv as dt2c

    e0  = 7 # kV
    det = findDetector("Oxford p4 05eV 2K")
    trs = [ epq.XRayTransitionSet(epq.Element.C, epq.XRayTransitionSet.K_FAMILY),
            epq.XRayTransitionSet(epq.Element.Si, epq.XRayTransitionSet.K_FAMILY)
          ]
    cFile = "/Path/C-std.msa"
    cSpc = wrap(readSpectrum(cFile, i=0, det=det))
    siFile = "/Path/Si-std.msa"
    siSpc = wrap(readSpectrum(siFile, i=0, det=det))
    unFile = "/Path/Unknown.msa"
    unSpc = wrap(readSpectrum(umFile, i=0, det=det))

    cStd  = {"El":element("C"),  "Spc":cSpc}
    siStd = {"El":element("Si"), "Spc":siSpc}

    stds = [cStd, siStd]

    res = dt2c.compKRs(unSpc, stds, trs, det, e0)
    print(res)
    """
    # Now set up the calc
    qa = epq.CompositionFromKRatios()
    ff=epq.FilterFit(det,epq.ToSI.keV(e0))
    l = len(stds)
    for i in range(0, l):
        st=stds[i]
        el=st["El"]
        sp=st["Spc"]
        ff.addReference(el,sp)
    krs=ff.getKRatios(unSpc)
    kr=[]
    n = len(trs)
    for i in range(0, n):
        tr=trs[i]
        k=round(checkNaN(krs.getKRatio(tr)), digits)
        kr.append(k)
    return kr

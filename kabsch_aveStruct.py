#!/usr/bin/env python


"""
This script takes in a pdb with multiple models, aligns them as specified using the Kabsch algorithm and calculates the average model, then spits out the RMSD between the average model and the original models.
Usage: ./kabsch_aveStruct.py [infile] [residue range] [atoms types]
[residue range] is the residue range of a given domain given in PyMOL format, i.e. low1-high1+low2-high2 etc., same goes for [atom types], eg. n+ca+c (write 'all' to include all atom types).

"""


__author__ = "Mateusz Maciejewski"
__email__ = "matt@mattmaciejewski.com"
__date__ = "December 2012"
__copyright__ = "Copyright (C) 2012 Mateusz Maciejewski"
__license__ = "Public Domain"


import sys
import re
import numpy
import datetime
from math import cos, sin, sqrt, acos, degrees, pi


if len(sys.argv) != 4:
    print >> sys.stderr, __doc__
    sys.exit(1)



class AtomCollector:

    """
    This class collects info on each atom in taret pdb.
    It creates a molecule object that stores atoms as attributes identified
    by atom numbers from pdb:
    .model - molecular model identifier
    .resid - type of residue the atom belongs to
    .chain - chain ID
    .residnum - sequence number of the residue
    .xcoor - x coordinate
    .ycoor - y coordinate
    .zcoor - z coordinate
    .occup - atom occupancy
    .tempfact - temperature factor
    .mass - atomic mass
    """

    def __init__(self):

        self.atom = {}
        self.atomNumber = {}
        self.resid = {}
        self.chain = {}
        self.residnum = {}
        self.xcoor = {}
        self.ycoor = {}
        self.zcoor = {}
        self.occup = {}
        self.tempfact = {}
        self.mass = {}
        self.model = {}
        self.atomCounter = 0
        whiteSpace = re.compile(r'\s+')

    def collect(self,fileToWorkOn):
                
        pattern1 = re.compile('^ATOM')
        pattern2 = re.compile('^MODEL')

        for line in fileToWorkOn.readlines():

            if pattern2.search(line):
                current_model = int(line.split()[1])
            
            if pattern1.search(line):
                
                    
                self.atomNumber[self.atomCounter] = line[6:11]
                self.atom[self.atomCounter] = line[12:16].strip()
                self.resid[self.atomCounter] = line[17:20].strip()
                self.chain[self.atomCounter] = line[21].strip()
                self.residnum[self.atomCounter] = int(line[22:26])
                self.xcoor[self.atomCounter] = float(line[30:38])
                self.ycoor[self.atomCounter] = float(line[38:46])
                self.zcoor[self.atomCounter] = float(line[46:54])
                self.occup[self.atomCounter] = line[54:60]
                self.tempfact[self.atomCounter] = line[60:66]

                try:
                    self.model[self.atomCounter] = current_model
                except NameError:
                    self.model[self.atomCounter] = 1
        

                if self.atom[self.atomCounter][0].lower() == 'c': self.mass[self.atomCounter] = 12.0
                elif self.atom[self.atomCounter][0].lower() == 'o': self.mass[self.atomCounter] = 16.0
                elif self.atom[self.atomCounter][0].lower() == 'n': self.mass[self.atomCounter] = 15.0
                elif self.atom[self.atomCounter][0].lower() == 's': self.mass[self.atomCounter] = 32.0
                elif self.atom[self.atomCounter][0].lower() == 'h': self.mass[self.atomCounter] = 1.0
                elif len(self.atom[self.atomCounter]) > 1 and self.atom[self.atomCounter][1].lower() == 'h': self.mass[self.atomCounter] = 1.0
                elif self.atom[self.atomCounter][0].lower() == 'd': self.mass[self.atomCounter] = 2.0
                elif self.atom[self.atomCounter][0].lower() == '2': self.mass[self.atomCounter] = 1.0
                elif self.atom[self.atomCounter][0].lower() == '4': self.mass[self.atomCounter] = 1.0
                elif self.atom[self.atomCounter][0].lower() == 'p': self.mass[self.atomCounter] = 31.0
                elif self.atom[self.atomCounter][0].lower() == 'z': self.mass[self.atomCounter] = 65.0
                elif self.atom[self.atomCounter][0].lower() == 'f': self.mass[self.atomCounter] = 55.9

            
                #self.mass[self.atomCounter] = 1.0 # uniform mass switch

                self.atomCounter += 1
        


def norm(vekt):
    normal = sqrt(vekt[0]**2+vekt[1]**2+vekt[2]**2)
    return vekt[0]/normal, vekt[1]/normal, vekt[2]/normal

def angle(a,b):
    return acos(a[0]*b[0]+a[1]*b[1]+a[2]*b[2])    

def rotate(rot1,rot2):
    rotation = numpy.matrix([(cos(rot1),(-1)*sin(rot1),0),
                           (sin(rot1),cos(rot1),0),
                           (0,0,1)])
    rotation = rotation * numpy.matrix([(cos(rot2),0,(-1)*sin(rot2)),
                                        (0,1,0),
                                        (sin(rot2),0,cos(rot2))])
    return rotation.T



def parseInput(resids, atoms):
    
    if resids != "all":
        residues = []
        if resids.find("+") != -1:
            resids = resids.split("+")

            for res in resids:
                if res.find("-") != -1:
                    residues += range(int(res.split("-")[0]),int(res.split("-")[1])+1)
                else:
                    residues += [int(res)]
        else:
            if resids.find("-") != -1:
                residues += range(int(resids.split("-")[0]),int(resids.split("-")[1])+1)
            else:
                residues += [int(resids)]
    else:
        residues = "all"

    if atoms != "all":
        atoms = atoms.split("+")

    return residues, atoms


def COMref(res, ats, ref_model, ref_molecule):

    COM_ref = [0,0,0]
    comdiv_ref = 0


    for entry in ref_molecule.atom.keys():
        if (int(ref_molecule.residnum[entry]) in res) and \
           (ref_molecule.atom[entry].lower() in ats or ats == "all") and \
           (ref_molecule.resid[entry] != "ANI") and ref_molecule.model[entry] == ref_model:

            COM_ref[0] += ref_molecule.xcoor[entry]*ref_molecule.mass[entry]
            COM_ref[1] += ref_molecule.ycoor[entry]*ref_molecule.mass[entry]
            COM_ref[2] += ref_molecule.zcoor[entry]*ref_molecule.mass[entry]
            comdiv_ref += ref_molecule.mass[entry]


    COM_ref[0] /= comdiv_ref
    COM_ref[1] /= comdiv_ref
    COM_ref[2] /= comdiv_ref


    for at in ref_molecule.atom.keys():

        ref_molecule.xcoor[at] -= COM_ref[0]
        ref_molecule.ycoor[at] -= COM_ref[1]
        ref_molecule.zcoor[at] -= COM_ref[2]

    return ref_molecule


def kabschAlign(res, ats, mobile_model, mobile_molecule, ref_model, ref_molecule):
    
    COM_mobile = [0,0,0]    #center of mass
    comdiv_mobile = 0

    atcoords_ref = []
    atcoords_mobile = []


    for entry in mobile_molecule.atom.keys():

        if (int(mobile_molecule.residnum[entry]) in res) and \
           (mobile_molecule.atom[entry].lower() in ats or ats == "all") and \
           (mobile_molecule.resid[entry] != "ANI") and mobile_molecule.model[entry] == mobile_model:

            COM_mobile[0] += mobile_molecule.xcoor[entry]*mobile_molecule.mass[entry]
            COM_mobile[1] += mobile_molecule.ycoor[entry]*mobile_molecule.mass[entry]
            COM_mobile[2] += mobile_molecule.zcoor[entry]*mobile_molecule.mass[entry]
            comdiv_mobile += mobile_molecule.mass[entry]

            atcoords_mobile.append([mobile_molecule.xcoor[entry], mobile_molecule.ycoor[entry], mobile_molecule.zcoor[entry]])


    COM_mobile[0] /= comdiv_mobile
    COM_mobile[1] /= comdiv_mobile
    COM_mobile[2] /= comdiv_mobile


    for entry in ref_molecule.atom.keys():
        if (int(ref_molecule.residnum[entry]) in res) and \
           (ref_molecule.atom[entry].lower() in ats or ats == "all") and \
           (ref_molecule.resid[entry] != "ANI") and ref_molecule.model[entry] == ref_model:

            atcoords_ref.append([ref_molecule.xcoor[entry], ref_molecule.ycoor[entry], ref_molecule.zcoor[entry]])


    L = len(atcoords_mobile)


    for i in range(len(atcoords_mobile)):
        atcoords_mobile[i][0] -= COM_mobile[0]
        atcoords_mobile[i][1] -= COM_mobile[1]
        atcoords_mobile[i][2] -= COM_mobile[2]
        


    #E0 = numpy.sum( numpy.sum(atcoords_ref * atcoords_ref, axis=0), axis=0) + numpy.sum( numpy.sum(atcoords_mobile * atcoords_mobile, axis=0), axis=0)

    V, S, Wt = numpy.linalg.svd( numpy.dot( numpy.transpose(atcoords_mobile), atcoords_ref))

    reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
 
    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]
 
    #RMSD = E0 - (2.0 * sum(S))
    #RMSD = numpy.sqrt(abs(RMSD / L))
 
    U = numpy.dot(V, Wt)
 
    for at in mobile_molecule.atom.keys():
        if mobile_molecule.model[at] == mobile_model:

            temp = [mobile_molecule.xcoor[at] - COM_mobile[0], mobile_molecule.ycoor[at] - COM_mobile[1], mobile_molecule.zcoor[at] - COM_mobile[2]]
            
            x, y, z = numpy.dot(temp, U)

            mobile_molecule.xcoor[at] = x
            mobile_molecule.ycoor[at] = y
            mobile_molecule.zcoor[at] = z
            
    return mobile_molecule


def RMSD(mol1, mod1, mol2, mod2, res, ats):

    RMSD = 0
    atnum = 0

    for i in mol1.atom.keys():
        if (int(mol1.residnum[i]) in res) and \
                (mol1.atom[i].lower() in ats or ats == "all") and \
                mol1.model[i] == mod1 and mol1.resid[i] != "ANI":

            for ii in mol2.atom.keys():
                if mol2.model[ii] == mod2 and \
                        mol2.atom[ii].replace(" ", "") == mol1.atom[i] and \
                        mol2.residnum[ii] == mol1.residnum[i]:
                        #mol2.atomNumber[ii] == mol1.atomNumber[i]:


                    dist = [mol1.xcoor[i] - mol2.xcoor[ii], mol1.ycoor[i] - mol2.ycoor[ii], mol1.zcoor[i] - mol2.zcoor[ii]]

                    RMSD += dist[0]**2 + dist[1]**2 + dist[2]**2
                    atnum += 1

    
    RMSD /= atnum
    RMSD = sqrt(RMSD)

    return RMSD


def main():


    fileToWorkOn = open(sys.argv[1],'r')
    id_mat = numpy.array([(1,0,0),(0,1,0),(0,0,1)])
    molecule = AtomCollector()
    molecule.collect(fileToWorkOn)
    fileToWorkOn.close()
    res, ats = parseInput(sys.argv[2], sys.argv[3])
    unique_models = list(set(molecule.model.values())) # this trick will remove the repetitions from values upon creation of the set, and then turn that into a list
    
    molecule = COMref(res, ats, 1, molecule)
        
    printRot = True
    calcAve = True
    
    for mod in unique_models:
        if mod != 1:
            molecule = kabschAlign(res, ats, mod, molecule, 1, molecule)
    

    if calcAve == True:

        aveFile = open(sys.argv[1].replace(".pdb", ".ave"), "w")

        for i in range(molecule.atomCounter):
            if molecule.model[i] == 1:
                ave_x = molecule.xcoor[i]; ave_y = molecule.ycoor[i]; ave_z = molecule.zcoor[i]

                q = 1

                for ii in range(molecule.atomCounter):
                    if molecule.atomNumber[ii] == molecule.atomNumber[i] and molecule.model[ii] != molecule.model[i]:

                        ave_x += molecule.xcoor[ii]
                        ave_y += molecule.ycoor[ii]
                        ave_z += molecule.zcoor[ii]
                        q += 1

                ave_x /= q
                ave_y /= q
                ave_z /= q

                print >> aveFile, "%-6s%5s %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6s%6s" % ("ATOM", molecule.atomNumber[i], molecule.atom[i], molecule.resid[i], molecule.chain[i], molecule.residnum[i], ave_x, ave_y, ave_z, molecule.occup[i], molecule.tempfact[i])

        print >> aveFile, "END"

        aveFile.close()


    calcRMSD = True # works only if the calcAve is True

    if calcRMSD == True:

        aveFile = open(sys.argv[1].replace(".pdb", ".ave"), "r")
        aveMol = AtomCollector()
        aveMol.collect(aveFile) # now the average file is loaded as a proper object
        aveFile.close()

        rms = []
        rms_ave = 0
        rms_sd = 0


        for i in unique_models:

            rms.append(RMSD(molecule, i, aveMol, 1, res, ats))
            rms_ave += rms[-1]

        rms_ave /= len(rms)


        for i in rms:
            rms_sd += (i - rms_ave)**2

        rms_sd /= len(rms)
        rms_sd = sqrt(rms_sd)

        print "RMSD from the fit of all models (on res %s and atoms %s) and their average: %.3f +/- %.3f" % (sys.argv[2], sys.argv[3], rms_ave, rms_sd)


    if printRot == True:

        rotFile = open(sys.argv[1].replace(".pdb", ".rot.pdb"), "w")

        for mod in unique_models:

            if mod == 1:

                print >> rotFile, """REMARK Produced by %s on %s,
REMARK a Kabsch algorithm-based structural alignment tool.
REMARK %d models aligned using %s atoms
REMARK and residues in %s range.
REMARK The RMSD of this alignment is %f +/- %f.
REMARK For information on the algorithm implementation visit
REMARK www.mattmaciejewski.com, or email me at matt@mattmaciejewski.com
""" %(sys.argv[0], datetime.datetime.now().strftime("%H:%m:%S, %d %h %Y"), len(unique_models), sys.argv[3], sys.argv[2], rms_ave, rms_sd)


            print >> rotFile, "MODEL      %3s" %mod
            
            for i in range(molecule.atomCounter):
                if molecule.model[i] == mod:
                    if len(molecule.atom[i]) < 4:
                        print >> rotFile, "%-6s%5s  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f%6s%6s" % ("ATOM", molecule.atomNumber[i], molecule.atom[i], molecule.resid[i], molecule.chain[i], molecule.residnum[i], molecule.xcoor[i], molecule.ycoor[i], molecule.zcoor[i], molecule.occup[i], molecule.tempfact[i])
                    else:
                        print >> rotFile, "%-6s%5s %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6s%6s" % ("ATOM", molecule.atomNumber[i], molecule.atom[i], molecule.resid[i], molecule.chain[i], molecule.residnum[i], molecule.xcoor[i], molecule.ycoor[i], molecule.zcoor[i], molecule.occup[i], molecule.tempfact[i])
            print >> rotFile, "ENDMDL"
        print >> rotFile, "END"
                
        rotFile.close()


if __name__ == "__main__":

    main()

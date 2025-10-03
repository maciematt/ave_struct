ave_struct
==========

This script takes in a pdb with multiple models, aligns them as specified using the Kabsch algorithm and calculates the average model, then spits out the RMSD between the average model and the original models.

I wrote this script while calculating structures of multidomain proteins using CYANA. I wanted to know the domain-wise rmsds and couldn't find any software that did that. At first I tried modifying Cameron Mura's script for PyMOL, but the PyMOL library is horribly slow at this task (the bottleneck being the calculation of the mean structure relative to which the rmsds of the structures at hand are to be calculated).

What ended up being a good solution was performing the same task using numpy to avoid the PyMOL overhead. The resulting script aligns models within a pdb on the portion of choice via the Kabsch algorithm. A lot of this code (and the inspiration to write it) was based on the Kabsch script written by Jason Vertrees, available on the PyMOL wiki.


Requires: numpy


Usage: `./ave_struct.py structure.pdb n+ca+c 15-52`

Where `structure.pdb` is the structure with multiple models to be aligned, `n+ca+c` is the PyMOL-formatted types of atoms used in the alignment (this can be any atom range, like `ca`, or `n+ca+c+o` for MOLMOL-style backbone heavy atoms), and `15-52` is the the PyMOL-formatted residue range used in the alignment (so it could be multiple ranges linked with a + sign, like `15-52+67-82`).


Output: 

`structure.ave` - the average structure file based on the minimal RMSD alignment of the models constituing `structure.pdb`.
`structure.rot.pdb` - a version of `structure.pdb` with all models aligned.
stdout - the RMSD of the alignment (averaged for the models of `structure.pdb`) relative to `structure.ave` is printed in stdout.

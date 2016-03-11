################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""
====================================
Pipeline pipeline_immunogenicity.py
====================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python


Overview
========

A pipeline to predict the immunogenicity i.e. MHC binding potential of a set of 
peptides. It uses tools from the IEDB epitope prediction resource http://www.iedb.org/.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

Input
-----

Input is a fasta file of peptide sequences. The fasta file must end in .fasta or .fasta.gz


Pipeline output
===============

The output is a table of peptides and the associated predicted IC50 values for
MHC binding. The IC50 is the minimum IC50 for the protein i.e. the 15mer that
has the strongest affinity.

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *

import sys
import os
import glob
import gzip
import itertools
import collections
import time
import optparse
import shutil
import sqlite3
from rpy2.robjects import r as R
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as Database
import CGAT.FastaIterator as FastaIterator
import numpy as np

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGATPipelines.Pipeline as P
P.getParameters( 
    [ "pipeline.ini" ] )

PARAMS = P.PARAMS

@follows(mkdir("fasta.dir"))
@split("*.fa.gz", "fasta.dir/*.fasta")
def splitFasta(infile, outfiles):
    '''
    split fasta file into separate files
    '''
    for fasta in FastaIterator.iterate(IOTools.openFile(infile[0])):
        filename = fasta.title.replace(" ", "_").split("_")[0]
        outf = IOTools.openFile(os.path.join(
            "fasta.dir", filename + ".fasta"), "w")
        outf.write(">%s\n%s\n" % (fasta.title, fasta.sequence))
    outf.close()

###################################################
###################################################
###################################################

@follows(mkdir("prediction.dir"))
@transform(splitFasta, regex("(\S+)/(\S+).fasta"), r"prediction.dir/\2.pred")
def predictEpitope(infile, outfile):
    '''
    make a prediction of IC50 of binding to MHC based on 
    IEDB resources
    '''
    allele=PARAMS.get("prediction_allele")
    method=PARAMS.get("prediction_method")
    tempfile=os.path.join(os.path.dirname(outfile), os.path.basename(infile)) + ".tmp"

    statement="""/usr/local/bin/python2.7 /home/nilott/apps/src/mhc_ii/mhc_II_binding.py \
                                       %(method)s \
                                       %(allele)s \
                                       %(infile)s \
                                       > %(tempfile)s""" % locals()
    os.system(statement)
    peptide2data = {}
    
    # store ic50 values
    ics = []
    res = IOTools.openFile(tempfile)
    header = res.readline()
    seq_id = P.snip(os.path.basename(infile), ".fasta")
    for line in res.readlines():
        data = line[:-1].split("\t")
        allele, seq_num, start, end, core_peptide, peptide, ic50 = data
        peptide2data[peptide] = [seq_num, start, end, core_peptide, peptide, ic50]
        ics.append(ic50)

    # get the peptide with the minimum predicted IC50
    min_ic50 = min(ics)
    outf = open(outfile, "w")
    outf.write("%s\tseq_id\n" % header[:-1])
    for peptide, dat in peptide2data.iteritems():
        ic50 = dat[-1]
        if ic50 != min_ic50: 
            continue
        else:
          outf.write("\t".join([peptide] + dat + [seq_id]) + "\n") 
    outf.close()
    os.unlink(tempfile)


###################################################
###################################################
###################################################

@merge(predictEpitope, "prediction.dir/agg_epitopes.tsv.gz")
def aggregatePredictions(infiles, outfile):
    '''
    concatenate the individual predictions
    '''
    infs = " ".join(infiles)
    statement = """awk 'FNR==1 && NR !=1 { while (/^allele/) getline; } 1 {print}' prediction.dir/*.pred | gzip > %s""" % outfile
    os.system(statement)



if __name__== "__main__":
    sys.exit( P.main(sys.argv) )    




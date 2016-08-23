'''
predictions2minimumic50.py
============================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Take a table of predicted IC50 values and output the sequence/epitope
with the minimum.

.. Overall purpose and function of the script>

Usage
-----

.. Example use case

Example::

   python predictions2minimumic50.py

Type::

   python predictions2minimumic50.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import os

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """
    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-n", "--seq-name", dest="seq_name", type="string",
                      help="supply sequence name")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    peptide2data = {}
    
    # store ic50 values
    ics = []
    res = sys.stdin
    header = res.readline()
    seq_id = options.seq_name
    for line in res:
        data = line.strip("\n").split("\t")
        allele, seq_num, start, end, core_peptide, peptide, ic50 = data
        peptide2data[peptide] = [allele, seq_num, start, end, core_peptide, peptide, float(ic50)]
        ics.append(float(ic50))
        
    # get the peptide with the minimum predicted IC50
    min_ic50 = min(ics)
    
    options.stdout.write("%s\tseq_id\n" % header[:-1])
    for peptide, dat in peptide2data.iteritems():
        ic50 = dat[-1]
        if ic50 != min_ic50: 
            continue
        else:
          options.stdout.write("\t".join(map(str,dat) + [seq_id]) + "\n") 

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

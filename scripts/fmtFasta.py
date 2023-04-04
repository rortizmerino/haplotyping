#!/usr/bin/python

# to deal with all sequence-related steps
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# regular expressions
import re

# to call system functions as sys.exit
import sys

# to get command line arguments
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile",
                    help="input file")
parser.add_argument("-o", "--outdir",
                    help="output directory")

# read arguments from the command line
args = parser.parse_args()

# check for no arguments and print usage
if len(sys.argv) == 1 :
 sys.exit("\n usage: python fmtFasta.py <--infile> <--outdir> <--help>\n") 

fileFields = re.split('/',args.infile)
outfile = f"{args.outdir}/{fileFields[-1]}"

lst2print = []

# read genomic fasta 
with open (args.infile) as input_handle:
 for record in SeqIO.parse(input_handle, "fasta"):
  
  # to write output fasta
  ID = record.id
  Desc = record.description
  NTs = record.seq

  new_record = SeqRecord(NTs, id=ID, description=Desc)
  lst2print.append(new_record)

with open(outfile, "w") as output_handle:
 SeqIO.write(lst2print, output_handle, "fasta")

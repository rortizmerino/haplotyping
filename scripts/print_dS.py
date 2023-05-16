#!/usr/bin/python

# to expand *
import glob

# to run external scripts / commands
import subprocess

# to iterate over files in directory
import os

# for regular expressions
import re

# to deal with all sequence-related steps
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# to call system functions as sys.exit
import sys

# file/path operations
from pathlib import Path

# to get command line arguments
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directory",
                    help="input directory e.g. Scer_vs_S1")

# read arguments from the command line
args = parser.parse_args()

# check for no arguments and print usage
if len(sys.argv) == 1 :
 sys.exit("\n usage: python3 print_dS.py <--directory> <--help>\n") 

# check for -directory
root_dir = args.directory
path = Path(root_dir)
if path.exists():
 msg = f'\n {root_dir} directory exists'
 print(msg)
else:
 msg = f'\n {root_dir} does not exist'
 print(msg)
 sys.exit(" check directory\n")

out_list = []

pathlist = Path(root_dir).glob('orthology/orth*/M*/*dS')
for path in pathlist:
 # because path is object not string
 path_in_str = str(path)

 # define here so its restarted on each file
 orths_list = []
 with open(path_in_str) as in_file:
  for line in in_file:
   # collapse multiple spaces to 1
   line = re.sub("\s+", " ", line)
   # ignore 1st line starting with spaces (paml orth count)
   if re.search("^(?! ).*", line):
    fields = line.split(" ")
    orths_list.append(fields[0])
   last_one = line

  dSvals = last_one.split(" ")

# all in 1 line
  # define here so its restarted on each file
#  out = ""
#  i = 0
#  for ID in reversed(orths_list):
#  if i == 0:
#    out = ID + f'\t'
#    i = i + 1
#   else:
#    out = out + ID + f'\t' + dSvals[i] + f'\t'
#    i = i + 1
  # save records in list
#  out_list.append(out)
# all in 1 line
 
  # define here so its restarted on each file
  out0 = ""
  out1 = ""
  i = 0
  for ID in reversed(orths_list):
   if i == 0:
    out0 = ID + f'\t'
    i = i + 1
   else:
    out1 = out0 + ID + f'\t' + dSvals[i]
    i = i + 1
    out_list.append(out1) 

#outfile = re.sub("/", "", root_dir)
fileflds = root_dir.split("/")
outfile = fileflds[2]
outfile = outfile + f'_dS.txt'
msg = f'\n printing into {outfile}\n'
print(msg)
with open(outfile, 'w') as f:
 for line in sorted(out_list):
  line = re.sub("-T1", "", line)
  f.write("%s\n" % line)

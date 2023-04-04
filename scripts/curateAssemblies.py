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

# to check if file exists
from pathlib import Path

# to get command line arguments
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directory",
                    help="input directory e.g. fltrdSortedPhased")
parser.add_argument("-i", "--infile",
                    help="input file e.g. phasedAssembliesStats.txt")

# read arguments from the command line
args = parser.parse_args()

# check for no arguments and print usage
if len(sys.argv) == 1 :
 sys.exit("\n usage: python3 curateAssemblies.py <--infile> <--directory> <--help>\n") 

# check for -infile
if args.infile:
 print(" infile set to: %s" % args.infile)
else:
 args.infile = "phasedAssembliesStats.txt"
 print(" infile set to default: %s" % args.infile)
infile_name = args.infile

# read lists to save cov, circ, repeat, mult, len, from Flye
all_IDs_list = []
all_IDs_hash = {}

with open(infile_name) as in_file:
 for line in in_file:
  # ignore header
  if re.search("^(?!#).*", line):
   fields = line.split("\t")
   ID = f"{fields[1]}"
   sampleFields = re.split('_|\t',fields[0])
   sample = sampleFields[0]

   tmp = f"{fields[3]} {fields[4]} {fields[5]} {fields[6]} {fields[2]}"
   all_IDs_hash[ID] = tmp
   all_IDs_list.append(ID)

path = args.directory
#ext = ('phunph.fa')
if re.search("canu", infile_name):
 ext = ('r2cat.canu.fasta')
if re.search("flye", infile_name):
 ext = ('r2cat.flye.fasta')

for file in os.listdir(path):
 if file.endswith(ext):

  fullfile = f"{path}/{file}"

  records = list(SeqIO.parse(fullfile, "fasta"))
  outfilename = file.replace(".fasta",".ids_n_stats.txt", ) # 2nd comma?
  fulloutfile = f"{path}/{outfilename}"

  outfile = open(f"{path}/fastANIlist.txt", "w")
  for i in range(len(records)):
   # to write fasta and run ANI
   ID = records[i].id
   NTs = records[i].seq
   new_record = SeqRecord(NTs, id=ID, description="")
   filename = f"{path}/{ID}.fasta"
   with open(filename, "w") as output_handle:
    SeqIO.write(new_record, output_handle, "fasta")

   # to write list
   outstr = f"{path}/{records[i].id}.fasta\n"
   outfile.write(outstr)
  # close file so it can be passed to subprocess
  outfile.close()

  # run ANI
  args = ["fastANI", "--queryList", f"{path}/fastANIlist.txt", "-r", "references/genomes/Scen.fa", "--fragLen", "1000", "--output", f"{path}/fastANI.txt", "--threads", "6"]
  p = subprocess.Popen(args)
  p.wait()

  # capture ANI output
  outfile2 = open(fulloutfile, "w")
  with open(f"{path}/fastANI.txt") as in_file:
   for line in in_file:
    # ignore header
    fields = line.split("\t")
    tmp = fields[0].replace(".fasta","", )
    contig = tmp.replace(f"{path}/","", )
    outstr = f"{contig} {all_IDs_hash[contig]} {fields[2]}\n"
    outfile2.write(outstr)
  outfile2.close()

 else:
  continue

exit()

# cleanup
# Get a list of all the file paths that ends with .txt from in specified directory
fileList = glob.glob('*.fa')
fileList.append("fastANIlist.txt")
fileList.append("fastANI.txt")
# Iterate over the list of filepaths & remove each file.
for filePath in fileList:
 try:
  os.remove(filePath)
 except:
  print("Error while deleting file : ", filePath)

# manually edit .ids.txt into curated.ids.txt
# e.g. keep only scaffolds from unph corresponding to missing chr regions in ph according to dotplot and cov
# special care about "repeat" and "mult" info from phasedAssembliesStats

# iterate over curated id txt files in curated directory

exit()

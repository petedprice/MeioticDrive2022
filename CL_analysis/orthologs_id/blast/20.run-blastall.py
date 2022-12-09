#!/usr/bin/python2.6
# -*- coding: utf-8 -*-

''' Run BLASTALL
This script automatizes a complete reciprocal besthit blast analysis.
A list of input genomes is blasted against each other and output is 
written in provided format. Name of outfile is db name.query name.bla.
Run with "6 qseqid sseqid pident length mismatch gapopen qstart 
qend sstart send evalue bitscore sseq"
'''
#==============================================================================
import argparse
import sys
import os
import time
from subprocess import Popen, list2cmdline
from itertools import product
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("genomes", type=str, nargs="+",
                    help="A list of input genomes in fasta format")
parser.add_argument("-e", "--evalue", type=float, default=1e-7,
                    help="E-value for blast [default 1e-7]")
parser.add_argument("-f", "--outfmt", type=str,
                    help="The output format produced by blast")
parser.add_argument("-p", "--processors", type=float,
                    help="Maximum number of processors [default: all]")
parser.add_argument("-id", "--identity", type=int, default=25,
                    help="min. percent identity [default: 25]")
parser.add_argument("-b", "--blast_type", type=str, default="blastp",
                    help="blastprogram {blastp} [default: blastp]")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def cpu_count():
    ''' Returns the number of CPUs in the system
    '''
    num = 1
    if sys.platform == 'win32':
        try:
            num = int(os.environ['NUMBER_OF_PROCESSORS'])
        except (ValueError, KeyError):
            pass
    elif sys.platform == 'darwin':
        try:
            num = int(os.popen('sysctl -n hw.ncpu').read())
        except ValueError:
            pass
    else:
        try:
            num = os.sysconf('SC_NPROCESSORS_ONLN')
        except (ValueError, OSError, AttributeError):
            pass

    return num

def exec_commands(cmds):
    ''' Exec commands in parallel in multiple process
    (as much as we have CPU). It uses cpu_count to check for the number
    of available processors and runs as many jobs in parallel. 
    '''
    if not cmds:
        return  # empty list

    def done(p):
        return p.poll() is not None

    def success(p):
        return p.returncode == 0

    def fail():
        sys.exit(1)
    # Check how many processors we have
    max_task = cpu_count()
    processes = []
    # While loop submits commands to precessors.
    while True:
        while cmds and len(processes) < max_task:
            task = cmds.pop()
            print list2cmdline(task)
            processes.append(Popen(task))

        for p in processes:
            if done(p):
                if success(p):
                    processes.remove(p)
                else:
                    fail()

        if not processes and not cmds:
            break
        else:
            time.sleep(0.05)

def exec_in_row(cmds):
    ''' Exec commands one after the other until finished. This is helpful
    if a program is already parallelized, but we have to submit many jobs
    after each other.'''
    if not cmds:
        return  # empty list

    def done(p):
        return p.poll() is not None

    def success(p):
        return p.returncode == 0

    def fail():
        sys.exit(1)

    for task in cmds:
        p = Popen(task)
        p.wait()

    if done(p):
            if success(p):
                print "done!"
            else:
                fail()

def set_seq_type(blast_type):
    '''Determines the sequence type based on the supplied blast command'''
    if blast_type == "blastp":
        return "prot"
    elif blast_type == "blastn":
        return "nucl"
    else:
        raise Exception("We cannot guess the sequence type!")

#==============================================================================
#Main==========================================================================
#==============================================================================

def main():
    '''Main runs only when the file itself is executed.
    If a part of the script is imported, main will not be executed.'''
    seq_type = set_seq_type(args.blast_type)
    
    # Lists for the created makeblastdb and blast commands
    makeblastdb = []
    makeblastrun = []

    # Creates the makeblastdb commands for each genome supplied.
    for g in args.genomes:
        makeblastdb.append(["makeblastdb", "-input_type", "fasta",
                                           "-dbtype", seq_type,
                                           "-in", g])
    for pair in product(args.genomes, repeat=2):
        if pair[0] != pair[1]:
            db = pair[0]
            query = pair[1]
            db_basename = os.path.basename(db)
            query_basename =os.path.basename(query)
            ''' Name of file is db then query'''
            blastout = db_basename+"."+query_basename+".bla"
            makeblastrun.append([args.blast_type, "-evalue", str(args.evalue),
                                                  "-task", args.blast_type,
                                                  "-outfmt", args.outfmt,
                                                  "-db", db,
                                                  "-query", query,
                                                  "-out", blastout])
    
    # Passes the makeblastdb commands (a list of lists) to the exec_commands
    # function and executes them in parallel.
    exec_commands(makeblastdb)
    
    # Passes the makeblastrun commands (a list of lists) to the exec_commands
    # function and executes them in parallel.
    exec_commands(makeblastrun)

if __name__ == "__main__":
    main()

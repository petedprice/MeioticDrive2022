import argparse
import sys
import os
from itertools import combinations
import cPickle as pickle


a = "data/blast_dbs/stalkie_longest.dros_longest.bla"
b = "data/blast_dbs/dros_longest.stalkie_longest.bla"
infiles=(a,b)
outfile = "del"

def find_pairs_ensembl(infiles):
    '''Finds files in a list that have the form a <--> b and b <--> a
    Where a and b are different species. Each pair will be used to determine
    the reciprocal_besthit later. Returns list of pairs
    '''
    pairs = []
    for a, b in combinations(infiles, 2):
        abase = os.path.basename(a).split(".")
        bbase = os.path.basename(b).split(".")
        a_db = abase[0].split("_")[0]
        a_query = abase[1].split("_")[0]
        b_db = bbase[0].split("_")[0]
        b_query = bbase[1].split("_")[0]
        if a_db == b_query and a_query == b_db:
            pairs.append((a,b))
    return pairs

def read_blastout(source):
    countdeleted = 0
    try:
        with open(source, "r") as infile:
            sbase = os.path.basename(source).split(".")
            s_db = sbase[0].split("_")[0]
            s_query = sbase[1].split("_")[0]
            count =0
            blasthits = {}
            blastlist = []
            currentscaffold = None
            for line in infile:
                count +=1
                qseqid = line.rstrip().split()[0]
                if currentscaffold is None:
                    currentscaffold = qseqid
                    blastlist.append(line)
                else:
                    if currentscaffold == qseqid:
                        blastlist.append(line)
                    else:
                        blasthits, countdel = read_blastout_bitscore(blastlist, blasthits, sbase)
                        countdeleted = countdeleted + countdel
                        blastlist = []
                        currentscaffold = qseqid
                        blastlist.append(line)
            #final loop for final scaffold
            blasthits, countdel = read_blastout_bitscore(blastlist, blasthits, sbase)
            countdeleted = countdeleted + countdel
            print source
            print "Number of lines =", count
            print "Number of blast hits filtered because no one tophit =", countdeleted
        return blasthits
    except IOError:
        print "!----ERROR----!"
        print "File %s does not exit!" % source
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(1)

def read_blastout_bitscore(blastlist, blasthits, sbase):
    ''' The blastout output format supported is:
    "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq"
    For each hit only the besthit is stored, rest is discarded.
    Minimum 30% identity, then pick tophit with greatest bitscore. If bitscores are identical
    pick the tophit with greatest %identity. If equal %identity then ortholog discarded.
    '''
    if len(blastlist) == 0:
        line = line.rstrip().split()
        s_db = sbase[0].split("_")[0]
        s_query = sbase[1].split("_")[0]
        qseqid = s_query+"_"+line[0]
        sseqid = s_db+"_"+line[1]
        pidentity = float(line[2])
        bitscore = float(line[11])
        blasthits[qseqid] = (sseqid, bitscore, pidentity)
    else:
        countdel = 0
        check = []
        for line in blastlist:
            line = line.rstrip().split()
            s_db = sbase[0].split("_")[0]
            s_query = sbase[1].split("_")[0]
            qseqid = s_query+"_"+line[0]
            sseqid = s_db+"_"+line[1]
            pidentity = float(line[2])
            bitscore = float(line[11])
            if pidentity > 30: # Check for min 30% identity
                if qseqid in blasthits:
                    if blasthits[qseqid][1] < bitscore:
                        blasthits[qseqid] = (sseqid, bitscore, pidentity)
                    elif blasthits[qseqid][1] == bitscore:
                        if blasthits[qseqid][2] < pidentity:
                            blasthits[qseqid] = (sseqid, bitscore, pidentity)
                        elif blasthits[qseqid][2] == pidentity:
                            info = [bitscore, pidentity]
                            check.append(info)
                else:
                    blasthits[qseqid] = (sseqid, bitscore, pidentity)
        #check that there arent two hits which have the highest bitscore and pidentity but are identical
        if float(len(check)) > 0:
            sortcheck = sorted(check, key=lambda t: t[0], reverse=True)
            bitscore = float(sortcheck[0][0])
            pidentity = float(sortcheck[0][1])
            if bitscore == blasthits[qseqid][1]:
                if pidentity == blasthits[qseqid][2]:
                    countdel += 1
                    del blasthits[qseqid]
    return blasthits, countdel

def reciprocal_besthit(a,b):
    '''Finds the tophit in each file and determines if they are the same.
    Identifies and returns list of 1-1 reciprocal orthologs.
    '''
    recip_besthit = []
    a_blastout = read_blastout(a)
    b_blastout = read_blastout(b)
    for g in a_blastout:
        query = g
        query_besthit = a_blastout[g][0]
        if query_besthit in b_blastout:
            target_besthit = b_blastout[query_besthit][0]
            if query == target_besthit:
                recip_besthit.append((query, query_besthit))
    return recip_besthit

pairs = find_pairs_ensembl(infiles)
print len(pairs)
print pairs
recip_besthit_all = []
for a, b in pairs:
    #pick tophit
    recip_besthit = reciprocal_besthit(a, b)
    #combine all tophit dictionaries into one
    recip_besthit_all.append(recip_besthit)
print "Number of pairs = ", len(recip_besthit_all)

print "Starting to pickle..."
with open(outfile+".pkl", "w") as outfile:
    pickle.dump(recip_besthit_all, outfile, -1)

print "Writing to file ..."
count = 0
with open(outfile+".txt", "w") as outtext:
    for gene in recip_besthit_all:
        for g in gene:
            count += 1
            outtext.write(str(g))
            outtext.write("\n")
print "Number of reciprocal tophits =", count
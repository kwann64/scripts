import sys
import csv
import random
import statistics

if len(sys.argv) < 2:
    print ('Usage: python3 mpileup_hets.py [basename]; calculates observed heterozygosit for each individual and each population based on a [basename].tped and a [basename].mpileup file from SAMtools; requires a "pops" no-header tsv file that has individual names in col-1 and pop names in col-2; email kevinwwann@gmail.com if you have issues')
    sys.exit(1)

# load in mpileup file and keep genotype column with no additional characters
genos = []
with open(sys.argv[1]+'.mpileup','r') as f1: # modify this to your .mpileup file path
    csv_reader = csv.reader(f1, delimiter='\t')
    for line in csv_reader:
        templist = []
        cleaned = []
        bases = 'ACGT'
        for x in range(4,len(line),3):
            templist.append(line[x])
        for val in templist:
            geno = ''
            for char in val.upper():
                if char in bases:
                    geno += char
            cleaned.append(geno)
        genos.append(cleaned)

# subsample two bases at random from sites with at least 2x coverage
genos_sub = []
for site in genos:
    templist = []
    for geno in site:
        if len(geno) > 1:
            templist.append(''.join(random.sample(geno,2)))
        else:
            templist.append('')
    genos_sub.append(templist)

# create dictionaries for {sample name:het count}, and for {mpileup position:sample name}
hets = {}
samp_IDs = {}
with open('output/'+sys.argv[1]+'.tfam','r') as f1: # modify this for your TFAM path
    csv_reader = csv.reader(f1, delimiter=' ')
    count = 0
    for line in csv_reader:
        hets[line[0]] = 0
        samp_IDs[count] = line[0]
        count += 1

# for each sample, loop through that position in the mpileup file and count the total number of sites with data and the number of het states
for x in range(0,len(hets)):
    het_count = 0
    total_count = 0
    count = 0
    for site in genos_sub:
        count += 1
        if len(site[x]) == 2:
            total_count += 1
            if site[x][0] != site[x][1]:
                het_count += 1
    hets[samp_IDs[x]] = het_count/total_count

# read in a 'pops' file with each sample and the ancestral group they belong to
indi_pops = []
with open("pops",'r') as myfile:
    csv_reader = csv.reader(myfile, delimiter='\t')
    for row in csv_reader:
        indi_pops.append(row)

# build dictionaries for {pop:[list of indis in pop]} and {pop:[list of het values for each indi in pop]}
popdict = {}
pophets = {}
for a in indi_pops:
    if a[1] not in list(popdict.keys()):
        popdict[a[1]] = []
        pophets[a[1]] = []

for a in indi_pops:
    if a[0] in hets.keys():
        popdict[a[1]].append(a[0])

for key in popdict.keys():
    for indi in popdict[key]:
        pophets[key].append(hets[indi])

# print each population's average heterozygosity and each individual's heterozygosity
for key in pophets.keys():
    print(key+': '+str(statistics.mean(pophets[key])))

for a in indi_pops:
    print(a[0]+','+a[1]+','+str(hets[a[0]]))

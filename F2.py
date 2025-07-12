import sys
import csv
import subprocess

if len(sys.argv) < 2:
    print ('Usage: python3 F2.py [basename]; calculates F2 sites shared within and between pre-defined groups; requires a "pops" no-header tsv file that has individual names in col-1 and pop names in col-2; email kevinwwann@gmail.com if you have issues')
    sys.exit(1)


######## READ IN AND CONVERT TPED TO WORKABLE LIST
tped = []
with open("../variants/output/"+sys.argv[1]+".tped",'r') as myfile: ## modify this to your TPED path
    csv_reader = csv.reader(myfile, delimiter='\t')
    for row in csv_reader:
        templist = []
        templist1 = []
        templist.append(row[4:])
        for x in range(0,len(templist[0])-1):
            if x % 2 != 0:
                templist1.append(templist[0][x])
        tped.append(templist1)

######## CALCULATE AND SAVE THE F2 SITES
baselist = ['A','C','G','T']
F2sites = []
F2counts = []
for line in tped:
    bases = {'A':0,'G':0,'C':0,'T':0}
    for state in line:
        if state in baselist:
            bases[state] += 1
    if bases['A'] == 2 or bases['C'] == 2 or bases['G'] == 2 or bases['T'] == 2:
        F2sites.append(line)
        F2counts.append(bases)

####### SET UP INDIVIDUAL AND POPULATION DICTIONARIES
indi_pops = []
with open("pops",'r') as myfile:
    csv_reader = csv.reader(myfile, delimiter='\t')
    for row in csv_reader:
        indi_pops.append(row)
allpops = []
for i in indi_pops:
    allpops.append(i[1])
pops = list(set(allpops))
pop_dict = dict.fromkeys(pops)
pop_counts = dict.fromkeys(pops)
pop_index = dict.fromkeys(pops)
for i in pop_dict:
    pop_dict[i] = []
    pop_counts[i] = []
    for x in range(0,len(pops)):
        pop_counts[i].append(0)
for i in range(0,len(indi_pops)):
    for x in pops:
        if indi_pops[i][1] == x:
            pop_dict[x].append(i)
for x in range(0,len(pops)):
    pop_index[pops[x]] = x
indis={}
for x in range(0,len(F2sites[0])):
    indis[x] = 0

####### IDENTIFY THE INDIVIDUALS INVOLVED IN EACH F2 SITE AND COUNT UP HITS BETWEEN EACH POP
for a,b in zip(F2sites,F2counts):
    if b['A'] == 2:
        q = 'A'
    elif b['C'] == 2:
        q = 'C'
    elif b['G'] == 2:
        q = 'G'
    else:
        q = 'T'

    hits = []
    for x in range(0,len(a)):
        if a[x] == q:
            indis[x] += 1
            hits.append(x)

    hitpops = []
    for i in pops:
        for x in hits:
            if x in pop_dict[i]:
                hitpops.append(i)

    pop_counts[hitpops[0]][pop_index[hitpops[1]]] += 1
    pop_counts[hitpops[1]][pop_index[hitpops[0]]] += 1

for i in pop_counts:
    for x in range(0,len(pop_counts[i])):
        pop_counts[i][x] = str(pop_counts[i][x])
    pop_counts[i].append(i)
    print(','.join(pop_counts[i]))

## HELPER FUNCTIONS ##

def preambule(family):
    s = "# STOCKHOLM 1.0\n#=GF ID    "+family
    return s

def structure_annotation(gapped_seq, ss_cons):
    ungapped_seq = ""
    ungapped_ss = ""

    pairing_dict = pairing_dictionary(ss_cons)

    for k in range(len(gapped_seq)):
        if gapped_seq[k]!='-':
            # in any case, we add the seq symbol
            ungapped_seq += gapped_seq[k]

            # depending on whether the position
            # is paired or not, and if the partner
            # is a gap if it is paired, we
            # change what we add to the structure
            add_ss_cons_char = k in pairing_dict.keys()
            add_ss_cons_char = add_ss_cons_char and  gapped_seq[pairing_dict[k]]!='-'
            add_ss_cons_char = add_ss_cons_char and compatible(gapped_seq[k],gapped_seq[pairing_dict[k]])
            # all 3 conditions must be verified:
            #    - the position is paired
            #    - the position it is paired with is not a gap
            #    - the position it is paired with is compatible
            if add_ss_cons_char:
                ungapped_ss += ss_cons[k]
            else:
                ungapped_ss += '.'

    return ungapped_seq, ungapped_ss
            
def compatible(x,y):
    return (x,y) in [('A','U'),
                     ('U','A'),
                     ('G','U'),
                     ('U','G'),
                     ('G','C'),
                     ('C','G')]

def pairing_dictionary(ss_cons):
    symbols = ['()','[]','{}','<>','Aa','Bb','Cc','Dd','Ee','Ff']

    # return dict
    d = {}

    stacks = {}
    # oc is for "opening/closing"
    for oc in symbols:
        stacks[oc] = []

    # classic stack structure to extract bps
    for k,c in enumerate(ss_cons):
        for oc in symbols:
            if c==oc[0]:
                stacks[oc].append(k)
            if c==oc[1]:
                j = stacks[oc].pop()
                d[k] = j
                d[j] = k

    return d

## END HELPER FUNCTIONS ##
## BEGINNING OF MAIN SCRIPT ##

ss_cons = {}

pk_list = open('pseudoknotted_families.txt','w')

# first pass on the file: extracting ss_cons for each family.
for line in open('Rfam.seed', encoding='ISO-8859-1').readlines():
    if line.startswith('#=GF AC'):
        family = line.split(' ')[-1].rstrip('\n')
        print(family)
    if line.startswith('#=GC SS_cons'):
        ss_cons_family = line.split(' ')[-1].rstrip('\n')
        ss_cons[family] = ss_cons_family
        if ss_cons_family.find('A') >=0:
            print(family, file=pk_list)


# second pass on the file: actually making annotated alignments per family
for line in open('Rfam.seed', encoding='ISO-8859-1').readlines():
    if line.startswith('#=GF AC'):
        family = line.split(' ')[-1].rstrip('\n')
        print(family)
        current_file = open('annotated_alignments/'+family+'.stk','w')
        print(preambule(family), file=current_file)

#    if family!='RF04300':
#        continue

    if line!='\n' and not line.startswith('#') and line!='//\n':
        gapped_seq = line.split(' ')[-1].rstrip('\n')

        # purpose is to keep the spaces after identifier
        identifier = line.split(gapped_seq)[0] 

        ungapped_seq, ungapped_ss = structure_annotation(gapped_seq,ss_cons[family])

        ss_annotation_start = '#=GR'+' '*(len(identifier)-4)

        print(identifier+ungapped_seq, file=current_file)
        print(ss_annotation_start+ungapped_ss, file=current_file)

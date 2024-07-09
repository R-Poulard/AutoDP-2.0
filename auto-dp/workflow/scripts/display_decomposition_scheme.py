from utils import read_td_lines
from colors import *
from DPTable import DPTable, compute_dp_table
colors = Set3

# extracting bags and tree from td file
root = open(snakemake.input.tdname).readlines()[0].split(' ')[1].rstrip('\n')
adj, index2bag = read_td_lines(open(snakemake.input.tdname).readlines())


# extracting helix extemities information 
cluster_extremities = {}
for helixline in open(snakemake.input.helix).readlines():
    label = helixline.split(' ')[0]
    extremities = [c.replace(' ','') for c in helixline.split('(')[1].split(')')[0].split(',')]
    cluster_extremities[label] = extremities
    
dp_tables = compute_dp_table(open(snakemake.input.tdname).readlines(),
                             open(snakemake.input.helix).readlines(), 
                             open(snakemake.output[0],'w'))

from autodp.minimal_expansion import MinimalExpansion

# setting options
inter_helix_gap = snakemake.config["inter_helix_gap"]
overarching = True # overarching 1 <--> n (extremities) arc.

# object init
minimal_expansion = MinimalExpansion()

# retrieve input
dbn = open(snakemake.input[0]).readlines()[0].rstrip('\n')

# converting to graph of minimal expansion
minimal_expansion.from_str(dbn, 
                           inter_helix_gap=inter_helix_gap, 
                           overarching=overarching)

# dumping to file
minimal_expansion.dump_to_gr(snakemake.output[0])

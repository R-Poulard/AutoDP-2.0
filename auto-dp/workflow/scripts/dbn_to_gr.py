from autodp.minimal_expansion import MinimalExpansion

# setting options
inter_helix_gap = True # a Boolean stating whether
                                                      # helices share extremities
                                                      # or not
overarching = True # overarching 1 <--> n (extremities) arc.

# object init
minimal_expansion = MinimalExpansion()

# retrieve input
dbn = open('resources/dbn_files/K5.dbn').readlines()[0].rstrip('\n')

# converting to graph of minimal expansion
minimal_expansion.from_str(dbn, 
                           inter_helix_gap=inter_helix_gap, 
                           overarching=overarching)

# dumping to file
minimal_expansion.dump_to_gr("results/gr_files/K5.gr")

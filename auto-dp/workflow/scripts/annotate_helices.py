from autodp.minimal_expansion import MinimalExpansion

inter_helix_gap=snakemake.config["inter_helix_gap"]
underlying_graph = MinimalExpansion()
dbn = open(snakemake.input.dbn).readlines()[0].rstrip('\n')
underlying_graph.from_str(dbn, inter_helix_gap=inter_helix_gap)

print(underlying_graph.helices)

f = open(snakemake.output[0], 'w')
for k, h in enumerate(underlying_graph.helices):
    print("H"+str(k),h, file=f)


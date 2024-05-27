from autodp.minimal_expansion import MinimalExpansion

inter_helix_gap=True
underlying_graph = MinimalExpansion()
dbn = open('resources/dbn_files/K5.dbn').readlines()[0].rstrip('\n')
underlying_graph.from_str(dbn, inter_helix_gap=inter_helix_gap)

print(underlying_graph.helices)

f = open("results/helix_annotations/K5.helix", 'w')
for k, h in enumerate(underlying_graph.helices):
    print("H"+str(k),h, file=f)


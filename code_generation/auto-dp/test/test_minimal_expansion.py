from minimal_expansion import MinimalExpansion

def test_H_type():

    min_exp_H = MinimalExpansion()

    min_exp_H.from_str('([)]')
    helices_should_be = set([(0,4,10,14),(5,9,15,19)])

    assert(helices_should_be==min_exp_H.helices)
    assert(0 in min_exp_H.adj[14])
    assert((0,1) in min_exp_H.edges)

    min_exp2 = MinimalExpansion()
    min_exp2.from_str('([)]', inter_helix_gap=False)

    helices_should_be = set([(0,4,8,12),(4,8,12,16)])

    assert(helices_should_be==min_exp2.helices)

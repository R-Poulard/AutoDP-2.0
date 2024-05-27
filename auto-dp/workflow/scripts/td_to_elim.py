bag_adj = {}
bag_content = {}

def add_to_adj(adj, key, vertex):
    try:
        adj[key].append(vertex)
    except KeyError:
        adj[key] = [vertex]

# readling file
after_b = False
for line in open(snakemake.input[0]).readlines():
    if line[0]=='b':
        after_b = True
        label = line.split(' ')[1]
        bag_content[label] = [v.rstrip('\n') for v in line.split(' ')[2:]]
    else:
        if after_b:
            u = line.split(' ')[0]
            v = line.split(' ')[1].rstrip('\n')
            add_to_adj(bag_adj, u, v)
            add_to_adj(bag_adj, v, u)

# to elim ordering:
root = '1' # arbitrary root

queue = [root]
bag_content['-1'] = []

# recursive function: vertices in bag are after orders for children
def build_order(parent, bag, bag_adj, bag_content):
    
    order = [u for u in bag_content[bag] if u not in bag_content[parent]]

    for child in bag_adj[bag]:
        if child!=parent:
            order = order + build_order(bag, child, bag_adj, bag_content)

    return order

def print_order_tree(parent, bag, bag_adj, bag_content, depth, f, positions, is_last, branching):

    print(positions, is_last, branching, end=" ")

    for k in range(depth):
        if k in positions:
            print('|', end="", file=f)
        else:
            print(' ', end="", file=f)
    
    if branching and is_last:
        positions.pop()

    print([u for u in bag_content[bag] if u not in bag_content[parent]],file=f)


    if len([c for c in bag_adj[bag] if c!=parent]) > 1:
        positions.append(depth)
        branching = True
    else:
        branching = False
    
    print([u for u in bag_content[bag] if u not in bag_content[parent]])

    for child in bag_adj[bag]:
        if child!=parent:
            if child!=[b for b in bag_adj[bag] if b!=parent][-1]:
                print_order_tree(bag, child, bag_adj, bag_content, depth+1, f, positions, False, branching)
            else:
                print_order_tree(bag, child, bag_adj, bag_content, depth+1, f, positions, True, branching)

f = open(snakemake.output[0], 'w')

order = build_order('-1', root, bag_adj, bag_content)

for vertex in order:
    if vertex==order[-1]:
        print(vertex, end="", file=f)
        continue
    print(vertex+'-', end="", file=f)
print("\n",file=f)
print_order_tree('-1', root, bag_adj, bag_content, 0, f, [], False, False)

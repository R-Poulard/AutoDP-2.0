import sys
from optparse import OptionParser
import os
import json

def convert_to_td1dot(pace_td_fname, input_dir, dataset_name):

    bag_content = {}

    for line in open(input_dir+pace_td_fname).readlines():
        if line[:2]=='b ':
            bag_number = line.split(' ')[1]
            content = line.split(' ')[2:]

            bag_content[bag_number] = content
    
    b_began = False

    adj = {}

    for line in open(input_dir+pace_td_fname).readlines():
        if line[:2]=='b ':
            b_began = True
        else:
            if b_began:
                bn1 = line.split(' ')[0]
                bn2 = line.split(' ')[1].rstrip('\n')

                try:
                    adj[bn1].append(bn2)
                except KeyError:
                    adj[bn1] = [bn2]

                try:
                    adj[bn2].append(bn1)
                except KeyError:
                    adj[bn2] = [bn1]


    root = min(adj.keys())
    
    visited = {k:False for k in adj.keys()}

    queue = [root]

    bag_new_number = {}
    counter = 0

    while len(queue) > 0:
        bn1 = queue.pop()
        visited[bn1] = True
        bag_new_number[bn1] = counter
        counter += 1
        for bn2 in adj[bn1]:
            if not visited[bn2]:
                queue.append(bn2)


    
    # writing into file
    outfile = open(snakemake.output[0],'w')

    outfile.write('graph G {\n')
    outfile.write('\n')
            

    for bag_number, content in sorted(bag_content.items(), key=lambda x:bag_new_number[x[0]]):
        outfile.write('\tbag'+str(bag_new_number[bag_number])+' [label="')

        for index in content:
            outfile.write(index+" ")

        outfile.write('"]\n')

    outfile.write("\n")

    for bn1 in sorted(bag_new_number.keys(), key=lambda x : bag_new_number[x]):
        for bn2 in sorted(adj[bn1], key=lambda x: bag_new_number[x]):
            if bag_new_number[bn2] > bag_new_number[bn1]:
                outfile.write('\t')
                outfile.write('bag'+str(bag_new_number[bn1]))
                outfile.write(' -- ')
                outfile.write('bag'+str(bag_new_number[bn2])+'\n')

    outfile.write('\n')
    outfile.write('}')

if __name__ == "__main__":
    
    pace_td_fname = snakemake.input[0].split('/')[-1]
    input_dir = snakemake.input[0].split(pace_td_fname)[0]
    dataset_name = snakemake.input[0].split('/')[2]

    convert_to_td1dot(pace_td_fname, input_dir, dataset_name)

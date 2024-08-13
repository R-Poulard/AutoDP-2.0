import re,sys
#check wether a folding got the right shaped
PSEUDO = sys.argv[1]  
DIRECTORY = sys.argv[2]

#H generation
if PSEUDO=="H":
    anchors = [[(0, 2)], [(1, 3)]]

if PSEUDO=="K":
    anchors = [[(0, 2)],[(1,4)],[(3,5)]]

if PSEUDO=="L":
    anchors = [[(0, 3)],[(1,4)],[(2,5)]]

if PSEUDO=="M":
    anchors = [[(0, 3)],[(1,5)],[(2,6)],[(4,7)]]

if PSEUDO=="C5":
    anchors = [[(0, 3)],[(1,8)],[(2,5)],[(4,7)],[(6,9)]]

if PSEUDO=="K4":
    anchors = [[(0, 4)],[(1,5)],[(2,6)],[(3,7)]]

if PSEUDO=="K5":
    anchors = [[(0, 5)],[(1,6)],[(2,7)],[(3,8)],[(4,9)]]

def get_fatgraph(folding):
    folding_list=[]
    for letter in "ABCDEFGH":
        #print("do ",letter)
        if letter not in folding:
            #print("not in it")
            if chr(ord(letter)+1) in folding:
                #print("next")
                return -1
            else:
                break
        else:
            big=[m.start() for m in re.finditer(letter,folding)]
            #print("big",big)
            little=[m.start() for m in re.finditer(letter.lower(),folding)]
            little.reverse()
            #print("little",little)
            if len(big)!=len(little):
                return -1
            else:
                folding_list.append(list(zip(big,little)))
    return folding_list

def is_equal(ref,fatgraph):
    try:
        #print(anchors,fatgraph)
        flat_list = [item for sublist in fatgraph for item in sorted(sublist[:1])]
    except Exception:
        print(anchors,fatgraph)
    
    #print(flat_list)
    fl = [item for sublist in flat_list for item in sublist]
    corresponding={b:a for a,b in  enumerate(sorted(fl))}
    #print(corresponding)
    reduced=[[(corresponding[a[0]],corresponding[a[1]])] for a in sorted(flat_list)]
    #print(reduced)
    if ref!=reduced:
        print(ref,fatgraph,reduced)
    return ref==reduced

def check_result(path):
    failed_lines=[]
    with open(path,"r") as f:
        f.readline()
        lines=f.readlines()
        for i in range(2,len(lines),2):
            
            if "End" in lines[i]:
                break
            #print(lines[i+1])
            if False and "Correct" not in lines[i]:
                failed_lines.append((lines[i],lines[i+1]))
            else:
                ft=get_fatgraph(lines[i+1])
                if ft==-1:
                    print(lines[i],lines[i+1])
                    failed_lines.append((lines[i],lines[i+1]))
                elif not is_equal(anchors,ft):
                    print(lines[i+1])
                    failed_lines.append((lines[i],lines[i+1]))
    return failed_lines

flines=check_result(DIRECTORY)
if len(flines)==0:
    print(" All lines are correct")
else:
    print(" Error on ",len(flines),"lines")
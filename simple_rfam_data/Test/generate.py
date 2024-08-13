import random
import os 

#os.chdir("..")

import math

lol=0
cant=100
nb=0
"""
for num in range(15977609,50000,-1):
    for i in range(2,num-1):
        pls=num%i
        
        print(num,i,pls)
        if(pls!=0):
            break
    if all(num%i!=0 for i in range(2,int(math.sqrt(num-1)))):
        if lol==0:
            print (num,",",end="",flush=True)
            nb+=1
            lol=20000
            cant-=1
            if cant==0 :
                break
        else:
            lol-=1
print("-----",nb)   """    

primer=[15977609 ,15646489 ,15315259 ,14984819 ,14654509 ,14324903 ,13995577 ,13666729 ,13337969 ,13009019 ,12681707 ,12355493 ,12029063 ,11703061 ,11379097 ,11054227 ,10729597 ,10407779 ,10084637 ,9761327 ,9440309 ,9119749 ,8798171 ,8480621 ,8161763 ,7843463 ,7526647 ,7209667 ,6894541 ,6581137 ,6267619 ,5954411 ,5643529 ,5332127 ,5023141 ,4716001 ,4410389 ,4103903 ,3800617 ,3498493 ,3197693 ,2899417 ,2603033 ,2308483 ,2016697 ,1728691 ,1443059 ,1161683 ,884951 ,614683 ,353333 ,108127 ,15977609 ,15646489 ,15315259 ,14984819 ,14654509 ,14324903 ,13995577 ,13666729 ,13337969 ,13009019 ,12681707 ,12355493 ,12029063 ,11703061 ,11379097 ,11054227 ,10729597 ,10407779 ,10084637 ,9761327 ,9440309 ,9119749 ,8798171 ,8480621 ,8161763 ,7843463 ,7526647 ,7209667 ,6894541 ,6581137 ,6267619 ,5954411 ,5643529 ,5332127 ,5023141 ,4716001 ,4410389 ,4103903 ,3800617 ,3498493 ,3197693 ,2899417 ,2603033 ,2308483 ,2016697 ,1728691 ,1443059 ,1161683 ,884951 ,614683 ,353333 ,108127 ]
for num in primer:
    if not all(num%i!=0 for i in range(2,num-1)):
        print(num,"false")
exit()
PSEUDO="H"
#H generation
if PSEUDO=="H":
    seq = "AGUC"
    anchors = [[(0, 2)], [(1, 3)]]
    score = -15

if PSEUDO=="K":
    #K generation
    seq = "ACUAGU"
    anchors = [[(0, 2)],[(1,4)],[(3,5)]]
    score = -20

if PSEUDO=="L":

    #L generation
    seq = "AUGUAC"
    anchors = [[(0, 3)],[(1,4)],[(2,5)]]
    score = -20

if PSEUDO=="M":
    #M generation
    seq = "GACCGUGC"
    anchors = [[(0, 3)],[(1,5)],[(2,6)],[(4,7)]]
    score = -35

if PSEUDO=="C5":
    #C5 generation
    seq = "AGCUAGGUCC"
    anchors = [[(0, 3)],[(1,8)],[(2,5)],[(4,7)],[(6,9)]]
    score = -40

if PSEUDO=="K4":
    #K4 generation
    seq = "AGCAUCGU"
    anchors = [[(0, 4)],[(1,5)],[(2,6)],[(3,7)]]
    score = -30

if PSEUDO=="K5":
    #K5 generation
    seq = "AGCAGUCGUC"
    anchors = [[(0, 5)],[(1,6)],[(2,7)],[(3,8)],[(4,9)]]
    score = -40

ss='.'*len(seq)
length = len(seq)
target_length = 50

ref=anchors
num_anchors = len(anchors)

def bigger(pair, where):
    a, b = pair
    if a >= where:
        a += 1
    if b >= where:
        b += 1
    return (a, b)

def update_anchors(anchors, where):
    return [[bigger(pair, where) for pair in anchor] for anchor in anchors]


def update_anchors_pairs(anchors,where):
    x,y=where
    #print("values = ",x,y)
    updated_anchors = []
    for anchor in anchors:
        updated_pairs = []
        for pair in anchor:
            a, b = pair
            #print(a,b," becomes",end=" ")
            if a >= x:
                a += 1
            if a>=y:
                a+=1
   
            if b >= x:
                b += 1
            if b>=y:
                b+=1
            #print(a,b)
            updated_pairs.append((a, b))
        updated_anchors.append(updated_pairs)
    return updated_anchors

def generate(anchors,seq,ss,target_length,score):
    while len(seq) < target_length:
        #print("Sequence:", seq)
        #print("Anchors:", anchors)
        choice = random.randint(0, 3)
        
        if choice <=2:  # Add 'X' (unlink) anywhere
            where = random.randint(0, len(seq))
            which = random.randint(0,3)
            print("ici",which)
            nuc=""
            if which==0:
                print("ici")
                nuc="A"
            elif which==1:
                nuc="G"
            elif which==2:
                nuc="C"
            elif which==3:
                nuc="U"
            print(nuc)
            seq = seq[:where] + nuc + seq[where:]
            ss = ss[:where] + "x" + ss[where:]
            anchors = update_anchors(anchors, where)
            score=score-1
        else:  # Add 2 characters (link) in an anchor
            helix = random.randint(0, num_anchors - 1)
            before = random.randint(0, 1)
            
            if before == 1:  # Exterior add
                #print("exterior")
                x_prev, y_prev = anchors[helix][0]
                x_new = x_prev
                y_new = y_prev + 2
                anchors = update_anchors_pairs(anchors, (x_new,y_new))
                anchors[helix].insert(0, (x_new, y_new))
                #print("insertion of ",x_new,y_new)
                which = random.randint(0, 3)
                #print(seq[:x_new],seq[x_new:y_new],seq[y_new:])
                if which == 0:
                    seq = seq[:x_new] + "A" + seq[x_new:y_new-1] + "U" + seq[y_new-1:]
                    ss = ss[:x_new] + "." + ss[x_new:y_new-1] + "." + ss[y_new-1:]
                    score -= 5
                elif which == 1:
                    seq = seq[:x_new] + "U" + seq[x_new:y_new-1] + "A" + seq[y_new-1:]
                    ss = ss[:x_new] + "." + ss[x_new:y_new-1] + "." + ss[y_new-1:]
                    score -= 5
                elif which == 2:
                    seq = seq[:x_new] + "G" + seq[x_new:y_new-1] + "C" + seq[y_new-1:]
                    ss = ss[:x_new] + "." + ss[x_new:y_new-1] + "." + ss[y_new-1:]
                    score -= 10
                elif which == 3:
                    seq = seq[:x_new] + "C" + seq[x_new:y_new-1] + "G" + seq[y_new-1:]
                    ss = ss[:x_new] + "." + ss[x_new:y_new-1] + "." + ss[y_new-1:]
                    score -= 10

            else:  # Interior add
                #print("interior")
                x_prev, y_prev = anchors[helix][-1]
                x_new = x_prev + 1
                y_new = y_prev +1
                print("avant1",anchors)
                anchors = update_anchors_pairs(anchors, (x_new,y_new))
                print("avant",anchors)
                anchors[helix].append((x_new, y_new))
                print("apres",anchors)
                #print("insertion of ",x_new,y_new)
                which = random.randint(0, 3)
                #print(seq[:x_new],seq[x_new:y_new-1],seq[y_new-1:])
                if which == 0:
                    seq = seq[:x_new] + "A" + seq[x_new:y_new-1] + "U" + seq[y_new-1:]
                    ss = ss[:x_new] + "." + ss[x_new:y_new-1] + "." + ss[y_new-1:]
                    score -= 5
                elif which == 1:
                    seq = seq[:x_new] + "U" + seq[x_new:y_new-1] + "A" + seq[y_new-1:]
                    ss = ss[:x_new] + "." + ss[x_new:y_new-1] + "." + ss[y_new-1:]
                    score -= 5
                elif which == 2:
                    seq = seq[:x_new] + "G" + seq[x_new:y_new-1] + "C" + seq[y_new-1:]
                    ss = ss[:x_new] + "." + ss[x_new:y_new-1] + "." + ss[y_new-1:]
                    score -= 10
                elif which == 3:
                    seq = seq[:x_new] + "C" + seq[x_new:y_new-1] + "G" + seq[y_new-1:]
                    ss = ss[:x_new] + "." + ss[x_new:y_new-1] + "." + ss[y_new-1:]
                    score -= 10
    return seq,ss,anchors,score

def intricated(x1,x2,y1,y2,pairing):
    for m in pairing:
        a,b=m
        if a<x2 and a>x1 and not( b>y2 and b<y1):
            return False
        if  b>y2 and b<y1 and not (a<x2 and a>x1):
            return False
    return True
        
def get_fatgraph(pairing):
    pairing=sorted(pairing)
    to_remove=[]
    for starting in pairing:
        x1,y1=starting
        if starting in to_remove:
            continue
        for y in pairing[pairing.index(starting)+1:]:
            x2,y2=y
            if x2>x1 and x2<y1 and y2>x1 and y2<y1 and intricated(x1,x2,y1,y2,pairing):
                to_remove.append(y)
            else:
                break
    return  [sublist for sublist in pairing if sublist not in to_remove]

def is_equal(ref,fatgraph):
    flat_list = [item for sublist in fatgraph for item in sublist]
    corresponding={b:a for a,b in  enumerate(sorted(flat_list))}
    reduced=[[(corresponding[a],corresponding[b])] for a,b in fatgraph]
    return ref==reduced



def create_file(path,length,nb):
    i=0
    with open(path+"_folding.txt","w") as f2:
        with open(path+".txt","w") as f:
            while(i<nb):
                print("i=",i)
                i+=1
                new_seq,new_ss,new_anchors,new_score=generate(anchors,seq,ss,length,score)
                print("new anchors",new_anchors)
                ft=get_fatgraph([item for sublist in new_anchors for item in sublist])
                print(len(new_seq))
                if len(ft)!=num_anchors or not is_equal(ref,ft):
                    print(new_seq)
                    print(new_anchors)
                    print(new_score)
                    print(ft)
                    exit()
                f.write(new_seq+"\n")
                f.write(new_ss+"\n")
                f.write(str(new_score)+"\n")
                print(new_anchors,file=f2,end="\n")

def check_result(path):
    failed_lines=[]
    with open(path,"r") as f:
        f.readline()
        lines=f.readlines()
        for i in range(0,len(lines),2):
            print(i)
            if "End" in lines[i]:
                break
            if "Correct" not in lines[i]:
                failed_lines.append((lines[i],lines[i+1]))
            else:
                folding=[]
                for i in lines[i+1].strip().split(","):
                    if i=="":
                        break
                    else:
                        m=i.split("->")
                        a,b=int(m[0]),int(m[1])
                        if(a==b):
                            failed_lines.append((lines[i],lines[i+1]))
                            continue
                        if a<0 or b<0:
                            failed_lines.append((lines[i],lines[i+1]))
                            continue
                        folding.append((min(a,b),max(a,b)))
                if len([item for sublist in folding for item in sublist])!=len(set([item for sublist in folding for item in sublist])):
                    failed_lines.append((lines[i],lines[i+1]))
                    continue
                if not is_equal(ref,get_fatgraph(folding)):
                    failed_lines.append((lines[i],lines[i+1]))
    return failed_lines


create_file("test_sequence/example_"+PSEUDO+"_big",75,50)
#print(check_result("logK4.txt"))
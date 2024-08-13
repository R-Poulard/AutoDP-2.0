#ce programme récupère les séquences ayant des pseudonoeuds de types H en un seul fichier

with open("../pseudoknotted_families.txt","r") as f:
    file_list=f.read().split('\n')

print(len(file_list),file_list)
retained_file=[]
seq=[]
for i in file_list:
    try:
        with open(i+".stk","r") as f:
            match=False
            text=f.readlines()
            for line in text:
                if "#=GR" in line:
                    line=line.replace("#=GR","").replace("\n","").replace(" ","")
                    if "A" in line and "B" in line and "C" not in line:
                        match=True
                        
                        retained_file.append(i)
                        break
                    else:
                        break
            if match is True:
                for line in text:
                    if '#' not in line:
                        seq.append([x for x in line.split(" ") if x!=""][1])
                        
    except:
        continue

with open("../sequence_H.txt","w") as f:
    f.write('\n'.join(seq))





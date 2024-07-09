import requests
from bs4 import BeautifulSoup

response=requests.get("https://rnavlab.utep.edu/dbresults?a_sequence=&a_keyword=&a_pkbno=&a_organism=&a_comment=&a_submittedby=&span_range=greater_than&a_span=&a_rnatype=&a_supportedby=&a_pseudoknot=H&loop_range=equal_to&a_loop=&stem_range=equal_to&a_stem=&a_continuous=&display=pkbn+ncbi_number+abbreviation+organism+rna_type+supported_by+classification+definition",verify=False)
soup=BeautifulSoup(response.text,'html.parser')

links=[]
a_tags = soup.find_all('a')
for tag in a_tags:
    if tag.get_text()!="view":
        continue
    links.append(tag.get('href'))

sequences=[]
dot_bracket=[]

i=0

for lk in links:
    response=requests.get("https://rnavlab.utep.edu"+lk,verify=False)
    i+=1
    print(response.text)
    print(lk,i,"/",len(links),end="    ")
    if response.status_code!=200:
        print("Unable to read "+lk+" (Code:"+response.status_code+")")
        continue
    r=response.text
    soup=BeautifulSoup(response.text,'html.parser')
    interest=soup.find('pre')
    
    seq=""
    db=""
    for dt in interest.get_text().split("\n"):
        dt=dt.strip()
        #print("line = ", dt)

        if(len(dt)<=1):
            continue
        if dt[0]=="$":
            #print("here")
            seq+=dt.split(" ")[2].split("=")[0]
        if dt[0]=="%":
            #print("there")
            db+=dt.split(" ")[2].replace(":",".").replace("[","B").replace("]","b").replace("(","A").replace(")","a")
    print("result=",seq,db)
    sequences.append(seq)
    dot_bracket.append(db)
    print("done")

pnb=[">"+lk.replace("/static/PKB_files/","") for lk in links]
with open("./test_sequence/sequence_H.txt","w") as f:
    for line in zip(pnb,sequences,dot_bracket):
        f.write("#"+line[0]+"\n")
        f.write(line[1]+"\n")
        f.write("."*len(line[1])+"\n")
        f.write("0\n")
        f.write("#"+line[2]+"\n")

with open("./res_benchmark_stack_biggnumbers.txt","r") as f:
    new_lines=""
    full_lines=[]
    full_lines.append("nb_elem;capacity;initial_capacity;biggest_chain;prct_chain;mean_chain_size;init_cap;time_exec;final_cap;disk_usage\n")
    for line in f.readlines():
      if("#-#" not in line and "/" in line):
         continue
      if "#-#" in line:
         if "#-#--------------------------" in line:
            continue
         if "#-#--TABLE STAT REPORT-----" in line:
            continue
         if "#-#nb of elem: "in line:
            line=line.replace("#-#nb of elem: ","")
         if "#-#capcaity: " in line:
            line=line.replace("#-#capcaity:","")
         if "#-#initial capacity " in line:
            line=line.replace("#-#initial capacity ","")
         if "#-#Biggest chain: " in line:
            line=line.replace("#-#Biggest chain: ","")
         if "#-#Biggest chain: " in line:
            line=line.replace("#-#Biggest chain: ","")
         if "#-# prct chaine: " in line:
            line=line.replace("#-# prct chaine: ","")
            line=line.split(" ")[0]
         if "#-# mean chaine size: " in line:
            line=line.replace("#-# mean chaine size: ","")

         new_lines+=line.strip()+";"
      elif ";" in line:
         new_lines+=line
         full_lines.append(new_lines)  
         print(new_lines)
         new_lines=""   
    with open("bch_normal_big_stack.csv","w") as f2:
       for i in full_lines:
          f2.write(i)
    
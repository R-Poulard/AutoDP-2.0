import_statements = ""
linkage_statements = ""
for shd in snakemake.input:
    shd=shd.replace("_folding.c","").replace("results/c_code/","")
    import_statements += f'\n#include "{shd}_folding.h"'
    linkage_statements += f'\n    setup_{shd}(fc,fc_aux, seq, con, bt);'
print(import_statements)
print(linkage_statements)
with open(snakemake.output[0], "w") as f:
    with open("resources/stack_Turner/auto-dp.c", "r") as template:
        template_content = template.read()
        updated_content = template_content.replace("##ADDS_ON_INCLUDE##", import_statements)
        updated_content = updated_content.replace("##LINKAGE##", linkage_statements)
        print(updated_content,file=f)

if True or not os.path.exists(DIRECTORY+"aux.c"):
    f = open("results/c_code/aux.c",'w')
    with open('resources/stack_Turner/aux.c','r') as f2:
        print("".join(f2.readlines()),file=f)
if True or not os.path.exists(DIRECTORY+"aux.h"):
    f = open("results/c_code/aux.h",'w')
    with open('resources/stack_Turner/aux.h','r') as f2:
        print("".join(f2.readlines()),file=f)
exit()
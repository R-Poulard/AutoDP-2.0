cd ~/pwb/Turner/
mkdir $1
mkdir $1/log_test/
mkdir $1/log_build/
list="H K K4 K5 L M C5"
for i in $list; do 
    echo "###########Generating $i############";
    /home/remipoul/miniconda3/bin/python /home/remipoul/pwb/auto-dp/workflow/scripts/produce_c_code2_Turner.py $i $1 > $1/log_build/$i.txt;
    echo "Testing $i:";
    gcc /home/remipoul/pwb/Turner/$1/$i.c -o /home/remipoul/pwb/Turner/$1/$i -lm -lRNA -g 2>> $1/log_build/$i.txt &&  "/home/remipoul/pwb/Turner/$1/"$i "/home/remipoul/pwb/test_sequence/example_${i}_little.txt" > $1/log_test/$i.txt;
    echo "  Result $? if 0 not failed";
    echo "Check Result validity:";
    /home/remipoul/miniconda3/bin/python /home/remipoul/pwb/check_fatgraph.py $i /home/remipoul/pwb/Turner/$1/log_test/$i.txt
done
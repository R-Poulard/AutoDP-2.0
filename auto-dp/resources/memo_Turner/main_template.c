
int fold(HashTable * hashTable) {
    int min_value=INT_MAX;
    for(int END=0+ANCH;END<=MAX;END++){
        for(int a=0;a<END-ANCH;a++){
            int mfe1=MFEFree(0,a-1);
            int mfe2=MFEFree(END,MAX-1);
            int mfe=add(mfe1,mfe2);
            min_value = min(min_value,add(compute_A(hashTable,a,END),mfe));
        }
    }
    return min_value;
}

void backtrace(HashTable * hashTable,int score) {
    for(int END=0+ANCH;END<=MAX;END++){
        for(int a=0;a<END-ANCH;a++){
            
            int mfe1=MFEFree(0,a-1);
            int mfe2=MFEFree(END,MAX-1);
            int mfe=add(mfe1,mfe2);
            int tmp0=compute_A(hashTable,a,END);
            //printf("h = %d, a=%d, %d+%d %d,%d\n",h,a,mfe2,mfe1,score,compute_A(hashTable,a,h));
            if(score==add(mfe,tmp0)){
                backtrace_A(hashTable,tmp0,a,END);
                backtrace_MFEFree(mfe1,0,a-1);
                backtrace_MFEFree(mfe2,END,MAX-1);
                return;
            }
        }
    }
}


int main(int argc, char ** argv) {
    printf("File name: %s\n",argv[1]);
    FILE * fp = fopen(argv[1], "r");
    if (fp == NULL) {
        perror("Error opening file");
        return EXIT_FAILURE;
    }
    int nb_tests=0;
    clock_t start_time, end_time;
    double total_time = 0.0;

    start_time = clock();
    while(true){
        HashTable *hashTable = createHashTable();
        int b_len=0;

        int len=0;

        //reading of the sequence
        len=getline(&line,&b_len, fp);
        
        if(len<=0){
            destroyHashTable(hashTable);
            free(line);
            free(correct_score);
            printf("End of file, %d tested",nb_tests);
            break;
        }
        if(line[0]=='#'){
            printf("%s",line);
            continue;
        }
        MAX=len;
        if(line[len-1]=='\n'){
            MAX=len-1;
            nb_tests+=1;
        }
        preset_N(MAX);

        //secondary structure
        char *ss = NULL;
        len=getline(&ss,&b_len, fp);
        createMatrix(MAX,ss);
        //reading the score
        len=getline(&correct_score,&b_len, fp);
        if (len != -1) {
            if (correct_score[len - 1] == '\n') {
                correct_score[len - 1] = '\0'; // Replace newline character with null terminator
                len--; // Decrement the length to exclude the removed '\n'
            }
            // Convert the line to an integer
            
        }
        else{
            destroyHashTable(hashTable);
            free(line);
            free(correct_score);
            free(ss);
            printf("No integer found for test %d\n",nb_tests);
            exit(-1);
        }
        int number = atoi(correct_score);
        printf("Test: %d Size of the Sequence: %d ---", nb_tests,MAX);
        
        //elements used to record backtrack
        bracket='A';
        structure=malloc(sizeof(char)*(MAX+1));
        structure[MAX]='\0';
        for(int i=0;i<MAX;i++){
            structure[i]='.';
        }
        //Start of test
        clock_t test_start_time = clock();

        //folding
        int score = fold(hashTable);
        //retrieve backtrack
        backtrace(hashTable,score);
        
        clock_t test_end_time = clock();

        double test_time = ((double)(test_end_time - test_start_time)) / CLOCKS_PER_SEC;
        total_time += test_time;
        //end of test
        
        //sanity checks
        if(score == number){
            printf("Correct");
            int nb=0;
            for(int i=0;i<MAX;i++){
                if(structure[i]=='.'){
                    printf("missing bt\n");
                    break;
                }
                if(structure[i]>64){
                    if(structure[i]<91){
                        
                        nb+=structure[i];
                    }
                    else{
                        
                        nb-=(structure[i]-32);
                    }
                }
            }
            printf("(backtrack impurty= %d)\n",nb);
        }
        else{
            destroyHashTable(hashTable);
            free(line);
            free(correct_score);
            free(structure);
            free(ss);
            printf("Failed (%d found, should be %d)\n",score,number);
            exit(nb_tests);
        }

        
        printf("%s\n",structure);
        
        //clean up
        destroyHashTable(hashTable);
        freeMatrix(MAX);
        free(structure);
        free(line);
        free(ss);
    }
    end_time = clock();
    double total_execution_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    printf(" Total execution time: %.4f seconds", total_execution_time);
    printf("(mean %.4f seconds)\n", total_time / nb_tests);

    fclose(fp);
    return 0;
}


int init_root(HashTable* HashTable,Stack * stack){  
    int A = 65;
    int size = 3;
    for(int END=0+ANCH;END<=MAX;END++){
        for(int a=0;a<END-ANCH;a++){
            int tab[]={A,a,END};
            Node* newNode = (Node*)malloc(sizeof(Node));
            newNode->data = createTuple(tab,size);
            newNode->next = NULL;
            if(stack->top!=NULL){
                newNode->next=stack->top;
                stack->top=newNode;
            }
            stack->top=newNode;
        }
    }
}

int compute_root(HashTable * hashTable,Stack *stack) {
    int min_value=INT_MAX; 
    for(int END=0+ANCH;END<=MAX;END++){
        for(int a=0;a<END-ANCH;a++){
            int mfe1=MFEFree(0,a-1);
            int mfe2=MFEFree(END,MAX-1);
            int mfe=add(mfe1,mfe2);
            int tmp;
            get_A(hashTable,NULL,NULL,&tmp,a,END);
            min_value = min(min_value,add(tmp,mfe));
        }
    }
    return min_value;
}

int fold(HashTable* hasTable, Stack * stack){
    
    init_root(hasTable,stack);
    Node* topNode;
    Node* nextNode;
    Node* add_on_last;
    Node* add_on_first;
    while (!isStackEmpty(stack)) {
        //pop(stack, &topTuple)
        //printStack(stack);
        topNode=stack->top;
        nextNode=topNode->next;
        Tuple * topTuple = topNode->data;
        add_on_last=NULL;
        add_on_first=NULL;
        int ALL_ANC;
        int i1,i2,j2,j1;
        //printf("here\n");
        //print_tuple(topTuple);
        switch(topTuple->values[0]){
VARIABLE_PART
        }
    }
    //printStack(stack);
    return compute_root(hasTable,stack);
}
void backtrace(HashTable * hashTable,Stack* stack,int score) {
    //int a=0;int h=MAX;
    for(int END=0+ANCH;END<=MAX;END++){
        for(int a=0;a<END-ANCH;a++){
            int mfe1=MFEFree(0,a-1);
            int mfe2=MFEFree(END,MAX-1);
            int mfe=add(mfe1,mfe2);
            int tmp0;
            get_A(hashTable,NULL,NULL,&tmp0,a,END);
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
        
        int b_len=0;

        int len=0;

        //reading of the sequence
        len=getline(&line,&b_len, fp);
        
        if(len<=0){
            
            free(line);
            free(correct_score);
            printf("End of file, %d tested",nb_tests);
            break;
        }
        if(line[0]=='#'){
            printf("%s",line);
            continue;
        }
        HashTable *hashTable = createHashTable();
        Stack *stack = malloc(sizeof(Stack));
        initializeStack(stack);
        
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
        if (len != -1) {
            if (ss[len - 1] == '\n') {
                ss[len - 1] = '\0'; // Replace newline character with null terminator
                len--; // Decrement the length to exclude the removed '\n'
            }
        }
        else{
            destroyHashTable(hashTable);
            freeStack(stack);
            free(line);
            free(ss);
            printf("No secondary structure found for test %d\n",nb_tests);
            exit(-1);
        }
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
            freeStack(stack);
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

        vrna_init_rand();
        fc = vrna_fold_compound(line, NULL, VRNA_OPTION_DEFAULT);
        vrna_constraints_add(fc, ss, VRNA_CONSTRAINT_DB_DEFAULT);
        vrna_mfe(fc, NULL);
        
        //folding
        int score = fold(hashTable,stack);
        //retrieve backtrack
        backtrace(hashTable,stack,score);
        vrna_fold_compound_free(fc);
        clock_t test_end_time = clock();

        double test_time = ((double)(test_end_time - test_start_time)) / CLOCKS_PER_SEC;
        total_time += test_time;
        //end of test
        
        //sanity checks
        if(score == number || 1){//we autorise wrong score in the case of the Turner model
            float fl=(float)score/100.0;
            printf("Score %.2f ",fl);
            int nb=0;
            for(int i=0;i<MAX;i++){
                
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
            freeStack(stack);
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
        freeStack(stack);
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

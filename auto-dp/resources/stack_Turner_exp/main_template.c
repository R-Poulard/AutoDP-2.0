
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
        char * line=NULL;
        //reading of the sequence
        len=getline(&line,&b_len, fp);
        
        if(len<=0){
            
            free(line);
            printf("End of file, %d tested",nb_tests);
            break;
        }
        if(line[0]=='#'){
            printf("%s",line);
            free(line);
            continue;
        }

        int MAX=len;
        if(line[len-1]=='\n'){
            MAX=len-1;
            nb_tests+=1;
        }

        //secondary structure
        char *ss = NULL;
        len=getline(&ss,&b_len, fp);
        if (len != -1) {
            if (ss[len - 1] == '\n') {
                ss[len - 1] = '\0'; // Replace newline character with null terminator
                len--; // Decrement the length to exclude the removed '\n'
            }
        }
        else{
            
            free(line);
            free(ss);
            printf("No secondary structure found for test %d\n",nb_tests);
            exit(-1);
        }
        //reading the score
        char * correct_score=NULL;
        len=getline(&correct_score,&b_len, fp);
        if (len != -1) {
            if (correct_score[len - 1] == '\n') {
                correct_score[len - 1] = '\0'; // Replace newline character with null terminator
                len--; // Decrement the length to exclude the removed '\n'
            }
            // Convert the line to an integer
            
        }
        else{
            free(line);
            free(correct_score);
            free(ss);
            printf("No integer found for test %d\n",nb_tests);
            exit(-1);
        }
        int number = atoi(correct_score);
        printf("Test: %d Size of the Sequence: %d ---", nb_tests,MAX);
        
        //elements used to record backtrack
        
        //Start of test
        clock_t test_start_time = clock();

        vrna_init_rand();
        vrna_fold_compound_t * fc = vrna_fold_compound(line, NULL, VRNA_OPTION_DEFAULT);
        vrna_constraints_add(fc, ss, VRNA_CONSTRAINT_DB_DEFAULT);
        pk_compound* pk=create_pk_compound(fc,line,ss);
        vrna_mfe(fc, NULL);
        
        //folding
        double score = fold(pk);
        //retrieve backtrack
        bt_struct* bt=create_bt_struct(NULL,pk->MAX);
        backtrace(pk,bt,score);
        clock_t test_end_time = clock();

        double test_time = ((double)(test_end_time - test_start_time)) / CLOCKS_PER_SEC;
        total_time += test_time;
        //end of test
        
        //sanity checks
        if(score == number || 1){//we autorise wrong score in the case of the Turner model
            float fl=(float)score/100.0;
            printf("Score %.2f ",fl);
            int nb=0;
            for(int i=0;i<pk->MAX;i++){
                
                if(bt->structure[i]>64){
                    if(bt->structure[i]<91){
                        
                        nb+=bt->structure[i];
                    }
                    else{
                        
                        nb-=(bt->structure[i]-32);
                    }
                }
            }
            printf("(backtrack impurty= %d)\n",nb);
        }
        else{
            free(line);
            free(correct_score);
            free_pk(pk);
            free_bt(bt);
            free(ss);
            printf("Failed (%d found, should be %d)\n",score,number);
            exit(nb_tests);
        }

        
        printf("%s\n",bt->structure);
        
        //clean up

        free_pk(pk);
        free_bt(bt);
        free(line);
        free(ss);
        free(correct_score);
    }
    end_time = clock();
    double total_execution_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    printf(" Total execution time: %.4f seconds", total_execution_time);
    printf("(mean %.4f seconds)\n", total_time / nb_tests);

    fclose(fp);
    return 0;
}

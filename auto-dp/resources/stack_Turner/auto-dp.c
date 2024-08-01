#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>
#include <ViennaRNA/mfe.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/structures.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/grammar.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/mfe.h"
#include "aux.h"
##ADDS_ON_INCLUDE##


int folding(const char* seq, const char* con){

    /* initialize random number generator */
    vrna_init_rand();
    /* Generate a random sequence of 50 nucleotides */
    vrna_fold_compound_t *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);
    vrna_constraints_add(fc, con, VRNA_CONSTRAINT_DB_DEFAULT);
    vrna_fold_compound_t *fc_aux = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);
    vrna_constraints_add(fc_aux, con, VRNA_CONSTRAINT_DB_DEFAULT);
    vrna_mfe(fc_aux,NULL);
    bt_struct * bt=create_bt_struct(NULL,strlen(seq));
    //displayMatrix(strlen(seq));
    /* allocate memory for MFE structure (length + 1) */
    char *structure = (char *)vrna_alloc(sizeof(char) * (strlen(seq) + 1));
    /* predict Minmum Free Energy and corresponding secondary structure */

##LINKAGE##

    vrna_gr_set_serialize_bp(fc,post_process_f,bt,NULL,NULL);

    float mfe = vrna_mfe(fc, structure);
    int nb=0;
    printf("Score [%6.2f] ",mfe);
    printf("%s\n", structure);

    /* cleanup memory */

    vrna_fold_compound_free(fc);
    vrna_fold_compound_free(fc_aux);
    free_bt(bt);
    free(structure);
    return mfe;
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
      char * line=NULL;
      //reading of the sequence
      len=getline(&line,&b_len, fp);
        
      if(len<=0){    
            free(line);
            printf("End of file, %d done",nb_tests);
            break;
      }
      if(line[0]=='#'){
            printf("%s",line);
            free(line);
            continue;
      }

      int MAX=len;
      if(line[len-1]=='\n'){
        line[len-1]='\0';
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
            printf("No secondary structure found for sequence %d\n",nb_tests);
            exit(-1);
      }
      
        
        printf("Sequence: %d Size of the Sequence: %d---", nb_tests,MAX);
        
        //elements used to record backtrack
        
        //Start of test
        clock_t test_start_time = clock();
        printf("%s\n",line);
        
        folding(line,ss);
        fflush(NULL);
        clock_t test_end_time = clock();

        double test_time = ((double)(test_end_time - test_start_time)) / CLOCKS_PER_SEC;
        total_time += test_time;
        //end of test
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
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

#include "H_folding.h"

int nb_tests=0;
int folding(const char* seq, const char* con){
    if(strlen(seq)>=100){
      return -1;
    }
    /* initialize random number generator */
    vrna_init_rand();
    /* Generate a random sequence of 50 nucleotides */
    vrna_fold_compound_t *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);
    vrna_constraints_add(fc, con, VRNA_CONSTRAINT_DB_DEFAULT);
    vrna_fold_compound_t *fc_aux = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);
    vrna_constraints_add(fc_aux, con, VRNA_CONSTRAINT_DB_DEFAULT);
    int mfe_normal=vrna_mfe(fc_aux,NULL);
    bt_struct * bt=create_bt_struct(NULL,strlen(seq));
    //displayMatrix(strlen(seq));
    /* allocate memory for MFE structure (length + 1) */
    char *structure = (char *)vrna_alloc(sizeof(char) * (strlen(seq) + 1));
    /* predict Minmum Free Energy and corresponding secondary structure */


    setup_H(fc,fc_aux, seq, con, bt);

    vrna_gr_set_serialize_bp(fc,post_process_f,bt,NULL,NULL);

    float mfe = vrna_mfe(fc, NULL);
    int nb=0;
    printf("%d,%6.2f;%6.2f\n",nb_tests,mfe,mfe_normal);
    //printf("%s\n", structure);

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
  //int nb_tests=0;
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
        line[len-1]='\0';
        MAX=len-1;
        nb_tests+=1;
      }

        //secondary structure
      char *ss = malloc(MAX);
      ss[MAX-1]='\0';
      for(int i=0;i<MAX;i++){
        ss[i]='.';
      }
        
        //printf("Sequence: %d Size of the Sequence: %d %d---", nb_tests,MAX,strlen(line));
        
        //elements used to record backtrack
        
        //Start of test
        clock_t test_start_time = clock();
        
        
        int lol=folding(line,ss);
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

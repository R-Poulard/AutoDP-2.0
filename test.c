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
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/mfe.h"
#include "H.c"

typedef struct {
    HashTable *ht;
    Stack *st;
    int** matrix;
    char *line;
}Custom_data;

int custom_ud_exp_prod_cb(vrna_fold_compound_t *fc,
                                         int i,
                                         int j,
                                         unsigned int loop_type,
                                         void *data) {
    //printf("%d %d\n",i,j);
  Custom_data *dt=(Custom_data*)data;
  int tmp;
  Node* first;
  Node* last;
  get_A(dt->ht,&first,&last,&tmp,i-1,j-1);
  printf("got %d %d,%d\n",tmp,i,j);
  return tmp;
}
void test(vrna_fold_compound_t *fc,void *data) {
  printf("ici\n");
  Custom_data *dt=(Custom_data*)data;
  printf("ici\n");
  dt->ht;
  printf("%d\n",dt->ht->size);
  dt->st;
  fold(dt->ht,dt->st);
  printf("finnnnnnnnn\n");
  //print_table(dt->ht);
  return;
}
void free_data(void *data){
  Custom_data *dt=(Custom_data*)data;
  freeStack(dt->st);
  destroyHashTable(dt->ht);
  freeMatrix(strlen(dt->line));
  free(dt->line);
}
int
main()
{
/* initialize random number generator */
vrna_init_rand();
/* Generate a random sequence of 50 nucleotides */
const char *seq = "CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUC";
const char *con = "......................................"; // Hard constraint dot-bracket notation

vrna_fold_compound_t *fc2 = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);
vrna_constraints_add(fc2, con, VRNA_CONSTRAINT_DB_DEFAULT);

Custom_data *data=malloc(sizeof(Custom_data));
data->ht = createHashTable();
data->st = malloc(sizeof(Stack));
initializeStack(data->st);
data->line=seq;
line=seq;
fc=fc2;
MAX=strlen(seq);
createMatrix(strlen(seq),con);
//displayMatrix(strlen(seq));
/* allocate memory for MFE structure (length + 1) */
char *structure = (char *)vrna_alloc(sizeof(char) * (strlen(seq) + 1));
/* predict Minmum Free Energy and corresponding secondary structure */
vrna_ud_set_data(fc2,data,free_data);
vrna_ud_set_prod_rule_cb(fc2, test,custom_ud_exp_prod_cb);

float mfe = vrna_mfe(fc2, structure);

printf("%s\n%s [ %6.2f ]\n", seq,structure, mfe);

/* cleanup memory */
//free(seq);
//free(structure);
}
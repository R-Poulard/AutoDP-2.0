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
#include "H.c"

int custom_ud_exp_prod_cb(vrna_fold_compound_t *fc,
                                         int i,
                                         int j,
                                         void *data) {
    //printf("%d %d\n",i,j);
  pk_compound *pk=(pk_compound*)data;
  int tmp;
  Node* first;
  Node* last;
  if(get_A(pk,&first,&last,&tmp,i-1,j-1)==0){
    init_root_inter(pk,i-1,j-1);
    fold(pk);
    get_A(pk,&first,&last,&tmp,i-1,j-1);
  }
  //printf("got %d %d,%d\n",tmp,i,j);
  return tmp;
}

void free_aux(void * data){
  pk_compound *pk=(pk_compound*)data;
  free_pk(data);
}

int
main()
{
/* initialize random number generator */
vrna_init_rand();
/* Generate a random sequence of 50 nucleotides */
const char *seq = "CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUC";
const char *con = "......................................"; // Hard constraint dot-bracket notation

vrna_fold_compound_t *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);
vrna_constraints_add(fc, con, VRNA_CONSTRAINT_DB_DEFAULT);
pk_compound * pk=create_pk_compound(fc,seq,con);

//displayMatrix(strlen(seq));
/* allocate memory for MFE structure (length + 1) */
char *structure = (char *)vrna_alloc(sizeof(char) * (strlen(seq) + 1));
/* predict Minmum Free Energy and corresponding secondary structure */
vrna_gr_set_data(fc,pk,free_aux);
vrna_gr_set_aux(fc,custom_ud_exp_prod_cb);
vrna_gr_set_aux_f(fc,custom_ud_exp_prod_cb);
vrna_gr_set_aux_c(fc,custom_ud_exp_prod_cb);
vrna_gr_set_aux_m(fc,custom_ud_exp_prod_cb);
vrna_gr_set_aux_m1(fc,custom_ud_exp_prod_cb);
float mfe = vrna_mfe(fc, structure);

printf("%s\n%s [ %6.2f ]\n", seq,structure, mfe);

/* cleanup memory */
//free(seq);
//free(structure);
}
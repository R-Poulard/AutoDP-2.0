#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>
#include <ViennaRNA/mfe.h>
int
main()
{
/* initialize random number generator */
vrna_init_rand();
/* Generate a random sequence of 50 nucleotides */
const char *seq = "CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG";
const char *con = "......................................................................."; // Hard constraint dot-bracket notation

/* Create a fold compound for the sequence */
vrna_fold_compound_t *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);
vrna_constraints_add(fc, con, VRNA_CONSTRAINT_DB_DEFAULT);

/* allocate memory for MFE structure (length + 1) */
char *structure = (char *)vrna_alloc(sizeof(char) * (strlen(seq) + 1));
/* predict Minmum Free Energy and corresponding secondary structure */
float mfe = vrna_mfe(fc, structure);

int * indx = fc->iindx;
int i=1;
int j=60;
int ij = indx[j++] + i++;
int * fML = fc->matrices->fML;
printf("fML[%d] = %d (ou %d)\n",ij,fML[ij],fc->matrices->fML[fc->iindx[j++]+(i++)]);

vrna_bp_stack_t * bp=(vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + strlen(seq) / 2)));
int b=bp[0].i;
int p,q,comp1,comp2;
int lol=vrna_BT_mb_loop_split(fc,&i,&j,&p,&q,&comp1,&comp2,bp,&b);
char *structure2 = (char *)vrna_alloc(sizeof(char) * (strlen(seq) + 1));


structure2 = vrna_db_from_bp_stack(bp, strlen(seq));

printf("%s\n%s [ %6.2f ]\n BT=%d\n %s\n", seq,structure, mfe,lol,structure2);

/* cleanup memory */
//free(seq);
//free(structure);

vrna_fold_compound_free(fc);
return 0;
}

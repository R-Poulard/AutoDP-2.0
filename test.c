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

void printing(vrna_bp_stack_t *bp,unsigned int    length){
  int   k, i, j, temp;
  //printf("iciii\n");
  if (bp) {
    printf("BP FOUND:\n");
    for (k = 1; k <= bp[0].i; k++) {
      i = bp[k].i;
      j = bp[k].j;
      if (i > length)
        i -= length;

      if (j > length)
        j -= length;

      if (i > j) {
        temp  = i;
        i     = j;
        j     = temp;
      }

      if (i == j) {
        /* Gquad bonds are marked as bp[i].i == bp[i].j */
        //structure[i - 1] = '+';
        printf("bp %d %d = +",i-1,i-1);
      } else {
        /* the following ones are regular base pairs */
        //structure[i - 1]  = '(';
        //structure[j - 1]  = ')';
        printf("bp %d %d = (",i-1,j-1);
        printf("bp %d %d = )",i-1,j-1);
      }
    }
    printf("-----------\n");
  }
}

int
main()
{
/* initialize random number generator */
vrna_init_rand();
/* Generate a random sequence of 50 nucleotides */
const char *seq = "CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG";
//const char *con = "......................................................................."; // Hard constraint dot-bracket notation
//char* seq="CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUAC";
/* Create a fold compound for the sequence */
vrna_fold_compound_t *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);
//vrna_constraints_add(fc, con, VRNA_CONSTRAINT_DB_DEFAULT);

/* allocate memory for MFE structure (length + 1) */
char *structure = (char *)vrna_alloc(sizeof(char) * (strlen(seq) + 1));
/* predict Minmum Free Energy and corresponding secondary structure */
float mfe = vrna_mfe(fc, structure);

int * indx = fc->jindx;
int i=4;
int j=40;
int ij = indx[j] + (i);
int * fML = fc->matrices->fML;
printf("fML[%d] = %d (%d)\n",ij,fML[ij],fc->matrices->fML[fc->jindx[j]+(i)]);

vrna_bp_stack_t * bp=(vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + strlen(seq) / 2)));

int p,q,comp1,comp2;
int lol=backtrack_MFE(fc,i,j,bp);
char *structure2 = (char *)vrna_alloc(sizeof(char) * (strlen(seq) + 1));
int b=bp[0].i;
structure2= vrna_db_from_bp_stack(bp,strlen(seq));
printing(bp, strlen(seq));
//exit(0);
printf("i=%d j=%d p=%d q=%d comp1=%d comp2=%d\n",i,j,p,q,comp1,comp2);
printf("%s\n%s [ %6.2f ]\nb_value= %d BT=%d\n%s\n", seq,structure, mfe,b,lol,structure2);
printf("%c\n",fc->sequence[0]);

/* cleanup memory */
//free(seq);
//free(structure);

vrna_fold_compound_free(fc);
return 0;
}

int backtrack_MFE(vrna_fold_compound_t *fc, int i,int j, vrna_bp_stack_t *bp_stack){
    int ret;
    int p,q,comp1,comp2;
    int s=0;
    int b=bp_stack[0].i;
    char * ss;
    sect              bt_stack[500]; /* stack of partial structures for backtracking */

    int BT=vrna_BT_mb_loop_split(fc, &i, &j, &p, &q, &comp1, &comp2, bp_stack, &b);
    if(BT!=1){
        return BT;
    }
    if (i > 0) {
        bt_stack[++s].i = i;
        bt_stack[s].j   = j;
        bt_stack[s].ml  = comp1;
    }

    if (p > 0) {
        bt_stack[++s].i = p;
        bt_stack[s].j   = q;
        bt_stack[s].ml  = comp2;
    }
    BT=my_backtrack(fc,bp_stack,bt_stack,s,NULL);

    free(ss);
    return BT;
}

int my_backtrack(vrna_fold_compound_t  *fc,
          vrna_bp_stack_t       *bp_stack,
          sect                  bt_stack[],
          int                   s,
          struct ms_helpers     *ms_dat)
{
  char          backtrack_type;
  int           i, j, ij, k, l, length, b, *my_c, *indx, noLP, *pscore, ret;
  vrna_param_t  *P;

  ret             = 1;
  b               = bp_stack[0].i;
  length          = fc->length;
  my_c            = fc->matrices->c;
  indx            = fc->jindx;
  P               = fc->params;
  noLP            = P->model_details.noLP;
  pscore          = fc->pscore;         /* covariance scores for comparative structure prediction */
  backtrack_type  = P->model_details.backtrack_type;

  if (s == 0) {
    bt_stack[++s].i = 1;
    bt_stack[s].j   = length;
    bt_stack[s].ml  = (backtrack_type == 'M') ? 1 : ((backtrack_type == 'C') ? 2 : 0);
  }

  while (s > 0) {
    int ml, cij;
    int canonical = 1;     /* (i,j) closes a canonical structure */

    /* pop one element from stack */
    i   = bt_stack[s].i;
    j   = bt_stack[s].j;
    ml  = bt_stack[s--].ml;

    switch (ml) {
      /* backtrack in f5 */
      case 0:
      {
        int p, q;
        if (vrna_BT_ext_loop_f5(fc, &j, &p, &q, bp_stack, &b)) {
          if (j > 0) {
            bt_stack[++s].i = 1;
            bt_stack[s].j   = j;
            bt_stack[s].ml  = 0;
          }

          if (p > 0) {
            i = p;
            j = q;
            goto repeat1;
          }

          continue;
        } else {
          vrna_message_warning("backtracking failed in f5, segment [%d,%d], e = %d\n",
                               i,
                               j,
                               fc->matrices->f5[j]);
          ret = 0;
          goto backtrack_exit;
        }
      }
      break;

      /* trace back in fML array */
      case 1:
      {
        int p, q, comp1, comp2;
        if (vrna_BT_mb_loop_split(fc, &i, &j, &p, &q, &comp1, &comp2, bp_stack, &b)) {
          if (i > 0) {
            bt_stack[++s].i = i;
            bt_stack[s].j   = j;
            bt_stack[s].ml  = comp1;
          }

          if (p > 0) {
            bt_stack[++s].i = p;
            bt_stack[s].j   = q;
            bt_stack[s].ml  = comp2;
          }

          continue;
        } else {
          vrna_message_warning("backtracking failed in fML, segment [%d,%d]\n", i, j);
          ret = 0;
          goto backtrack_exit;
        }
      }
      break;

      /* backtrack in c */
      case 2:
        bp_stack[++b].i = i;
        bp_stack[b].j   = j;
        goto repeat1;

        break;

      /* backtrack in fms5 */
      case 5:
      {
        unsigned int strand = j;

        if (BT_fms5_split(fc, strand, &i, &k, &l, ms_dat)) {
          if (k > 0) {
            bt_stack[++s].i = i;
            bt_stack[s].j   = k;
            bt_stack[s].ml  = 2;

            if (k < fc->strand_end[strand]) {
              bt_stack[++s].i = l;
              bt_stack[s].j   = strand;
              bt_stack[s].ml  = 5;
            }
          } else if (i > 0) {
            bt_stack[++s].i = i;
            bt_stack[s].j   = strand;
            bt_stack[s].ml  = 5;
          }

          continue;
        } else {
          vrna_message_warning("backtracking failed in fsm5[%d][%d] (%d:%d)\n",
                               strand,
                               i,
                               fc->strand_start[strand],
                               fc->strand_end[strand]);
          ret = 0;
          goto backtrack_exit;
        }
      }
      break;

      /* backtrack in fms3 */
      case 6:
      {
        unsigned int strand = i;

        if (BT_fms3_split(fc, strand, &j, &k, &l, ms_dat)) {
          if (k > 0) {
            bt_stack[++s].i = k;
            bt_stack[s].j   = j;
            bt_stack[s].ml  = 2;

            if (k > fc->strand_start[strand]) {
              bt_stack[++s].i = strand;
              bt_stack[s].j   = l;
              bt_stack[s].ml  = 6;
            }
          } else if (j > 0) {
            bt_stack[++s].i = strand;
            bt_stack[s].j   = j;
            bt_stack[s].ml  = 6;
          }

          continue;
        } else {
          vrna_message_warning("backtracking failed in fsm3[%d][%d] (%d:%d)\n",
                               strand,
                               j,
                               fc->strand_start[strand],
                               fc->strand_end[strand]);
          ret = 0;
          goto backtrack_exit;
        }
      }
      break;

      default:
        ret = 0;
        goto backtrack_exit;
    }

repeat1:

    /*----- begin of "repeat:" -----*/
    ij = indx[j] + i;

    if (canonical)
      cij = my_c[ij];

    if (noLP) {
      if (vrna_BT_stack(fc, &i, &j, &cij, bp_stack, &b)) {
        canonical = 0;
        goto repeat1;
      }
    }

    canonical = 1;

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
      cij += pscore[indx[j] + i];

    if (vrna_BT_hp_loop(fc, i, j, cij, bp_stack, &b))
      continue;

    if (vrna_BT_int_loop(fc, &i, &j, cij, bp_stack, &b)) {
      if (i < 0)
        continue;
      else
        goto repeat1;
    }

    /* (i.j) must close a multi-loop */
    int comp1, comp2;

    if (vrna_BT_mb_loop(fc, &i, &j, &k, cij, &comp1, &comp2)) {
      bt_stack[++s].i = i;
      bt_stack[s].j   = k;
      bt_stack[s].ml  = comp1;
      bt_stack[++s].i = k + 1;
      bt_stack[s].j   = j;
      bt_stack[s].ml  = comp2;
    } else {
      vrna_message_warning("backtracking failed in repeat, segment [%d,%d]\n", i, j);
      ret = 0;
      goto backtrack_exit;
    }

    /* end of repeat: --------------------------------------------------*/
  } /* end of infinite while loop */

backtrack_exit:

  bp_stack[0].i = b;    /* save the total number of base pairs */

  return ret;
}
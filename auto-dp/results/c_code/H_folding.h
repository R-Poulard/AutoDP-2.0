#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <string.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/loop_energies.h>
#include "aux.h"

int setup_H(vrna_fold_compound_t *fc,vrna_fold_compound_t *fc_aux,const char* seq, const char* con,bt_struct* BT);


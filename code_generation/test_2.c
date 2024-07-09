#include <stdio.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/default.h>
#include <ViennaRNA/constraints.h>
#include <ViennaRNA/constraints_hard.h>
#include <ViennaRNA/constraints_soft.h>
#include <ViennaRNA/energy_par.h>
#include <ViennaRNA/eval.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/unstructured_domains.h>
#include <string.h>
#include <H.c>

typedef struct {
    HashTable *ht;
    Stack *st;
    int** matrix;
    char *line;
}Custom_data;
/* Define your custom scoring callback function */
 int custom_ud_exp_prod_cb(vrna_fold_compound_t *fc,
                                         int i,
                                         int j,
                                         unsigned int loop_type,
                                         void *data) {
    if(i==1 && j==5){
        return -100;
    }
    return 0;
}

int main(void) {
    const char *sequence = "GCGCUUCGCCAGCUCUGUACGGCUA";
    const char *ss = ".........................";
    vrna_fold_compound_t *fc;


    // Initialize the fold compound
    fc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_EVAL_ONLY);
    HashTable *hashTable = createHashTable();
    Stack *stack = malloc(sizeof(Stack));
    initializeStack(stack);
    Custom_data *data=malloc(sizeof(Custom_data));
    data->st=stack;
    data->ht=hashTable;
    data->line=sequence;
    createMatrix(strlen(sequence),ss);
    displayMatrix(strlen(sequence));
    // Initialize unstructured domains
    //domains_up = vrna_ud_init(fc);

    // Set custom unstructured domain callback
    //vrna_ud_set_prod_rule_cb(fc, NULL,custom_ud_exp_prod_cb);

    // Perform MFE folding
    char *structure = (char *)vrna_alloc(sizeof(char) * (26 + 1));
    double mfe = vrna_mfe(fc, structure);

    printf("MFE: %6.2f\n%s\n%s\n", mfe, sequence, structure);

    // Cleanup
    free(structure);
    //vrna_ud_destroy(domains_up);
    vrna_fold_compound_free(fc);

    return 0;
}

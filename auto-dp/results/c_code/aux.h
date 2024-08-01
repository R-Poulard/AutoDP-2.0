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

#ifndef AUX_H
#define AUX_H

#define LOAD_FACTOR_THRESHOLD 0.7

PUBLIC typedef struct{
    char * structure;
    char bracket;
    int MAX;
    vrna_bts_t bt_stack;
    vrna_bps_t bp_stack;
} bt_struct;


bt_struct * create_bt_struct(char* structure,int MAX);
void free_bt(bt_struct* bt);


//PAIRING FUNCTIONS
int min(int a, int b);
int max(int a, int b);

int add(int a,int b);
int mult(int a,int b);

double random_double(double b);

void createMatrix(int *** matrix, int N,const char * ss);
// Function to display the matrix
void displayMatrix(int** matrix,int N);

// Function to free the memory allocated for the matrix
void freeMatrix(int ** matrix,int N);


int matching(char x, char y);


char * post_process_f(vrna_fold_compound_t *fc,vrna_bps_t bp_stack,void *data);

#endif // AUX_H

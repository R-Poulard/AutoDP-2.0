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

#define LOAD_FACTOR_THRESHOLD 0.7

PUBLIC typedef struct{
    char * structure;
    char bracket;
    int MAX;
    vrna_bts_t bt_stack;
    vrna_bps_t bp_stack;
} bt_struct;


bt_struct * create_bt_struct(char* structure,int MAX){
    bt_struct* bt=malloc(sizeof(bt_struct));
    bt->bracket='A';
    bt->MAX=MAX;
    bt->bt_stack=NULL;
    bt->bp_stack=NULL;
    if(structure==NULL){
        bt->structure=malloc(sizeof(char)*(MAX+1));
        bt->structure[MAX]='\0';
        for(int i=0;i<MAX;i++){
            bt->structure[i]='.';
        }
    }
    else{
        bt->structure=structure;
    }
    return bt;
}
void free_bt(bt_struct* bt){
    free(bt->structure);
    free(bt);
}


//PAIRING FUNCTIONS
int min(int a, int b) { if (a<b) {return a;} else {return b;}};
int max(int a, int b) { if (a<b) {return b;} else {return a;}};

int add(int a,int b){
    if(a==INF || b==INF){
        return INF;
    }

    return a+b;
}

int mult(int a,int b){
    if(a==INF || b==INF){
        return INF;
    }

    return a*b;
}

double random_double(double b) {
    // Generate a random double between 0.0 and 1.0
    double random = (double)rand() / RAND_MAX;
    // Scale it to the desired range [0.0, b]
    return random * b;
}

void createMatrix(int *** matrix, int N,const char * ss) {

    // Allocate memory for the matrix
    *matrix = (int **)malloc(N * sizeof(int *));
    for (int i = 0; i < N; i++) {
        (*matrix)[i] = (int *)malloc(N * sizeof(int));
    }

    // Fill the matrix
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if(i==j){
                if(ss[i]=='|'){
                    (*matrix)[i][j]=0;
                }
                else{
                    (*matrix)[i][j]=1;
                }
                continue;
            }
            (*matrix)[i][j] = 0;
            (*matrix)[j][i] = 0;
            if ( ss[i]=='x' || ss[j]=='x') {
                (*matrix)[i][j] = 0;
                (*matrix)[j][i] = 0;
                continue;
            }
            if( (ss[i]=='(' && ss[j]==')') || (ss[j]=='(' && ss[i]==')')){
                (*matrix)[i][j]= 1;
                (*matrix)[j][i] = 1;
                continue;
            }
            if( ss[i]=='>'){
                (*matrix)[i][j]= 0;
                (*matrix)[j][i]= 0;
                continue;
            }
            if( (ss[i]=='.' || ss[i]=='<') && (ss[j]=='.' || ss[j]=='>')){
                (*matrix)[i][j]= 1;
                (*matrix)[j][i]= 1;
                continue;
            }
            if(ss[i]=='['){
                if(ss[j]==']' && j>i){
                    int matching=0;
                    for(int tmp=i+1;tmp<j;tmp++){
                        if(ss[tmp]==']'){
                            matching--;
                        }
                        if(ss[tmp]=='['){
                            matching++;
                        }
                    }
                    if(matching==0){
                        (*matrix)[i][j]=1;
                        (*matrix)[j][i]=1;
                    }
                    else{
                        (*matrix)[i][j]=0;
                        (*matrix)[j][i]=0;
                    }
                }
                else{
                    (*matrix)[i][j] = 0;
                    (*matrix)[j][i] = 0;
                }
                continue;
            }
        }
    }
}

// Function to display the matrix
void displayMatrix(int** matrix,int N) {
    printf("Matrix:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

// Function to free the memory allocated for the matrix
void freeMatrix(int ** matrix,int N) {
    for (int i = 0; i < N; i++) {
        free(matrix[i]);
    }
    free(matrix);
}


int matching(char x, char y){
    if(x!='X' && x!='G' && x!='C' && x!='A' && x!='U' && x!='T'){
        printf("x is not a recognizable charactere\n");
        exit(-1);
    }
    if(y!='X' && y!='G' && y!='C' && y!='A' && y!='U' && y!='T'){
        printf("y is not a recognizable charactere\n");
        exit(-1);
    }
    //printf("%c,%c\n",x,y);
    if (x=='G' && y=='C') { return -10; }
    if (x=='C' && y=='G') { return -10; }
    if (x=='G' && y=='U') { return -5; }
    if (x=='U' && y=='G') { return -5; }
    if (x=='A' && y=='U') { return -5; }
    if (x=='U' && y=='A') { return -5; }
    return INF;
}


char * post_process_f(vrna_fold_compound_t *fc,vrna_bps_t bp_stack,void *data){
  bt_struct *bt=(bt_struct*)data;
  char * ss=vrna_db_from_bps(bp_stack,bt->MAX);
  for(int i=0;i<bt->MAX;i++){
    if(bt->structure[i]!='.'){
      ss[i]=bt->structure[i];
    }
  }
  return ss;
}
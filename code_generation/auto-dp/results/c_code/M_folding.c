#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdbool.h>

int min(int a, int b) { if (a<b) {return a;} else {return b;}};

int index_B(int a,int e,int f,int h);
int index_C(int b,int d,int f,int h);
int index_CLIQUE(int i, int j, int k, int l);

void init_fill_B();
void init_fill_C();
void init_fill_CLIQUE();

int compute_A();
int compute_B(int a,int e,int f,int h);
int compute_C(int b,int d,int f,int h);
int fold();

int n;

char * line = NULL;

double * B;
double * C;
double * CLIQUE;
double * CLIQUE2;

int main(int argc, char ** argv) {
    char * line = NULL;
    size_t len = 0;
    FILE * fp = fopen(argv[1], "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);
    getline(&line, &len, fp); 
    printf("%s", line);
    double * B = malloc(n*n*n*n*sizeof(double));
    double * C = malloc(n*n*n*n*sizeof(double));
    double * CLIQUE = malloc(n*n*n*n*sizeof(double));
    init_fill_B();
    init_fill_C();
    init_fill_CLIQUE();
    int score = fold();
    char * structure = NULL;
    structure = backtrace();
    free(B);
    free(C);
    free(CLIQUE);
}
int bp_score(char x, char y) {
    if (x=='G' && y=='C') { return 10; }
    if (x=='C' && y=='G') { return 10; }
    if (x=='G' && y=='U') { return 5; }
    if (x=='U' && y=='G') { return 5; }
    if (x=='A' && y=='U') { return 5; }
    if (x=='U' && y=='A') { return 5; }
    return 0;
}
int fold() {
    compute_A();
    }

int compute_CLIQUE2(int i, int j, int k, int l);

int compute_CLIQUE(int i, int j, int k, int l) {
    //C_boxtimes in article"
    if (CLIQUE[index_CLIQUE(i,j,k,l)] > INT_MIN) { 
        return CLIQUE[index_CLIQUE(i,j,k,l)];
    }
    
    int min_value = INT_MAX;

    if (k < l) { 
        min_value = min(min_value, compute_CLIQUE2(i,j,k,l-1)); 
    }
    if (i < j) {
        min_value = min(min_value, compute_CLIQUE(i+1,j,k,l));
    }
    if (i < j && k < l) {
        min_value = min(min_value, 
                    compute_CLIQUE(i+1, j, k, l-1)+bp_score(line[i],line[l]));
    }
    if (k==l) { 
        min_value = min(min_value, bp_score(line[i], line[l]));
    }

    CLIQUE[index_CLIQUE(i,j,k,l)] = min_value;
    return min_value;
}

int compute_CLIQUE2(int i, int j, int k, int l) {
    //C_boxtimes' in article:w

    if (CLIQUE2[index_CLIQUE(i,j,k,l)] > INT_MIN) { 
        return CLIQUE2[index_CLIQUE(i,j,k,l)];
    }
    int min_value = INT_MAX;

    if (k < l) { 
        min_value = min(min_value, compute_CLIQUE2(i,j,k,l-1)); 
    }
    if (i < j && k < l) {
        min_value = min(min_value, 
                    compute_CLIQUE(i+1, j, k, l-1)+bp_score(line[i],line[l]));
    }
    if (k==l) { 
        min_value = min(min_value, bp_score(line[i], line[l]));
    }

    CLIQUE2[index_CLIQUE(i,j,k,l)] = min_value;
    return min_value;
}
void init_fill_B() {
    for (int a=0;a<n;a++) {
        for (int e=a;e<n;e++) {
            for (int f=e;f<n;f++) {
                for (int h=f;h<n;h++) {
                    B[index_B(a,e,f,h)] = INT_MIN;
                }
            }
        }
    }
}

int index_B(int a,int e,int f,int h)  {
    return n*n*n*a+n*n*e+n*f+h;
}

void init_fill_C() {
    for (int b=0;b<n;b++) {
        for (int d=b;d<n;d++) {
            for (int f=d;f<n;f++) {
                for (int h=f;h<n;h++) {
                    C[index_C(b,d,f,h)] = INT_MIN;
                }
            }
        }
    }
}

int index_C(int b,int d,int f,int h)  {
    return n*n*n*b+n*n*d+n*f+h;
}

void init_fill_CLIQUE() {
    for (int i=0; i < n;i++) {
        for (int j=i; j < n;j++) {
            for (int k=j; k < n;k++) {
                for (int l=k; l < n;l++) {
                    CLIQUE[i,j,k,l] = INT_MIN;
                }
            }
        }
    }
}
int index_CLIQUE(int i, int j, int k, int l) {
    return n*n*n*i+n*n*j+n*k+l;
}

int compute_A() {
    
    int min_value = INT_MAX;
    
    for (int a=0;a<n;a++) {
        for (int e=a;e<n;e++) {
            for (int f=e;f<n;f++) {
                for (int h=f;h<n;h++) {
                    for (int i=h;i<n;i++) {
                        min_value = min(min_value, compute_CLIQUE(e,f,h,i)+compute_B(e,h,f,a));
                    }
                }
            }
        }
    }

    return min_value;
}


int compute_B(int a,int e,int f,int h) {
    if (B[index_B(a,e,f,h)] > INT_MIN) {
        return B[index_B(a,e,f,h)];
    }
    
    int min_value = INT_MAX;
    
    for (int d=a;d<h;d++) {
        for (int b=d;b<n;b++) {
            min_value = min(min_value, compute_C(h,b,f,d)+compute_CLIQUE(a,b,d,e));
        }
    }

    B[index_B(a,e,f,h)] = min_value;
    return min_value;
}

int compute_C(int b,int d,int f,int h) {
    if (C[index_C(b,d,f,h)] > INT_MIN) {
        return C[index_C(b,d,f,h)];
    }
    
    int min_value = INT_MAX;
    
    for (int g=f;g<h;g++) {
        for (int c=g;c<n;c++) {
            min_value = min(min_value, compute_CLIQUE(c,d,g,h)+compute_CLIQUE(b,c,f,g));
        }
    }

    C[index_C(b,d,f,h)] = min_value;
    return min_value;
}




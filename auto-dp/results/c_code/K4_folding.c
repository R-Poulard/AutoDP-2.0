#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdbool.h>

int min(int a, int b) { if (a<b) {return a;} else {return b;}};

int index_B(int a,int d,int e,int h);
int index_C(int a,int c,int e,int g);
int index_CLIQUE(int i, int j, int k, int l);

void init_fill_B();
void init_fill_C();
void init_fill_CLIQUE();

int compute_A();
int compute_B(int a,int d,int e,int h);
int compute_C(int a,int c,int e,int g);
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
        for (int d=a;d<n;d++) {
            for (int e=d;e<n;e++) {
                for (int h=e;h<n;h++) {
                    B[index_B(a,d,e,h)] = INT_MIN;
                }
            }
        }
    }
}

int index_B(int a,int d,int e,int h)  {
    return n*n*n*a+n*n*d+n*e+h;
}

void init_fill_C() {
    for (int a=0;a<n;a++) {
        for (int c=a;c<n;c++) {
            for (int e=c;e<n;e++) {
                for (int g=e;g<n;g++) {
                    C[index_C(a,c,e,g)] = INT_MIN;
                }
            }
        }
    }
}

int index_C(int a,int c,int e,int g)  {
    return n*n*n*a+n*n*c+n*e+g;
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
        for (int d=a;d<n;d++) {
            for (int e=d;e<n;e++) {
                for (int h=e;h<n;h++) {
                    for (int i=h;i<n;i++) {
                        min_value = min(min_value, compute_B(a,d,h,e)+compute_CLIQUE(d,e,h,i));
                    }
                }
            }
        }
    }

    return min_value;
}

int compute_B(int a,int d,int e,int h) {
    if (B[index_B(a,d,e,h)] > INT_MIN) {
        return B[index_B(a,d,e,h)];
    }
    
    int min_value = INT_MAX;
    
    for (int g=e;g<h;g++) {
        for (int c=g;c<n;c++) {
            min_value = min(min_value, compute_C(g,a,c,e)+compute_CLIQUE(c,d,g,h));
        }
    }

    B[index_B(a,d,e,h)] = min_value;
    return min_value;
}

int compute_C(int a,int c,int e,int g) {
    if (C[index_C(a,c,e,g)] > INT_MIN) {
        return C[index_C(a,c,e,g)];
    }
    
    int min_value = INT_MAX;
    
    for (int f=e;f<g;f++) {
        for (int b=f;b<c;b++) {
            min_value = min(min_value, compute_CLIQUE(a,b,e,f)+compute_CLIQUE(b,c,f,g));
        }
    }

    C[index_C(a,c,e,g)] = min_value;
    return min_value;
}





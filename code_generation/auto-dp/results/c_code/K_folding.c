#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdbool.h>

int min(int a, int b) { if (a<b) {return a;} else {return b;}};

int index_B(int d,int g,int a,int d2);
int index_C(int a,int d,int e,int f);
int index_CLIQUE(int i, int j, int k, int l);

void init_fill_B();
void init_fill_C();
void init_fill_CLIQUE();

int compute_A();
int compute_B(int d,int g,int a,int d2);
int compute_B2(int d,int g,int a,int d2);
int compute_C(int a,int d,int e,int f);
int compute_C2(int a,int d,int e,int f);
int fold();

int n;

char * line = NULL;

double * B;
double * B2;
double * C;
double * C2;
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
double * B2 = malloc(n*n*n*n*sizeof(double));
    double * C = malloc(n*n*n*n*sizeof(double));
double * C2 = malloc(n*n*n*n*sizeof(double));
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
    for (int d=0;d<n;d++) {
        for (int g=d;g<n;g++) {
            for (int a=g;a<n;a++) {
                for (int d=a;d<n;d++) {
                    B[index_B(d,g,a,d)] = INT_MIN;
                }
            }
        }
    }
}

int index_B(int d,int g,int a,int d2)  {
    return n*n*n*d+n*n*g+n*a+d2;
}

void init_fill_C() {
    for (int a=0;a<n;a++) {
        for (int d=a;d<n;d++) {
            for (int e=d;e<n;e++) {
                for (int f=e;f<n;f++) {
                    C[index_C(a,d,e,f)] = INT_MIN;
                }
            }
        }
    }
}

int index_C(int a,int d,int e,int f)  {
    return n*n*n*a+n*n*d+n*e+f;
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
            for (int g=d;g<n;g++) {
                min_value = min(min_value, compute_B(d,g,a,d));
            }
        }
    }

    return min_value;
}

int compute_B(int g, int d, int a,int d2) {
    if (B[index_B(g,d, a,d2)] > INT_MIN) {
        return B[index_B(g,d, a,d2)];
    }

    int min_value = INT_MAX;
    bool eq_some_const1 = (d-1!=a) && (d-1!=d2);
    bool eq_some_const2 = (d-1!=a) && (d-1!=d2);

    if (g+1 < d) {
        if (!eq_some_const1) {
            min_value = min(min_value, compute_B(g+1,d,a,d2));
        }
    }
    if (d-1 > g) {
        if (!eq_some_const2) {
            min_value = min(min_value, compute_B2(g, d-1, a,d2));
        }
    }
    if (!eq_some_const1 && !eq_some_const2) {
        min_value = min(min_value, compute_B(g+1, d-1, a,d2)+bp_score(line[g], 
                                                                                      line[d]));
    }

    min_value = min(min_value, compute_C(a,d2,13,25));

    B[index_B(g,d,a,d2)] = min_value;
    return min_value;
} 

int compute_B2(int g, int d, int a,int d2) {
    if (B2[index_B(g,d, a,d2)] > INT_MIN) {
        return B2[index_B(g,d, a,d2)];
    }

    int min_value = INT_MAX;
    bool eq_some_const1 = (d-1!=a) && (d-1!=d2);
    bool eq_some_const2 = (d-1!=a) && (d-1!=d2);

    if (d-1 > g && !eq_some_const2) {
        min_value = min(min_value, compute_B2(g, d-1, a,d2));
    }
    if (!eq_some_const1 && !eq_some_const2) {
        min_value = min(min_value, compute_B(g+1, d-1, a,d2)+bp_score(line[g], 
                                                                                      line[d]));
    }

    B2[index_B(g,d,a,d2)] = min_value;
    return min_value;
}


int compute_C(int a, int d, int f,int e) {
    if (C[index_C(a,d, f,e)] > INT_MIN) {
        return C[index_C(a,d, f,e)];
    }

    int min_value = INT_MAX;
    bool eq_some_const1 = (d-1!=f) && (d-1!=e);
    bool eq_some_const2 = (d-1!=f) && (d-1!=e);

    if (a+1 < d) {
        if (!eq_some_const1) {
            min_value = min(min_value, compute_C(a+1,d,f,e));
        }
    }
    if (d-1 > a) {
        if (!eq_some_const2) {
            min_value = min(min_value, compute_C2(a, d-1, f,e));
        }
    }
    if (!eq_some_const1 && !eq_some_const2) {
        min_value = min(min_value, compute_C(a+1, d-1, f,e)+bp_score(line[a], 
                                                                                      line[d]));
    }

    min_value = min(min_value, compute_CLIQUE(1,13,e,f));

    C[index_C(a,d,f,e)] = min_value;
    return min_value;
} 

int compute_C2(int a, int d, int f,int e) {
    if (C2[index_C(a,d, f,e)] > INT_MIN) {
        return C2[index_C(a,d, f,e)];
    }

    int min_value = INT_MAX;
    bool eq_some_const1 = (d-1!=f) && (d-1!=e);
    bool eq_some_const2 = (d-1!=f) && (d-1!=e);

    if (d-1 > a && !eq_some_const2) {
        min_value = min(min_value, compute_C2(a, d-1, f,e));
    }
    if (!eq_some_const1 && !eq_some_const2) {
        min_value = min(min_value, compute_C(a+1, d-1, f,e)+bp_score(line[a], 
                                                                                      line[d]));
    }

    C2[index_C(a,d,f,e)] = min_value;
    return min_value;
}



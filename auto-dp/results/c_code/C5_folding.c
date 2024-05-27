#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdbool.h>

int min(int a, int b) { if (a<b) {return a;} else {return b;}};

int index_B(int a,int m,int n,int r);
int index_C(int a,int l,int n,int r);
int index_D(int a,int h,int k,int l,int p,int r);
int index_E(int b,int g,int k,int l,int p,int r);
int index_F(int b,int e,int p,int r);
int index_G(int c,int e,int p,int r);
int index_H(int c,int e,int q,int r);
int index_I(int b,int e,int g,int k,int l);
int index_J(int a,int h,int k,int n,int p);
int index_K(int i,int k,int n,int p);
int index_L(int i,int k,int o,int p);
int index_M(int m,int n,int r,int t);
int index_CLIQUE(int i, int j, int k, int l);

void init_fill_B();
void init_fill_C();
void init_fill_D();
void init_fill_E();
void init_fill_F();
void init_fill_G();
void init_fill_H();
void init_fill_I();
void init_fill_J();
void init_fill_K();
void init_fill_L();
void init_fill_M();
void init_fill_CLIQUE();

int compute_A();
int compute_B(int a,int m,int n,int r);
int compute_C(int a,int l,int n,int r);
int compute_D(int a,int h,int k,int l,int p,int r);
int compute_D2(int a,int h,int k,int l,int p,int r);
int compute_E(int b,int g,int k,int l,int p,int r);
int compute_F(int b,int e,int p,int r);
int compute_G(int c,int e,int p,int r);
int compute_H(int c,int e,int q,int r);
int compute_I(int b,int e,int g,int k,int l);
int compute_J(int a,int h,int k,int n,int p);
int compute_K(int i,int k,int n,int p);
int compute_L(int i,int k,int o,int p);
int compute_M(int m,int n,int r,int t);
int fold();

int n;

char * line = NULL;

double * B;
double * C;
double * D;
double * D2;
double * E;
double * F;
double * G;
double * H;
double * I;
double * J;
double * K;
double * L;
double * M;
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
    double * D = malloc(n*n*n*n*n*n*sizeof(double));
double * D2 = malloc(n*n*n*n*n*n*sizeof(double));
    double * E = malloc(n*n*n*n*n*n*sizeof(double));
    double * F = malloc(n*n*n*n*sizeof(double));
    double * G = malloc(n*n*n*n*sizeof(double));
    double * H = malloc(n*n*n*n*sizeof(double));
    double * I = malloc(n*n*n*n*n*sizeof(double));
    double * J = malloc(n*n*n*n*n*sizeof(double));
    double * K = malloc(n*n*n*n*sizeof(double));
    double * L = malloc(n*n*n*n*sizeof(double));
    double * M = malloc(n*n*n*n*sizeof(double));
    double * CLIQUE = malloc(n*n*n*n*sizeof(double));
    init_fill_B();
    init_fill_C();
    init_fill_D();
    init_fill_E();
    init_fill_F();
    init_fill_G();
    init_fill_H();
    init_fill_I();
    init_fill_J();
    init_fill_K();
    init_fill_L();
    init_fill_M();
    init_fill_CLIQUE();
    int score = fold();
    char * structure = NULL;
    structure = backtrace();
    free(B);
    free(C);
    free(D);
    free(E);
    free(F);
    free(G);
    free(H);
    free(I);
    free(J);
    free(K);
    free(L);
    free(M);
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
        for (int m=a;m<n;m++) {
            for (int n=m;n<n;n++) {
                for (int r=n;r<n;r++) {
                    B[index_B(a,m,n,r)] = INT_MIN;
                }
            }
        }
    }
}

int index_B(int a,int m,int n,int r)  {
    return n*n*n*a+n*n*m+n*n+r;
}

void init_fill_C() {
    for (int a=0;a<n;a++) {
        for (int l=a;l<n;l++) {
            for (int n=l;n<n;n++) {
                for (int r=n;r<n;r++) {
                    C[index_C(a,l,n,r)] = INT_MIN;
                }
            }
        }
    }
}

int index_C(int a,int l,int n,int r)  {
    return n*n*n*a+n*n*l+n*n+r;
}

void init_fill_D() {
    for (int a=0;a<n;a++) {
        for (int h=a;h<n;h++) {
            for (int k=h;k<n;k++) {
                for (int l=k;l<n;l++) {
                    for (int p=l;p<n;p++) {
                        for (int r=p;r<n;r++) {
                            D[index_D(a,h,k,l,p,r)] = INT_MIN;
                        }
                    }
                }
            }
        }
    }
}

int index_D(int a,int h,int k,int l,int p,int r)  {
    return n*n*n*n*n*a+n*n*n*n*h+n*n*n*k+n*n*l+n*p+r;
}

void init_fill_E() {
    for (int b=0;b<n;b++) {
        for (int g=b;g<n;g++) {
            for (int k=g;k<n;k++) {
                for (int l=k;l<n;l++) {
                    for (int p=l;p<n;p++) {
                        for (int r=p;r<n;r++) {
                            E[index_E(b,g,k,l,p,r)] = INT_MIN;
                        }
                    }
                }
            }
        }
    }
}

int index_E(int b,int g,int k,int l,int p,int r)  {
    return n*n*n*n*n*b+n*n*n*n*g+n*n*n*k+n*n*l+n*p+r;
}

void init_fill_F() {
    for (int b=0;b<n;b++) {
        for (int e=b;e<n;e++) {
            for (int p=e;p<n;p++) {
                for (int r=p;r<n;r++) {
                    F[index_F(b,e,p,r)] = INT_MIN;
                }
            }
        }
    }
}

int index_F(int b,int e,int p,int r)  {
    return n*n*n*b+n*n*e+n*p+r;
}

void init_fill_G() {
    for (int c=0;c<n;c++) {
        for (int e=c;e<n;e++) {
            for (int p=e;p<n;p++) {
                for (int r=p;r<n;r++) {
                    G[index_G(c,e,p,r)] = INT_MIN;
                }
            }
        }
    }
}

int index_G(int c,int e,int p,int r)  {
    return n*n*n*c+n*n*e+n*p+r;
}

void init_fill_H() {
    for (int c=0;c<n;c++) {
        for (int e=c;e<n;e++) {
            for (int q=e;q<n;q++) {
                for (int r=q;r<n;r++) {
                    H[index_H(c,e,q,r)] = INT_MIN;
                }
            }
        }
    }
}

int index_H(int c,int e,int q,int r)  {
    return n*n*n*c+n*n*e+n*q+r;
}

void init_fill_I() {
    for (int b=0;b<n;b++) {
        for (int e=b;e<n;e++) {
            for (int g=e;g<n;g++) {
                for (int k=g;k<n;k++) {
                    for (int l=k;l<n;l++) {
                        I[index_I(b,e,g,k,l)] = INT_MIN;
                    }
                }
            }
        }
    }
}

int index_I(int b,int e,int g,int k,int l)  {
    return n*n*n*n*b+n*n*n*e+n*n*g+n*k+l;
}

void init_fill_J() {
    for (int a=0;a<n;a++) {
        for (int h=a;h<n;h++) {
            for (int k=h;k<n;k++) {
                for (int n=k;n<n;n++) {
                    for (int p=n;p<n;p++) {
                        J[index_J(a,h,k,n,p)] = INT_MIN;
                    }
                }
            }
        }
    }
}

int index_J(int a,int h,int k,int n,int p)  {
    return n*n*n*n*a+n*n*n*h+n*n*k+n*n+p;
}

void init_fill_K() {
    for (int i=0;i<n;i++) {
        for (int k=i;k<n;k++) {
            for (int n=k;n<n;n++) {
                for (int p=n;p<n;p++) {
                    K[index_K(i,k,n,p)] = INT_MIN;
                }
            }
        }
    }
}

int index_K(int i,int k,int n,int p)  {
    return n*n*n*i+n*n*k+n*n+p;
}

void init_fill_L() {
    for (int i=0;i<n;i++) {
        for (int k=i;k<n;k++) {
            for (int o=k;o<n;o++) {
                for (int p=o;p<n;p++) {
                    L[index_L(i,k,o,p)] = INT_MIN;
                }
            }
        }
    }
}

int index_L(int i,int k,int o,int p)  {
    return n*n*n*i+n*n*k+n*o+p;
}

void init_fill_M() {
    for (int m=0;m<n;m++) {
        for (int n=m;n<n;n++) {
            for (int r=n;r<n;r++) {
                for (int t=r;t<n;t++) {
                    M[index_M(m,n,r,t)] = INT_MIN;
                }
            }
        }
    }
}

int index_M(int m,int n,int r,int t)  {
    return n*n*n*m+n*n*n+n*r+t;
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
        for (int m=a;m<n;m++) {
            for (int n=m;n<n;n++) {
                for (int r=n;r<n;r++) {
                    for (int t=r;t<n;t++) {
                        min_value = min(min_value, compute_B(m,a,r,n)+compute_M(m,t,r,n));
                    }
                }
            }
        }
    }

    return min_value;
}

int compute_B(int a,int m,int n,int r) {
    if (B[index_B(a,m,n,r)] > INT_MIN) {
        return B[index_B(a,m,n,r)];
    }
    
    int min_value = INT_MAX;
    
    for (int l=a;l<r;l++) {
        min_value = min(min_value, compute_C(l,a,r,n));
    }

    B[index_B(a,m,n,r)] = min_value;
    return min_value;
}

int compute_C(int a,int l,int n,int r) {
    if (C[index_C(a,l,n,r)] > INT_MIN) {
        return C[index_C(a,l,n,r)];
    }
    
    int min_value = INT_MAX;
    
    for (int h=a;h<r;h++) {
        for (int k=h;k<r;k++) {
            for (int p=k;p<r;p++) {
                min_value = min(min_value, compute_D(a,h,k,l,p,r)+compute_J(n,h,p,k,a));
            }
        }
    }

    C[index_C(a,l,n,r)] = min_value;
    return min_value;
}

int compute_D(int a, int h, int p,int l,int k,int r) {
    if (D[index_D(a,h, p,l,k,r)] > INT_MIN) {
        return D[index_D(a,h, p,l,k,r)];
    }

    int min_value = INT_MAX;
    bool eq_some_const1 = (h-1!=p) && (h-1!=l) && (h-1!=k) && (h-1!=r);
    bool eq_some_const2 = (h-1!=p) && (h-1!=l) && (h-1!=k) && (h-1!=r);

    if (a+1 < h) {
        if (!eq_some_const1) {
            min_value = min(min_value, compute_D(a+1,h,p,l,k,r));
        }
    }
    if (h-1 > a) {
        if (!eq_some_const2) {
            min_value = min(min_value, compute_D2(a, h-1, p,l,k,r));
        }
    }
    if (!eq_some_const1 && !eq_some_const2) {
        min_value = min(min_value, compute_D(a+1, h-1, p,l,k,r)+bp_score(line[a], 
                                                                                      line[h]));
    }

    min_value = min(min_value, compute_E(r,p,1,k,l,20));

    D[index_D(a,h,p,l,k,r)] = min_value;
    return min_value;
} 

int compute_D2(int a, int h, int p,int l,int k,int r) {
    if (D2[index_D(a,h, p,l,k,r)] > INT_MIN) {
        return D2[index_D(a,h, p,l,k,r)];
    }

    int min_value = INT_MAX;
    bool eq_some_const1 = (h-1!=p) && (h-1!=l) && (h-1!=k) && (h-1!=r);
    bool eq_some_const2 = (h-1!=p) && (h-1!=l) && (h-1!=k) && (h-1!=r);

    if (h-1 > a && !eq_some_const2) {
        min_value = min(min_value, compute_D2(a, h-1, p,l,k,r));
    }
    if (!eq_some_const1 && !eq_some_const2) {
        min_value = min(min_value, compute_D(a+1, h-1, p,l,k,r)+bp_score(line[a], 
                                                                                      line[h]));
    }

    D2[index_D(a,h,p,l,k,r)] = min_value;
    return min_value;
}


int compute_E(int b,int g,int k,int l,int p,int r) {
    if (E[index_E(b,g,k,l,p,r)] > INT_MIN) {
        return E[index_E(b,g,k,l,p,r)];
    }
    
    int min_value = INT_MAX;
    
    for (int e=0;e<r;e++) {
        min_value = min(min_value, compute_F(p,b,e,r)+compute_I(b,e,k,l,g));
    }

    E[index_E(b,g,k,l,p,r)] = min_value;
    return min_value;
}

int compute_F(int b,int e,int p,int r) {
    if (F[index_F(b,e,p,r)] > INT_MIN) {
        return F[index_F(b,e,p,r)];
    }
    
    int min_value = INT_MAX;
    
    for (int c=r;c<n;c++) {
        min_value = min(min_value, compute_G(p,e,c,r));
    }

    F[index_F(b,e,p,r)] = min_value;
    return min_value;
}

int compute_G(int c,int e,int p,int r) {
    if (G[index_G(c,e,p,r)] > INT_MIN) {
        return G[index_G(c,e,p,r)];
    }
    
    int min_value = INT_MAX;
    
    for (int q=p;q<r;q++) {
        min_value = min(min_value, compute_H(q,e,c,r));
    }

    G[index_G(c,e,p,r)] = min_value;
    return min_value;
}

int compute_H(int c,int e,int q,int r) {
    if (H[index_H(c,e,q,r)] > INT_MIN) {
        return H[index_H(c,e,q,r)];
    }
    
    int min_value = INT_MAX;
    
    for (int d=0;d<r;d++) {
        min_value = min(min_value, compute_CLIQUE(c,d,q,r));
    }

    H[index_H(c,e,q,r)] = min_value;
    return min_value;
}


int compute_I(int b,int e,int g,int k,int l) {
    if (I[index_I(b,e,g,k,l)] > INT_MIN) {
        return I[index_I(b,e,g,k,l)];
    }
    
    int min_value = INT_MAX;
    
    for (int f=e;f<l;f++) {
        min_value = min(min_value, compute_CLIQUE(e,f,k,l));
    }

    I[index_I(b,e,g,k,l)] = min_value;
    return min_value;
}


int compute_J(int a,int h,int k,int n,int p) {
    if (J[index_J(a,h,k,n,p)] > INT_MIN) {
        return J[index_J(a,h,k,n,p)];
    }
    
    int min_value = INT_MAX;
    
    for (int i=h;i<p;i++) {
        min_value = min(min_value, compute_K(p,i,k,n));
    }

    J[index_J(a,h,k,n,p)] = min_value;
    return min_value;
}

int compute_K(int i,int k,int n,int p) {
    if (K[index_K(i,k,n,p)] > INT_MIN) {
        return K[index_K(i,k,n,p)];
    }
    
    int min_value = INT_MAX;
    
    for (int o=n;o<p;o++) {
        min_value = min(min_value, compute_L(p,o,i,k));
    }

    K[index_K(i,k,n,p)] = min_value;
    return min_value;
}

int compute_L(int i,int k,int o,int p) {
    if (L[index_L(i,k,o,p)] > INT_MIN) {
        return L[index_L(i,k,o,p)];
    }
    
    int min_value = INT_MAX;
    
    for (int j=i;j<p;j++) {
        min_value = min(min_value, compute_CLIQUE(i,j,o,p));
    }

    L[index_L(i,k,o,p)] = min_value;
    return min_value;
}


int compute_M(int m,int n,int r,int t) {
    if (M[index_M(m,n,r,t)] > INT_MIN) {
        return M[index_M(m,n,r,t)];
    }
    
    int min_value = INT_MAX;
    
    for (int s=r;s<t;s++) {
        min_value = min(min_value, compute_CLIQUE(m,n,s,t));
    }

    M[index_M(m,n,r,t)] = min_value;
    return min_value;
}


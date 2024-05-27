#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdbool.h>

int min(int a, int b) { if (a<b) {return a;} else {return b;}};

int index_B(int a,int i,int j,int r);
int index_C(int a,int h,int j,int r);
int index_D(int a,int g,int j,int q);
int index_E(int c,int g,int l,int q);
int index_F(int c,int d,int l,int n);
int index_G(int d,int g,int n,int q);
int index_H(int d,int g,int n,int p);
int index_I(int e,int p,int g,int n);
int index_J(int a,int l,int c,int j);
int index_K(int i,int j,int r,int t);
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
void init_fill_CLIQUE();

int compute_A();
int compute_B(int a,int i,int j,int r);
int compute_C(int a,int h,int j,int r);
int compute_D(int a,int g,int j,int q);
int compute_E(int c,int g,int l,int q);
int compute_F(int c,int d,int l,int n);
int compute_G(int d,int g,int n,int q);
int compute_H(int d,int g,int n,int p);
int compute_I(int e,int p,int g,int n);
int compute_I2(int e,int p,int g,int n);
int compute_J(int a,int l,int c,int j);
int compute_J2(int a,int l,int c,int j);
int compute_K(int i,int j,int r,int t);
int fold();

int n;

char * line = NULL;

double * B;
double * C;
double * D;
double * E;
double * F;
double * G;
double * H;
double * I;
double * I2;
double * J;
double * J2;
double * K;
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
    double * D = malloc(n*n*n*n*sizeof(double));
    double * E = malloc(n*n*n*n*sizeof(double));
    double * F = malloc(n*n*n*n*sizeof(double));
    double * G = malloc(n*n*n*n*sizeof(double));
    double * H = malloc(n*n*n*n*sizeof(double));
    double * I = malloc(n*n*n*n*sizeof(double));
double * I2 = malloc(n*n*n*n*sizeof(double));
    double * J = malloc(n*n*n*n*sizeof(double));
double * J2 = malloc(n*n*n*n*sizeof(double));
    double * K = malloc(n*n*n*n*sizeof(double));
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
        for (int i=a;i<n;i++) {
            for (int j=i;j<n;j++) {
                for (int r=j;r<n;r++) {
                    B[index_B(a,i,j,r)] = INT_MIN;
                }
            }
        }
    }
}

int index_B(int a,int i,int j,int r)  {
    return n*n*n*a+n*n*i+n*j+r;
}

void init_fill_C() {
    for (int a=0;a<n;a++) {
        for (int h=a;h<n;h++) {
            for (int j=h;j<n;j++) {
                for (int r=j;r<n;r++) {
                    C[index_C(a,h,j,r)] = INT_MIN;
                }
            }
        }
    }
}

int index_C(int a,int h,int j,int r)  {
    return n*n*n*a+n*n*h+n*j+r;
}

void init_fill_D() {
    for (int a=0;a<n;a++) {
        for (int g=a;g<n;g++) {
            for (int j=g;j<n;j++) {
                for (int q=j;q<n;q++) {
                    D[index_D(a,g,j,q)] = INT_MIN;
                }
            }
        }
    }
}

int index_D(int a,int g,int j,int q)  {
    return n*n*n*a+n*n*g+n*j+q;
}

void init_fill_E() {
    for (int c=0;c<n;c++) {
        for (int g=c;g<n;g++) {
            for (int l=g;l<n;l++) {
                for (int q=l;q<n;q++) {
                    E[index_E(c,g,l,q)] = INT_MIN;
                }
            }
        }
    }
}

int index_E(int c,int g,int l,int q)  {
    return n*n*n*c+n*n*g+n*l+q;
}

void init_fill_F() {
    for (int c=0;c<n;c++) {
        for (int d=c;d<n;d++) {
            for (int l=d;l<n;l++) {
                for (int n=l;n<n;n++) {
                    F[index_F(c,d,l,n)] = INT_MIN;
                }
            }
        }
    }
}

int index_F(int c,int d,int l,int n)  {
    return n*n*n*c+n*n*d+n*l+n;
}

void init_fill_G() {
    for (int d=0;d<n;d++) {
        for (int g=d;g<n;g++) {
            for (int n=g;n<n;n++) {
                for (int q=n;q<n;q++) {
                    G[index_G(d,g,n,q)] = INT_MIN;
                }
            }
        }
    }
}

int index_G(int d,int g,int n,int q)  {
    return n*n*n*d+n*n*g+n*n+q;
}

void init_fill_H() {
    for (int d=0;d<n;d++) {
        for (int g=d;g<n;g++) {
            for (int n=g;n<n;n++) {
                for (int p=n;p<n;p++) {
                    H[index_H(d,g,n,p)] = INT_MIN;
                }
            }
        }
    }
}

int index_H(int d,int g,int n,int p)  {
    return n*n*n*d+n*n*g+n*n+p;
}

void init_fill_I() {
    for (int e=0;e<n;e++) {
        for (int p=e;p<n;p++) {
            for (int g=p;g<n;g++) {
                for (int n=g;n<n;n++) {
                    I[index_I(e,p,g,n)] = INT_MIN;
                }
            }
        }
    }
}

int index_I(int e,int p,int g,int n)  {
    return n*n*n*e+n*n*p+n*g+n;
}

void init_fill_J() {
    for (int a=0;a<n;a++) {
        for (int l=a;l<n;l++) {
            for (int c=l;c<n;c++) {
                for (int j=c;j<n;j++) {
                    J[index_J(a,l,c,j)] = INT_MIN;
                }
            }
        }
    }
}

int index_J(int a,int l,int c,int j)  {
    return n*n*n*a+n*n*l+n*c+j;
}

void init_fill_K() {
    for (int i=0;i<n;i++) {
        for (int j=i;j<n;j++) {
            for (int r=j;r<n;r++) {
                for (int t=r;t<n;t++) {
                    K[index_K(i,j,r,t)] = INT_MIN;
                }
            }
        }
    }
}

int index_K(int i,int j,int r,int t)  {
    return n*n*n*i+n*n*j+n*r+t;
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
        for (int i=a;i<n;i++) {
            for (int j=i;j<n;j++) {
                for (int r=j;r<n;r++) {
                    for (int t=r;t<n;t++) {
                        min_value = min(min_value, compute_B(a,i,j,r)+compute_K(i,j,t,r));
                    }
                }
            }
        }
    }

    return min_value;
}

int compute_B(int a,int i,int j,int r) {
    if (B[index_B(a,i,j,r)] > INT_MIN) {
        return B[index_B(a,i,j,r)];
    }
    
    int min_value = INT_MAX;
    
    for (int h=a;h<r;h++) {
        min_value = min(min_value, compute_C(a,h,j,r));
    }

    B[index_B(a,i,j,r)] = min_value;
    return min_value;
}

int compute_C(int a,int h,int j,int r) {
    if (C[index_C(a,h,j,r)] > INT_MIN) {
        return C[index_C(a,h,j,r)];
    }
    
    int min_value = INT_MAX;
    
    for (int g=a;g<r;g++) {
        for (int q=g;q<r;q++) {
            min_value = min(min_value, compute_D(a,q,g,j)+compute_CLIQUE(g,h,q,r));
        }
    }

    C[index_C(a,h,j,r)] = min_value;
    return min_value;
}

int compute_D(int a,int g,int j,int q) {
    if (D[index_D(a,g,j,q)] > INT_MIN) {
        return D[index_D(a,g,j,q)];
    }
    
    int min_value = INT_MAX;
    
    for (int l=j;l<q;l++) {
        for (int c=l;c<n;c++) {
            min_value = min(min_value, compute_E(l,q,g,c)+compute_J(a,l,c,j));
        }
    }

    D[index_D(a,g,j,q)] = min_value;
    return min_value;
}

int compute_E(int c,int g,int l,int q) {
    if (E[index_E(c,g,l,q)] > INT_MIN) {
        return E[index_E(c,g,l,q)];
    }
    
    int min_value = INT_MAX;
    
    for (int d=0;d<q;d++) {
        for (int n=d;n<q;n++) {
            min_value = min(min_value, compute_F(l,d,n,c)+compute_G(q,g,n,d));
        }
    }

    E[index_E(c,g,l,q)] = min_value;
    return min_value;
}

int compute_F(int c,int d,int l,int n) {
    if (F[index_F(c,d,l,n)] > INT_MIN) {
        return F[index_F(c,d,l,n)];
    }
    
    int min_value = INT_MAX;
    
    for (int m=l;m<n;m++) {
        min_value = min(min_value, compute_CLIQUE(c,d,m,n));
    }

    F[index_F(c,d,l,n)] = min_value;
    return min_value;
}


int compute_G(int d,int g,int n,int q) {
    if (G[index_G(d,g,n,q)] > INT_MIN) {
        return G[index_G(d,g,n,q)];
    }
    
    int min_value = INT_MAX;
    
    for (int p=n;p<q;p++) {
        min_value = min(min_value, compute_H(g,p,n,d));
    }

    G[index_G(d,g,n,q)] = min_value;
    return min_value;
}

int compute_H(int d,int g,int n,int p) {
    if (H[index_H(d,g,n,p)] > INT_MIN) {
        return H[index_H(d,g,n,p)];
    }
    
    int min_value = INT_MAX;
    
    for (int e=d;e<p;e++) {
        min_value = min(min_value, compute_I(e,p,g,n));
    }

    H[index_H(d,g,n,p)] = min_value;
    return min_value;
}

int compute_I(int e, int p, int g,int n) {
    if (I[index_I(e,p, g,n)] > INT_MIN) {
        return I[index_I(e,p, g,n)];
    }

    int min_value = INT_MAX;
    bool eq_some_const1 = (p-1!=g) && (p-1!=n);
    bool eq_some_const2 = (p-1!=g) && (p-1!=n);

    if (e+1 < p) {
        if (!eq_some_const1) {
            min_value = min(min_value, compute_I(e+1,p,g,n));
        }
    }
    if (p-1 > e) {
        if (!eq_some_const2) {
            min_value = min(min_value, compute_I2(e, p-1, g,n));
        }
    }
    if (!eq_some_const1 && !eq_some_const2) {
        min_value = min(min_value, compute_I(e+1, p-1, g,n)+bp_score(line[e], 
                                                                                      line[p]));
    }

    min_value = min(min_value, INT_MAX // no children);

    I[index_I(e,p,g,n)] = min_value;
    return min_value;
} 

int compute_I2(int e, int p, int g,int n) {
    if (I2[index_I(e,p, g,n)] > INT_MIN) {
        return I2[index_I(e,p, g,n)];
    }

    int min_value = INT_MAX;
    bool eq_some_const1 = (p-1!=g) && (p-1!=n);
    bool eq_some_const2 = (p-1!=g) && (p-1!=n);

    if (p-1 > e && !eq_some_const2) {
        min_value = min(min_value, compute_I2(e, p-1, g,n));
    }
    if (!eq_some_const1 && !eq_some_const2) {
        min_value = min(min_value, compute_I(e+1, p-1, g,n)+bp_score(line[e], 
                                                                                      line[p]));
    }

    I2[index_I(e,p,g,n)] = min_value;
    return min_value;
}


int compute_J(int a, int l, int j,int c) {
    if (J[index_J(a,l, j,c)] > INT_MIN) {
        return J[index_J(a,l, j,c)];
    }

    int min_value = INT_MAX;
    bool eq_some_const1 = (l-1!=j) && (l-1!=c);
    bool eq_some_const2 = (l-1!=j) && (l-1!=c);

    if (a+1 < l) {
        if (!eq_some_const1) {
            min_value = min(min_value, compute_J(a+1,l,j,c));
        }
    }
    if (l-1 > a) {
        if (!eq_some_const2) {
            min_value = min(min_value, compute_J2(a, l-1, j,c));
        }
    }
    if (!eq_some_const1 && !eq_some_const2) {
        min_value = min(min_value, compute_J(a+1, l-1, j,c)+bp_score(line[a], 
                                                                                      line[l]));
    }

    min_value = min(min_value, INT_MAX // no children);

    J[index_J(a,l,j,c)] = min_value;
    return min_value;
} 

int compute_J2(int a, int l, int j,int c) {
    if (J2[index_J(a,l, j,c)] > INT_MIN) {
        return J2[index_J(a,l, j,c)];
    }

    int min_value = INT_MAX;
    bool eq_some_const1 = (l-1!=j) && (l-1!=c);
    bool eq_some_const2 = (l-1!=j) && (l-1!=c);

    if (l-1 > a && !eq_some_const2) {
        min_value = min(min_value, compute_J2(a, l-1, j,c));
    }
    if (!eq_some_const1 && !eq_some_const2) {
        min_value = min(min_value, compute_J(a+1, l-1, j,c)+bp_score(line[a], 
                                                                                      line[l]));
    }

    J2[index_J(a,l,j,c)] = min_value;
    return min_value;
}



int compute_K(int i,int j,int r,int t) {
    if (K[index_K(i,j,r,t)] > INT_MIN) {
        return K[index_K(i,j,r,t)];
    }
    
    int min_value = INT_MAX;
    
    for (int s=r;s<t;s++) {
        min_value = min(min_value, compute_CLIQUE(i,j,s,t));
    }

    K[index_K(i,j,r,t)] = min_value;
    return min_value;
}


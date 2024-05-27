#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <float.h>

//HASHING TABLE FUNCTION
#define LOAD_FACTOR_THRESHOLD 0.7
#define TUPLE_SIZE 6

int TABLE_SIZE = 64019;
int N = 10;
int THETA_PAIRING=1;
int THETA = 100;

void preset_tablesize(int size) {
    if (size > 0) {
        TABLE_SIZE = size;
    }
    else {
        printf("Size negative or null\n");
    }
}

void preset_N(int n) {
    if (n > 0) {
        N = n;
    }
    else {
        printf("N negative or null\n");
    }
}

typedef struct Tuple {
    int values[TUPLE_SIZE];
} Tuple;

Tuple *createTuple(int val[],int size) {
    Tuple *tuple = (Tuple *)malloc(sizeof(Tuple));

    // Initialize all values in tuple to 0
    for (int i = 0; i < TUPLE_SIZE; i++) {
        tuple->values[i] = 0;
    }

    // Copy values from val array, up to its length or TUPLE_SIZE
    for (int i = 0; i < size; i++) {
        tuple->values[i] = val[i];
    }

    return tuple;
}

int compare_tuple(Tuple *tp1, int val[], int size) {
    int i;
    // Compare elements until the end of val or tp1->values
    for (i = 0; i < size ; i++) {
        if (tp1->values[i] != val[i]) {
            return 0; // Tuples are different
        }
    }
    // Check if remaining elements of tp1->values are all zeros
    for (; i < TUPLE_SIZE; i++) {
        if (tp1->values[i] != 0) {
            return 0; // Tuples are different
        }
    }
    return 1; // Tuples are identical
}

// Define a structure for each node in the hash table
typedef struct HashNode {
    struct Tuple* key;
    float value;
    struct HashNode* next;
} HashNode;

// Define the hash table structure
typedef struct {
    HashNode **buckets;
    int size;
    int capacity;
} HashTable;

// Hash function
int hash(int keys[], HashTable *hashtable,int size) {
    int key = 0;
    for (int i = 0; i < size; i++) {
        key += keys[i] * pow(N, TUPLE_SIZE-i);
    }
    return ( hashtable->capacity + (key % hashtable->capacity) ) % hashtable->capacity;
}


// Function to initialize a new hash table
HashTable *createHashTable() {
    HashTable *hashTable = (HashTable *)malloc(sizeof(HashTable));
    hashTable->buckets = (HashNode **)calloc(TABLE_SIZE, sizeof(HashNode *));
    hashTable->size = 0;
    hashTable->capacity = TABLE_SIZE;
    return hashTable;
}

// Function to free memory used by the hash table 
void destroyHashTable(HashTable *hashTable) {
    for (int i = 0; i < hashTable->capacity; i++) {
        HashNode *current = hashTable->buckets[i];
        while (current != NULL) {
            HashNode *temp = current;
            free(temp->key);
            current = current->next;
            free(temp);
        }
    }
    free(hashTable->buckets);
    free(hashTable);
}

// Function to insert a key-value pair into the hash table
void insert(HashTable *hashTable, int keys[], int size,float value) {
    // Combine the numbers into a single key
    int index = hash(keys, hashTable,size);
    // Create a new node
    HashNode *newNode = (HashNode *)malloc(sizeof(HashNode));
    newNode->key = createTuple(keys,size);
    newNode->value = value;
    newNode->next = NULL;

    // Check for collision
    if (hashTable->buckets[index] == NULL) {
        hashTable->buckets[index] = newNode;
    }
    else {
        // Collision detected, add node to the end of the linked list
        HashNode *current = hashTable->buckets[index];
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = newNode;
    }

    hashTable->size++;

    // Check if resize is needed
    if ((double)hashTable->size / hashTable->capacity > LOAD_FACTOR_THRESHOLD) {
        // Resize the hash table
        //printf("%d resize \n", hashTable->capacity);
        int new_capacity = hashTable->capacity * 2;
        HashNode **new_buckets = (HashNode **)calloc(new_capacity, sizeof(HashNode *));
        for (int i = 0; i < hashTable->capacity; i++) {
            HashNode *current = hashTable->buckets[i];
            while (current != NULL) {
                HashNode *temp = current;
                current = current->next;
                int new_index = hash(temp->key->values, hashTable,size);
                temp->next = new_buckets[new_index];
                new_buckets[new_index] = temp;
            }
        }
        free(hashTable->buckets);
        hashTable->buckets = new_buckets;
        hashTable->capacity = new_capacity;
    }
}

// Function to retrieve a value from the hash table given a key
bool get(HashTable *hashTable, int keys[], int size, float *value) {
    
    int index = hash(keys, hashTable,size);
    
    // Traverse the linked list at the index
    HashNode *current = hashTable->buckets[index];
    while (current != NULL) {
        if (compare_tuple(current->key, keys,size)) {
            *value = current->value;
            return true; // Key found
        }
        current = current->next;
    }
    return false; // Key not found
}

void print_tuple(Tuple *tpl){
    printf("tp=");
    for (int i = 0; i < TUPLE_SIZE; i++) {
        printf("%d,",tpl->values[i]);
    }

}

void print_table(HashTable *hashTable){
    
    int index = 0;
    
    // Traverse the linked list at the index
     for (int i = 0; i < hashTable->capacity; i++) {
        HashNode *current = hashTable->buckets[i];
        while (current != NULL) {
            print_tuple(current->key);
            printf(": %f",current->value);
                
            printf("||");
            current = current->next;
        }
        printf("\n");
    }
}
//SECONDARY STRUCTURE FUNCTION
int **matrix = NULL; 

void createMatrix(int N,char * ss) {

    // Allocate memory for the matrix
    matrix = (int **)malloc(N * sizeof(int *));
    for (int i = 0; i < N; i++) {
        matrix[i] = (int *)malloc(N * sizeof(int));
    }

    // Fill the matrix
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if(i==j){
                if(ss[i]=='|'){
                    matrix[i][j]=0;
                }
                else{
                    matrix[i][j]=1;
                }
                continue;
            }
            matrix[i][j] = 0;
            matrix[j][i] = 0;
            if ( ss[i]=='x' || ss[j]=='x') {
                matrix[i][j] = 0;
                matrix[j][i] = 0;
                continue;
            }
            if( (ss[i]=='(' && ss[j]==')') || (ss[j]=='(' && ss[i]==')')){
                matrix[i][j]= 1;
                matrix[j][i] = 1;
                continue;
            }
            if( ss[i]=='>'){
                matrix[i][j]= 0;
                matrix[j][i]= 0;
                continue;
            }
            if( (ss[i]=='.' || ss[i]=='<') && (ss[j]=='.' || ss[j]=='>')){
                matrix[i][j]= 1;
                matrix[j][i]= 1;
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
                        matrix[i][j]=1;
                        matrix[j][i]=1;
                    }
                    else{
                        matrix[i][j]=0;
                        matrix[j][i]=0;
                    }
                }
                else{
                    matrix[i][j] = 0;
                    matrix[j][i] = 0;
                }
                continue;
            }
        }
    }
}

// Function to display the matrix
void displayMatrix(int N) {
    printf("Matrix:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

// Function to free the memory allocated for the matrix
void freeMatrix(int N) {
    for (int i = 0; i < N; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

//PAIRING FUNCTIONS
int min(int a, int b) { if (a<b) {return a;} else {return b;}};
int max(int a, int b) { if (a<b) {return b;} else {return a;}};

float min_float(float a, float b) { if (a<b) {return a;} else {return b;}};
float max_float(float a, float b) { if (a<b) {return b;} else {return a;}};

char * line = NULL;
int MAX;
char * correct_score = NULL;
char * structure=NULL;
char bracket='A';

float matching(char x, char y){
    if(x!='X' && x!='G' && x!='C' && x!='A' && x!='U' && x!='T'){
        printf("x is not a recognizable charactere\n");
        exit(-1);
    }
    if(y!='X' && y!='G' && y!='C' && y!='A' && y!='U' && y!='T'){
        printf("y is not a recognizable charactere\n");
        exit(-1);
    }
    //printf("%c,%c\n",x,y);
    if (x=='G' && y=='C') { return -10.0; }
    if (x=='C' && y=='G') { return -10.0; }
    if (x=='G' && y=='U') { return -5.0; }
    if (x=='U' && y=='G') { return -5.0; }
    if (x=='A' && y=='U') { return -5.0; }
    if (x=='U' && y=='A') { return -5.0; }
    return FLT_MAX;
}

int evaluate(int x, int y) {
    if(abs(x-y)<THETA_PAIRING || matrix[x][y]==0 || matching(line[x],line[y])==FLT_MAX){
        return 0;
    }
    return 1;
}

float bp_score(int x, int y) {
    if(abs(x-y)<THETA_PAIRING || matrix[x][y]==0){
        return FLT_MAX;
    }
    return matching(line[x],line[y]);
}

float add(float a,float b){
    if(a==FLT_MAX || b==FLT_MAX){
        return FLT_MAX;
    }

    return a+b;
}

float MFEFree(float a,float b){
    if(b<a){
        return 0.0;
    }
    return -(b-a)-1.0;
}


void backtrace_MFEFree(float score, int a,int b){
    for(int i=a;i<=b;i++){
        if(score==MFEFree(i,b)){
            for(int y=i;y<=b;y++){
                 structure[y]=bracket-17;
            }
            return;
        }
    }
}

float INTB(int a,int b,int c,int d){
    float bp=bp_score(a,b);
    float mfe=add(MFEFree(c,a-1),MFEFree(b+1,d));
    return add(bp,mfe);
}

void backtrace_INTB(float score,int a,int b,int c,int d){
    structure[a]=bracket;
    structure[b]=bracket+32;

    backtrace_MFEFree(MFEFree(c,a-1),c,a-1);
    backtrace_MFEFree(MFEFree(b+1,d),b+1,d);
    return;
}




//declarations
float compute_CLIQUE0(HashTable *hashTable, int i, int j, int k, int l);

float compute_CLIQUE1(HashTable *hashTable, int i, int j, int k, int l);

void backtrace_CLIQUE0(HashTable *hashTable, float score, int i,int j, int k, int l);

void backtrace_CLIQUE1(HashTable *hashTable, float score, int i,int j, int k, int l);

float compute_K(HashTable *hashTable,int i,int j,int r,int t) ;

void backtrace_K(HashTable *hashTable,float score,int i,int j,int r,int t) ;

float compute_J0(HashTable *hashTable,int a, int l, int c,int j);

float compute_J1(HashTable *hashTable,int a, int l, int c,int j);

void backtrace_J0(HashTable *hashTable, float score, int a, int l, int c,int j) ;

void backtrace_J1(HashTable *hashTable, float score, int a, int l, int c,int j) ;

float compute_I0(HashTable *hashTable,int e, int p, int g,int n);

float compute_I1(HashTable *hashTable,int e, int p, int g,int n);

void backtrace_I0(HashTable *hashTable, float score, int e, int p, int g,int n) ;

void backtrace_I1(HashTable *hashTable, float score, int e, int p, int g,int n) ;

float compute_H(HashTable *hashTable,int d,int g,int n,int p) ;

void backtrace_H(HashTable *hashTable,float score,int d,int g,int n,int p) ;

float compute_G(HashTable *hashTable,int d,int g,int n,int q) ;

void backtrace_G(HashTable *hashTable,float score,int d,int g,int n,int q) ;

float compute_F(HashTable *hashTable,int c,int d,int l,int n) ;

void backtrace_F(HashTable *hashTable,float score,int c,int d,int l,int n) ;

float compute_E(HashTable *hashTable,int c,int g,int l,int q) ;

void backtrace_E(HashTable *hashTable,float score,int c,int g,int l,int q) ;

float compute_D(HashTable *hashTable,int a,int g,int j,int q) ;

void backtrace_D(HashTable *hashTable,float score,int a,int g,int j,int q) ;

float compute_C(HashTable *hashTable,int a,int h,int j,int r) ;

void backtrace_C(HashTable *hashTable,float score,int a,int h,int j,int r) ;

float compute_B(HashTable *hashTable,int a,int i,int j,int r) ;

void backtrace_B(HashTable *hashTable,float score,int a,int i,int j,int r) ;

float compute_A(HashTable *hashTable,int a,int t) ;

void backtrace_A(HashTable *hashTable,float score,int a,int t) ;


float fold(HashTable * hashTable) {
    float min_value=FLT_MAX;
    for(int t=0+9;t<=MAX;t++){
        for(int a=0;a<t-9;a++){
            float mfe1=MFEFree(0,a-1);
            float mfe2=MFEFree(t,MAX-1);
            float mfe=add(mfe1,mfe2);
            min_value = min_float(min_value,add(compute_A(hashTable,a,t),mfe));
        }
    }
    return min_value;
}

void backtrace(HashTable * hashTable,float score) {
    for(int t=0+9;t<=MAX;t++){
        for(int a=0;a<t-9;a++){
            
            float mfe1=MFEFree(0,a-1);
            float mfe2=MFEFree(t,MAX-1);
            float mfe=add(mfe1,mfe2);
            float tmp0=compute_A(hashTable,a,t);
            //printf("h = %d, a=%d, %d+%d %d,%d\n",h,a,mfe2,mfe1,score,compute_A(hashTable,a,h));
            if(score==add(mfe,tmp0)){
                backtrace_A(hashTable,tmp0,a,t);
                backtrace_MFEFree(mfe1,0,a-1);
                backtrace_MFEFree(mfe2,t,MAX-1);
                return;
            }
        }
    }
}


int main(int argc, char ** argv) {
    printf("File name: %s\n",argv[1]);
    FILE * fp = fopen(argv[1], "r");
    if (fp == NULL) {
        perror("Error opening file");
        return EXIT_FAILURE;
    }
    int nb_tests=0;
    clock_t start_time, end_time;
    double total_time = 0.0;

    start_time = clock();
    while(true){
        HashTable *hashTable = createHashTable();
        int b_len=0;

        int len=0;

        //reading of the sequence
        len=getline(&line,&b_len, fp);
        
        if(len<=0){
            destroyHashTable(hashTable);
            free(line);
            free(correct_score);
            printf("End of file, %d tested",nb_tests);
            break;
        }
        if(line[0]=='#'){
            printf("%s",line);
            continue;
        }
        MAX=len;
        if(line[len-1]=='\n'){
            MAX=len-1;
            nb_tests+=1;
        }
        preset_N(MAX);

        //secondary structure
        char *ss = NULL;
        len=getline(&ss,&b_len, fp);
        createMatrix(MAX,ss);
        //reading the score
        len=getline(&correct_score,&b_len, fp);
        if (len != -1) {
            if (correct_score[len - 1] == '\n') {
                correct_score[len - 1] = '\0'; // Replace newline character with null terminator
                len--; // Decrement the length to exclude the removed '\n'
            }
            // Convert the line to an integer
            
        }
        else{
            destroyHashTable(hashTable);
            free(line);
            free(correct_score);
            free(ss);
            printf("No integer found for test %d\n",nb_tests);
            exit(-1);
        }
        float number = (float)atoi(correct_score);
        printf("Test: %d Size of the Sequence: %d ---", nb_tests,MAX);
        
        //elements used to record backtrack
        bracket='A';
        structure=malloc(sizeof(char)*(MAX+1));
        structure[MAX]='\0';
        for(int i=0;i<MAX;i++){
            structure[i]='.';
        }
        //Start of test
        clock_t test_start_time = clock();

        //folding
        float score = fold(hashTable);
        //retrieve backtrack
        backtrace(hashTable,score);
        
        clock_t test_end_time = clock();

        double test_time = ((double)(test_end_time - test_start_time)) / CLOCKS_PER_SEC;
        total_time += test_time;
        //end of test
        
        //sanity checks
        if(score == number){
            printf("Correct");
            int nb=0;
            for(int i=0;i<MAX;i++){
                if(structure[i]=='.'){
                    printf("missing bt\n");
                    break;
                }
                if(structure[i]>64){
                    if(structure[i]<91){
                        
                        nb+=structure[i];
                    }
                    else{
                        
                        nb-=(structure[i]-32);
                    }
                }
            }
            printf("(backtrack impurty= %d)\n",nb);
        }
        else{
            destroyHashTable(hashTable);
            free(line);
            free(correct_score);
            free(structure);
            free(ss);
            printf("Failed (%d found, should be %d)\n",score,number);
            exit(nb_tests);
        }

        
        printf("%s\n",structure);
        
        //clean up
        destroyHashTable(hashTable);
        freeMatrix(MAX);
        free(structure);
        free(line);
        free(ss);
    }
    end_time = clock();
    double total_execution_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    printf(" Total execution time: %.4f seconds", total_execution_time);
    printf("(mean %.4f seconds)\n", total_time / nb_tests);

    fclose(fp);
    return 0;
}

float compute_CLIQUE0(HashTable *hashTable,int i, int i2, int j2, int j){
    int CLIQUE=-1;
    float value;
    int tab[] = {CLIQUE,i,i2,j2,j};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        return FLT_MAX;
    }

    if (get(hashTable,tab,size,&value)) { 
        return value;
    }
    
    float min_value = add(INTB(i,j,i,j),compute_CLIQUE1(hashTable,i+1,i2,j2,j-1));

    insert(hashTable,tab,size,min_value);
    return min_value;
}

float compute_CLIQUE1(HashTable *hashTable,int i, int i2, int j2, int j) {
    int CLIQUE=1;
    float value;
    int tab[] = {CLIQUE,i,i2,j2,j};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        return FLT_MAX;
    }

    if (get(hashTable,tab,size,&value)) { 
        return value;
    }
    
    float min_value = 0.0;
    float tmp;
    if (j == j2 && i <= i2) { 
        min_value = min_float(min_value, MFEFree(i,i2)); 
    }
    for(int k=i;k<=i2;k++){
        for(int l=max(j2,k+1);l<=j;l++){
            if(evaluate(k,l) && l-k+1<=THETA && (k<i2 || l==j2)){
                tmp=INTB(k,l,i,j);
                if(tmp!=INT_MAX){
                    min_value = min_float(min_value,add(compute_CLIQUE1(hashTable,k+1,i2,j2,l-1),tmp));
                }
            }
        }
    }
    
    insert(hashTable,tab,size,min_value);
    return min_value;
}

float compute_K(HashTable *hashTable,int i,int j,int r,int t) {
    int K = 75;
    float value;
    int tab[]={K,i,j,r,t};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value = FLT_MAX;
    
    for (int s=r;s<t;s++) {
      if(!evaluate(i,t-1)||!evaluate(j-1,s)){continue;}
      float mfe0 = MFEFree(r,s-1);

      float tmp0= compute_CLIQUE0( hashTable,i,j-1,s,t-1);
      if(tmp0==FLT_MAX){continue;}

      min_value = min_float(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

float compute_J0(HashTable *hashTable,int a, int l, int c,int j){
    int J0 = 74;
    float value;

    int tab[]={J0,a,l,c,j};
    int size = 5;

    if(a<0 || l<0 || l>=MAX || a>=MAX || a>l){
        return FLT_MAX;
    }
    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value=add(INTB(a,l,a,l),compute_J1(hashTable,a+1,l-1,c,j));

    insert(hashTable,tab,size,min_value);
    return min_value;
}
float compute_J1(HashTable *hashTable,int a, int l, int c,int j) {
    int J1 = -74;
    float value;

    int tab[]={J1,a,l,c,j};
    int size = 5;

    if(a<0 || l<0 || l>=MAX || a>=MAX || a>l){
        return FLT_MAX;
    }
    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value=FLT_MAX;
    float mfe1=MFEFree(a,c-1);
    float mfe2=MFEFree(j,l);
    
    min_value=add(0,add(mfe1,mfe2));

    loop:
    for(int tmp1=a;tmp1<=min(l-1,c-1);tmp1++){
        for(int tmp2=max(j,tmp1+5);tmp2<=l;tmp2++){
            if(evaluate(tmp1,tmp2) && tmp1-tmp2+1<=THETA){
                float tmp3=INTB(tmp1,tmp2,a,l);

                if(tmp3!=INT_MAX){
                    min_value=min_float(min_value,
                    add(compute_J1(hashTable,tmp1+1,tmp2-1,c,j),tmp3));
                }
            }
        }    
    }
    insert(hashTable,tab,size,min_value);
    return min_value;
}
float compute_I0(HashTable *hashTable,int e, int p, int g,int n){
    int I0 = 73;
    float value;

    int tab[]={I0,e,p,g,n};
    int size = 5;

    if(e<0 || p<0 || p>=MAX || e>=MAX || e>p){
        return FLT_MAX;
    }
    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value=add(INTB(e,p,e,p),compute_I1(hashTable,e+1,p-1,g,n));

    insert(hashTable,tab,size,min_value);
    return min_value;
}
float compute_I1(HashTable *hashTable,int e, int p, int g,int n) {
    int I1 = -73;
    float value;

    int tab[]={I1,e,p,g,n};
    int size = 5;

    if(e<0 || p<0 || p>=MAX || e>=MAX || e>p){
        return FLT_MAX;
    }
    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value=FLT_MAX;
    float mfe1=MFEFree(e,g-1);
    float mfe2=MFEFree(n,p);
    
    min_value=add(0,add(mfe1,mfe2));

    loop:
    for(int tmp1=e;tmp1<=min(p-1,g-1);tmp1++){
        for(int tmp2=max(n,tmp1+5);tmp2<=p;tmp2++){
            if(evaluate(tmp1,tmp2) && tmp1-tmp2+1<=THETA){
                float tmp3=INTB(tmp1,tmp2,e,p);

                if(tmp3!=INT_MAX){
                    min_value=min_float(min_value,
                    add(compute_I1(hashTable,tmp1+1,tmp2-1,g,n),tmp3));
                }
            }
        }    
    }
    insert(hashTable,tab,size,min_value);
    return min_value;
}
float compute_H(HashTable *hashTable,int d,int g,int n,int p) {
    int H = 72;
    float value;
    int tab[]={H,d,g,n,p};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value = FLT_MAX;
    
    for (int e=d;e<g;e++) {
      if(!evaluate(e,p-1)){continue;}
      float mfe0 = MFEFree(d,e-1);

      float tmp0= compute_I0( hashTable,e,p-1,g,n);
      if(tmp0==FLT_MAX){continue;}

      min_value = min_float(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

float compute_G(HashTable *hashTable,int d,int g,int n,int q) {
    int G = 71;
    float value;
    int tab[]={G,d,g,n,q};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value = FLT_MAX;
    
    for (int p=n+1;p<q+1;p++) {

      float mfe0 = MFEFree(p,q-1);

      float tmp0= compute_H( hashTable,d,g,n,p);
      if(tmp0==FLT_MAX){continue;}

      min_value = min_float(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

float compute_F(HashTable *hashTable,int c,int d,int l,int n) {
    int F = 70;
    float value;
    int tab[]={F,c,d,l,n};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value = FLT_MAX;
    
    for (int m=l;m<n;m++) {
      if(!evaluate(c,n-1)||!evaluate(d-1,m)){continue;}
      float mfe0 = MFEFree(l,m-1);

      float tmp0= compute_CLIQUE0( hashTable,c,d-1,m,n-1);
      if(tmp0==FLT_MAX){continue;}

      min_value = min_float(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

float compute_E(HashTable *hashTable,int c,int g,int l,int q) {
    int E = 69;
    float value;
    int tab[]={E,c,g,l,q};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value = FLT_MAX;
    
    for (int d=c+1;d<g;d++) {
        for (int n=l+1;n<q;n++) {


          float tmp0= compute_G( hashTable,d,g,n,q);
          if(tmp0==FLT_MAX){continue;}
          float tmp1= compute_F( hashTable,c,d,l,n);
          if(tmp1==FLT_MAX){continue;}

          min_value = min_float(min_value,add(add(tmp0,tmp1),0));
        }
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

float compute_D(HashTable *hashTable,int a,int g,int j,int q) {
    int D = 68;
    float value;
    int tab[]={D,a,g,j,q};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value = FLT_MAX;
    
    for (int c=a+1;c<g-1;c++) {
        for (int l=j+1;l<q-1;l++) {
          if(!evaluate(a,l-1)){continue;}

          float tmp0= compute_E( hashTable,c,g,l,q);
          if(tmp0==FLT_MAX){continue;}
          float tmp1= compute_J0( hashTable,a,l-1,c,j);
          if(tmp1==FLT_MAX){continue;}

          min_value = min_float(min_value,add(add(tmp0,tmp1),0));
        }
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

float compute_C(HashTable *hashTable,int a,int h,int j,int r) {
    int C = 67;
    float value;
    int tab[]={C,a,h,j,r};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value = FLT_MAX;
    
    for (int g=a+3;g<h;g++) {
        for (int q=j+3;q<r;q++) {
          if(!evaluate(g,r-1)||!evaluate(h-1,q)){continue;}

          float tmp0= compute_D( hashTable,a,g,j,q);
          if(tmp0==FLT_MAX){continue;}
          float tmp1= compute_CLIQUE0( hashTable,g,h-1,q,r-1);
          if(tmp1==FLT_MAX){continue;}

          min_value = min_float(min_value,add(add(tmp0,tmp1),0));
        }
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

float compute_B(HashTable *hashTable,int a,int i,int j,int r) {
    int B = 66;
    float value;
    int tab[]={B,a,i,j,r};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value = FLT_MAX;
    
    for (int h=a+4;h<i+1;h++) {

      float mfe0 = MFEFree(h,i-1);

      float tmp0= compute_C( hashTable,a,h,j,r);
      if(tmp0==FLT_MAX){continue;}

      min_value = min_float(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

float compute_A(HashTable *hashTable,int a,int t) {
    int A = 65;
    float value;
    int tab[]={A,a,t};
    int size = 3;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value = FLT_MAX;
    
    for (int i=a+4;i<t-5;i++) {
        for (int j=i+1;j<t-4;j++) {
            for (int r=j+4;r<t;r++) {


              float tmp0= compute_B( hashTable,a,i,j,r);
              if(tmp0==FLT_MAX){continue;}
              float tmp1= compute_K( hashTable,i,j,r,t);
              if(tmp1==FLT_MAX){continue;}

              min_value = min_float(min_value,add(add(tmp0,tmp1),0));
            }
        }
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}
void backtrace_CLIQUE0(HashTable *hashTable,float score,int i, int i2, int j2, int j) {
    
   if(i2>=j2 || i2<i || j2>j){
        return ;
    }
    float tmp=compute_CLIQUE1(hashTable,i+1,i2,j2,j-1);
    float sc=INTB(i,j,i,j);
    if(score==add(sc,tmp)){
        backtrace_INTB(sc,i,j,i,j);
        backtrace_CLIQUE1(hashTable,tmp,i+1,i2,j2,j-1);
    }

    return;
}

void backtrace_CLIQUE1(HashTable *hashTable,float score,int i, int i2, int j2, int j) {
    
   if(i2>=j2 || i2<i || j2>j){
        return ;
    }
    float tmp;
    if (j == j2 && i <= i2 && score==MFEFree(i,i2)) { 
        backtrace_MFEFree(score,i,i2);
    }
    for(int k=i;k<=i2;k++){
        for(int l=max(j2,k+1);l<=j;l++){
            if(evaluate(k,l) && l-k+1<=THETA && (k<i2 || l==j2)){
                tmp=INTB(k,l,i,j);
                if(tmp!=FLT_MAX){
                    float tmp0=compute_CLIQUE1(hashTable,k+1,i2,j2,l-1);
                    if(score==add(tmp0,tmp)){
                        backtrace_INTB(tmp,k,l,i,j); 
                        backtrace_CLIQUE1(hashTable,tmp0,k+1,i2,j2,l-1);                       
                        return;
                    }
                }
            }
        }
    }
    return;
}

void backtrace_K(HashTable *hashTable,float score,int i,int j,int r,int t) {

    for (int s=r;s<t;s++) {
      if(!evaluate(i,t-1)||!evaluate(j-1,s)){continue;}
      float mfe0 = MFEFree(r,s-1);

      float tmp0= compute_CLIQUE0( hashTable,i,j-1,s,t-1);
      if(tmp0==FLT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,r,s-1);

        backtrace_CLIQUE0(hashTable,tmp0, i,j-1,s,t-1);
        bracket+=1;

        return;
      }
    }


}

void backtrace_J0(HashTable *hashTable, float score, int a, int l, int c,int j) {
    
    if(a<0 || l<0 || l>=MAX || a>=MAX || a>l){
        return;
    }
    float tmp=compute_J1(hashTable,a+1,l-1,c,j);
    float sc=INTB(a,l,a,l);
    if(score==add(sc,tmp)){
        backtrace_INTB(sc,a,l,a,l);
        backtrace_J1(hashTable,tmp,a+1,l-1,c,j);
    }
    return;
}
void backtrace_J1(HashTable *hashTable, float score, int a, int l, int c,int j) {
    
    if(a<0 || l<0 || l>=MAX || a>=MAX || a>l){
        return;
    }

    float mfe1=MFEFree(a,c-1);
    float mfe2=MFEFree(j,l);
    

    if(score==add(0,add(mfe1,mfe2))){
        
        backtrace_MFEFree(mfe1,a,c-1);
        backtrace_MFEFree(mfe2,j,l);

        return;
    }
    
    loop:
    for(int tmp1=a;tmp1<=min(l-1,c-1);tmp1++){
        for(int tmp2=max(j,tmp1+5);tmp2<=l;tmp2++){
            if(evaluate(tmp1,tmp2) && tmp1-tmp2+1<=THETA){
                float tmp3=INTB(tmp1,tmp2,a,l);
                if(tmp3!=INT_MAX){
                    float tmp0=compute_J1(hashTable,tmp1+1,tmp2-1,c,j);
                    if(score==add(tmp0,tmp3)){
                        backtrace_INTB(tmp3,tmp1,tmp2,a,l);
                        backtrace_J1(hashTable,tmp0,tmp1+1,tmp2-1,c,j);
                        return;
                    }
                }
            }
        }    
    }
    return;
} 
void backtrace_I0(HashTable *hashTable, float score, int e, int p, int g,int n) {
    
    if(e<0 || p<0 || p>=MAX || e>=MAX || e>p){
        return;
    }
    float tmp=compute_I1(hashTable,e+1,p-1,g,n);
    float sc=INTB(e,p,e,p);
    if(score==add(sc,tmp)){
        backtrace_INTB(sc,e,p,e,p);
        backtrace_I1(hashTable,tmp,e+1,p-1,g,n);
    }
    return;
}
void backtrace_I1(HashTable *hashTable, float score, int e, int p, int g,int n) {
    
    if(e<0 || p<0 || p>=MAX || e>=MAX || e>p){
        return;
    }

    float mfe1=MFEFree(e,g-1);
    float mfe2=MFEFree(n,p);
    

    if(score==add(0,add(mfe1,mfe2))){
        
        backtrace_MFEFree(mfe1,e,g-1);
        backtrace_MFEFree(mfe2,n,p);

        return;
    }
    
    loop:
    for(int tmp1=e;tmp1<=min(p-1,g-1);tmp1++){
        for(int tmp2=max(n,tmp1+5);tmp2<=p;tmp2++){
            if(evaluate(tmp1,tmp2) && tmp1-tmp2+1<=THETA){
                float tmp3=INTB(tmp1,tmp2,e,p);
                if(tmp3!=INT_MAX){
                    float tmp0=compute_I1(hashTable,tmp1+1,tmp2-1,g,n);
                    if(score==add(tmp0,tmp3)){
                        backtrace_INTB(tmp3,tmp1,tmp2,e,p);
                        backtrace_I1(hashTable,tmp0,tmp1+1,tmp2-1,g,n);
                        return;
                    }
                }
            }
        }    
    }
    return;
} 
void backtrace_H(HashTable *hashTable,float score,int d,int g,int n,int p) {

    for (int e=d;e<g;e++) {
      if(!evaluate(e,p-1)){continue;}
      float mfe0 = MFEFree(d,e-1);

      float tmp0= compute_I0( hashTable,e,p-1,g,n);
      if(tmp0==FLT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,d,e-1);

        backtrace_I0(hashTable,tmp0, e,p-1,g,n);
        bracket+=1;

        return;
      }
    }


}

void backtrace_G(HashTable *hashTable,float score,int d,int g,int n,int q) {

    for (int p=n+1;p<q+1;p++) {

      float mfe0 = MFEFree(p,q-1);

      float tmp0= compute_H( hashTable,d,g,n,p);
      if(tmp0==FLT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,p,q-1);

        backtrace_H(hashTable,tmp0, d,g,n,p);

        return;
      }
    }


}

void backtrace_F(HashTable *hashTable,float score,int c,int d,int l,int n) {

    for (int m=l;m<n;m++) {
      if(!evaluate(c,n-1)||!evaluate(d-1,m)){continue;}
      float mfe0 = MFEFree(l,m-1);

      float tmp0= compute_CLIQUE0( hashTable,c,d-1,m,n-1);
      if(tmp0==FLT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,l,m-1);

        backtrace_CLIQUE0(hashTable,tmp0, c,d-1,m,n-1);
        bracket+=1;

        return;
      }
    }


}

void backtrace_E(HashTable *hashTable,float score,int c,int g,int l,int q) {

    for (int d=c+1;d<g;d++) {
        for (int n=l+1;n<q;n++) {


          float tmp0= compute_G( hashTable,d,g,n,q);
          if(tmp0==FLT_MAX){continue;}
          float tmp1= compute_F( hashTable,c,d,l,n);
          if(tmp1==FLT_MAX){continue;}

          if(score==add(add(tmp0,tmp1),0)){

            backtrace_G(hashTable,tmp0, d,g,n,q);
            backtrace_F(hashTable,tmp1, c,d,l,n);

            return;
          }
        }
    }


}

void backtrace_D(HashTable *hashTable,float score,int a,int g,int j,int q) {

    for (int c=a+1;c<g-1;c++) {
        for (int l=j+1;l<q-1;l++) {
          if(!evaluate(a,l-1)){continue;}

          float tmp0= compute_E( hashTable,c,g,l,q);
          if(tmp0==FLT_MAX){continue;}
          float tmp1= compute_J0( hashTable,a,l-1,c,j);
          if(tmp1==FLT_MAX){continue;}

          if(score==add(add(tmp0,tmp1),0)){

            backtrace_E(hashTable,tmp0, c,g,l,q);
            backtrace_J0(hashTable,tmp1, a,l-1,c,j);
            bracket+=1;

            return;
          }
        }
    }


}

void backtrace_C(HashTable *hashTable,float score,int a,int h,int j,int r) {

    for (int g=a+3;g<h;g++) {
        for (int q=j+3;q<r;q++) {
          if(!evaluate(g,r-1)||!evaluate(h-1,q)){continue;}

          float tmp0= compute_D( hashTable,a,g,j,q);
          if(tmp0==FLT_MAX){continue;}
          float tmp1= compute_CLIQUE0( hashTable,g,h-1,q,r-1);
          if(tmp1==FLT_MAX){continue;}

          if(score==add(add(tmp0,tmp1),0)){

            backtrace_D(hashTable,tmp0, a,g,j,q);
            backtrace_CLIQUE0(hashTable,tmp1, g,h-1,q,r-1);
            bracket+=1;

            return;
          }
        }
    }


}

void backtrace_B(HashTable *hashTable,float score,int a,int i,int j,int r) {

    for (int h=a+4;h<i+1;h++) {

      float mfe0 = MFEFree(h,i-1);

      float tmp0= compute_C( hashTable,a,h,j,r);
      if(tmp0==FLT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,h,i-1);

        backtrace_C(hashTable,tmp0, a,h,j,r);

        return;
      }
    }


}

void backtrace_A(HashTable *hashTable,float score,int a,int t) {

    for (int i=a+4;i<t-5;i++) {
        for (int j=i+1;j<t-4;j++) {
            for (int r=j+4;r<t;r++) {


              float tmp0= compute_B( hashTable,a,i,j,r);
              if(tmp0==FLT_MAX){continue;}
              float tmp1= compute_K( hashTable,i,j,r,t);
              if(tmp1==FLT_MAX){continue;}

              if(score==add(add(tmp0,tmp1),0)){

                backtrace_B(hashTable,tmp0, a,i,j,r);
                backtrace_K(hashTable,tmp1, i,j,r,t);

                return;
              }
            }
        }
    }


}

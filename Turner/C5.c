#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <time.h>

//HASHING TABLE FUNCTION
#define LOAD_FACTOR_THRESHOLD 0.7
#define TUPLE_SIZE 8

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
    int value;
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
void insert(HashTable *hashTable, int keys[], int size,int value) {
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
bool get(HashTable *hashTable, int keys[], int size, int *value) {
    
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
            printf(": %d",current->value);
                
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

char * line = NULL;
int MAX;
char * correct_score = NULL;
char * structure=NULL;
char bracket='A';

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
    return INT_MAX;
}

int evaluate(int x, int y) {
    if(abs(x-y)<THETA_PAIRING || matrix[x][y]==0 || matching(line[x],line[y])==INT_MAX){
        return 0;
    }
    return 1;
}

int bp_score(int x, int y) {
    if(abs(x-y)<THETA_PAIRING || matrix[x][y]==0){
        return INT_MAX;
    }
    return matching(line[x],line[y]);
}

int add(int a,int b){
    if(a==INT_MAX || b==INT_MAX){
        return INT_MAX;
    }

    return a+b;
}

int MFEFree(int a,int b){
    if(b<a){
        return 0;
    }
    return -(b-a)-1;
}


void backtrace_MFEFree(int score, int a,int b){
    for(int i=a;i<=b;i++){
        if(score==MFEFree(i,b)){
            for(int y=i;y<=b;y++){
                 structure[y]=bracket-17;
            }
            return;
        }
    }
}

int INTB(int a,int b,int c,int d){
    int bp=bp_score(a,b);
    int mfe=add(MFEFree(c,a-1),MFEFree(b+1,d));
    return add(bp,mfe);
}

void backtrace_INTB(int score,int a,int b,int c,int d){
    structure[a]=bracket;
    structure[b]=bracket+32;

    backtrace_MFEFree(MFEFree(c,a-1),c,a-1);
    backtrace_MFEFree(MFEFree(b+1,d),b+1,d);
    return;
}




//declarations
int compute_CLIQUE0(HashTable *hashTable, int i, int j, int k, int l);

int compute_CLIQUE1(HashTable *hashTable, int i, int j, int k, int l);

void backtrace_CLIQUE0(HashTable *hashTable, int score, int i,int j, int k, int l);

void backtrace_CLIQUE1(HashTable *hashTable, int score, int i,int j, int k, int l);

int compute_M(HashTable *hashTable,int m,int n,int r,int t) ;

void backtrace_M(HashTable *hashTable,int score,int m,int n,int r,int t) ;

int compute_L(HashTable *hashTable,int i,int k,int o,int p) ;

void backtrace_L(HashTable *hashTable,int score,int i,int k,int o,int p) ;

int compute_K(HashTable *hashTable,int i,int k,int n,int p) ;

void backtrace_K(HashTable *hashTable,int score,int i,int k,int n,int p) ;

int compute_J(HashTable *hashTable,int a,int h,int k,int n,int p) ;

void backtrace_J(HashTable *hashTable,int score,int a,int h,int k,int n,int p) ;

int compute_I(HashTable *hashTable,int b,int e,int g,int k,int l) ;

void backtrace_I(HashTable *hashTable,int score,int b,int e,int g,int k,int l) ;

int compute_H(HashTable *hashTable,int c,int e,int q,int r) ;

void backtrace_H(HashTable *hashTable,int score,int c,int e,int q,int r) ;

int compute_G(HashTable *hashTable,int c,int e,int p,int r) ;

void backtrace_G(HashTable *hashTable,int score,int c,int e,int p,int r) ;

int compute_F(HashTable *hashTable,int b,int e,int p,int r) ;

void backtrace_F(HashTable *hashTable,int score,int b,int e,int p,int r) ;

int compute_E(HashTable *hashTable,int b,int g,int k,int l,int p,int r) ;

void backtrace_E(HashTable *hashTable,int score,int b,int g,int k,int l,int p,int r) ;

int compute_D0(HashTable *hashTable,int a, int h, int k,int l,int p,int r);

int compute_D1(HashTable *hashTable,int a, int h, int k,int l,int p,int r);

void backtrace_D0(HashTable *hashTable, int score, int a, int h, int k,int l,int p,int r) ;

void backtrace_D1(HashTable *hashTable, int score, int a, int h, int k,int l,int p,int r) ;

int compute_C(HashTable *hashTable,int a,int l,int n,int r) ;

void backtrace_C(HashTable *hashTable,int score,int a,int l,int n,int r) ;

int compute_B(HashTable *hashTable,int a,int m,int n,int r) ;

void backtrace_B(HashTable *hashTable,int score,int a,int m,int n,int r) ;

int compute_A(HashTable *hashTable,int a,int t) ;

void backtrace_A(HashTable *hashTable,int score,int a,int t) ;


int fold(HashTable * hashTable) {
    int min_value=INT_MAX;
    for(int t=0+9;t<=MAX;t++){
        for(int a=0;a<t-9;a++){
            int mfe1=MFEFree(0,a-1);
            int mfe2=MFEFree(t,MAX-1);
            int mfe=add(mfe1,mfe2);
            min_value = min(min_value,add(compute_A(hashTable,a,t),mfe));
        }
    }
    return min_value;
}

void backtrace(HashTable * hashTable,int score) {
    for(int t=0+9;t<=MAX;t++){
        for(int a=0;a<t-9;a++){
            
            int mfe1=MFEFree(0,a-1);
            int mfe2=MFEFree(t,MAX-1);
            int mfe=add(mfe1,mfe2);
            int tmp0=compute_A(hashTable,a,t);
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
        int number = atoi(correct_score);
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
        int score = fold(hashTable);
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

int compute_CLIQUE0(HashTable *hashTable,int i, int i2, int j2, int j){
    int CLIQUE=-1;
    int value;
    int tab[] = {CLIQUE,i,i2,j2,j};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        return INT_MAX;
    }

    if (get(hashTable,tab,size,&value)) { 
        return value;
    }
    
    int min_value = add(INTB(i,j,i,j),compute_CLIQUE1(hashTable,i+1,i2,j2,j-1));

    insert(hashTable,tab,size,min_value);
    return min_value;
}

int compute_CLIQUE1(HashTable *hashTable,int i, int i2, int j2, int j) {
    int CLIQUE=1;
    int value;
    int tab[] = {CLIQUE,i,i2,j2,j};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        return INT_MAX;
    }

    if (get(hashTable,tab,size,&value)) { 
        return value;
    }
    
    int min_value = 0;
    int tmp;
    if (j == j2 && i <= i2) { 
        min_value = min(min_value, MFEFree(i,i2)); 
    }
    for(int k=i;k<=i2;k++){
        for(int l=max(j2,k+1);l<=j;l++){
            if(evaluate(k,l) && l-k+1<=THETA && (k<i2 || l==j2)){
                tmp=INTB(k,l,i,j);
                if(tmp!=INT_MAX){
                    min_value = min(min_value,add(compute_CLIQUE1(hashTable,k+1,i2,j2,l-1),tmp));
                }
            }
        }
    }
    
    insert(hashTable,tab,size,min_value);
    return min_value;
}

int compute_M(HashTable *hashTable,int m,int n,int r,int t) {
    int M = 77;
    int value;
    int tab[]={M,m,n,r,t};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int s=r;s<t;s++) {
      if(!evaluate(m,t-1)||!evaluate(n-1,s)){continue;}
      int mfe0 = MFEFree(r,s-1);

      int tmp0= compute_CLIQUE0( hashTable,m,n-1,s,t-1);
      if(tmp0==INT_MAX){continue;}

      min_value = min(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

int compute_L(HashTable *hashTable,int i,int k,int o,int p) {
    int L = 76;
    int value;
    int tab[]={L,i,k,o,p};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int j=i+1;j<k+1;j++) {
      if(!evaluate(i,p-1)||!evaluate(j-1,o)){continue;}
      int mfe0 = MFEFree(j,k-1);

      int tmp0= compute_CLIQUE0( hashTable,i,j-1,o,p-1);
      if(tmp0==INT_MAX){continue;}

      min_value = min(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

int compute_K(HashTable *hashTable,int i,int k,int n,int p) {
    int K = 75;
    int value;
    int tab[]={K,i,k,n,p};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int o=n;o<p;o++) {

      int mfe0 = MFEFree(n,o-1);

      int tmp0= compute_L( hashTable,i,k,o,p);
      if(tmp0==INT_MAX){continue;}

      min_value = min(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

int compute_J(HashTable *hashTable,int a,int h,int k,int n,int p) {
    int J = 74;
    int value;
    int tab[]={J,a,h,k,n,p};
    int size = 6;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int i=h;i<k;i++) {

      int mfe0 = MFEFree(h,i-1);

      int tmp0= compute_K( hashTable,i,k,n,p);
      if(tmp0==INT_MAX){continue;}

      min_value = min(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

int compute_I(HashTable *hashTable,int b,int e,int g,int k,int l) {
    int I = 73;
    int value;
    int tab[]={I,b,e,g,k,l};
    int size = 6;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int f=e+1;f<g+1;f++) {
      if(!evaluate(e,l-1)||!evaluate(f-1,k)){continue;}
      int mfe0 = MFEFree(f,g-1);

      int tmp0= compute_CLIQUE0( hashTable,e,f-1,k,l-1);
      if(tmp0==INT_MAX){continue;}

      min_value = min(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

int compute_H(HashTable *hashTable,int c,int e,int q,int r) {
    int H = 72;
    int value;
    int tab[]={H,c,e,q,r};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int d=c+1;d<e+1;d++) {
      if(!evaluate(c,r-1)||!evaluate(d-1,q)){continue;}
      int mfe0 = MFEFree(d,e-1);

      int tmp0= compute_CLIQUE0( hashTable,c,d-1,q,r-1);
      if(tmp0==INT_MAX){continue;}

      min_value = min(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

int compute_G(HashTable *hashTable,int c,int e,int p,int r) {
    int G = 71;
    int value;
    int tab[]={G,c,e,p,r};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int q=p;q<r;q++) {

      int mfe0 = MFEFree(p,q-1);

      int tmp0= compute_H( hashTable,c,e,q,r);
      if(tmp0==INT_MAX){continue;}

      min_value = min(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

int compute_F(HashTable *hashTable,int b,int e,int p,int r) {
    int F = 70;
    int value;
    int tab[]={F,b,e,p,r};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int c=b;c<e;c++) {

      int mfe0 = MFEFree(b,c-1);

      int tmp0= compute_G( hashTable,c,e,p,r);
      if(tmp0==INT_MAX){continue;}

      min_value = min(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

int compute_E(HashTable *hashTable,int b,int g,int k,int l,int p,int r) {
    int E = 69;
    int value;
    int tab[]={E,b,g,k,l,p,r};
    int size = 7;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int e=b+1;e<g;e++) {


      int tmp0= compute_F( hashTable,b,e,p,r);
      if(tmp0==INT_MAX){continue;}
      int tmp1= compute_I( hashTable,b,e,g,k,l);
      if(tmp1==INT_MAX){continue;}

      min_value = min(min_value,add(add(tmp0,tmp1),0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

int compute_D0(HashTable *hashTable,int a, int h, int k,int l,int p,int r){
    int D0 = 68;
    int value;

    int tab[]={D0,a,h,k,l,p,r};
    int size = 7;

    if(a<0 || h<0 || h>=MAX || a>=MAX || a>h){
        return INT_MAX;
    }
    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value=add(INTB(a,h,a,h),compute_D1(hashTable,a+1,h-1,k,l,p,r));

    insert(hashTable,tab,size,min_value);
    return min_value;
}
int compute_D1(HashTable *hashTable,int a, int h, int k,int l,int p,int r) {
    int D1 = -68;
    int value;

    int tab[]={D1,a,h,k,l,p,r};
    int size = 7;

    if(a<0 || h<0 || h>=MAX || a>=MAX || a>h){
        return INT_MAX;
    }
    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value=INT_MAX;
    
    int tmp0= compute_E( hashTable,a,h+1,k,l,p,r);
    if(tmp0==INT_MAX){ goto loop;}

    min_value=add(tmp0,0);

    loop:
    for(int tmp1=a;tmp1<=min(h-1,INT_MAX);tmp1++){
        for(int tmp2=max(INT_MIN,tmp1+3);tmp2<=h;tmp2++){
            if(evaluate(tmp1,tmp2) && tmp1-tmp2+1<=THETA){
                int tmp3=INTB(tmp1,tmp2,a,h);

                if(tmp3!=INT_MAX){
                    min_value=min(min_value,
                    add(compute_D1(hashTable,tmp1+1,tmp2-1,k,l,p,r),tmp3));
                }
            }
        }    
    }
    insert(hashTable,tab,size,min_value);
    return min_value;
}
int compute_C(HashTable *hashTable,int a,int l,int n,int r) {
    int C = 67;
    int value;
    int tab[]={C,a,l,n,r};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int h=a+4;h<l-1;h++) {
        for (int k=h+1;k<l;k++) {
            for (int p=n+1;p<r;p++) {
              if(!evaluate(a,h-1)){continue;}

              int tmp0= compute_D0( hashTable,a,h-1,k,l,p,r);
              if(tmp0==INT_MAX){continue;}
              int tmp1= compute_J( hashTable,a,h,k,n,p);
              if(tmp1==INT_MAX){continue;}

              min_value = min(min_value,add(add(tmp0,tmp1),0));
            }
        }
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

int compute_B(HashTable *hashTable,int a,int m,int n,int r) {
    int B = 66;
    int value;
    int tab[]={B,a,m,n,r};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int l=a+6;l<m+1;l++) {

      int mfe0 = MFEFree(l,m-1);

      int tmp0= compute_C( hashTable,a,l,n,r);
      if(tmp0==INT_MAX){continue;}

      min_value = min(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

int compute_A(HashTable *hashTable,int a,int t) {
    int A = 65;
    int value;
    int tab[]={A,a,t};
    int size = 3;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int m=a+6;m<t-3;m++) {
        for (int n=m+1;n<t-2;n++) {
            for (int r=n+2;r<t;r++) {


              int tmp0= compute_B( hashTable,a,m,n,r);
              if(tmp0==INT_MAX){continue;}
              int tmp1= compute_M( hashTable,m,n,r,t);
              if(tmp1==INT_MAX){continue;}

              min_value = min(min_value,add(add(tmp0,tmp1),0));
            }
        }
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}
void backtrace_CLIQUE0(HashTable *hashTable,int score,int i, int i2, int j2, int j) {
    
   if(i2>=j2 || i2<i || j2>j){
        return ;
    }
    int tmp=compute_CLIQUE1(hashTable,i+1,i2,j2,j-1);
    int sc=INTB(i,j,i,j);
    if(score==add(sc,tmp)){
        backtrace_INTB(sc,i,j,i,j);
        backtrace_CLIQUE1(hashTable,tmp,i+1,i2,j2,j-1);
    }

    return;
}

void backtrace_CLIQUE1(HashTable *hashTable,int score,int i, int i2, int j2, int j) {
    
   if(i2>=j2 || i2<i || j2>j){
        return ;
    }
    int tmp;
    if (j == j2 && i <= i2 && score==MFEFree(i,i2)) { 
        backtrace_MFEFree(score,i,i2);
    }
    for(int k=i;k<=i2;k++){
        for(int l=max(j2,k+1);l<=j;l++){
            if(evaluate(k,l) && l-k+1<=THETA && (k<i2 || l==j2)){
                tmp=INTB(k,l,i,j);
                if(tmp!=INT_MAX){
                    int tmp0=compute_CLIQUE1(hashTable,k+1,i2,j2,l-1);
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

void backtrace_M(HashTable *hashTable,int score,int m,int n,int r,int t) {

    for (int s=r;s<t;s++) {
      if(!evaluate(m,t-1)||!evaluate(n-1,s)){continue;}
      int mfe0 = MFEFree(r,s-1);

      int tmp0= compute_CLIQUE0( hashTable,m,n-1,s,t-1);
      if(tmp0==INT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,r,s-1);

        backtrace_CLIQUE0(hashTable,tmp0, m,n-1,s,t-1);
        bracket+=1;

        return;
      }
    }


}

void backtrace_L(HashTable *hashTable,int score,int i,int k,int o,int p) {

    for (int j=i+1;j<k+1;j++) {
      if(!evaluate(i,p-1)||!evaluate(j-1,o)){continue;}
      int mfe0 = MFEFree(j,k-1);

      int tmp0= compute_CLIQUE0( hashTable,i,j-1,o,p-1);
      if(tmp0==INT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,j,k-1);

        backtrace_CLIQUE0(hashTable,tmp0, i,j-1,o,p-1);
        bracket+=1;

        return;
      }
    }


}

void backtrace_K(HashTable *hashTable,int score,int i,int k,int n,int p) {

    for (int o=n;o<p;o++) {

      int mfe0 = MFEFree(n,o-1);

      int tmp0= compute_L( hashTable,i,k,o,p);
      if(tmp0==INT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,n,o-1);

        backtrace_L(hashTable,tmp0, i,k,o,p);

        return;
      }
    }


}

void backtrace_J(HashTable *hashTable,int score,int a,int h,int k,int n,int p) {

    for (int i=h;i<k;i++) {

      int mfe0 = MFEFree(h,i-1);

      int tmp0= compute_K( hashTable,i,k,n,p);
      if(tmp0==INT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,h,i-1);

        backtrace_K(hashTable,tmp0, i,k,n,p);

        return;
      }
    }


}

void backtrace_I(HashTable *hashTable,int score,int b,int e,int g,int k,int l) {

    for (int f=e+1;f<g+1;f++) {
      if(!evaluate(e,l-1)||!evaluate(f-1,k)){continue;}
      int mfe0 = MFEFree(f,g-1);

      int tmp0= compute_CLIQUE0( hashTable,e,f-1,k,l-1);
      if(tmp0==INT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,f,g-1);

        backtrace_CLIQUE0(hashTable,tmp0, e,f-1,k,l-1);
        bracket+=1;

        return;
      }
    }


}

void backtrace_H(HashTable *hashTable,int score,int c,int e,int q,int r) {

    for (int d=c+1;d<e+1;d++) {
      if(!evaluate(c,r-1)||!evaluate(d-1,q)){continue;}
      int mfe0 = MFEFree(d,e-1);

      int tmp0= compute_CLIQUE0( hashTable,c,d-1,q,r-1);
      if(tmp0==INT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,d,e-1);

        backtrace_CLIQUE0(hashTable,tmp0, c,d-1,q,r-1);
        bracket+=1;

        return;
      }
    }


}

void backtrace_G(HashTable *hashTable,int score,int c,int e,int p,int r) {

    for (int q=p;q<r;q++) {

      int mfe0 = MFEFree(p,q-1);

      int tmp0= compute_H( hashTable,c,e,q,r);
      if(tmp0==INT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,p,q-1);

        backtrace_H(hashTable,tmp0, c,e,q,r);

        return;
      }
    }


}

void backtrace_F(HashTable *hashTable,int score,int b,int e,int p,int r) {

    for (int c=b;c<e;c++) {

      int mfe0 = MFEFree(b,c-1);

      int tmp0= compute_G( hashTable,c,e,p,r);
      if(tmp0==INT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,b,c-1);

        backtrace_G(hashTable,tmp0, c,e,p,r);

        return;
      }
    }


}

void backtrace_E(HashTable *hashTable,int score,int b,int g,int k,int l,int p,int r) {

    for (int e=b+1;e<g;e++) {


      int tmp0= compute_F( hashTable,b,e,p,r);
      if(tmp0==INT_MAX){continue;}
      int tmp1= compute_I( hashTable,b,e,g,k,l);
      if(tmp1==INT_MAX){continue;}

      if(score==add(add(tmp0,tmp1),0)){

        backtrace_F(hashTable,tmp0, b,e,p,r);
        backtrace_I(hashTable,tmp1, b,e,g,k,l);

        return;
      }
    }


}

void backtrace_D0(HashTable *hashTable, int score, int a, int h, int k,int l,int p,int r) {
    
    if(a<0 || h<0 || h>=MAX || a>=MAX || a>h){
        return;
    }
    int tmp=compute_D1(hashTable,a+1,h-1,k,l,p,r);
    int sc=INTB(a,h,a,h);
    if(score==add(sc,tmp)){
        backtrace_INTB(sc,a,h,a,h);
        backtrace_D1(hashTable,tmp,a+1,h-1,k,l,p,r);
    }
    return;
}
void backtrace_D1(HashTable *hashTable, int score, int a, int h, int k,int l,int p,int r) {
    
    if(a<0 || h<0 || h>=MAX || a>=MAX || a>h){
        return;
    }

    
    int tmp0= compute_E( hashTable,a,h+1,k,l,p,r);
    if(tmp0==INT_MAX){ goto loop;}


    if(score==add(tmp0,0)){
               bracket+=1;
       backtrace_E(hashTable,tmp0, a,h+1,k,l,p,r);

        
        return;
    }
    
    loop:
    for(int tmp1=a;tmp1<=min(h-1,INT_MAX);tmp1++){
        for(int tmp2=max(INT_MIN,tmp1+3);tmp2<=h;tmp2++){
            if(evaluate(tmp1,tmp2) && tmp1-tmp2+1<=THETA){
                int tmp3=INTB(tmp1,tmp2,a,h);
                if(tmp3!=INT_MAX){
                    int tmp0=compute_D1(hashTable,tmp1+1,tmp2-1,k,l,p,r);
                    if(score==add(tmp0,tmp3)){
                        backtrace_INTB(tmp3,tmp1,tmp2,a,h);
                        backtrace_D1(hashTable,tmp0,tmp1+1,tmp2-1,k,l,p,r);
                        return;
                    }
                }
            }
        }    
    }
    return;
} 
void backtrace_C(HashTable *hashTable,int score,int a,int l,int n,int r) {

    for (int h=a+4;h<l-1;h++) {
        for (int k=h+1;k<l;k++) {
            for (int p=n+1;p<r;p++) {
              if(!evaluate(a,h-1)){continue;}

              int tmp0= compute_D0( hashTable,a,h-1,k,l,p,r);
              if(tmp0==INT_MAX){continue;}
              int tmp1= compute_J( hashTable,a,h,k,n,p);
              if(tmp1==INT_MAX){continue;}

              if(score==add(add(tmp0,tmp1),0)){

                backtrace_D0(hashTable,tmp0, a,h-1,k,l,p,r);
                backtrace_J(hashTable,tmp1, a,h,k,n,p);

                return;
              }
            }
        }
    }


}

void backtrace_B(HashTable *hashTable,int score,int a,int m,int n,int r) {

    for (int l=a+6;l<m+1;l++) {

      int mfe0 = MFEFree(l,m-1);

      int tmp0= compute_C( hashTable,a,l,n,r);
      if(tmp0==INT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,l,m-1);

        backtrace_C(hashTable,tmp0, a,l,n,r);

        return;
      }
    }


}

void backtrace_A(HashTable *hashTable,int score,int a,int t) {

    for (int m=a+6;m<t-3;m++) {
        for (int n=m+1;n<t-2;n++) {
            for (int r=n+2;r<t;r++) {


              int tmp0= compute_B( hashTable,a,m,n,r);
              if(tmp0==INT_MAX){continue;}
              int tmp1= compute_M( hashTable,m,n,r,t);
              if(tmp1==INT_MAX){continue;}

              if(score==add(add(tmp0,tmp1),0)){

                backtrace_B(hashTable,tmp0, a,m,n,r);
                backtrace_M(hashTable,tmp1, m,n,r,t);

                return;
              }
            }
        }
    }


}
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

float compute_G(HashTable *hashTable,int f,int h,int k,int l) ;

void backtrace_G(HashTable *hashTable,float score,int f,int h,int k,int l) ;

float compute_F(HashTable *hashTable,int c,int e,int i,int j) ;

void backtrace_F(HashTable *hashTable,float score,int c,int e,int i,int j) ;

float compute_E(HashTable *hashTable,int b,int e,int i,int j) ;

void backtrace_E(HashTable *hashTable,float score,int b,int e,int i,int j) ;

float compute_D0(HashTable *hashTable,int a, int f, int i,int j);

float compute_D1(HashTable *hashTable,int a, int f, int i,int j);

void backtrace_D0(HashTable *hashTable, float score, int a, int f, int i,int j) ;

void backtrace_D1(HashTable *hashTable, float score, int a, int f, int i,int j) ;

float compute_C(HashTable *hashTable,int a,int f,int i,int k) ;

void backtrace_C(HashTable *hashTable,float score,int a,int f,int i,int k) ;

float compute_B(HashTable *hashTable,int a,int f,int h,int k) ;

void backtrace_B(HashTable *hashTable,float score,int a,int f,int h,int k) ;

float compute_A(HashTable *hashTable,int a,int l) ;

void backtrace_A(HashTable *hashTable,float score,int a,int l) ;


float fold(HashTable * hashTable) {
    float min_value=FLT_MAX;
    for(int l=0+5;l<=MAX;l++){
        for(int a=0;a<l-5;a++){
            float mfe1=MFEFree(0,a-1);
            float mfe2=MFEFree(l,MAX-1);
            float mfe=add(mfe1,mfe2);
            min_value = min_float(min_value,add(compute_A(hashTable,a,l),mfe));
        }
    }
    return min_value;
}

void backtrace(HashTable * hashTable,float score) {
    for(int l=0+5;l<=MAX;l++){
        for(int a=0;a<l-5;a++){
            
            float mfe1=MFEFree(0,a-1);
            float mfe2=MFEFree(l,MAX-1);
            float mfe=add(mfe1,mfe2);
            float tmp0=compute_A(hashTable,a,l);
            //printf("h = %d, a=%d, %d+%d %d,%d\n",h,a,mfe2,mfe1,score,compute_A(hashTable,a,h));
            if(score==add(mfe,tmp0)){
                backtrace_A(hashTable,tmp0,a,l);
                backtrace_MFEFree(mfe1,0,a-1);
                backtrace_MFEFree(mfe2,l,MAX-1);
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

float compute_G(HashTable *hashTable,int f,int h,int k,int l) {
    int G = 71;
    float value;
    int tab[]={G,f,h,k,l};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value = FLT_MAX;
    
    for (int g=f;g<h;g++) {
      if(!evaluate(g,l-1)||!evaluate(h-1,k)){continue;}
      float mfe0 = MFEFree(f,g-1);

      float tmp0= compute_CLIQUE0( hashTable,g,h-1,k,l-1);
      if(tmp0==FLT_MAX){continue;}

      min_value = min_float(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

float compute_F(HashTable *hashTable,int c,int e,int i,int j) {
    int F = 70;
    float value;
    int tab[]={F,c,e,i,j};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value = FLT_MAX;
    
    for (int d=c+1;d<e+1;d++) {
      if(!evaluate(c,j-1)||!evaluate(d-1,i)){continue;}
      float mfe0 = MFEFree(d,e-1);

      float tmp0= compute_CLIQUE0( hashTable,c,d-1,i,j-1);
      if(tmp0==FLT_MAX){continue;}

      min_value = min_float(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

float compute_E(HashTable *hashTable,int b,int e,int i,int j) {
    int E = 69;
    float value;
    int tab[]={E,b,e,i,j};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value = FLT_MAX;
    
    for (int c=b;c<e;c++) {

      float mfe0 = MFEFree(b,c-1);

      float tmp0= compute_F( hashTable,c,e,i,j);
      if(tmp0==FLT_MAX){continue;}

      min_value = min_float(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

float compute_D0(HashTable *hashTable,int a, int f, int i,int j){
    int D0 = 68;
    float value;

    int tab[]={D0,a,f,i,j};
    int size = 5;

    if(a<0 || f<0 || f>=MAX || a>=MAX || a>f){
        return FLT_MAX;
    }
    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value=add(INTB(a,f,a,f),compute_D1(hashTable,a+1,f-1,i,j));

    insert(hashTable,tab,size,min_value);
    return min_value;
}
float compute_D1(HashTable *hashTable,int a, int f, int i,int j) {
    int D1 = -68;
    float value;

    int tab[]={D1,a,f,i,j};
    int size = 5;

    if(a<0 || f<0 || f>=MAX || a>=MAX || a>f){
        return FLT_MAX;
    }
    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value=FLT_MAX;
    
    float tmp0= compute_E( hashTable,a,f+1,i,j);
    if(tmp0==FLT_MAX){ goto loop;}

    min_value=add(tmp0,0);

    loop:
    for(int tmp1=a;tmp1<=min(f-1,FLT_MAX);tmp1++){
        for(int tmp2=max(FLT_MIN,tmp1+2);tmp2<=f;tmp2++){
            if(evaluate(tmp1,tmp2) && tmp1-tmp2+1<=THETA){
                float tmp3=INTB(tmp1,tmp2,a,f);

                if(tmp3!=INT_MAX){
                    min_value=min_float(min_value,
                    add(compute_D1(hashTable,tmp1+1,tmp2-1,i,j),tmp3));
                }
            }
        }    
    }
    insert(hashTable,tab,size,min_value);
    return min_value;
}
float compute_C(HashTable *hashTable,int a,int f,int i,int k) {
    int C = 67;
    float value;
    int tab[]={C,a,f,i,k};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value = FLT_MAX;
    
    for (int j=i+1;j<k+1;j++) {
      if(!evaluate(a,f-1)){continue;}
      float mfe0 = MFEFree(j,k-1);

      float tmp0= compute_D0( hashTable,a,f-1,i,j);
      if(tmp0==FLT_MAX){continue;}

      min_value = min_float(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

float compute_B(HashTable *hashTable,int a,int f,int h,int k) {
    int B = 66;
    float value;
    int tab[]={B,a,f,h,k};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value = FLT_MAX;
    
    for (int i=h;i<k;i++) {

      float mfe0 = MFEFree(h,i-1);

      float tmp0= compute_C( hashTable,a,f,i,k);
      if(tmp0==FLT_MAX){continue;}

      min_value = min_float(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

float compute_A(HashTable *hashTable,int a,int l) {
    int A = 65;
    float value;
    int tab[]={A,a,l};
    int size = 3;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    float min_value = FLT_MAX;
    
    for (int f=a+3;f<l-2;f++) {
        for (int h=f+1;h<l-1;h++) {
            for (int k=h+1;k<l;k++) {


              float tmp0= compute_B( hashTable,a,f,h,k);
              if(tmp0==FLT_MAX){continue;}
              float tmp1= compute_G( hashTable,f,h,k,l);
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

void backtrace_G(HashTable *hashTable,float score,int f,int h,int k,int l) {

    for (int g=f;g<h;g++) {
      if(!evaluate(g,l-1)||!evaluate(h-1,k)){continue;}
      float mfe0 = MFEFree(f,g-1);

      float tmp0= compute_CLIQUE0( hashTable,g,h-1,k,l-1);
      if(tmp0==FLT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,f,g-1);

        backtrace_CLIQUE0(hashTable,tmp0, g,h-1,k,l-1);
        bracket+=1;

        return;
      }
    }


}

void backtrace_F(HashTable *hashTable,float score,int c,int e,int i,int j) {

    for (int d=c+1;d<e+1;d++) {
      if(!evaluate(c,j-1)||!evaluate(d-1,i)){continue;}
      float mfe0 = MFEFree(d,e-1);

      float tmp0= compute_CLIQUE0( hashTable,c,d-1,i,j-1);
      if(tmp0==FLT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,d,e-1);

        backtrace_CLIQUE0(hashTable,tmp0, c,d-1,i,j-1);
        bracket+=1;

        return;
      }
    }


}

void backtrace_E(HashTable *hashTable,float score,int b,int e,int i,int j) {

    for (int c=b;c<e;c++) {

      float mfe0 = MFEFree(b,c-1);

      float tmp0= compute_F( hashTable,c,e,i,j);
      if(tmp0==FLT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,b,c-1);

        backtrace_F(hashTable,tmp0, c,e,i,j);

        return;
      }
    }


}

void backtrace_D0(HashTable *hashTable, float score, int a, int f, int i,int j) {
    
    if(a<0 || f<0 || f>=MAX || a>=MAX || a>f){
        return;
    }
    float tmp=compute_D1(hashTable,a+1,f-1,i,j);
    float sc=INTB(a,f,a,f);
    if(score==add(sc,tmp)){
        backtrace_INTB(sc,a,f,a,f);
        backtrace_D1(hashTable,tmp,a+1,f-1,i,j);
    }
    return;
}
void backtrace_D1(HashTable *hashTable, float score, int a, int f, int i,int j) {
    
    if(a<0 || f<0 || f>=MAX || a>=MAX || a>f){
        return;
    }

    
    float tmp0= compute_E( hashTable,a,f+1,i,j);
    if(tmp0==FLT_MAX){ goto loop;}


    if(score==add(tmp0,0)){
               bracket+=1;
       backtrace_E(hashTable,tmp0, a,f+1,i,j);

        
        return;
    }
    
    loop:
    for(int tmp1=a;tmp1<=min(f-1,FLT_MAX);tmp1++){
        for(int tmp2=max(FLT_MIN,tmp1+2);tmp2<=f;tmp2++){
            if(evaluate(tmp1,tmp2) && tmp1-tmp2+1<=THETA){
                float tmp3=INTB(tmp1,tmp2,a,f);
                if(tmp3!=INT_MAX){
                    float tmp0=compute_D1(hashTable,tmp1+1,tmp2-1,i,j);
                    if(score==add(tmp0,tmp3)){
                        backtrace_INTB(tmp3,tmp1,tmp2,a,f);
                        backtrace_D1(hashTable,tmp0,tmp1+1,tmp2-1,i,j);
                        return;
                    }
                }
            }
        }    
    }
    return;
} 
void backtrace_C(HashTable *hashTable,float score,int a,int f,int i,int k) {

    for (int j=i+1;j<k+1;j++) {
      if(!evaluate(a,f-1)){continue;}
      float mfe0 = MFEFree(j,k-1);

      float tmp0= compute_D0( hashTable,a,f-1,i,j);
      if(tmp0==FLT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,j,k-1);

        backtrace_D0(hashTable,tmp0, a,f-1,i,j);

        return;
      }
    }


}

void backtrace_B(HashTable *hashTable,float score,int a,int f,int h,int k) {

    for (int i=h;i<k;i++) {

      float mfe0 = MFEFree(h,i-1);

      float tmp0= compute_C( hashTable,a,f,i,k);
      if(tmp0==FLT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,h,i-1);

        backtrace_C(hashTable,tmp0, a,f,i,k);

        return;
      }
    }


}

void backtrace_A(HashTable *hashTable,float score,int a,int l) {

    for (int f=a+3;f<l-2;f++) {
        for (int h=f+1;h<l-1;h++) {
            for (int k=h+1;k<l;k++) {


              float tmp0= compute_B( hashTable,a,f,h,k);
              if(tmp0==FLT_MAX){continue;}
              float tmp1= compute_G( hashTable,f,h,k,l);
              if(tmp1==FLT_MAX){continue;}

              if(score==add(add(tmp0,tmp1),0)){

                backtrace_B(hashTable,tmp0, a,f,h,k);
                backtrace_G(hashTable,tmp1, f,h,k,l);

                return;
              }
            }
        }
    }


}

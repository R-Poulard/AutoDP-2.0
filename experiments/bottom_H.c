#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>

#define LOAD_FACTOR_THRESHOLD 0.7
#define TUPLE_SIZE 5

int TABLE_SIZE = 64019;
int N = 10;
int THETA=1;

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
//
int min(int a, int b) { if (a<b) {return a;} else {return b;}};

char * line = NULL;
int MAX;
char * correct_score = NULL;
int * structure=NULL;
int index;

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
    if (x=='G' && y=='U') { return -1; }
    if (x=='U' && y=='G') { return -1; }
    if (x=='A' && y=='U') { return -5; }
    if (x=='U' && y=='A') { return -5; }
    return INT_MAX;
}

int bp_score(int x, int y) {
    if(abs(x-y)<THETA){
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


//declarations
int compute_CLIQUE(HashTable *hashTable, int i, int j, int k, int l);

int compute_CLIQUE2(HashTable *hashTable, int i, int j, int k, int l);

void backtrace_CLIQUE(HashTable *hashTable, int score, int i,int j, int k, int l);

void backtrace_CLIQUE2(HashTable *hashTable,int score,int i,int j,int k,int l);

void compute_C2(HashTable *hashTable,int a, int f, int c,int d) ;

void compute_C(HashTable *hashTable,int a, int f, int c,int d) ;

void backtrace_C2(HashTable *hashTable,int score,int a, int f, int c,int d) ;

void backtrace_C(HashTable *hashTable,int score,int a, int f, int c,int d) ;

void compute_B(HashTable *hashTable,int c,int d,int f,int h) ;

int backtrace_B(HashTable *hashTable,int score,int c,int d,int f,int h) ;

int compute_A(HashTable *hashTable,int START,int h) ;

int backtrace_A(HashTable *hashTable,int score,int START,int h) ;



void compute_CLIQUE_aux1(HashTable *hashTable, int c,int d, int g, int h){

}
void compute_C2_aux(HashTable *hashTable,int a,int f,int c,int d){
    //borns has to be check out (+1 >= ect.. + which one)
    //is it good? not sure...
    for(int a2=c+1;a2>=a;a--){
        for(int f2=d;f2<f;f++){
            compute_C(hashTable,a2,f2,c,d);
            compute_C2(hashTable,a2,f2,c,d);
        }
    }
}

int fold(HashTable * hashTable) {
    int START=0;
    int h=MAX;
    
    for (int a=START;a<h-3;a++) {
        for (int c=a+1;c<h-2;c++) {
            for (int d=c+1;d<h-1;d++) {
                for (int f=d+1;f<h;f++) {
                    compute_C2_aux(hashTable,a,f-1,c,d);
                    for (int g=f;g<h;g++) {
                        compute_CLIQUE_aux1(hashTable,c,d-1,g,h-1);
                    }
                    compute_B(hashTable,c,d,f,h);
                }
            }
        }
    }
    return compute_A(hashTable,0,MAX);
}

int backtrace(HashTable * hashTable,int score) {
    backtrace_A(hashTable,score,0,MAX);
}


int main(int argc, char ** argv) {
    printf("File name: %s\n",argv[1]);
    FILE * fp = fopen(argv[1], "r");
    int nb_tests=0;
    while(true){
        HashTable *hashTable = createHashTable();
        size_t b_len=0;

        int len=0;
        if (fp == NULL)
            exit(EXIT_FAILURE);
        len=getline(&line,&b_len, fp);
        
        if(len<=0){
            printf("End of file, %d tested",nb_tests);
            exit(-1);
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
        structure=malloc(sizeof(int)*MAX*2);
        index=0;
        len=getline(&correct_score,&b_len, fp);
        if (len != -1) {
            if (correct_score[len - 1] == '\n') {
                correct_score[len - 1] = '\0'; // Replace newline character with null terminator
                len--; // Decrement the length to exclude the removed '\n'
            }
            // Convert the line to an integer
            
        }
        else{
            printf("No integer found for test %d\n",nb_tests);
            exit(-1);
        }
        int number = atoi(correct_score);
        printf("Test: %d Size of the Sequence: %d ---", nb_tests,MAX);
        int score = fold(hashTable);
        if(score == number){
            printf("Correct\n");
        }
        else{
            printf("Failed (%d found, should be %d)\n",score,number);
            //exit(nb_tests);
        }
        int btscore=0;
        backtrace(hashTable,score);
        for(int i=0;i<index;i+=2){
            printf("%d->%d,",structure[i],structure[i+1]);
            btscore+=bp_score(structure[i],structure[i+1]);
        }
        if(btscore!=score){
            printf("Issue with backtrack (%d instead of %d)",btscore,score);
        }
        printf("\n");
        //print_table(hashTable);
        destroyHashTable(hashTable);
        free(structure);
        free(line);
    }
    return 0;
}



int compute_CLIQUE(HashTable *hashTable,int i, int j, int k, int l) {
    int CLIQUE=1;
    int value;
    int tab[] = {CLIQUE,i,j,k,l};
    int size= 5;

    if(j>=k || j<i || k>l){
        return INT_MAX;
    }

    if (get(hashTable,tab,size,&value)) { 
        return value;
    }
    
    int min_value = INT_MAX;

    if (k < l) { 
        min_value = min(min_value, compute_CLIQUE2(hashTable,i,j,k,l-1)); 
    }
    if (i < j) {
        min_value = min(min_value, compute_CLIQUE(hashTable,i+1,j,k,l));
    }
    if (i < j && k < l) {
        min_value = min(min_value, add(compute_CLIQUE(hashTable,i+1, j, k, l-1),bp_score(i,l)));
    }
    if (k==l) { 
        min_value = min(min_value, bp_score(i,l));
    }
    
    insert(hashTable,tab,size,min_value);
    return min_value;
}

int compute_CLIQUE2(HashTable *hashTable,int i, int j, int k, int l) {
    int CLIQUE2=-1;
    int value;
    int tab[] = {CLIQUE2,i,j,k,l};
    int size= 5;

    if(j>=k || j<i || k>l){
        return INT_MAX;
    }

    if (get(hashTable,tab,size,&value)) {
        return value;
    }

    int min_value = INT_MAX;

    if (k < l) { 
        min_value = min(min_value, compute_CLIQUE2(hashTable,i,j,k,l-1)); 
    }
    if (i < j && k < l) {
        min_value = min(min_value,add( 
                    compute_CLIQUE(hashTable,i+1, j, k, l-1),bp_score(i,l)));
    }
    
    if (k==l) { 
        min_value = min(min_value, bp_score(i, l));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

void compute_C(HashTable *hashTable,int a, int f, int c,int d) {
    int C = 67;
    int C2=-67;
    int value;

    int tab[]={C,a,f,c,d};
    int size = 5;

    if(a<0 || f<0 || f>=MAX || a>=MAX || a>f){
        return INT_MAX;
    }

    int min_value = INT_MAX;
    bool eq_some_const1 = (f>c-1) && (f>d-1);
    bool eq_some_const2 = (a<c) && (a<d);

    if (a+1 < f && eq_some_const2) {
            int tab_C1[]={C,a+1,f,c,d};
            int size_C1=5;
            int value_C1;
            get(hashTable,tab_C1,size_C1,&value_C1);
            min_value = min(min_value, value_C1);
    }
    if (f-1 > a && eq_some_const1) {
            int tab_C2[]={C2,a,f-1,c,d};
            int size_C2=5;
            int value_C2;
            get(hashTable,tab_C2,size_C2,&value_C2);
            min_value = min(min_value, value_C2);
    }
    if (eq_some_const1 && eq_some_const2) {
        int tab_C3[]={C,a+1,f-1,c,d};
        int size_C3=5;
        int value_C3;
        get(hashTable,tab_C3,size_C3,&value_C3);
        min_value = min(min_value, add(value_C3,bp_score(a,f)));
    }

    min_value = min(min_value, 0);

    insert(hashTable,tab,size,min_value);
} 

void compute_C2(HashTable *hashTable,int a, int f, int c,int d) {
    int C2 = -67;
    int C = 67;
    int value;
    int tab[]={C2,a,f,c,d};
    int size = 5;

    if(a<0 || f<0 || f>=MAX || a>=MAX || a>f){
        return INT_MAX;
    }

    int min_value = INT_MAX;
    bool eq_some_const1 = (f>c-1) && (f>d-1);
    bool eq_some_const2 = (a<c) && (a<d);

    if (f-1 > a && eq_some_const1) {
            int tab_C2[]={C2,a,f-1,c,d};
            int size_C2=5;
            int value_C2;
            get(hashTable,tab_C2,size_C2,&value_C2);
            min_value = min(min_value, value_C2);
    }
    if (eq_some_const1 && eq_some_const2) {
        int tab_C3[]={C,a+1,f-1,c,d};
        int size_C3=5;
        int value_C3;
        get(hashTable,tab_C3,size_C3,&value_C3);
        min_value = min(min_value, add(value_C3,bp_score(a,f)));
    }

    min_value = min(min_value, 0);

    insert(hashTable,tab,size,min_value);
}

void compute_B(HashTable *hashTable,int c,int d,int f,int h) {
    int B = 66;
    int value;
    int tab[]={B,c,d,f,h};
    int size = 5;

    int min_value = INT_MAX;
    //modulable
    int CLIQUE2=-1;
    int size_CLIQUE2=5;
    //modulable
    for (int g=f;g<h;g++) {
        //modulable
        int value_CLIQUE2;
        int tab_CLIQUE2[] = {CLIQUE2,c,d-1,g,h-1};
        get(hashTable,tab_CLIQUE2,size_CLIQUE2,&value_CLIQUE2);
        //modulable
        min_value = min(min_value, value_CLIQUE2);
    }

    insert(hashTable,tab,size,min_value);
}

int compute_A(HashTable *hashTable,int START,int h) {
    int A = 65;
    int value;
    int tab[]={A,h};
    int size = 2;

    
    int min_value = INT_MAX;
    //modulable
    int B=66;
    int size_B=5;

    int C2=-67;
    int size_C2=6;
    //modulable

    for (int a=START;a<h-3;a++) {
        for (int c=a+1;c<h-2;c++) {
            for (int d=c+1;d<h-1;d++) {
                for (int f=d+1;f<h;f++) {

                    int value_B;
                    int tab_B[] = {B,c,d,f,h};
                    get(hashTable,tab_B,size_B,&value_B);

                    int value_C2;
                    int tab_C2[] = {C2,a,f-1,c,d};
                    get(hashTable,tab_C2,size_C2,&value_C2);

                    min_value = min(min_value, add(value_B,value_C2));
                }
            }
        }
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

void backtrace_CLIQUE(HashTable *hashTable, int score, int i,int j, int k, int l){

    if(j>=k || j<i || k>l){
        return;
    }

    if (k < l) { 
        if(score==compute_CLIQUE2(hashTable,i,j,k,l-1)){
            backtrace_CLIQUE2(hashTable,score,i,j,k,l-1);
            return;
        }
    }
    if (i < j) {
        if(score==compute_CLIQUE(hashTable,i+1,j,k,l)){
            backtrace_CLIQUE(hashTable,score,i+1,j,k,l);
            return;
        }
    }
    if (i < j && k < l) {
        if(score == add(compute_CLIQUE(hashTable,i+1, j, k, l-1),bp_score(i,l))){
            backtrace_CLIQUE(hashTable,compute_CLIQUE(hashTable,i+1, j, k, l-1),i+1,j,k,l-1);
            structure[index]=i;
            index+=1;
            structure[index]=l;
            index+=1;
            return;
        }
    }
    if (k==l) { 
        if(score== bp_score(i, l)){
            structure[index]=i;
            index+=1;
            structure[index]=l;
            index+=1;
            return;
        }
    }
}

void backtrace_CLIQUE2(HashTable *hashTable,int score,int i,int j,int k,int l){

    if(j>=k || j<i || k>l){
        return;
    }

    if (k < l) { 
        if(score==compute_CLIQUE2(hashTable,i,j,k,l-1)){
            backtrace_CLIQUE2(hashTable,score,i,j,k,l-1);
            return;
        }
    }
    if (i < j && k < l) {
        if(score == add(compute_CLIQUE(hashTable,i+1, j, k, l-1),bp_score(i,l))){
            backtrace_CLIQUE(hashTable,compute_CLIQUE(hashTable,i+1, j, k, l-1),i+1,j,k,l-1);
            structure[index]=i;
            index+=1;
            structure[index]=l;
            index+=1;
            return;
        }
    }
    
    if (k==l) {
        if(score== bp_score(i, l)){
            structure[index]=i;
            index+=1;
            structure[index]=l;
            index+=1;
            return;
        }
    }

}

void backtrace_C(HashTable *hashTable,int score,int a, int f, int c,int d) {

    if(a<0 || f<0 || f>=MAX || a>=MAX || a>f){
        return;
    }

    bool eq_some_const1 = (f>c-1) && (f>d-1);
    bool eq_some_const2 = (a<c) && (a<d);

    if (a+1 < f && eq_some_const2) {
        if(score==compute_C(hashTable,a+1, f, c,d)){
                backtrace_C(hashTable,compute_C(hashTable,a+1, f, c,d),a+1, f, c,d);
                return;
        }
    }
    if (f-1 > a && eq_some_const1) {
            if(score==compute_C2(hashTable,a, f-1, c,d)){
            backtrace_C2(hashTable,compute_C2(hashTable,a, f-1, c,d),a, f-1, c,d);
            return;
        }
    }
    if (eq_some_const1 && eq_some_const2) {
        if(score==add(compute_C(hashTable,a+1, f-1, c,d),bp_score(a,f))){
            backtrace_C(hashTable,compute_C(hashTable,a+1, f-1, c,d),a+1, f-1, c,d);
            structure[index]=a;//be careful about incrementation
            index+=1;
            structure[index]=f;
            index+=1;
            return;
        }
        
    }
    if(score==0){
        
        return;
    }

    return;
} 

void backtrace_C2(HashTable *hashTable,int score, int a, int f, int c,int d) {
    
    if(a<0 || f<0 || f>=MAX || a>=MAX || a>f){
        return;
    }

    bool eq_some_const1 = (f>c-1) && (f>d-1);
    bool eq_some_const2 = (a<c) && (a<d);
        if (f-1 > a && eq_some_const1) {
            if(score==compute_C2(hashTable,a, f-1, c,d)){
            backtrace_C2(hashTable,compute_C2(hashTable,a, f-1, c,d),a, f-1, c,d);
            return;
        }
    }

     if (eq_some_const1 && eq_some_const2) {
        if(score==add(compute_C(hashTable,a+1, f-1, c,d),bp_score(a,f))){
            backtrace_C(hashTable,compute_C(hashTable,a+1, f-1, c,d),a+1, f-1, c,d);
            structure[index]=a;//be careful about incrementation
            index+=1;
            structure[index]=f;
            index+=1;
            return;
        }
        
    }
    return;
}

int backtrace_B(HashTable *hashTable,int score,int c,int d,int f,int h) {

    for (int g=f;g<h;g++) {
         if(score==compute_CLIQUE2( hashTable,c,d-1,g,h-1)){
    backtrace_CLIQUE2(hashTable,compute_CLIQUE2( hashTable,c,d-1,g,h-1), c,d-1,g,h-1);

            return;
         }
    }


}

int backtrace_A(HashTable *hashTable,int score,int START,int h) {

    for (int a=START;a<h-3;a++) {
        for (int c=a+1;c<h-2;c++) {
            for (int d=c+1;d<h-1;d++) {
                for (int f=d+1;f<h;f++) {
                     if(score==add(compute_B( hashTable,c,d,f,h),compute_C2( hashTable,a,f-1,c,d))){
                backtrace_B(hashTable,compute_B( hashTable,c,d,f,h), c,d,f,h);
                backtrace_C2(hashTable,compute_C2( hashTable,a,f-1,c,d), a,f-1,c,d);

                        return;
                     }
                }
            }
        }
    }


}

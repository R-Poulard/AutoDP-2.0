#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>

#define LOAD_FACTOR_THRESHOLD 0.7
#define TUPLE_SIZE 6

int TABLE_SIZE = 4001;
int N = 10;

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

Tuple *createTuple(int val[]) {
    Tuple *tuple = (Tuple *)malloc(sizeof(Tuple));

    // Initialize all values in tuple to 0
    for (int i = 0; i < TUPLE_SIZE; i++) {
        tuple->values[i] = 0;
    }

    // Copy values from val array, up to its length or TUPLE_SIZE
    for (int i = 0; i < TUPLE_SIZE; i++) {
        if (val[i] != '\0') {
            tuple->values[i] = val[i];
        }
        else {
            break; // Stop copying if val[i] is '\0' indicating end of values
        }
    }

    return tuple;
}

int compare_tuple(Tuple *tp1, int val[]) {
    int i;
    // Compare elements until the end of val or tp1->values
    for (i = 0; i < TUPLE_SIZE && val[i] != '\0'; i++) {
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
int hash(int keys[], HashTable *hashtable) {
    int key = 0;
    for (int i = 0; i < TUPLE_SIZE; i++) {
        if (keys[i] != '\0') {
            break;
        }
        key += keys[i] * pow(N, TUPLE_SIZE-i);
    }
    return key % hashtable->capacity;
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
void insert(HashTable *hashTable, int keys[], int value) {
    // Combine the numbers into a single key
    int index = hash(keys, hashTable);

    // Create a new node
    HashNode *newNode = (HashNode *)malloc(sizeof(HashNode));
    newNode->key = createTuple(keys);
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
        printf("%d resize \n", hashTable->capacity);
        int new_capacity = hashTable->capacity * 2;
        HashNode **new_buckets = (HashNode **)calloc(new_capacity, sizeof(HashNode *));
        for (int i = 0; i < hashTable->capacity; i++) {
            HashNode *current = hashTable->buckets[i];
            while (current != NULL) {
                HashNode *temp = current;
                current = current->next;
                int new_index = hash(temp->key->values, hashTable);
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
bool get(HashTable *hashTable, int keys[], int *value) {
    int index = hash(keys, hashTable);
    
    // Traverse the linked list at the index
    HashNode *current = hashTable->buckets[index];
    while (current != NULL) {
        if (compare_tuple(current->key, keys)) {
            *value = current->value;
            return true; // Key found
        }
        current = current->next;
    }
    return false; // Key not found
}

//
int min(int a, int b) { if (a<b) {return a;} else {return b;}};

char * line = NULL;
int MAX;



//declarations
int compute_CLIQUE(HashTable *hashTable, int i, int j, int k, int l);

int compute_CLIQUE2(HashTable *hashTable, int i, int j, int k, int l);

int compute_M(HashTable *hashTable,int m,int n,int r,int t) ;

int compute_L(HashTable *hashTable,int i,int k,int o,int p) ;

int compute_K(HashTable *hashTable,int i,int k,int n,int p) ;

int compute_J(HashTable *hashTable,int a,int h,int k,int n,int p) ;

int compute_I(HashTable *hashTable,int b,int e,int g,int k,int l) ;

int compute_H(HashTable *hashTable,int c,int e,int q,int r) ;

int compute_G(HashTable *hashTable,int c,int e,int p,int r) ;

int compute_F(HashTable *hashTable,int b,int e,int p,int r) ;

int compute_E(HashTable *hashTable,int b,int g,int k,int l,int p,int r) ;

int compute_D2(HashTable *hashTable,int h, int a, int l,int p,int k,int r) ;

int compute_D(HashTable *hashTable,int h, int a, int l,int p,int k,int r) ;

int compute_C(HashTable *hashTable,int a,int l,int n,int r) ;

int compute_B(HashTable *hashTable,int a,int m,int n,int r) ;

int compute_A(HashTable *hashTable) ;

int bp_score(char x, char y) {
   if (x=='G' && y=='C') { return 10; }
    if (x=='C' && y=='G') { return 10; }
    if (x=='G' && y=='U') { return 5; }
    if (x=='U' && y=='G') { return 5; }
    if (x=='A' && y=='U') { return 5; }
    if (x=='U' && y=='A') { return 5; }
    return 0;
}

int fold(HashTable * hashTable) {
    compute_A(hashTable);
}


int main(int argc, char ** argv) {
    printf("File name: %s\n",argv[1]);
    HashTable *hashTable = createHashTable();
    size_t b_len=0;
    FILE * fp = fopen(argv[1], "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);
    int len=getline(&line,&b_len, fp); 
    
    preset_N(len);
    MAX=len;
    printf("Sequence: %s\n Size of the Sequence: %d\n", line,len);
    int score = fold(hashTable);
    printf("%d\n",score);
    //char * structure = NULL; not implemented yet
    //structure = backtrace();
    destroyHashTable(hashTable);
    return score;
}


int compute_CLIQUE(HashTable *hashTable,int i, int j, int k, int l) {
    int CLIQUE=1;
    
    int value;
    int tab[] = {CLIQUE,i,j,k,l};
    if (get(hashTable,tab,&value)) { 
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
        min_value = min(min_value, 
                    compute_CLIQUE(hashTable,i+1, j, k, l-1)+bp_score(line[i],line[l]));
    }
    if (k==l) { 
        min_value = min(min_value, bp_score(line[i], line[l]));
    }
    
    insert(hashTable,tab,min_value);
    return min_value;
}

int compute_CLIQUE2(HashTable *hashTable,int i, int j, int k, int l) {
    int CLIQUE2=-1;

    int value;
    int tab[] = {CLIQUE2,i,j,k,l};
    if (get(hashTable,tab,&value)) { 
        return value;
    }

    int min_value = INT_MAX;

    if (k < l) { 
        min_value = min(min_value, compute_CLIQUE2(hashTable,i,j,k,l-1)); 
    }
    if (i < j && k < l) {
        min_value = min(min_value, 
                    compute_CLIQUE(hashTable,i+1, j, k, l-1)+bp_score(line[i],line[l]));
    }
    if (k==l) { 
        min_value = min(min_value, bp_score(line[i], line[l]));
    }

    insert(hashTable,tab,min_value);
    return min_value;
}

int compute_M(HashTable *hashTable,int m,int n,int r,int t) {
    int M = 77;
    int value;
    int tab[]={M,m,n,r,t};
    if (get(hashTable,tab,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int s=r;s<t;s++) {
        min_value = min(min_value, compute_CLIQUE( hashTable,m,n,s,t));
    }

    insert(hashTable,tab,min_value);
    return min_value;
}

int compute_L(HashTable *hashTable,int i,int k,int o,int p) {
    int L = 76;
    int value;
    int tab[]={L,i,k,o,p};
    if (get(hashTable,tab,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int j=i;j<p;j++) {
        min_value = min(min_value, compute_CLIQUE( hashTable,i,j,o,p));
    }

    insert(hashTable,tab,min_value);
    return min_value;
}

int compute_K(HashTable *hashTable,int i,int k,int n,int p) {
    int K = 75;
    int value;
    int tab[]={K,i,k,n,p};
    if (get(hashTable,tab,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int o=n;o<p;o++) {
        min_value = min(min_value, compute_L( hashTable,i,k,p,o));
    }

    insert(hashTable,tab,min_value);
    return min_value;
}

int compute_J(HashTable *hashTable,int a,int h,int k,int n,int p) {
    int J = 74;
    int value;
    int tab[]={J,a,h,k,n,p};
    if (get(hashTable,tab,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int i=h;i<p;i++) {
        min_value = min(min_value, compute_K( hashTable,n,i,p,k));
    }

    insert(hashTable,tab,min_value);
    return min_value;
}

int compute_I(HashTable *hashTable,int b,int e,int g,int k,int l) {
    int I = 73;
    int value;
    int tab[]={I,b,e,g,k,l};
    if (get(hashTable,tab,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int f=e;f<l;f++) {
        min_value = min(min_value, compute_CLIQUE( hashTable,e,f,k,l));
    }

    insert(hashTable,tab,min_value);
    return min_value;
}

int compute_H(HashTable *hashTable,int c,int e,int q,int r) {
    int H = 72;
    int value;
    int tab[]={H,c,e,q,r};
    if (get(hashTable,tab,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int d=0;d<r;d++) {
        min_value = min(min_value, compute_CLIQUE( hashTable,c,d,q,r));
    }

    insert(hashTable,tab,min_value);
    return min_value;
}

int compute_G(HashTable *hashTable,int c,int e,int p,int r) {
    int G = 71;
    int value;
    int tab[]={G,c,e,p,r};
    if (get(hashTable,tab,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int q=p;q<r;q++) {
        min_value = min(min_value, compute_H( hashTable,c,q,e,r));
    }

    insert(hashTable,tab,min_value);
    return min_value;
}

int compute_F(HashTable *hashTable,int b,int e,int p,int r) {
    int F = 70;
    int value;
    int tab[]={F,b,e,p,r};
    if (get(hashTable,tab,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int c=r;c<MAX;c++) {
        min_value = min(min_value, compute_G( hashTable,c,p,e,r));
    }

    insert(hashTable,tab,min_value);
    return min_value;
}

int compute_E(HashTable *hashTable,int b,int g,int k,int l,int p,int r) {
    int E = 69;
    int value;
    int tab[]={E,b,g,k,l,p,r};
    if (get(hashTable,tab,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int e=0;e<r;e++) {
        min_value = min(min_value, compute_F( hashTable,p,b,e,r)+compute_I( hashTable,k,g,l,e,b));
    }

    insert(hashTable,tab,min_value);
    return min_value;
}

int compute_D(HashTable *hashTable,int h, int a, int l,int p,int k,int r) {
    int D = 68;
    int value;

    int tab[]={D,l,p,k,r};
    if (get(hashTable,tab,&value)){
        return value;
    }

    int min_value = INT_MAX;
    bool eq_some_const1 = (a-1!=l) && (a-1!=p) && (a-1!=k) && (a-1!=r);
    bool eq_some_const2 = (a-1!=l) && (a-1!=p) && (a-1!=k) && (a-1!=r);

    if (h+1 < a) {
        if (!eq_some_const1) {
            min_value = min(min_value, compute_D(hashTable,h+1,a,l,p,k,r));
        }
    }
    if (a-1 > h) {
        if (!eq_some_const2) {
            min_value = min(min_value, compute_D2(hashTable,h, a-1, l,p,k,r));
        }
    }
    if (!eq_some_const1 && !eq_some_const2) {
        min_value = min(min_value, compute_D(hashTable,h+1, a-1, l,p,k,r)+bp_score(line[h],line[a]));
    }

    min_value = min(min_value, compute_E( hashTable,k,20,l,r,p,1));

    insert(hashTable,tab,min_value);
    return min_value;
} 

int compute_D2(HashTable *hashTable,int h, int a, int l,int p,int k,int r) {
    int D2 = -68;
    int value;

    int tab[]={D2,h,a,l,p,k,r};
    if (get(hashTable,tab,&value)){
        return value;
    }

    int min_value = INT_MAX;
    bool eq_some_const1 = (a-1!=l) && (a-1!=p) && (a-1!=k) && (a-1!=r);
    bool eq_some_const2 = (a-1!=l) && (a-1!=p) && (a-1!=k) && (a-1!=r);

    if (a-1 > h && !eq_some_const2) {
        min_value = min(min_value, compute_D2(hashTable,h, a-1, l,p,k,r));
    }
    if (!eq_some_const1 && !eq_some_const2) {
        min_value = min(min_value, compute_D(hashTable,h+1, a-1, l,p,k,r)+bp_score(line[h],line[a]));
    }

    insert(hashTable,tab,min_value);
    return min_value;
}

int compute_C(HashTable *hashTable,int a,int l,int n,int r) {
    int C = 67;
    int value;
    int tab[]={C,a,l,n,r};
    if (get(hashTable,tab,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int h=a;h<r;h++) {
        for (int k=h;k<r;k++) {
            for (int p=k;p<r;p++) {
                min_value = min(min_value, compute_D( hashTable,a,h,k,l,p,r)+compute_J( hashTable,k,h,a,p,n));
            }
        }
    }

    insert(hashTable,tab,min_value);
    return min_value;
}

int compute_B(HashTable *hashTable,int a,int m,int n,int r) {
    int B = 66;
    int value;
    int tab[]={B,a,m,n,r};
    if (get(hashTable,tab,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int l=a;l<r;l++) {
        min_value = min(min_value, compute_C( hashTable,l,n,a,r));
    }

    insert(hashTable,tab,min_value);
    return min_value;
}

int compute_A(HashTable *hashTable) {
    int min_value = INT_MAX;
    
    for (int a=0;a<MAX;a++) {
        for (int m=a;m<MAX;m++) {
            for (int n=m;n<MAX;n++) {
                for (int r=n;r<MAX;r++) {
                    for (int t=r;t<MAX;t++) {
                        min_value = min(min_value, compute_B( hashTable,a,n,m,r)+compute_M( hashTable,t,n,m,r));
                    }
                }
            }
        }
    }

    return min_value;
}

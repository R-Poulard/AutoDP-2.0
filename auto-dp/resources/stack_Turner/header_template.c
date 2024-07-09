#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <string.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/loop_energies.h>

//HASHING TABLE FUNCTION
#define LOAD_FACTOR_THRESHOLD 0.7
#define TUPLE_SIZE ##MODIFY_SIZE_TUPLE##

int TABLE_SIZE = 64019;
int N = 10;


typedef struct Tuple {
    int values[TUPLE_SIZE];
} Tuple;

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

typedef struct Node {
    Tuple* data;
    struct Node* next;
} Node;

// Define the Stack structure
PUBLIC typedef struct {
    Node* top;
} Stack;

PUBLIC typedef struct {
    HashTable* hashtable;
    Stack* Stack;
    vrna_fold_compound_t* fc;
    const char* line;
    int** matrix;
    int MAX;
    int THETA;
    int THETA_PAIRING;
} pk_compound;

PUBLIC typedef struct{
    char * structure;
    char bracket;
} bt_struct;


//PAIRING FUNCTIONS
int min(int a, int b) { if (a<b) {return a;} else {return b;}};
int max(int a, int b) { if (a<b) {return b;} else {return a;}};

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

// Hash function
int hash(int keys[], int capacity,int size) {
    int key = 0;
    for (int i = 0; i < size; i++) {
        key += keys[i] * (int)pow(N, TUPLE_SIZE-i-1);
    }
    return abs(key % capacity);
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
    int index = hash(keys, hashTable->capacity,size);
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
                int new_index = hash(temp->key->values, new_capacity,size);
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
    
    int index = hash(keys, hashTable->capacity,size);

    // Traverse the linked list at the index
    HashNode *current = hashTable->buckets[index];
    while (current != NULL) {
        if (compare_tuple(current->key, keys,size)) {
            
            *value = current->value;
            
            //printf("\nvalue in hash +> %d",current->value);
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

void print_table_stat(HashTable *hashTable){
    printf("#-#--TABLE STAT REPORT-----\n#-#nb of elem: %d\n#-#capcaity: %d\n#-#initial capacity %d\n",hashTable->size,hashTable->capacity,TABLE_SIZE);
    
    int nb_chaines=0;
    int biggest_chaines=0;
    int sum_chaines=0;
    for (int i = 0; i < hashTable->capacity; i++) {
        HashNode *current = hashTable->buckets[i];
        int chaine_size=0;
        if(current!=NULL && current->next!=NULL){
            nb_chaines++;
        }
        while (current != NULL) {      
            chaine_size++;
            current = current->next;
        }
        biggest_chaines=max(biggest_chaines,chaine_size);
        sum_chaines+=chaine_size;
    }
    printf("#-#Biggest chain: %d\n#-# prct chaine: %.3f (%d/%d)\n#-# mean chaine size: %.3f\n",biggest_chaines,(double)nb_chaines/hashTable->size,nb_chaines,hashTable->size,sum_chaines/nb_chaines);
    printf("#-#--------------------------\n");
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

void createMatrix(int *** matrix, int N,const char * ss) {

    // Allocate memory for the matrix
    *matrix = (int **)malloc(N * sizeof(int *));
    for (int i = 0; i < N; i++) {
        (*matrix)[i] = (int *)malloc(N * sizeof(int));
    }

    // Fill the matrix
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if(i==j){
                if(ss[i]=='|'){
                    (*matrix)[i][j]=0;
                }
                else{
                    (*matrix)[i][j]=1;
                }
                continue;
            }
            (*matrix)[i][j] = 0;
            (*matrix)[j][i] = 0;
            if ( ss[i]=='x' || ss[j]=='x') {
                (*matrix)[i][j] = 0;
                (*matrix)[j][i] = 0;
                continue;
            }
            if( (ss[i]=='(' && ss[j]==')') || (ss[j]=='(' && ss[i]==')')){
                (*matrix)[i][j]= 1;
                (*matrix)[j][i] = 1;
                continue;
            }
            if( ss[i]=='>'){
                (*matrix)[i][j]= 0;
                (*matrix)[j][i]= 0;
                continue;
            }
            if( (ss[i]=='.' || ss[i]=='<') && (ss[j]=='.' || ss[j]=='>')){
                (*matrix)[i][j]= 1;
                (*matrix)[j][i]= 1;
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
                        (*matrix)[i][j]=1;
                        (*matrix)[j][i]=1;
                    }
                    else{
                        (*matrix)[i][j]=0;
                        (*matrix)[j][i]=0;
                    }
                }
                else{
                    (*matrix)[i][j] = 0;
                    (*matrix)[j][i] = 0;
                }
                continue;
            }
        }
    }
}

// Function to display the matrix
void displayMatrix(int** matrix,int N) {
    printf("Matrix:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

// Function to free the memory allocated for the matrix
void freeMatrix(int ** matrix,int N) {
    for (int i = 0; i < N; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

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

int evaluate(pk_compound* pk,int x, int y) {
    int a=abs(x-y)<pk->THETA_PAIRING;
    int b=pk->matrix[x][y];
    int c=matching(pk->line[x],pk->line[y]);

    if(abs(x-y)<pk->THETA_PAIRING || pk->matrix[x][y]==0 || matching(pk->line[x],pk->line[y])==INT_MAX){
        return 0;
    }
    return 1;
}


int add(int a,int b){
    if(a==INT_MAX || b==INT_MAX){
        return INT_MAX;
    }

    return a+b;
}

int write_structure(pk_compound* pk,bt_struct* bt,vrna_bp_stack_t *bp){
    int i,j,temp;
    for (int k = 1; k <= bp[0].i; k++) {
      i = bp[k].i;
      j = bp[k].j;
      if (i > pk->MAX)
        i -= pk->MAX;

      if (j > pk->MAX)
        j -= pk->MAX;

      if (i > j) {
        temp  = i;
        i     = j;
        j     = temp;
      }

      if (i == j) {
        /* Gquad bonds are marked as bp[i].i == bp[i].j */
        bt->structure[i - 1] = '+';
      } else {
        /* the following ones are regular base pairs */
        bt->structure[i - 1]  = '(';
        bt->structure[j - 1]  = ')';
      }
    }
}

int MFEFree(pk_compound* pk,int a,int b){
    if(a>b){
        return 0;
    }
    a=a++;
    b=b++;

    int * indx = pk->fc->jindx;
    int ij = indx[b] + a;
    int * fML =pk->fc->matrices->fML;
    int up_energy=pk->fc->params->MLbase*(b-a+1);
    return min(up_energy,fML[ij]);
}

int compute_BT(vrna_fold_compound_t *fc, int i,int j, vrna_bp_stack_t *bp_stack){
    int ret;
    int p,q,comp1,comp2;
    int s=0;
    int b=bp_stack[0].i;
    sect bt_stack[500]; /* stack of partial structures for backtracking */
    //printf("a=%d,b=%d\n",i,j);
    int BT=vrna_BT_mb_loop_split(fc, &i, &j, &p, &q, &comp1, &comp2, bp_stack, &b);
    if(BT!=1){
        return BT;
    }
    if (i > 0) {
        bt_stack[++s].i = i;
        bt_stack[s].j   = j;
        bt_stack[s].ml  = comp1;
    }

    if (p > 0) {
        bt_stack[++s].i = p;
        bt_stack[s].j   = q;
        bt_stack[s].ml  = comp2;
    }
    BT=vrna_backtrack_from_intervals(fc,bp_stack,bt_stack,s);
    return BT;
}

void backtrace_MFEFree(pk_compound* pk,bt_struct* bt,int score, int a,int b){
    int up_energy=pk->fc->params->MLbase*(b-a+1);
    if(a>b || score==up_energy){
        return;
    }
    a++;
    b++;
    vrna_bp_stack_t * bp=(vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * ( 1 + (b-a) / 2)));
    int score_bt=compute_BT(pk->fc,a,b,bp);
    write_structure(pk,bt,bp);
    free(bp);
}

int INTB(pk_compound* pk,int a,int b,int c,int d){
    if(d==-1 && c==-1){
        return 0;
    }    
    //printf("\na=%d,b=%d,c=%d,d=%d  ",a,b,c,d);
    a++;b++;
    int energy_full = vrna_eval_int_loop(pk->fc,c,d,a,b);
    //printf("=%d\n",energy_full);
    return energy_full;
}

void add_pairing(pk_compound* pk,bt_struct* bt, int score,int a,int b,int c,int d){
    bt->structure[a]=bt->bracket;
    bt->structure[b]=tolower(bt->bracket);
    return;
}

void freeNode(Node * node){
    free(node->data);
    free(node);
}

// Function to initialize the stack
void initializeStack(Stack* stack) {
    stack->top = NULL;
}

// Function to check if the stack is empty
bool isStackEmpty(Stack* stack) {
    return stack->top == NULL;
}

// Function to push a tuple onto the stack
void push(Stack* stack,int tab[],int size) {
    Node* newNode = (Node*)malloc(sizeof(Node));
    if (newNode == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    newNode->data = createTuple(tab,size);
    newNode->next = stack->top;
    stack->top = newNode;
}

// Function to pop a tuple from the stack
bool pop(Stack* stack, Tuple* item) {
    if (isStackEmpty(stack)) {
        //printf("Stack is empty. Cannot pop item.\n");
        return false;
    }
    Node* temp = stack->top;
    item = temp->data;
    stack->top = stack->top->next;
    free(temp);
    return true;
}


// Function to free the stack
void freeStack(Stack* stack) {
    Node* current = stack->top;
    Node* next;
    while (current != NULL) {
        next = current->next;
        free(current);
        current = next;
    }
    stack->top = NULL;
    free(stack);
}
void printStack(Stack* stack) {
    Node* current = stack->top;
    
    Node* next;
    printf("\n----Stack print-------\n");
    while (current != NULL) {
        next = current->next;
        print_tuple(current->data);
        printf("\n");
        current = next;
    }
    printf("--------------------\n");
}


pk_compound* create_pk_compound(vrna_fold_compound_t *fc, char* seq,const char* ss){
    pk_compound *pk=malloc(sizeof(pk_compound));
    pk->fc=fc;
    pk->hashtable=createHashTable();
    pk->Stack=malloc(sizeof(Stack));
    initializeStack(pk->Stack);
    pk->line=seq;
    pk->MAX=strlen(seq)-1;
    createMatrix(&(pk->matrix),pk->MAX,ss);
    preset_N(pk->MAX);
    pk->THETA=100;
    pk->THETA_PAIRING=1;
}
void free_pk(pk_compound* pk){
    destroyHashTable(pk->hashtable);
    freeMatrix(pk->matrix,pk->MAX);
    freeStack(pk->Stack);
    vrna_fold_compound_free(pk->fc);
    free(pk);
}

bt_struct * create_bt_struct(char* structure,int MAX){
    bt_struct* bt=malloc(sizeof(bt_struct));
    bt->bracket='A';
    if(structure==NULL){
        bt->structure=malloc(sizeof(char)*(MAX+1));
        bt->structure[MAX]='\0';
        for(int i=0;i<MAX;i++){
            bt->structure[i]='.';
        }
    }
    else{
        bt->structure=structure;
    }
    return bt;
}
void free_bt(bt_struct* bt){
    free(bt->structure);
    free(bt);
}
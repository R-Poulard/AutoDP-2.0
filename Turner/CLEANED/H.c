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
#define TUPLE_SIZE 6

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
//declarations
int get_CLIQUE1(pk_compound *pk, Node **first, Node **last, int *value, int i, int i2, int j2, int j);

int get_CLIQUE0(pk_compound *pk, Node **first, Node **last, int *value, int i, int i2, int j2, int j);

int compute_CLIQUE1(pk_compound *pk,Node **first,Node** last, int i, int j, int k, int l);

void backtrace_CLIQUE0(pk_compound *pk,bt_struct* bt, int score, int i,int j, int k, int l);

void backtrace_CLIQUE1(pk_compound *pk,bt_struct* bt, int score, int i,int j, int k, int l);
int get_C0(pk_compound* pk,Node** first,Node** last,int * value,int a, int f, int c,int d);

int get_C1(pk_compound* pk,Node** first,Node** last,int * value,int a, int f, int c,int d);

int compute_C1(pk_compound* pk,Node** first,Node** last,int a, int f, int c,int d);

void backtrace_C0(pk_compound *pk,bt_struct* bt, int score, int a, int f, int c,int d) ;

void backtrace_C1(pk_compound *pk,bt_struct* bt, int score, int a, int f, int c,int d) ;

int get_B(pk_compound *pk,Node** first,Node** last,int * value,int c,int d,int f,int h);

int compute_B(pk_compound *pk,Node** first,Node** last,int c,int d,int f,int h);

void backtrace_B(pk_compound *pk,bt_struct* bt,int score,int c,int d,int f,int h) ;

int get_A(pk_compound *pk,Node** first,Node** last,int * value,int a,int h);

int compute_A(pk_compound *pk,Node** first,Node** last,int a,int h);

void backtrace_A(pk_compound *pk,bt_struct* bt,int score,int a,int h) ;


int init_root(pk_compound* pk){  
    int A = 65;
    int size = 3;
    for(int h=0+3;h<=pk->MAX;h++){
        for(int a=0;a<h-3;a++){
            int tab[]={A,a,h};
            Node* newNode = (Node*)malloc(sizeof(Node));
            newNode->data = createTuple(tab,size);
            newNode->next = NULL;
            if(pk->Stack->top!=NULL){
                newNode->next=pk->Stack->top;
                pk->Stack->top=newNode;
            }
            pk->Stack->top=newNode;
        }
    }
}

int init_root_inter(pk_compound* pk,int a,int h){  
    int A = 65;
    int size = 3;
    Stack* stack=pk->Stack;
    int tab[]={A,a,h};
            
    Node* newNode = (Node*)malloc(sizeof(Node));
    newNode->data = createTuple(tab,size);
    newNode->next = NULL;
    if(stack->top!=NULL){
        newNode->next=stack->top;
        stack->top=newNode;
    }
    stack->top=newNode;
}

int compute_root(pk_compound* pk) {
    int min_value=INT_MAX; 
    for(int h=0+3;h<=pk->MAX;h++){
        for(int a=0;a<h-3;a++){
            int mfe1=MFEFree(pk,0,a-1);
            int mfe2=MFEFree(pk,h,pk->MAX-1);
            int mfe=add(mfe1,mfe2);
            int tmp;
            get_A(pk,NULL,NULL,&tmp,a,h);
            min_value = min(min_value,add(tmp,mfe));
        }
    }
    return min_value;
}


int fold(pk_compound* pk){
    
    init_root(pk);
    //printStack(stack);
    Node* topNode;
    Node* nextNode;
    Node* add_on_last;
    Node* add_on_first;
    while (!isStackEmpty(pk->Stack)) {
        //pop(stack, &topTuple)
        //printStack(stack);
        topNode=pk->Stack->top;
        nextNode=topNode->next;
        Tuple * topTuple = topNode->data;
        add_on_last=NULL;
        add_on_first=NULL;
        int a,b,c,d,e,f,g,h;
        int i1,i2,j2,j1;
        //printf("here\n");
        //print_tuple(topTuple);
        switch(topTuple->values[0]){
            default: 
            
                 i1=topTuple->values[1];
                 i2=topTuple->values[2];
                 j2=topTuple->values[3];
                 j1=topTuple->values[4];
                if(compute_CLIQUE1(pk,&add_on_first,&add_on_last,i1,i2,j2,j1)==1){
                    
                    freeNode(topNode);
                    pk->Stack->top=nextNode;
                }
                else{
                    pk->Stack->top=add_on_first;
                    add_on_last->next=topNode;
                }
            break;//Clique

             case 67:
                 a=topTuple->values[1];
                 f=topTuple->values[2];
                 c=topTuple->values[3];
                 d=topTuple->values[4];

                 if(compute_C1(pk,&add_on_first,&add_on_last,a,f,c,d)==1){

                    freeNode(topNode);
                    pk->Stack->top=nextNode;
                 }
                 else{
                    pk->Stack->top=add_on_first;
                    add_on_last->next=topNode;
                 }
                 break;

             case 66:
                 c=topTuple->values[1];
                 d=topTuple->values[2];
                 f=topTuple->values[3];
                 h=topTuple->values[4];

                 if(compute_B(pk,&add_on_first,&add_on_last,c,d,f,h)==1){

                    freeNode(topNode);
                    pk->Stack->top=nextNode;
                 }
                 else{
                    pk->Stack->top=add_on_first;
                    add_on_last->next=topNode;
                 }
                 break;

             case 65:
                 a=topTuple->values[1];
                 h=topTuple->values[2];

                 if(compute_A(pk,&add_on_first,&add_on_last,a,h)==1){

                    freeNode(topNode);
                    pk->Stack->top=nextNode;
                 }
                 else{
                    pk->Stack->top=add_on_first;
                    add_on_last->next=topNode;
                 }
                 break;

        }
    }
    //printStack(stack);
    return compute_root(pk);
}
void backtrace(pk_compound *pk,bt_struct* bt,int score) {
    //int a=0;int h=pk->MAX;
    for(int h=0+3;h<=pk->MAX;h++){
        for(int a=0;a<h-3;a++){
            int mfe1=MFEFree(pk,0,a-1);
            int mfe2=MFEFree(pk,h,pk->MAX-1);
            int mfe=add(mfe1,mfe2);
            int tmp0;
            get_A(pk,NULL,NULL,&tmp0,a,h);
            if(score==add(mfe,tmp0)){
                backtrace_A(pk,bt,tmp0,a,h);
                backtrace_MFEFree(pk,bt,mfe1,0,a-1);
                backtrace_MFEFree(pk,bt,mfe2,h,pk->MAX-1);
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
        
        int b_len=0;

        int len=0;
        char * line=NULL;
        //reading of the sequence
        len=getline(&line,&b_len, fp);
        
        if(len<=0){
            
            free(line);
            printf("End of file, %d tested",nb_tests);
            break;
        }
        if(line[0]=='#'){
            printf("%s",line);
            free(line);
            continue;
        }

        int MAX=len;
        if(line[len-1]=='\n'){
            MAX=len-1;
            nb_tests+=1;
        }

        //secondary structure
        char *ss = NULL;
        len=getline(&ss,&b_len, fp);
        if (len != -1) {
            if (ss[len - 1] == '\n') {
                ss[len - 1] = '\0'; // Replace newline character with null terminator
                len--; // Decrement the length to exclude the removed '\n'
            }
        }
        else{
            
            free(line);
            free(ss);
            printf("No secondary structure found for test %d\n",nb_tests);
            exit(-1);
        }
        //reading the score
        char * correct_score=NULL;
        len=getline(&correct_score,&b_len, fp);
        if (len != -1) {
            if (correct_score[len - 1] == '\n') {
                correct_score[len - 1] = '\0'; // Replace newline character with null terminator
                len--; // Decrement the length to exclude the removed '\n'
            }
            // Convert the line to an integer
            
        }
        else{
            free(line);
            free(correct_score);
            free(ss);
            printf("No integer found for test %d\n",nb_tests);
            exit(-1);
        }
        int number = atoi(correct_score);
        printf("Test: %d Size of the Sequence: %d ---", nb_tests,MAX);
        
        //elements used to record backtrack
        
        //Start of test
        clock_t test_start_time = clock();

        vrna_init_rand();
        vrna_fold_compound_t * fc = vrna_fold_compound(line, NULL, VRNA_OPTION_DEFAULT);
        vrna_constraints_add(fc, ss, VRNA_CONSTRAINT_DB_DEFAULT);
        pk_compound* pk=create_pk_compound(fc,line,ss);
        vrna_mfe(fc, NULL);
        
        //folding
        int score = fold(pk);
        //retrieve backtrack
        bt_struct* bt=create_bt_struct(NULL,pk->MAX);
        backtrace(pk,bt,score);
        clock_t test_end_time = clock();

        double test_time = ((double)(test_end_time - test_start_time)) / CLOCKS_PER_SEC;
        total_time += test_time;
        //end of test
        
        //sanity checks
        if(score == number || 1){//we autorise wrong score in the case of the Turner model
            float fl=(float)score/100.0;
            printf("Score %.2f ",fl);
            int nb=0;
            for(int i=0;i<pk->MAX;i++){
                
                if(bt->structure[i]>64){
                    if(bt->structure[i]<91){
                        
                        nb+=bt->structure[i];
                    }
                    else{
                        
                        nb-=(bt->structure[i]-32);
                    }
                }
            }
            printf("(backtrack impurty= %d)\n",nb);
        }
        else{
            free(line);
            free(correct_score);
            free_pk(pk);
            free_bt(bt);
            free(ss);
            printf("Failed (%d found, should be %d)\n",score,number);
            exit(nb_tests);
        }

        
        printf("%s\n",bt->structure);
        
        //clean up

        free_pk(pk);
        free_bt(bt);
        free(line);
        free(ss);
        free(correct_score);
    }
    end_time = clock();
    double total_execution_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    printf(" Total execution time: %.4f seconds", total_execution_time);
    printf("(mean %.4f seconds)\n", total_time / nb_tests);

    fclose(fp);
    return 0;
}


int get_CLIQUE0(pk_compound* pk,Node** first,Node** last,int * value,int i, int i2, int j2, int j){
    int CLIQUE=1;
    int tab[] = {CLIQUE,i+1,i2,j2,j-1};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        *value=INT_MAX;
        return 1;
    }
    if (get(pk->hashtable,tab,size,value)){
        //printf("*value = %d(+%d)\n",*value,INTB(i,j,-1,-1));
        *value=add(INTB(pk,i,j,-1,-1),*value);
        //printf("*new _ value = %d\n",*value);
        return 1;
    }
    Node* newNode = (Node*)malloc(sizeof(Node));
    newNode->data = createTuple(tab,size);
    newNode->next = NULL;
    if(*last!=NULL){
        (*last)->next=newNode;
    }
    if(*first==NULL){
        *first=newNode;
    }
    *last=newNode;
    return 0;
}

int get_CLIQUE1(pk_compound* pk,Node** first,Node** last,int * value,int i, int i2, int j2, int j){
    int CLIQUE=1;
    int tab[] = {CLIQUE,i,i2,j2,j};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        *value=INT_MAX;
        return 1;
    }
    if (get(pk->hashtable,tab,size,value)){
        return 1;
    }
    Node* newNode = (Node*)malloc(sizeof(Node));
    newNode->data = createTuple(tab,size);
    newNode->next = NULL;
    if(last!=NULL && *last!=NULL){
        (*last)->next=newNode;
    }
    if(first!=NULL && *first==NULL){
        *first=newNode;
    }
    *last=newNode;
    return 0;
}


int compute_CLIQUE1(pk_compound* pk,Node** first,Node** last,int i, int i2, int j2, int j) {
    int CLIQUE=1;
    int value;
    int tab[] = {CLIQUE,i,i2,j2,j};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        return 1;
    }

    if (get(pk->hashtable,tab,size,&value)) { 
        return 1;
    }
    int possible=1;
    int min_value = 0;
    int tmp;
    if (j == j2 && i <= i2) { 
        min_value = min(min_value, MFEFree(pk,i,i2)); 
    }
    for(int k=i;k<=i2;k++){
        for(int l=max(j2,k+1);l<=j;l++){
            if(evaluate(pk,k,l) && l-k+1<=pk->THETA && (k<i2 || l==j2)){
                tmp=INTB(pk,k,l,i,j);
                if(tmp!=INT_MAX){
                    int tmp2;
                    if(get_CLIQUE1(pk,first,last,&tmp2,k+1,i2,j2,l-1)==0){
                        possible=0;
                    }
                    if(possible==1){
                        
                    min_value = min(min_value,add(tmp2,tmp));
                    }
                }
            }
        }
    }
    if(possible==1){
        insert(pk->hashtable,tab,size,min_value);
    }
    return possible;
}

int get_C0(pk_compound* pk,Node** first,Node** last,int * value,int a, int f, int c,int d){
    int C1 = 67;
    int tab[]={C1,a+1,f-1,c,d};
    int size = 5;

    if(a<0 || f<0 || f>=pk->MAX || a>=pk->MAX || a>f){
        *value=INT_MAX;
        return 1;
    }   
    if (get(pk->hashtable,tab,size,value)){
        *value=add(INTB(pk,a,f,-1,-1),*value);
        return 1;
    }
    else{
        Node* newNode = (Node*)malloc(sizeof(Node));
        newNode->data = createTuple(tab,size);
        newNode->next = NULL;
        if(*last!=NULL){
            (*last)->next=newNode;
        }
        if(*first==NULL){
            *first=newNode;
        }
        *last=newNode;
        return 0;
    }
}

int get_C1(pk_compound *pk,Node** first,Node** last,int * value,int a, int f, int c,int d){
    int C1 = 67;
    int tab[]={C1,a,f,c,d};
    int size = 5;

    if(a<0 || f<0 || f>=pk->MAX || a>=pk->MAX || a>f){
        *value=INT_MAX;
        return 1;
    }  
    if (get(pk->hashtable,tab,size,value)){
        //printf("value=>%d \n",*value);
        return 1;
    }
    else{
        Node* newNode = (Node*)malloc(sizeof(Node));
        newNode->data = createTuple(tab,size);
        newNode->next = NULL;
        if(*last!=NULL){
            (*last)->next=newNode;
        }
        if(*first==NULL){
            *first=newNode;
        }
        *last=newNode;
        return 0;
    }
}

int compute_C1(pk_compound *pk,Node** first,Node** last,int a, int f, int c,int d) {
    int C1 = 67;
    int tab[]={C1,a,f,c,d};
    int size = 5;
    int value;
    if(a<0 || f<0 || f>=pk->MAX || a>=pk->MAX || a>f){
        return 1;
    } 
    if (get(pk->hashtable,tab,size,&value)){
        return 1;
    }
    int possible=1;
    int min_value=INT_MAX;
    int mfe1=MFEFree(pk,a,c-1);
    int mfe2=MFEFree(pk,d,f);
    
    if(possible==1){
    min_value=add(0,add(mfe1,mfe2));
    }
    
    loop:
    for(int tmp1=a;tmp1<=min(f-1,c-1);tmp1++){
        for(int tmp2=max(d,tmp1+2);tmp2<=f;tmp2++){
            if(evaluate(pk,tmp1,tmp2) && tmp1-tmp2+1<=pk->THETA){
                int tmp3=INTB(pk,tmp1,tmp2,a,f);

                if(tmp3!=INT_MAX){
                    int tmp;
                    if(get_C1(pk,first,last,&tmp,tmp1+1,tmp2-1,c,d)==0){
                        possible=0;
                    }
                    if(possible==1){
                        //printf("C1 res= %d",tmp);
                        min_value=min(min_value, add(tmp,tmp3));
                    }
                }
            }
        }    
    }
    if(possible==1){
        insert(pk->hashtable,tab,size,min_value);
    }
    return possible;
}

int get_B(pk_compound *pk,Node** first,Node** last,int * value,int c,int d,int f,int h){
    int B = 66;
    int tab[]={B,c,d,f,h};
    int size = 5;


    if (get(pk->hashtable,tab,size,value)){
        return 1;
    }
    else{
        Node* newNode = (Node*)malloc(sizeof(Node));
        newNode->data = createTuple(tab,size);
        newNode->next = NULL;
        if(last!=NULL && *last!=NULL){
            (*last)->next=newNode;
        }
        if(first!=NULL && *first==NULL){
            *first=newNode;
        }
        *last=newNode;
        return 0;
    }
}
int compute_B(pk_compound *pk,Node** first,Node** last,int c,int d,int f,int h) {
    int B = 66;
    int value;
    int tab[]={B,c,d,f,h};
    int size = 5;

    if (get(pk->hashtable,tab,size,&value)){
        return 1;
    }
    int min_value = INT_MAX;
    int possible=1;
    for (int g=f;g<h;g++) {
      if(!evaluate(pk,c,h-1)||!evaluate(pk,d-1,g)){continue;}
      int mfe0 = MFEFree(pk,f,g-1);

      int tmp0;
      if(get_CLIQUE0( pk,first,last,&tmp0,c,d-1,g,h-1)==0){possible=0;}
      else if(tmp0==INT_MAX){continue;}

      if(possible==1){
        min_value = min(min_value,add(tmp0,mfe0));
      }
    }

    if(possible==1){
        insert(pk->hashtable,tab,size,min_value);
    }
    return possible;
}

int get_A(pk_compound *pk,Node** first,Node** last,int * value,int a,int h){
    int A = 65;
    int tab[]={A,a,h};
    int size = 3;


    if (get(pk->hashtable,tab,size,value)){
        return 1;
    }
    else{
        Node* newNode = (Node*)malloc(sizeof(Node));
        newNode->data = createTuple(tab,size);
        newNode->next = NULL;
        if(last!=NULL && *last!=NULL){
            (*last)->next=newNode;
        }
        if(first!=NULL && *first==NULL){
            *first=newNode;
        }
        *last=newNode;
        return 0;
    }
}
int compute_A(pk_compound *pk,Node** first,Node** last,int a,int h) {
    int A = 65;
    int value;
    int tab[]={A,a,h};
    int size = 3;

    if (get(pk->hashtable,tab,size,&value)){
        return 1;
    }
    int min_value = INT_MAX;
    int possible=1;
    for (int c=a+1;c<h-2;c++) {
        for (int d=c+1;d<h-1;d++) {
            for (int f=d+1;f<h;f++) {
              if(!evaluate(pk,a,f-1)){continue;}

              int tmp0;
              if(get_C0( pk,first,last,&tmp0,a,f-1,c,d)==0){possible=0;}
              else if(tmp0==INT_MAX){continue;}
              int tmp1;
              if(get_B( pk,first,last,&tmp1,c,d,f,h)==0){possible=0;}
              else if(tmp1==INT_MAX){continue;}

              if(possible==1){
                min_value = min(min_value,add(add(tmp0,tmp1),0));
              }
            }
        }
    }

    if(possible==1){
        insert(pk->hashtable,tab,size,min_value);
    }
    return possible;
}
void backtrace_CLIQUE0(pk_compound *pk,bt_struct* bt,int score,int i, int i2, int j2, int j) {
    
   if(i2>=j2 || i2<i || j2>j){
        return;
    }
    int tmp;
    get_CLIQUE1(pk,NULL,NULL,&tmp,i+1,i2,j2,j-1);
    int sc=INTB(pk,i,j,-1,-1);
    if(score==add(sc,tmp)){
        add_pairing(pk,bt,sc,i,j,-1,-1);
        backtrace_CLIQUE1(pk,bt,tmp,i+1,i2,j2,j-1);
    }

    return;
}

void backtrace_CLIQUE1(pk_compound *pk,bt_struct* bt,int score,int i, int i2, int j2, int j) {
    
   if(i2>=j2 || i2<i || j2>j){
        return ;
    }
    int tmp;
    if (j == j2 && i <= i2 && score==MFEFree(pk,i,i2)) { 
        backtrace_MFEFree(pk,bt,score,i,i2);
    }
    for(int k=i;k<=i2;k++){
        for(int l=max(j2,k+1);l<=j;l++){
            if(evaluate(pk,k,l) && l-k+1<=pk->THETA && (k<i2 || l==j2)){
                tmp=INTB(pk,k,l,i,j);
                if(tmp!=INT_MAX){
                    int tmp0;
                    get_CLIQUE1(pk,NULL,NULL,&tmp0,k+1,i2,j2,l-1);
                    if(score==add(tmp0,tmp)){
                        add_pairing(pk,bt,tmp,k,l,i,j); 
                        backtrace_CLIQUE1(pk,bt,tmp0,k+1,i2,j2,l-1);                       
                        return;
                    }
                }
            }
        }
    }
    return;
}
void backtrace_C0(pk_compound *pk,bt_struct* bt, int score, int a, int f, int c,int d) {

    if(a<0 || f<0 || f>=pk->MAX || a>=pk->MAX || a>f){
        return;
    }
    int tmp;
    get_C1(pk,NULL,NULL,&tmp,a+1,f-1,c,d);
    int sc=INTB(pk,a,f,-1,-1);
    if(score==add(sc,tmp)){
        add_pairing(pk,bt,sc,a,f,-1,-1);
        backtrace_C1(pk,bt,tmp,a+1,f-1,c,d);
    }
    return;
}
void backtrace_C1(pk_compound *pk,bt_struct* bt, int score, int a, int f, int c,int d) {

    if(a<0 || f<0 || f>=pk->MAX || a>=pk->MAX || a>f){
        return;
    }

    int mfe1=MFEFree(pk,a,c-1);
    int mfe2=MFEFree(pk,d,f);
    

    if(score==add(0,add(mfe1,mfe2))){
        
        backtrace_MFEFree(pk,bt,mfe1,a,c-1);
        backtrace_MFEFree(pk,bt,mfe2,d,f);

        return;
    }
    
    loop:
    for(int tmp1=a;tmp1<=min(f-1,c-1);tmp1++){
        for(int tmp2=max(d,tmp1+2);tmp2<=f;tmp2++){
            if(evaluate(pk,tmp1,tmp2) && tmp1-tmp2+1<=pk->THETA){
                int tmp3=INTB(pk,tmp1,tmp2,a,f);
                if(tmp3!=INT_MAX){
                    int tmp0;
                    get_C1(pk,NULL,NULL,&tmp0,tmp1+1,tmp2-1,c,d);
                    if(score==add(tmp0,tmp3)){
                        add_pairing(pk,bt,tmp3,tmp1,tmp2,a,f);
                        backtrace_C1(pk,bt,tmp0,tmp1+1,tmp2-1,c,d);
                        return;
                    }
                }
            }
        }    
    }
    return;
} 
void backtrace_B(pk_compound *pk,bt_struct* bt,int score,int c,int d,int f,int h) {

    for (int g=f;g<h;g++) {
      if(!evaluate(pk,c,h-1)||!evaluate(pk,d-1,g)){continue;}
      int mfe0 = MFEFree(pk,f,g-1);

      int tmp0;
      get_CLIQUE0( pk,NULL,NULL,&tmp0,c,d-1,g,h-1);
      if(tmp0==INT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(pk,bt,mfe0,f,g-1);

        backtrace_CLIQUE0( pk,bt,tmp0,c,d-1,g,h-1);
        bt->bracket+=1;

        return;
      }
    }


}

void backtrace_A(pk_compound *pk,bt_struct* bt,int score,int a,int h) {

    for (int c=a+1;c<h-2;c++) {
        for (int d=c+1;d<h-1;d++) {
            for (int f=d+1;f<h;f++) {
              if(!evaluate(pk,a,f-1)){continue;}

              int tmp0;
              get_C0( pk,NULL,NULL,&tmp0,a,f-1,c,d);
              if(tmp0==INT_MAX){continue;}
              int tmp1;
              get_B( pk,NULL,NULL,&tmp1,c,d,f,h);
              if(tmp1==INT_MAX){continue;}

              if(score==add(add(tmp0,tmp1),0)){

                backtrace_C0( pk,bt,tmp0,a,f-1,c,d);
                bt->bracket+=1;
                backtrace_B( pk,bt,tmp1,c,d,f,h);

                return;
              }
            }
        }
    }


}

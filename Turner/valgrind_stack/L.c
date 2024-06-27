#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/loop_energies.h>

//HASHING TABLE FUNCTION
#define LOAD_FACTOR_THRESHOLD 0.7
#define TUPLE_SIZE 5

int TABLE_SIZE = 64019;
int N = 10;
int THETA_PAIRING=1;
int THETA = 100;

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

char * line = NULL;
int MAX;
char * correct_score = NULL;
char * structure=NULL;
char bracket='A';
vrna_fold_compound_t *fc ;

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

int write_structure(vrna_bp_stack_t *bp,int length){
    int i,j,temp;
    for (int k = 1; k <= bp[0].i; k++) {
      i = bp[k].i;
      j = bp[k].j;
      if (i > length)
        i -= length;

      if (j > length)
        j -= length;

      if (i > j) {
        temp  = i;
        i     = j;
        j     = temp;
      }

      if (i == j) {
        /* Gquad bonds are marked as bp[i].i == bp[i].j */
        structure[i - 1] = '+';
      } else {
        /* the following ones are regular base pairs */
        structure[i - 1]  = '(';
        structure[j - 1]  = ')';
      }
    }
}

int MFEFree(int a,int b){
    if(a>b){
        return 0;
    }
    a=a++;
    b=b++;
    int * indx = fc->jindx;
    int ij = indx[b] + a;
    int * fML =fc->matrices->fML;
    return min(0,fML[ij]);
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

void backtrace_MFEFree(int score, int a,int b){
    if(a>b || score==0){
        return;
    }
    a++;
    b++;
    vrna_bp_stack_t * bp=(vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * ( 1 + (b-a) / 2)));
    int bt=compute_BT(fc,a,b,bp);
    write_structure(bp,strlen(line));
    free(bp);
}

int INTB(int a,int b,int c,int d){
    if(d==-1 && c==-1){
        return 0;
    }    
    a++;b++;

    int energy_full = vrna_eval_int_loop(fc,c,d,a,b);
    
    return energy_full;
}

void backtrace_INTB(int score,int a,int b,int c,int d){
    structure[a]=bracket;
    structure[b]=bracket+32;
    return;
}



typedef struct Node {
    Tuple* data;
    struct Node* next;
} Node;

void freeNode(Node * node){
    free(node->data);
    free(node);
}
// Define the Stack structure
typedef struct {
    Node* top;
} Stack;

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

//declarations
int get_CLIQUE1(HashTable *hashTable, Node **first, Node **last, int *value, int i, int i2, int j2, int j);

int get_CLIQUE0(HashTable *hashTable, Node **first, Node **last, int *value, int i, int i2, int j2, int j);

int compute_CLIQUE1(HashTable *hashTable,Node **first,Node** last, int i, int j, int k, int l);

void backtrace_CLIQUE0(HashTable *hashTable, int score, int i,int j, int k, int l);

void backtrace_CLIQUE1(HashTable *hashTable, int score, int i,int j, int k, int l);
int get_G(HashTable *hashTable,Node** first,Node** last,int * value,int d,int f,int k,int l);

int compute_G(HashTable *hashTable,Node** first,Node** last,int d,int f,int k,int l);

void backtrace_G(HashTable *hashTable,int score,int d,int f,int k,int l) ;

int get_F(HashTable *hashTable,Node** first,Node** last,int * value,int c,int d,int h,int j);

int compute_F(HashTable *hashTable,Node** first,Node** last,int c,int d,int h,int j);

void backtrace_F(HashTable *hashTable,int score,int c,int d,int h,int j) ;

int get_E(HashTable *hashTable,Node** first,Node** last,int * value,int c,int d,int h,int k);

int compute_E(HashTable *hashTable,Node** first,Node** last,int c,int d,int h,int k);

void backtrace_E(HashTable *hashTable,int score,int c,int d,int h,int k) ;

int get_D(HashTable *hashTable,Node** first,Node** last,int * value,int b,int d,int h,int k);

int compute_D(HashTable *hashTable,Node** first,Node** last,int b,int d,int h,int k);

void backtrace_D(HashTable *hashTable,int score,int b,int d,int h,int k) ;

int get_C(HashTable *hashTable,Node** first,Node** last,int * value,int a,int d,int g,int k);

int compute_C(HashTable *hashTable,Node** first,Node** last,int a,int d,int g,int k);

void backtrace_C(HashTable *hashTable,int score,int a,int d,int g,int k) ;

int get_B(HashTable *hashTable,Node** first,Node** last,int * value,int a,int d,int f,int k);

int compute_B(HashTable *hashTable,Node** first,Node** last,int a,int d,int f,int k);

void backtrace_B(HashTable *hashTable,int score,int a,int d,int f,int k) ;

int get_A(HashTable *hashTable,Node** first,Node** last,int * value,int a,int l);

int compute_A(HashTable *hashTable,Node** first,Node** last,int a,int l);

void backtrace_A(HashTable *hashTable,int score,int a,int l) ;


int init_root(HashTable* HashTable,Stack * stack){  
    int A = 65;
    int size = 3;
    for(int l=0+5;l<=MAX;l++){
        for(int a=0;a<l-5;a++){
            int tab[]={A,a,l};
            Node* newNode = (Node*)malloc(sizeof(Node));
            newNode->data = createTuple(tab,size);
            newNode->next = NULL;
            if(stack->top!=NULL){
                newNode->next=stack->top;
                stack->top=newNode;
            }
            stack->top=newNode;
        }
    }
}

int compute_root(HashTable * hashTable,Stack *stack) {
    int min_value=INT_MAX; 
    for(int l=0+5;l<=MAX;l++){
        for(int a=0;a<l-5;a++){
            int mfe1=MFEFree(0,a-1);
            int mfe2=MFEFree(l,MAX-1);
            int mfe=add(mfe1,mfe2);
            int tmp;
            get_A(hashTable,NULL,NULL,&tmp,a,l);
            min_value = min(min_value,add(tmp,mfe));
        }
    }
    return min_value;
}

int fold(HashTable* hasTable, Stack * stack){
    
    init_root(hasTable,stack);
    Node* topNode;
    Node* nextNode;
    Node* add_on_last;
    Node* add_on_first;
    while (!isStackEmpty(stack)) {
        //pop(stack, &topTuple)
        //printStack(stack);
        topNode=stack->top;
        nextNode=topNode->next;
        Tuple * topTuple = topNode->data;
        add_on_last=NULL;
        add_on_first=NULL;
        int a,b,c,d,e,f,g,h,i,j,k,l;
        int i1,i2,j2,j1;
        //printf("here\n");
        //print_tuple(topTuple);
        switch(topTuple->values[0]){
            default: 
            
                 i1=topTuple->values[1];
                 i2=topTuple->values[2];
                 j2=topTuple->values[3];
                 j1=topTuple->values[4];
                if(compute_CLIQUE1(hasTable,&add_on_first,&add_on_last,i1,i2,j2,j1)==1){
                    
                    freeNode(topNode);
                    stack->top=nextNode;
                }
                else{
                    stack->top=add_on_first;
                    add_on_last->next=topNode;
                }
            break;//Clique

             case 71:
                 d=topTuple->values[1];
                 f=topTuple->values[2];
                 k=topTuple->values[3];
                 l=topTuple->values[4];

                 if(compute_G(hasTable,&add_on_first,&add_on_last,d,f,k,l)==1){

                    freeNode(topNode);
                    stack->top=nextNode;
                 }
                 else{
                    stack->top=add_on_first;
                    add_on_last->next=topNode;
                 }
                 break;

             case 70:
                 c=topTuple->values[1];
                 d=topTuple->values[2];
                 h=topTuple->values[3];
                 j=topTuple->values[4];

                 if(compute_F(hasTable,&add_on_first,&add_on_last,c,d,h,j)==1){

                    freeNode(topNode);
                    stack->top=nextNode;
                 }
                 else{
                    stack->top=add_on_first;
                    add_on_last->next=topNode;
                 }
                 break;

             case 69:
                 c=topTuple->values[1];
                 d=topTuple->values[2];
                 h=topTuple->values[3];
                 k=topTuple->values[4];

                 if(compute_E(hasTable,&add_on_first,&add_on_last,c,d,h,k)==1){

                    freeNode(topNode);
                    stack->top=nextNode;
                 }
                 else{
                    stack->top=add_on_first;
                    add_on_last->next=topNode;
                 }
                 break;

             case 68:
                 b=topTuple->values[1];
                 d=topTuple->values[2];
                 h=topTuple->values[3];
                 k=topTuple->values[4];

                 if(compute_D(hasTable,&add_on_first,&add_on_last,b,d,h,k)==1){

                    freeNode(topNode);
                    stack->top=nextNode;
                 }
                 else{
                    stack->top=add_on_first;
                    add_on_last->next=topNode;
                 }
                 break;

             case 67:
                 a=topTuple->values[1];
                 d=topTuple->values[2];
                 g=topTuple->values[3];
                 k=topTuple->values[4];

                 if(compute_C(hasTable,&add_on_first,&add_on_last,a,d,g,k)==1){

                    freeNode(topNode);
                    stack->top=nextNode;
                 }
                 else{
                    stack->top=add_on_first;
                    add_on_last->next=topNode;
                 }
                 break;

             case 66:
                 a=topTuple->values[1];
                 d=topTuple->values[2];
                 f=topTuple->values[3];
                 k=topTuple->values[4];

                 if(compute_B(hasTable,&add_on_first,&add_on_last,a,d,f,k)==1){

                    freeNode(topNode);
                    stack->top=nextNode;
                 }
                 else{
                    stack->top=add_on_first;
                    add_on_last->next=topNode;
                 }
                 break;

             case 65:
                 a=topTuple->values[1];
                 l=topTuple->values[2];

                 if(compute_A(hasTable,&add_on_first,&add_on_last,a,l)==1){

                    freeNode(topNode);
                    stack->top=nextNode;
                 }
                 else{
                    stack->top=add_on_first;
                    add_on_last->next=topNode;
                 }
                 break;

        }
    }
    //printStack(stack);
    return compute_root(hasTable,stack);
}
void backtrace(HashTable * hashTable,Stack* stack,int score) {
    //int a=0;int h=MAX;
    for(int l=0+5;l<=MAX;l++){
        for(int a=0;a<l-5;a++){
            int mfe1=MFEFree(0,a-1);
            int mfe2=MFEFree(l,MAX-1);
            int mfe=add(mfe1,mfe2);
            int tmp0;
            get_A(hashTable,NULL,NULL,&tmp0,a,l);
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
        
        int b_len=0;

        int len=0;

        //reading of the sequence
        len=getline(&line,&b_len, fp);
        
        if(len<=0){
            
            free(line);
            free(correct_score);
            printf("End of file, %d tested",nb_tests);
            break;
        }
        if(line[0]=='#'){
            printf("%s",line);
            continue;
        }
        HashTable *hashTable = createHashTable();
        Stack *stack = malloc(sizeof(Stack));
        initializeStack(stack);
        
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
        if (len != -1) {
            if (ss[len - 1] == '\n') {
                ss[len - 1] = '\0'; // Replace newline character with null terminator
                len--; // Decrement the length to exclude the removed '\n'
            }
        }
        else{
            destroyHashTable(hashTable);
            freeStack(stack);
            free(line);
            free(ss);
            printf("No secondary structure found for test %d\n",nb_tests);
            exit(-1);
        }
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
            freeStack(stack);
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

        vrna_init_rand();
        fc = vrna_fold_compound(line, NULL, VRNA_OPTION_DEFAULT);
        vrna_constraints_add(fc, ss, VRNA_CONSTRAINT_DB_DEFAULT);
        vrna_mfe(fc, NULL);
        
        //folding
        int score = fold(hashTable,stack);
        //retrieve backtrack
        backtrace(hashTable,stack,score);
        vrna_fold_compound_free(fc);
        clock_t test_end_time = clock();

        double test_time = ((double)(test_end_time - test_start_time)) / CLOCKS_PER_SEC;
        total_time += test_time;
        //end of test
        
        //sanity checks
        if(score == number || 1){//we autorise wrong score in the case of the Turner model
            float fl=(float)score/100.0;
            printf("Score %.2f ",fl);
            int nb=0;
            for(int i=0;i<MAX;i++){
                
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
            freeStack(stack);
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
        freeStack(stack);
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

int get_CLIQUE0(HashTable *hashTable,Node** first,Node** last,int * value,int i, int i2, int j2, int j){
    int CLIQUE=1;
    int tab[] = {CLIQUE,i+1,i2,j2,j-1};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        *value=INT_MAX;
        return 1;
    }
    if (get(hashTable,tab,size,value)){
        //printf("*value = %d(+%d)\n",*value,INTB(i,j,-1,-1));
        *value=add(INTB(i,j,-1,-1),*value);
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

int get_CLIQUE1(HashTable *hashTable,Node** first,Node** last,int * value,int i, int i2, int j2, int j){
    int CLIQUE=1;
    int tab[] = {CLIQUE,i,i2,j2,j};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        *value=INT_MAX;
        return 1;
    }
    if (get(hashTable,tab,size,value)){
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


int compute_CLIQUE1(HashTable *hashTable,Node** first,Node** last,int i, int i2, int j2, int j) {
    int CLIQUE=1;
    int value;
    int tab[] = {CLIQUE,i,i2,j2,j};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        return 1;
    }

    if (get(hashTable,tab,size,&value)) { 
        return 1;
    }
    int possible=1;
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
                    int tmp2;
                    if(get_CLIQUE1(hashTable,first,last,&tmp2,k+1,i2,j2,l-1)==0){
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
        insert(hashTable,tab,size,min_value);
    }
    return possible;
}

int get_G(HashTable *hashTable,Node** first,Node** last,int * value,int d,int f,int k,int l){
    int G = 71;
    int tab[]={G,d,f,k,l};
    int size = 5;


    if (get(hashTable,tab,size,value)){
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
int compute_G(HashTable *hashTable,Node** first,Node** last,int d,int f,int k,int l) {
    int G = 71;
    int value;
    int tab[]={G,d,f,k,l};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return 1;
    }
    int min_value = INT_MAX;
    int possible=1;
    for (int e=d;e<f;e++) {
      if(!evaluate(e,l-1)||!evaluate(f-1,k)){continue;}
      int mfe0 = MFEFree(d,e-1);

      int tmp0;
      if(get_CLIQUE0( hashTable,first,last,&tmp0,e,f-1,k,l-1)==0){possible=0;}
      else if(tmp0==INT_MAX){continue;}

      if(possible==1){
        min_value = min(min_value,add(tmp0,mfe0));
      }
    }

    if(possible==1){
        insert(hashTable,tab,size,min_value);
    }
    return possible;
}

int get_F(HashTable *hashTable,Node** first,Node** last,int * value,int c,int d,int h,int j){
    int F = 70;
    int tab[]={F,c,d,h,j};
    int size = 5;


    if (get(hashTable,tab,size,value)){
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
int compute_F(HashTable *hashTable,Node** first,Node** last,int c,int d,int h,int j) {
    int F = 70;
    int value;
    int tab[]={F,c,d,h,j};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return 1;
    }
    int min_value = INT_MAX;
    int possible=1;
    for (int i=h;i<j;i++) {
      if(!evaluate(c,j-1)||!evaluate(d-1,i)){continue;}
      int mfe0 = MFEFree(h,i-1);

      int tmp0;
      if(get_CLIQUE0( hashTable,first,last,&tmp0,c,d-1,i,j-1)==0){possible=0;}
      else if(tmp0==INT_MAX){continue;}

      if(possible==1){
        min_value = min(min_value,add(tmp0,mfe0));
      }
    }

    if(possible==1){
        insert(hashTable,tab,size,min_value);
    }
    return possible;
}

int get_E(HashTable *hashTable,Node** first,Node** last,int * value,int c,int d,int h,int k){
    int E = 69;
    int tab[]={E,c,d,h,k};
    int size = 5;


    if (get(hashTable,tab,size,value)){
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
int compute_E(HashTable *hashTable,Node** first,Node** last,int c,int d,int h,int k) {
    int E = 69;
    int value;
    int tab[]={E,c,d,h,k};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return 1;
    }
    int min_value = INT_MAX;
    int possible=1;
    for (int j=h+1;j<k+1;j++) {

      int mfe0 = MFEFree(j,k-1);

      int tmp0;
      if(get_F( hashTable,first,last,&tmp0,c,d,h,j)==0){possible=0;}
      else if(tmp0==INT_MAX){continue;}

      if(possible==1){
        min_value = min(min_value,add(tmp0,mfe0));
      }
    }

    if(possible==1){
        insert(hashTable,tab,size,min_value);
    }
    return possible;
}

int get_D(HashTable *hashTable,Node** first,Node** last,int * value,int b,int d,int h,int k){
    int D = 68;
    int tab[]={D,b,d,h,k};
    int size = 5;


    if (get(hashTable,tab,size,value)){
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
int compute_D(HashTable *hashTable,Node** first,Node** last,int b,int d,int h,int k) {
    int D = 68;
    int value;
    int tab[]={D,b,d,h,k};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return 1;
    }
    int min_value = INT_MAX;
    int possible=1;
    for (int c=b;c<d;c++) {

      int mfe0 = MFEFree(b,c-1);

      int tmp0;
      if(get_E( hashTable,first,last,&tmp0,c,d,h,k)==0){possible=0;}
      else if(tmp0==INT_MAX){continue;}

      if(possible==1){
        min_value = min(min_value,add(tmp0,mfe0));
      }
    }

    if(possible==1){
        insert(hashTable,tab,size,min_value);
    }
    return possible;
}

int get_C(HashTable *hashTable,Node** first,Node** last,int * value,int a,int d,int g,int k){
    int C = 67;
    int tab[]={C,a,d,g,k};
    int size = 5;


    if (get(hashTable,tab,size,value)){
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
int compute_C(HashTable *hashTable,Node** first,Node** last,int a,int d,int g,int k) {
    int C = 67;
    int value;
    int tab[]={C,a,d,g,k};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return 1;
    }
    int min_value = INT_MAX;
    int possible=1;
    for (int b=a+1;b<d;b++) {
        for (int h=g+1;h<k;h++) {
          if(!evaluate(a,h-1)||!evaluate(b-1,g)){continue;}

          int tmp0;
          if(get_D( hashTable,first,last,&tmp0,b,d,h,k)==0){possible=0;}
          else if(tmp0==INT_MAX){continue;}
          int tmp1;
          if(get_CLIQUE0( hashTable,first,last,&tmp1,a,b-1,g,h-1)==0){possible=0;}
          else if(tmp1==INT_MAX){continue;}

          if(possible==1){
            min_value = min(min_value,add(add(tmp0,tmp1),0));
          }
        }
    }

    if(possible==1){
        insert(hashTable,tab,size,min_value);
    }
    return possible;
}

int get_B(HashTable *hashTable,Node** first,Node** last,int * value,int a,int d,int f,int k){
    int B = 66;
    int tab[]={B,a,d,f,k};
    int size = 5;


    if (get(hashTable,tab,size,value)){
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
int compute_B(HashTable *hashTable,Node** first,Node** last,int a,int d,int f,int k) {
    int B = 66;
    int value;
    int tab[]={B,a,d,f,k};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return 1;
    }
    int min_value = INT_MAX;
    int possible=1;
    for (int g=f;g<k-1;g++) {

      int mfe0 = MFEFree(f,g-1);

      int tmp0;
      if(get_C( hashTable,first,last,&tmp0,a,d,g,k)==0){possible=0;}
      else if(tmp0==INT_MAX){continue;}

      if(possible==1){
        min_value = min(min_value,add(tmp0,mfe0));
      }
    }

    if(possible==1){
        insert(hashTable,tab,size,min_value);
    }
    return possible;
}

int get_A(HashTable *hashTable,Node** first,Node** last,int * value,int a,int l){
    int A = 65;
    int tab[]={A,a,l};
    int size = 3;


    if (get(hashTable,tab,size,value)){
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
int compute_A(HashTable *hashTable,Node** first,Node** last,int a,int l) {
    int A = 65;
    int value;
    int tab[]={A,a,l};
    int size = 3;

    if (get(hashTable,tab,size,&value)){
        return 1;
    }
    int min_value = INT_MAX;
    int possible=1;
    for (int d=a+2;d<l-3;d++) {
        for (int f=d+1;f<l-2;f++) {
            for (int k=f+2;k<l;k++) {


              int tmp0;
              if(get_B( hashTable,first,last,&tmp0,a,d,f,k)==0){possible=0;}
              else if(tmp0==INT_MAX){continue;}
              int tmp1;
              if(get_G( hashTable,first,last,&tmp1,d,f,k,l)==0){possible=0;}
              else if(tmp1==INT_MAX){continue;}

              if(possible==1){
                min_value = min(min_value,add(add(tmp0,tmp1),0));
              }
            }
        }
    }

    if(possible==1){
        insert(hashTable,tab,size,min_value);
    }
    return possible;
}
void backtrace_CLIQUE0(HashTable *hashTable,int score,int i, int i2, int j2, int j) {
    
   if(i2>=j2 || i2<i || j2>j){
        return;
    }
    int tmp;
    get_CLIQUE1(hashTable,NULL,NULL,&tmp,i+1,i2,j2,j-1);
    int sc=INTB(i,j,-1,-1);
    if(score==add(sc,tmp)){
        backtrace_INTB(sc,i,j,-1,-1);
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
                    int tmp0;
                    get_CLIQUE1(hashTable,NULL,NULL,&tmp0,k+1,i2,j2,l-1);
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

void backtrace_G(HashTable *hashTable,int score,int d,int f,int k,int l) {

    for (int e=d;e<f;e++) {
      if(!evaluate(e,l-1)||!evaluate(f-1,k)){continue;}
      int mfe0 = MFEFree(d,e-1);

      int tmp0;
      get_CLIQUE0( hashTable,NULL,NULL,&tmp0,e,f-1,k,l-1);
      if(tmp0==INT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,d,e-1);

        backtrace_CLIQUE0( hashTable,tmp0,e,f-1,k,l-1);
        bracket+=1;

        return;
      }
    }


}

void backtrace_F(HashTable *hashTable,int score,int c,int d,int h,int j) {

    for (int i=h;i<j;i++) {
      if(!evaluate(c,j-1)||!evaluate(d-1,i)){continue;}
      int mfe0 = MFEFree(h,i-1);

      int tmp0;
      get_CLIQUE0( hashTable,NULL,NULL,&tmp0,c,d-1,i,j-1);
      if(tmp0==INT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,h,i-1);

        backtrace_CLIQUE0( hashTable,tmp0,c,d-1,i,j-1);
        bracket+=1;

        return;
      }
    }


}

void backtrace_E(HashTable *hashTable,int score,int c,int d,int h,int k) {

    for (int j=h+1;j<k+1;j++) {

      int mfe0 = MFEFree(j,k-1);

      int tmp0;
      get_F( hashTable,NULL,NULL,&tmp0,c,d,h,j);
      if(tmp0==INT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,j,k-1);

        backtrace_F( hashTable,tmp0,c,d,h,j);

        return;
      }
    }


}

void backtrace_D(HashTable *hashTable,int score,int b,int d,int h,int k) {

    for (int c=b;c<d;c++) {

      int mfe0 = MFEFree(b,c-1);

      int tmp0;
      get_E( hashTable,NULL,NULL,&tmp0,c,d,h,k);
      if(tmp0==INT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,b,c-1);

        backtrace_E( hashTable,tmp0,c,d,h,k);

        return;
      }
    }


}

void backtrace_C(HashTable *hashTable,int score,int a,int d,int g,int k) {

    for (int b=a+1;b<d;b++) {
        for (int h=g+1;h<k;h++) {
          if(!evaluate(a,h-1)||!evaluate(b-1,g)){continue;}

          int tmp0;
          get_D( hashTable,NULL,NULL,&tmp0,b,d,h,k);
          if(tmp0==INT_MAX){continue;}
          int tmp1;
          get_CLIQUE0( hashTable,NULL,NULL,&tmp1,a,b-1,g,h-1);
          if(tmp1==INT_MAX){continue;}

          if(score==add(add(tmp0,tmp1),0)){

            backtrace_D( hashTable,tmp0,b,d,h,k);
            backtrace_CLIQUE0( hashTable,tmp1,a,b-1,g,h-1);
            bracket+=1;

            return;
          }
        }
    }


}

void backtrace_B(HashTable *hashTable,int score,int a,int d,int f,int k) {

    for (int g=f;g<k-1;g++) {

      int mfe0 = MFEFree(f,g-1);

      int tmp0;
      get_C( hashTable,NULL,NULL,&tmp0,a,d,g,k);
      if(tmp0==INT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,f,g-1);

        backtrace_C( hashTable,tmp0,a,d,g,k);

        return;
      }
    }


}

void backtrace_A(HashTable *hashTable,int score,int a,int l) {

    for (int d=a+2;d<l-3;d++) {
        for (int f=d+1;f<l-2;f++) {
            for (int k=f+2;k<l;k++) {


              int tmp0;
              get_B( hashTable,NULL,NULL,&tmp0,a,d,f,k);
              if(tmp0==INT_MAX){continue;}
              int tmp1;
              get_G( hashTable,NULL,NULL,&tmp1,d,f,k,l);
              if(tmp1==INT_MAX){continue;}

              if(score==add(add(tmp0,tmp1),0)){

                backtrace_B( hashTable,tmp0,a,d,f,k);
                backtrace_G( hashTable,tmp1,d,f,k,l);

                return;
              }
            }
        }
    }


}
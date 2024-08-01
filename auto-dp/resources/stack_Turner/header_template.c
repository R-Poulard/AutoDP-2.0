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
#include "aux.h"

//HASHING TABLE FUNCTION
#define TUPLE_SIZE_##NAME## ##MODIFY_SIZE_TUPLE##

PRIVATE int TABLE_SIZE = 64019;
PRIVATE int N = 10;


typedef struct Tuple {
    int values[TUPLE_SIZE_##NAME##];
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
typedef struct {
    Node* top;
} Stack;

typedef struct {
    HashTable* hashtable;
    Stack* Stack;
    vrna_fold_compound_t* fc;
    vrna_fold_compound_t* fc_aux;
    const char* line;
    int** matrix;
    int MAX;
    int THETA;
    int THETA_PAIRING;
} pk_compound;

PRIVATE void preset_tablesize(int size) {
    if (size > 0) {
        TABLE_SIZE = size;
    }
    else {
        printf("Size negative or null\n");
    }
}

PRIVATE void preset_N(int n) {
    if (n > 0) {
        N = n;
    }
    else {
        printf("N negative or null\n");
    }
}

PRIVATE Tuple *createTuple(int val[],int size) {
    Tuple *tuple = (Tuple *)malloc(sizeof(Tuple));

    // Initialize all values in tuple to 0
    for (int i = 0; i < TUPLE_SIZE_##NAME##; i++) {
        tuple->values[i] = 0;
    }

    // Copy values from val array, up to its length or TUPLE_SIZE
    for (int i = 0; i < size; i++) {
        tuple->values[i] = val[i];
    }

    return tuple;
}

PRIVATE int compare_tuple(Tuple *tp1, int val[], int size) {
    int i;
    // Compare elements until the end of val or tp1->values
    for (i = 0; i < size ; i++) {
        if (tp1->values[i] != val[i]) {
            return 0; // Tuples are different
        }
    }
    // Check if remaining elements of tp1->values are all zeros
    for (; i < TUPLE_SIZE_##NAME##; i++) {
        if (tp1->values[i] != 0) {
            return 0; // Tuples are different
        }
    }
    return 1; // Tuples are identical
}

// Hash function
PRIVATE int hash(int keys[], int capacity,int size) {
    int key = 0;
    for (int i = 0; i < size; i++) {
        key += keys[i] * (int)pow(N, TUPLE_SIZE_##NAME##-i-1);
    }
    return abs(key % capacity);
}


// Function to initialize a new hash table
PRIVATE HashTable *createHashTable() {
    HashTable *hashTable = (HashTable *)malloc(sizeof(HashTable));
    hashTable->buckets = (HashNode **)calloc(TABLE_SIZE, sizeof(HashNode *));
    hashTable->size = 0;
    hashTable->capacity = TABLE_SIZE;
    return hashTable;
}

// Function to free memory used by the hash table 
PRIVATE void destroyHashTable(HashTable *hashTable) {
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
PRIVATE void insert(HashTable *hashTable, int keys[], int size,int value) {
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
PRIVATE bool get(HashTable *hashTable, int keys[], int size, int *value) {
    
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

PRIVATE void print_tuple(Tuple *tpl){
    printf("tp=");
    for (int i = 0; i < TUPLE_SIZE_##NAME##; i++) {
        printf("%d,",tpl->values[i]);
    }

}

PRIVATE void print_table_stat(HashTable *hashTable){
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

PRIVATE void print_table(HashTable *hashTable){
    
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

PRIVATE int evaluate(pk_compound* pk,int x, int y) {
    int a=abs(x-y)<pk->THETA_PAIRING;
    int b=pk->matrix[x][y];
    int c=matching(pk->line[x],pk->line[y]);

    if(abs(x-y)<pk->THETA_PAIRING || pk->matrix[x][y]==0 || matching(pk->line[x],pk->line[y])==INF){
        return 0;
    }
    return 1;
}

PRIVATE int write_structure(pk_compound* pk,bt_struct* bt,vrna_bp_stack_t *bp){
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

PRIVATE int MFEFree(pk_compound* pk,int a,int b){
    return pk->fc_aux->params->MLbase;
    if(a>b){
        return 0;
    }
    a=a++;
    b=b++;

    int * indx = pk->fc_aux->jindx;
    int ij = indx[b] + a;
    int * fML =pk->fc_aux->matrices->fML;
    int up_energy=pk->fc_aux->params->MLbase;
    return min(up_energy,fML[ij]);
}

PRIVATE int compute_BT(pk_compound *pk, int i,int j, vrna_bp_stack_t *bp_stack){
    int ret;
    int p,q,comp1,comp2;
    int s=0;
    int b=bp_stack[0].i;
    sect bt_stack[500]; /* stack of partial structures for backtracking */

    int BT=vrna_BT_mb_loop_split(pk->fc_aux, &i, &j, &p, &q, &comp1, &comp2, bp_stack, &b);
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
    BT=vrna_backtrack_from_intervals(pk->fc_aux,bp_stack,bt_stack,s);
    return BT;
}

PRIVATE int backtrace_MFEFree(pk_compound*pk,bt_struct* bt,int score, int a,int b){
    int up_energy=pk->fc_aux->params->MLbase;
    if(a>b || score==up_energy){
        return 1;
    }
    a++;
    b++;
    vrna_bp_stack_t * bp=(vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * ( 1 + (b-a) / 2)));
    int cmp_bt=compute_BT(pk,a,b,bp);
    write_structure(pk,bt,bp);
    free(bp);
    return 1;
}

PRIVATE int INTB(pk_compound* pk,int a,int b,int c,int d){
    
    if(d==-1 && c==-1){
        return 0;
    } 
    int * fML =pk->fc_aux->matrices->fML;
    int * indx = pk->fc_aux->jindx;   

    a++;b++;

    int right_index = indx[d] + b;
    int right_ml=fML[right_index];
    int left_index = indx[c] + a;
    int left_ml=fML[left_index];

    left_ml=-5;
    right_ml=-5;


    int open=pk->fc_aux->params->MLclosing;//penality for creating multiloop
    
    //case C (multiloop with helix on each side)
    
    int free_energy= open+pk->fc_aux->params->MLintern[3]+right_ml+left_ml;

    //case A (interior loop)
    if(c-a+1<pk->THETA && d-b+1<pk->THETA){
        
         free_energy = min(free_energy,vrna_eval_int_loop(pk->fc_aux,c,d,a,b));
        
    }
    //case B (multiloop only on left loop)
    if(d-b+1<pk->THETA){
        int unpaired=pk->fc_aux->params->MLbase;
        free_energy = min(free_energy,unpaired+open+pk->fc_aux->params->MLintern[2]+left_ml);

    }
    //case B (multiloop loop)
    if(c-a+1<pk->THETA){
        int unpaired=pk->fc_aux->params->MLbase;
        free_energy = min(free_energy,unpaired+open+pk->fc_aux->params->MLintern[2]+right_ml);
        
    }
    

    return free_energy;
}

PRIVATE void add_pairing(pk_compound* pk,bt_struct* bt, int score,int a,int b,int c,int d){
    bt->structure[a]=bt->bracket;
    bt->structure[b]=tolower(bt->bracket);
    return;
}

PRIVATE int backtrace_INTB(pk_compound* pk,bt_struct* bt, int score,int a,int b,int c,int d){
    
   
    add_pairing(pk,bt,score,a,b,c,d);
    //return 1;
    if(d==-1 && c==-1){
        return 1;
    } 
    int * fML =pk->fc_aux->matrices->fML;
    int * indx = pk->fc_aux->jindx;   
    
    a++;b++;
    
    int right_index = indx[d] + b;
    int right_ml=fML[right_index];
    int left_index = indx[c] + a;
    int left_ml=fML[left_index];
    
    left_ml=-5;
    right_ml=-5;
 

    int open=pk->fc_aux->params->MLclosing;//penality for creating multiloop

    //case C (multiloop with helix on each side)
    
    if(score==open+pk->fc_aux->params->MLintern[3]+right_ml+left_ml){
        
        int res1= backtrace_MFEFree(pk,bt,left_ml,a--,c--);
        int res2= backtrace_MFEFree(pk,bt,right_ml,b--,d--);
        return res1 && res2;
    }
    
    //case A (interior loop)
    if(c-a+1<pk->THETA && d-b+1<pk->THETA){
        
        if(score==vrna_eval_int_loop(pk->fc,c,d,a,b)){
            
            return 1;
        }
    }
    //case B (multiloop only on left loop)
    if(d-b+1<pk->THETA){
        int unpaired=pk->fc_aux->params->MLbase;
        if(score==unpaired+open+pk->fc_aux->params->MLintern[2]+left_ml){
            
            return backtrace_MFEFree(pk,bt,left_ml,a--,c--);
        }
    }
    //case B (multiloop loop)
    if(c-a+1<pk->THETA){
        int unpaired=pk->fc_aux->params->MLbase;
        if(score==unpaired+open+pk->fc_aux->params->MLintern[2]+right_ml){
            
            return backtrace_MFEFree(pk,bt,right_ml,b--,d--);
        }
    }

    return 0;
}



PRIVATE void freeNode(Node * node){
    free(node->data);
    free(node);
}

// Function to initialize the stack
PRIVATE void initializeStack(Stack* stack) {
    stack->top = NULL;
}

// Function to check if the stack is empty
PRIVATE bool isStackEmpty(Stack* stack) {
    return stack->top == NULL;
}

// Function to push a tuple onto the stack
PRIVATE void push(Stack* stack,int tab[],int size) {
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
PRIVATE bool pop(Stack* stack, Tuple* item) {
    if (isStackEmpty(stack)) {
        
        return false;
    }
    Node* temp = stack->top;
    item = temp->data;
    stack->top = stack->top->next;
    free(temp);
    return true;
}


// Function to free the stack
PRIVATE void freeStack(Stack* stack) {
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
PRIVATE void printStack(Stack* stack) {
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


PRIVATE pk_compound* create_pk_compound(vrna_fold_compound_t *fc,vrna_fold_compound_t *fc_aux, const char* seq,const char* ss){
    pk_compound *pk=malloc(sizeof(pk_compound));
    pk->fc=fc;
    pk->fc_aux=fc_aux;
    pk->hashtable=createHashTable();
    pk->Stack=malloc(sizeof(Stack));
    initializeStack(pk->Stack);
    pk->line=seq;
    pk->MAX=strlen(seq)-1;
    createMatrix(&(pk->matrix),pk->MAX,ss);
    preset_N(pk->MAX);
    pk->THETA=MAXLOOP;
    pk->THETA_PAIRING=TURN;
}
PRIVATE void free_pk(pk_compound* pk){
    destroyHashTable(pk->hashtable);
    freeMatrix(pk->matrix,pk->MAX);
    freeStack(pk->Stack);
    free(pk);
}

##DECLARATIONS##

PRIVATE int init_root(pk_compound* pk){  
    int A = 65;
    int size = 3;
    for(int END=0+ANCH;END<=pk->MAX;END++){
        for(int a=0;a<END-ANCH;a++){
            int tab[]={A,a,END};
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

PRIVATE int init_root_inter(pk_compound* pk,int a,int h){  
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

PRIVATE int compute_root(pk_compound* pk) {
    int min_value=INF; 
    for(int END=0+ANCH;END<=pk->MAX;END++){
        for(int a=0;a<END-ANCH;a++){
            int mfe1=MFEFree(pk,0,a-1);
            int mfe2=MFEFree(pk,END,pk->MAX-1);
            int mfe=add(mfe1,mfe2);
            int tmp;
            get_A(pk,NULL,NULL,&tmp,a,END);
            min_value = min(min_value,add(tmp,mfe));
        }
    }
    return min_value;
}


PRIVATE int fold(pk_compound* pk){
    
    init_root(pk);
    Node* topNode;
    Node* nextNode;
    Node* add_on_last;
    Node* add_on_first;
    while (!isStackEmpty(pk->Stack)) {
        topNode=pk->Stack->top;
        nextNode=topNode->next;
        Tuple * topTuple = topNode->data;
        add_on_last=NULL;
        add_on_first=NULL;
        int ALL_ANC;
        int i1,i2,j2,j1;
        switch(topTuple->values[0]){
VARIABLE_PART
        }
    }
    return compute_root(pk);
}
PRIVATE void backtrace(pk_compound *pk,bt_struct* bt,int score) {
    
    for(int END=0+ANCH;END<=pk->MAX;END++){
        for(int a=0;a<END-ANCH;a++){
            int mfe1=MFEFree(pk,0,a-1);
            int mfe2=MFEFree(pk,END,pk->MAX-1);
            int mfe=add(mfe1,mfe2);
            int tmp0;
            get_A(pk,NULL,NULL,&tmp0,a,END);
            if(score==add(mfe,tmp0)){
                backtrace_A(pk,bt,tmp0,a,END);
                backtrace_MFEFree(pk,bt,mfe1,0,a-1);
                backtrace_MFEFree(pk,bt,mfe2,END,pk->MAX-1);
                return;
            }
       }
    }
}

PRIVATE int backtrace_on_interval(pk_compound *pk,bt_struct* bt,int score,int a,int END) {
    int tmp0;
    get_A(pk,NULL,NULL,&tmp0,a,END);
    
    if(score==tmp0){
        if(score==INF){
            return 1;
        }
        int res1=backtrace_A(pk,bt,tmp0,a,END);
        
        return res1;
    }
    
    return 0;
}

PRIVATE bt_struct *bt;

PRIVATE int scoring(vrna_fold_compound_t *fc,
                                         int i,
                                         int j,
                                         void *data) {
  pk_compound *pk=(pk_compound*)data;
  int tmp;
  Node* first;
  Node* last;
  if(get_A(pk,&first,&last,&tmp,i-1,j-1)==0){
    init_root_inter(pk,i-1,j-1);
    fold(pk);
    get_A(pk,&first,&last,&tmp,i-1,j-1);
  }
  if(tmp==INF){
    return INF+1;
  }
  return tmp;
}

PRIVATE int scoring_f(vrna_fold_compound_t *fc,
                                         int i,
                                         int j,
                                         void *data) {
  pk_compound *pk=(pk_compound*)data;
  int min=INF;
  int free_energy;
  Node* first;
  Node* last;
  int tmp;
  
  for(int tmpi=0;tmpi<j-1;tmpi++){
    
    if(get_A(pk,&first,&last,&tmp,tmpi,j-1)==0){
      init_root_inter(pk,tmpi,j-1);
      fold(pk);
      get_A(pk,&first,&last,&tmp,tmpi,j-1);
    }
    free_energy=MFEFree(pk,0,tmpi-1);
    
    if(min>add(tmp,free_energy)){
      min=add(tmp,free_energy);
    }
  }
  if(min==INF){
    
    return INF+1;
  }
  
  return min;
}

//requiered beceause i and j are forming the base pair
PRIVATE int scoring_c(vrna_fold_compound_t *fc,
                                         int i,
                                         int j,
                                         void *data) {
  pk_compound *pk=(pk_compound*)data;
  int tmp;
  Node* first;
  Node* last;
  if(get_A(pk,&first,&last,&tmp,i,j-2)==0){
    init_root_inter(pk,i,j-2);
    fold(pk);
    get_A(pk,&first,&last,&tmp,i,j-2);
  }
  if(tmp==INF){
    return INF+1;
  }
  return tmp;
}


PRIVATE int custom_BT_c(vrna_fold_compound_t *gc,unsigned int ui,unsigned int uj,vrna_bps_t bp_stack,vrna_bts_t bt_stack,void *data){
  
  
  pk_compound *pk=(pk_compound*)data;
  int tmp;
  Node* first;
  Node* last;
  bt->bt_stack=bt_stack;
  bt->bp_stack=bp_stack;
  int i=(int)ui;
  int j=(int)uj;
  
  get_A(pk,&first,&last,&tmp,i,j-2);
  
  if(tmp==gc->matrices->c[gc->jindx[j] + i]){
    return backtrace_on_interval(pk,bt,tmp,i,j-2);
  }
  else{
    
    return 0;
  }
}


PRIVATE int custom_BT_fML(vrna_fold_compound_t *gc,unsigned int ui,unsigned int uj,vrna_bps_t bp_stack,vrna_bts_t bt_stack,void *data){
  pk_compound *pk=(pk_compound*)data;
  int tmp;
  Node* first;
  Node* last;
  bt->bt_stack=bt_stack;
  bt->bp_stack=bp_stack;
  int i=(int)ui;
  int j=(int)uj;
  get_A(pk,&first,&last,&tmp,i-1,j-1);
  if(tmp==gc->matrices->fML[gc->jindx[j] + i]){
  return backtrace_on_interval(pk,bt,tmp,i-1,j-1);
  }
  else{
    return 0;
  }
}

PRIVATE int custom_BT_fM1(vrna_fold_compound_t *gc,unsigned int ui,unsigned int uj,vrna_bps_t bp_stack,vrna_bts_t bt_stack,void *data){
  pk_compound *pk=(pk_compound*)data;
  int tmp;
  Node* first;
  Node* last;
  bt->bt_stack=bt_stack;
  bt->bp_stack=bp_stack;
  int i=(int)ui;
  int j=(int)uj;
  get_A(pk,&first,&last,&tmp,i-1,j-1);
  if(tmp==gc->matrices->fM1[gc->jindx[j] + i]){
  return backtrace_on_interval(pk,bt,tmp,i-1,j-1);
  }
  else{
    return 0;
  }
}

PRIVATE int custom_BT_f(vrna_fold_compound_t *gc,unsigned int ui,unsigned int uj,vrna_bps_t bp_stack,vrna_bts_t bt_stack,void *data){
  
  pk_compound *pk=(pk_compound*)data;
  int tmp;
  Node* first;
  Node* last;
  
  int min_free_energy=INF;
  int min_tmp=INF;
  int free_energy;

  int best_i=0;
  int j=(int) uj;
  int right_score=gc->matrices->f5[j];

  for(;j>=2;j--){
    for(int i=0;i<j-1;i++){
        get_A(pk,&first,&last,&tmp,i,j-1);
        free_energy=MFEFree(pk,0,i-1);
        
        if(right_score==add(tmp,free_energy)){
        best_i=i;
        min_free_energy=free_energy;
        min_tmp=tmp;
        bt->bt_stack=bt_stack;
        bt->bp_stack=bp_stack;
        int res1=backtrace_on_interval(pk,bt,min_tmp,best_i,j-1);
        
        int res2=backtrace_MFEFree(pk,bt,min_free_energy,0,best_i-1);
        return res1 && res2;
        }
    }
  }
  return 0;
}

PRIVATE void free_aux(void * data){
  pk_compound *pk=(pk_compound*)data;
  free_pk(data);
}

int setup_##NAME##(vrna_fold_compound_t *fc,vrna_fold_compound_t *fc_aux,const char* seq, const char* con,bt_struct* BT){
    bt=BT;
    pk_compound * pk=create_pk_compound(fc,fc_aux,seq,con);

    vrna_gr_add_aux_f(fc,scoring_f,custom_BT_f,pk,NULL,free_aux);
    vrna_gr_add_aux_c(fc,scoring_c,custom_BT_c,pk,NULL,NULL);
    vrna_gr_add_aux_m(fc,scoring,custom_BT_fML,pk,NULL,NULL);
    vrna_gr_add_aux_m1(fc,scoring,custom_BT_fM1,pk,NULL,NULL);
}
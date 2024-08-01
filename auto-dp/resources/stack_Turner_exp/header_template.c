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
    double value;
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
PRIVATE void insert(HashTable *hashTable, int keys[], int size,double value) {
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
PRIVATE bool get(HashTable *hashTable, int keys[], int size, double *value) {
    
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

PRIVATE double MFEFree(pk_compound* pk,int a,int b){
    if(a>b){
        return 1.0;
    }
    return 1.0;
    a=a++;
    b=b++;

    int * indx = pk->fc->jindx;
    int ij = indx[b] + a;
    int * fML =pk->fc->matrices->fML;
    int up_energy=pk->fc->params->MLbase*(b-a+1);
    return min(up_energy,fML[ij]);
}


PRIVATE int compute_BT(vrna_fold_compound_t *fc, int i,int j,bt_struct* bt, vrna_bp_stack_t* bp_stack){
    int ret;
    int p,q,comp1,comp2;
    int s=0;
    int b=0;
    int BT=1;
    if(bt->bp_stack==NULL){
        b=0;  
        vrna_bps_t bp=vrna_bps_init((4 * ( 1 + (j-i) / 2)));
        BT=vrna_bt_mb_loop_split(fc, &i, &j, &p, &q, &comp1, &comp2, bp);  
    }
    else{
        b=vrna_bps_at(bt->bp_stack,0).i;
        BT=vrna_bt_mb_loop_split(fc, &i, &j, &p, &q, &comp1, &comp2, bt->bp_stack);
    }
    /* stack of partial structures for backtracking */
    //printf("a=%d,b=%d\n",i,j);
    
    if(BT!=1){
        return BT;
    }
    if(bt->bt_stack==NULL){
        sect bt_stack[500];
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
    }
    else{
       if (i > 0) {
            vrna_sect_t *vt=(vrna_sect_t*)vrna_alloc(sizeof(vrna_sect_t)) ;
            vt->i=i;
            vt->j=j;
            vt->ml=comp1; 
            vrna_bts_push(bt->bt_stack,*vt);
        }

        if (p > 0) {
            vrna_sect_t *vt2=(vrna_sect_t*)vrna_alloc(sizeof(vrna_sect_t)) ;
            vt2->i=p;
            vt2->j=q;
            vt2->ml=comp2; 
            vrna_bts_push(bt->bt_stack,*vt2);
        }  
    }
    return BT;
}

PRIVATE int backtrace_MFEFree(pk_compound* pk,bt_struct* bt,int score, int a,int b){

    int up_energy=pk->fc->params->MLbase*(b-a+1);
    if(a>b || score==up_energy){
        return 1;
    }
    return 1;
    a++;
    b++;
    vrna_bp_stack_t *bp=NULL;
    if(bt->bp_stack==NULL){
        bp=(vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + (a-b) / 2)));
    }
    int score_bt=compute_BT(pk->fc,a,b,bt,bp);
    if(bt->bt_stack==NULL){
        write_structure(pk,bt,bp);
    }
    free(bp);
    return 1;
}

PRIVATE double INTB(pk_compound* pk,int a,int b,int c,int d){
    if(d==-1 && c==-1){
        return 1.0;
    } 
    a++;b++;
    
    return vrna_exp_E_interior_loop(pk->fc,c,d,a,b);
}

PRIVATE void add_pairing(pk_compound* pk,bt_struct* bt, int score,int a,int b,int c,int d){
    bt->structure[a]=bt->bracket;
    bt->structure[b]=tolower(bt->bracket);
    return;
}

PRIVATE int backtrace_INTB(pk_compound* pk,bt_struct* bt, int score,int a,int b,int c,int d){
    //printf("pairing %d,%d,%d,%d %d\n",a,b,c,d,score);
    add_pairing(pk,bt,score,a,b,c,d);
    if(d==-1 && c==-1){
        return 1;
    } 
    return 1;
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


PRIVATE pk_compound* create_pk_compound(vrna_fold_compound_t *fc, const char* seq,const char* ss){
    pk_compound *pk=malloc(sizeof(pk_compound));
    pk->fc=fc;
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
    double min_value=0.0; 
    for(int END=0+ANCH;END<=pk->MAX;END++){
        for(int a=0;a<END-ANCH;a++){
            double mfe1=MFEFree(pk,0,a-1);
            double mfe2=MFEFree(pk,END,pk->MAX-1);
            double mfe=add(mfe1,mfe2);
            double tmp;
            get_A(pk,NULL,NULL,&tmp,a,END);
            min_value = min_value + tmp * mfe;
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
PRIVATE void backtrace(pk_compound *pk,bt_struct* bt,double score) {
    
    for(int END=0+ANCH;END<=pk->MAX;END++){
        for(int a=0;a<END-ANCH;a++){
            double mfe1=MFEFree(pk,0,a-1);
            double mfe2=MFEFree(pk,END,pk->MAX-1);
            double mfe= mfe1 *mfe2;
            double tmp0;
            get_A(pk,NULL,NULL,&tmp0,a,END);
            if( (score-= (mfe*tmp0)) <=0.0 ){
                backtrace_A(pk,bt,tmp0,a,END);
                backtrace_MFEFree(pk,bt,mfe1,0,a-1);
                backtrace_MFEFree(pk,bt,mfe2,END,pk->MAX-1);
                return;
            }
       }
    }
}

PRIVATE int backtrace_on_interval(pk_compound *pk,bt_struct* bt,int score,int a,int END) {
    double tmp0;
    get_A(pk,NULL,NULL,&tmp0,a,END);
    if(score==0.0){
        return 0;
    }
    score=random_double(score);
    if((score-=tmp0)<=0.0){
        int res1=backtrace_A(pk,bt,tmp0,a,END);
        return res1;
    }
    
    return 0;
}

PRIVATE bt_struct *bt;

PRIVATE double scoring(vrna_fold_compound_t *fc,
                                         int i,
                                         int j,
                                         void *data) {
  pk_compound *pk=(pk_compound*)data;
  double tmp;
  Node* first;
  Node* last;
  if(get_A(pk,&first,&last,&tmp,i-1,j-1)==0){
    init_root_inter(pk,i-1,j-1);
    fold(pk);
    get_A(pk,&first,&last,&tmp,i-1,j-1);
  }
  
  return tmp;
}

PRIVATE int scoring_f(vrna_fold_compound_t *fc,
                                         int i,
                                         int j,
                                         void *data) {
  pk_compound *pk=(pk_compound*)data;
  double min=0.0;
  double free_energy;
  Node* first;
  Node* last;
  double tmp;
  //printf("-------f %d,%d\n",1,j);
  for(int tmpi=0;tmpi<j-1;tmpi++){
    
    if(get_A(pk,&first,&last,&tmp,tmpi,j-1)==0){
      init_root_inter(pk,tmpi,j-1);
      fold(pk);
      get_A(pk,&first,&last,&tmp,tmpi,j-1);
    }
    free_energy=MFEFree(pk,0,tmpi-1);
    //printf("%d %d %d\n",tmpi,free_energy,tmp);
    
      min=min + tmp * free_energy;
    }
  //printf("RETURN %d\n",min);
  return min;
}

//requiered beceause i and j are forming the base pair
PRIVATE int scoring_c(vrna_fold_compound_t *fc,
                                         int i,
                                         int j,
                                         void *data) {
  pk_compound *pk=(pk_compound*)data;
  double tmp;
  Node* first;
  Node* last;
  if(get_A(pk,&first,&last,&tmp,i,j-2)==0){
    init_root_inter(pk,i,j-2);
    fold(pk);
    get_A(pk,&first,&last,&tmp,i,j-2);
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
  //printf("custom_BT_c %d,%d, %d == %d\n",i,j,gc->matrices->c[gc->jindx[j] + i],tmp);
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
  printf("custom_BT_f %d,%d, %d \n",1,j,gc->matrices->f5[j]);
  for(;j>=2;j--){
    for(int i=0;i<j-1;i++){
        get_A(pk,&first,&last,&tmp,i,j-1);
        free_energy=MFEFree(pk,0,i-1);
        printf("%d %d %d\n",i,free_energy,tmp);
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

int setup_##NAME##(vrna_fold_compound_t *fc,const char* seq, const char* con,bt_struct* BT){
    bt=BT;
    pk_compound * pk=create_pk_compound(fc,seq,con);

    vrna_gr_add_aux_exp_f(fc,scoring_f,custom_BT_f,pk,NULL,free_aux);
    vrna_gr_add_aux_exp_c(fc,scoring_c,custom_BT_c,pk,NULL,NULL);
    vrna_gr_add_aux_exp_m(fc,scoring,custom_BT_fML,pk,NULL,NULL);
    vrna_gr_add_aux_exp_m1(fc,scoring,custom_BT_fM1,pk,NULL,NULL);
}
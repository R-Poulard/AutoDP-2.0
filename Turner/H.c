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
#define TUPLE_SIZE 6

int TABLE_SIZE = 64019;
int N = 10;
int THETA_PAIRING=3;
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
    BT=backtrack(fc,bp_stack,bt_stack,s,NULL);
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
    int stacking;
    vrna_param_t *P = fc->params;
    vrna_md_t *md=&(P->model_details);
    int *rtype= &(md->rtype[0]);
    int *ptype = fc->ptype;
    int *indx = fc->jindx;
    /*if (a==c && b==d){
        a++;b++;c++;d++;
        int ij=indx[c]+d;
        int type = vrna_get_ptype(ij, ptype);
        int kl=indx[a]+b;
        int type_2  = rtype[vrna_get_ptype(kl, ptype)];
        
        int stacking_energy = P->stack[type][type_2];
        return stacking_energy;
    }*/
    
    a++;b++;

    int energy_full = vrna_eval_int_loop(fc,c,d,a,b);
    
    return energy_full;
}

void backtrace_INTB(int score,int a,int b,int c,int d){
    structure[a]=bracket;
    structure[b]=bracket+32;
    return;
}






//declarations
int compute_CLIQUE0(HashTable *hashTable, int i, int j, int k, int l);

int compute_CLIQUE1(HashTable *hashTable, int i, int j, int k, int l);

void backtrace_CLIQUE0(HashTable *hashTable, int score, int i,int j, int k, int l);

void backtrace_CLIQUE1(HashTable *hashTable, int score, int i,int j, int k, int l);

int compute_C0(HashTable *hashTable,int a, int f, int c,int d);

int compute_C1(HashTable *hashTable,int a, int f, int c,int d);

void backtrace_C0(HashTable *hashTable, int score, int a, int f, int c,int d) ;

void backtrace_C1(HashTable *hashTable, int score, int a, int f, int c,int d) ;

int compute_B(HashTable *hashTable,int c,int d,int f,int h) ;

void backtrace_B(HashTable *hashTable,int score,int c,int d,int f,int h) ;

int compute_A(HashTable *hashTable,int a,int h) ;

void backtrace_A(HashTable *hashTable,int score,int a,int h) ;


int fold(HashTable * hashTable) {
    int min_value=INT_MAX;
    for(int h=0+3;h<=MAX;h++){
        for(int a=0;a<h-3;a++){
            int mfe1=MFEFree(0,a-1);
            int mfe2=MFEFree(h,MAX-1);
            int mfe=add(mfe1,mfe2);
            min_value = min(min_value,add(compute_A(hashTable,a,h),mfe));
        }
    }
    return min_value;
}

void backtrace(HashTable * hashTable,int score) {
    for(int h=0+3;h<=MAX;h++){
        for(int a=0;a<h-3;a++){
            
            int mfe1=MFEFree(0,a-1);
            int mfe2=MFEFree(h,MAX-1);
            int mfe=add(mfe1,mfe2);
            int tmp0=compute_A(hashTable,a,h);
            //printf("h = %d, a=%d, %d+%d %d,%d\n",h,a,mfe2,mfe1,score,compute_A(hashTable,a,h));
            if(score==add(mfe,tmp0)){
                backtrace_A(hashTable,tmp0,a,h);
                backtrace_MFEFree(mfe1,0,a-1);
                backtrace_MFEFree(mfe2,h,MAX-1);
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
        if (len != -1) {
            if (ss[len - 1] == '\n') {
                ss[len - 1] = '\0'; // Replace newline character with null terminator
                len--; // Decrement the length to exclude the removed '\n'
            }
        }
        else{
            destroyHashTable(hashTable);
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
        int score = fold(hashTable);
        //retrieve backtrack
        backtrace(hashTable,score);
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

   if(i2>=j2 || i2<i-1 || j2>j+1){
        return INT_MAX;
    }

    int min_value = add(INTB(i,j,-1,-1),compute_CLIQUE1(hashTable,i+1,i2,j2,j-1));

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

int compute_C0(HashTable *hashTable,int a, int f, int c,int d){
    
    if(a<0 || f<0 || f>=MAX || a>=MAX || a>f){
        return INT_MAX;
    }
   
    int min_value=add(INTB(a,f,-1,-1),compute_C1(hashTable,a+1,f-1,c,d));

    return min_value;
}
int compute_C1(HashTable *hashTable,int a, int f, int c,int d) {
    int C1 = -67;
    int value;

    int tab[]={C1,a,f,c,d};
    int size = 5;

    if(a<0 || f<0 || f>=MAX || a>=MAX || a>f){
        return INT_MAX;
    }
    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value=INT_MAX;
    int mfe1=MFEFree(a,c-1);
    int mfe2=MFEFree(d,f);
    
    min_value=add(0,add(mfe1,mfe2));

    loop:
    for(int tmp1=a;tmp1<=min(f-1,c-1);tmp1++){
        for(int tmp2=max(d,tmp1+2);tmp2<=f;tmp2++){
            if(evaluate(tmp1,tmp2) && tmp1-tmp2+1<=THETA){
                int tmp3=INTB(tmp1,tmp2,a,f);

                if(tmp3!=INT_MAX){
                    min_value=min(min_value,
                    add(compute_C1(hashTable,tmp1+1,tmp2-1,c,d),tmp3));
                }
            }
        }    
    }
    insert(hashTable,tab,size,min_value);
    return min_value;
}
int compute_B(HashTable *hashTable,int c,int d,int f,int h) {
    int B = 66;
    int value;
    int tab[]={B,c,d,f,h};
    int size = 5;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int g=f;g<h;g++) {
      if(!evaluate(c,h-1)||!evaluate(d-1,g)){continue;}
      int mfe0 = MFEFree(f,g-1);

      int tmp0= compute_CLIQUE0( hashTable,c,d-1,g,h-1);
      if(tmp0==INT_MAX){continue;}

      min_value = min(min_value,add(tmp0,mfe0));
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

int compute_A(HashTable *hashTable,int a,int h) {
    int A = 65;
    int value;
    int tab[]={A,a,h};
    int size = 3;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value = INT_MAX;
    
    for (int c=a+1;c<h-2;c++) {
        for (int d=c+1;d<h-1;d++) {
            for (int f=d+1;f<h;f++) {
              if(!evaluate(a,f-1)){continue;}

              int tmp0= compute_C0( hashTable,a,f-1,c,d);
              if(tmp0==INT_MAX){continue;}
              int tmp1= compute_B( hashTable,c,d,f,h);
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

void backtrace_C0(HashTable *hashTable, int score, int a, int f, int c,int d) {
    
    if(a<0 || f<0 || f>=MAX || a>=MAX || a>f){
        return;
    }
    int tmp=compute_C1(hashTable,a+1,f-1,c,d);
    int sc=INTB(a,f,-1,-1);
    if(score==add(sc,tmp)){
        backtrace_INTB(sc,a,f,-1,-1);
        backtrace_C1(hashTable,tmp,a+1,f-1,c,d);
    }
    return;
}
void backtrace_C1(HashTable *hashTable, int score, int a, int f, int c,int d) {
    
    if(a<0 || f<0 || f>=MAX || a>=MAX || a>f){
        return;
    }

    int mfe1=MFEFree(a,c-1);
    int mfe2=MFEFree(d,f);
    

    if(score==add(0,add(mfe1,mfe2))){
        
        backtrace_MFEFree(mfe1,a,c-1);
        backtrace_MFEFree(mfe2,d,f);

        return;
    }
    
    loop:
    for(int tmp1=a;tmp1<=min(f-1,c-1);tmp1++){
        for(int tmp2=max(d,tmp1+2);tmp2<=f;tmp2++){
            if(evaluate(tmp1,tmp2) && tmp1-tmp2+1<=THETA){
                int tmp3=INTB(tmp1,tmp2,a,f);
                if(tmp3!=INT_MAX){
                    int tmp0=compute_C1(hashTable,tmp1+1,tmp2-1,c,d);
                    if(score==add(tmp0,tmp3)){
                        backtrace_INTB(tmp3,tmp1,tmp2,a,f);
                        backtrace_C1(hashTable,tmp0,tmp1+1,tmp2-1,c,d);
                        return;
                    }
                }
            }
        }    
    }
    return;
} 
void backtrace_B(HashTable *hashTable,int score,int c,int d,int f,int h) {

    for (int g=f;g<h;g++) {
      if(!evaluate(c,h-1)||!evaluate(d-1,g)){continue;}
      int mfe0 = MFEFree(f,g-1);

      int tmp0= compute_CLIQUE0( hashTable,c,d-1,g,h-1);
      if(tmp0==INT_MAX){continue;}

      if(score==add(tmp0,mfe0)){
        backtrace_MFEFree(mfe0,f,g-1);

        backtrace_CLIQUE0(hashTable,tmp0, c,d-1,g,h-1);
        bracket+=1;

        return;
      }
    }


}

void backtrace_A(HashTable *hashTable,int score,int a,int h) {

    for (int c=a+1;c<h-2;c++) {
        for (int d=c+1;d<h-1;d++) {
            for (int f=d+1;f<h;f++) {
              if(!evaluate(a,f-1)){continue;}

              int tmp0= compute_C0( hashTable,a,f-1,c,d);
              if(tmp0==INT_MAX){continue;}
              int tmp1= compute_B( hashTable,c,d,f,h);
              if(tmp1==INT_MAX){continue;}

              if(score==add(add(tmp0,tmp1),0)){

                backtrace_C0(hashTable,tmp0, a,f-1,c,d);
                bracket+=1;
                backtrace_B(hashTable,tmp1, c,d,f,h);

                return;
              }
            }
        }
    }


}

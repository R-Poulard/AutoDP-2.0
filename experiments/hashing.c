#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>

#define LOAD_FACTOR_THRESHOLD 0.7
#define TUPLE_SIZE 4

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
    printf("compare tuple start\n");
    for (i = 0; i < size ; i++) {
        if (tp1->values[i] != val[i]) {
            printf("compare tuple stop1\n");
            return 0; // Tuples are different
        }
    }
    // Check if remaining elements of tp1->values are all zeros
    for (; i < TUPLE_SIZE; i++) {
        if (tp1->values[i] != 0) {
            printf("compare tuple stop2\n");
            return 0; // Tuples are different
        }
    }
    printf("compare tuple stop\n");
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
    printf("ICII\n");
    for (int i = 0; i < size; i++) {
        printf("i");
        key += keys[i] * pow(N, TUPLE_SIZE-i);
        printf("lol");
    }
    printf("la\n");
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
        printf("%d resize \n", hashTable->capacity);
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
    printf("hash retrieved");
    printf("%d\n",index);
    // Traverse the linked list at the index
    HashNode *current = hashTable->buckets[index];
    printf("%d\n",index);
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
            if(current->value!=INT_MAX){
                print_tuple(current->key);
                printf(": %d",current->value);
                
                printf("\n");
            }
            current = current->next;
        }
    }
}
//
int min(int a, int b) { if (a<b) {return a;} else {return b;}};

char * line = NULL;
int MAX;

int main() {
    HashTable *hashTable = createHashTable();

    // Insert some key-value pairs
    int search1[] = { 1, 2, 3, 4 };
    insert(hashTable, search1,4, 100);

    int search2[] = { 1, 10 };
    insert(hashTable, search2,2, 200);

    int search3[] = { 40, 50, 100 };
    insert(hashTable, search3,3, 300);

    int search4[] = { 1, 2, 3, 7 };
    insert(hashTable, search4,4, 400);

    // Retrieve values
    int value;
    printf("debut des pb");
    int search5[] = { -1, -10,13};
    if (get(hashTable, search5,3, &value)) {
        printf("Value for key (1, 2, 3, 4, 5): %d\n", value);
    }
    else {
        printf("Key (1, 2, 3, 4, 5) not found\n");
    }

    int search6[] = { 1, 11 };
    if (get(hashTable, search6,2, &value)) {
        printf("Value for key (5, 6, 734, 8, 9): %d\n", value);
    }
    else {
        printf("Key (5, 6, 7, 8, 9) not found\n");
    }

    destroyHashTable(hashTable);
    printf("We're good\n");
    return 0;
}

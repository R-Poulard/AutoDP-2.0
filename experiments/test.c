#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
/*
CLIQUE=0;//do not keep, won t run 

int compute_CLIQUE(int i, int j, int k, int l) {
    CLIQUE=1
    
    int value;
    if (get(hashTable,i,j,k,l,CLIQUE,&value)) { 
        return value;
    }
    
    int min_value = INT_MAX;

    if (k < l) { 
        min_value = min(min_value, compute_CLIQUE2(i,j,k,l-1)); 
    }
    if (i < j) {
        min_value = min(min_value, compute_CLIQUE(i+1,j,k,l));
    }
    if (i < j && k < l) {
        min_value = min(min_value, 
                    compute_CLIQUE(i+1, j, k, l-1)+bp_score(line[i],line[l]));
    }
    if (k==l) { 
        min_value = min(min_value, bp_score(line[i], line[l]));
    }
    
    insert(hashTable,i,j,k,l,CLIQUE,min_value);
    return min_value;
}

int compute_CLIQUE2(int i, int j, int k, int l) {
    CLIQUE2=-1


    if (get(hashTable,i,j,k,l,CLIQUE2,&value)) { 
        return value;
    }

    int min_value = INT_MAX;

    if (k < l) { 
        min_value = min(min_value, compute_CLIQUE2(i,j,k,l-1)); 
    }
    if (i < j && k < l) {
        min_value = min(min_value, 
                    compute_CLIQUE(i+1, j, k, l-1)+bp_score(line[i],line[l]));
    }
    if (k==l) { 
        min_value = min(min_value, bp_score(line[i], line[l]));
    }

    insert(hashTable,i,j,k,l,CLIQUE2,min_value);
    return min_value;
}
*/
void printBinary(char ch) {
    // Iterate through each bit from left to right
    for (int i = 7; i >= 0; i--) {
        // Use bitwise AND to check if the i-th bit is set
        if (ch & (1 << i)) {
            printf("1");
        } else {
            printf("0");
        }
    }
    printf("\n");
}
int TUPLE_SIZE=5;
int N=10;
int hash(int keys[],int size) {
    int key = 0;
    for (int i = 0; i < size; i++) {
        key += keys[i] * pow(N, TUPLE_SIZE-i);
    }
    return key;
}

typedef struct Tuple {
    int values[5];
} Tuple;

Tuple *createTuple(int val[],int size) {
    Tuple *tuple = (Tuple *)malloc(sizeof(Tuple));

    // Initialize all values in tuple to 0
    for (int i = 0; i < TUPLE_SIZE; i++) {
        tuple->values[i] = 0;
    }

    // Copy values from val array, up to its length or TUPLE_SIZE
    for (int i = 0; i < size; i++) {
        if (val[i] != '\0') {
            tuple->values[i] = val[i];
        }
        else {
            break; // Stop copying if val[i] is '\0' indicating end of values
        }
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

int main(int argc, char ** argv) {
    int tp[]={1,2,3}; 
    Tuple *tp1=createTuple(tp,3);
    for(int i=0;i<TUPLE_SIZE;i++){
        printf("%d  ",tp1->values[i]);
    }
    int tb[]={1,2,3}; // Example character
    printf(" cmp %d\n",compare_tuple(tp1,tb,3));
    exit(0);
    int tab[]={1,2,3,4,5}; // Example character
    printf("%d\n",hash(tab,5));

    int tab2[]={10,2,14,25,5}; // Example character
    printf("%d\n",hash(tab2,5));

    int tab3[]={1,22,33,42,50}; // Example character
    printf("%d\n",hash(tab3,5));
    return 0;
}
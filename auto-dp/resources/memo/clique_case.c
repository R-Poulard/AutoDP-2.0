
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
    int evaluate;
    if (k < l) { 
        evaluate=up_score(l,l);
        if(evaluate!=INT_MAX){
            min_value = min(min_value, add(compute_CLIQUE2(hashTable,i,j,k,l-1),evaluate)); 
        }
    }
    if (i < j ) {
        evaluate=up_score(i,i);
        if(evaluate!=INT_MAX){
            min_value = min(min_value,add(compute_CLIQUE(hashTable,i+1,j,k,l),evaluate));
        }
    }
    if (i < j && k < l) {
        evaluate=bp_score(i,l);
        if(evaluate!=INT_MAX){
            min_value = min(min_value, add(compute_CLIQUE(hashTable,i+1, j, k, l-1),evaluate));
        }
    }
    if (k==l) {
        evaluate=bp_score(i,l);
        if(evaluate!=INT_MAX){ 
            min_value = min(min_value,evaluate);
        }
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

    int evaluate;
    if (k < l) { 
        evaluate=up_score(l,l);
        if(evaluate!=INT_MAX){
            min_value = min(min_value, add(compute_CLIQUE2(hashTable,i,j,k,l-1),evaluate)); 
        }
    }
    if (i < j && k < l) {
        evaluate=bp_score(i,l);
        if(evaluate!=INT_MAX){
            min_value = min(min_value, add(compute_CLIQUE(hashTable,i+1, j, k, l-1),evaluate));
        }
    }
    if (k==l) {
        evaluate=bp_score(i,l);
        if(evaluate!=INT_MAX){ 
            min_value = min(min_value,evaluate);
        }
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}
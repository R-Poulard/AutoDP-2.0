int compute_CLIQUE0(HashTable *hashTable,int i, int i2, int j2, int j){
    int CLIQUE=-1;
    int value;
    int tab[] = {CLIQUE,i,i2,j2,j};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        return INT_MAX;
    }

    if (get(hashTable,tab,size,&value)) { 
        return value;
    }
    
    int min_value = add(INTB(i,j,i,j),compute_CLIQUE1(hashTable,i+1,i2,j2,j-1));

    insert(hashTable,tab,size,min_value);
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

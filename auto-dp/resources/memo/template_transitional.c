int compute_MAINNAME(HashTable *hashTable,INT_INDICES) {
    // index_start
    int MAINNAME = INT_VALUE_NAME;
    int value;
    int tab[]={MAINNAME,INDICES};
    int size = SIZE_ARRAY;

    if (get(hashTable,tab,size,&value)){
        return value;
    }
    // index_end    
    int min_value = INT_MAX;
    
FOR_LOOP_NEW_VARIABLES_OPEN
CHILDREN_SUM
  INDENTmin_value = min(min_value, CHILDREN_MAX);
FOR_LOOP_NEW_VARIABLES_CLOSE
    // index_start
    insert(hashTable,tab,size,min_value);
    // index_end
    return min_value;
}

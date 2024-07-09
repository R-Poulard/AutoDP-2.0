int compute_MAINNAME(HashTable *hashTable,int V1, int V2, CONST_INT) {
    int MAINNAME = INT_VALUE_NAME;
    int value;

    int tab[]={MAINNAME,V1,V2,CONST};
    int size = SIZE_ARRAY;

    if(V1<0 || V2<0 || V2>=MAX || V1>=MAX || V1>V2){
        return INT_MAX;
    }
    if (get(hashTable,tab,size,&value)){
        return value;
    }

    int min_value = INT_MAX;
    bool eq_some_const1 = V1_INC1_CONST_COMP;
    bool eq_some_const2 = V2_INC2_CONST_COMP;
    int evaluate;

    if (eq_some_const2) {
        evaluate=up_score(V1,V1);
        if (evaluate!=INT_MAX){
            min_value = min(min_value, add(compute_MAINNAME(hashTable,V1INC1,V2,CONST),evaluate));
        }
    }
    if (eq_some_const1) {
        evaluate=up_score(V2,V2);
        if (evaluate!=INT_MAX){
            min_value = min(min_value, add(compute_MAINNAME2(hashTable,V1, V2INC2, CONST),evaluate));
        }
    }
    if (eq_some_const1 && eq_some_const2) {
        evaluate=bp_score(V1,V2);
        if (evaluate!=INT_MAX){
            min_value = min(min_value, add(compute_MAINNAME(hashTable,V1INC1, V2INC2, CONST),evaluate));
        }
    }
    CHILDREN_SUM
    min_value = min(min_value, CHILDREN_MAX);

    insert(hashTable,tab,size,min_value);
    return min_value;
} 

int compute_MAINNAME2(HashTable *hashTable,int V1, int V2, CONST_INT) {
    int MAINNAME2 = INT_VALUE_SECOND_NAME;
    int value;

    int tab[]={MAINNAME2,V1,V2,CONST};
    int size = SIZE_ARRAY;

    if(V1<0 || V2<0 || V2>=MAX || V1>=MAX || V1>V2){
        return INT_MAX;
    }

    if (get(hashTable,tab,size,&value)){
        return value;
    }

    int min_value = INT_MAX;
    bool eq_some_const1 = V1_INC1_CONST_COMP;
    bool eq_some_const2 = V2_INC2_CONST_COMP;
    int evaluate;
    
    if (eq_some_const1) {
        evaluate=up_score(V2,V2);
        if (evaluate!=INT_MAX){
            min_value = min(min_value, add(compute_MAINNAME2(hashTable,V1, V2INC2, CONST),evaluate));
        }
    }
    if (eq_some_const1 && eq_some_const2) {
        evaluate=bp_score(V1,V2);
        if (evaluate!=INT_MAX){
            min_value = min(min_value, add(compute_MAINNAME(hashTable,V1INC1, V2INC2, CONST),evaluate));
        }
    }

    insert(hashTable,tab,size,min_value);
    return min_value;
}

int compute_MAINNAME(int V1, int V2, CONST_INT) {
    int MAINNAME = INT_VALUE_NAME;
    int value;
    if (get(hashTable,CONST,MAINNAME,&value)){
        return value;
    }
    /*
    if (MAINNAME[index_MAINNAME(V1,V2, CONST)] > INT_MIN) {
        return MAINNAME[index_MAINNAME(V1,V2, CONST)];
    }*/

    int min_value = INT_MAX;
    bool eq_some_const1 = V1_INC1_CONST_COMP;
    bool eq_some_const2 = V2_INC2_CONST_COMP;

    if (V1INC1 < V2) {
        if (!eq_some_const1) {
            min_value = min(min_value, compute_MAINNAME(V1INC1,V2,CONST));
        }
    }
    if (V2INC2 > V1) {
        if (!eq_some_const2) {
            min_value = min(min_value, compute_MAINNAME2(V1, V2INC2, CONST));
        }
    }
    if (!eq_some_const1 && !eq_some_const2) {
        min_value = min(min_value, compute_MAINNAME(V1INC1, V2INC2, CONST)+bp_score(line[V1], 
                                                                                      line[V2]));
    }

    min_value = min(min_value, CHILDREN_SUM);

    insert(hashTable,CONST,MAINNAME,min_value);
    //MAINNAME[index_MAINNAME(V1,V2,CONST)] = min_value;
    return min_value;
} 

int compute_MAINNAME2(int V1, int V2, CONST_INT) {
    int MAINNAME2 = INT_VALUE_SECOND_NAME;
    int value;
    if (get(hashTable,CONST,MAINNAME2,&value)){
        return value;
    }
    /*
    if (MAINNAME2[index_MAINNAME(V1,V2, CONST)] > INT_MIN) {
        return MAINNAME2[index_MAINNAME(V1,V2, CONST)];
    }
    */
    int min_value = INT_MAX;
    bool eq_some_const1 = V1_INC1_CONST_COMP;
    bool eq_some_const2 = V2_INC2_CONST_COMP;

    if (V2INC2 > V1 && !eq_some_const2) {
        min_value = min(min_value, compute_MAINNAME2(V1, V2INC2, CONST));
    }
    if (!eq_some_const1 && !eq_some_const2) {
        min_value = min(min_value, compute_MAINNAME(V1INC1, V2INC2, CONST)+bp_score(line[V1], 
                                                                                      line[V2]));
    }

    insert(hashTable,CONST,MAINNAME2,min_value);
    //MAINNAME2[index_MAINNAME(V1,V2,CONST)] = min_value;
    return min_value;
}

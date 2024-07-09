void backtrace_MAINNAME(HashTable *hashTable,int score,int V1, int V2, CONST_INT) {

    if(V1<0 || V2<0 || V2>=MAX || V1>=MAX || V1>V2){
        return;
    }

    bool eq_some_const1 = V1_INC1_CONST_COMP;
    bool eq_some_const2 = V2_INC2_CONST_COMP;
    int evaluate;
    int tmp;
    if (eq_some_const2) {
        evaluate=up_score(V1,V1);
        if(evaluate!=INT_MAX){
            tmp=compute_MAINNAME(hashTable,V1INC1, V2, CONST);
            if( score==add(tmp,evaluate)){
                backtrace_MAINNAME(hashTable,tmp,V1INC1, V2, CONST);
                return;
            }
        }
    }
    if (eq_some_const1) {
        evaluate=up_score(V2,V2);
        if(evaluate!=INT_MAX){
            tmp=compute_MAINNAME2(hashTable,V1, V2INC2, CONST);
            if(score==add(tmp,evaluate)){
                backtrace_MAINNAME2(hashTable,tmp,V1, V2INC2, CONST);
                return;
            }
        }
    }
    if (eq_some_const1 && eq_some_const2) {
        evaluate=bp_score(V1,V2);
        if(evaluate!=INT_MAX){
            tmp=compute_MAINNAME(hashTable,V1INC1, V2INC2, CONST);
            if(score==add(tmp,evaluate)){
                backtrace_MAINNAME(hashTable,tmp,V1INC1, V2INC2, CONST);
                structure[V1]=bracket;
                structure[V2]=bracket+32;
                return;
            }
        }
        
    }
    CHILDREN_SUM
    if(score==CHILDREN_MAX){
        CHILDREN_BACKTRACE
        return;
    }

    return;
} 

void backtrace_MAINNAME2(HashTable *hashTable,int score, int V1, int V2, CONST_INT) {
    
    if(V1<0 || V2<0 || V2>=MAX || V1>=MAX || V1>V2){
        return;
    }

    bool eq_some_const1 = V1_INC1_CONST_COMP;
    bool eq_some_const2 = V2_INC2_CONST_COMP;
    int evaluate;
    int tmp;
    if (V2INC2 > V1 && eq_some_const1) {
        evaluate=up_score(V2,V2);
        if(evaluate!=INT_MAX){
            tmp=compute_MAINNAME2(hashTable,V1, V2INC2, CONST);
            if(score==add(tmp,evaluate)){
                backtrace_MAINNAME2(hashTable,tmp,V1, V2INC2, CONST);
                return;
            }
        }
    }

    if (eq_some_const1 && eq_some_const2) {
        evaluate=bp_score(V1,V2);
        if(evaluate!=INT_MAX){
            tmp=compute_MAINNAME(hashTable,V1INC1, V2INC2, CONST);
            if(score==add(tmp,evaluate)){
                backtrace_MAINNAME(hashTable,tmp,V1INC1, V2INC2, CONST);
                structure[V1]=bracket;
                structure[V2]=bracket+32;
                return;
            }
        }
        
    }

    return;
}

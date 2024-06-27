void backtrace_MAINNAME0(HashTable *hashTable, int score, int V1, int V2, CONST_INT) {

    if(V1<0 || V2<0 || V2>=MAX || V1>=MAX || V1>V2){
        return;
    }
    int tmp;
    get_MAINNAME1(hashTable,NULL,NULL,&tmp,V1+1,V2-1,CONST);
    int sc=INTB(V1,V2,-1,-1);
    if(score==add(sc,tmp)){
        backtrace_INTB(sc,V1,V2,-1,-1);
        backtrace_MAINNAME1(hashTable,tmp,V1+1,V2-1,CONST);
    }
    return;
}
void backtrace_MAINNAME1(HashTable *hashTable, int score, int V1, int V2, CONST_INT) {

    if(V1<0 || V2<0 || V2>=MAX || V1>=MAX || V1>V2){
        return;
    }

    CONST_SUM
    CHILDREN_SUM

    if(score==add(CHILDREN_MAX,CONST_MAX)){
        CHILDREN_BACKTRACE
        CONST_BACKTRACE
        return;
    }
    
    loop:
    for(int tmp1=V1;tmp1<=min(V2-1,BOUNDIADJI);tmp1++){
        for(int tmp2=max(BOUNDJADJJ,tmp1+NB_HELIX);tmp2<=V2;tmp2++){
            if(evaluate(tmp1,tmp2) && tmp1-tmp2+1<=THETA){
                int tmp3=INTB(tmp1,tmp2,V1,V2);
                if(tmp3!=INT_MAX){
                    int tmp0;
                    get_MAINNAME1(hashTable,NULL,NULL,&tmp0,tmp1+1,tmp2-1,CONST);
                    if(score==add(tmp0,tmp3)){
                        backtrace_INTB(tmp3,tmp1,tmp2,V1,V2);
                        backtrace_MAINNAME1(hashTable,tmp0,tmp1+1,tmp2-1,CONST);
                        return;
                    }
                }
            }
        }    
    }
    return;
} 
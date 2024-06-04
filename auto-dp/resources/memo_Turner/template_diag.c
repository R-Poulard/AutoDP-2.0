int compute_MAINNAME0(HashTable *hashTable,int V1, int V2, CONST_INT){
    
    if(V1<0 || V2<0 || V2>=MAX || V1>=MAX || V1>V2){
        return INT_MAX;
    }
   
    int min_value=add(INTB(V1,V2,-1,-1),compute_MAINNAME1(hashTable,V1+1,V2-1,CONST));

    return min_value;
}
int compute_MAINNAME1(HashTable *hashTable,int V1, int V2, CONST_INT) {
    int MAINNAME1 = -INT_VALUE_NAME;
    int value;

    int tab[]={MAINNAME1,V1,V2,CONST};
    int size = SIZE_ARRAY;

    if(V1<0 || V2<0 || V2>=MAX || V1>=MAX || V1>V2){
        return INT_MAX;
    }
    if (get(hashTable,tab,size,&value)){
        return value;
    }
    int min_value=INT_MAX;
    CONST_SUM
    CHILDREN_SUM
    min_value=add(CHILDREN_MAX,CONST_MAX);

    loop:
    for(int tmp1=V1;tmp1<=min(V2INC2,BOUNDIADJI);tmp1++){
        for(int tmp2=max(BOUNDJADJJ,tmp1+NB_HELIX);tmp2<=V2;tmp2++){
            if(evaluate(tmp1,tmp2) && tmp1-tmp2+1<=THETA){
                int tmp3=INTB(tmp1,tmp2,V1,V2);

                if(tmp3!=INT_MAX){
                    min_value=min(min_value,
                    add(compute_MAINNAME1(hashTable,tmp1+1,tmp2-1,CONST),tmp3));
                }
            }
        }    
    }
    insert(hashTable,tab,size,min_value);
    return min_value;
}
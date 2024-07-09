int get_MAINNAME0(pk_compound* pk,Node** first,Node** last,int * value,int V1, int V2, CONST_INT){
    int MAINNAME1 = INT_VALUE_NAME;
    int tab[]={MAINNAME1,V1+1,V2-1,CONST};
    int size = SIZE_ARRAY;

    if(V1<0 || V2<0 || V2>=pk->MAX || V1>=pk->MAX || V1>V2){
        *value=INT_MAX;
        return 1;
    }   
    if (get(pk->hashtable,tab,size,value)){
        *value=add(INTB(pk,V1,V2,-1,-1),*value);
        return 1;
    }
    else{
        Node* newNode = (Node*)malloc(sizeof(Node));
        newNode->data = createTuple(tab,size);
        newNode->next = NULL;
        if(*last!=NULL){
            (*last)->next=newNode;
        }
        if(*first==NULL){
            *first=newNode;
        }
        *last=newNode;
        return 0;
    }
}

int get_MAINNAME1(pk_compound *pk,Node** first,Node** last,int * value,int V1, int V2, CONST_INT){
    int MAINNAME1 = INT_VALUE_NAME;
    int tab[]={MAINNAME1,V1,V2,CONST};
    int size = SIZE_ARRAY;

    if(V1<0 || V2<0 || V2>=pk->MAX || V1>=pk->MAX || V1>V2){
        *value=INT_MAX;
        return 1;
    }  
    if (get(pk->hashtable,tab,size,value)){
        //printf("value=>%d \n",*value);
        return 1;
    }
    else{
        Node* newNode = (Node*)malloc(sizeof(Node));
        newNode->data = createTuple(tab,size);
        newNode->next = NULL;
        if(*last!=NULL){
            (*last)->next=newNode;
        }
        if(*first==NULL){
            *first=newNode;
        }
        *last=newNode;
        return 0;
    }
}

int compute_MAINNAME1(pk_compound *pk,Node** first,Node** last,int V1, int V2, CONST_INT) {
    int MAINNAME1 = INT_VALUE_NAME;
    int tab[]={MAINNAME1,V1,V2,CONST};
    int size = SIZE_ARRAY;
    int value;
    if(V1<0 || V2<0 || V2>=pk->MAX || V1>=pk->MAX || V1>V2){
        return 1;
    } 
    if (get(pk->hashtable,tab,size,&value)){
        return 1;
    }
    int possible=1;
    int min_value=INT_MAX;
    CONST_SUM
    CHILDREN_SUM
    if(possible==1){
    min_value=add(CHILDREN_MAX,CONST_MAX);
    }
    
    loop:
    for(int tmp1=V1;tmp1<=min(V2INC2,BOUNDIADJI);tmp1++){
        for(int tmp2=max(BOUNDJADJJ,tmp1+NB_HELIX);tmp2<=V2;tmp2++){
            if(evaluate(pk,tmp1,tmp2) && tmp1-tmp2+1<=pk->THETA){
                int tmp3=INTB(pk,tmp1,tmp2,V1,V2);

                if(tmp3!=INT_MAX){
                    int tmp;
                    if(get_MAINNAME1(pk,first,last,&tmp,tmp1+1,tmp2-1,CONST)==0){
                        possible=0;
                    }
                    if(possible==1){
                        //printf("C1 res= %d",tmp);
                        min_value=min(min_value, add(tmp,tmp3));
                    }
                }
            }
        }    
    }
    if(possible==1){
        insert(pk->hashtable,tab,size,min_value);
    }
    return possible;
}

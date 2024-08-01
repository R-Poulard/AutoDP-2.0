int get_MAINNAME0(HashTable *hashTable,Node** first,Node** last,double * value,int V1, int V2, CONST_INT){
    int MAINNAME1 = INT_VALUE_NAME;
    int tab[]={MAINNAME1,V1-1,V2+1,CONST};
    int size = SIZE_ARRAY;

    if(V1<0 || V2<0 || V2>=MAX || V1>=MAX || V1>V2){
        *value=0.;
        return 1;
    }   
    if (get(hashTable,tab,size,value)){
        *value=INTB(V1,V2,-1,-1)+(*value);
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

int get_MAINNAME1(HashTable *hashTable,Node** first,Node** last,double * value,int V1, int V2, CONST_INT){
    int MAINNAME1 = INT_VALUE_NAME;
    int tab[]={MAINNAME1,V1,V2,CONST};
    int size = SIZE_ARRAY;

    if(V1<-1 || V2<-1 || V2>=MAX || V1>=MAX || V1>V2){
        *value=0.;
        return 1;
    }  
    if (get(hashTable,tab,size,value)){
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

int compute_MAINNAME1(HashTable *hashTable,Node** first,Node** last,int V1, int V2, CONST_INT) {
    int MAINNAME1 = INT_VALUE_NAME;
    int tab[]={MAINNAME1,V1,V2,CONST};
    int size = SIZE_ARRAY;
    double value;
    
    if(V1<-1 || V2<-1 || V2>=MAX || V1>=MAX || V1>V2){
        *value=0.;
        return 1;
    } 
    if (get(hashTable,tab,size,&value)){
        return 1;
    }
    int possible=1;
    double min_value=0.;
    CONST_SUM
    CHILDREN_SUM
    min_value= CHILDREN_MAX * CONST_MAX;

    
    loop:
    for(int tmp1=max(0,BOUNDIADJI);tmp1<=V1;tmp1++){
        for(int tmp2=V2;tmp2<=min(MAX-1,BOUNDJADJJ);tmp2++){
            if(evaluate(tmp1,tmp2) && tmp1-tmp2+1<=THETA){
                int tmp3=INTB(V1,V2,tmp1,tmp2);

                if(tmp3!=INT_MAX){
                    int tmp;
                    if(get_MAINNAME1(hashTable,first,last,&tmp,tmp1-1,tmp2+1,CONST)==0){
                        possible=0;
                    }
                    if(possible==1){
                        //printf("C1 res= %d",tmp);
                        min_value= min_value + tmp*tmp3);
                    }
                }
            }
        }    
    }
    if(possible==1){
        insert(hashTable,tab,size,min_value);
    }
    return possible;
}

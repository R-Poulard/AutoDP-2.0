int get_CLIQUE0(HashTable *hashTable,Node** first,Node** last,int * value,int i, int i2, int j2, int j){
    int CLIQUE=1;
    int tab[] = {CLIQUE,i+1,i2,j2,j-1};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        *value=INT_MAX;
        return 1;
    }
    if (get(hashTable,tab,size,value)){
        //printf("*value = %d(+%d)\n",*value,INTB(i,j,-1,-1));
        *value=add(INTB(i,j,-1,-1),*value);
        //printf("*new _ value = %d\n",*value);
        return 1;
    }
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

int get_CLIQUE1(HashTable *hashTable,Node** first,Node** last,int * value,int i, int i2, int j2, int j){
    int CLIQUE=1;
    int tab[] = {CLIQUE,i,i2,j2,j};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        *value=INT_MAX;
        return 1;
    }
    if (get(hashTable,tab,size,value)){
        return 1;
    }
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


int compute_CLIQUE1(HashTable *hashTable,Node** first,Node** last,int i, int i2, int j2, int j) {
    int CLIQUE=1;
    int value;
    int tab[] = {CLIQUE,i,i2,j2,j};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        return 1;
    }

    if (get(hashTable,tab,size,&value)) { 
        return 1;
    }
    int possible=1;
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
                    int tmp2;
                    if(get_CLIQUE1(hashTable,first,last,&tmp2,k+1,i2,j2,l-1)==0){
                        possible=0;
                    }
                    if(possible==1){
                        
                    min_value = min(min_value,add(tmp2,tmp));
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

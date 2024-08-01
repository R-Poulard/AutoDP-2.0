
PRIVATE int get_CLIQUE0(pk_compound* pk,Node** first,Node** last,double * value,int i, int i2, int j2, int j){
    int CLIQUE=1;
    int tab[] = {CLIQUE,i+1,i2,j2,j-1};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        *value=0.;
        return 1;
    }
    if (get(pk->hashtable,tab,size,value)){
        //printf("*value = %d(+%d)\n",*value,INTB(i,j,-1,-1));
        *value=INTB(pk,i,j,-1,-1)*(*value);
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

PRIVATE int get_CLIQUE1(pk_compound* pk,Node** first,Node** last,double * value,int i, int i2, int j2, int j){
    int CLIQUE=1;
    int tab[] = {CLIQUE,i,i2,j2,j};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        *value=0.;
        return 1;
    }
    if (get(pk->hashtable,tab,size,value)){
        return 1;
    }
    Node* newNode = (Node*)malloc(sizeof(Node));
    newNode->data = createTuple(tab,size);
    newNode->next = NULL;
    if(last!=NULL && *last!=NULL){
        (*last)->next=newNode;
    }
    if(first!=NULL && *first==NULL){
        *first=newNode;
    }
    *last=newNode;
    return 0;
}


PRIVATE int compute_CLIQUE1(pk_compound* pk,Node** first,Node** last,int i, int i2, int j2, int j) {
    int CLIQUE=1;
    double value;
    int tab[] = {CLIQUE,i,i2,j2,j};
    int size= 5;

   if(i2>=j2 || i2<i-1 || j2>j+1){
        return 1;
    }

    if (get(pk->hashtable,tab,size,&value)) { 
        return 1;
    }
    int possible=1;
    double min_value = 0.;
    double tmp;
    if (j == j2 && i <= i2) { 
        min_value = min_value + MFEFree(pk,i,i2); 
    }
    for(int k=i;k<=i2;k++){
        for(int l=max(j2,k+1);l<=j;l++){
            if(evaluate(pk,k,l) && (k<i2 || l==j2)){
                tmp=INTB(pk,k,l,i,j);
                if(tmp!=0.0){
                    double tmp2;
                    if(get_CLIQUE1(pk,first,last,&tmp2,k+1,i2,j2,l-1)==0){
                        possible=0;
                    }
                    if(possible==1){
                        
                    min_value = min_value + tmp2*tmp;
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

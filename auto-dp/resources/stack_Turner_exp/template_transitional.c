PRIVATE int get_MAINNAME(pk_compound *pk,Node** first,Node** last,double * value,INT_INDICES){
    int MAINNAME = INT_VALUE_NAME;
    int tab[]={MAINNAME,INDICES};
    int size = SIZE_ARRAY;


    if (get(pk->hashtable,tab,size,value)){
        //printf("found\n");
        return 1;
    }
    else{
        //printf("not found\n");
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
}
PRIVATE int compute_MAINNAME(pk_compound *pk,Node** first,Node** last,INT_INDICES) {
    // index_start
    int MAINNAME = INT_VALUE_NAME;
    double value;
    int tab[]={MAINNAME,INDICES};
    int size = SIZE_ARRAY;

    if (get(pk->hashtable,tab,size,&value)){
        return 1;
    }
    // index_end    
    double min_value = 0.;
    int possible=1;
FOR_LOOP_NEW_VARIABLES_OPEN
CONDITIONS
MFE_SUM
CHILDREN_SUM
  INDENTif(possible==1){
  INDENT  min_value = min_value + CHILDREN_MAX + MFE_MAX;
  INDENT}
FOR_LOOP_NEW_VARIABLES_CLOSE
    if(possible==1){
        insert(pk->hashtable,tab,size,min_value);
    }
    return possible;
}

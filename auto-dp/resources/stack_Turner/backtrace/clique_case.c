PRIVATE int backtrace_CLIQUE0(pk_compound *pk,bt_struct* bt,int score,int i, int i2, int j2, int j) {
    
   if(i2>=j2 || i2<i || j2>j){
        return 1;
    }
    int tmp;
    get_CLIQUE1(pk,NULL,NULL,&tmp,i+1,i2,j2,j-1);
    int sc=INTB(pk,i,j,-1,-1);
    if(score==add(sc,tmp)){
        int res1=backtrace_INTB(pk,bt,sc,i,j,-1,-1);
        return res1 && backtrace_CLIQUE1(pk,bt,tmp,i+1,i2,j2,j-1);
    }

    return 0;
}

PRIVATE int backtrace_CLIQUE1(pk_compound *pk,bt_struct* bt,int score,int i, int i2, int j2, int j) {
    
   if(i2>=j2 || i2<i || j2>j){
        return 1;
    }
    int tmp;
    
    if (j == j2 && i <= i2 && score==MFEFree(pk,i,i2)) { 
        return backtrace_MFEFree(pk,bt,score,i,i2);
    }
    for(int k=i;k<=i2;k++){
        for(int l=max(j2,k+1);l<=j;l++){
            if(evaluate(pk,k,l) && (k<i2 || l==j2)){
                tmp=INTB(pk,k,l,i,j);
                if(tmp!=INF){
                    int tmp0;
                    get_CLIQUE1(pk,NULL,NULL,&tmp0,k+1,i2,j2,l-1);
                    if(score==add(tmp0,tmp)){
                        int res1=backtrace_INTB(pk,bt,tmp,k,l,i,j); 
                        int res2=backtrace_CLIQUE1(pk,bt,tmp0,k+1,i2,j2,l-1);                       
                        return res1 && res2;
                    }
                }
            }
        }
    }
    return 0;
}
PRIVATE int backtrace_CLIQUE0(pk_compound *pk,bt_struct* bt,double score,int i, int i2, int j2, int j) {
    
   if(i2>=j2 || i2<i || j2>j){
        return 1;
    }
    double tmp;
    get_CLIQUE1(pk,NULL,NULL,&tmp,i+1,i2,j2,j-1);
    double sc=INTB(pk,i,j,-1,-1);
    score=random_double(score);
    if((score-=(sc*tmp)) <=0.0){
        int res1=backtrace_INTB(pk,bt,sc,i,j,-1,-1);
        return res1 && backtrace_CLIQUE1(pk,bt,tmp,i+1,i2,j2,j-1);
    }

    return 0;
}

PRIVATE int backtrace_CLIQUE1(pk_compound *pk,bt_struct* bt,double score,int i, int i2, int j2, int j) {
    
   if(i2>=j2 || i2<i || j2>j){
        return 1;
    }
    double tmp;
    score=random_double(score);
    if (j == j2 && i <= i2 ) { 
        score -= MFEFree(pk,i,i2);
        if(score<=0.0){
        return backtrace_MFEFree(pk,bt,MFEFree(pk,i,i2),i,i2);
        }
    }
    for(int k=i;k<=i2;k++){
        for(int l=max(j2,k+1);l<=j;l++){
            if(evaluate(pk,k,l) && (k<i2 || l==j2)){
                tmp=INTB(pk,k,l,i,j);
                if(tmp!=0.0){
                    double tmp0;
                    get_CLIQUE1(pk,NULL,NULL,&tmp0,k+1,i2,j2,l-1);
                    if( (score-=(tmp0*tmp))<=0.0 ){
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
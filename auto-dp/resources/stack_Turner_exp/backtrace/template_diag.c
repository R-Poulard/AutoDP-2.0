PRIVATE int backtrace_MAINNAME0(pk_compound *pk,bt_struct* bt, double score, int V1, int V2, CONST_INT) {

    if(V1<0 || V2<0 || V2>=pk->MAX || V1>=pk->MAX || V1>V2){
        return 1;
    }
    double tmp;
    get_MAINNAME1(pk,NULL,NULL,&tmp,V1+1,V2-1,CONST);
    double sc=INTB(pk,V1,V2,-1,-1);
    score=random_double(score);
    if((score-=(sc*tmp)) <=0.0){
        int res1=backtrace_INTB(pk,bt,sc,V1,V2,-1,-1);
        return res1 && backtrace_MAINNAME1(pk,bt,tmp,V1+1,V2-1,CONST);
    }
    return 0;
}
PRIVATE int backtrace_MAINNAME1(pk_compound *pk,bt_struct* bt, double score, int V1, int V2, CONST_INT) {

    if(V1<0 || V2<0 || V2>=pk->MAX || V1>=pk->MAX || V1>V2){
        return 1;
    }

    CONST_SUM
    CHILDREN_SUM
    score=random_double(score);
    if( (score-=(CHILDREN_MAX*CONST_MAX))<=0.0){
        CHILDREN_BACKTRACE
        CONST_BACKTRACE
        return RETURN_CHECK;
    }
    
    loop:
    for(int tmp1=V1;tmp1<=min(V2-1,BOUNDIADJI);tmp1++){
        for(int tmp2=max(BOUNDJADJJ,tmp1+NB_HELIX);tmp2<=V2;tmp2++){
            if(evaluate(pk,tmp1,tmp2)){
                double tmp3=INTB(pk,tmp1,tmp2,V1,V2);
                if(tmp3!=INF){
                    double tmp0;
                    get_MAINNAME1(pk,NULL,NULL,&tmp0,tmp1+1,tmp2-1,CONST);
                    if((score-=(tmp0 * tmp3))<=0.0 ){
                        int res1 = backtrace_INTB(pk,bt,tmp3,tmp1,tmp2,V1,V2);
                        return res1 && backtrace_MAINNAME1(pk,bt,tmp0,tmp1+1,tmp2-1,CONST);
                        
                    }
                }
            }
        }    
    }
    return 0;
} 
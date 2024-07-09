void backtrace_CLIQUE(HashTable *hashTable, int score, int i,int j, int k, int l){

    if(j>=k || j<i || k>l){
        return;
    }
    int evaluate;
    int tmp;
    if (k < l) { 
        evaluate=up_score(l,l);
        if(evaluate!=INT_MAX){
            tmp=compute_CLIQUE2(hashTable,i,j,k,l-1);
            if(score==add(tmp,evaluate)){
                backtrace_CLIQUE2(hashTable,tmp,i,j,k,l-1);
                return;
            }
        }
    }
    if (i < j) {
        evaluate=up_score(i,i);
        if(evaluate!=INT_MAX){
            tmp=compute_CLIQUE(hashTable,i+1,j,k,l);
            if(score==add(tmp,evaluate)){
                backtrace_CLIQUE(hashTable,tmp,i+1,j,k,l);
                return;
            }
        }
    }
    if (i < j && k < l) {
        evaluate=bp_score(i,l);
        if(evaluate!=INT_MAX){
            tmp=compute_CLIQUE(hashTable,i+1, j, k, l-1);
            if(score == add(tmp,evaluate)){
                backtrace_CLIQUE(hashTable,tmp,i+1,j,k,l-1);
                structure[i]=bracket;
                structure[l]=bracket+32;
                return;
            }
        }
    }
    if (k==l) { 
        if(score== bp_score(i, l)){
            structure[i]=bracket;
            structure[l]=bracket+32;
            return;
        }
    }
}

void backtrace_CLIQUE2(HashTable *hashTable,int score,int i,int j,int k,int l){

    if(j>=k || j<i || k>l){
        return;
    }
    int evaluate;
    int tmp;
    if (k < l) { 
        evaluate=up_score(l,l);
        if(evaluate!=INT_MAX){
            tmp=compute_CLIQUE2(hashTable,i,j,k,l-1);
            if(score==add(tmp,evaluate)){
                backtrace_CLIQUE2(hashTable,tmp,i,j,k,l-1);
                return;
            }
        }
    }

    if (i < j && k < l) {
        evaluate=bp_score(i,l);
        if(evaluate!=INT_MAX){
            tmp=compute_CLIQUE(hashTable,i+1, j, k, l-1);
            if(score == add(tmp,evaluate)){
                backtrace_CLIQUE(hashTable,tmp,i+1,j,k,l-1);
                structure[i]=bracket;
                structure[l]=bracket+32;
                return;
            }
        }
    }

    if (k==l) { 
        if(score== bp_score(i, l)){
            structure[i]=bracket;
            structure[l]=bracket+32;
            return;
        }
    }

}

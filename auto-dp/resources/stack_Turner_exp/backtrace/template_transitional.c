PRIVATE int backtrace_MAINNAME(pk_compound *pk,bt_struct* bt,double score,INT_INDICES) {
    // index_start
score=random_double(score);
FOR_LOOP_NEW_VARIABLES_OPEN
CONDITIONS
MFE_SUM
CHILDREN_SUM
INDENT  if((score-= CHILDREN_MAX * MFE_MAX )<=0.0 ){
MFE_BACKTRACE
CHILDREN_BACKTRACE
INDENT    return RETURN_CHECK;
INDENT  }
FOR_LOOP_NEW_VARIABLES_CLOSE
return 0;
}

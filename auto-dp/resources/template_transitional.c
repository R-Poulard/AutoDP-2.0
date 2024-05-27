int compute_MAINNAME(INT_INDICES) {
    // index_start
    if (MAINNAME[index_MAINNAME(INDICES)] > INT_MIN) {
        return MAINNAME[index_MAINNAME(INDICES)];
    }
    // index_end
    
    int min_value = INT_MAX;
    
FOR_LOOP_NEW_VARIABLES_OPEN
    INDENTmin_value = min(min_value, CHILDREN_SUM);
FOR_LOOP_NEW_VARIABLES_CLOSE
    // index_start
    MAINNAME[index_MAINNAME(INDICES)] = min_value;
    // index_end
    return min_value;
}

int bp_score(char x, char y) {
    if (x=='G' && y=='C') { return 10; }
    if (x=='C' && y=='G') { return 10; }
    if (x=='G' && y=='U') { return 5; }
    if (x=='U' && y=='G') { return 5; }
    if (x=='A' && y=='U') { return 5; }
    if (x=='U' && y=='A') { return 5; }
    return 0;
}

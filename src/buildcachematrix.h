
SEXP buildcachematrix(SEXP R_dagdim, SEXP R_dagbanned, SEXP R_dagretained, SEXP R_maxparents, SEXP R_whichnodes);

int isretained(int curnode, int **dagdefn,int **retained_m, int numNodes);


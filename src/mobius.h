
double qprime_prior(int numNodes, int numparents, int priorchoice);

void getAlphaMax(double *ptr_score,int *ptr_nodeid,int **parents,double *alpha,int *alphalookup, int numNodes, int numRows);

int issubset(int **parents,int curcomb,int newcomb,int numNodes);

void getAlphaMaxSingle(int *Sprime, int mynode, double *ptr_score,int *ptr_nodeid,int **parents,double *alpha,int *alphalookup,
		       int numNodes, int numRows, double *myalpha, int *myalphaparentindex, int *start, int *end);

int issubset1D(int *Sprime,int *parents,int numNodes);


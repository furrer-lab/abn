
int generateRandomDAG(gsl_rng *r,network *dag,cache *nodecache,int *store,cycle *cyclestore, network *dagretain, int verbose);
void copynetwork(network *dagsrc,network *dagdest);
void hillSingleIteration(network *dag,cache *nodecache,int *betterscore,network *dag_orig, cycle *cyclestore, network *dag_new, int verbose);
void swapnode(int i,int j,network *dag,cache *nodecache);
int improvedscore(network *dag, network *dag_orig);
int checknumofArcs(network *dag,int maxtotalarcs);

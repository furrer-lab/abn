#include <gsl/gsl_rng.h>
void init_hascycle(cycle *cyclestore, network *dag);
unsigned int hascycle(cycle *cyclestore,network *dag);
void get_numincomingedges(unsigned int *incomingedges,unsigned int **graph, unsigned int numnodes);
void checkandfixcycle(cycle *cyclestore,network *dag, gsl_rng *r, network *dagretain, int verbose);
void droplinks(network *dag, int **retaingraph, unsigned int editnode);




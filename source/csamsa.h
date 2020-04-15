#include <stdio.h>

struct _treenode;
struct _linkedblock;

FILE *debugfile;

int numberofseqs; // number of sequences
char **texts; // list of texts that will originate the tree
int *textsizes; // list of the number of characters of the texts
char **descs; // array of the descriptions of each one of the sequences
int *rotations; // list of rotations for all sequences

// global variables needed for suffix tree functions, so they are not always passed as arguments
char *text; // text of the current sequence
int textsize; // number of characters of the text of the current sequence
int currentseq; // index of the current sequence

void exitMessage(char *msg);

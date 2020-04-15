struct _bordernode;

struct _bordernode *firstbordernode;
struct _bordernode *lastbordernode;
int *startposarray;
int *endposarray;

char CharAt(int pos, int seq);
void PrepareTreeForAlignment();
void RunAlignment();
void SaveAlignment(char *outputfilename);

/*********************************************/
/* GENERALIZED CYCLIC SUFFIX TREES (UKKONEN) */
/*********************************************/

struct _linkedblock;
struct _bordernode;
struct _linkedpos;

// TODO: keep isleaf/leaffrom? (when a rotation of one sequence is a prefix of a rotation of another sequence)
// TODO: cleanup id, isleaf (numbranches==0), leaffrom (numbranches==0&&labelfrom), rotation (=(endpos-depth)%textsizes[labelfrom]))
// TODO: add method GetRotation and structure "suffixfromseqs"
// Suffix Tree Node Structure
typedef struct _treenode {
	int id;							// node identifier
	int startpos;					// position in the text of the first character in the node's label
	int endpos;						// position in the text of the last character in the node's label
	int depth;						// size of the text spelled by the path from the root to this node
	//int count;					// number of times that this suffix occurs in the text
	int isleaf;						// indicates if this is a final node or not
	int rotation;					// number of characters that this text is rotated
	struct _treenode *backlink;		// pointer to the previous above node that links to this node
	struct _treenode *suffixlink;	// pointer to the node whose path label is the suffix of this node's label
	struct _treenode **branches;	// branches of the node
	struct _linkedblock *fromblock;	// pointer to the corresponding block in the largest common unique blocks list
	int numbranches;				// number of branches of the node
	int labelfrom;					// sequence to which the indexes of the label of the node belong
	int leaffrom;					// sequence from which this is a final node
	unsigned long long int fromseqs;	// sequences to which this node belongs
	struct _linkedpos **positions;
	struct _bordernode *bordernode;
} treenode;

treenode *root; // root node of the tree
int steps; // number of steps in the construction of the suffix tree

treenode * buildGeneralizedTree();
void freeTreeNode(treenode *node);
int nodeFromSeq(treenode *node, int i);
int nodeFromAllSeqs(treenode *node);
void markNodeFromSeq(treenode *node, int i);
int getLabelSize(treenode *node);
char *getNodeLabel(treenode *node);
int **getNodePositions(treenode *node, int *numpositions);
treenode * getSuffixNode(treenode *node, int dist, int *pos);
treenode * followChar(treenode *node, int *pos, int charpos);
treenode * followText(char *chars);
treenode * splitNode(treenode *node, int pos);
void printTreeNodeInfo(treenode *currentnode, int pos, int seq, int i);

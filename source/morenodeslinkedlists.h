struct _treenode;

// Linked Position List Structure
typedef struct _linkedpos {
	int k;						// position
	struct _linkedpos *next;	// next position structure
} linkedpos;

// TODO: remove the 'hidden' tag
// Double Linked Bordernode List Structure
typedef struct _bordernode {
	struct _treenode *treenode;		// corresponding node from the suffix tree
	int size;						// size of the substring associated with this node
	struct _linkedpos **positions;	// positions on each sequence where this substring occurs
	int *activeposcount;			// number of positions on each sequence that are lower than the current ending position (in endposarray[i])
	char hidden;					// 't' if this node is hidden and 'f' otherwise
	struct _bordernode *hiddennode;	// pointer to list of hidden nodes inside this node
	struct _linkedpos **hiddenpos;	// pointer to list of hidden positions inside this node (used when node has more than one active position for each sequence)
	struct _bordernode *next;		// pointer to next border node in the list
	struct _bordernode *prev;		// pointer to previous border node in the list
	struct _bordernode *suffix;		// pointer to the border node that corresponds to the longest suffix of this substring
} bordernode;

void ResetBitArrays(struct _treenode *node);
void MarkUsedNodes();
void DeleteUnusedNodes(struct _treenode *node);
void CollectBorderNodes(struct _treenode *node);
void MarkSuffixNodes();
void PrintBorderNodesList(int onlyactivenodes);
void PrintTreeNodeInfoByLabel(char *nodelabel);
int UpdateActiveBorderNodes();
bordernode *NewBorderNode(struct _treenode *treenode);
void ReSortBorderNode(bordernode *node);
void HideFirstPositions(bordernode *node);
void DeleteBorderNodesList();

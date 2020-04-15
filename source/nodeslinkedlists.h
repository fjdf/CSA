struct _treenode;

// Linked Treenode List Structure
typedef struct _linkedblock {
	struct _treenode *item;			// corresponding node from the suffix tree
	int size;						// size of the block chain starting on this block (sum of the sizes/depths of all the blocks)
	int totalsize;					// size of the block chain including intervals (sum of all depths and all intervals in the middle)
	int interval;					// interval between this and the next block (average of all sequences)
	int *positions;					// list of positions of this block on each sequence
	struct _linkedblock *nextblock;	// next block in the block chain
	struct _linkedblock *prev;		// previous similar structure in the linked list
	struct _linkedblock *next;		// next similar structure in the linked list
} linkedblock;

linkedblock *createItem(struct _treenode *item);
linkedblock *addItem(linkedblock *list, struct _treenode *item);
linkedblock *insertSortedItem(linkedblock *list, struct _treenode *item);
linkedblock *sortList(linkedblock *list);
linkedblock *deleteItem(linkedblock *list);
linkedblock *deleteItemAndAdvance(linkedblock *list);
void clearList(linkedblock *list);
char *blockLabel(linkedblock *block);
void printBlocksList(linkedblock *list);
void printChainsList(linkedblock *list);

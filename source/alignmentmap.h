
// Linked Alignment Map Segment List Structure
typedef struct _alignmapsegment {
	int *positions;					// positions on each sequence were the gap to be aligned begins
	int size;						// length of the block correspondent to this segment
	int mingapsize;					// minimum size among all sequences of the gap between this segment and the next
	int maxgapsize;					// maximum size among all sequences of the gap between this segment and the next
	char **alignedstrings;			// strings correspondent to all the aligned sequences
	struct _alignmapsegment *next;	// pointer to next alignment map segment
} alignmapsegment;

// Linked Heaviest Increasing Subsequence Item List Structure
typedef struct _chainitem {
	int *positions;					// position on each sequence were this block occurs
	int size;						// length of this block
	int weight;						// sum of the sizes of all the blocks that constitute the H.I.S. starting in this block
	struct _chainitem *backtrack;	// previous block in the H.I.S.
	struct _chainitem *next;		// next item in the linked list of items
	struct _chainitem *prev;		// previous item in the linked list of items
} chainitem;


alignmapsegment *firstsegment; // fake first alignment segment starting before the begginging of every sequence
alignmapsegment *lastsegment; // fake last alignment segment starting after the ending of every sequence
alignmapsegment *startsegment; // first alignment segment in the current gap being processed
alignmapsegment *endsegment; // last alignment segment in the current gap being processed
chainitem *chain; // linked list of heaviest increasing subsequence chain items

alignmapsegment *NewAlignmentMapSegment(chainitem *chainitem);
void UpdateSegmentGapSizes(alignmapsegment *segment);
void CalculateHeaviestIncreasingSubsequence();
int SetAlignmentMapSegments();
void DeleteAlignmentMap();
void PrintChain();
void PrintAlignmentMap(int onlyactivesegments);
void PrintAlignmentSegmentStrings(alignmapsegment *segment);

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/timeb.h>
#include <limits.h>
#include "csamsa.h"
#include "gencycsuffixtrees.h"
#include "nodeslinkedlists.h"
#include "graphics.h"
#include "alignment.h"
#include "console.h"
#include "tools.h"

#include "bitmap.h"

#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif

#define PAUSE_AT_EXIT 1

#define MAXNUMBEROFSEQS 64
#define DEBUG_FILENAME "debug.txt"
#define POSITIONS_FILENAME "-positions.txt"
#define IMAGEMAP_FILENAME "-imagemap.txt"
#define ROTATIONS_FILENAME "-Rotated.fasta"
#define ALIGNMENT_FILENAME "-Aligned.fasta"
#define BLOCKSINFO_FILENAME "-Blocks.csv"
#define BLOCKSIMAGE_FILENAME "-Blocks.bmp"
#define CIRCULARIMAGE_FILENAME "-CircularAlignment.bmp"

int minblocksize; // minimum block size
int maxblocksize; // maximum block size
int maxinterval; // maximum allowed interval length between two consecutive blocks in a chain
struct _linkedblock *blockslist; // list of nodes of maximum depth that belong to all the sequences
struct _linkedblock **positionslist; // list of the final positions on each sequence for the maximum depth nodes
int *countslist; // number of final positions for each list
int mcscount; // number of maximum unique common sequences
char *inputfilename; // name of the input file containing the original sequences


// Joins the basename of the input file with another string
char* newOutputFilename(char *extra){
	char *resultfilename;
	int i, n;
	n=(int)strlen(inputfilename);
	for(i=(n-1);i>0;i--){
		if(inputfilename[i]=='.') break;
	}
	if(i==0) i=n;
	n=(int)strlen(extra);
	resultfilename=(char *)calloc((i+n+1),sizeof(char));
	strncpy(resultfilename,inputfilename,i);
	resultfilename[i]='\0';
	strcat(resultfilename,extra);
	return resultfilename;
}

void exitMessage(char *msg){
	printf("\n> ERROR: %s\n",msg);
	if(PAUSE_AT_EXIT) getchar();
	exit(0);
}


// TODO: *** only insert the nodes and sort them in the end (find largest depth, array of pointers with that size, group nodes with same depth, link them together)
// follows all branches of the tree (recursively) and collects in the blockslist array the deepest nodes that belong to all sequences
int collectNodes(treenode *node){
	int i, continues;
	if(!nodeFromAllSeqs(node)) return 0;
	continues=0;
	for(i=0;i<(node->numbranches);i++)
		continues = continues | collectNodes(node->branches[i]);
	if(!continues){
		blockslist=insertSortedItem(blockslist,node);
		//blockslist=addItem(blockslist,node);
		mcscount++;
	}
	return 1;
}

// TODO: implement more sofisticated way of getting rid of the suffixes (fill suffixfrom structure when building tree)
// removes from the blockslist array all the nodes that are suffixes of another node
void removeSuffixNodes(){
	linkedblock *node,*searchnode;
	treenode *suffix;
	int labelsize,pos;
	if(blockslist==NULL || (blockslist->item)==root) return;
	node=blockslist; // the blockslist is sorted by decreasing depth
	while(node!=NULL){ // for each block, get the corresponding suffixes and delete them if they appear ahead in the blocks list
		labelsize=((node->item->endpos)-(node->item->startpos)+1);
		pos=labelsize; // the position inside the tree node is the end of the label
		suffix=getSuffixNode(node->item,labelsize,&pos); // get first suffix
		searchnode=(node->next);
		while(searchnode!=NULL && suffix!=root){
			//printf("#%d[%d-%d](%d)@%d - %d/%d\n",(suffix->id),(suffix->startpos),(suffix->endpos),(suffix->isleaf),(suffix->depth),pos,((suffix->endpos)-(suffix->startpos)+1));
			while(searchnode!=NULL && (suffix->depth)<(searchnode->item->depth)) // advance in the blocks list until the blocks have the same depth as the suffix
				searchnode=(searchnode->next);
			while(searchnode!=NULL && (suffix->depth)==(searchnode->item->depth)){ // check if any of the blocks at that depth have that suffix, and delete it if so
				if(suffix==(searchnode->item)) { searchnode=deleteItem(searchnode); mcscount--; break; }
				searchnode=(searchnode->next);
			}
			labelsize=((suffix->endpos)-(suffix->startpos)+1);
			suffix=getSuffixNode(suffix,labelsize,&pos); // get next suffix
		} // suffixes loop
		node=(node->next);
	} // blocks list loop
}


// TODO: implement non-recursive collectPositions and remove positionslist and countslist variables
// collects the positions and the number of positions of a specific node in the positionslist and countslist arrays
void collectPositions(treenode *node){
	int i,n;
	if((node->isleaf)){
		n=node->leaffrom;
		positionslist[n]=insertSortedItem(positionslist[n],node);
		countslist[n]++;
	}
	n=node->numbranches;
	for(i=0;i<n;i++) collectPositions(node->branches[i]);
}

// TODO: *** reformulate method (sort blocks by position in shortest sequence, connect chains + 1st block again, sort chains by size)
// TODO: remove "totalsize" and only keep size
// TODO: check if all blocks in same order on all sequences gives problems
// NOTE: in the end, the blocks that are at the beggining of a chain have totalsize>0 and all other blocks have totalsize=-1
// NOTE: starting blocks have size set to the sum of all blocks in the chain, and the remaining blocks have size set to their depth/length
// NOTE: size refers to the sum of lengths/depths excluding intervals/gaps, and totalsize includes intervals/gaps length
// links the blocks that appear in the same order in all the sequences
void collectNodeChains(){
	treenode *currentnode,*foundnode;
	linkedblock *block,*prevblock,*currblock;
	int i,k,n,pos,count,interval;
	if(blockslist==NULL) return;
	block=blockslist;
	while(block!=NULL){ // mark all valid nodes in the tree
		block->totalsize=0; // block->size=0 already too
		block->item->fromblock=block;
		block=(block->next);
	}
	for(k=0;k<numberofseqs;k++){ // follow the whole sequence in the tree and link the blocks in the order that they appear in the sequence
		text=texts[k];
		textsize=textsizes[k];
		n=textsize;
		currentnode=root;
		pos=1;
		prevblock=NULL;
		for(i=0;i<n;i++){
			//printTreeNodeInfo(currentnode,pos,k,i);
			foundnode=followChar(currentnode,&pos,i); // follow char in tree
			if(!nodeFromAllSeqs(foundnode)){ // if it does not belong to all sequences, analyze last scanned common node
				block=(currentnode->fromblock);
				if(block!=NULL){ // if valid block
					if(prevblock!=NULL){ // if is not first block found
						if((prevblock->size)==0){ // if is still valid
							if((prevblock->nextblock)==NULL) prevblock->nextblock=block; // first time linking to a block
							else if((prevblock->nextblock)!=block){ // already linked to a block (by one of previous seqs) but not the same
								prevblock->nextblock=NULL; // invalidate previously defined block chain
								prevblock->size=-1;
							}
						}
					} else n+=(i-(currentnode->depth)); // set where to stop following chars (before first block) to catch blocks in the "wrap around" part
					prevblock=block; // current block is the next previous block
				}
				if(currentnode!=root){
					currentnode=currentnode->suffixlink; // follow suffix link
					pos=getLabelSize(currentnode); // go to end of suffix node
					i--; // follow the same character again
				}
			} else { // if it belongs to all sequences, continue to next node below in the tree
				currentnode=foundnode;
				pos=getLabelSize(currentnode); // go to end of node
				i+=(pos-1); // skip the number of chars inside the node
			}
		} // end textsize loop
	} // end sequences loop

	i=0;
	block=blockslist; // at the beggining, block->totalsize=0 for all blocks
	while(block!=NULL){
		if((block->totalsize)==-1){ // if we are starting with a block that already appeared in the middle of another chain
			block=(block->next);
			continue;
		}
		block->size=(block->item->depth); // all blocks have their size set to their length
		prevblock=block;
		currblock=(block->nextblock);
		while(currblock!=NULL){ // process all blocks in the chain starting in this one
			interval=INT_MAX; // check if the intervals between them on every sequence do not exceed the maximum bound
			for(k=0;k<numberofseqs;k++){ // calculate shortest interval between the 2 blocks on all sequences
				count=0;
				if((currblock->positions[k])<(prevblock->positions[k])) count+=textsizes[k]; // if the current block has a position lower than the previous block, update it to go over the size of the text
				count+=((currblock->positions[k])-((prevblock->positions[k])+(prevblock->item->depth))); // distance between end of previous block and beginning of current block
				if(count<interval) interval=count;
			}
			if(interval>maxinterval){ // if the interval is larger than the maximum size allowed, break this chain here (initialized with INT_MAX)
				prevblock->nextblock=NULL;
				break;
			}
			if(currblock->totalsize>0){ // if this block was already processed before as the beggining of a chain, now it is in the middle of a larger one
				block->size+=(currblock->size); // update the size of the chain (in the first block) with the size of this existing chain
				block->totalsize+=(currblock->totalsize); // update the size including the intervals too
				prevblock->interval=interval;
				block->totalsize+=interval; // add this intervals length to the total size of the chain
				currblock->size=(currblock->item->depth); // reset the size of this block/chain to the size of the first block
				currblock->totalsize=-1; // mark this block as already processed
				mcscount--; // one less chain
				break; // stop following the chain until the end (it was already processed)
			}
			currblock->size=(currblock->item->depth); // this block was not processed before, so, set its size
			block->size+=(currblock->size); // update the size of the chain (in the first block) with the size of this block
			prevblock->interval=interval;
			block->totalsize+=interval; // add this intervals length to the total size of the chain
			currblock->totalsize=-1; // mark this block as already processed
			mcscount--; // this block is part of an existing chain, so, it does not count
			prevblock=currblock;
			currblock=(currblock->nextblock); // process next block in the chain
		}
		block->totalsize+=block->size; // sum blocks length and intervals length
		block=(block->next); // next block in the list of all blocks
		i++;
	}
	blockslist=sortList(blockslist);
}

// TODO: non-recursive collectPositions (fixed array, if slot already filled then stop)
// collects the positions and removes the nodes that occur more than once on at least one sequence
void removeNonUniqueNodes(){
	linkedblock *node;
	int i,k;
	positionslist=(linkedblock **)calloc(numberofseqs,sizeof(linkedblock *));
	countslist=(int *)calloc(numberofseqs,sizeof(int));
	if(blockslist==NULL) return;
	node=blockslist;
	while(node!=NULL){ // finds the nodes that appears only once on each sequence and stores their positions
		collectPositions((node->item)); // fills positionslist and countslist arrays for this node
		for(k=0;k<numberofseqs;k++){ // check if this subsequence appears in each sequence only once
			if(countslist[k]>1) break;
		}
		if(k==numberofseqs){ // if check is ok, fill positions array
			node->positions=(int *)calloc(numberofseqs,sizeof(int));
			for(i=0;i<numberofseqs;i++) node->positions[i]=positionslist[i]->item->rotation;
			node=(node->next);
		} else {
			if(node==blockslist) blockslist=(node->next); // if we are deleting the first node of the list, do not lose the list pointer
			node=deleteItemAndAdvance(node); // if node is not unique, delete it and advance to next node
			mcscount--;
		}
		for(i=0;i<numberofseqs;i++) { // clear positionslist and countslist arrays
			clearList(positionslist[i]);
			positionslist[i]=NULL;
			countslist[i]=0;
		}
	}
}

// identifies the start positions in all the sequences of the largest unique common subsequence
void getRotations(){
	linkedblock *node;
	int i;
	if(blockslist==NULL) return;
	node=blockslist; // the largest chain is at the beginning of the sorted list
	rotations=(int *)calloc(numberofseqs,sizeof(int));
	for(i=0;i<numberofseqs;i++) rotations[i]=node->positions[i];
}

// TODO: ??? after calculating rotations update all tree nodes with the new positions
// analyzes the tree in order to calculate the common subsequences and the correct rotation
void analyzeTree(){
	blockslist=NULL;
	mcscount=0;
	printf("> Collecting maximum common subsequences... ");fflush(stdout);
	collectNodes(root); // collect subsequences that belong to all sequences
	if(mcscount==0){
		exitMessage("No common subsequences found");
	}
	printf("%d nodes found\n",mcscount);
	#ifdef DEBUG
	printBlocksList(blockslist);
	#endif
	printf("> Removing suffixes... ");fflush(stdout);
	removeSuffixNodes(); // remove subsequences that are suffixes of other subsequence
	printf("%d nodes left\n",mcscount);
	#ifdef DEBUG
	printBlocksList(blockslist);
	#endif
	printf("> Removing repeats... ");fflush(stdout);
	removeNonUniqueNodes();
	if(mcscount==0){
		exitMessage("No unique subsequences found");
	}
	printf("%d nodes left\n",mcscount);
	#ifdef DEBUG
	printBlocksList(blockslist);
	#endif
	printf("> Connecting block chains... ");fflush(stdout);
	collectNodeChains();
	printf("%d chains found\n",mcscount);
	#ifdef DEBUG
	printChainsList(blockslist);
	#endif
	getRotations(); // collect rotations for each sequence for the largest unique common subsequence
	if(rotations==NULL){
		exitMessage("No unique common subsequences found");
	}
}

void createImageAndShowResults(char *imagefilename, int showchainsonly, int withrotation){
	FILE *originalblocksfile;
	FILE *datafile;
	linkedblock *node;
	int ntoprint=20;
	int charstoprint=100;
	int i,size,ndrawn;
	char *string;
	int *rgbarray, *rotatedpositions;
	int *nullrotations=NULL;
	char *tmpfilename;
	if(blockslist==NULL) return;
	tmpfilename=newOutputFilename(POSITIONS_FILENAME);
	datafile=fopen(tmpfilename,"w");
	if(datafile==NULL) {
		exitMessage("Can't write data file");
	}
	free(tmpfilename);
	fprintf(datafile,"%d\n",numberofseqs);
	tmpfilename=newOutputFilename(BLOCKSINFO_FILENAME);
	originalblocksfile=fopen(tmpfilename,"w");
	if(originalblocksfile==NULL) {
		exitMessage("Can't write original blocks file");
	}
	free(tmpfilename);
	fprintf(originalblocksfile,"Length,Sequence");
	for(i=0;i<numberofseqs;i++)	fprintf(originalblocksfile,",Position_%d",(i+1));
	fprintf(originalblocksfile,"\n");
	tmpfilename=newOutputFilename(IMAGEMAP_FILENAME);
	rgbarray=(int *)calloc(3,sizeof(int));
	rotatedpositions=(int *)calloc(numberofseqs,sizeof(int));
	if(withrotation){
		initializeBlocks(numberofseqs,textsizes,rotations,tmpfilename); // true rotations
	} else {
		nullrotations=(int *)calloc(numberofseqs,sizeof(int));
		for(i=0;i<numberofseqs;i++) nullrotations[i]=0;
		initializeBlocks(numberofseqs,textsizes,nullrotations,tmpfilename); // all rotations set to zero
	}
	free(tmpfilename);
	node=blockslist;
	mcscount=0;
	ndrawn=0;
	printf("> Length, sequence and rotations for the first %d longest block chains:\n",ntoprint);
	while(node!=NULL){
		if(showchainsonly) size=(node->totalsize); // size of the chain starting at this block (sum of all individual blocks and intervals)
		else size=(node->item->depth); // size of this single block (not including other blocks in chain)
		if(size>0){ // only draw this block in the image if we want all individual blocks or if this is the beggining of a chain (totalsize!=-1)
			if(size>=minblocksize && size<=maxblocksize){ // check if the block has the correct size to be shown in image
				for(i=0;i<numberofseqs;i++)
					rotatedpositions[i]=drawBlockRotated(node->positions[i],size,i);
				getRGBColor(rgbarray);
				fprintf(datafile,"%d %d %d %d",rgbarray[0],rgbarray[1],rgbarray[2],size);
				for(i=0;i<numberofseqs;i++) fprintf(datafile," %d",rotatedpositions[i]);
				fprintf(datafile,"\n");
				connectBlocks();
				ndrawn++;
			}
		}
		size=(node->totalsize);
		if(size==-1){ // draw all single blocks in the image, but only print info about the chains
			node=(node->next);
			continue;
		}
		string=blockLabel(node);
		if(mcscount<ntoprint){
			printf(":: (%d) ",(node->size)); // only print size of blocks without interval/gap lengths
			if(((int)strlen(string))<charstoprint) printf("%s",string);
			else{
				for(i=0;i<charstoprint;i++) printf("%c",string[i]);
				printf("...");
			}
			/*
			printf(" [");
			for(i=0;i<numberofseqs;i++)
				printf("%d,",(node->positions[i]));
			ConsoleMoveCursorPosition(-1,0);
			printf("]");
			*/
			printf("\n");
			fflush(stdout);
		}
		fprintf(originalblocksfile,"%d,%s",size,string); // print totalsize including gaps to file
		for(i=0;i<numberofseqs;i++)
			fprintf(originalblocksfile,",%d",(node->positions[i]));
		fprintf(originalblocksfile,"\n");
		free(string);
		mcscount++; // total number of chains
		node=(node->next); // next node in list
	}
	if(mcscount>ntoprint) printf(":: ... (%d total)\n",mcscount);
	fclose(datafile);
	fclose(originalblocksfile);
	drawLabels(descs);
	string=(char *)calloc(255,sizeof(char));
	if(maxblocksize==INT_MAX && minblocksize==1) sprintf(string,"%d chain blocks",mcscount); // default values
	else if(maxblocksize==INT_MAX && minblocksize!=1) sprintf(string,"%d %s with size >=%d of a total of %d block chains",ndrawn,((showchainsonly)?"chains":"blocks"),minblocksize,mcscount);
	else if(maxblocksize!=INT_MAX && minblocksize==1) sprintf(string,"%d %s with size <=%d of a total of %d block chains",ndrawn,((showchainsonly)?"chains":"blocks"),maxblocksize,mcscount);
	else sprintf(string,"%d %s with size >=%d and <=%d of a total of %d block chains",ndrawn,((showchainsonly)?"chains":"blocks"),minblocksize,maxblocksize,mcscount);
	drawBottomLabel(string);
	free(rotatedpositions);
	free(rgbarray);
	free(string);
	if(nullrotations!=NULL) free(nullrotations);
	finalizeGraphics(imagefilename);
}

void saveRotatedSequences(char *outputfilename){
	FILE *file;
	int i,rot;
	file=fopen(outputfilename,"w");
	if(file==NULL) {
		exitMessage("Can't write rotated sequences file");
	}
	for(i=0;i<numberofseqs;i++){
		rot=rotations[i];
		fprintf(file,">%s @ %d\n",descs[i],rotations[i]);
		fputs((texts[i]+rot),file);
		fwrite(texts[i],sizeof(char),(size_t)rot,file);
		fprintf(file,"\n");
	}
	fclose(file);
}

void LoadSequences(){
	FILE *file;
	char c;
	int k,totalsize,desclen,seqlen;
	long int seqstart,seqend,seqsize;
	fpos_t startpos;
	printf("> Loading sequences from file <%s> ... ",inputfilename);
	if((file=fopen(inputfilename,"r"))==NULL){
		exitMessage("Sequence file not found");
	} else {
		seqstart=ftell(file);
		fseek(file,0L,SEEK_END);
		seqend=ftell(file);
		seqsize=(seqend-seqstart);
		rewind(file);
		printf("(%ld bytes)\n",seqsize);
	}
	texts=(char **)calloc(MAXNUMBEROFSEQS,sizeof(char *));
	descs=(char **)calloc(MAXNUMBEROFSEQS,sizeof(char *));
	textsizes=(int *)calloc(MAXNUMBEROFSEQS,sizeof(int));
	numberofseqs=0;
	totalsize=0;
	c=fgetc(file);
	while(c!=EOF && c!='>') c=fgetc(file);
	if(c==EOF){
		exitMessage("No sequences in file");
	}
	while(1){
		while(c!=EOF && c!='>') c=fgetc(file);
		if(c==EOF) break;
		desclen=0;
		fgetpos(file,&startpos);
		while((c=fgetc(file))!=EOF && c!='\n' && c!='\r') desclen++;
		descs[numberofseqs]=(char *)calloc((desclen+1),sizeof(char));
		fsetpos(file,&startpos);
		k=0;
		while((c=fgetc(file))!=EOF && c!='\n' && c!='\r') descs[numberofseqs][k++]=c;
		printf("# %02d [",(numberofseqs+1));
		k=0;
		while(k<40 && k<desclen) printf("%c",descs[numberofseqs][k++]);
		while(k<40) { printf(" "); k++; }
		printf("] ");
		fflush(stdout);
		seqlen=0;
		fgetpos(file,&startpos);
		while((c=fgetc(file))!=EOF && c!='>') seqlen++;
		texts[numberofseqs]=(char *)calloc(((int)seqlen+1),sizeof(char));
		fsetpos(file,&startpos);
		k=0;
		while((c=fgetc(file))!=EOF && c!='>'){
			if(c=='\n' || c=='\r' || c=='\0' || c=='-' || c==' ') continue;
			if(c=='A' || c=='C' || c=='G' || c=='T') texts[numberofseqs][k++]=c;
			else if(c=='a' || c=='c' || c=='g' || c=='t') texts[numberofseqs][k++]=(char)(c-32);
			else if(c=='R' || c=='Y' || c=='S' || c=='W' || c=='K' || c=='M' ||
					c=='D' || c=='H' || c=='B' || c=='V' || c=='N') texts[numberofseqs][k++]=c;
			else if(c=='r' || c=='y' || c=='s' || c=='w' || c=='k' || c=='m' ||
					c=='d' || c=='h' || c=='b' || c=='v' || c=='n') texts[numberofseqs][k++]=(char)(c-32);
			else break;
		}
		if(k==0){
			printf("EMPTY\n");
			free(descs[numberofseqs]);
			free(texts[numberofseqs]);
			continue;
		}
		if(c!=EOF && c!='>'){
			printf("INVALID_CHARS\n");
			free(descs[numberofseqs]);
			free(texts[numberofseqs]);
			continue;
		}
		printf("OK (%d characters)\n",k);
		fflush(stdout);
		textsizes[numberofseqs]=k;
		totalsize+=k;
		numberofseqs++;
		if(numberofseqs==MAXNUMBEROFSEQS){
			printf("> WARNING: Current version only supports up to %d sequences\n",MAXNUMBEROFSEQS);
			break;
		}
	}
	fclose(file);
	if(numberofseqs<2){
		exitMessage("Not enough valid sequences found");
	}
	printf("> %d sequences successfully loaded\n",numberofseqs);
}

/********/
/* MAIN */
/********/
int main(int argc, char *argv[]) {
	int i;
	treenode *tree;
	//struct timeb starttb,endtb;
	//double timetb;
	char *rotationsfilename;
	char *alignmentfilename;
	char *blocksimagefilename;
	char *circularimagefilename;
	char mode;
	ConsoleSetTextColor(COLOR_BRIGHT,COLOR_RED,COLOR_WHITE);
	printf("[ Multiple Circular Sequence Aligner v1.11 ]");
	ConsoleResetTextColor();
	printf("\n");
	mode='\0';
	if(argc==2){
		mode='N'; // normal mode
		inputfilename=argv[1];
	}
	if(argc==3){
		mode=argv[1][0]; // get mode option char
		if(mode>='a' && mode<='z') mode=(char)(mode-32); // convert to uppercase if needed
		if(mode!='R' && mode!='A' && mode!='I' && mode!='C' && mode!='S' && mode!='M') mode='\0'; // check valid mode chars
		inputfilename=argv[2];
	}
	if( mode=='\0' ){
		printf("> USAGE:\n");
		printf("\t[ Rotate+Align+Image ]\t%s <multi-fasta-file>\n",argv[0]);
		printf("\t[   Rotation only    ]\t%s R <multi-fasta-file>\n",argv[0]);
		printf("\t[   Alignment only   ]\t%s A <multi-fasta-file>\n",argv[0]);
		printf("\t[Alignment Image only]\t%s I <multi-fasta-file>\n",argv[0]);
		printf("> TOOLS:\n");
		printf("\t[  Clean FASTA file  ]\t%s C <multi-fasta-file>\n",argv[0]);
		printf("\t[Get Alignment Score ]\t%s S <multi-fasta-file>\n",argv[0]);
		printf("\t[Convert FASTA to MSF]\t%s M <multi-fasta-file>\n",argv[0]);
		/*
		printf("> OPTIONS:\n");
		printf("\t-M <max-interval>   : maximum interval between two blocks in a chain (optional;default=16)\n");
		printf("\t-S <min-block-size> : shortest block size to show in blocks image    (optional;default= 0)\n");
		printf("\t-W <max-block-size> : widest block size to show in blocks image      (optional;default= %c)\n",((char)236));
		printf("\t-L                  : linear alignment image instead of circular     (optional;default= -)\n");
		printf("\t-F                  : alignment in MSF format instead of FASTA       (optional;default= -)\n");
		printf("\t-A                  : use alternative scoring function               (optional;default= -)\n");
		*/
		printf("> Done!\n");
		return 0;
	}

	// initialize and parse other arguments
	minblocksize=10;
	maxblocksize=INT_MAX;
	maxinterval=INT_MAX;
	rotations=NULL;
	rotationsfilename=NULL;
	alignmentfilename=NULL;
	blocksimagefilename=NULL;
	circularimagefilename=NULL;
	debugfile=NULL;
	tree=NULL;

	#if defined(DEBUG) || defined(DEBUGDP) || defined(DEBUGSORT)
	debugfile=fopen(DEBUG_FILENAME,"w");
	if(debugfile==NULL) {
		exitMessage("Can't write debug file");
	}
	#endif

	// load the sequences and build the cyclic suffix tree
	if( mode=='N' || mode=='R' || mode=='A' ){
		LoadSequences();
		printf("> Building generalized%s suffix tree",(mode=='A')?"":" cyclic");
		fflush(stdout);
		steps=0;
		//ftime(&starttb);
		tree=buildGeneralizedTree();
		printf("\n");
		//checkSuffixTree();
		//ftime(&endtb);
		//timetb=((endtb.time) + (endtb.millitm)/1000.0) - ((starttb.time) + (starttb.millitm)/1000.0);
		//printf("> Statistics: %d sequences ; %d chars ; %d blocks ; %.3lf seconds\n",numberofseqs,totalsize,mcscount,timetb);
	}

	// find the best rotation for the sequences
	if( mode=='N' || mode=='R' ){
		rotationsfilename=newOutputFilename(ROTATIONS_FILENAME);
		blocksimagefilename=newOutputFilename(BLOCKSIMAGE_FILENAME);
		analyzeTree();
		saveRotatedSequences(rotationsfilename);
		createImageAndShowResults(blocksimagefilename,1,1);
	}

	// do the multiple sequence alignment
	if( mode=='N' || mode=='A' ){
		if(mode=='A') rotations=(int *)calloc(numberofseqs,sizeof(int));
		alignmentfilename=newOutputFilename(ALIGNMENT_FILENAME);
		PrepareTreeForAlignment();
		printf("> Running multiple sequence alignment...\n");
		fflush(stdout);
		RunAlignment();
		SaveAlignment(alignmentfilename);
		TestAlignmentFileOutput(((rotationsfilename!=NULL)?(rotationsfilename):(inputfilename)),alignmentfilename);
	}

	// draw the circular alignment plot
	if( mode=='N' || mode=='I' ){
		circularimagefilename=newOutputFilename(CIRCULARIMAGE_FILENAME);
		DrawCircularAlignmentPlot(((alignmentfilename!=NULL)?(alignmentfilename):(inputfilename)),circularimagefilename);
	}

	// clean fasta file
	if(mode=='C'){
		CleanDNAFastaFile(inputfilename);
	}

	// calculate alignment score
	if(mode=='S'){
		CalculateSumOfPairsScore(inputfilename);
	}

	// convert fasta file to msf file
	if(mode=='M'){
		ConvertFastaToMsf(inputfilename);
	}

	if(debugfile!=NULL) fclose(debugfile);

	printf("> Done!\n");
	fflush(stdout);
	if(PAUSE_AT_EXIT) getchar();

	if(rotationsfilename!=NULL) free(rotationsfilename);
	if(alignmentfilename!=NULL) free(alignmentfilename);
	if(blocksimagefilename!=NULL) free(blocksimagefilename);
	if(circularimagefilename!=NULL) free(circularimagefilename);

	if(root!=NULL) freeTreeNode(tree);
	if(texts!=NULL || descs!=NULL){
		for(i=0;i<numberofseqs;i++){
			free(texts[i]);
			free(descs[i]);
		}
		free(descs);
		free(texts);
	}
	if(textsizes!=NULL) free(textsizes);
	//free(masks);
	if(positionslist!=NULL){
		for(i=0;i<numberofseqs;i++) clearList(positionslist[i]);
		free(positionslist);
	}
	if(countslist!=NULL) free(countslist);
	if(rotations!=NULL) free(rotations);
	if(blockslist!=NULL) clearList(blockslist);
	return 0;
}

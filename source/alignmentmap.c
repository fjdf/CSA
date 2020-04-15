#include <stdlib.h>
#include <limits.h>
#include "csamsa.h"
#include "alignmentmap.h"
#include "alignment.h"
#include "morenodeslinkedlists.h"

// Creates a new chain item starting at the first position of each sequence in the bordernode
chainitem *NewChainItem(bordernode *bordernode){
	chainitem *newchainitem;
	int i,pos,actualsize,newsize,auxsize;
	newchainitem=(chainitem *)malloc(sizeof(chainitem));
	newchainitem->positions=(int *)malloc(numberofseqs*sizeof(int));
	actualsize=bordernode->size;
	newsize=actualsize;
	for(i=0;i<numberofseqs;i++){ // check if the ending position of this chain item on each sequence does not go beyond the current ending positions
		pos=bordernode->positions[i]->k;
		newchainitem->positions[i]=pos;
		auxsize=(pos+actualsize);
		if(auxsize>=endposarray[i]){
			auxsize=(endposarray[i]-pos);
			if(auxsize<newsize) newsize=auxsize; // update the chain item size if needed
		}
	}
	newchainitem->size=newsize;
	newchainitem->weight=newsize;
	newchainitem->backtrack=NULL;
	newchainitem->next=NULL;
	newchainitem->prev=NULL;
	return newchainitem;
}

// Deletes all the items in the chain
void DeleteChain(){
	chainitem *chainitem,*nextchainitem;
	chainitem=chain;
	while(chainitem!=NULL){
		nextchainitem=chainitem->next;
		if(chainitem->positions!=NULL) free(chainitem->positions);
		free(chainitem);
		chainitem=nextchainitem;
	}
	chain=NULL;
}

void PrintChain(){
	int i,n;
	chainitem *chainitem, *backitem;
	chainitem=chain;
	n=0;
	fprintf(debugfile,"\n");
	while(chainitem!=NULL){
		fprintf(debugfile,"->[p=(");
		for(i=0;i<numberofseqs;i++){
			if(i!=0) fprintf(debugfile,",");
			fprintf(debugfile,"%d",(chainitem->positions[i])); // first position
		}
		fprintf(debugfile,")-(");
		for(i=0;i<numberofseqs;i++){
			if(i!=0) fprintf(debugfile,",");
			fprintf(debugfile,"%d",((chainitem->positions[i])+(chainitem->size)-1)); // last position
		}
		fprintf(debugfile,"),w=%d,s=%d,\"",(chainitem->weight),(chainitem->size));
		for(i=0;i<(chainitem->size);i++){
			fprintf(debugfile,"%c",CharAt((chainitem->positions[0]+i),0)); // string
		}
		fprintf(debugfile,"\"]");
		backitem=(chainitem->backtrack); // backtrack node
		if(backitem!=NULL){
			fprintf(debugfile," ^[p=(");
			for(i=0;i<numberofseqs;i++){
				if(i!=0) fprintf(debugfile,",");
				fprintf(debugfile,"%d",(backitem->positions[i]));
			}
			fprintf(debugfile,")]");
		}
		fprintf(debugfile,"\n");
		chainitem=chainitem->next;
		n++;
	}
	fprintf(debugfile,"(%d chain items)\n\n",n);
	fflush(debugfile);
}

// Check if the chainitem1 is lower than chainitem2 in all sequences (and that they do not overlap)
int LowerThan(chainitem *chainitem1, chainitem *chainitem2){
	int i;
	for(i=0;i<numberofseqs;i++){ // the ending pos of chainitem1 should always be lower than the starting pos of chainitem2: (pos1+size)<=(pos2)
		if(((chainitem1->positions[i])+(chainitem1->size))>(chainitem2->positions[i]))
			return 0;
	}
	return 1;
}

// Check if the chainitem1 is higher than chainitem2 in all sequences (and that they do not overlap)
int GreaterThan(chainitem *chainitem1, chainitem *chainitem2){
	int i;
	for(i=0;i<numberofseqs;i++){   // the starting pos of chainitem1 should always be higher than the ending pos of chainitem2: (pos1)>=(pos2+size)
		if((chainitem1->positions[i])<((chainitem2->positions[i])+(chainitem2->size)))
			return 0;
	}
	return 1;
}

// Chains all the active nodes from the border nodes list, creating a weight-sorted chain corresponding to the Heaviest Increasing Subsequence
// NOTE: the chain is always sorted in decreasing weight order (the first item is the highest/heaviest)
void CalculateHeaviestIncreasingSubsequence(){
	bordernode *bordernode, *nextnode;
	chainitem *newitem,*currentitem,*nextitem,*previtem;
	chain=NULL;
	bordernode=firstbordernode->next;
	while(bordernode!=NULL && (bordernode->positions[0]->k)<endposarray[0]){ // only process active border nodes
		newitem=NewChainItem(bordernode); // create a new chain item for each node
		currentitem=NULL;
		nextitem=chain; // find the first (heaviest) item in the chain that is lower (in positions) than this new item
		while(nextitem!=NULL && !GreaterThan(newitem,nextitem)){ // stop when the new item is greater (in positions) than the next item in the chain
			currentitem=nextitem;
			nextitem=currentitem->next;
		}
		if(nextitem!=NULL){
			newitem->weight+=nextitem->weight; // add the weight of the next (lower in positions) item ahead
			newitem->backtrack=nextitem; // link it to backtrack later
		}
		previtem=currentitem;
		currentitem=nextitem; // find (backwards) the right place in the weight-sorted chain to place the new item
		while(previtem!=NULL && (newitem->weight)>=(previtem->weight)){ // while the weight of the new item is higher than the weight of the current chain item, go backwards
			currentitem=previtem;
			previtem=currentitem->prev;
		}
		if(previtem==NULL) chain=newitem;
		else previtem->next=newitem;
		newitem->prev=previtem;
		if(currentitem!=NULL) currentitem->prev=newitem;
		newitem->next=currentitem;
		nextnode=bordernode->next;
		if(bordernode->activeposcount[0]>1){ // if the node has more positions available, hide the first positions (already processed) and re-sort (by position) the node in the nodes list
			HideFirstPositions(bordernode);
			ReSortBorderNode(bordernode);
			if((bordernode->next)==nextnode) nextnode=bordernode; // if after sorting, the position of the node did not change, process this same node again
		}
		bordernode=nextnode;
	}
}

// Creates a new alignment map segment with the same size and positions of a chain item
// NOTE: the list of positions is moved from the chain item to the alignment map segment
alignmapsegment *NewAlignmentMapSegment(chainitem *chainitem){
	alignmapsegment *newalignmapsegment;
	newalignmapsegment=(alignmapsegment *)malloc(sizeof(alignmapsegment));
	newalignmapsegment->positions=NULL;
	newalignmapsegment->mingapsize=INT_MAX;
	newalignmapsegment->maxgapsize=INT_MAX;
	newalignmapsegment->alignedstrings=NULL;
	newalignmapsegment->next=NULL;
	if(chainitem!=NULL){
		newalignmapsegment->positions=chainitem->positions;
		chainitem->positions=NULL;
		newalignmapsegment->size=chainitem->size;
	} else {
		newalignmapsegment->positions=(int *)calloc(numberofseqs,sizeof(int));
		newalignmapsegment->size=0;
	}
	return newalignmapsegment;
}

// Deletes all alignment map segments from the alignment map list
void DeleteAlignmentMap(){
	alignmapsegment *segment,*nextsegment;
	int i;
	segment=firstsegment;
	while(segment!=NULL){
		nextsegment=segment->next;
		if(segment->positions!=NULL) free(segment->positions);
		if(segment->alignedstrings!=NULL){
			for(i=0;i<numberofseqs;i++){
				if(segment->alignedstrings[i]!=NULL) free(segment->alignedstrings[i]);
			}
			free(segment->alignedstrings);
		}
		free(segment);
		segment=nextsegment;
	}
	firstsegment=NULL;
	lastsegment=NULL;
}

void PrintAlignmentMap(int onlyactivesegments){
	int i, count, activecount;
	alignmapsegment *segment;
	segment=firstsegment;
	count=0;
	activecount=0;
	fprintf(debugfile,"\n");
	while(segment!=NULL){
		if(onlyactivesegments){
			for(i=0;i<numberofseqs;i++){
				if( (segment->positions[i])<startposarray[i] || (segment->positions[i])>=endposarray[i] ) break;
			}
			if(i==numberofseqs){ // it is an active segment
				fprintf(debugfile,"*");
				activecount++;
			}
		}
		if(segment==startsegment || segment==endsegment) fprintf(debugfile,"=>");
		else fprintf(debugfile,"->");
		fprintf(debugfile,"[p=(");
		for(i=0;i<numberofseqs;i++){
			if(i!=0) fprintf(debugfile,",");
			fprintf(debugfile,"%d",(segment->positions[i]));
		}
		fprintf(debugfile,")-(");
		for(i=0;i<numberofseqs;i++){
			if(i!=0) fprintf(debugfile,",");
			fprintf(debugfile,"%d",((segment->positions[i])+(segment->size)-1));
		}
		fprintf(debugfile,"),s=%d,g=(%d-%d),\"",(segment->size),(segment->mingapsize),(segment->maxgapsize));
		for(i=0;i<(segment->size);i++){
			fprintf(debugfile,"%c",CharAt((segment->positions[0]+i),0));
		}
		fprintf(debugfile,"\"]\n");
		segment=segment->next;
		count++;
	}
	if(onlyactivesegments) fprintf(debugfile,"(%d active alignment map segments)\n\n",activecount);
	else fprintf(debugfile,"(%d alignment map segments)\n\n",count);
	fflush(debugfile);
}

void PrintAlignmentSegmentStrings(alignmapsegment *segment){
	int i;
	for(i=0;i<numberofseqs;i++){
		fprintf(debugfile,"%d:\"",i);
		if((segment->alignedstrings[i])!=NULL) fprintf(debugfile,"%s",segment->alignedstrings[i]);
		fprintf(debugfile,"\"\n");
	}
	fflush(debugfile);
}

// Calculate minimum and maximum gap sizes between this and the next segment among all sequences
void UpdateSegmentGapSizes(alignmapsegment *segment){
	int i,gapsize,min,max,startpos,endpos,textsize;
	min=INT_MAX;
	max=INT_MIN;
	for(i=0;i<numberofseqs;i++){
		textsize=textsizes[i];
		startpos=((segment->positions[i])+(segment->size));
		endpos=(segment->next->positions[i]);
		gapsize=(endpos-startpos);
		if(gapsize<0) gapsize+=textsize;
		if(gapsize<min) min=gapsize;
		if(gapsize>max) max=gapsize;
	}
	segment->mingapsize=min;
	segment->maxgapsize=max;
}

// Converts each of the Heaviest Increasing Subsequence chain items to alignment map segments
// NOTE: the new segments are linked from left to right starting at the endsegment and ending at the startsegment
int SetAlignmentMapSegments(){
	alignmapsegment *currentsegment,*newsegment;
	chainitem *chainitem;
	int i,gapsize,sizesum,min,max,averagemin,averagemax,count;
	int textsize,startpos,endpos;
	currentsegment=endsegment; // currentsegment always points to the next valid segment on the right
	chainitem=chain; // the first element of the H.I.S. is the first element of the chain
	count=0;
	while(chainitem!=NULL){ // create a new alignment map segment for each H.I.S. chain item
		newsegment=NewAlignmentMapSegment(chainitem);
		newsegment->next=currentsegment; // the H.I.S. chain items are sorted from the rightmost position (close to end) to the leftmost position (close to beginning)
		UpdateSegmentGapSizes(newsegment); // update the gap sizes in each sequence between this segment (new) and the next one (current)
		sizesum=0;
		for(i=0;i<numberofseqs;i++){ // calculate the sum of the sizes of the gaps in all sequences
			textsize=textsizes[i];
			startpos=((newsegment->positions[i])+(newsegment->size)); // first position of the gap
			endpos=(currentsegment->positions[i]); // position next to the last position of the gap
			gapsize=(endpos-startpos);
			if(gapsize<0) gapsize+=textsize; // when the end position is rotated before the starting position
			sizesum+=gapsize;
		}
		min=(newsegment->mingapsize);
		max=(newsegment->maxgapsize);
		averagemin=(sizesum-min)/(numberofseqs-1); // average of the gap sizes without the sequence with the shortest gap size
		averagemax=(sizesum-max)/(numberofseqs-1); // average of the gap sizes without the sequence with the largest gap size
		if( (min<(averagemin/2)) || (max>((averagemax*3)/2)) ){ // if the min is less than 1/2 of the average or the max is higher than 1+1/2, discard this segment
			#ifdef DEBUG
			fprintf(debugfile,"DiscardSegment: [p=(");
			for(i=0;i<numberofseqs;i++){
				if(i!=0) fprintf(debugfile,",");
				fprintf(debugfile,"%d",(newsegment->positions[i]));
			}
			fprintf(debugfile,"),s=%d] min=%d (avg=%d) max=%d (avg=%d)\n",(newsegment->size),min,averagemin,max,averagemax);
			#endif
			free(newsegment->positions);
			free(newsegment);
			newsegment=NULL;
			currentsegment=currentsegment; // use the same gap-ending segment as before
		} else {
			#ifdef DEBUG
			fprintf(debugfile,"SetSegment: [p=(");
			for(i=0;i<numberofseqs;i++){
				if(i!=0) fprintf(debugfile,",");
				fprintf(debugfile,"%d",(newsegment->positions[i]));
			}
			fprintf(debugfile,"),s=%d]\n",(newsegment->size));
			#endif
			currentsegment=newsegment; // if it is valid, use this new segment as the next gap-ending segment
			count++;
		}
		chainitem=chainitem->backtrack; // process the previous (to the left) H.I.S. chain item
	}
	startsegment->next=currentsegment; // link the start segment (on the left) to the last valid segment (on the right)
	UpdateSegmentGapSizes(startsegment); // updates the gap sizes in the first segment
	DeleteChain(); // deletes all the no longer needed chain items
	return count;
}


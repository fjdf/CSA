#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "csamsa.h"
#include "alignment.h"
#include "alignmentmap.h"
#include "gencycsuffixtrees.h"
#include "morenodeslinkedlists.h"
#include "dynamicprogramming.h"
#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif


char CharAt(int pos, int seq){
	int i=rotations[seq]+pos;
	if(i>=textsizes[seq]) i-=textsizes[seq];
	return texts[seq][i];
}

// Updates the arrays of starting and ending positions of all the sequences for the gap between the current start and end segments
// NOTE: startpos is the first position that is valid and endpos is the first position that is not valid
void UpdateLimitsPositionArrays(){
	int i;
	for(i=0;i<numberofseqs;i++){
		startposarray[i]=(startsegment->positions[i])+(startsegment->size);
		endposarray[i]=endsegment->positions[i];
	}
}

/*
void FillSegmentGaps(){
	alignmapsegment *segment;
	segment=firstsegment;
	while(segment!=lastsegment){
		ProgressiveDP(segment);
		segment=segment->next;
	}
}
*/

// Initialize first and last border nodes and alignment segments
void InitializeAlignmentVariables(){
	linkedpos *newpositem;
	int i;
	firstbordernode=NewBorderNode(root);
	for(i=0;i<numberofseqs;i++){
		newpositem=(linkedpos *)malloc(sizeof(linkedpos));
		newpositem->k=-1;
		newpositem->next=NULL;
		firstbordernode->positions[i]=newpositem;
	}
	lastbordernode=firstbordernode;
	startposarray=(int *)calloc(numberofseqs,sizeof(int));
	endposarray=(int *)calloc(numberofseqs,sizeof(int));
	firstsegment=NewAlignmentMapSegment(NULL);
	lastsegment=NewAlignmentMapSegment(NULL);
	firstsegment->next=lastsegment;
	firstsegment->size=1;
	for(i=0;i<numberofseqs;i++){
		firstsegment->positions[i]=-1;
		lastsegment->positions[i]=textsizes[i];
	}
	UpdateSegmentGapSizes(firstsegment); // the first segment will have a mingapsize of the size of the smallest sequence and a maxgapsize of the size of the largest sequence
}


// Modify the generalized cyclic suffix tree to mark the deepest nodes that belong to all the sequences (border nodes) and store their rotated positions
void PrepareTreeForAlignment(){
	printf("> Preparing tree for alignment");fflush(stdout);
	InitializeAlignmentVariables();
	printf(".");fflush(stdout);
	ResetBitArrays(root);
	printf(".");fflush(stdout);
	MarkUsedNodes();
	printf(".");fflush(stdout);
	DeleteUnusedNodes(root);
	printf(".");fflush(stdout);
	CollectBorderNodes(root);
	printf(".");fflush(stdout);
	MarkSuffixNodes();
	printf(" ok\n");fflush(stdout);
	#ifdef DEBUG
	PrintBorderNodesList(0);
	#endif
}

// TODO: if the a DP string of a sequence is NULL, check if we should print the sequences or all gaps
// Saves the complete aligment (segment string and DP strings) to file
void SaveAlignment(char *outputfilename){
	alignmapsegment *segment;
	FILE *file;
	char **strings;
	int i,rot,end,size,length,startpos,endpos,gapsize,segmentcount,alignlength;
	segmentcount=0;
	alignlength=0;
	file=fopen(outputfilename,"w");
	if(file==NULL) {printf("> ERROR: Can't write output file\n");return;}
	for(i=0;i<numberofseqs;i++){
		if(rotations!=NULL){ // print description including rotation if it exists
			rot=rotations[i];
			fprintf(file,">%s @ %d\n",descs[i],rot);
		} else {
			rot=0;
			fprintf(file,">%s\n",descs[i]);
		}
		length=textsizes[i];
		end=(length-rot-1); // position in the rotated seqs correspondent to the last character in the original seqs
		segment=firstsegment;
		segmentcount=0;
		alignlength=0;
		while(segment!=lastsegment){ // process all alignment segments
			if(segment!=firstsegment){ // the first segment is a fake one with no valid segment strings
				size=(segment->size); // print the string correspondent to the segment
				alignlength+=size;
				startpos=(segment->positions[i]); // these positions all correspond to rotated positions
				endpos=(startpos+size-1);
				if((startpos<=end) && (endpos<=end)){ // if all characters of the segment string occur before the end of the original seq
					startpos=(startpos+rot);
					fwrite((texts[i]+startpos),sizeof(char),(size_t)size,file);
				}
				else if((startpos>end) && (endpos>end)){ // if all characters of the segment string occur after the end of the original seq
					startpos=(startpos+rot-length);
					fwrite((texts[i]+startpos),sizeof(char),(size_t)size,file);
				} else { // if the characters run over the rotation cutting point
					startpos=(startpos+rot);
					size=(length-startpos);
					fwrite((texts[i]+startpos),sizeof(char),(size_t)size,file);
					size=(startpos+(segment->size)-length);
					fwrite((texts[i]),sizeof(char),(size_t)size,file);
				}
			}
			strings=segment->alignedstrings;
			if(strings!=NULL){ // print the strings correspondent to the DP between this segment and the next one
				if(strings[i]!=NULL){
					fputs(strings[i],file);
					alignlength+=(int)strlen(strings[i]);
				}
				else { // if no strings were set for this gap
					startpos=((segment->positions[i])+(segment->size)); // first position to print
					endpos=(segment->next->positions[i]); // position after the last position to print
					gapsize=(endpos-startpos);
					alignlength+=gapsize;
					startpos+=rot;
					endpos+=rot; // update the positions in the rotated seqs to the positions in the original seqs
					if(startpos>=length) startpos-=length;
					if(endpos>length) endpos-=length; // if it is equal, it always fall on the next first case
					if(startpos<endpos){
						fwrite((texts[i]+startpos),sizeof(char),(size_t)(gapsize),file);
					} else { // if the ending position is before the starting position
						fwrite((texts[i]+startpos),sizeof(char),(size_t)(length-startpos),file); // from startpos to original end (length-1)
						fwrite((texts[i]),sizeof(char),(size_t)(endpos),file); // from original beginning at (0) to (endpos-1)
					}
				}
			}
			segment=segment->next;
			segmentcount++;
		}
		fprintf(file,"\n");
	}
	printf("> Alignment size: %d (%d alignment segments)\n",alignlength,segmentcount);
	fflush(stdout);
	fclose(file);
	DeleteAlignmentMap();
}

// Aligns all the sequences by using iteratively calculating the Heaviest Increasing Subsequences and performing Dynamic Programming
void RunAlignment(){
	int count;
	startsegment=firstsegment; // start with the fake first segment that points to beginning of all sequences
	endsegment=lastsegment; // end with the fake last segment that points to end of all sequences
	/*
	// to align the whole sequences using only DP (without any anchor segments)
	UpdateLimitsPositionArrays();
	ProgressiveDP(startsegment);
	return;
	*/
	while(startsegment!=lastsegment){ // process the segments until the last one at the end of all sequences
		endsegment=startsegment->next; // process the gaps between this segment and the next one
		#ifdef DEBUG
		PrintAlignmentMap(1);
		#endif
		if((startsegment->mingapsize)==0){ // if there are no gaps between this segment and the next one
			startsegment=startsegment->next; // set the next starting segment
			continue;
		}
		UpdateLimitsPositionArrays(); // set the valid starting and ending positions of each sequence for this segment/gap
		count=UpdateActiveBorderNodes(); // deletes no longer needed nodes and hides the ones with no valid positions (for now)
		#ifdef DEBUG
		PrintBorderNodesList(1);
		#endif
		if(count>0){ // if there are active nodes for this gap/segment
			CalculateHeaviestIncreasingSubsequence(); // find the best chain of nodes/positions for the alignment in this gap
			#ifdef DEBUG
			PrintChain();
			#endif
			count=SetAlignmentMapSegments(); // converts the blocks in the chain to new fixed alignment anchors/segments
		}
		if(count==0){ // if no valid nodes/positions or no valid segments exist, perform DP in this segment/gap
			ProgressiveDP(startsegment);
			#ifdef DEBUG
			PrintAlignmentSegmentStrings(startsegment);
			#endif
			startsegment=startsegment->next; // set the next starting segment 
			continue;
		} // if new segments were created, keep the same starting segment
	}
	DeleteBorderNodesList(); // delete no longer needed border nodes
	free(startposarray);
	free(endposarray);
	startposarray=NULL;
	endposarray=NULL;
}

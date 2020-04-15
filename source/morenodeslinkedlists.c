#include <stdio.h>
#include <stdlib.h>
#include "csamsa.h"
#include "alignment.h"
#include "gencycsuffixtrees.h"
#include "morenodeslinkedlists.h"

// Adds a border node structure to the tree node and updates the linked list of border nodes with this new one
bordernode *NewBorderNode(treenode *treenode){
	bordernode *newbordernode;
	newbordernode=(bordernode *)malloc(sizeof(bordernode));
	newbordernode->treenode=treenode;
	newbordernode->size=treenode->depth;
	newbordernode->positions=(linkedpos **)calloc(numberofseqs,sizeof(linkedpos *));
	newbordernode->activeposcount=(int *)calloc(numberofseqs,sizeof(int));
	newbordernode->hidden='f';
	newbordernode->hiddennode=NULL;
	newbordernode->hiddenpos=NULL;
	newbordernode->next=NULL;
	newbordernode->prev=NULL;
	newbordernode->suffix=NULL;
	treenode->bordernode=newbordernode;
	if(lastbordernode!=NULL){
		lastbordernode->next=newbordernode;
		newbordernode->prev=lastbordernode;
		lastbordernode=newbordernode;
	}
	return newbordernode;
}

// Deletes the border node from the tree node and from the linked list of border nodes
void DeleteBorderNode(bordernode *node){
	bordernode *nextnode,*prevnode;
	linkedpos *positem,*nextpositem;
	int i;
	#ifdef DEBUG
	char *nodelabel;
	nodelabel=getNodeLabel(node->treenode);
	fprintf(debugfile,"DeleteNode: [{%d}(%d)\"%s\"]\n",node->treenode->id,node->treenode->depth,nodelabel);
	free(nodelabel);
	#endif
	for(i=0;i<numberofseqs;i++){
		positem=node->positions[i];
		while(positem!=NULL){
			nextpositem=positem->next;
			free(positem);
			positem=nextpositem;
		}
	}
	free(node->positions);
	free(node->activeposcount);
	if(node->hiddenpos!=NULL){
		for(i=0;i<numberofseqs;i++){
			positem=node->hiddenpos[i];
			while(positem!=NULL){
				nextpositem=positem->next;
				free(positem);
				positem=nextpositem;
			}
		}
		free(node->hiddenpos);
	}
	nextnode=node->next;
	prevnode=node->prev;
	if(prevnode!=NULL) prevnode->next=nextnode;
	else firstbordernode=nextnode;
	if(nextnode!=NULL) nextnode->prev=prevnode;
	else lastbordernode=prevnode;
	node->treenode->bordernode=NULL;
	free(node);
}

// Deletes all the nodes from the border nodes list
void DeleteBorderNodesList(){
	bordernode *node,*nextnode;
	node=firstbordernode;
	while(node!=NULL){
		nextnode=node->next;
		DeleteBorderNode(node);
		node=nextnode;
	}
	firstbordernode=NULL;
	lastbordernode=NULL;
}

// Deletes the first position of the sequence in the positions list of the border node
void DeleteFirstPosition(bordernode *node, int seq){
	linkedpos *nextpos;
	#ifdef DEBUG
	char *nodelabel;
	nodelabel=getNodeLabel(node->treenode);
	fprintf(debugfile,"DeletePosition: [{%d}(%d)\"%s\"] %d:(%d)\n",node->treenode->id,node->treenode->depth,nodelabel,seq,node->positions[seq]->k);
	free(nodelabel);
	#endif
	if(node->positions[seq]==NULL) return;
	nextpos=node->positions[seq]->next;
	free(node->positions[seq]);
	node->activeposcount[seq]--;
	node->positions[seq]=nextpos;
}

// Hides the node inside the previous node
// NOTE: if several hidden nodes exist, they are stored in the correct order and hiddennode always points to the higher node
void HideBorderNode(bordernode *nodetohide){
	bordernode *storagenode;
	#ifdef DEBUG
	int i;
	char *nodelabel;
	nodelabel=getNodeLabel(nodetohide->treenode);
	fprintf(debugfile,"HideNode: [{%d}(%d)\"%s\"] (",nodetohide->treenode->id,nodetohide->treenode->depth,nodelabel);
	free(nodelabel);
	for(i=0;i<numberofseqs;i++){
		if(i!=0) fprintf(debugfile,":");
		fprintf(debugfile,"%d",nodetohide->activeposcount[i]);
	}
	fprintf(debugfile,")\n");
	#endif
	if((nodetohide->hidden)=='t') return;
	storagenode=nodetohide->prev;
	storagenode->next=nodetohide->next;
	if((nodetohide->next)!=NULL) nodetohide->next->prev=storagenode;
	nodetohide->next=NULL;
	nodetohide->prev=storagenode->hiddennode;
	if(storagenode->hiddennode!=NULL) storagenode->hiddennode->next=nodetohide;
	storagenode->hiddennode=nodetohide;
	nodetohide->hidden='t';
}

// Restores all hidden nodes inside the given node
void UnHideBorderNodes(bordernode *node){
	bordernode *auxnode;
	auxnode=node->hiddennode;
	if(auxnode==NULL) return;
	auxnode->hidden='f';
	auxnode->next=node->next;
	if(node->next!=NULL) node->next->prev=auxnode;
	while(auxnode->prev!=NULL){
		auxnode=auxnode->prev;
		auxnode->hidden='f';
	}
	auxnode->prev=node;
	node->next=auxnode;
	node->hiddennode=NULL;
}

// Hides the first position of all sequences in the node
// NOTE: the hidden positions list is stored in reverse order (from higher to lower)
// NOTE: only called checking first that it has more than one position for each sequence
void HideFirstPositions(bordernode *node){
	linkedpos *postohide,*nextvisiblepos,*nexthiddenpos;
	int i;
	#ifdef DEBUG
	char *nodelabel;
	nodelabel=getNodeLabel(node->treenode);
	fprintf(debugfile,"HidePositions: [{%d}(%d)\"%s\"] (",node->treenode->id,node->treenode->depth,nodelabel);
	free(nodelabel);
	for(i=0;i<numberofseqs;i++){
		if(i!=0) fprintf(debugfile,",");
		fprintf(debugfile,"%d,",node->positions[i]->k);
	}
	fprintf(debugfile,")\n");
	#endif
	if(node->hiddenpos==NULL) node->hiddenpos=(linkedpos **)calloc(numberofseqs,sizeof(linkedpos *));
	for(i=0;i<numberofseqs;i++){
		postohide=node->positions[i];
		nextvisiblepos=postohide->next;
		nexthiddenpos=node->hiddenpos[i];
		node->positions[i]=nextvisiblepos;
		node->hiddenpos[i]=postohide;
		postohide->next=nexthiddenpos;
		node->activeposcount[i]--;
	}
}

// TODO: simplify the inner loop because the positions are already stored in reverse order (and no new positions were added to the node)
// Restores the hidden positions inside the node
void UnHidePositions(bordernode *node){
	linkedpos *postounhide,*nexthiddenpos,*currentpos,*prevpos;
	int i;
	if(node->hiddennode==NULL) return;
	for(i=0;i<numberofseqs;i++){
		postounhide=node->hiddenpos[i];
		while(postounhide!=NULL){
			nexthiddenpos=postounhide->next;
			prevpos=NULL;
			currentpos=node->positions[i];
			while(currentpos!=NULL && (currentpos->k)<(postounhide->k)){
				prevpos=currentpos;
				currentpos=currentpos->next;
			}
			if(prevpos!=NULL) prevpos->next=postounhide;
			else node->positions[i]=postounhide;
			postounhide->next=currentpos;
			node->activeposcount[i]++;
			postounhide=nexthiddenpos;
		}
		node->hiddenpos[i]=NULL;
	}
}

// TODO: add comments and descriptions to code
void HideSuffixPositions(bordernode *node, bordernode *suffixnode){
	linkedpos *positem,*suffixpositem,*prevpositem;
	int i,j,k,n,diff,posval;
	diff=(node->size)-(suffixnode->size);
	n=node->activeposcount[0];
	if(n>(suffixnode->activeposcount[0])) return;
	if(suffixnode->hiddenpos==NULL) suffixnode->hiddenpos=(linkedpos **)calloc(numberofseqs,sizeof(linkedpos *));
	for(j=0;j<n;j++){
		for(i=0;i<numberofseqs;i++){
			k=1;
			positem=node->positions[i];
			while(k<j){
				positem=positem->next;
				k++;
			}
			posval=(positem->k)+diff;
			if(posval<endposarray[i]){
				prevpositem=NULL;
				suffixpositem=suffixnode->positions[i];
				while(suffixpositem!=NULL && (suffixpositem->k)!=posval){
					prevpositem=suffixpositem;
					suffixpositem=suffixpositem->next;
				}
				if(suffixpositem!=NULL){
					if(prevpositem!=NULL) prevpositem->next=suffixpositem->next;
					else suffixnode->positions[i]=suffixpositem->next;
					suffixpositem->next=suffixnode->hiddenpos[i];
					suffixnode->hiddenpos[i]=suffixpositem;
					suffixnode->activeposcount[i]--;
				}
			}
		}
	}
}

// TODO: add comments and descriptions to code
// TODO: find if function is usefull somewhere or not (currently not called anywhere)
void HideSuffixNodes(bordernode *node){
	bordernode *suffixnode;
	int i;
	suffixnode=node->suffix;
	while(suffixnode!=NULL){
		if((suffixnode->hidden)=='t') return;
		if((suffixnode->activeposcount[0])==(node->activeposcount[0])) HideBorderNode(suffixnode);
		else{
			HideSuffixPositions(node,suffixnode);
			for(i=0;i<numberofseqs;i++){
				if(suffixnode->activeposcount[i]<1){
					HideBorderNode(suffixnode);
					break;
				}
			}
		}
		suffixnode=suffixnode->suffix;
	}
}

// Move positions from linked list in tree node to linked list in border node in a sorted way (lowest position to highest position)
void AddPositions(treenode *tnode, bordernode *bnode){
	linkedpos **fromarray,**toarray;
	linkedpos *sourceitem,*destitem,*nextsourceitem,*prevdestitem;
	int i;
	fromarray=tnode->positions;
	toarray=bnode->positions;
	for(i=0;i<numberofseqs;i++){
		sourceitem=fromarray[i];
		while(sourceitem!=NULL){
			nextsourceitem=sourceitem->next;
			destitem=toarray[i];
			if(destitem==NULL || ((destitem->k)>(sourceitem->k))){ // if bnode positions do not exist or are larger than tnode position
				sourceitem->next=destitem;
				toarray[i]=sourceitem;
				sourceitem=nextsourceitem;
				continue;
			}
			prevdestitem=destitem;
			destitem=destitem->next;
			while(destitem!=NULL && ((destitem->k)<(sourceitem->k))){
				prevdestitem=destitem;
				destitem=destitem->next;
			}
			prevdestitem->next=sourceitem;
			sourceitem->next=destitem;
			sourceitem=nextsourceitem;
		}
	}
	free(tnode->positions);
	tnode->positions=NULL;
}

// Collect in the border node all the positions from the tree node and from all branches bellow
void CollectBorderPositions(treenode *tnode, bordernode *bnode){
	int i,n;
	if((tnode->positions)!=NULL) AddPositions(tnode,bnode);
	n=tnode->numbranches;
	for(i=0;i<n;i++) CollectBorderPositions(tnode->branches[i],bnode);
}


// TODO: prevent positions to be added to the root node
// Collect positions for all border nodes (nodes that belong to all sequences but at least one of its branches does not)
void CollectBorderNodes(treenode *node){
	bordernode *bordernode;
	treenode *branch;
	int i, n;
	if(!nodeFromAllSeqs(node)) return;
	bordernode=node->bordernode;
	if((node->positions)!=NULL){ // if the tree node has positions, collect them too, and not only the positions from the branches
		if(bordernode==NULL) bordernode=NewBorderNode(node);
		AddPositions(node,bordernode);
	}
	n=(node->numbranches);
	for(i=0;i<n;i++){
		branch=node->branches[i];
		if(!nodeFromAllSeqs(branch)){
			if(bordernode==NULL) bordernode=NewBorderNode(node); // the tree node is updated with the new border node inside the function
			CollectBorderPositions(branch,bordernode);
		}
		else CollectBorderNodes(branch);
	}
	if(bordernode!=NULL){
		for(i=0;i<numberofseqs;i++)	if(bordernode->positions[i]==NULL) break;
		if(i!=numberofseqs) DeleteBorderNode(bordernode);  // delete the border node if not all sequences have positions in the positions list
	}
}

// Print the information of the node corresponding to the given label
void PrintTreeNodeInfoByLabel(char *nodelabel){
	treenode *currentnode, *branch;
	linkedpos *pos;
	int k, n, i, m, ls, textsize;
	char *text, *fullnodelabel;
	currentnode=followText(nodelabel);
	if(currentnode==NULL){
		fprintf(debugfile,"\nNode with label \"%s\" not found!\n\n",nodelabel);
		return;
	}
	fullnodelabel=getNodeLabel(currentnode);
	ls=((currentnode->endpos)-(currentnode->startpos)+1);
	fprintf(debugfile,"\n[nid=%d,\"%s\",d=%d,l=%d,s=",(currentnode->id),fullnodelabel,(currentnode->depth),ls);
	for(k=0;k<numberofseqs;k++) fprintf(debugfile,"%c",(nodeFromSeq(currentnode,k)?'1':'0'));
	fprintf(debugfile,",p=(");
	if(currentnode->positions!=NULL){ // print all positions inside node
		for(k=0;k<numberofseqs;k++){
			if(k!=0) fprintf(debugfile,";");
			pos=(currentnode->positions[k]);
			while(pos!=NULL){
				fprintf(debugfile,"%d",(pos->k));
				pos=(pos->next);
				if(pos!=NULL) fprintf(debugfile,",");
			}
		}
	}
	fprintf(debugfile,")]");
	if(currentnode->isleaf) fprintf(debugfile,"*");
	fprintf(debugfile,"\n");
	n=(currentnode->numbranches);
	for(k=0;k<n;k++){ // print info for all branches
		branch=(currentnode->branches[k]);
		text=texts[(branch->labelfrom)];
		textsize=textsizes[(branch->labelfrom)];
		ls=((branch->endpos)-(branch->startpos)+1);
		fprintf(debugfile,"\t|->[%d:\"...",(branch->id));
		m=(branch->endpos);
		for(i=(branch->startpos);i<=m;i++) fprintf(debugfile,"%c",((i<textsize)?text[i]:text[(i-textsize)]));
		fprintf(debugfile,"\":");
		for(m=0;m<numberofseqs;m++) fprintf(debugfile,"%c",(nodeFromSeq(branch,m)?'1':'0'));
		fprintf(debugfile,":(");
		if(branch->positions!=NULL){ // print all positions inside branch
			for(i=0;i<numberofseqs;i++){
				if(i!=0) fprintf(debugfile,";");
				pos=(branch->positions[i]);
				while(pos!=NULL){
					fprintf(debugfile,"%d",(pos->k));
					pos=(pos->next);
					if(pos!=NULL) fprintf(debugfile,",");
				}
			}
		}
		fprintf(debugfile,")]");
		if(branch->isleaf) fprintf(debugfile,"*");
		fprintf(debugfile,"\n");
	}
	fprintf(debugfile,"\n");
	fflush(debugfile);
	free(fullnodelabel);
}

// TODO: check if it is ok to have suffix links to border nodes that are more than one suffix link jump away
// Links border nodes to their correspondent suffix border nodes
void MarkSuffixNodes(){
	treenode *treenode;
	bordernode *bordernode,*suffixbordernode;
	bordernode=firstbordernode->next;
	while(bordernode!=NULL){
		treenode=bordernode->treenode->suffixlink;
		while(treenode!=root){
			suffixbordernode=treenode->bordernode;
			if(suffixbordernode!=NULL){
				bordernode->suffix=suffixbordernode;
				break;
			}
			treenode=treenode->suffixlink;
		}
		bordernode=bordernode->next;
	}
}

// Sort all the valid border nodes (before the current ending position) by their lower position in the first sequence
void SortBorderNodes(){
	bordernode *currentnode,*prevnode,*nextnode,*backsearchnode,*forwardsearchnode,*followingnode;
	#ifdef DEBUGSORT
	int i,n=0;
	currentnode=firstbordernode->next;
	while(currentnode!=NULL && (currentnode->positions[0]->k)<endposarray[0]){
		currentnode=currentnode->next;
		n++;
	}
	i=0;
	printf("\n");
	#endif
	currentnode=firstbordernode->next;
	while(currentnode!=NULL && (currentnode->positions[0]->k)<endposarray[0]){
		#ifdef DEBUGSORT
		i++;
		printf("\033[30Dnode %d of %d",i,n);
		fflush(stdout);
		#endif
		prevnode=currentnode->prev;
		if((currentnode->positions[0]->k)<(prevnode->positions[0]->k)){ // if the current node is lower than the previous one
			backsearchnode=currentnode->prev;
			while(backsearchnode!=NULL && (backsearchnode->positions[0]->k)>(currentnode->positions[0]->k))
				backsearchnode=backsearchnode->prev; // find the first node backwards that is higher than the current one
			followingnode=backsearchnode->next;
			backsearchnode->next=currentnode;
			currentnode->prev=backsearchnode;
			forwardsearchnode=currentnode;
			while((forwardsearchnode->next)!=NULL // find the furthest node ahead of the current node that is still lower than the node that is next to the position where the current node will be placed
				&& (forwardsearchnode->next->positions[0]->k)>(forwardsearchnode->positions[0]->k)
				&& (forwardsearchnode->next->positions[0]->k)<(followingnode->positions[0]->k))
				forwardsearchnode=forwardsearchnode->next;
			nextnode=forwardsearchnode->next;
			forwardsearchnode->next=followingnode;
			followingnode->prev=forwardsearchnode;
			prevnode->next=nextnode;
			if(nextnode!=NULL) nextnode->prev=prevnode;
		} else { // if the current node is in the correct position
			nextnode=currentnode->next;
		}
		currentnode=nextnode; // go to the next node
	}
}

// Finds the correct position of the node in the sorted border nodes list after its first positions have been hidden
void ReSortBorderNode(bordernode *node){
	bordernode *currentnode,*prevnode,*nextnode;
	if(node->next==NULL || (node->next->positions[0]->k)>(node->positions[0]->k)) return; // the new positions are always higher than the old hidden ones
	currentnode=(node->next);
	while((currentnode->next)!=NULL && (currentnode->next->positions[0]->k)<(node->positions[0]->k))
		currentnode=currentnode->next; // find the node that will be behind this one
	prevnode=node->prev;
	nextnode=node->next;
	if(prevnode!=NULL) prevnode->next=nextnode;
	if(nextnode!=NULL) nextnode->prev=prevnode;
	nextnode=currentnode->next;
	currentnode->next=node;
	node->prev=currentnode;
	if(nextnode!=NULL) nextnode->prev=node;
	node->next=nextnode;
}

// Deletes all the positions/nodes that are behind the current starting position and hides the nodes that do not have valid positions
int UpdateActiveBorderNodes(){
	bordernode *bordernode,*nextbordernode;
	linkedpos *positem,*nextpositem;
	int i,count,activenodescount;
	bordernode=firstbordernode->next; // the first node has always to be present to store hidden nodes (e.g. when hidding the 2nd node)
	while(bordernode!=NULL && (bordernode->positions[0]->k)<endposarray[0]){ // delete positions behind current start position
		if(bordernode->hiddennode!=NULL) UnHideBorderNodes(bordernode); // restore all previously hidden nodes
		if(bordernode->hiddenpos!=NULL) UnHidePositions(bordernode); // restores all previously hidden positions (hidden by CalculateHeaviestIncreasingSubsequence if the node had more than one active position for each sequence)
		nextbordernode=bordernode->next;
		for(i=0;i<numberofseqs;i++){
			positem=bordernode->positions[i]; // the first time the function is called the start positions are all -1 and the end positions are the sizes of the sequences
			while(positem!=NULL && (positem->k)<startposarray[i]){ // for each sequence delete the positions that are already behind the current start position
				nextpositem=positem->next;
				DeleteFirstPosition(bordernode,i);
				positem=nextpositem;
			}
			if((bordernode->positions[i])==NULL){ // delete the border node if one of the sequences runs out of positions
				DeleteBorderNode(bordernode);
				break;
			}
		}
		bordernode=nextbordernode;
	}
	SortBorderNodes(); // sort border nodes by their (first) position in the first sequence
	activenodescount=0;
	bordernode=firstbordernode->next;
	while(bordernode!=NULL && (bordernode->positions[0]->k)<endposarray[0]){
		activenodescount++;
		for(i=0;i<numberofseqs;i++){ // count the number of valid positions for each sequence
			count=0;
			positem=bordernode->positions[i];
			while(positem!=NULL && (positem->k)<endposarray[i]){
				count++;
				positem=positem->next;
			}
			if(count==0) break;
			bordernode->activeposcount[i]=count;
		}
		nextbordernode=bordernode->next;
		if(i!=numberofseqs){ // if at least one sequence does not have any valid positions (for now), hide the node
			HideBorderNode(bordernode);
			activenodescount--;
			bordernode=nextbordernode;
			continue;
		}
		count=bordernode->activeposcount[0];
		for(i=1;i<numberofseqs;i++){ // if the node does not have the same number of positions on all sequences, hide it (because when calculating the HIS, all the first positions of all sequences are hidden at the same time)
			if((bordernode->activeposcount[i])!=count){
				HideBorderNode(bordernode);
				activenodescount--;
				break;
			}
		}
		bordernode=nextbordernode;
	}
	return activenodescount;
}

// Clear the information on this node and all nodes bellow, of which sequences the node belongs to (for latter marking only the deepest nodes that correspond to the rotation and all its suffixes on all sequences)
void ResetBitArrays(treenode *node){
	int i,n;
	node->fromseqs=0UL;
	n=(node->numbranches);
	for(i=0;i<n;i++) ResetBitArrays(node->branches[i]);
}

// Update the linked list of positions in this tree node with a new position
void InsertPosition(treenode *node, int seq, int pos){
	linkedpos *newpos;
	newpos=(linkedpos *)malloc(sizeof(linkedpos));
	newpos->k=pos;
	newpos->next=NULL;
	if(node->positions==NULL){
		node->positions=(linkedpos **)calloc((size_t)numberofseqs,sizeof(linkedpos *));
		node->positions[seq]=newpos;
		return;
	}
	newpos->next=node->positions[seq];
	node->positions[seq]=newpos;
}

// Adds the positions of the rotation and all its suffixes of all the sequences to the corresponding nodes in the tree (as if the tree was directly built with the linear rotated sequences)
void MarkUsedNodes(){
	treenode *deepestnode,*prevnode,*currentnode,*backnode,*splitnode,*prevsplitnode;
	int i,j,k,pos,labelsize;
	for(j=0;j<numberofseqs;j++){
		text=texts[j];
		textsize=textsizes[j];
		currentseq=j; // needed to be set for splitNode
		markNodeFromSeq(root,j);
		currentnode=root;
		pos=1;
		i=rotations[j];
		k=0;
		while(k<textsize){ // find deepest node starting at the rotation and mark all tree nodes belonging to this sequence
			currentnode=followChar(currentnode,&pos,i); // get correspondent branch bellow
			markNodeFromSeq(currentnode,j);
			pos=getLabelSize(currentnode); // go to end of node
			i+=pos; // advance as many chars as the size of the node label
			k+=pos;
		} // in the last node (leaf), pos already points to the end of the node (label size) where the rotation ends
		InsertPosition(currentnode,j,0); // the deepest node starts at position 0 of the (rotated) sequence
		k=1; // the suffix of the node will start at position 1
		deepestnode=currentnode;
		prevnode=deepestnode;
		prevsplitnode=NULL;
		labelsize=getLabelSize(currentnode);
		currentnode=getSuffixNode(currentnode,labelsize,&pos); // the deepest node is already explicit
		while(currentnode!=root){ // make all nodes corresponding to the final position and all its suffixes explicit, until the root
			labelsize=getLabelSize(currentnode);
			if(pos!=labelsize){ // split if needed
				splitnode=splitNode(currentnode,pos); // new split node is created between backnode above and currentnode below
				if(prevsplitnode!=NULL) prevsplitnode->suffixlink=splitnode;
				prevsplitnode=splitnode;
			} else { // if the node is already explicit
				splitnode=currentnode;
				if(prevsplitnode!=NULL){
					prevsplitnode->suffixlink=currentnode;
					prevsplitnode=NULL;
				}
			}
			backnode=splitnode;
			while(!nodeFromSeq(backnode,j)){ // mark all nodes above as belonging to this sequence
				markNodeFromSeq(backnode,j);
				backnode=backnode->backlink;
			}
			if(prevnode==deepestnode){ // if this is the second deepest node, link it to the deepest (not split because it is explicit)
				deepestnode->suffixlink=splitnode;
				deepestnode->isleaf=0; // to follow the suffix link to the correct node
			}
			InsertPosition(splitnode,j,k);
			k++; // advance position in (rotated) sequence
			prevnode=currentnode;
			if(splitnode->isleaf) pos++; // if this node is a leaf/rotation from another sequence, fix the number of chars to go up when following the suffix link
			currentnode=getSuffixNode(currentnode,labelsize,&pos); // go to next suffix
		}
	}
}

// Deletes all unnecessary tree nodes that do not correspond to a suffix of a rotation of any sequence
void DeleteUnusedNodes(treenode *node){
	int i,j,n,m;
	n=(node->numbranches);
	for(i=0;i<n;i++){
		if((node->branches[i]->fromseqs)==0UL){
			freeTreeNode(node->branches[i]);
			(node->branches[i])=NULL;
			(node->numbranches)--;
		}
		else DeleteUnusedNodes(node->branches[i]);
	}
	m=(node->numbranches);
	if(m==n) return;
	for(i=0;i<m;i++){ // join all remaining branches in the beginning of the array so no NULL entries appear in the middle
		if((node->branches[i])==NULL){
			j=i+1;
			while((node->branches[j])==NULL) j++;
			(node->branches[i])=(node->branches[j]);
			(node->branches[j])=NULL;
		}
	}
}

/*
void TestSort(){
	bordernode *bordernode;
	linkedpos *positem;
	int i,count,prevnodeval,nodeval,prevposval,posval,notsortednodes,notsortedpos;
	notsortednodes=0;
	notsortedpos=0;
	count=0;
	prevnodeval=-1;
	bordernode=firstbordernode->next;
	while(bordernode!=NULL){
		count++;
		nodeval=bordernode->positions[0]->k;
		if(nodeval<prevnodeval) notsortednodes++;
		prevnodeval=nodeval;
		for(i=0;i<numberofseqs;i++){
			positem=bordernode->positions[i];
			prevposval=-1;
			while(positem!=NULL){
				posval=positem->k;
				if(posval<prevposval) notsortedpos++;
				prevposval=posval;
				positem=positem->next;
			}
		}
		bordernode=bordernode->next;
	}
	printf("[%d nodes: (n:%d,p:%d) not sorted]\n",count,notsortednodes,notsortedpos);
}
*/

void PrintBorderNodesList(int onlyactivenodes){
	treenode *tnode;
	bordernode *bnode;
	linkedpos *pos;
	char *nodelabel;
	int i,k,count,activecount;
	fprintf(debugfile,"\n");
	if(onlyactivenodes){
		fprintf(debugfile,"start=(");
		for(i=0;i<numberofseqs;i++){
			if(i!=0) fprintf(debugfile,",");
			fprintf(debugfile,"%5d",startposarray[i]);
		}
		fprintf(debugfile,")\nend  =(");
		for(i=0;i<numberofseqs;i++){
			if(i!=0) fprintf(debugfile,",");
			fprintf(debugfile,"%5d",endposarray[i]);
		}
		fprintf(debugfile,")\n");
	}
	bnode=firstbordernode;
	count=0;
	activecount=0;
	while(bnode!=NULL){
		if(onlyactivenodes){
			for(i=0;i<numberofseqs;i++){
				pos=(bnode->positions[i]);
				if( (pos->k)<startposarray[i] || (pos->k)>=endposarray[i] ) break;
			}
			if(i==numberofseqs){ // it is an active node
				fprintf(debugfile,"*");
				activecount++;
			}
		}
		tnode=bnode->treenode;
		nodelabel=getNodeLabel(tnode);
		fprintf(debugfile,"->[{%d}(%d)\"%s\"]",tnode->id,tnode->depth,nodelabel);
		free(nodelabel);
		if(bnode->positions!=NULL){
			for(i=0;i<numberofseqs;i++){
				k=0;
				pos=(bnode->positions[i]);
				while(pos!=NULL){
					pos=(pos->next);
					k++;
				}
				fprintf(debugfile," (%d:",k);
				pos=(bnode->positions[i]);
				while(pos!=NULL){
					fprintf(debugfile,"%d",pos->k);
					pos=(pos->next);
					if(pos!=NULL) fprintf(debugfile,",");
				}
				fprintf(debugfile,")");
			}
			fprintf(debugfile,"\n");
		}
		bnode=bnode->next;
		count++;
	}
	if(onlyactivenodes) fprintf(debugfile,"(%d active border nodes)\n\n",activecount);
	else fprintf(debugfile,"(%d border nodes)\n\n",count);
	fflush(debugfile);
}

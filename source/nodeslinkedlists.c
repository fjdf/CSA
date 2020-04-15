#include <stdio.h>
#include <stdlib.h>
#include "csamsa.h"
#include "nodeslinkedlists.h"
#include "gencycsuffixtrees.h"
#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif

linkedblock *createItem(treenode *item){
	linkedblock *newitem;
	newitem=(linkedblock *)malloc(sizeof(linkedblock));
	newitem->item=item;
	newitem->size=0;
	newitem->totalsize=0;
	newitem->interval=0;
	newitem->positions=NULL;
	newitem->nextblock=NULL;
	newitem->prev=NULL;
	newitem->next=NULL;
	return newitem;
}

// adds a new block to the beggining of a list of blocks (and returns a pointer to that new block)
linkedblock *addItem(linkedblock *list, treenode *item){
	linkedblock *newitem;
	newitem=createItem(item);
	newitem->next=list;
	if(list!=NULL) list->prev=newitem;
	return newitem;
}

// inserts item in list by decreasing depth
linkedblock *insertSortedItem(linkedblock *list, treenode *item){
	linkedblock *newitem,*current;
	newitem=createItem(item);
	if(list==NULL) return newitem;
	if((item->depth)>=(list->item->depth)){
		newitem->next=list;
		list->prev=newitem;
		return newitem;
	}
	current=list;
	while((current->next)!=NULL && (item->depth)<(current->next->item->depth))
		current=current->next;
	newitem->next=current->next;
	newitem->prev=current;
	if((current->next)!=NULL) current->next->prev=newitem;
	current->next=newitem;
	return list;
}

// TODO: faster algorithm
// Sort the linked list of blocks by their size
linkedblock *sortList(linkedblock *list){
	linkedblock *block,*tmpblock,*maxblock;
	block=list;
	while(block!=NULL){
		maxblock=block;
		tmpblock=block->next;
		while(tmpblock!=NULL){ // find the largest block/chain from here to the end of the list
			if((tmpblock->size)>(maxblock->size)) maxblock=tmpblock;
			tmpblock=tmpblock->next;
		}
		if(maxblock!=block){ // and insert it before the current position
			if((maxblock->prev)!=NULL) (maxblock->prev)->next=maxblock->next;
			if((maxblock->next)!=NULL) (maxblock->next)->prev=maxblock->prev;
			if((block->prev)!=NULL)	(block->prev)->next=maxblock;
			maxblock->prev=block->prev;
			block->prev=maxblock;
			maxblock->next=block;
		} else block=block->next; // the maximum is now before the current position, so lets do the search again to check if the same current block is now the maximum
	}
	block=list; // if items were moved to the beggining of the list, then find the new list starting point
	while((block->prev)!=NULL) block=block->prev;
	return block;
}

// deletes the item in the current position and returns the item in the previous position
// (if it does not exist, the item in the following position is returned instead)
linkedblock *deleteItem(linkedblock *list){
	linkedblock* updatedlist;
	updatedlist=NULL;
	if((list->next)!=NULL) {
		list->next->prev=list->prev;
		updatedlist=list->next;
	}
	if((list->prev)!=NULL) {
		list->prev->next=list->next;
		updatedlist=list->prev;
	}
	if((list->positions)!=NULL) free(list->positions);
	free(list);
	return updatedlist;
}

linkedblock *deleteItemAndAdvance(linkedblock *list){
	linkedblock *nextitem;
	if(list==NULL) return NULL;
	nextitem=(list->next);
	deleteItem(list); // if deletion occurs on the first position, the next item is already returned
	return nextitem;
}

void clearList(linkedblock *list){
	linkedblock *current,*temp;
	if(list==NULL) return;
	current=list->prev;
	while(current!=NULL){
		temp=current->prev;
		if((current->positions)!=NULL) free(current->positions);
		free(current);
		current=temp;
	}
	current=list->next;
	while(current!=NULL){
		temp=current->next;
		if((current->positions)!=NULL) free(current->positions);
		free(current);
		current=temp;
	}
	if((list->positions)!=NULL) free(list->positions);
	free(list);
}

// TODO: check if the interval can be lower than zero
// Outputs a string containing all the substrings of all the blocks in this chain, separated by as many gaps as the corresponding intervals
char *blockLabel(linkedblock *block){
	int i,n,labelsize,seq,pos,labelpos,nodetextsize;
	char *nodetext,*label,*numgapsstring;
	linkedblock *tmpblock;
	treenode *node,*current;
	numgapsstring=(char *)calloc(15,sizeof(char));
	labelsize=(block->totalsize);
	label=(char *)calloc((labelsize+1),sizeof(char));
	labelpos=0;
	/*
	tmpblock=block;
	while(tmpblock!=NULL){
		node=tmpblock->item;
		n=(node->endpos)-(node->depth)+1;
		printf("->(id:%d;depth:%d;from:%d@[%d-%d];int:%d;rot:%d;pos:%d)",
			node->id,node->depth,node->labelfrom,node->startpos,node->endpos,tmpblock->interval,node->rotation,n);
		tmpblock=tmpblock->nextblock;
	}
	printf("\n");
	*/
	tmpblock=block;
	while(tmpblock!=NULL){
		node=tmpblock->item;
		n=node->depth;
		current=node;
		labelpos+=(n-1);
		while(current!=root){
			seq=current->labelfrom;
			nodetext=texts[seq];
			nodetextsize=textsizes[seq];
			for(i=(current->endpos);i>=(current->startpos);i--){
				pos=i;
				if(pos>=nodetextsize) pos-=nodetextsize;
				label[labelpos]=nodetext[pos];
				labelpos--;
			}
			current=current->backlink;
		}
		labelpos+=(n+1);
		n=tmpblock->interval;
		if(n<0){
			labelpos+=n; // ???
		} else {
			if(n>7){
				sprintf(numgapsstring,"-(%d)-",n);
				n=0;
				while(numgapsstring[n]!='\0') n++;
				for(i=0;i<n;i++){
					label[labelpos]=numgapsstring[i];
					labelpos++;
				}
			} else {
				for(i=0;i<n;i++){
					label[labelpos]='-';
					labelpos++;
				}
			}
		}
		tmpblock=tmpblock->nextblock;
	}
	free(numgapsstring);
	label[labelpos]='\0';
	return label;
}

void printBlocksList(linkedblock *list){
	linkedblock *current;
	char *nodelabel;
	int i, j, n, *numpos, **pos;
	if(list==NULL) return;
	numpos=(int *)calloc(numberofseqs,sizeof(int));
	current=list;
	while((current->prev)!=NULL) current=current->prev;
	n=0;
	while(current!=NULL){
		nodelabel=getNodeLabel(current->item);
		fprintf(debugfile,"->[{%d}(%d)\"%s\"]",current->item->id,current->item->depth,nodelabel);
		free(nodelabel);
		pos=getNodePositions(current->item,numpos);
		for(i=0;i<numberofseqs;i++){
			if(numpos[i]!=1) break;
		}
		if(i==numberofseqs) fprintf(debugfile,"*");
		for(i=0;i<numberofseqs;i++){
			fprintf(debugfile,"(");
			for(j=0;j<numpos[i];j++){
				fprintf(debugfile,"%d",pos[i][j]);
				if(j!=(numpos[i]-1)) fprintf(debugfile,",");
			}
			fprintf(debugfile,")");
		}
		for(i=0;i<numberofseqs;i++) free(pos[i]);
		free(pos);
		fprintf(debugfile,"\n");
		current=current->next;
		n++;
	}
	fprintf(debugfile,"(%d nodes)\n\n",n);
	fflush(debugfile);
	free(numpos);
}

void printChainsList(linkedblock *list){
	linkedblock *current, *chainnode;
	char *nodelabel;
	int i,n,numblocks,sumsizes;
	if(list==NULL) return;
	current=list;
	while((current->prev)!=NULL) current=current->prev;
	n=0;
	while(current!=NULL){
		sumsizes=0;
		numblocks=0;
		chainnode=current;
		while(chainnode!=NULL){
			numblocks++;
			sumsizes+=chainnode->item->depth;
			chainnode=chainnode->nextblock;
		}
		nodelabel=getNodeLabel(current->item);
		if(current->totalsize>0){
			fprintf(debugfile,"=>[%d](%d)",numblocks,sumsizes);
			n++;
		}
		fprintf(debugfile,"->[{%d}(%d)\"%s\"(%d;%d)]",current->item->id,current->item->depth,nodelabel,current->size,current->totalsize);
		free(nodelabel);
		if(current->positions!=NULL){
			fprintf(debugfile,"(");
			for(i=0;i<numberofseqs;i++){
				fprintf(debugfile,"%d",current->positions[i]);
				if(i!=(numberofseqs-1)) fprintf(debugfile,",");
			}
			fprintf(debugfile,")\n");
		}
		chainnode=current->nextblock;
		while(chainnode!=NULL){
			nodelabel=getNodeLabel(chainnode->item);
			fprintf(debugfile,"\t->[{%d}(%d)\"%s\"(%d;%d)]",chainnode->item->id,chainnode->item->depth,nodelabel,chainnode->size,chainnode->totalsize);
			free(nodelabel);
			if(chainnode->positions!=NULL){
				fprintf(debugfile,"(");
				for(i=0;i<numberofseqs;i++){
					fprintf(debugfile,"%d",chainnode->positions[i]);
					if(i!=(numberofseqs-1)) fprintf(debugfile,",");
				}
				fprintf(debugfile,")\n");
			}
			chainnode=chainnode->nextblock;
		}
		current=current->next;
	}
	fprintf(debugfile,"(%d node chains)\n\n",n);
	fflush(debugfile);
}

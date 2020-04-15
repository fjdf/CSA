#include "csamsa.h"
#include "gencycsuffixtrees.h"
#include "morenodeslinkedlists.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif

#define ALPHABETSIZE 5

int infinity; // last position of the label of the "open leaves"
int numberofnodes; // total number of nodes in the tree
int alphabetsize; // number of characters of the alphabet
unsigned long long int *masks; // list of the marks of each one of the sequences
unsigned long long int allseqsmask; // mark of the nodes that belong to all sequences at once


// marca um no' como pertencente a uma sequencia
void markNodeFromSeq(treenode *node, int i){
	node->fromseqs = (node->fromseqs) | masks[i];
}

// verifica se o no' pertence a uma sequencia
int nodeFromSeq(treenode *node, int i){
	//return ((node->fromseqs) & masks[i]);
	if( (node->fromseqs) & masks[i] ) return 1;
	return 0;
}

// verifica se o no' pertence a todas as sequencias
int nodeFromAllSeqs(treenode *node){
	//return ( (node->fromseqs) == allseqsmask );
	if( (~(node->fromseqs)) & allseqsmask ) return 0;
	return 1;
}

// devolve o numero de sequencias 'as quais o no' pertence
int getNumberOfSeqs(treenode *node){
	int i,k=0;
	for(i=0;i<numberofseqs;i++)	if( (node->fromseqs) & masks[i] ) k++;
	return k;
}

// inicializa a lista das marcas para cada sequencia e a marca para todas as sequencias
void initializeMasks(){
	int i;
	unsigned long long int mask = 1u;
	allseqsmask = 0u;
	masks=(unsigned long long int *)calloc(numberofseqs,sizeof(unsigned long long int));
	for(i=0;i<numberofseqs;i++){
		masks[i]=mask;
		allseqsmask = ( allseqsmask | mask );
		mask = mask << 1;
	}
}

// mostra os bits do identificador de cada sequencia
void printMasks(){
	int i,j,sizeofmask,bit;
	unsigned long long int bitmask;
	sizeofmask=8*(int)sizeof(allseqsmask);
	for(i=0;i<numberofseqs;i++){
		printf("%02d (%dbits) ",i,sizeofmask);
		bitmask = 1u;
		bitmask = bitmask << (sizeofmask-1);
		for(j=0;j<sizeofmask;j++){
			bit=0;
			if(( masks[i] & bitmask )) bit=1;
			printf("%d",bit);
			bitmask = bitmask >> 1;
		}
		printf("\n");
	}
	fflush(stdout);
}

// Checks if the given node is a leaf node
int isLeaf(treenode *node){
	return ((node->numbranches)==0);
}

// Returns the size of the node label
int getLabelSize(treenode *node){
	return ((node->endpos)-(node->startpos)+1);
}

// Returns the position in the sequence where the suffix represented by this node begins
int getPosition(treenode *node){
	return ((node->endpos)-(node->depth)+1);
}

// Returns the label of the node (from the root to the end of the node)
char *getNodeLabel(treenode *node){
	int i,n,seq,pos,labelpos,nodetextsize;
	char *nodetext,*label;
	treenode *current;
	n=node->depth;
	label=(char *)calloc((n+1),sizeof(char));
	current=node;
	labelpos=(n-1);
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
	return label;
}

// TODO: get positions from sequences were one rotation is a prefix of a rotation of another sequences (leaf node followed by deeper nodes if we continue down in the tree)
// Returns the positions in each sequence were this substring associated with this node occurs
int **getNodePositions(treenode *node, int *numpositions){
	treenode *currentnode, *backnode;
	int i,n,**positions;
	if(numpositions==NULL) numpositions=(int *)malloc(numberofseqs*sizeof(int));
	for(i=0;i<numberofseqs;i++) numpositions[i]=0;
	currentnode=node;
	while(1){
		while(!isLeaf(currentnode)){ // starting from this node, find the deepest leaf node
			/*
			i=0;
			while((currentnode->branches[i])==NULL) i++; // follow the first/leftmost branch
			currentnode=(currentnode->branches[i]);
			*/
			currentnode=(currentnode->branches[0]);
		}
		for(i=0;i<numberofseqs;i++){ // increment number of positions
			if(nodeFromSeq(currentnode,i)) numpositions[i]++;
		}
		backnode=currentnode;
		/*
		i=alphabetsize;
		*/
		n=i; // just initializing variables to enter loop
		while(i==n && backnode!=node){ // find the first branch at the right or above that has not been processed
			currentnode=backnode;
			backnode=(currentnode->backlink);
			n=(backnode->numbranches);
			for(i=0;i<n;i++){ // find which branch of the upper node is the current node
				if((backnode->branches[i])==currentnode){
					i++; // proceed to the next branch (if it is the last one, we will go up)
					break;
				}
			}
			/*
			while(i<alphabetsize && (backnode->branches[i])==NULL) i++; // skip non-existent branches
			*/
		}
		if(i==n) break; // if we have returned to the initial node, we are done
		currentnode=(backnode->branches[i]); // continue down by checking this branch
	}
	positions=(int **)calloc(numberofseqs,sizeof(int *));
	for(i=0;i<numberofseqs;i++) positions[i]=(int *)calloc(numpositions[i],sizeof(int));
	currentnode=node;
	while(1){
		while(!isLeaf(currentnode)) currentnode=(currentnode->branches[0]); // starting from this node, find the deepest leftmost leaf node
		for(i=0;i<numberofseqs;i++){
			if(nodeFromSeq(currentnode,i)){
				positions[i][0]=getPosition(currentnode); // store position
				positions[i]++; // increment pointer
			}
		}
		backnode=currentnode;
		n=i;
		while(i==n && backnode!=node){ // find the first branch at the right or above that has not been processed
			currentnode=backnode;
			backnode=(currentnode->backlink);
			n=(backnode->numbranches);
			for(i=0;i<n;i++){ // find which branch of the upper node is the current node
				if((backnode->branches[i])==currentnode){
					i++; // proceed to the next branch (if it is the last one, we will go up)
					break;
				}
			}
		}
		if(i==n) break; // if we have returned to the initial node, we are done
		currentnode=(backnode->branches[i]); // continue down by checking this branch
	}
	for(i=0;i<numberofseqs;i++) positions[i]-=numpositions[i]; // return pointers to their original first position
	return positions;
}

// Creates a new node
treenode * newNode(int i, int j, int depth){
	treenode *newnode=(treenode *)malloc(sizeof(treenode));
	if(newnode==NULL) exitMessage("Not enough memory to create suffix tree!");
	newnode->startpos=i;
	newnode->endpos=j;
	//newnode->count=1;
	newnode->depth=depth;
	newnode->backlink=NULL;
	newnode->suffixlink=NULL;
	newnode->branches=NULL;
	newnode->fromblock=NULL;
	newnode->numbranches=0;
	newnode->id=numberofnodes;
	newnode->isleaf=1;
	newnode->rotation=(infinity-textsize+1);
	newnode->labelfrom=currentseq;
	newnode->leaffrom=currentseq;
	newnode->fromseqs=masks[currentseq];
	newnode->positions=NULL;
	newnode->bordernode=NULL;
	numberofnodes++;
	return newnode;
}


// FALTA: tirar verificacao
// Adds a new branch to an existing node
// NOTE: the branches list in the node is initialized with the alphabet size and new branches are added sequently from position 0 to (alphabetsize-1), no matter what the transition char is
treenode * addBranch(treenode *node, int i, int j){
	treenode *newbranch;
	newbranch=newNode(i,j,(node->depth)+1+(j-i));
	if((node->numbranches)>=alphabetsize) return NULL;
	if(node->branches==NULL)
		node->branches=(treenode **)calloc((alphabetsize),sizeof(treenode *));
	node->branches[node->numbranches]=newbranch;
	node->numbranches++;
	newbranch->backlink=node;
	return newbranch;
}


// FALTA: tirar verificacao
// Splits the node in the given position
// NOTE: the new split node is placed between the current node below and the back node above, and returned from the function
treenode * splitNode(treenode *node, int pos){
	int branchpos;
	treenode *newnode, *backnode;
	if( (node->backlink)==NULL || pos<1 || pos>((node->endpos)-(node->startpos)) ) return NULL; // verificacao
	backnode=node->backlink;
	newnode=newNode(node->startpos,(node->startpos)+pos-1,(backnode->depth)+pos); // novo no' entre o no' anterior e o actual
	newnode->isleaf=0; // atributos do no'
	newnode->rotation=node->rotation;
	newnode->labelfrom=node->labelfrom;
	newnode->leaffrom=-1;
	newnode->fromseqs=node->fromseqs;
	branchpos=0;
	while(backnode->branches[branchpos]!=node) branchpos++; // procurar posicao no no' anterior onde o no' actual esta' ligado
	backnode->branches[branchpos]=newnode; // substitui ligacao do no' anterior para o no' actual por ligacao para o novo no'
	newnode->backlink=backnode; // liga novo no' ao no' anterior
	newnode->branches=(treenode **)malloc((alphabetsize)*sizeof(treenode *)); // inicializa ramos no novo no'
	newnode->branches[0]=node; // liga novo no' ao no' actual
	newnode->numbranches=1;
	node->backlink=newnode; // actualiza apontador do no' actual para o novo no' anterior
	node->startpos=(node->startpos)+pos; // actualiza primeira posicao da label do no' actual
	return newnode;
}

// Follows the node path corresponding to spelling the text from the root and returns the final node
treenode * followText(char *chars){
	int i,j,k,n,textsize,textpos,seq,labelsize;
	char *text,c;
	treenode *node,*branch;
	branch=NULL; // initialize some values to prevent compile warnings
	text=NULL;
	textsize=0;
	textpos=0;
	node=root;
	labelsize=0;
	i=0;
	while(1){
		for(j=1;j<labelsize;j++){ // follow chars inside the node
			c=chars[i];
			if(c=='\0') return node;
			if(textpos>=textsize) textpos-=textsize;
			if(text[textpos]!=c) return NULL;
			textpos++;
			i++;
		}
		c=chars[i];
		if(c=='\0') return node;
		n=node->numbranches;
		for(k=0;k<n;k++){ // find correct branch
			branch=node->branches[k];
			seq=branch->labelfrom;
			text=texts[seq];
			textsize=textsizes[seq];
			textpos=branch->startpos;
			if(textpos>=textsize) textpos-=textsize;
			if(text[textpos]==c) break;
		}
		if(k==n) return NULL;
		node=branch;
		labelsize=(node->endpos)-(node->startpos)+1;
		textpos++;
		i++;
	}
	return NULL;
}

// Checks if the given text is a rotation of a previous sequence
int isRotation(char *text){
	treenode *node;
	node=followText(text);
	if(node==NULL) return 0;
	if(node->isleaf) return ((node->leaffrom)+1);
	return 0;
}

// Returns the node and the position corresponding to the transition by the given character
// NOTE: variables text and textsize are needed to be set for the current sequence
treenode * followChar(treenode *node, int *pos, int charpos){
	treenode *branch;
	char *nodetext;
	char textchar,treechar;
	int i,n,k,seq,labelsize,nodetextsize;
	labelsize=(node->endpos)-(node->startpos)+1;
	if( (*pos)<1 || (*pos)>labelsize ) return NULL;
	if(charpos<textsize) textchar=text[charpos];
	else textchar=text[(charpos-textsize)];
	if(textchar!='A' && textchar!='C' && textchar!='G' && textchar!='T') textchar='-';
	if((*pos)==labelsize){
		n=node->numbranches;
		for(i=0;i<n;i++){
			branch=node->branches[i];
			seq=branch->labelfrom;
			nodetext=texts[seq];
			nodetextsize=textsizes[seq];
			k=branch->startpos;
			if(k>=nodetextsize) k-=nodetextsize;
			treechar=nodetext[k];
			if(treechar!='A' && treechar!='C' && treechar!='G' && treechar!='T') treechar='-';
			if(treechar==textchar){
				(*pos)=1;
				return branch;
			}
		}
		return NULL;
	}
	seq=node->labelfrom;
	nodetext=texts[seq];
	nodetextsize=textsizes[seq];
	k=(node->startpos)+(*pos);
	if(k>=nodetextsize) k-=nodetextsize;
	treechar=nodetext[k];
	if(treechar!='A' && treechar!='C' && treechar!='G' && treechar!='T') treechar='-';
	if(treechar==textchar){
		(*pos)=(*pos)+1;
		return node;
	}
	return NULL;
}


// Returns the node and the position that matches the suffix of the given node/position
// NOTE: position pos goes from 1 to the size of the label
// NOTE: dist is always equal to the label size (passed as argument only to avoid repeated calculation in the calling functions)
treenode * getSuffixNode(treenode *node, int dist, int *pos){
	treenode *backnode;
	int togoup,labelsize;
	backnode=node->suffixlink; // follow the suffix link
	togoup=dist-(*pos); // number of characters to go up in the chain of nodes (dist is always the label size)
	if(node->isleaf) togoup++; // if it is a suffix rotation link then go up one more character
	while( ( labelsize=((backnode->endpos)-(backnode->startpos)+1) ) <= togoup ){
		togoup=togoup-labelsize; // decrement the number of characters in the node's label
		backnode=backnode->backlink; // go up to the next node behind
	}
	(*pos)=labelsize-togoup; // new position in the found node
	return backnode;
}

// Deletes the sequence information from all variables
void discardSequence(int seq){
	int i;
	unsigned long long int maskclear;
	free(texts[seq]);
	free(descs[seq]);
	for(i=(seq+1);i<numberofseqs;i++){
		texts[(i-1)]=texts[i];
		descs[(i-1)]=descs[i];
		textsizes[(i-1)]=textsizes[i];
	}
	maskclear = 1u;
	maskclear = maskclear << (numberofseqs-1);
	maskclear = ~ maskclear;
	allseqsmask = allseqsmask & maskclear;
	numberofseqs--;
	if(numberofseqs<2) exitMessage("The program needs at least 2 sequences to run");
}

// TODO: fix issue of two identical input sequences with different rotations
// TODO: give alphabetsize as input
// Builds the generalized suffix tree of all the given sequences
treenode * buildGeneralizedTree(){
	treenode *currentnode, *foundnode, *newnode, *prevnewnode, *splitnode, *prevsplitnode, *backnode;
	int ii,i,j,n,pos,labelsize,nodedepth;
	int *oldids;
	alphabetsize=ALPHABETSIZE;
	//steps=0;
	oldids=(int *)malloc(numberofseqs*sizeof(int));
	for(i=0;i<numberofseqs;i++) oldids[i]=i;
	numberofnodes=0;
	currentseq=0;
	initializeMasks();
	//printMasks();
	root=newNode(-1,-1,0);
	root->leaffrom=-1;

	for(j=0;j<numberofseqs;j++){ // loop for all the texts

	//if((rotseq=isRotation(texts[j]))) printf("WARNING: Sequence %d is a rotation of sequence %d\n",(j+1),(rotseq));

	currentseq=j;
	text=texts[j];
	textsize=textsizes[j];
	infinity=(textsize-1);
	n=textsize;
	foundnode=NULL;
	newnode=NULL;
	prevnewnode=NULL;
	splitnode=NULL;
	prevsplitnode=NULL;
	currentnode=root;
	markNodeFromSeq(root,j);
	pos=1;
	for(ii=0;ii<(2*n);ii++){ // loop for all chars of the text twice
		i=ii;
		if(ii>=n) i-=n; // if we passed over the end of the sequence, then cycle from the beggining
		while( ( foundnode=followChar(currentnode,&pos,i) ) == NULL ){ // if the transition does not exist yet
			steps++;
			// lets create the transition for this char
			// but first we need to check if the current node needs to be split or not
			labelsize=(currentnode->endpos)-(currentnode->startpos)+1; // size of the node label
			if( pos != labelsize ){ // if we are not at the end of the node label
				splitnode=splitNode(currentnode,pos); // split the node in two nodes
				newnode=addBranch(splitnode,ii,infinity); // add the new char transition
				if(prevsplitnode!=NULL) // link the last split node to this one
					prevsplitnode->suffixlink=splitnode;
				prevsplitnode=splitnode;
			} else { // if we are at the end of the node label
				splitnode=currentnode; // this node was already split at this position
				newnode=addBranch(currentnode,ii,infinity); // add the new char transition
				if(prevsplitnode!=NULL){ // link the last split node to this one
					prevsplitnode->suffixlink=currentnode;
					prevsplitnode=NULL; // this node was not split
				}
			}
			backnode=splitnode; // mark all nodes backward from this node as belonging to this sequence
			while(!nodeFromSeq(backnode,j)){
				markNodeFromSeq(backnode,j);
				backnode=backnode->backlink;
			}
			if(prevnewnode!=NULL) // link the previous created node to the new created node
				prevnewnode->suffixlink=newnode;
			prevnewnode=newnode;
			infinity++; // increase the length of the next created node
			if(currentnode==root){ // if we created the new node at the root, then stop
				newnode->suffixlink=root;
				pos=1;
				break;
			}
			// follow the suffix link and eventually climb to higher nodes in the tree
			currentnode=getSuffixNode(currentnode,labelsize,&pos);
		}
		if(foundnode!=NULL){ // if the char transition already existed
			if(prevsplitnode!=NULL){ // link the last split node to the current node
				prevsplitnode->suffixlink=currentnode;
				prevsplitnode=NULL; // this node was not split
			}
			currentnode=foundnode; // continue from this node forward
			nodedepth=((foundnode->backlink)->depth)+pos;
			if(nodedepth==n){
				newnode=foundnode;
				labelsize=(foundnode->endpos)-(foundnode->startpos)+1;
				if(labelsize!=pos) newnode=splitNode(foundnode,pos);
				if(prevnewnode!=NULL) prevnewnode->suffixlink=newnode;
				prevnewnode=newnode;
				backnode=newnode;
				while(!nodeFromSeq(backnode,j)){
					markNodeFromSeq(backnode,j);
					backnode=backnode->backlink;
				}
				currentnode=getSuffixNode(foundnode,labelsize,&pos);
				if(!(newnode->isleaf)){
					newnode->isleaf=1;
					newnode->rotation=(infinity-textsize+1);
					newnode->leaffrom=currentseq;
					steps++;
				} else if(((newnode->leaffrom)!=currentseq)){
					printf("> WARNING: Discarding seq. %d because it is an identical rotation of seq. %d\n",oldids[currentseq]+1,oldids[(newnode->leaffrom)]+1);
					for(i=j;i<numberofseqs;i++) oldids[i]++;
					discardSequence(currentseq);
					j--;
					break;
				}
				infinity++;
			}
		}
	} // end of the loop for all the chars of the text

	printf(".");fflush(stdout);
	} // end of the loop for all the texts

	free(oldids);
	return root;
}


// imprime a estrutura da arvore de sufixos de todos os nos abaixo do no'
void printSuffixTreeNode(treenode *node, int pos){
	int i,n,slid,blid,seq,textsize;
	char *lf,*text,*seqs;
	if(pos==0){
		seqs=(char *)calloc((numberofseqs+1),sizeof(char));
		for(i=0;i<numberofseqs;i++) seqs[i]=(48+(i+1)%10);
		printf("[NID](DEP)(NON)(LF)(ROT){<-BKL}{->SXL}|%s|\n",seqs);
		free(seqs);
	}
	if((node->backlink)!=NULL) blid=node->backlink->id;
	else blid=-1;
	if((node->suffixlink)!=NULL) slid=node->suffixlink->id;
	else slid=-1;
	//if(node->isleaf) lf='X'; else lf=' ';
	lf=(char *)calloc(3,sizeof(char));
	if((node->leaffrom)==-1) strcpy(lf,"  ");
	else sprintf(lf,"%2d",(node->leaffrom));
	seqs=(char *)calloc((numberofseqs+1),sizeof(char));
	for(i=0;i<numberofseqs;i++){
		if(nodeFromSeq(node,i)) seqs[i]='X';
		else seqs[i]='-';
	}
	printf("[%3d](%3d)(%3d)(%s)(%3d){<-%3d}{->%3d}|%s|",
		node->id,node->depth,node->numbranches,lf,node->rotation,blid,slid,seqs);
	free(seqs);
	free(lf);
	for(i=0;i<pos;i++)	printf("\t");
	if(node!=root){
		n=(node->endpos);
		seq=(node->labelfrom);
		text=texts[seq];
		textsize=textsizes[seq];
		for(i=(node->startpos);i<=n;i++) printf("%c",text[(i%textsize)]);
	}
	printf("\n");
	fflush(stdout);
	n=node->numbranches;
	for(i=0;i<n;i++){
		printSuffixTreeNode(node->branches[i],pos+1);
	}
}


// imprime a estrutura da arvore de sufixos
void printSuffixTree(){
	printSuffixTreeNode(root,0);
}


void printTreeNodeInfo(treenode *currentnode, int pos, int seq, int i){
	treenode *branch;
	int k, n, m, textsize;
	char *text;
	n=(currentnode->numbranches);
	text=texts[seq];
	textsize=textsizes[seq];
	printf("\ntext[%d][%d]='%c' ; pos=%d ; nid=%d\t(d=%d,l=%d,s=",seq,i,((i<textsize)?text[i]:text[(i-textsize)]),pos,(currentnode->id),(currentnode->depth),(currentnode->endpos)-(currentnode->startpos)+1);
	for(k=0;k<numberofseqs;k++) printf("%c",(nodeFromSeq(currentnode,k)?'1':'0'));
	printf(")\n ->");
	for(k=0;k<n;k++){
		branch=(currentnode->branches[k]);
		text=texts[(branch->labelfrom)];
		textsize=textsizes[(branch->labelfrom)];
		m=(branch->startpos);
		printf("\t%d:'%c':",(branch->id),((m<textsize)?text[m]:text[(m-textsize)]));
		for(m=0;m<numberofseqs;m++) printf("%c",(nodeFromSeq(branch,m)?'1':'0'));
	}
	fflush(stdout);
}

// Releases the memory allocated for this node and all child nodes bellow
void freeTreeNode(treenode *node){
	int i;
	linkedpos *pos,*nextpos;
	if((node->branches)!=NULL){
		while((node->numbranches)>0){
			freeTreeNode(node->branches[(node->numbranches)-1]);
			(node->numbranches)--;
		}
		free(node->branches);
	}
	if((node->positions)!=NULL){
		for(i=0;i<numberofseqs;i++){
			pos=(node->positions[i]);
			while(pos!=NULL){
				nextpos=(pos->next);
				free(pos);
				pos=nextpos;
			}
		}
		free(node->positions);
	}
	free(node);
}


// devolve o numero total de nos da arvore
int getNumberOfTreeNodes(){
	return numberofnodes;
}


// devolve o tamanho da memoria ocupada pela arvore
int getTreeSize(){
	return numberofnodes*sizeof(treenode);
}

// TODO: acabar codigo
// actualiza todos os contadores dos nos da arvore
void updateNodeCounts(){
	return;
}

// TODO: remove this function because it is the same as getNodeLabel
// Returns the string spelled by the node path from the root to the given node
char * getNodeText(treenode *node, int *size){
	int i,j,n,textsize,textpos,seq,labelsize;
	char *text,*string;
	treenode *tmpnode;
	if(node==NULL) return NULL;
	n=(node->depth);
	(*size)=n;
	string=(char *)calloc((n+1),sizeof(char));
	tmpnode=node;
	i=(n-1);
	while(tmpnode!=root){
		seq=tmpnode->labelfrom;
		text=texts[seq];
		textsize=textsizes[seq];
		textpos=tmpnode->endpos;
		if(textpos>=textsize) textpos-=textsize;
		labelsize=(node->endpos)-(node->startpos)+1;
		for(j=0;j<labelsize;j++){
			string[i]=text[textpos];
			i--;
			textpos--;
			if(textpos<0) textpos+=textsize;
		}
		tmpnode=tmpnode->backlink;
	}
	return string;
}


// Checks if the suffix tree is well built
int checkSuffixTree(){
	int k,i,j,m,pos,nodeid,labelsize;
	treenode *currentnode, *foundnode, *prevendnode, *suffixnode;
	treenode **blockslist;
	nodeid=0;
	blockslist=(treenode **)calloc(numberofnodes,sizeof(treenode *));
	printf("> Checking cyclic suffix tree consistency");
	fflush(stdout);
	for(k=0;k<numberofseqs;k++){
		currentseq=k;
		text=texts[k];
		textsize=textsizes[k];
		prevendnode=NULL;
		for(i=0;i<=textsize;i++){
			foundnode=NULL;
			currentnode=root;
			pos=1;
			for(j=0;j<textsize;j++){
				m=(i+j);
				//if(m>=textsize) m-=textsize;
				foundnode=followChar(currentnode,&pos,m);
				if(foundnode==NULL){ printf("\nERROR: Forward transition not found [Seq:%02d;Rot:%d;Pos:%d;Node:%d]\n",currentseq,i,j,currentnode->id); return 0; }
				if(foundnode!=currentnode){
					if(!nodeFromSeq(foundnode,currentseq)){ printf("\nERROR: Sequence mark not found [Seq:%02d;Rot:%d;Pos:%d;Node:%d]\n",currentseq,i,j,currentnode->id); return 0; }
					if((foundnode->backlink)!=currentnode){ printf("\nERROR: Backward link does not match [Seq:%02d;Rot:%d;Pos:%d;Node:%d]\n",currentseq,i,j,foundnode->id); return 0; }
					nodeid=foundnode->id;
					if(nodeid!=0){
						suffixnode=foundnode->suffixlink;
						if(suffixnode==NULL){ printf("\nERROR: Suffix link not found [Seq:%02d;Rot:%d;Pos:%d;Node:%d]\n",currentseq,i,j,nodeid); return 0; }
						if(!nodeFromSeq(suffixnode,currentseq)){ printf("\nERROR: Sequence mark not found at suffix [Seq:%02d;Rot:%d;Pos:%d;Node:%d]\n",currentseq,i,j,nodeid); return 0; }
						if(!(foundnode->isleaf) && (suffixnode->depth)!=((foundnode->depth)-1)){ printf("\nERROR: Incorrect suffix link depth [Seq:%02d;Rot:%d;Pos:%d;Node:%d]\n",currentseq,i,j,nodeid); return 0; }
						if((foundnode->isleaf) && (suffixnode->depth)!=((foundnode->depth))){ printf("\nERROR: Incorrect suffix link depth [Seq:%02d;Rot:%d;Pos:%d;Node:%d]\n",currentseq,i,j,nodeid); return 0; }
						labelsize=getLabelSize(foundnode);
						if((foundnode->depth)!=((currentnode->depth)+labelsize)){ printf("\nERROR: Incorrect back link depth [Seq:%02d;Rot:%d;Pos:%d;Node:%d]\n",currentseq,i,j,nodeid); return 0; }
					}
					if(blockslist[nodeid]==NULL) blockslist[nodeid]=foundnode;
					currentnode=foundnode;
				}
			}
			if(prevendnode!=NULL){
				suffixnode=(prevendnode->suffixlink);
				if(suffixnode!=currentnode){ printf("\nERROR: Suffix link does not match at final node [Seq:%02d;Rot:%d;Pos:%d;Node:%d]\n",currentseq,i,j,prevendnode->id); return 0; }
			}
			if((currentnode->depth)!=textsize){ printf("\nERROR: Incorrect final node depth [Seq:%02d;Rot:%d;Pos:%d;Node:%d]\n",currentseq,i,j,nodeid); return 0; }
			if(!(currentnode->isleaf)){ printf("\nERROR: Final mark not found [Seq:%02d;Rot:%d;Pos:%d;Node:%d]\n",currentseq,i,j,currentnode->id); return 0; }
			if((currentnode->leaffrom)!=currentseq){ printf("\nERROR: Incorrect final mark [Seq:%02d;Rot:%d;Pos:%d;Node:%d]\n",currentseq,i,j,currentnode->id); return 0; }
			if((currentnode->rotation)!=i && i!=textsize){ printf("\nERROR: Incorrect rotation [Seq:%02d;Rot:%d;Pos:%d;Node:%d]\n",currentseq,i,j,currentnode->id); return 0; }
			prevendnode=currentnode;
		}
		printf(".");
		fflush(stdout);
	}
	for(i=1;i<numberofnodes;i++){
		if(blockslist[i]==NULL){ printf("\nERROR: Node [%d] not referenced.\n",i); return 0; }
	}
	free(blockslist);
	printf(" OK\n");
	return 1;
}


// corre alguns testes com strings pre-definidas
void runTests(){
	int i;
	treenode *tree;
	char * testarray[] = {"AAAB","AABA","ABAB","ABAC","AABBC","AAAAAA","AABABC","ABCABD","ABCDAE","CABACA","BDABCABD",
						"CABACADA","ABCABDABC","BDABCABDA","ABCABDABCA","ABCABDABCAB","ABCDABCEABF","BCDABCEABCF",
						"BFBCEABCDABK","ABCABDABEABEF","ABCABDABCEABDE","BEABCDFABCDABK","DGABCABDABEABEFDBH"};
	texts=testarray;
	numberofseqs=23;
	textsizes=(int *)calloc(numberofseqs,sizeof(int));
	for(i=0;i<numberofseqs;i++) textsizes[i]=(int)strlen(texts[i]);
	tree=buildGeneralizedTree();
	printSuffixTree();
	if(checkSuffixTree()) printf("OK\n");
	free(textsizes);
	freeTreeNode(tree);
	/*
	int testresultok,i,n=23;
	for(i=0;i<n;i++){
		printf("> Test %2d/%2d (\"%s\"): ",(i+1),n,texts[i]);
		tree=buildTree(texts[i],NULL);
		testresultok=checkSuffixTree();
		if(testresultok) printf("OK\n");
		freeTreeNode(tree);
	}
	*/
}

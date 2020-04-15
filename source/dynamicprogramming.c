#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "dynamicprogramming.h"
#include "csamsa.h"
#include "alignment.h"
#include "alignmentmap.h"

#define ALPHABET "ACGT-"
#define ALPHABETSIZE 4
#define GAPCODE 4
#define NONGAPCODE 5

// TODO: change DP direction selection order in ProgressiveDP() function
/**/
#define MATCHSCORE (+1)
#define DOUBLEGAPSCORE (0)
#define MISMATCHSCORE (-1)
#define INDELSCORE (-1)
/**/
/*
#define MATCHSCORE (+1)
#define DOUBLEGAPSCORE (-1)
#define MISMATCHSCORE (0)
#define INDELSCORE (-1)
*/

int **dpmatrix=NULL; // dynamic programming matrix
int **scorevector=NULL; // array/matrix of the size of the current consensus size containing the counts of each char on each column (NOTE: starts at position 1 ?)
int *orderedseqs=NULL; // sorted sequences by the order they will be aligned
int *seqlengths=NULL; // length of each of the ordered sequences that are going to be aligned
int prevnrows=-1;
int prevncols=-1;
/*
int **matchscore=NULL;
int *scores=NULL;
int prevndpcols=-1;
int prevndprows=-1;
int ndpcols=0;
int ndprows=0;
int dprowsize=0;
int firstseq=0;
*/

// Returns the score of aligning char code a against char code b
int Score(int a, int b){
	if(a==GAPCODE && b==GAPCODE) return 0; // gap aligned against gap does nothing
	if(a==b) return 1; // match score = 1
	return -1; // mismatch = insertion = deletion = -1
	/*
	if(a==b) return 1; // match score = 1
	return 0; // mismatch = insertion = deletion = 0
	*/
}

// Retrieves the index code of the character in the given position of the given rotated sequence
int CharCodeFromSeq(int pos, int seq){
	char c;
	int i,code;
	i=rotations[seq]+pos;
	if(i>=textsizes[seq]) i-=textsizes[seq];
	c=texts[seq][i];
	switch(c){
		case 'A' : code=0; break;
		case 'C' : code=1; break;
		case 'G' : code=2; break;
		case 'T' : code=3; break;
		default : code=-1; break;
	}
	return code;
}

// Retrieves the index code of the character
int GetCharCode(char c){
	int code;
	switch(c){
		case 'A' : code=0; break;
		case 'C' : code=1; break;
		case 'G' : code=2; break;
		case 'T' : code=3; break;
		case '-' : code=4; break;
		default : code=-1; break;
	}
	return code;
}

// TODO: print chars, scorevector counts and directions
void PrintMatrix(int **matrix, char **dirs, int vtextstart, int vtextid, int **hscores, int nrows, int ncols){
	int j,k;
	char c;
	printf("\n");
	/*
	for(k=0;k<=ncols;k++){
		if(k==0) printf("   ");
		printf("| %3d ",k);
		if(k==ncols) printf("|\n");
	}
	*/
	for(j=0;j<(ALPHABETSIZE+1);j++){
		switch(j){
			case 0: c='A'; break;
			case 1: c='C'; break;
			case 2: c='G'; break;
			case 3: c='T'; break;
			case 4: c='-'; break;
		}
		printf(" %c ",c);
		for(k=0;k<=ncols;k++){
			if(k==0) printf("|     ");
			else printf("| %2d ",hscores[k][j]);
		}
		printf(" |\n");
	}
	for(k=0;k<=ncols;k++){
		if(k==0) printf("    ");
		printf("-----");
		if(k==ncols) printf("-\n");
	}
	for(j=0;j<=nrows;j++){
		for(k=0;k<=ncols;k++){
			//if(k==0) printf(" %2d |",j);
			if(k==0){
				if(j==0) printf("   | ");
				else printf(" %c | ",CharAt(vtextstart+j-1,vtextid));
			}
			switch(dirs[j][k]){
				case 'D': c='\\'; break;
				case 'L': c=(char)27; break;
				case 'U': c=(char)24; break;
				default: c=' '; break;
			}
			printf("%c% 3d ",c,matrix[j][k]);
		}
		printf("|\n");
	}
}

/*
void BuildTwoWaysMatchScoreMatrix(alignmapsegment *segment){
	int **forwardscore=NULL;
	int **reversescore=NULL;
	int **currentscore=NULL;
	int *mrow,*frow,*rrow;
	int istart,iend,fr;
	int i,j;
	int n,m,k,pos,endpos,startpos,charcode,score;
	if(forwardscore==NULL){
		forwardscore=(int **)calloc(ndpcols,sizeof(int));
		for(j=0;j<ndpcols;j++) forwardscore[j]=(int *)calloc(ALPHABETSIZE,sizeof(int));
	}
	if(reversescore==NULL){
		reversescore=(int **)calloc(ndpcols,sizeof(int));
		for(j=0;j<ndpcols;j++) reversescore[j]=(int *)calloc(ALPHABETSIZE,sizeof(int));
	}
	for(i=0;i<ndpcols;i++){
		mrow=matchscore[i];
		frow=forwardscore[i];
		rrow=reversescore[i];
		for(j=0;j<ALPHABETSIZE;j++){
			frow[j]=mrow[j];
			rrow[j]=mrow[j];
		}
	}
	fr=0;
	istart=firstseq;
	iend=(numberofseqs-1);
	currentscore=forwardscore;
	i=istart;
	while(1){
		n=orderedseqs[i];
		if(seqlengths[i]==0) continue;
		startpos=(segment->positions[n])+(segment->size);
		endpos=(segment->next->positions[n]);
		pos=startpos;
		ndprows=seqlengths[i];
		dprowsize=(ndpcols-ndprows+1);
		for(j=1;j<=ndprows;j++){
			charcode=CharCodeFromSeq(pos,n);
			for(k=1;k<=dprowsize;k++){
				score=currentscore[(j+k-2)][charcode];
				score+=dpmatrix[j-1][k];
				if((dpmatrix[j][k-1])>score) dpmatrix[j][k]=dpmatrix[j][k-1];
				else dpmatrix[j][k]=score;
			}
			pos++;
		}
		j=ndprows;
		k=dprowsize;
		m=(ndpcols-1);
		pos=(endpos-1);
		while(m>=0){
			if((dpmatrix[j][k])==(dpmatrix[j][k-1])){
				k--;
			} else {
				charcode=CharCodeFromSeq(pos,n);
				currentscore[m][charcode]++;
				pos--;
				j--;
			}
			m--;
		}
		if(fr==0){
			i++;
			if(i>iend){
				fr=1;
				istart=(numberofseqs-1);
				iend=firstseq;
				currentscore=reversescore;
				i=istart;
			}
		} else {
			i--;
			if(i<iend) break;
		}
	}
	for(i=0;i<ndpcols;i++){
		mrow=matchscore[i];
		frow=forwardscore[i];
		rrow=reversescore[i];
		for(j=0;j<ALPHABETSIZE;j++){
			mrow[j]=(frow[j]+rrow[j]);
		}
		free(forwardscore[i]);
		free(reversescore[i]);
	}
	free(forwardscore);
	free(reversescore);
	forwardscore=NULL;
	reversescore=NULL;
}

int NextBestSeqForDP(alignmapsegment *segment){
	int i,j,maxval,maxpos;
	int n,k,pos,endpos,startpos,charcode,score;
	if(scores==NULL) scores=(int *)calloc(numberofseqs,sizeof(int));
	for(i=firstseq;i<numberofseqs;i++){
		if(scores[i]==-1) continue;
		n=orderedseqs[i];
		if(seqlengths[i]==0) continue;
		startpos=(segment->positions[n])+(segment->size);
		endpos=(segment->next->positions[n]);
		pos=startpos;
		ndprows=seqlengths[i];
		dprowsize=(ndpcols-ndprows+1);
		for(j=1;j<=ndprows;j++){
			charcode=CharCodeFromSeq(pos,n);
			for(k=1;k<=dprowsize;k++){
				score=matchscore[(j+k-2)][charcode];
				score+=dpmatrix[j-1][k];
				if((dpmatrix[j][k-1])>score) dpmatrix[j][k]=dpmatrix[j][k-1];
				else dpmatrix[j][k]=score;
			}
			pos++;
		}
		scores[i]=dpmatrix[ndprows][dprowsize];
	}
	maxpos=firstseq;
	maxval=scores[firstseq];
	//printf("[");
	for(i=(firstseq+1);i<numberofseqs;i++){
		//printf("%d,",scores[i]);
		if(scores[i]>maxval){
			maxpos=i;
			maxval=scores[i];
		}
	}
	//printf("](%d@%d)\n",maxval,(maxpos-firstseq));
	scores[maxpos]=-1;
	return maxpos;
}
*/

//TODO: set prev* variables global
//TODO: get max length and alloc dpmatrix only once
// Sorts the sequences to be aligned by their increasing length
void SortSequencesForDP(alignmapsegment *segment){
	int i,j,ii,jj,n;
	int nrows,ncols;
	int startposi,startposj,posi,posj,charcodei,charcodej;
	int left,up,diag,score;
	int maxi,maxj,maxscore;
	int min,minpos,aux;
	int **pairwisescoresmatrix;
	if(orderedseqs==NULL) orderedseqs=(int *)calloc(numberofseqs,sizeof(int));
	if(seqlengths==NULL) seqlengths=(int *)calloc(numberofseqs,sizeof(int));
	for(i=0;i<numberofseqs;i++){ // initialize arrays
		orderedseqs[i]=i;
		seqlengths[i]=(segment->next->positions[i])-((segment->positions[i])+(segment->size));
	}
	for(i=0;i<(numberofseqs-1);i++){ // sort sequences by their length (from shortest to longest)
		min=seqlengths[i];
		minpos=i;
		for(j=i+1;j<numberofseqs;j++){ // find shortest sequence ahead
			if(seqlengths[j]<min){
				min=seqlengths[j];
				minpos=j;
			}
		}
		if(minpos!=i){ // switch sequences in both orderedseqs and seqlengths
			aux=orderedseqs[i];
			orderedseqs[i]=orderedseqs[minpos];
			orderedseqs[minpos]=aux;
			aux=seqlengths[i];
			seqlengths[i]=seqlengths[minpos];
			seqlengths[minpos]=aux;
		}
	}
	return;
	/**/
	// initialize pairwise scores matrix
	pairwisescoresmatrix=(int **)calloc(numberofseqs,sizeof(int *));
	for(i=0;i<numberofseqs;i++){
		pairwisescoresmatrix[i]=(int *)malloc(numberofseqs*sizeof(int));
		pairwisescoresmatrix[i][i]=0;
	}
	dpmatrix=NULL;
	prevnrows=-1;
	prevncols=-1;
	// calculate all pairwise scores
	for(i=0;i<(numberofseqs-1);i++){
		nrows=seqlengths[i];
		for(j=i+1;j<numberofseqs;j++){
			ncols=seqlengths[j];
			// realloc dpmatrix if needed
			if(nrows>prevnrows || ncols>prevncols){
				if(dpmatrix!=NULL){
					for(ii=0;ii<=prevnrows;ii++) free(dpmatrix[ii]);
					free(dpmatrix);
					dpmatrix=NULL;
				}
				prevnrows=nrows;
				dpmatrix=(int **)malloc((nrows+1)*sizeof(int));
				for(ii=0;ii<=nrows;ii++){
					dpmatrix[ii]=(int *)malloc((ncols+1)*sizeof(int));
					dpmatrix[ii][0]=ii*(-1);
				}
				for(jj=0;jj<=ncols;jj++){
					dpmatrix[0][jj]=jj*(-1);
				}
			}
			startposi=(segment->positions[i])+(segment->size);
			startposj=(segment->positions[j])+(segment->size);
			posi=startposi;
			for(ii=1;ii<=nrows;ii++){
				charcodei=CharCodeFromSeq(posi,i);
				posj=startposj;
				for(jj=1;jj<=ncols;jj++){
					charcodej=CharCodeFromSeq(posj,j);
					if(charcodei==charcodej) score=1;
					else score=-1;
					diag=dpmatrix[ii-1][jj-1]+score;
					up=dpmatrix[ii-1][jj]-1;
					left=dpmatrix[ii][jj-1]-1;
					if(up>left && up>diag){
						dpmatrix[ii][jj]=up;
						continue;
					}
					if(left>up && left>diag){
						dpmatrix[ii][jj]=left;
						continue;
					}
					dpmatrix[ii][jj]=diag;
					posj++;
				}
				posi++;
			}
			score=dpmatrix[nrows][ncols];
			pairwisescoresmatrix[i][j]=score;
			pairwisescoresmatrix[j][i]=score;
			prevncols=ncols;
		}
		if(nrows>prevnrows) prevnrows=nrows;
		if(ncols>prevncols) prevncols=ncols;
		printf(":");fflush(stdout);
	}
	// set best first two sequences
	maxi=-1;
	maxj=-1;
	maxscore=INT_MIN;
	for(i=0;i<(numberofseqs-1);i++){
		for(j=i+1;j<numberofseqs;j++){
			score=pairwisescoresmatrix[i][j];
			if(score>maxscore){
				maxi=i;
				maxj=j;
				maxscore=score;
			}
		}
	}
	orderedseqs[0]=maxi;
	orderedseqs[1]=maxj;
	pairwisescoresmatrix[maxi][maxi]=INT_MIN;
	pairwisescoresmatrix[maxj][maxj]=INT_MIN;
	for(j=0;j<numberofseqs;j++){
		if(j!=maxi && j!=maxj){
				pairwisescoresmatrix[j][j]=pairwisescoresmatrix[j][maxi];
				pairwisescoresmatrix[j][j]+=pairwisescoresmatrix[j][maxj];
		}
	}
	// set following sequences by score
	for(i=2;i<numberofseqs;i++){
		maxscore=INT_MIN;
		maxj=-1;
		for(j=0;j<numberofseqs;j++){
			score=pairwisescoresmatrix[j][j];
			if(score>maxscore){
				maxj=j;
				maxscore=score;
			}
		}
		orderedseqs[i]=maxj;
		pairwisescoresmatrix[maxj][maxj]=INT_MIN;
		for(j=0;j<numberofseqs;j++){
			if(pairwisescoresmatrix[j][j]!=INT_MIN)
				pairwisescoresmatrix[j][j]+=pairwisescoresmatrix[j][maxj];
		}
	}
	// update sequences lengths
	for(i=0;i<numberofseqs;i++){
		n=orderedseqs[i];
		seqlengths[i]=(segment->next->positions[n])-((segment->positions[n])+(segment->size));
	}
	// clean variables
	if(pairwisescoresmatrix!=NULL){
		for(i=0;i<numberofseqs;i++) free(pairwisescoresmatrix[i]);
		free(pairwisescoresmatrix);
		pairwisescoresmatrix=NULL;
	}
	if(dpmatrix!=NULL){
		for(ii=0;ii<=nrows;ii++) free(dpmatrix[ii]);
		free(dpmatrix);
		dpmatrix=NULL;
	}
	/**/
}

/*
void ReSortNextSequence(int i, int ncols,alignmapsegment *segment){
	int nrows,ii,j,k,n;
	int pos,startpos,charcode,diag,up,left,score,rowgapscore,colgapscore;
	int max,maxpos,aux;
	int *scores;
	if(i<=1 || i==(numberofseqs-1)) return;
	nrows=seqlengths[i];
	for(j=(i+1);j<numberofseqs;j++)
		if(seqlengths[j]>nrows) nrows=seqlengths[j];
	if(ncols>prevncols || nrows>prevnrows){
		if(dpmatrix!=NULL){
			for(j=0;j<=prevnrows;j++) free(dpmatrix[j]);
			free(dpmatrix);
			dpmatrix=NULL;
		}
		prevnrows=nrows;
		rowgapscore=i;
		colgapscore=1;
		dpmatrix=(int **)malloc((nrows+1)*sizeof(int *));
		for(j=0;j<=nrows;j++){
			dpmatrix[j]=(int *)malloc((ncols+1)*sizeof(int));
			dpmatrix[j][0]=j*(-rowgapscore);
		}
		colgapscore=0;
		for(j=1;j<=ncols;j++){
			colgapscore+=(i-(scorevector[j][GAPCODE]));
			dpmatrix[0][j]=(-colgapscore);
		}
	}
	if(nrows>prevnrows) prevnrows=nrows;
	if(ncols>prevncols) prevncols=ncols;
	scores=(int *)malloc(numberofseqs*sizeof(int));
	for(ii=i;ii<numberofseqs;ii++){
		n=orderedseqs[ii];
		nrows=seqlengths[ii];
		startpos=(segment->positions[n])+(segment->size);
		pos=startpos;
		for(j=1;j<=nrows;j++){
			charcode=CharCodeFromSeq(pos,n);
			for(k=1;k<=ncols;k++){
				score=(2*(scorevector[k][charcode])-i);
				rowgapscore=i;
				colgapscore=(i-(scorevector[k][GAPCODE]));
				diag=dpmatrix[j-1][k-1]+score;
				up=dpmatrix[j-1][k]-rowgapscore;
				left=dpmatrix[j][k-1]-colgapscore;
				if(up>left && up>diag){
					dpmatrix[j][k]=up;
					continue;
				}
				if(left>up && left>diag){
					dpmatrix[j][k]=left;
					continue;
				}
				dpmatrix[j][k]=diag;
			}
			pos++;
		}
		scores[ii]=dpmatrix[nrows][ncols];
		printf(":");fflush(stdout);
	}
	max=scores[i];
	maxpos=i;
	for(j=(i+1);j<numberofseqs;j++){
		if(scores[j]>max){
			max=scores[j];
			maxpos=j;
		}
	}
	if(maxpos!=i){
		aux=orderedseqs[i];
		orderedseqs[i]=orderedseqs[maxpos];
		orderedseqs[maxpos]=aux;
		aux=seqlengths[i];
		seqlengths[i]=seqlengths[maxpos];
		seqlengths[maxpos]=aux;
	}
	free(scores);
}
*/

void PrintStrings(char **strings, int startpos, int size, int dir, int *validseqs, int numvalidseqs, int *seqstoshift, int numseqstoshift, int shiftsize){
	int n,k,i,shift,endpos,firstpostoshift,lastpostoshift;
	char *string;
	if(shiftsize<0) shiftsize=(-shiftsize); // absolute value of the shift
	if(dir<0){ // define correct starting and ending positions according to direction
		startpos=(startpos-size+1);
		endpos=startpos;
		dir=(-1);
	} else {
		startpos=startpos;
		endpos=(startpos+size-1);
		dir=(+1);
	}
	fprintf(debugfile,"dir=%+d ; shiftsize=%d ; startpos=%d ; endpos=%d\n",dir,shiftsize,startpos,endpos);
	for(n=0;n<numvalidseqs;n++){
		k=validseqs[n];
		string=strings[k];
		fprintf(debugfile,"%2d",k);
		for(i=0;i<numseqstoshift;i++){
			if(k==seqstoshift[i]) break;
		}
		if(i!=numseqstoshift){ // if this is a sequence to shift
			if(dir>0){ // find the first gap ahead (all chars before that are to be shifted to the right)
				firstpostoshift=startpos;
				lastpostoshift=startpos;
				while((lastpostoshift<=endpos) && (string[lastpostoshift]!='\0') && (string[lastpostoshift]!='-')) lastpostoshift++;
				lastpostoshift--; // the last checked position in not valid
				lastpostoshift+=shiftsize; // add gaps at the beggining
			} else { // find the first gap behind (all chars before that are to be shifted to the left)
				lastpostoshift=endpos;
				firstpostoshift=endpos;
				while((firstpostoshift>=startpos) && (firstpostoshift>=0) && (string[firstpostoshift]!='-')) firstpostoshift--;
				firstpostoshift++; // the last checked position in not valid
				firstpostoshift-=shiftsize; // add gaps at the end
			}
			shift=(dir*shiftsize);
			fprintf(debugfile,"*\t");
		} else {
			firstpostoshift=-1;
			lastpostoshift=-1;
			shift=0;
			fprintf(debugfile," \t");
		}
		fprintf(debugfile,"...%c\t",((startpos>0)?(string[(startpos-1)]):('|'))); // previous character
		for(i=startpos;i<=endpos;i++){
			if((shift!=0) && (i>=firstpostoshift) && (i<=lastpostoshift)){
				if((i-shift)<startpos || (i-shift)>endpos) fprintf(debugfile,"%c   ",'-'); // print gap char
				else fprintf(debugfile,"%c   ",string[(i-shift)]); // print shifted char
			}
			else fprintf(debugfile,"%c   ",string[i]); // print un-shifted char
		}
		fprintf(debugfile,"\t%c...",((string[(endpos+1)]!='\0')?(string[(endpos+1)]):('|'))); // next character
		fprintf(debugfile,"\n");
	}
	//fprintf(debugfile,"\n");
}

// TODO: fix score calculation
int GetCountsVectorColumnScore(int **vector, int col){
	int n,totalcount,score;
	totalcount=0;
	for(n=0;n<=ALPHABETSIZE;n++) totalcount += (vector[col][n]);
	score=0;
	/**/
	for(n=0;n<ALPHABETSIZE;n++){
		score += MATCHSCORE * (vector[col][n]) * (vector[col][n]-1); // #matches = nc*(nc-1)
		score += MISMATCHSCORE * (vector[col][n]) * (totalcount-(vector[col][n])-(vector[col][GAPCODE])); // #mismatches = nc*(N-(nc+ng))
		score += INDELSCORE * (vector[col][n]) * (vector[col][GAPCODE]); // #deletions = nc*ng
	}
	score += DOUBLEGAPSCORE * (vector[col][GAPCODE]) * (vector[col][GAPCODE]-1); // #doublegaps = ng*(ng-1)
	score += INDELSCORE * (vector[col][GAPCODE]) * (totalcount-(vector[col][GAPCODE])); // #insertions = ng*(N-ng)
	/**/
	/*
	for(n=0;n<ALPHABETSIZE;n++){
		if(vector[col][n]!=0)
			score += vector[col][n]*( MATCHSCORE*(vector[col][n]-1) + MISMATCHSCORE*(totalcount-vector[col][n]) ); // = (#)*(MATCH*(score[char]-1)+MISMATCH*(N-score[char]))
	}
	if(vector[col][GAPCODE]!=0){
		score += vector[col][GAPCODE]*( DOUBLEGAPSCORE*(vector[col][GAPCODE]-1) + INDELSCORE*(totalcount-vector[col][GAPCODE]) ); // = (#)*(MATCH*(score[gaps]-1)+MISMATCH*(N-score[gaps]))
	}
	*/
	return score;
}

void PrintCountsVector(int **vector, int size, int startpos, int dir){
	int i,n,totalscore;
	fprintf(debugfile,"\t\t\t");
	for(i=0;i<size;i++) fprintf(debugfile,"%-3d ",(startpos==0)?(startpos+i):(startpos+i-1));
	fprintf(debugfile,"\n");
	for(n=0;n<=ALPHABETSIZE;n++){
		fprintf(debugfile,"%c\t\t\t",ALPHABET[n]);
		if(dir>0) for(i=0;i<=(size-1);i++) fprintf(debugfile,"%02d  ",vector[(startpos+i)][n]); // forward direction
		else { // reverse direction
			if(startpos==0) for(i=(size-1);i>=0;i--) fprintf(debugfile,"%02d  ",vector[i][n]); // if it is a reverse temporary vector, print in reverse
			else for(i=(size-1);i>=0;i--) fprintf(debugfile,"%02d  ",vector[(startpos-i)][n]);
		}
		fprintf(debugfile,"\n");
	}
	totalscore=0;
	fprintf(debugfile,"%c\t\t\t",'=');
	if(dir>0){
		for(i=0;i<=(size-1);i++){
			n=GetCountsVectorColumnScore(vector,(startpos+i));
			fprintf(debugfile,"%-3d ",n);
			totalscore+=n;
		}
	} else {
		for(i=(size-1);i>=0;i--){
			if(startpos==0) n=GetCountsVectorColumnScore(vector,i);
			else n=GetCountsVectorColumnScore(vector,(startpos-i));
			fprintf(debugfile,"%-3d ",n);
			totalscore+=n;
		}
	}
	fprintf(debugfile,"= %d\n\n",totalscore);
}

// TODO: **** do not calculate scores of the empty columns at the beginning
// TODO: **** fix score (two times):  = moving[char]*(MATCH*(score[char]-1)+MISMATCH*(N-(score[char]+score[gaps]))+INDEL*(score[gaps]))
// TODO: *** remove empty columns in all the maxaffectedpos area
// TODO: *** check if signal operations are ok
// TODO: *** skip un-necessary columns after a check ?
// TODO: use static variables
// TODO: !!! implement non-linear merge method
void DeleteGappedColumns(int *usableseqs, char **alignedstrings, int numseqs, int *consize, int maxnongaps){
	int i,j,k,n,m,col,mingaps,ii,jj;
	int testedcols, deletedcols, scoregain, colscore;
	int ntoshift, postofarthestgap, minnextgaps, maxposaffected;
	int charcode, currentscore, bestscore, bestshift;
	int looplimit, dirsignal;
	int bestmaxposaffected;
	static int tempsvsize = 0;
	static int auxsvcolsize = 0;
	static int **staticsv = NULL;
	static int **movingsv = NULL;
	static int **workingsv = NULL;
	static int **bestworkingsv = NULL;
	static int **auxsvcol = NULL;
	static int *shiftedscore = NULL;
	static int currentmaxnongaps = 0;
	static int *nposaffected = NULL;
	static int *seqstoshift = NULL;
	static int *postonextgap = NULL;
	static int *nnextgaps = NULL;
	static int *bestnposaffected = NULL;
	bestmaxposaffected=0;
	testedcols=0;
	deletedcols=0;
	scoregain=0;
	mingaps=(numseqs-maxnongaps);
	if(maxnongaps>currentmaxnongaps){
		nposaffected=(int *)realloc(nposaffected,maxnongaps*sizeof(int));
		seqstoshift=(int *)realloc(seqstoshift,maxnongaps*sizeof(int));
		postonextgap=(int *)realloc(postonextgap,maxnongaps*sizeof(int));
		nnextgaps=(int *)realloc(nnextgaps,maxnongaps*sizeof(int));
		bestnposaffected=(int *)realloc(bestnposaffected,maxnongaps*sizeof(int));
		currentmaxnongaps=maxnongaps;
	}
	for(col=1;col<=(*consize);col++){ // check all columns of the aligned strings
		if(scorevector[col][GAPCODE]<mingaps) continue; // if the column has enought aligned chars, leave it alone (scorevector starts at 1)
		testedcols++;
		ntoshift=0;
		for(i=0;i<numseqs;i++){ // mark sequences that we are going to shift
			ii=usableseqs[i];
			if((alignedstrings[ii][(col-1)])!='-'){ // sequences that do not have gaps on this column (alignedstrings starts at 0)
				seqstoshift[ntoshift]=ii;
				ntoshift++;
			}
		}
		if(ntoshift==0){
			printf("!"); // in this column all chars are gaps
			continue;
		}
		bestscore=0; // best score among all the shifts
		bestshift=0; // number of positions to shift to achieve the best score
		looplimit=((*consize)+1); // end of sequences (how far can we go to the right)
		dirsignal=+1; // start checking to the right
		while(1){ // testing all possible shifts in both directions and find the one with the best score
			postofarthestgap=0;
			minnextgaps=(*consize);
			for(k=0;k<ntoshift;k++){
				i=seqstoshift[k];
				j=col;
				postonextgap[k]=0;
				while((j!=looplimit)&&((alignedstrings[i][(j-1)])!='-')){ // find end of non-gapped block
					postonextgap[k]++;
					j=j+dirsignal;
				}
				if(j==looplimit) break; // if the block is at the end of the sequence, we cannot shift it
				if(postonextgap[k]>postofarthestgap) postofarthestgap=postonextgap[k]; // maximum run of consecutive non-gap chars
				nnextgaps[k]=0;
				while((j!=looplimit)&&((alignedstrings[i][(j-1)])=='-')){ // count number of gaps that follow the block
					nnextgaps[k]++;
					j=j+dirsignal;
				}
				if(nnextgaps[k]<minnextgaps) minnextgaps=nnextgaps[k]; // keep the minimum number of gaps that we can shift
			}
			if(k!=ntoshift){ // if one of the blocks was at the end
				if(dirsignal==-1) break;
				looplimit=0; // try the opposite direction
				dirsignal=-1;
				continue;
			}
			for(k=0;k<ntoshift;k++) nposaffected[k]=(postonextgap[k]+minnextgaps); // number of positions on each sequence that will be affected by the shifting of the block
			maxposaffected=(postofarthestgap+minnextgaps);
			if(maxposaffected>tempsvsize){ // realloc the arrays if needed
				staticsv=(int **)realloc(staticsv,maxposaffected*sizeof(int *));
				movingsv=(int **)realloc(movingsv,maxposaffected*sizeof(int *));
				workingsv=(int **)realloc(workingsv,maxposaffected*sizeof(int *));
				bestworkingsv=(int **)realloc(bestworkingsv,maxposaffected*sizeof(int *));
				shiftedscore=(int *)realloc(shiftedscore,maxposaffected*sizeof(int));
				for(k=tempsvsize;k<maxposaffected;k++){
					staticsv[k]=(int *)malloc((ALPHABETSIZE+1)*sizeof(int));
					movingsv[k]=(int *)malloc((ALPHABETSIZE+1)*sizeof(int));
					workingsv[k]=(int *)malloc((ALPHABETSIZE+1)*sizeof(int));
					bestworkingsv[k]=(int *)malloc((ALPHABETSIZE+1)*sizeof(int));
				}
				tempsvsize=maxposaffected; // current size of the arrays
			}
			currentscore=0; // calculate the score of the moving vector in its current original position
			for(j=0;j<maxposaffected;j++){ // process all affected positions
				jj=(col+dirsignal*j); // column index in the main full scorevector
				for(n=0;n<(ALPHABETSIZE+1);n++){ // copy the scorevector counts to the current working arrays
					staticsv[j][n]=scorevector[jj][n]; // the static vector will contain the counts for all other sequences except the shifting ones
					movingsv[j][n]=0; // the moving vector will contain only the counts for the shifting sequences
				}
				for(k=0;k<ntoshift;k++){ // process all shift-able sequences
					if(j<nposaffected[k]){ // if this column still contains the block (and the minimum gaps in front) of the current sequence
						i=seqstoshift[k];
						charcode=GetCharCode(alignedstrings[i][(jj-1)]); // get the char of this seq in this column
						movingsv[j][charcode]++; // add it to the moving vector
						staticsv[j][charcode]--; // remove it from the static vector
					}
				}
				colscore=0; // score of the current column
				for(n=0;n<ALPHABETSIZE;n++){
					if(movingsv[j][n]!=0)
						colscore += movingsv[j][n]*( MATCHSCORE*(scorevector[jj][n]-1) + MISMATCHSCORE*(numseqs-(scorevector[jj][n]+scorevector[jj][GAPCODE])) +  INDELSCORE*(scorevector[jj][GAPCODE]) ); // = moving*(MATCH*(score-1)+MISMATCH*(N-score))
				}
				if(movingsv[j][GAPCODE]!=0)
						colscore += movingsv[j][GAPCODE]*( DOUBLEGAPSCORE*(scorevector[jj][GAPCODE]-1) + INDELSCORE*(numseqs-scorevector[jj][GAPCODE]) ); // = moving*(MATCH*(score-1)+MISMATCH*(N-score))
				currentscore+=colscore;
			} // end of calculating the current score
			for(i=1;i<=minnextgaps;i++){ // simulate the shift of the moving vector each one of the minnextgaps positions
				shiftedscore[(i-1)]=0; // array with the scores resulting from shifting the moving vector
				for(k=0;k<ntoshift;k++){
					j=(nposaffected[k]-1);
					movingsv[j][GAPCODE]--; // delete the gap at the end of each sequence
					nposaffected[k]--;
				}
				for(j=0;j<maxposaffected;j++){ // process all affected columns
					colscore=0;
					if(j<i){ // the places where the blocks were before shifting are now filled with gaps (at least on all shifted sequences)
						for(n=0;n<ALPHABETSIZE;n++) workingsv[j][n]=0;
						workingsv[j][GAPCODE]=(staticsv[j][GAPCODE]+ntoshift); // add a gap to the beginning of all seqs that moved
						if(workingsv[j][GAPCODE]==numseqs) continue; // do not calculate the score of totally empty (all gaps) columns because they will be removed
						colscore += ntoshift*( DOUBLEGAPSCORE*(workingsv[j][GAPCODE]-1) + INDELSCORE*(numseqs-workingsv[j][GAPCODE]) ); // = (#moved)*(MATCH*(score[gaps]-1)+MISMATCH*(N-score[gaps]))
						shiftedscore[(i-1)]+=colscore;
						continue;
					}
					for(n=0;n<ALPHABETSIZE;n++) workingsv[j][n]=(staticsv[j][n]+movingsv[(j-i)][n]); // chars that were already there plus chars that were shifted to this column now
					workingsv[j][GAPCODE]=(staticsv[j][GAPCODE]+movingsv[(j-i)][GAPCODE]);
					if(workingsv[j][GAPCODE]==numseqs) continue; // do not calculate the score of totally empty (all gaps) columns because they will be removed
					for(n=0;n<ALPHABETSIZE;n++){ // update score of the shifted sequences
						if(movingsv[(j-i)][n]!=0)
							colscore += movingsv[(j-i)][n]*( MATCHSCORE*(workingsv[j][n]-1) + MISMATCHSCORE*(numseqs-(workingsv[j][n]+workingsv[j][GAPCODE])) + INDELSCORE*(workingsv[j][GAPCODE]) ); // = (#moved)*(MATCH*(score[char]-1)+MISMATCH*(N-score[char]))
					}
					if(movingsv[(j-i)][GAPCODE]!=0)
							colscore += movingsv[(j-i)][GAPCODE]*( DOUBLEGAPSCORE*(workingsv[j][GAPCODE]-1) + INDELSCORE*(numseqs-workingsv[j][GAPCODE]) ); // = (#movedgaps)*(MATCH*(score[gaps]-1)+MISMATCH*(N-score[gaps]))
					shiftedscore[(i-1)]+=colscore; // update score of this shift
				}
				shiftedscore[(i-1)]-=currentscore; // decrease score of the un-shifted blocks (if it is the same as the currentscore, it will be equal to zero)
				if(shiftedscore[(i-1)]>=bestscore){ // bestscore is equal to zero initially
					bestshift=dirsignal*i; // number (and direction) of positions to shift
					bestscore=shiftedscore[(i-1)];
				}
			} // end of calculating the score for all shifts
			if(bestshift!=0 && (bestshift*dirsignal)>0){ // if we found a shift with a better (or equal) score and it was in this direction
				bestmaxposaffected=maxposaffected;
				i=bestshift*dirsignal; // number of positions to shift (absolute value)
				n=(minnextgaps-i); // number of gaps still remaining at the end of the sequences after the best shift
				for(k=0;k<ntoshift;k++){
					m=postonextgap[k]; // position of the first gap (they were all removed while calculating the scores of all shifts)
					for(j=0;j<n;j++){
						movingsv[m][GAPCODE]++; // re-add all removed valid gaps
						m++;
					}
					bestnposaffected[k]=(postonextgap[k]+i); // number of columns that will be affected
				}
				for(j=0;j<maxposaffected;j++){ // save best working score vector
					if(j<i){ // if it is in the positions at the beginning of the shifted blocks
						for(n=0;n<(ALPHABETSIZE+1);n++)
							bestworkingsv[j][n]=staticsv[j][n];
						bestworkingsv[j][GAPCODE]+=ntoshift; // add the new gaps
						continue;
					}
					for(n=0;n<(ALPHABETSIZE+1);n++) // sum chars in the rest of the columns
						bestworkingsv[j][n]=(staticsv[j][n]+movingsv[(j-i)][n]);
				}
			}
			if(dirsignal==-1) break; // if we have already tested shifting to the left
			looplimit=0; // start of sequences (how far can we go to the left)
			dirsignal=-1; // check the shifts in the reverse direction
		} // end of testing all possible shifts in both directions and finding the one with the best score
		if(bestshift==0) continue; // no shift with a better score was found, so proceed to next candidate columns
		dirsignal=+1;
		if(bestshift<0){
			dirsignal=-1;
			bestshift=-bestshift; // absolute value of shift
		}
		#ifdef DEBUGDP
		fprintf(debugfile,"ORIGINAL:\n");
		PrintStrings(alignedstrings,(col-1),bestmaxposaffected,dirsignal,usableseqs,numseqs,seqstoshift,ntoshift,0);
		PrintCountsVector(scorevector,bestmaxposaffected,col,dirsignal);
		fprintf(debugfile,"SHIFTED:\n");
		PrintStrings(alignedstrings,(col-1),bestmaxposaffected,dirsignal,usableseqs,numseqs,seqstoshift,ntoshift,bestshift);
		PrintCountsVector(bestworkingsv,bestmaxposaffected,0,dirsignal);
		#endif
		for(j=0;j<bestmaxposaffected;j++){ // replace char counts with counts of best shift
			for(n=0;n<(ALPHABETSIZE+1);n++)
				scorevector[(col+dirsignal*j)][n]=bestworkingsv[j][n];
		}
		for(k=0;k<ntoshift;k++){ // move string chars to their shifted positions
			i=seqstoshift[k];
			m=(dirsignal*bestshift); // shift amount and direction
			for(j=(bestnposaffected[k]-1);j>=0;j--){ // process all affected columns from farthest to closest
				n=(col+dirsignal*j); // current column
				if(j<bestshift){ // add gaps in the beginning
					alignedstrings[i][(n-1)]='-';
					continue;
				}
				alignedstrings[i][(n-1)]=alignedstrings[i][(n-m-1)]; // copy char from original column to shifted column
			}
		}
		n=(*consize);
		m=0;
		for(j=col;j<=n;j++){ // get number of consecutive columns to the right of the current column entirely filled with gaps
			if(scorevector[j][GAPCODE]!=numseqs) break;
			m++;
		}
		k=0;
		for(j=(col-1);j>=1;j--){ // get number of consecutive columns to the left of the current column entirely filled with gaps
			if(scorevector[j][GAPCODE]!=numseqs) break;
			k++;
		}
		m=(m+k); // total number of empty columns to be removed
		if(m>auxsvcolsize){
			auxsvcol=(int **)realloc(auxsvcol,m*sizeof(int *));
			auxsvcolsize=m;
		}
		for(j=0;j<m;j++){
			auxsvcol[j]=scorevector[(col-k)+j]; // save empty column pointers to be later placed at the end of the vector
		}
		for(j=(col-k);j<=(n-m);j++){ // process all columns from the left-most empty column to the end of the original vector
			scorevector[j]=scorevector[(j+m)]; // move vector columns to their shifted positions
			for(i=0;i<numseqs;i++){
				ii=usableseqs[i];
				alignedstrings[ii][(j-1)]=alignedstrings[ii][(j+m-1)]; // move the chars in the strings
			}
		}
		for(j=0;j<m;j++){ // place the empty vector columns at the end of the vector
			auxsvcol[j][GAPCODE]=0;
			scorevector[(n-j)]=auxsvcol[j];
			for(i=0;i<numseqs;i++){
				ii=usableseqs[i];
				alignedstrings[ii][(n-j-1)]='\0'; // fill the end positions of the strings with terminator symbols
			}
		}
		(*consize)=((*consize)-m); // update the consensus size of the strings
		col=(col-(k+1)); // if the current column was shifted backwards, check it again and continue from there
		#ifdef DEBUGDP
		fprintf(debugfile,"FINAL:\n");
		PrintStrings(alignedstrings,(col+1-1),bestmaxposaffected,dirsignal,usableseqs,numseqs,seqstoshift,ntoshift,0);
		PrintCountsVector(scorevector,bestmaxposaffected,(col+1),dirsignal);
		#endif
		scoregain+=bestscore;
		deletedcols++;
	} // end of checking all columns for shiftable candidates
	//printf(" %d/%d %d\n",deletedcols,testedcols,scoregain);
	fflush(stdout);
}

// TODO: test different alignment orders (increasing/decreasing length) to see which one gets a better score
// TODO: use realloc instead of repeated free+malloc
// TODO: use: rowgapscore=INDELSCORE*(i)+DOUBLEGAPSCORE*(i-1); ?
// TODO: ?free scorevector variable at end
// TODO: ?test performance with scorevector by linked int lists
void ProgressiveDP(alignmapsegment *segment){
	int i,j,k,l,p,n,m,ncols,nrows;
	int score,charcode;
	int consensussize,prevconsensussize,rowgapscore,colgapscore;
	int **newscorevector;
	int up,left,diag;
	char **strings,**newstrings,*string;
	int pos,startpos,endpos;
	char **dpdirs;
	//int ti,tj,tk,tn,tcharcode,*tmpcount;
	if((segment->maxgapsize)==0) return; // nothing to do if there are no gaps between this segment and the next one
	printf("[(%-4d-%4d)",segment->mingapsize,segment->maxgapsize);fflush(stdout);
	SortSequencesForDP(segment); // sort the sequences by length
	scorevector=NULL;
	strings=NULL;
	dpmatrix=NULL;
	dpdirs=NULL;
	strings=(char **)calloc(numberofseqs,sizeof(char *)); // contains the alignment strings of the (unsorted) sequences
	prevconsensussize=0;
	prevnrows=0;
	consensussize=seqlengths[0]; // initialize variables for first string
	nrows=0;
	ncols=consensussize;
	if(scorevector==NULL){ // ncols x alphabetsize
		scorevector=(int **)calloc((ncols+1),sizeof(int *));
		for(j=1;j<=ncols;j++) scorevector[j]=(int *)calloc((ALPHABETSIZE+1),sizeof(int)); // position 0 of scorevector is not used
	}
	n=orderedseqs[0]; // first ordered sequence
	strings[n]=(char *)malloc((ncols+1)*sizeof(char));
	strings[n][ncols]='\0';
	startpos=(segment->positions[n])+(segment->size);
	pos=startpos;
	string=strings[n];
	for(m=1;m<=ncols;m++){
		string[(m-1)]=CharAt(pos,n); // positions in scorevector start at 1 but on string start at 0
		charcode=CharCodeFromSeq(pos,n);
		scorevector[m][charcode]++;
		pos++;
	}
	for(i=1;i<numberofseqs;i++){ // start at the second ordered sequence
		//ReSortNextSequence(i,consensussize,segment);
		ncols=consensussize;
		nrows=seqlengths[i];
		n=orderedseqs[i];  // current sequence
		if(seqlengths[i]==0){ // if the gap has zero length fill it with the gap symbol
			strings[n]=(char *)malloc((ncols+1)*sizeof(char));
			strings[n][ncols]='\0';
			string=strings[n];
			for(m=(ncols-1);m>=0;m--) string[m]='-';
			continue;
		}
		if(consensussize!=prevconsensussize || nrows>prevnrows){ // if the number of columns or the number of rows is larger than before, realloc the dpmatrix
			if(dpmatrix!=NULL){ // free previous dp matrix
				for(j=0;j<=prevnrows;j++) free(dpmatrix[j]);
				free(dpmatrix);
				dpmatrix=NULL;
			}
			rowgapscore = INDELSCORE*(i); // aligning a char of the current sequence against a gap in each one of the previous i sequences
			dpmatrix=(int **)malloc((nrows+1)*sizeof(int *));
			for(j=0;j<=nrows;j++){ // initialize all rows and all positions of first column
				dpmatrix[j]=(int *)malloc((ncols+1)*sizeof(int));
				dpmatrix[j][0]=j*(rowgapscore); // first column: align char of current seq against i gaps of previous seqs
			}
			colgapscore=0;
			for(j=1;j<=ncols;j++){ // initialize all positions of first row
				colgapscore += DOUBLEGAPSCORE*(scorevector[j][GAPCODE]) + INDELSCORE*(i-(scorevector[j][GAPCODE])); // aligning a gap in the current sequence against all the previously i aligned chars in this column
				dpmatrix[0][j]=(colgapscore); // first row: align gap on current seq against non-gaps of previous i seqs
			}
			if(dpdirs!=NULL){ // free previous dp directions matrix
				for(j=0;j<=prevnrows;j++) free(dpdirs[j]);
				free(dpdirs);
				dpdirs=NULL;
			}
			dpdirs=(char **)malloc((nrows+1)*sizeof(char *));
			for(j=0;j<=nrows;j++){
				dpdirs[j]=(char *)malloc((ncols+1)*sizeof(char));
				dpdirs[j][0]='U'; // first col
			}
			dpdirs[0][0]='D';
			for(j=1;j<=ncols;j++) dpdirs[0][j]='L'; // first row
			prevnrows=nrows;
		}
		startpos=(segment->positions[n])+(segment->size);
		pos=startpos;
		for(j=1;j<=nrows;j++){ // go through all rows
			charcode=CharCodeFromSeq(pos,n);
			for(k=1;k<=ncols;k++){ // go through all columns of each row
				score = MATCHSCORE*(scorevector[k][charcode]) + INDELSCORE*(scorevector[k][GAPCODE]) + MISMATCHSCORE*(i-(scorevector[k][charcode]+scorevector[k][GAPCODE])); // score of aligning the current char in the current seq agains the already aligned char in this column
				rowgapscore = INDELSCORE*(i); // vertical: insert gap in all i previous seq.s, aligned against char of current seq.
				colgapscore = DOUBLEGAPSCORE*(scorevector[k][GAPCODE]) + INDELSCORE*(i-(scorevector[k][GAPCODE])); // horizontal: insert gap in current seq., aligned against all i chars of previous seq.s
				diag=dpmatrix[j-1][k-1]+score; // cell from the upper-left diagonal
				up=dpmatrix[j-1][k]+rowgapscore; // cell from above
				left=dpmatrix[j][k-1]+colgapscore; // cell from the left
				/*
				if(up>diag && up>left){ // prefer insertions in the previous seqs than insertions in the current seq or mismatches
					dpmatrix[j][k]=up;
					dpdirs[j][k]='U';
					continue;
				}
				if(left>diag && left>=up){ // prefer insertions in the current seq than mismatches
					dpmatrix[j][k]=left;
					dpdirs[j][k]='L';
					continue;
				}
				dpmatrix[j][k]=diag;
				dpdirs[j][k]='D';
				*/
				/**/
				if(diag>=up && diag>=left){ // prefer mismatches than insertion or deletions
					dpmatrix[j][k]=diag;
					dpdirs[j][k]='D';
					continue;
				}
				if(left>=up){ // prefer insertions in the current single sequence than in all the previously aligned sequences
					dpmatrix[j][k]=left;
					dpdirs[j][k]='L';
					continue;
				}
				dpmatrix[j][k]=up;
				dpdirs[j][k]='U';
				/**/
			}
			pos++;
		}
		//PrintMatrix(dpmatrix,dpdirs,startpos,n,scorevector,nrows,ncols);
		//printf(" -> ");fflush(stdout);
		// new consensus size
		prevconsensussize=consensussize;
		consensussize=0;
		j=nrows;
		k=ncols;
		while(j>0 && k>0){ // get the new consensus size
			if((dpdirs[j][k])=='D'){
				j--;
				k--;
			} else if((dpdirs[j][k])=='L'){
				k--;
			} else if((dpdirs[j][k])=='U'){
				j--;
			}
			consensussize++;
		}
		if(j>0)	consensussize+=j;
		if(k>0)	consensussize+=k;
		if(consensussize!=prevconsensussize){ // if the consensus size increased, expand the arrays
			newscorevector=(int **)calloc((consensussize+1),sizeof(int *));
			for(j=1;j<=consensussize;j++) newscorevector[j]=(int *)calloc((ALPHABETSIZE+1),sizeof(int));
			newstrings=(char **)calloc(numberofseqs,sizeof(char *));
			for(j=0;j<i;j++){
				p=orderedseqs[j];
				newstrings[p]=(char *)malloc((consensussize+1)*sizeof(char));
				newstrings[p][consensussize]='\0';
			}
		} else { // if the consensus size did not change, keep the previous arrays
			newscorevector=scorevector;
			newstrings=strings;
		}
		newstrings[n]=(char *)malloc((consensussize+1)*sizeof(char)); // string for current seq
		newstrings[n][consensussize]='\0';
		string=newstrings[n];
		j=nrows;
		k=ncols;
		m=(consensussize-1); // current position in the new alignment string
		endpos=(segment->next->positions[n]);
		pos=(endpos-1); // current position in the current sequence
		charcode=CharCodeFromSeq(pos,n); // code of the character at the current position
		while(j>0 && k>0){ // backtrack and fill alignment strings
			if((dpdirs[j][k])=='D'){
				if(consensussize!=prevconsensussize){ // if we have new expanded arrays, copy the contents of the old ones to the new ones
					for(l=0;l<(ALPHABETSIZE+1);l++) newscorevector[(m+1)][l]=scorevector[k][l]; // m starts at 0 , k starts at 1 , scorevectors starts at 1
					for(l=0;l<i;l++){
						p=orderedseqs[l]; // m starts at 0 , k starts at 1 , strings starts at 0
						newstrings[p][m]=strings[p][(k-1)]; // copy from the old column k to the new column m 
					}
				}
				string[m]=CharAt(pos,n);
				newscorevector[(m+1)][charcode]++; // update the scores vector with the newly added char in this column
				pos--;
				charcode=CharCodeFromSeq(pos,n);
				j--;
				k--;
			} else if((dpdirs[j][k])=='L'){
				if(consensussize!=prevconsensussize){
					for(l=0;l<(ALPHABETSIZE+1);l++) newscorevector[(m+1)][l]=scorevector[k][l];
					for(l=0;l<i;l++){
						p=orderedseqs[l];
						newstrings[p][m]=strings[p][(k-1)];
					}
				}
				string[m]='-';
				newscorevector[(m+1)][GAPCODE]++; // one more gap in this column
				k--;
			} else if((dpdirs[j][k])=='U'){
				if(consensussize!=prevconsensussize){
					newscorevector[(m+1)][GAPCODE]=0; // create a new column in the middle of all previously aligned seqs arrays
					for(l=0;l<i;l++){
						p=orderedseqs[l];
						newstrings[p][m]='-';
						newscorevector[(m+1)][GAPCODE]++;
					}
				}
				string[m]=CharAt(pos,n);
				newscorevector[(m+1)][charcode]++;
				pos--;
				charcode=CharCodeFromSeq(pos,n);
				j--;
			}
			m--; // new column behind
		}
		while(j>0){ // if we are not at the first row yet, fill up
			for(l=0;l<i;l++){ // add the missing gaps to the previous seqs
				p=orderedseqs[l];
				newstrings[p][m]='-';
				newscorevector[(m+1)][GAPCODE]++;
			}
			string[m]=CharAt(pos,n);
			newscorevector[(m+1)][charcode]++;
			pos--;
			charcode=CharCodeFromSeq(pos,n);
			j--;
			m--;
		}
		while(k>0){ // if we are not at the first column yet, fill left
			for(l=0;l<(ALPHABETSIZE+1);l++) newscorevector[(m+1)][l]=scorevector[k][l];
			for(l=0;l<i;l++){ // copy the contents of the previous seqs exactly
				p=orderedseqs[l];
				newstrings[p][m]=strings[p][(k-1)];
			}
			string[m]='-'; // add gaps to the current sequence
			newscorevector[(m+1)][GAPCODE]++;
			k--;
			m--;
		}
		if(consensussize!=prevconsensussize){ // free previous arrays if needed
			for(j=1;j<=ncols;j++) free(scorevector[j]);
			free(scorevector);
			scorevector=NULL;
			for(j=0;j<i;j++){
				p=orderedseqs[j];
				free(strings[p]);
			}
			free(strings);
			strings=NULL;
		}
		scorevector=newscorevector;
		strings=newstrings;
		newscorevector=NULL;
		newstrings=NULL;
		if(nrows>prevnrows) prevnrows=nrows; // needs the check because it is re-initialized with the length of the next sequence and it may be lower
		if(ncols>prevncols) prevncols=ncols;
		printf(".");fflush(stdout);
		if(i>1) DeleteGappedColumns(orderedseqs,strings,(i+1),&consensussize,(i+1)/2);
	}
	printf("->%4d]\n",consensussize);
	segment->alignedstrings=strings;
	if(dpmatrix!=NULL){
		for(j=0;j<=nrows;j++) free(dpmatrix[j]);
		free(dpmatrix);
		dpmatrix=NULL;
	}
	if(dpdirs!=NULL){
		for(j=0;j<=nrows;j++) free(dpdirs[j]);
		free(dpdirs);
		dpdirs=NULL;
	}
}

/*
// TODO: reuse and realloc dpmatrix
void RunDynamicProgrammingAlignment(alignmapsegment *segment){
	int i,j,k,n,m,pos,startpos,endpos;
	int maxndprows,maxdprowsize;
	int max,maxpos,aux;
	int charcode,score;
	char *string;
	int ii;
	ndpcols=segment->maxgapsize;
	if(ndpcols==0) return;
	//segment->alignedstrings=(char **)calloc(numberofseqs,sizeof(char *));
	//if((segment->maxgapsize)==(segment->mingapsize)) return;
	if(orderedseqs==NULL) orderedseqs=(int *)calloc(numberofseqs,sizeof(int));
	if(seqlengths==NULL) seqlengths=(int *)calloc(numberofseqs,sizeof(int));
	for(i=0;i<numberofseqs;i++){
		orderedseqs[i]=i;
		seqlengths[i]=(segment->next->positions[i])-((segment->positions[i])+(segment->size));
	}
	for(i=0;i<(numberofseqs-1);i++){
		max=seqlengths[i];
		maxpos=i;
		for(j=i+1;j<numberofseqs;j++){
			if(seqlengths[j]<max){
				max=seqlengths[j];
				maxpos=j;
			}
		}
		if(maxpos!=i){
			aux=orderedseqs[i];
			orderedseqs[i]=orderedseqs[maxpos];
			orderedseqs[maxpos]=aux;
			aux=seqlengths[i];
			seqlengths[i]=seqlengths[maxpos];
			seqlengths[maxpos]=aux;
		}
	}

	ProgressiveDP(segment);
	return;

	if(matchscore==NULL){
		matchscore=(int **)calloc(ndpcols,sizeof(int));
		for(j=0;j<ndpcols;j++) matchscore[j]=(int *)calloc(ALPHABETSIZE,sizeof(int));
	}
	for(i=0;i<numberofseqs;i++){
		if(seqlengths[i]!=ndpcols) break;
		n=orderedseqs[i];
		startpos=(segment->positions[n])+(segment->size);
		pos=startpos;
		for(j=0;j<ndpcols;j++){
			charcode=CharCodeFromSeq(pos,n);
			matchscore[j][charcode]++;
			pos++;
		}
	}
	if((i!=numberofseqs) && (seqlengths[i]!=0)){
		maxndprows=seqlengths[i];
		maxdprowsize=(seqlengths[0]-seqlengths[(numberofseqs-1)]+1);
		if(dpmatrix==NULL){
			dpmatrix=(int **)malloc((maxndprows+1)*sizeof(int));
			for(j=0;j<=maxndprows;j++){
				dpmatrix[j]=(int *)malloc((maxdprowsize+1)*sizeof(int));
				dpmatrix[j][0]=-1;
			}
			for(j=0;j<=maxdprowsize;j++) dpmatrix[0][j]=0;
		}
	}
	firstseq=i;
	//BuildTwoWaysMatchScoreMatrix(segment);
	while(i<numberofseqs){
		ii=i;
		//ii=NextBestSeqForDP(segment);
		n=orderedseqs[ii];
		segment->alignedstrings[n]=(char *)malloc((ndpcols+1)*sizeof(char));
		segment->alignedstrings[n][ndpcols]='\0';
		if(seqlengths[ii]==0){
			string=segment->alignedstrings[n];
			m=(ndpcols-1);
			while(m>=0){
				string[m]='-';
				m--;
			}
			i++;
			continue;
		}
		startpos=(segment->positions[n])+(segment->size);
		endpos=(segment->next->positions[n]);
		pos=startpos;
		ndprows=seqlengths[ii];
		dprowsize=(ndpcols-ndprows+1);
		for(j=1;j<=ndprows;j++){
			charcode=CharCodeFromSeq(pos,n);
			for(k=1;k<=dprowsize;k++){
				score=matchscore[(j+k-2)][charcode];
				score+=dpmatrix[j-1][k];
				if((dpmatrix[j][k-1])>score) dpmatrix[j][k]=dpmatrix[j][k-1];
				else dpmatrix[j][k]=score;
			}
			pos++;
		}
		j=ndprows;
		k=dprowsize;
		m=(ndpcols-1);
		pos=(endpos-1);
		string=segment->alignedstrings[n];
		while(m>=0){
			if((dpmatrix[j][k])==(dpmatrix[j][k-1])){
				string[m]='-';
				k--;
			} else {
				string[m]=CharAt(pos,n);
				charcode=CharCodeFromSeq(pos,n);
				matchscore[m][charcode]++;
				pos--;
				j--;
			}
			m--;
		}
		i++;
	}
	if(dpmatrix!=NULL){
		for(i=0;i<=maxndprows;i++) free(dpmatrix[i]);
		free(dpmatrix);
		dpmatrix=NULL;
	}
	if(matchscore!=NULL){
		for(i=0;i<ndpcols;i++) free(matchscore[i]);
		free(matchscore);
		matchscore=NULL;
	}
	if(scores!=NULL){
		free(scores);
		scores=NULL;
	}
}
*/

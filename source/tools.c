#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/timeb.h>

#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif

// Cleans the FASTA file by removing all the gaps, newlines and non-DNA characters
void CleanDNAFastaFile(char *fastafilename){
	struct timeb starttb,endtb;
	double timetb;
	char c;
	FILE *inputfile,*outputfile;
	char *outputfilename;
	int numberofseqs,desclen;
	long int seqlen,nvalidchars,ninvalidchars,nspecialchars,nextrachars,seqstart,seqend,seqsize;
	printf("> Loading sequences from file <%s> ... ",fastafilename);
	if((inputfile=fopen(fastafilename,"r"))==NULL) {
		printf("\n> ERROR: Sequence file not found\n");
		return;
	} else { // get file size
		seqstart=ftell(inputfile);
		fseek(inputfile,0L,SEEK_END);
		seqend=ftell(inputfile);
		seqsize=(seqend-seqstart);
		rewind(inputfile);
		printf("(%ld bytes)\n",seqsize);
		fflush(stdout);
	}
	c=fgetc(inputfile); // check if it is a valid fasta file
	if(c!='>') {
		printf("> ERROR: Invalid FASTA file\n");
		return;
	}
	outputfilename=(char *)calloc(((int)strlen(fastafilename)+6+1),sizeof(char));
	strcpy(outputfilename,"Clean-"); // append this suffix to the file name
	strcat(outputfilename,fastafilename);
	outputfile=fopen(outputfilename,"w");
	if(outputfile==NULL) {
		printf("> ERROR: Can't write output file\n");
		return;
	}
	numberofseqs=0;
	ftime(&starttb);
	while(1){
		fputc('>',outputfile);
		printf("  [");
		desclen=0;
		while((c=fgetc(inputfile))!=EOF && c!='\n' && c!='\r') { // get sequence description
			fputc(c,outputfile);
			if(desclen<20) putchar(c);
			desclen++;
		}
		fputc('\n',outputfile);
		while(desclen<20) {
			putchar(' ');
			desclen++;
		}
		printf("] ");
		fflush(stdout);
		seqlen=0;
		nvalidchars=0;
		ninvalidchars=0;
		nspecialchars=0;
		nextrachars=0;
		while((c=fgetc(inputfile))!=EOF && c!='>') { // get new chars until the end of the file or the beginning of a new sequence
			seqlen++;
			if(c=='A' || c=='C' || c=='G' || c=='T') { // basic DNA letters
				fputc(c,outputfile);
				nvalidchars++;
				continue;
			}
			if(c=='a' || c=='c' || c=='g' || c=='t') { // low case letters
				fputc((char)(c-32),outputfile);
				nvalidchars++;
				continue;
			}
			if(c=='\n' || c=='\r' || c=='\0' || c==' ') { // blank characters
				nspecialchars++;
				continue;
			}
			if(c=='R' || c=='Y' || c=='S' || c=='W' || c=='K' || c=='M' ||
					c=='D' || c=='H' || c=='B' || c=='V' || c=='N') { // extended IUPAC characters
				nextrachars++;
				continue;
			}
			if(c=='r' || c=='y' || c=='s' || c=='w' || c=='k' || c=='m' ||
					c=='d' || c=='h' || c=='b' || c=='v' || c=='n') { // low case characters
				nextrachars++;
				continue;
			}
			ninvalidchars++; // invalid character
		}
		fputc('\n',outputfile);
		printf("(%ld chars: %ldV %ldS %ldX %ldI)\n",seqlen,nvalidchars,nspecialchars,nextrachars,ninvalidchars);
		fflush(stdout);
		numberofseqs++;
		if(c==EOF) break;
	}
	ftime(&endtb);
	timetb=((endtb.time) + (endtb.millitm)/1000.0) - ((starttb.time) + (starttb.millitm)/1000.0);
	printf("> %d sequence(s) processed in %.3lf seconds\n",numberofseqs,timetb);
	printf("> Saving sequences to file <%s> ... ",outputfilename);
	fflush(stdout);
	seqend=ftell(outputfile);
	rewind(outputfile);
	seqstart=ftell(outputfile);
	seqsize=(seqend-seqstart);
	printf("(%ld bytes)\n",seqsize);
	fclose(inputfile);
	free(outputfilename);
	if(fclose(outputfile)==EOF) {
		printf("> ERROR: Can't save output file\n");
		return;
	}
	//printf("> Done!\n");
}

// Checks if the final aligned strings match their original sequences
void TestAlignmentFileOutput(char *filename1, char *filename2){
	FILE *file1,*file2;
	char *desc1,*desc2;
	char c1,c2;
	int i1,i2,n1,n2,ns1,ns2;
	int retval=0;
	retval++;
	if((file1=fopen(filename1,"r"))==NULL) return;
	if((file2=fopen(filename2,"r"))==NULL) return;
	printf("> Checking integrity of aligned sequences... ");
	fflush(stdout);
	desc1=(char *)calloc((64+1),sizeof(char));
	desc2=(char *)calloc((64+1),sizeof(char));
	ns1=0; // number of sequences
	ns2=0;
	i1=0; // number of chars
	i2=0;
	n1=0; // number of valid chars
	n2=0;
	c1=fgetc(file1);
	c2=fgetc(file2);
	while(c1!=EOF || c2!=EOF){ // while at least one of the sequences has chars left
		if(c1=='>'){ // new seq
			//if(ns1>0) printf("[<%s>:%d:\"%s\":(%d/%d)]\n",filename1,ns1,desc1,n1,(i1-1)); // print info of previous seq
			retval=fscanf(file1,"%64[^\n]s",desc1); // get seq description
			ns1++; // one more seq
			i1=0; // reset char counters
			n1=0;
		}
		if(c2=='>'){ // new seq
			//if(ns2>0) printf("[<%s>:%d:\"%s\":(%d/%d)]\n",filename2,ns2,desc2,n2,(i2-1)); // print info of previous seq
			retval=fscanf(file2,"%64[^\n]s",desc2); // get seq description
			ns2++; // one more seq
			i2=0; // reset char counters
			n2=0;
		}
		if(c1!=c2){ // report error if chars in both sequences do not match
			//printf("> ERROR: Mismatch found at: '%c'@[%d:%d]=!='%c'@[%d:%d]\n",c1,ns1,i1,c2,ns2,i2);
			printf("ERROR at: '%c'@[%d:%d]=!='%c'@[%d:%d]\n",c1,ns1,i1,c2,ns2,i2);
			break;
		} else { // one more valid non gap char
			if(c1!='>') n1++;
			if(c2!='>') n2++;
		}
		c1=fgetc(file1); // get next char
		if(c1!='\n') i1++;
		while(c1=='-' || c1=='\n'){ // skip gap chars and newlines
			c1=fgetc(file1);
			if(c1!='\n') i1++;
		}
		c2=fgetc(file2); // get next char
		if(c2!='\n') i2++;
		while(c2=='-' || c2=='\n'){ // skip gap chars and newlines
			c2=fgetc(file2);
			if(c2!='\n') i2++;
		}
	}
	if(c1==EOF && c2==EOF){ // if the loop ended correctly
		//printf("[<%s>:%d:\"%s\":(%d/%d)]\n",filename1,ns1,desc1,n1,(i1-1)); // print info for last seq
		//printf("[<%s>:%d:\"%s\":(%d/%d)]\n",filename2,ns2,desc2,n2,(i2-1));
		//printf("> Aligned sequences match their original sequences!\n");
		printf("OK\n");
		fflush(stdout);
	}
	fclose(file1);
	fclose(file2);
	free(desc1);
	free(desc2);
}

// Calculate the alignment score (sum of pairs score)
int CalculateSumOfPairsScore(char *fastafilename){
	struct timeb starttb,endtb;
	double timetb;
	char c;
	FILE *inputfile;
	fpos_t *positions;
	char *description,*chars;
	int *sizes;
	int i,j,n,numberofseqs,size,ngaps,nconscols,score;
	long int filestart,fileend,filesize;
	char *retval=NULL;
	retval++;
	printf("> Opening file <%s> ... ",fastafilename);
	if((inputfile=fopen(fastafilename,"r"))==NULL) {
		printf("\n> ERROR: Sequence file not found\n");
		return -1;
	} else { // get file size
		filestart=ftell(inputfile);
		fseek(inputfile,0L,SEEK_END);
		fileend=ftell(inputfile);
		filesize=(fileend-filestart);
		rewind(inputfile);
		printf("(%ld bytes)\n",filesize);
		fflush(stdout);
	}
	c=fgetc(inputfile);
	if(c==EOF || c!='>') {
		printf("> ERROR: No sequences or no descriptions in file\n");
		return -1;
	}
	description=(char *)calloc(255,sizeof(char));
	numberofseqs=0;
	while(c!=EOF){ // count number of sequences
		if(c=='>'){
			numberofseqs++;
			while(c!='\n' && c!=EOF) c=fgetc(inputfile); // skip rest of description
		}
		c=fgetc(inputfile);
	}
	if(numberofseqs<2) {
		printf("> ERROR: Not enough sequences in file\n");
		return -1;
	}
	rewind(inputfile);
	sizes=(int *)malloc(numberofseqs*sizeof(int)); // length of each sequence
	positions=(fpos_t *)malloc(numberofseqs*sizeof(fpos_t)); // starting position of each sequence inside the file
	for(i=0;i<numberofseqs;i++){
		retval=fgets(description,255,inputfile); // get description
		fgetpos(inputfile,&positions[i]); // save start position
		size=0;
		while((c=fgetc(inputfile))!='>' && c!='=' && c!=EOF) if(c!='\r' && c!='\n') size++; // count valid chars
		if(c=='=') while(c!='>' && c!=EOF) c=fgetc(inputfile); // skip extra info line in case of XMFA format
		sizes[i]=size;
	}
	size=sizes[0];
	for(i=0;i<numberofseqs;i++){ // check if all sequences have the same size
		if(sizes[i]!=size){
			printf("> ERROR: Consensus sizes are not consistent\n");
			return -1;
		}
	}
	chars=(char *)malloc(numberofseqs*sizeof(char)); // stores an alignment column, i.e., the chars at the same position on each sequence
	ngaps=0;
	nconscols=0;
	score=0;
	ftime(&starttb);
	for(n=0;n<size;n++){ // process all positions
		for(i=0;i<numberofseqs;i++){ // get the char at this position in each sequence
			fsetpos(inputfile,&positions[i]); // restore saved position for this sequence
			c=fgetc(inputfile); // get next valid char
			while(c=='\r' || c=='\n') c=fgetc(inputfile); // skip newlines
			fgetpos(inputfile,&positions[i]); // saved current position
			chars[i]=c;
			if(c=='-') ngaps++; // count total number of gaps in the alignment
		}
		c=chars[0];
		for(i=1;i<numberofseqs;i++){
			if(c!=chars[i]) break;
		}
		if(i==numberofseqs) nconscols++; // check if the column is totally conserved (the same char is on all sequences)
		for(i=0;i<=(numberofseqs-2);i++){ // for this position, compare the char on each sequence with that chars on all following sequences
			for(j=(i+1);j<=(numberofseqs-1);j++){
				if(chars[i]=='-' && chars[j]=='-') continue; // a gap aligned against another gap does not have any score
				if(chars[i]==chars[j]) score++; // if the char is the same
				else score--; // if the chars are different
			}
		}
	}
	ftime(&endtb);
	timetb=((endtb.time) + (endtb.millitm)/1000.0) - ((starttb.time) + (starttb.millitm)/1000.0);
	printf("> %d sequence(s) processed in %.3lf seconds\n",numberofseqs,timetb);
	printf("> Statistics:\nConsensus size = %d\nAverage gaps per sequence = %d\nNumber of conserved columns = %d\nSum-of-Pairs score = %d\n",size,(ngaps/numberofseqs),nconscols,score);
	fclose(inputfile);
	free(positions);
	free(sizes);
	free(description);
	free(chars);
	//printf("> Done!\n");
	return score;
}

/*
void CompareToReferenceAlignment(char *testfilename, char *referencefilename){
	struct timeb starttb,endtb;
	double timetb;
	char tc, rc;
	FILE *testfile,*reffile;
	fpos_t *tpositions, *rpositions;
	char *tchars,*rchars;
	int *sizes;
	int i,j,n,numberofseqs,size,ngaps,nconscols,score;
	long int filestart,fileend,filesize;
	printf("> Opening test alignment file <%s> ... ",testfilename); // open and check test file
	if((testfile=fopen(testfilename,"r"))==NULL) {
		printf("\n> ERROR: Sequence file not found\n");
		return;
	} else { // get file size
		filestart=ftell(testfile);
		fseek(testfile,0L,SEEK_END);
		fileend=ftell(testfile);
		filesize=(fileend-filestart);
		rewind(testfile);
		printf("(%ld bytes)\n",filesize);
		fflush(stdout);
	}
	tc=fgetc(testfile);
	if(tc==EOF || tc!='>') {
		printf("> ERROR: No sequences or no descriptions in file\n");
		return;
	}
	numberofseqs=0;
	while(tc!=EOF){ // count number of sequences
		if(tc=='>'){
			numberofseqs++;
			while(tc!='\n' && tc!=EOF) tc=fgetc(testfile); // skip rest of description
		}
		tc=fgetc(testfile);
	}
	if(numberofseqs<2) {
		printf("> ERROR: Not enough sequences in file\n");
		return;
	}
	n=numberofseqs; // store number of sequences in test file to later compare with ref file
	printf("> Opening reference alignment file <%s> ... ",referencefilename); // open and check reference file
	if((reffile=fopen(referencefilename,"r"))==NULL) {
		printf("\n> ERROR: Sequence file not found\n");
		return;
	} else { // get file size
		filestart=ftell(reffile);
		fseek(reffile,0L,SEEK_END);
		fileend=ftell(reffile);
		filesize=(fileend-filestart);
		rewind(reffile);
		printf("(%ld bytes)\n",filesize);
		fflush(stdout);
	}
	rc=fgetc(reffile);
	if(rc==EOF || rc!='>') {
		printf("> ERROR: No sequences or no descriptions in file\n");
		return;
	}
	numberofseqs=0;
	while(rc!=EOF){ // count number of sequences
		if(rc=='>'){
			numberofseqs++;
			while(rc!='\n' && rc!=EOF) rc=fgetc(reffile); // skip rest of description
		}
		rc=fgetc(reffile);
	}
	if(numberofseqs<2) {
		printf("> ERROR: Not enough sequences in file\n");
		return;
	}
	if(numberofseqs!=n){ // compare number of sequences in both files
		printf("> ERROR: Number of sequences do not match in both files\n");
		return;
	}
	rewind(testfile);
	description=(char *)calloc(255,sizeof(char));
	sizes=(int *)malloc(numberofseqs*sizeof(int)); // length of each sequence
	positions=(fpos_t *)malloc(numberofseqs*sizeof(fpos_t)); // starting position of each sequence inside the file
	for(i=0;i<numberofseqs;i++){
		fgets(description,255,inputfile); // get description
		fgetpos(inputfile,&positions[i]); // save start position
		size=0;
		while((c=fgetc(inputfile))!='>' && c!=EOF) if(c!='\r' && c!='\n') size++; // count valid chars
		sizes[i]=size;
	}
	size=sizes[0];
	for(i=0;i<numberofseqs;i++){ // check if all sequences have the same size
		if(sizes[i]!=size){
			printf("> ERROR: Consensus sizes are not consistent\n");
			return -1;
		}
	}
	chars=(char *)malloc(numberofseqs*sizeof(char)); // stores an alignment column, i.e., the chars at the same position on each sequence
	ngaps=0;
	nconscols=0;
	score=0;
	ftime(&starttb);
	for(n=0;n<size;n++){ // process all positions
		for(i=0;i<numberofseqs;i++){ // get the char at this position in each sequence
			fsetpos(inputfile,&positions[i]); // restore saved position for this sequence
			c=fgetc(inputfile); // get next valid char
			while(c=='\r' || c=='\n') c=fgetc(inputfile); // skip newlines
			fgetpos(inputfile,&positions[i]); // saved current position
			chars[i]=c;
			if(c=='-') ngaps++; // count total number of gaps in the alignment
		}
		c=chars[0];
		for(i=1;i<numberofseqs;i++){
			if(c!=chars[i]) break;
		}
		if(i==numberofseqs) nconscols++; // check if the column is totally conserved (the same char is on all sequences)
		for(i=0;i<=(numberofseqs-2);i++){ // for this position, compare the char on each sequence with that chars on all following sequences
			for(j=(i+1);j<=(numberofseqs-1);j++){
				if(chars[i]=='-' && chars[j]=='-') continue; // a gap aligned against another gap does not have any score
				if(chars[i]==chars[j]) score++; // if the char is the same
				else score--; // if the chars are different
			}
		}
	}
	ftime(&endtb);
	timetb=((endtb.time) + (endtb.millitm)/1000.0) - ((starttb.time) + (starttb.millitm)/1000.0);
	printf("> %d sequence(s) processed in %.3lf seconds\n",numberofseqs,timetb);
	printf("> Statistics:\nConsensus size = %d\nAverage gaps per sequence = %d\nNumber of conserved columns = %d\nSum-of-Pairs score = %d\n",size,(ngaps/numberofseqs),nconscols,score);
	fclose(inputfile);
	free(positions);
	free(sizes);
	free(description);
	free(chars);
	//printf("> Done!\n");
	return score;
}
*/

// Converts a file in the FASTA format to the MSF format
void ConvertFastaToMsf(char *fastafilename){
	FILE *fastafile,*msffile;
	char *msffilename;
	char **descriptions;
	fpos_t *currentpos;
	char c;
	int numseqs,*seqsizes;
	int s,i,j,m,n;
	int alignmentsize;
	printf("> Opening FASTA file <%s>... ",fastafilename);
	fflush(stdout);
	n=(int)strlen(fastafilename);
	for(i=(n-1);i>0;i--){ // get base name
		if(fastafilename[i]=='.') break;
	}
	if(i==0) i=n;
	msffilename=(char *)calloc((i+4+1),sizeof(char));
	strncpy(msffilename,fastafilename,i);
	msffilename[i]='\0';
	strcat(msffilename,".msf"); // create new filename
	if((fastafile=fopen(fastafilename,"r"))==NULL){
		printf("> ERROR: Failed to open FASTA file\n");
		return;
	}
	if((msffile=fopen(msffilename,"w"))==NULL){
		printf("> ERROR: Failed to create MSF file\n");
		return;
	}
	numseqs=0;
	c='\n';
	while(c!=EOF){ // count number of sequences
		c=fgetc(fastafile);
		if(c=='>'){
			numseqs++;
			while(c!='\n' && c!=EOF) c=fgetc(fastafile);
		}
	}
	if(numseqs==0){
		printf("> ERROR: No sequences found in FASTA file\n");
		return;
	}
	descriptions=(char **)malloc(numseqs*sizeof(char *));
	seqsizes=(int *)malloc(numseqs*sizeof(int));
	currentpos=(fpos_t *)malloc(numseqs*sizeof(fpos_t));
	rewind(fastafile);
	for(s=0;s<numseqs;s++){
		while(c!='>') c=fgetc(fastafile); // advance to beggining of sequence if not already at it
		descriptions[s]=(char *)malloc(11*sizeof(char));
		i=0;
		c=fgetc(fastafile);
		while(i<10 && c!='\n'){ // get description with no spaces and a maximum of 10 chars
			if(c!=' '){
				descriptions[s][i]=c;
				i++;
			}
			c=fgetc(fastafile);
		}
		descriptions[s][i]='\0';
		while(c!='\n') c=fgetc(fastafile); // skip rest of description
		seqsizes[s]=0;
		fgetpos(fastafile,&(currentpos[s])); // save start position of alignment chars
		c=fgetc(fastafile);
		while(c!='>' && c!='=' && c!=EOF){ // get number of valid chars (only alphabet letters and gaps)
			if( (c>='A' && c<='Z') || (c>='a' && c<='z') || c=='-' ) seqsizes[s]++;
			c=fgetc(fastafile);
		} // stop at new sequence, at end of file or at extra info line in case of XMFA format
		if(s!=0 && seqsizes[s]!=seqsizes[(s-1)]){
			printf("> ERROR: Sequences alignment sizes do not match\n");
			return;
		}
	}
	alignmentsize=seqsizes[0];
	printf("(%d aligned sequences of size %d)\n",numseqs,alignmentsize);
	printf("> Saving alignments to MSF file <%s>... ",msffilename);
	fflush(stdout);
	fprintf(msffile,"!!NA_MULTIPLE_ALIGNMENT 1.0\n\n"); // print MSF header
	fprintf(msffile," %s \tMSF: %d \tType: N \tCheck: 0 \t..\n\n",msffilename,alignmentsize);
	for(s=0;s<numseqs;s++){
		fprintf(msffile," Name: %s oo\tLen: %d \tCheck: 0 \tWeight: 1.00 \n",descriptions[s],seqsizes[s]);
	}
	fprintf(msffile,"\n//\n\n");
	n=0;
	while(n<alignmentsize){ // print sequences
		for(s=0;s<numseqs;s++){ // all sequences
			m=n;
			fsetpos(fastafile,&(currentpos[s])); // set current position of sequence
			fprintf(msffile,"%s \t",descriptions[s]);
			for(i=0;i<5;i++){ // 5 columns
				for(j=0;j<10;j++){ // of 10 chars each
					while(1){ // get next valid char from sequence
						c=fgetc(fastafile);
						if((c>='A' && c<='Z') || (c>='a' && c<='z')) break;
						if(c=='-'){ // convert dashes to dots
							c='.';
							break;
						}
					}
					fprintf(msffile,"%c",c); // write char
					m++;
					if(m>=alignmentsize) break;
				}
				if(m>=alignmentsize) break;
				fprintf(msffile," "); // separate columns by space
			}
			fgetpos(fastafile,&(currentpos[s])); // save current position of sequence
			fprintf(msffile,"\n"); // next sequence after 5 columns
		}
		n+=50;
		fprintf(msffile,"\n"); // double new line after all 5 columns of all sequences
	}
	fclose(fastafile);
	if(fclose(msffile)==EOF){
		printf("> ERROR: Failed to save MSF file\n");
		return;
	}
	printf("OK\n");
	fflush(stdout);
	free(msffilename);
	for(s=0;s<numseqs;s++) free(descriptions[s]);
	free(descriptions);
	free(currentpos);
	free(seqsizes);
}

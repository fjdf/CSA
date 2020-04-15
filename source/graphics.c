#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "bitmap.h"
#include "graphics.h"
#include "csamsa.h"
#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif

#define ALPHABETSIZE 32

// *****************
// VARIAVEIS GLOBAIS
// *****************

static int **alphabet = NULL; // alfabeto
static int alphabetsize = 96; // tamanho do alfabeto
//static int charheight = 6; // tamanho dos caracteres
//static int charwidth = 5; // tamanho dos caracteres
static uint8_t black,white,red,green,blue,yellow,purple,cyan,grey; // cores
static int bmpwidth = 320; // largura da imagem
static int bmpheight = 240; // altura da imagem

static int numberofsequences = 0;
//static int numberofboxes = 0;
static int boxstartx = 0;
static int *boxstarty = NULL;
static int boxheight = 0;
static double pixelsperunit = 0.0;
static int *boxmiddlepoint = NULL;
static int numberofcolors = 0;
static uint8_t currentcolor = 0;
static double *numcolorspersizeunit = NULL;
//static uint8_t *boxcolor = NULL;
static int *sequencesize = NULL;
static int *startposition = NULL;
static int headersize = 0;
static int lowheadersize = 0;
static int interseqsize = 0;
static FILE *mapfile = NULL;

// ******************
// FUN�ES AUXILIARES
// ******************

// m�imo
int xymax(int x, int y){
	//return ( (x>y) ? x : y );
	if(x>y) return x;
	return y;
}

// m�imo
int xymin(int x, int y){
	//return ( (x<y) ? x : y );
	if(x<y) return x;
	return y;
}

// m�ulo positivo
// ex.: posmod(-12,10)=2
int posmod(int x, int y){
	int modres=x%y;
	if(modres<0) modres+=y;
	return modres;
}

// m�ulo negativo
// ex.: posmod(-12,10)=8
int negmod(int x, int y){
	int modres=x%y;
	if(modres>0) modres-=y;
	if(modres<0) modres=-modres;
	return modres;
}

// dist�cia entre dois pontos
int distance(int x, int y){
	return (int)floor(sqrt(pow((double)x,2.0)+pow((double)y,2.0)));
}

// troca o contedo de duas vari�eis
void switchint(int *i1, int *i2){
	int aux;
	aux=(*i1);
	(*i1)=(*i2);
	(*i2)=aux;
}

// compara dois valores inteiros (para usar no 'qsort')
int comp(void *a, void *b) {
	return ( ( (*(int *)a) > (*(int *)b) ) ? 1 : 0 );
}

// converte um double para um int (arredonda para baixo)
int double2int(double d){
	return (int)floor(d);
}

// converte um double para um int (arredonda para cima)
int double2intup(double d){
	//if(d<0) return -(int)floor(-d);
	return (int)ceil(d);
}

// converte uma lista de inteiros para uma lista de doubles
double * intlist2doublelist(int *intlist, int n){
	int i=0;
	double *doublelist;
	doublelist=(double *)malloc(n*sizeof(double));
	for(i=0;i<n;i++) doublelist[i]=((double)(intlist[i]));
	return doublelist;
}


// *********************
// PRIMITIVAS DE DESENHO
// *********************

// inicializa as vari�eis das cores
void initializeColors(){
	black=getColorFromPalette(0,0,0);
	white=getColorFromPalette(255,255,255);
	red=getColorFromPalette(255,0,0);
	green=getColorFromPalette(0,255,0);
	blue=getColorFromPalette(0,0,255);
	yellow=getColorFromPalette(255,255,0);
	purple=getColorFromPalette(255,0,255);
	cyan=getColorFromPalette(0,255,255);
	grey=getColorFromPalette(128,128,128);
}

// desenha uma linha horizontal (da esquerda para a direita)
int drawHLine(int x, int y, int size, uint8_t color){
	int i,pos,limit;
	pos=x;
	limit=bmpwidth;
	for(i=1;(i<=size && pos<limit);i++)
		drawPoint(pos++,y,color);	
	return i;
}

// desenha uma linha vertical (de cima para baixo / Y menor para Y maior)
int drawVLine(int x, int y, int size, uint8_t color){
	int i,pos,limit;
	pos=y;
	limit=bmpheight;
	for(i=1;(i<=size && pos<limit);i++)
		drawPoint(x,pos++,color);
	return i;
}

// desenha uma linha diagonal (da esquerda para a dreita)
int drawDLine(int x, int y, double m, int size, uint8_t color){
	int i,xpos,ypos,xlimit,ylimit,prevypos,step;
	xpos=x;
	ypos=y;
	xlimit=bmpwidth;
	ylimit=bmpheight;
	prevypos=y;
	if(m>0) step=1;
	else step=-1;
	for(i=0;(i<=size && xpos<=xlimit && ypos<=ylimit);i++){
		ypos=(int)floor(m*(double)i)+y;
		drawPoint(xpos,ypos,color);
		while(prevypos!=ypos){
			drawPoint(xpos-1,prevypos,color);
			prevypos=prevypos+step;
		}
		prevypos=ypos;
		xpos++;
	}
	return i;
}

// TODO: move code to bitmap.c
// FALTA: alterar VLine e HLine e DLine para suportarem recta de sentido contr�io
// desenha uma linha (geral)
int drawLine(int x1, int y1, int x2, int y2, uint8_t color){
	int xdif,ydif,size;
	int xx1,xx2,yy1,yy2,xswitch,yswitch;
	double m,signal;
	xx1=x1;
	xx2=x2;
	yy1=y1;
	yy2=y2;
	xswitch=0;
	yswitch=0;
	if(x2<x1){
		xx2=x1;
		xx1=x2;
		xswitch=1;
	}
	if(y2<y1){
		yy2=y1;
		yy1=y2;
		yswitch=1;
	}
	xdif=xx2-xx1;
	ydif=yy2-yy1;
	if(xdif==0) return drawVLine(xx1,yy1,ydif,color);
	if(ydif==0) return drawHLine(xx1,yy1,xdif,color);
	if(xswitch==yswitch) signal=1.0;
	else signal=-1.0;
	if(xswitch==1) yy1=y2;
	else yy1=y1;
	m=((double)ydif/(double)xdif)*signal;
	size=xdif;
	return drawDLine(xx1,yy1,m,size,color);
}

// TODO: define border length
// desenha um rect�gulo vazio
void drawRectangle(int x, int y, int w, int h, uint8_t color){
	int i,j;
	for(j=0;j<h;j++){
		if(j==0 || j==(h-1)) for(i=0;i<w;i++) drawPoint(x+i,y+j,color);
		else { drawPoint(x,y+j,color); drawPoint(x+(w-1),y+j,color); }
	}
}

// desenha um rect�gulo preenchido
void drawFilledRectangle(int x, int y, int w, int h, uint8_t color){
	int i,j;
	for(i=0;i<w;i++)
		for(j=0;j<h;j++)
			drawPoint(x+i,y+j,color);
}

// desenha um c�culo
void drawCircle(int x, int y, int r, uint8_t color){
	int xx,yy,i;
	for(xx=1;xx<r;xx++){
		yy=(int)floor(sqrt(pow((double)r,2.0)-pow((double)xx,2.0)));
		//drawPoint(x+xx,y+yy,color);
		//drawPoint(x+xx,y-yy,color);
		//drawPoint(x-xx,y+yy,color);
		//drawPoint(x-xx,y-yy,color);
		for(i=-xx;i<=xx;i++) drawPoint(x+i,y+yy,color);
		for(i=-xx;i<=xx;i++) drawPoint(x+i,y-yy,color);
		for(i=-yy;i<=yy;i++) drawPoint(x+xx,y+i,color);
		for(i=-yy;i<=yy;i++) drawPoint(x-xx,y+i,color);
	}
	drawVLine(x,y-r+1,2*r-1,color);
}

// desenha uma estrela
void drawStar(int x, int y, int r, uint8_t color){
	int i;
	for(i=0;i<r;i++){
		drawPoint(x+i,y+i,color);
		drawPoint(x-i,y-i,color);
		drawPoint(x-i,y+i,color);
		drawPoint(x+i,y-i,color);
		drawPoint(x+i,y,color);
		drawPoint(x-i,y,color);
		drawPoint(x,y+i,color);
		drawPoint(x,y-i,color);
	}
}


// ********************
//  PRIMITIVAS DE TEXTO
// ********************

// TODO: convert and store alphabet characters as integers
// cria um caracter de 6x5 pixels
int *newChar(int i1,int i2,int i3,int i4,int i5,int i6,int i7,int i8,int i9,int i10,int i11,int i12,int i13,int i14,int i15,int i16,int i17,int i18,int i19,int i20,int i21,int i22,int i23,int i24,int i25,int i26,int i27,int i28,int i29,int i30){
	int *digit=(int *)malloc(30*sizeof(int));
	digit[0]=i1;
	digit[1]=i2;
	digit[2]=i3;
	digit[3]=i4;
	digit[4]=i5;
	digit[5]=i6;
	digit[6]=i7;
	digit[7]=i8;
	digit[8]=i9;
	digit[9]=i10;
	digit[10]=i11;
	digit[11]=i12;
	digit[12]=i13;
	digit[13]=i14;
	digit[14]=i15;
	digit[15]=i16;
	digit[16]=i17;
	digit[17]=i18;
	digit[18]=i19;
	digit[19]=i20;
	digit[20]=i21;
	digit[21]=i22;
	digit[22]=i23;
	digit[23]=i24;
	digit[24]=i25;
	digit[25]=i26;
	digit[26]=i27;
	digit[27]=i28;
	digit[28]=i29;
	digit[29]=i30;
	return digit;
}

// inicializa o conjunto dos 96 caracteres a usar:
void initializeAlphabet(){
	int **digits;
	if(alphabet!=NULL) return;
	digits=(int **)malloc(96*sizeof(int *));
	
	digits[0]=newChar // ' '
	(0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0);
	digits[1]=newChar // '!'
	(0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,0,0,0
	,0,0,1,0,0);
	digits[2]=newChar // '"'
	(0,1,0,1,0
	,0,1,0,1,0
	,0,1,0,1,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0);
	digits[3]=newChar // '#'
	(0,1,0,1,0
	,0,1,0,1,0
	,1,1,1,1,1
	,0,1,0,1,0
	,1,1,1,1,1
	,0,1,0,1,0);
	digits[4]=newChar // '$'
	(0,0,1,0,0
	,0,1,1,1,1
	,1,0,1,0,0
	,0,1,1,1,0
	,0,0,1,0,1
	,1,1,1,1,0);
	digits[5]=newChar // '%'
	(1,1,0,0,0
	,1,1,0,0,1
	,0,0,0,1,0
	,0,0,1,0,0
	,0,1,0,1,1
	,1,0,0,1,1);
	digits[6]=newChar // '&'
	(0,1,0,0,0
	,1,0,1,0,0
	,0,1,0,0,0
	,0,1,1,0,1
	,1,0,0,1,0
	,0,1,1,0,1);
	digits[7]=newChar // '''
	(0,0,1,0,0
	,0,0,1,0,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0);
	digits[8]=newChar // '('
	(0,0,1,1,0
	,0,1,0,0,0
	,0,1,0,0,0
	,0,1,0,0,0
	,0,1,0,0,0
	,0,0,1,1,0);
	digits[9]=newChar // ')'
	(0,1,1,0,0
	,0,0,0,1,0
	,0,0,0,1,0
	,0,0,0,1,0
	,0,0,0,1,0
	,0,1,1,0,0);
	digits[10]=newChar // '*'
	(0,0,0,0,0
	,0,0,1,0,0
	,1,0,1,0,1
	,0,1,1,1,0
	,0,1,1,1,0
	,1,0,1,0,1);
	digits[11]=newChar // '+'
	(0,0,0,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,1,1,1,1,1
	,0,0,1,0,0
	,0,0,1,0,0);	
	digits[12]=newChar // ','
	(0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,0,1,0,0
	,0,1,0,0,0);
	digits[13]=newChar // '-'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,1,1,1,0
	,0,0,0,0,0
	,0,0,0,0,0);
	digits[14]=newChar // '.'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,1,1,0,0
	,0,1,1,0,0);
	digits[15]=newChar // '/'
	(0,0,0,0,1
	,0,0,0,1,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,1,0,0,0
	,1,0,0,0,0);
	digits[16]=newChar // '0'
	(0,1,1,1,0
	,1,0,0,1,1
	,1,0,1,0,1
	,1,0,1,0,1
	,1,1,0,0,1
	,0,1,1,1,0);
	digits[17]=newChar // '1'
	(0,0,1,0,0
	,0,1,1,0,0
	,1,0,1,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,1,1,1,1,1);
	digits[18]=newChar // '2'
	(1,1,1,1,0
	,0,0,0,0,1
	,0,0,0,0,1
	,0,1,1,1,0
	,1,0,0,0,0
	,1,1,1,1,1);
	digits[19]=newChar // '3'
	(1,1,1,1,0
	,0,0,0,0,1
	,0,1,1,1,0
	,0,0,0,0,1
	,0,0,0,0,1
	,1,1,1,1,0);
	digits[20]=newChar // '4'
	(0,0,0,1,0
	,0,0,1,1,0
	,0,1,0,1,0
	,1,1,1,1,1
	,0,0,0,1,0
	,0,0,0,1,0);
	digits[21]=newChar // '5'
	(1,1,1,1,1
	,1,0,0,0,0
	,1,1,1,1,0
	,0,0,0,0,1
	,0,0,0,0,1
	,1,1,1,1,0);
	digits[22]=newChar // '6'
	(0,1,1,1,1
	,1,0,0,0,0
	,1,1,1,1,0
	,1,0,0,0,1
	,1,0,0,0,1
	,0,1,1,1,0);
	digits[23]=newChar // '7'
	(1,1,1,1,1
	,0,0,0,0,1
	,0,0,0,1,0
	,0,0,1,0,0
	,0,1,0,0,0
	,1,0,0,0,0);
	digits[24]=newChar // '8'
	(0,1,1,1,0
	,1,0,0,0,1
	,0,1,1,1,0
	,1,0,0,0,1
	,1,0,0,0,1
	,0,1,1,1,0);
	digits[25]=newChar // '9'
	(0,1,1,1,0
	,1,0,0,0,1
	,1,0,0,0,1
	,0,1,1,1,1
	,0,0,0,0,1
	,1,1,1,1,0);
	digits[26]=newChar // ':'
	(0,0,0,0,0
	,0,1,1,0,0
	,0,1,1,0,0
	,0,0,0,0,0
	,0,1,1,0,0
	,0,1,1,0,0);
	digits[27]=newChar // ';'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,0,1,0,0
	,0,0,0,0,0
	,0,0,1,0,0
	,0,1,0,0,0);
	digits[28]=newChar // '<'
	(0,0,0,0,0
	,0,0,1,1,0
	,0,1,0,0,0
	,1,0,0,0,0
	,0,1,0,0,0
	,0,0,1,1,0);
	digits[29]=newChar // '='
	(0,0,0,0,0
	,0,0,0,0,0
	,0,1,1,1,0
	,0,0,0,0,0
	,0,1,1,1,0
	,0,0,0,0,0);
	digits[30]=newChar // '>'
	(0,0,0,0,0
	,0,1,1,0,0
	,0,0,0,1,0
	,0,0,0,0,1
	,0,0,0,1,0
	,0,1,1,0,0);
	digits[31]=newChar // '?'
	(0,1,1,1,0
	,1,0,0,0,1
	,0,0,0,1,0
	,0,0,1,0,0
	,0,0,0,0,0
	,0,0,1,0,0);

	digits[32]=newChar // '@'
	(0,1,1,1,0
	,1,0,0,1,1
	,1,0,1,0,1
	,1,0,1,1,1
	,1,0,0,0,0
	,0,1,1,1,1);
	digits[33]=newChar // 'A'
	(0,0,1,0,0
	,0,1,0,1,0
	,1,0,0,0,1
	,1,1,1,1,1
	,1,0,0,0,1
	,1,0,0,0,1);
	digits[34]=newChar // 'B'
	(1,1,1,1,0
	,1,0,0,0,1
	,1,1,1,1,0
	,1,0,0,0,1
	,1,0,0,0,1
	,1,1,1,1,0);
	digits[35]=newChar // 'C'
	(0,1,1,1,1
	,1,0,0,0,0
	,1,0,0,0,0
	,1,0,0,0,0
	,1,0,0,0,0
	,0,1,1,1,1);
	digits[36]=newChar // 'D'
	(1,1,1,1,0
	,1,0,0,0,1
	,1,0,0,0,1
	,1,0,0,0,1
	,1,0,0,0,1
	,1,1,1,1,0);
	digits[37]=newChar // 'E'
	(1,1,1,1,1
	,1,0,0,0,0
	,1,1,1,1,0
	,1,0,0,0,0
	,1,0,0,0,0
	,1,1,1,1,1);
	digits[38]=newChar // 'F'
	(1,1,1,1,1
	,1,0,0,0,0
	,1,1,1,1,0
	,1,0,0,0,0
	,1,0,0,0,0
	,1,0,0,0,0);
	digits[39]=newChar // 'G'
	(0,1,1,1,1
	,1,0,0,0,0
	,1,0,0,0,0
	,1,0,1,1,1
	,1,0,0,0,1
	,0,1,1,1,0);
	digits[40]=newChar // 'H'
	(1,0,0,0,1
	,1,0,0,0,1
	,1,1,1,1,1
	,1,0,0,0,1
	,1,0,0,0,1
	,1,0,0,0,1);
	digits[41]=newChar // 'I'
	(1,1,1,1,1
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,1,1,1,1,1);
	digits[42]=newChar // 'J'
	(1,1,1,1,1
	,0,0,0,1,0
	,0,0,0,1,0
	,0,0,0,1,0
	,1,0,0,1,0
	,0,1,1,0,0);
	digits[43]=newChar // 'K'
	(1,0,0,1,1
	,1,0,1,0,0
	,1,1,0,0,0
	,1,1,0,0,0
	,1,0,1,0,0
	,1,0,0,1,1);	
	digits[44]=newChar // 'L'
	(1,0,0,0,0
	,1,0,0,0,0
	,1,0,0,0,0
	,1,0,0,0,0
	,1,0,0,0,0
	,1,1,1,1,1);
	digits[45]=newChar // 'M'
	(1,0,0,0,1
	,1,1,0,1,1
	,1,0,1,0,1
	,1,0,0,0,1
	,1,0,0,0,1
	,1,0,0,0,1);
	digits[46]=newChar // 'N'
	(1,0,0,0,1
	,1,1,0,0,1
	,1,0,1,0,1
	,1,0,1,0,1
	,1,0,0,1,1
	,1,0,0,0,1);
	digits[47]=newChar // 'O'
	(0,1,1,1,0
	,1,0,0,0,1
	,1,0,0,0,1
	,1,0,0,0,1
	,1,0,0,0,1
	,0,1,1,1,0);
	digits[48]=newChar // 'P'
	(1,1,1,1,0
	,1,0,0,0,1
	,1,0,0,0,1
	,1,1,1,1,0
	,1,0,0,0,0
	,1,0,0,0,0);
	digits[49]=newChar // 'Q'
	(0,1,1,1,0
	,1,0,0,0,1
	,1,0,0,0,1
	,1,0,1,0,1
	,1,0,0,1,0
	,0,1,1,0,1);
	digits[50]=newChar // 'R'
	(1,1,1,1,0
	,1,0,0,0,1
	,1,0,0,0,1
	,1,1,1,1,0
	,1,0,0,1,0
	,1,0,0,0,1);
	digits[51]=newChar // 'S'
	(0,1,1,1,1
	,1,0,0,0,0
	,0,1,1,1,0
	,0,0,0,0,1
	,0,0,0,0,1
	,1,1,1,1,0);
	digits[52]=newChar // 'T'
	(1,1,1,1,1
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,0,0);
	digits[53]=newChar // 'U'
	(1,0,0,0,1
	,1,0,0,0,1
	,1,0,0,0,1
	,1,0,0,0,1
	,1,0,0,0,1
	,0,1,1,1,0);
	digits[54]=newChar // 'V'
	(1,0,0,0,1
	,1,0,0,0,1
	,1,0,0,0,1
	,0,1,0,1,0
	,0,1,0,1,0
	,0,0,1,0,0);
	digits[55]=newChar // 'W'
	(1,0,0,0,1
	,1,0,0,0,1
	,1,0,0,0,1
	,1,0,1,0,1
	,1,1,0,1,1
	,1,0,0,0,1);
	digits[56]=newChar // 'X'
	(1,0,0,0,1
	,0,1,0,1,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,1,0,1,0
	,1,0,0,0,1);
	digits[57]=newChar // 'Y'
	(1,0,0,0,1
	,1,0,0,0,1
	,0,1,0,1,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,0,0);
	digits[58]=newChar // 'Z'
	(1,1,1,1,1
	,0,0,0,1,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,1,0,0,0
	,1,1,1,1,1);
	digits[59]=newChar // '['
	(0,1,1,1,0
	,0,1,0,0,0
	,0,1,0,0,0
	,0,1,0,0,0
	,0,1,0,0,0
	,0,1,1,1,0);
	digits[60]=newChar // '\'
	(1,0,0,0,0
	,0,1,0,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,0,1,0
	,0,0,0,0,1);
	digits[61]=newChar // ']'
	(0,1,1,1,0
	,0,0,0,1,0
	,0,0,0,1,0
	,0,0,0,1,0
	,0,0,0,1,0
	,0,1,1,1,0);
	digits[62]=newChar // '^'
	(0,0,1,0,0
	,0,1,0,1,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0);
	digits[63]=newChar // '_'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,1,1,1,0);

	digits[64]=newChar // '`'
	(0,1,0,0,0
	,0,0,1,0,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0
	,0,0,0,0,0);
	digits[65]=newChar // 'a'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,1,1,1,0
	,0,0,0,1,0
	,0,1,1,1,0
	,0,1,1,1,1);
	digits[66]=newChar // 'b'
	(0,0,0,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,1,1
	,0,0,1,0,1
	,0,1,1,1,1);
	digits[67]=newChar // 'c'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,0,1,1,1
	,0,1,0,0,0
	,0,1,0,0,0
	,0,0,1,1,1);
	digits[68]=newChar // 'd'
	(0,0,0,0,0
	,0,0,0,1,0
	,0,0,0,1,0
	,0,1,1,1,0
	,0,1,0,1,0
	,0,1,1,1,1);
	digits[69]=newChar // 'e'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,0,1,1,0
	,0,1,1,1,1
	,0,1,0,0,0
	,0,0,1,1,1);
	digits[70]=newChar // 'f'
	(0,0,0,0,0
	,0,0,1,1,0
	,0,0,1,0,1
	,0,0,1,0,0
	,0,1,1,1,0
	,0,0,1,0,0);
	digits[71]=newChar // 'g'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,0,1,1,1
	,0,1,0,0,1
	,0,0,1,1,1
	,0,1,1,1,1);
	digits[72]=newChar // 'h'
	(0,0,0,0,0
	,0,1,0,0,0
	,0,1,0,0,0
	,0,1,1,1,0
	,0,1,0,0,1
	,0,1,0,0,1);
	digits[73]=newChar // 'i'
	(0,0,0,0,0
	,0,0,1,0,0
	,0,0,0,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,1,0);
	digits[74]=newChar // 'j'
	(0,0,0,0,0
	,0,0,0,1,0
	,0,0,0,0,0
	,0,0,0,1,0
	,0,0,0,1,0
	,0,0,1,1,0);
	digits[75]=newChar // 'k'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,1,0,1,1
	,0,1,1,0,0
	,0,1,1,0,0
	,0,1,0,1,1);
	digits[76]=newChar // 'l'
	(0,0,0,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,1,0);
	digits[77]=newChar // 'm'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,0,1,1,1
	,0,1,1,1,1
	,0,1,0,1,1
	,0,1,0,0,1);
	digits[78]=newChar // 'n'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,0,1,1,1
	,0,1,0,0,1
	,0,1,0,0,1
	,0,1,0,0,1);
	digits[79]=newChar // 'o'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,0,1,1,0
	,0,1,0,0,1
	,0,1,0,0,1
	,0,0,1,1,0);
	digits[80]=newChar // 'p'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,1,1,1,1
	,0,0,1,0,1
	,0,0,1,1,1
	,0,0,1,0,0);
	digits[81]=newChar // 'q'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,1,1,1,1
	,0,1,0,1,0
	,0,1,1,1,0
	,0,0,0,1,0);
	digits[82]=newChar // 'r'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,1,0,1,1
	,0,0,1,0,1
	,0,0,1,0,0
	,0,0,1,0,0);
	digits[83]=newChar // 's'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,0,1,1,1
	,0,1,1,0,0
	,0,0,0,1,1
	,0,1,1,1,0);
	digits[84]=newChar // 't'
	(0,0,0,0,0
	,0,0,1,0,0
	,0,1,1,1,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,1,0);
	digits[85]=newChar // 'u'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,1,0,0,1
	,0,1,0,0,1
	,0,1,0,0,1
	,0,0,1,1,0);
	digits[86]=newChar // 'v'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,1,0,0,1
	,0,1,0,0,1
	,0,0,1,0,1
	,0,0,0,1,0);
	digits[87]=newChar // 'w'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,1,0,0,1
	,0,1,0,0,1
	,0,1,1,1,1
	,0,0,1,1,0);
	digits[88]=newChar // 'x'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,1,0,0,1
	,0,0,1,1,0
	,0,0,1,1,0
	,0,1,0,0,1);
	digits[89]=newChar // 'y'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,1,0,0,1
	,0,0,1,1,0
	,0,0,0,1,0
	,0,1,1,1,0);
	digits[90]=newChar // 'z'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,1,1,1,1
	,0,0,0,1,0
	,0,0,1,0,0
	,0,1,1,1,1);
	digits[91]=newChar // '{'
	(0,0,1,1,0
	,0,1,0,0,0
	,0,1,0,0,0
	,1,1,0,0,0
	,0,1,0,0,0
	,0,0,1,1,0);
	digits[92]=newChar // '|'
	(0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,0,0
	,0,0,1,0,0);
	digits[93]=newChar // '}'
	(0,1,1,0,0
	,0,0,0,1,0
	,0,0,0,1,0
	,0,0,0,1,1
	,0,0,0,1,0
	,0,1,1,0,0);
	digits[94]=newChar // '~'
	(0,0,0,0,0
	,0,0,0,0,0
	,0,0,1,0,1
	,0,1,0,1,0
	,0,0,0,0,0
	,0,0,0,0,0);
	digits[95]=newChar // ''
	(1,1,1,1,1
	,1,1,1,1,1
	,1,1,1,1,1
	,1,1,1,1,1
	,1,1,1,1,1
	,1,1,1,1,1);

	alphabetsize=96;
	alphabet=digits;
}

// liberta a mem�ia alocada pelas vari�eis do alfabeto
void freeAlphabet(){
	int i,n=alphabetsize;
	if(alphabet==NULL) return;
	for(i=0;i<n;i++) free(alphabet[i]);
	free(alphabet);
	alphabet=NULL;
}

// converte o formato dos caracteres para formato bin�io com 30 bits (1 bit por pixel)
// ex: digit[2]==1 -> 000000000000000000000000000100 (base 2) -> 00000004 (base 10)
int *convertAlphabetToBinary(){
	int i,j,bindigit = 0;
	int *digit=NULL;
	int *bindigits=(int *)malloc(alphabetsize*sizeof(int));
	for(i=0;i<alphabetsize;i++){
		digit=(int *)alphabet[i];
		bindigit=0;
		for(j=29;j>=0;j--){
			bindigit = bindigit | (int)digit[j];
			bindigit = bindigit << 1;
		}
		bindigits[i]=bindigit;
		printf("@[%d] = %d = %u = %#x\n",i,bindigits[i],bindigit,bindigit);
	}
	return bindigits;
}

// devolve o �dice do caracter ASCII dado no alfabeto a usar
int alphabetCode(char c){
	int charcode=(int)c;
	if(charcode>=32 && charcode<=126) return (charcode-32);
	return (alphabetsize-1);
}

// desenha o caracter na posi�o indicada da imagem
void drawChar(int n, int x, int y){
	int i,j,k;
	int *digit;
	int **digits;
	digits=alphabet;
	if(digits==NULL) return;
	digit=digits[n];
	for(i=0;i<6;i++)
		for(j=0;j<5;j++){
			k=5*i+j;
			if(digit[k]!=0)	drawPoint(x+j,y+i,black);
		}
}

// mostra todos os caracteres do alfabeto no topo da imagem
void printAlphabet(){
	int i,n=alphabetsize;
	int x=2,y=2;
	for(i=0;i<n;i++) drawChar(i,x+6*i,y);
}

// imprime os caracteres ASCII e o seu respectivo c�igo
void printAsciiChars(){
	int i;
	for(i=32;i<=255;i++){ // '0' -> 48
		printf("%c => %d\n",(char)i,(int)i);
	}
}

// devolve o nmero de algarismos de um nmero
int digitCount(int n){
	int left,count;
	if(n==0) return 1;
	left=n;
	count=0;
	if(n<0){
		left=-n;
		count=1;
	}
	while(left!=0){
		left=left/10;
		count++;
	}
	return count;
}

// desenha o nmero na posi�o indicada da imagem
void drawNumber(int n, int x, int y){
	int left,i,k;
	int *digits;
	k=digitCount(n);
	if(n==0) {
		drawChar(0+16,x,y);
		return;
	}
	digits=(int *)malloc(k*sizeof(int));
	left=n;
	if(n<0){
		left=-n;
		digits[k-1]=10;
	}
	i=0;
	while(left!=0){
		digits[i]=left%10;
		left=left/10;
		i++;
	}
	for(i=0;i<k;i++) drawChar(digits[k-i-1]+16,x+i*5+i,y);
	free(digits);
}

// desenha o nmero decimal na posi�o indicada da imagem
void drawDecimalNumber(double n, int x, int y){
        int intpart,decpart;
	int intpartsize;
	intpart=(int)floor(n);
	decpart=(int)floor((n-(double)intpart)*10);
	drawNumber(intpart,x,y);
	intpartsize=digitCount(intpart);
	drawText(".",x+6*intpartsize,y);
	drawNumber(decpart,x+6*(intpartsize+1),y);
}

// desenha o nmero alinhado �esquerda (termina na posi�o dada)
void drawNumberAtLeft(int n, int x, int y){
	drawNumber(n,x-digitCount(n)*6,y);
}

// desenha o nmero alinhado ao centro (da posi�o dada)
void drawNumberAtCenter(int n, int x, int y){
	drawNumber(n,x-digitCount(n)*3,y);
}

// desenha o texto na posi�o indicada da imagem
void drawText(char *text, int x, int y){
	int i;
	int charpos;
	char currentchar;
	if(alphabet==NULL) initializeAlphabet();
	currentchar=text[0];
	for(i=0;currentchar!='\0';i++){
		charpos=alphabetCode(currentchar);
		drawChar(charpos,x+6*i,y);
		currentchar=text[i+1];
	}
}

// TODO: draw over white background
// TODO: check box margins and add "..." as needed
// desenha o texto na posi�o indicada da imagem
void drawTextNoBackground(char *text, int x, int y){
	int i,n;
	int charpos;
	char currentchar;
	if(alphabet==NULL) initializeAlphabet();
	n=(int)strlen(text);
	drawFilledRectangle(x-1,y-1,6*n+2,6+2,white);
	currentchar=text[0];
	for(i=0;currentchar!='\0';i++){
		charpos=alphabetCode(currentchar);
		drawChar(charpos,x+6*i,y);
		currentchar=text[i+1];
	}
}

// desenha o texto alinhado �esquerda (termina na posi�o dada)
void drawTextAtLeft(char *text, int x, int y){
	drawText(text,x-(int)strlen(text)*(5+1),y);
}

// desenha o texto alinhado ao centro
void drawTextAtCenter(char *text, int x, int y){
	int n=((int)strlen(text))*(5+1);
	drawText(text,x-(n/2),y);
}

// draws text with a circular shape
void drawTextAtCircleTop(char *text, int xc, int yc, int r){
	int i,n,m,x,y;
	if(alphabet==NULL) initializeAlphabet();
	m=(int)strlen(text);
	n=(m*(5+1));
	if(n>(2*r)) {m=((2*r)/(5+1));n=(m*(5+1));}
	x=(-(n/2));
	for(i=0;i<m;i++){
		y=(int)sqrt((double)(r*r-x*x));
		//drawFilledRectangle((xc+x-1),(yc-y-1),5+2,6+2,white);
		drawChar(alphabetCode(text[i]),(xc+x),(yc-y));
		x+=6;
	}
}


// FALTA: melhorar com strcats
// desenha uma legenda centrada no topo da imagem
void drawLabel(int l, int p, int i, char *seqtext){
	int n;
	char *labeltext;
	n=(int)strlen(seqtext);
	labeltext=(char *)malloc((20+n)*sizeof(char));
	if(i<=0 && l<=0) sprintf(labeltext,"*=%d ",p);
	else if(l<=0) sprintf(labeltext,"*=%d XI=%d ",p,i);
	else if(p<=0) sprintf(labeltext,"L=%d XI=%d ",l,i);
	else if(i<=0) sprintf(labeltext,"L=%d *=%d ",l,p);
	else sprintf(labeltext,"L=%d *=%d XI=%d ",l,p,i);
	if(seqtext!=NULL && n!=0) strcat(labeltext,seqtext);
	n=((int)strlen(labeltext))*(5+1);
	drawText(labeltext,(bmpwidth/2)-(n/2),10);
	free(labeltext);
}


// ********
// GRAPHICS
// ********

// initializes the bitmap, colors and alphabet variables
void initializeGraphics(int w, int h){
	initializeBitmap(w,h,0);
	initializeColors();
	initializeAlphabet();
	bmpwidth=getBitmapWidth();
	bmpheight=getBitmapHeight();
}

// saves the bitmap file and frees initialized variables
void finalizeGraphics(char *filename){
	saveBitmap(filename);
	free(boxstarty);
	free(boxmiddlepoint);
	free(numcolorspersizeunit);
}

// draws the outline of the sequence
void drawEmptyBlock(int startpos, int length, int boxid, uint8_t color){
	drawRectangle(boxstartx+double2int(startpos*pixelsperunit),boxstarty[boxid],double2int(length*pixelsperunit)+3,boxheight,color); //  with +3 pixels for overflow at right
}

// draws the inside of the sequence
void drawFilledBlock(int startpos, int length, int boxid, uint8_t color){
	drawFilledRectangle(boxstartx+double2int(startpos*pixelsperunit)+1,boxstarty[boxid]+1,double2int(length*pixelsperunit)-2+3,boxheight-2,color); //  with +3 pixels for overflow at right
}

// TODO: use boxcolor;
// draws a coloured box inside a sequence in its original position
void drawBlock(int startpos, int length, int boxid){
	//uint8_t color=currentcolor;
	if(currentcolor==0){
		currentcolor=(startpos/(sequencesize[boxid]/numberofcolors));
		if(currentcolor<3) currentcolor+=3;
	}
	drawFilledRectangle(boxstartx+double2int(startpos*pixelsperunit)+1,boxstarty[boxid]+1,double2intup(length*pixelsperunit),boxheight-2,currentcolor);
	boxmiddlepoint[boxid]=boxstartx+double2intup((startpos+length/2)*pixelsperunit);
}

// returns the red, green and blue components of the color of the last block
void getRGBColor(int *rgbarray){
	if(rgbarray==NULL) return;
	rgbarray[0]=getColorComponent(currentcolor,'r');
	rgbarray[1]=getColorComponent(currentcolor,'g');
	rgbarray[2]=getColorComponent(currentcolor,'b');
}

// TODO: use boxcolor;
// draws a coloured box inside a sequence in its rotated position
// returns the new position already rotated
int drawBlockRotated(int startpos, int length, int boxid){
	//uint8_t color=currentcolor;
	int newstartpos=startpos;
	if(startpos<startposition[boxid]) newstartpos+=sequencesize[boxid];
	newstartpos-=startposition[boxid];
	if(currentcolor==0){
		//currentcolor=(startpos/(sequencesize[boxid]/numberofcolors));
		currentcolor=(uint8_t)double2int(((double)startpos)*numcolorspersizeunit[boxid]);
		if(currentcolor<3) currentcolor+=3;
	}
	drawFilledRectangle(boxstartx+double2int(newstartpos*pixelsperunit)+1,boxstarty[boxid]+1,double2intup(length*pixelsperunit),boxheight-2,currentcolor);
	boxmiddlepoint[boxid]=boxstartx+double2intup((newstartpos+length/2)*pixelsperunit);
	return newstartpos;
}

// connects the current block on all sequences with a line
void connectBlocks(){
	int i;
	if(boxmiddlepoint[0]==0) return;
	for(i=0;i<(numberofsequences-1);i++){
		drawLine(boxmiddlepoint[i],boxstarty[i]+boxheight,boxmiddlepoint[(i+1)],boxstarty[(i+1)],currentcolor);
		boxmiddlepoint[i]=0;
	}
	boxmiddlepoint[i]=0;
	//currentcolor++;
	//if(currentcolor>=255) currentcolor=2;
	currentcolor=0;
}

// draws the box and the description for each sequence
void drawLabels(char **labels){
	int i;
	for(i=0;i<numberofsequences;i++){
		drawTextNoBackground(labels[i],(boxstartx+2),(boxstarty[i]+boxheight/2-3));
		drawEmptyBlock(0,sequencesize[i],i,black); // draw black frame around each sequence box
	}
}

// writes a text centered at the bottom of the image
void drawBottomLabel(char *label){
	drawTextAtCenter(label,bmpwidth/2,bmpheight-lowheadersize);
}

// TODO: dynamically change minimumboxlength
// initialize variables
void initializeBlocks(int n, int *sizes, int *positions, char *imagemapfilename){
	int i,maxsize,markstep,marksize,spaceformarks,numberofmarks,markseqpos,markxpos;
	numberofsequences=n;
	sequencesize=sizes;
	startposition=positions;
	//boxstartx=bmpwidth/20;
	boxstartx=32;
	headersize=25;
	lowheadersize=15;
	interseqsize=50;
	//boxheight=(bmpheight-2*headersize)/(n+2*(n-1));
	boxheight=32;
	bmpheight=2*headersize+n*boxheight+(n-1)*interseqsize;
	initializeGraphics(640,bmpheight);
	//currentcolor=2;
	currentcolor=0;
	numberofcolors=getBitmapNumberOfColors()-3; // from bitmap pallette initialization
	maxsize=sizes[0];
	for(i=1;i<n;i++) if(sizes[i]>maxsize) maxsize=sizes[i];
	pixelsperunit=((double)(bmpwidth-2*boxstartx))/((double)maxsize);
	boxmiddlepoint=(int *)calloc(n,sizeof(int));
	boxstarty=(int *)calloc(n,sizeof(int));
	numcolorspersizeunit=(double *)calloc(n,sizeof(double));
	for(i=0;i<n;i++){
		boxstarty[i]=headersize+i*boxheight+i*interseqsize;
		//boxstarty[i]=headersize+boxheight*(3*i);
		drawFilledBlock(0,sequencesize[i],i,2); // draw grey background of sequence boxes
		numcolorspersizeunit[i]=((double)numberofcolors)/((double)sequencesize[i]);
	}
	drawLine(boxstartx,lowheadersize,(bmpwidth-boxstartx),lowheadersize,black);
	//drawNumberAtCenter(maxsize,(bmpwidth-boxstartx),headersize-10);
	drawLine((bmpwidth-boxstartx),lowheadersize-3,(bmpwidth-boxstartx),lowheadersize+4,black);
	marksize=digitCount(maxsize)*6;
	spaceformarks=(bmpwidth-2*boxstartx);
	markstep=1;
	while(1){
		numberofmarks=maxsize/markstep+1;
		if((numberofmarks*marksize)<spaceformarks) break;
		markstep=markstep*5;
		numberofmarks=maxsize/markstep+1;
		if((numberofmarks*marksize)<spaceformarks) break;
		markstep=markstep*2;
	}
	for(i=0;i<numberofmarks;i++){
		markseqpos=(i*markstep);
		markxpos=boxstartx+double2int(markseqpos*pixelsperunit);
		drawNumberAtCenter(markseqpos,markxpos,lowheadersize-10);
		drawLine(markxpos,lowheadersize-3,markxpos,lowheadersize+4,black);
	}
	if((mapfile=fopen(imagemapfilename,"w"))==NULL){
		exitMessage("Can't write image map file");
	}
	fprintf(mapfile,"%d %d\n",(numberofmarks+1),numberofsequences);
	for(i=0;i<numberofmarks;i++) fprintf(mapfile,"%d ",(boxstartx+double2int((i*markstep)*pixelsperunit)));
	fprintf(mapfile,"%d",(boxstartx+double2int(maxsize*pixelsperunit)));
	fprintf(mapfile,"\n");
	for(i=0;i<numberofmarks;i++) fprintf(mapfile,"%d ",(i*markstep));
	fprintf(mapfile,"%d",maxsize);
	fprintf(mapfile,"\n");
	for(i=0;i<numberofsequences;i++) fprintf(mapfile,"%d ",boxstarty[i]);
	fprintf(mapfile,"\n");
	for(i=0;i<numberofsequences;i++) fprintf(mapfile,"%d ",(boxstarty[i]+boxheight));
	fprintf(mapfile,"\n");
	fclose(mapfile);
}

void DrawCircularAlignmentPlot(char *filename, char *imagefilename){
	FILE *file;
	int numseqs,seqsize,*seqsizes,**counts,i,j,k;
	char c,*label;
	fpos_t *seqstarts;
	int emptycentersize,bandsize,bandgapsize,*circlesizes;
	int xc,yc,npoints,pos,n,x,y,r,conservation,gaps;
	double **posperpoint,dx,dy;
	uint8_t color,conscolor,notconscolor,gapcolor,addcolor;
	int linesize,labelsize;
	double markinterval;
	if((file=fopen(filename,"r"))==NULL) {printf("> ERROR: Can't open alignment file\n");return;}
	numseqs=0;
	while((c=fgetc(file))!=EOF) if(c=='>') numseqs++;
	if(numseqs<2) {printf("> ERROR: Not enough sequences in file\n");return;}
	seqsizes=(int *)calloc(numseqs,sizeof(int));
	seqstarts=(fpos_t *)calloc(numseqs,sizeof(fpos_t));
	rewind(file);
	i=0;
	c=fgetc(file);
	while(1){
		if(c!='>') break;
		while(c!='\n' && c!=EOF) c=fgetc(file);
		fgetpos(file,&(seqstarts[i]));
		while(c!='>' && c!=EOF){
			if(c=='A'||c=='C'||c=='G'||c=='T'||c=='-'||c=='a'||c=='c'||c=='g'||c=='t'||c=='N'||c=='n') seqsizes[i]++;
			else if (c!='\n') seqsizes[i]++;
			c=fgetc(file);
		}
		i++;
	}
	seqsize=seqsizes[0];
	for(i=1;i<numseqs;i++) if(seqsizes[i]!=seqsize) {printf("> ERROR: Consensus sizes don't match\n");return;}
	free(seqsizes);
	counts=(int **)calloc(seqsize,sizeof(int *));
	for(i=0;i<seqsize;i++) counts[i]=(int *)calloc(5,sizeof(int));
	for(i=0;i<numseqs;i++){
		fsetpos(file,&(seqstarts[i]));
		j=0;
		while(j<seqsize){
			c=fgetc(file);
			if(((int)c)>90) c=(char)(c-32); // convert to uppercase if needed
			switch(c){
				case '-':
					counts[j][0]++;
					break;
				case 'A':
					counts[j][1]++;
					break;
				case 'C':
					counts[j][2]++;
					break;
				case 'G':
					counts[j][3]++;
					break;
				case 'T':
					counts[j][4]++;
					break;
				case '\n':
					j--;
					break;
				default:
					break;
			}
			j++;
		}
	}
	bandsize=5; // size of sequence plot
	bandgapsize=(2*bandsize); // size of interval between sequences
	emptycentersize=(bandgapsize*numseqs); // hollow middle
	if(emptycentersize<50) emptycentersize=100; // expand hollow middle if too small
	j=(2*(emptycentersize+numseqs*(bandsize+bandgapsize)+bandsize)+1); // diameter of all sequences
	j+=2*(6*digitCount(seqsize)+6); // add space to draw position numbers on both sides
	initializeBitmap(j,j,2);
	xc=((j+1)/2); // center coordinates
	yc=((j+1)/2);
	circlesizes=(int *)malloc(numseqs*sizeof(int));
	posperpoint=(double **)malloc(bandsize*sizeof(double *));
	for(k=0;k<bandsize;k++){ // draw one single line circle at a time to fill the entire band
		posperpoint[k]=(double *)malloc(numseqs*sizeof(double));
		for(i=0;i<numseqs;i++){
			circlesizes[i]=(emptycentersize+(numseqs-i)*(bandsize+bandgapsize)); // circle radius of each sequence (from outter to inner)
			circlesizes[i]-=k;
			r=circlesizes[i];
			npoints=0;
			y=r;
			for(x=0;x<=y;x++){
				y=(int)floor(sqrt((double)(r*r-x*x)));
				npoints++;
			}
			for(y=x;y>=0;y--){
				x=(int)floor(sqrt((double)(r*r-y*y)));
				npoints++;
			}
			npoints=(4*(npoints)-1); // perimeter of the circle
			//printf("\n(%d:%d)",i,npoints);
			if(npoints>seqsize) {printf("> ERROR: Sequence length is too short to draw correct circular plot.\n");return;}
			posperpoint[k][i]=(((double)seqsize)/((double)npoints)); // each pixel will store this many sequence positions
			npoints=1;
			j=0; // process the whole circle, a quarter of circle at a time: 90',45',315',225',135',90'
			y=(-r);
			fsetpos(file,&(seqstarts[i]));
			// top,right
			for(x=1;x<=(-y);x++){
				dy=-sqrt((double)(r*r-x*x));
				conservation=0;
				gaps=0;
				pos=(int)floor(npoints*posperpoint[k][i]); // get chars until this position
				n=(pos-j); // number of chars that this pixels represents
				while(j<pos){
					c=fgetc(file);
					if(((int)c)>90) c=(char)(c-32);
					switch(c){
						case '-':
							gaps++;
							break;
						case 'A':
							conservation+=counts[j][1];
							break;
						case 'C':
							conservation+=counts[j][2];
							break;
						case 'G':
							conservation+=counts[j][3];
							break;
						case 'T':
							conservation+=counts[j][4];
							break;
						case 'N':
							break;
						default:
							j--;
							break;
					}
					j++;
				} // gap color in blue, conservation color ranging from green to red
				conscolor=(uint8_t)floor(((double)(conservation*255))/((double)((numseqs)*n))); // if it had full conservation, the count would be numchars*numseqs
				gapcolor=(uint8_t)floor(((double)(gaps*255))/((double)n)); // if it was all gaps, the count would be numchars
				notconscolor=(uint8_t)(255-(conscolor+gapcolor)); // what is not conserved and is not gaps, is non-conserved
				if(conscolor>=notconscolor && conscolor>=gapcolor) addcolor=(255-conscolor); // make at least one of the colors maximum, so we do not get grey tones
				else if(notconscolor>=gapcolor) addcolor=(255-notconscolor);
				else addcolor=(255-gapcolor);
				color=getColorFromPalette( (uint8_t)(conscolor+addcolor) , (uint8_t)(notconscolor+addcolor) , (uint8_t)(gapcolor) ); // red, green, blue
				y=(int)ceil(dy); // draw both rounded up and rounded down points to avoid empty pixels in the middle
				drawPoint((xc+x),(yc+y),color);
				y=(int)floor(dy);
				drawPoint((xc+x),(yc+y),color);
				npoints++;
			}
			// right
			for(y=-(x+1);y<=x;y++){
				dx=sqrt((double)(r*r-y*y));
				conservation=0;
				gaps=0;
				pos=(int)floor(npoints*posperpoint[k][i]);
				n=(pos-j);
				while(j<pos){
					c=fgetc(file);
					if(((int)c)>90) c=(char)(c-32);
					switch(c){
						case '-':
							gaps++;
							break;
						case 'A':
							conservation+=counts[j][1];
							break;
						case 'C':
							conservation+=counts[j][2];
							break;
						case 'G':
							conservation+=counts[j][3];
							break;
						case 'T':
							conservation+=counts[j][4];
							break;
						case 'N':
							break;
						default:
							j--;
							break;
					}
					j++;
				}
				conscolor=(uint8_t)floor(((double)(conservation*255))/((double)((numseqs)*n)));
				gapcolor=(uint8_t)floor(((double)(gaps*255))/((double)n));
				notconscolor=(uint8_t)(255-(conscolor+gapcolor));
				if(conscolor>=notconscolor && conscolor>=gapcolor) addcolor=(255-conscolor);
				else if(notconscolor>=gapcolor) addcolor=(255-notconscolor);
				else addcolor=(255-gapcolor);
				color=getColorFromPalette( (uint8_t)(conscolor+addcolor) , (uint8_t)(notconscolor+addcolor) , (uint8_t)(gapcolor) );
				x=(int)ceil(dx);
				drawPoint((xc+x),(yc+y),color);
				x=(int)floor(dx);
				drawPoint((xc+x),(yc+y),color);
				npoints++;
			}
			// down
			for(x=(y-1);(-x)<=y;x--){
				dy=sqrt((double)(r*r-x*x));
				conservation=0;
				gaps=0;
				pos=(int)floor(npoints*posperpoint[k][i]);
				n=(pos-j);
				while(j<pos){
					c=fgetc(file);
					if(((int)c)>90) c=(char)(c-32);
					switch(c){
						case '-':
							gaps++;
							break;
						case 'A':
							conservation+=counts[j][1];
							break;
						case 'C':
							conservation+=counts[j][2];
							break;
						case 'G':
							conservation+=counts[j][3];
							break;
						case 'T':
							conservation+=counts[j][4];
							break;
						case 'N':
							break;
						default:
							j--;
							break;
					}
					j++;
				}
				conscolor=(uint8_t)floor(((double)(conservation*255))/((double)((numseqs)*n)));
				gapcolor=(uint8_t)floor(((double)(gaps*255))/((double)n));
				notconscolor=(uint8_t)(255-(conscolor+gapcolor));
				if(conscolor>=notconscolor && conscolor>=gapcolor) addcolor=(255-conscolor);
				else if(notconscolor>=gapcolor) addcolor=(255-notconscolor);
				else addcolor=(255-gapcolor);
				color=getColorFromPalette( (uint8_t)(conscolor+addcolor) , (uint8_t)(notconscolor+addcolor) , (uint8_t)(gapcolor) );
				y=(int)ceil(dy);
				drawPoint((xc+x),(yc+y),color);
				y=(int)floor(dy);
				drawPoint((xc+x),(yc+y),color);
				npoints++;
			}
			// left
			for(y=-(x-1);(-y)<=(-x);y--){
				dx=-sqrt((double)(r*r-y*y));
				conservation=0;
				gaps=0;
				pos=(int)floor(npoints*posperpoint[k][i]);
				n=(pos-j);
				while(j<pos){
					c=fgetc(file);
					if(((int)c)>90) c=(char)(c-32);
					switch(c){
						case '-':
							gaps++;
							break;
						case 'A':
							conservation+=counts[j][1];
							break;
						case 'C':
							conservation+=counts[j][2];
							break;
						case 'G':
							conservation+=counts[j][3];
							break;
						case 'T':
							conservation+=counts[j][4];
							break;
						case 'N':
							break;
						default:
							j--;
							break;
					}
					j++;
				}
				conscolor=(uint8_t)floor(((double)(conservation*255))/((double)((numseqs)*n)));
				gapcolor=(uint8_t)floor(((double)(gaps*255))/((double)n));
				notconscolor=(uint8_t)(255-(conscolor+gapcolor));
				if(conscolor>=notconscolor && conscolor>=gapcolor) addcolor=(255-conscolor);
				else if(notconscolor>=gapcolor) addcolor=(255-notconscolor);
				else addcolor=(255-gapcolor);
				color=getColorFromPalette( (uint8_t)(conscolor+addcolor) , (uint8_t)(notconscolor+addcolor) , (uint8_t)(gapcolor) );
				x=(int)ceil(dx);
				drawPoint((xc+x),(yc+y),color);
				x=(int)floor(dx);
				drawPoint((xc+x),(yc+y),color);
				npoints++;
			}
			// top,left
			for(x=(y+1);x<0;x++){
				dy=-sqrt((double)(r*r-x*x));
				conservation=0;
				gaps=0;
				pos=(int)floor(npoints*posperpoint[k][i]);
				n=(pos-j);
				while(j<pos){
					if(j==seqsize) break;
					c=fgetc(file);
					if(((int)c)>90) c=(char)(c-32);
					switch(c){
						case '-':
							gaps++;
							break;
						case 'A':
							conservation+=counts[j][1];
							break;
						case 'C':
							conservation+=counts[j][2];
							break;
						case 'G':
							conservation+=counts[j][3];
							break;
						case 'T':
							conservation+=counts[j][4];
							break;
						case '\n':
							j--;
							break;
						default:
							break;
					}
					j++;
				}
				conscolor=(uint8_t)floor(((double)(conservation*255))/((double)((numseqs)*n)));
				gapcolor=(uint8_t)floor(((double)(gaps*255))/((double)n));
				notconscolor=(uint8_t)(255-(conscolor+gapcolor));
				if(conscolor>=notconscolor && conscolor>=gapcolor) addcolor=(255-conscolor);
				else if(notconscolor>=gapcolor) addcolor=(255-notconscolor);
				else addcolor=(255-gapcolor);
				color=getColorFromPalette( (uint8_t)(conscolor+addcolor) , (uint8_t)(notconscolor+addcolor) , (uint8_t)(gapcolor) );
				y=(int)ceil(dy);
				drawPoint((xc+x),(yc+y),color);
				y=(int)floor(dy);
				drawPoint((xc+x),(yc+y),color);
				npoints++;
			}
			//printf("(%d:%d)",i,npoints);
		}
	}
	initializeColors();
	initializeAlphabet();
	bmpwidth=getBitmapWidth();
	bmpheight=getBitmapHeight();
	rewind(file);
	label=(char *)calloc(65,sizeof(char));
	c='\n';
	for(i=0;i<numseqs;i++){
		circlesizes[i]=(emptycentersize+(numseqs-i)*(bandsize+bandgapsize)); // radius of the outter line of the band of the i-th sequence
		drawVLine(xc,(yc-circlesizes[i]),(bandsize+1),grey);
		while(c!='>') c=getc(file);
		labelsize=0;
		while((c=getc(file))!='\n' && labelsize<64) label[labelsize++]=c;
		label[labelsize]='\0';
		drawTextAtCircleTop(label,xc,yc,(circlesizes[i]-bandsize-1)); // draw sequence description bellow the corresponding band
	}
	free(label);
	linesize=5;
	markinterval=(((double)seqsize)/8.0); // draw 8 position marks around the circle
	drawVLine(xc,(yc-circlesizes[0]-linesize),linesize,black);
	drawNumberAtCenter(0,(xc+1),(yc-circlesizes[0]-linesize-6-1));
	drawNumberAtCenter(seqsize,(xc+1),(yc-circlesizes[0]-linesize-6-1-6-1));
	drawVLine(xc,(yc+circlesizes[0]+1),linesize,black);
	drawNumberAtCenter(double2int(4.0*markinterval),xc,(yc+circlesizes[0]+linesize+1+1));
	drawHLine((xc-circlesizes[0]-linesize),yc,linesize,black);
	drawNumberAtLeft(double2int(6.0*markinterval),(xc-circlesizes[0]-linesize),(yc-3));
	drawHLine((xc+circlesizes[0]+1),yc,linesize,black);
	drawNumber(double2int(2.0*markinterval),(xc+circlesizes[0]+linesize+1+1),(yc-3));
	r=circlesizes[0]; // draw the diagonal marks
	y=r;
	for(x=0;x<=y;x++) y=(int)floor(sqrt((double)(r*r-x*x)));
	drawDLine((xc+x),(yc-y),-1.0,linesize,black);
	drawNumber(double2int(1.0*markinterval),(xc+x+linesize+1+1),(yc-y-linesize-6));
	drawDLine((xc+x),(yc+y),1.0,linesize,black);
	drawNumber(double2int(3.0*markinterval),(xc+x+linesize+1+1),(yc+y+linesize));
	drawDLine((xc-x-linesize),(yc+y+linesize),-1.0,linesize,black);
	drawNumberAtLeft(double2int(5.0*markinterval),(xc-x-linesize),(yc+y+linesize));
	drawDLine((xc-x-linesize),(yc-y-linesize),1.0,linesize,black);
	drawNumberAtLeft(double2int(7.0*markinterval),(xc-x-linesize),(yc-y-linesize-6));
	x=(bmpwidth-1-(5+1));
	y=(bmpheight-1-6*(6+1));
	n=(12*(5+1));
	drawTextAtLeft("Conservation",x,y); // draw bottom-right legend
	y+=(6+1);
	markinterval=(((double)(255))/((double)(n/2-1)));
	for(i=0;i<(n/2);i++){
		color=getColorFromPalette( (uint8_t)255 , (uint8_t)double2int(i*markinterval) , (uint8_t)0 );
		drawVLine((x-12*(5+1)+i),y,(6+1),color);
	}
	for(i=(n/2);i<n;i++){
		color=getColorFromPalette( (uint8_t)double2int(((n-1)-i)*markinterval) , (uint8_t)255 , (uint8_t)0 );
		drawVLine((x-12*(5+1)+i),y,(6+1),color);
	}
	drawTextAtLeft("+          -",x,y);
	y+=(2*(6+1));
	drawTextAtLeft("GapFrequency",x,y);
	y+=(6+1);
	for(i=0;i<n;i++){
		color=getColorFromPalette( (uint8_t)double2int((i/2)*markinterval) , (uint8_t)double2int((i/2)*markinterval) , (uint8_t)255 );
		drawVLine((x-12*(5+1)+i),y,(6+1),color);
	}
	drawTextAtLeft("+          -",x,y);
	/*
	if(rotations!=NULL){
			for(i=0;i<numseqs;i++){
				circlesizes[i]=(emptycentersize+(numseqs-i)*(bandsize+bandgapsize)); // radius of the outter line of the band of the i-th sequence

			}
	}
	*/
	saveBitmap(imagefilename);
	free(circlesizes);
	for(i=0;i<bandsize;i++) free(posperpoint[i]);
	free(posperpoint);
	for(i=0;i<seqsize;i++) free(counts[i]);
	free(counts);
	free(seqstarts);
	fclose(file);
}

void TestColorCircles(){
	int xc,yc,r,size,x,y,i,j,s;
	int ncolors,color,npoints,point;
	double m,ncolorsperpoint;
	initializeGraphics(320,240);
	xc=(bmpwidth/2);
	yc=(bmpheight/2);
	if(bmpwidth<bmpheight) size=(bmpwidth/6);
	else size=(bmpheight/6);
	r=(2*size);
	ncolors=(getBitmapNumberOfColors()-3);
	npoints=(int)floor(4.0*M_SQRT2*((double)r));
	ncolorsperpoint=(((double)ncolors)/npoints);
	point=0;
	y=r;
	// up,right
	for(x=0;x!=y;x++){
		y=(int)floor(sqrt((double)(r*r-x*x))+0.5);
		m=((double)x)/((double)y);
		s=0;
		i=x;
		color=(3+(int)floor(((double)point)*ncolorsperpoint));
		for(j=y;(s<size);j--){
			drawPoint((xc+i),(yc-j),color);
			i=(int)floor(m*j+0.5);
			s=(int)floor(sqrt((double)((x-i)*(x-i)+(y-j)*(y-j))));
		}
		point++;
	}
	// right,up
	for(y=(x-1);y>0;y--){
		x=(int)floor(sqrt((double)(r*r-y*y))+0.5);
		m=((double)y)/((double)x);
		s=0;
		j=y;
		color=(3+(int)floor(((double)point)*ncolorsperpoint));
		for(i=x;(s<size);i--){
			drawPoint((xc+i),(yc-j),color);
			j=(int)floor(m*i+0.5);
			s=(int)floor(sqrt((double)((x-i)*(x-i)+(y-j)*(y-j))));
		}
		point++;
	}
	// right,down
	for(y=0;y!=x;y++){
		x=(int)floor(sqrt((double)(r*r-y*y))+0.5);
		m=((double)y)/((double)x);
		s=0;
		j=y;
		color=(3+(int)floor(((double)point)*ncolorsperpoint));
		for(i=x;(s<size);i--){
			drawPoint((xc+i),(yc+j),color);
			j=(int)floor(m*i+0.5);
			s=(int)floor(sqrt((double)((x-i)*(x-i)+(y-j)*(y-j))));
		}
		point++;
	}
	// down,right
	for(x=(y-1);x>0;x--){
		y=(int)floor(sqrt((double)(r*r-x*x))+0.5);
		m=((double)x)/((double)y);
		s=0;
		i=x;
		color=(3+(int)floor(((double)point)*ncolorsperpoint));
		for(j=y;(s<size);j--){
			drawPoint((xc+i),(yc+j),color);
			i=(int)floor(m*j+0.5);
			s=(int)floor(sqrt((double)((x-i)*(x-i)+(y-j)*(y-j))));
		}
		point++;
	}
	// down,left
	for(x=0;x!=y;x++){
		y=(int)floor(sqrt((double)(r*r-x*x))+0.5);
		m=((double)x)/((double)y);
		s=0;
		i=x;
		color=(3+(int)floor(((double)point)*ncolorsperpoint));
		for(j=y;(s<size);j--){
			drawPoint((xc-i),(yc+j),color);
			i=(int)floor(m*j+0.5);
			s=(int)floor(sqrt((double)((x-i)*(x-i)+(y-j)*(y-j))));
		}
		point++;
	}
	// left,down
	for(y=(x-1);y>0;y--){
		x=(int)floor(sqrt((double)(r*r-y*y))+0.5);
		m=((double)y)/((double)x);
		s=0;
		j=y;
		color=(3+(int)floor(((double)point)*ncolorsperpoint));
		for(i=x;(s<size);i--){
			drawPoint((xc-i),(yc+j),color);
			j=(int)floor(m*i+0.5);
			s=(int)floor(sqrt((double)((x-i)*(x-i)+(y-j)*(y-j))));
		}
		point++;
	}
	// left,up
	for(y=0;y!=x;y++){
		x=(int)floor(sqrt((double)(r*r-y*y))+0.5);
		m=((double)y)/((double)x);
		s=0;
		j=y;
		color=(3+(int)floor(((double)point)*ncolorsperpoint));
		for(i=x;(s<size);i--){
			drawPoint((xc-i),(yc-j),color);
			j=(int)floor(m*i+0.5);
			s=(int)floor(sqrt((double)((x-i)*(x-i)+(y-j)*(y-j))));
		}
		point++;
	}
	// up,left
	for(x=(y-1);x>0;x--){
		y=(int)floor(sqrt((double)(r*r-x*x))+0.5);
		m=((double)x)/((double)y);
		s=0;
		i=x;
		color=(3+(int)floor(((double)point)*ncolorsperpoint));
		for(j=y;(s<size);j--){
			drawPoint((xc-i),(yc-j),color);
			i=(int)floor(m*j+0.5);
			s=(int)floor(sqrt((double)((x-i)*(x-i)+(y-j)*(y-j))));
		}
		point++;
	}
	saveBitmap("circletest.bmp");
}


// desenha linhas diagonais de teste
void testLines(){
	drawLine(bmpwidth/2-16,bmpheight,bmpwidth/2+16,0,green);// D/L->U/R
	drawLine(bmpwidth/2+16,bmpheight,bmpwidth/2-16,0,red);// D/R->U/L
	drawLine(bmpwidth/2-32,0,bmpwidth/2+32,bmpheight,cyan);// U/L->D/R
	drawLine(bmpwidth/2+32,0,bmpwidth/2-32,bmpheight,yellow);// U/R->D/L
	drawLine(0,bmpheight/2+16,bmpwidth,bmpheight/2-16,green);// D/L->U/R
	drawLine(bmpwidth,bmpheight/2+16,0,bmpheight/2-16,red);// D/R->U/L
	drawLine(0,bmpheight/2-32,bmpwidth,bmpheight/2+32,cyan);// U/L->D/R
	drawLine(bmpwidth,bmpheight/2-32,0,bmpheight/2+32,yellow);// U/R->D/L
}

// ****
// MAIN
// ****

// FALTA: gr�icos com mltiplas fun�es
// FALTA: gr�ico de representa�o das frequ�cias relativas das letras
// FALTA: eixos com nmeros decimais
// FALTA: c�culo autom�ico da largura e altura da imagem
/*
int main(int argc, char *argv[]){
	char *file;
	int i, n, *xlist, *ylist;

	file=(char *)malloc(255*sizeof(char));
	strcpy(file,"plot.bmp");
	
	initializeGraphics(640,480);
	testBitmap(1);
	saveBitmap(file);
	return 1;
	//initializeGraphics(400,320);
	//initializeGraphics(32,32);

	//printAlphabet();
	//drawText("N=0 XI=ACTGU L=0 P=0",2,2);
	
	if(argc>3){
		n=atoi(argv[1]);
		printf("n=%d\n",n);
		xlist=(int *)malloc(n*sizeof(int));
		ylist=(int *)malloc(n*sizeof(int));
		for(i=0;(i<n && (2+2*i+1)<=argc); i++){
			xlist[i]=atoi(argv[2+2*i]);
			ylist[i]=atoi(argv[2+2*i+1]);
			printf("[%d]=(%d,%d)\n",i,xlist[i],ylist[i]);
		}
		n=i;
		printf("n=%d\n",n);
		//drawPlot(intlist2doublelist(xlist,n),intlist2doublelist(ylist,n),n,1);
	}
	else{
		//testHistogram(4);
		//testPlot(-32,64);
		//testPlot(-6,13);
		//testLines();
		//testBitmap(0);
	}
	saveBitmap(file);
	//showFileHexData(file);
	//freeBitmap();
	//system("pause");
	return 1;
}
*/

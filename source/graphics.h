// ********
// GRAPHICS
// ********

// inicializa as variaveis da imagem, do alfabeto e das cores
void initializeGraphics(int w, int h);

// desenha uma legenda centrada no topo da imagem
void drawLabel(int l, int p, int i, char *seqtext);

// desenha o texto na posicao indicada da imagem
void drawText(char *text, int x, int y);

// desenha o numero na posicao indicada da imagem
void drawNumber(int n, int x, int y);

// desenha o numero alinhado a esquerda (termina na posicao dada)
void drawNumberAtLeft(int n, int x, int y);

// desenha o numero alinhado ao centro (da posicao dada)
void drawNumberAtCenter(int n, int x, int y);

// liberta a memoria alocada pelas variaveis do alfabeto
void freeAlphabet();

// devolve o numero de algarismos de um numero
int digitCount(int n);

void initializeBlocks(int n, int *sizes, int *positions, char *imagemapfilename);
void drawBlock(int startpos, int length, int boxid);
int drawBlockRotated(int startpos, int length, int boxid);
void connectBlocks();
void drawLabels(char **labels);
void finalizeGraphics(char *filename);
void drawBottomLabel(char *label);
void getRGBColor(int *rgbarray);
void DrawCircularAlignmentPlot(char *filename, char *imagefilename);

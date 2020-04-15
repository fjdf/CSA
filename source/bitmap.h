#include <stdint.h>

// mostra na imagem todas as cores dispon�veis na palete
void testBitmap(int method);

// comprime a BitmapData com o algoritmo Run Length Encoding
void compressBitmapData();

// mostra o conte�do bin�rio do ficheiro
int showFileHexData(char *filename);
int showBitmapInfo(char *filename);

// inicializa vari�veis do bitmap
void initializeBitmap(int w, int h, int colorset);

// salva o Bitmap para um ficheiro
int saveBitmap(char *filename);

// liberta a mem�ria alocada pelo Bitmap
void freeBitmap();

// devolve a altura da imagem em pixels
int getBitmapHeight();

// devolve a largura da imagem em pixels
int getBitmapWidth();

int getBitmapNumberOfColors();
int getColorComponent(uint8_t colorpos, char rgbchar);

// devolve a posi��o na palette da cor mais pr�xima da cor RGB dada
uint8_t getColorFromPalette(uint8_t r, uint8_t g, uint8_t b);

// devolve a posi��o na BitmapData das coordenadas (x,y) na imagem
int dataPosition(int x, int y);

// define a cor do pixel na posi��o (x,y)
void drawPoint(int x, int y, uint8_t colorpos);

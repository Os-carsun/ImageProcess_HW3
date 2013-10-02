#include <stdio.h>
#include <stdlib.h>

#define TRUE 1
#define FALSE 0
typedef unsigned char BYTE;

//�w�q�Ϥ���T����
typedef	struct STRUCT_IMAGE{
	int width,height;
	int depth;
	void* pixels;
}ImageData;
//�̤j255
#define PIXELMAX 255
typedef struct STRUCT_PIXEL{
	int r,g,b;
}Pixel;

typedef unsigned long       DWORD;
typedef int                 BOOL;
typedef unsigned short      WORD;
typedef unsigned long       LONG;

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

//�w�q���Y���O(��KŪ�g
typedef struct tagBITMAPFILEHEADER {
        WORD    bfType;
        DWORD   bfSize;
        WORD    bfReserved1;
        WORD    bfReserved2;
        DWORD   bfOffBits;
} BITMAPFILEHEADER, *PBITMAPFILEHEADER;
//�w�q���Y���O(��KŪ�g
typedef struct tagBITMAPCOREHEADER {
        DWORD   bcSize;
        WORD    bcWidth;
        WORD    bcHeight;
        WORD    bcPlanes;
        WORD    bcBitCount;
} BITMAPCOREHEADER, *PBITMAPCOREHEADER;
//�w�q���Y���O(��KŪ�g
typedef struct tagBITMAPINFOHEADER{
        DWORD      biSize;
        LONG       biWidth;
        LONG       biHeight;
        WORD       biPlanes;
        WORD       biBitCount;
        DWORD      biCompression;
        DWORD      biSizeImage;
        LONG       biXPelsPerMeter;
        LONG       biYPelsPerMeter;
        DWORD      biClrUsed;
        DWORD      biClrImportant;
} BITMAPINFOHEADER, *PBITMAPINFOHEADER;
//�̤j��256��
#define MAXCOLORS 256



ImageData* createImage(int width,int height,int depth);
void disposeImage(ImageData *img);
int readBMPfile(char *fname,ImageData **img);
int writeBMPfile(char *fname,ImageData *img);
int getPixel(ImageData *img,int x,int y,Pixel *pix);
int setPixel(ImageData *img,int x,int y,Pixel *pix);


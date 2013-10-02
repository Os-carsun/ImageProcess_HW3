#include "image.h"

ImageData* createImage(int width,int height,int depth)
{
//�إߤ@�i�s���Ϥ�
  ImageData *newimg;//�ΨӦs���s�Ϥ������c����
  int byte_per_pixel;//�s�C��pixel���h��BYTE

  if(width<0 || height<0) return NULL;//�Y���e���t�ƫh������(�^��NULL
  if(depth!=8 && depth!=24) return NULL;//�u��O�Ƕ��Υ��m
  
  newimg=(ImageData*)malloc(sizeof(ImageData));//���t�@�q�Ŷ��Ӧs�Ϥ�(���Ъ��Ŷ�
  if(newimg==NULL) return NULL;//�Y���t�Ŷ����Ѧ^��NULL

  byte_per_pixel=depth/8;//�ھ�depth�Ӥ��R�C��pixel��BYTE�j�p(���m��3,�Ƕ���2

  newimg->pixels=malloc(sizeof(BYTE)*byte_per_pixel*width*height);//�����Ϥ��u�����Ŷ�(���c�����Х�����W�bmalloc�@��
																  //�j�p�� �C��pixel��byte�� * �e * �� �� byte 
  if(newimg->pixels==NULL) {//�T�wnewimg->pixels���Q���t�Ŷ�
    free(newimg);
    return NULL;
  }

  newimg->width=width;  //�s�U�Ϥ����e��
  newimg->height=height;//�s�U�Ϥ�������
  newimg->depth=depth;  //�s�Udepth
  return newimg;//����V�Ϥ������Ц^��
}


void disposeImage(ImageData *img)
{
  if(img->pixels!=NULL) //����img->pixels ���O�����m�٦^�h
	  free(img->pixels);
  free(img);//�A��img���O�����m�٦^�h
  //���Ǥ��i�f,���M�|memory leakage
  return;
}


int getPixel(ImageData *img,int x,int y,Pixel *pix)
{
  //Ū���Ϥ����@��Pixel(��(x,y)��
  int ret=1;
  int adr;  
  int dep,val;
  BYTE *pixels;

  if(img==NULL) return -1;//�S�Ϥ��N�^��-1
  if(img->pixels==NULL) return -1;//�Ϥ����Spixel�^��-1

  if(x<0) {//X�p��0X����0
    x=0;
    ret=0;
  }
  if(x >= img->width ) {//X �W�L�Ϥ����e��X�ܬ��Ϥ��e��-1
    x=img->width -1;
    ret=0;
  }
  if(y<0) {//Y�p��0Y����0
    y=0;
    ret=0;
  }
  if(y >= img->height ) {//X �W�L�Ϥ�������X�ܬ��Ϫ���-1
    y=img->height -1;
    ret=0;
  }
  dep=img->depth;//�o��depth (���m�ΦǶ�
  adr=x + y*img->width;//�p���ePixel��memory address
  pixels=(BYTE*)img->pixels;//���o��Pixel
  if(dep==8) {//�Ƕ����p  
    val=pixels[adr];
    pix->r=val;//���⪺��
    pix->g=val;//��⪺��
    pix->b=val;//�Ŧ⪺��
  }
  else if(dep==24) {//���m���p
    pixels+=(adr*3);
    pix->r=(*pixels);//���⪺��
    pixels++;
    pix->g=(*pixels);//��⪺��
    pixels++;
    pix->b=(*pixels);//�Ŧ⪺��
  }
  else {
    return -1;
  }
  return ret;
}


int correctValue(int val,int max)
{
  if(val<0) return 0;
  if(val>max) return max;
  return val;
}


int setPixel(ImageData *img,int x,int y,Pixel *pix)
{//�NPixel���Ϥ��̭�
  int adr;  
  int dep;
  int val;
  BYTE *pixels;

  if(img==NULL) return -1;//�S�Ϥ��N�^��-1
  if(img->pixels==NULL) return -1;//�Ϥ����Spixel�^��-1

  if(x<0 || x >= img->width || y<0 || y >= img->height ) {
    return 0;
  }
  dep=img->depth;
  adr=x + y*img->width;
  pixels=(BYTE*)img->pixels;
  if(dep==8) {//�Ƕ����p
    pixels[adr]=correctValue(pix->r,PIXELMAX);
  }
  else if(dep==24) {//���m���p

    pixels+=(adr*3);
    (*pixels)=correctValue(pix->r,PIXELMAX);
    pixels++;
    (*pixels)=correctValue(pix->g,PIXELMAX);
    pixels++;
    (*pixels)=correctValue(pix->b,PIXELMAX);
  }
  else {
    return -1;
  }
  return 1;
}

int fwriteWORD(WORD val,FILE *fp)
{	//�g�J�@�� WORD�]2BYTE�^�j�p����
	int i,c;

	c=val;
	for(i=0;i<2;i++) {
		fputc(c%256,fp);//256�O�]��1��BYTE = 8��BIT = 256
		c/=256;
	}
	return TRUE;
}


int fwriteDWORD(DWORD val,FILE *fp)
{	//�g�J�@�� DOUBLE WORD�]4BYTE�^�j�p����
	int i,c;

	c=val;
	for(i=0;i<4;i++) {
		fputc(c%256,fp);//256�O�]��1��BYTE = 8��BIT = 256
		c/=256;
	}
	return TRUE;
}


int freadWORD(WORD *res,FILE *fp)
{	//Ū�J�@�� WORD�]2BYTE�^�j�p����
	int i,c;
	int val[2];

	for(i=0;i<2;i++) {
		c=fgetc(fp);
		if(c==EOF) return FALSE;
		val[i]=c;
	}
	*res=val[1]*256+val[0];//256�O�]��1��BYTE = 8��BIT = 256
	return TRUE;
}


int  freadDWORD(DWORD *res,FILE *fp)
{	//Ū�J�@�� DOUBLE WORD�]2BYTE�^�j�p����
	int i,c;
	int val[4];
	DWORD tmp=0;

	for(i=0;i<4;i++) {
		c=fgetc(fp);
		if(c==EOF) return FALSE;
		val[i]=c;
	}
	tmp=0;
	for(i=3;i>=0;i--) {
		tmp*=256;//256�O�]��1��BYTE = 8��BIT = 256
		tmp+=val[i];
	}
	*res=tmp;
	return TRUE;
}


static BOOL	IsWinDIB(BITMAPINFOHEADER* pBIH)
{	//�T�{�o��Bmp(Bitmap��)�O���ODIB�榡 �O���ܦ^��TRUE �Ϥ��h��FALSE
	//�@�ӧ��㪺DIB�Ѩⳡ���զ��G�@��BITMAPINFO���c�M�@�Ӧs�x�����}�C���Ʋ�
	if (((BITMAPCOREHEADER*)pBIH)->bcSize == sizeof(BITMAPCOREHEADER)) {
		return FALSE;
	}
	return TRUE;
}


int countOfDIBColorEntries(int iBitCount)
{
	int	iColors;
	//�����C��
	switch (iBitCount) {
	case 1://���� �Ƕ�
		iColors	= 2;
		break;
	case 4://�m��
		iColors	= 16;
		break;
	case 8://���m
		iColors	= 256;
		break;
	default:
		iColors	= 0;
		break;
	}

	return iColors;
}


int getDIBxmax(int mx,int dep)
{
	switch(dep) {
	case 32:
		return mx*4;
	case 24:
		//return mx;
		return ((mx*3)+3)/4*4;
		break;
	case 16:
		return (mx+1)/2*2;
		break;
	case 8:
		return (mx+3)/4*4;
		break;
	case 4:
		return (((mx+1)/2)+3)/4*4;
		break;
	case 1:
		return (((mx+7)/8)+3)/4*4;
	}
	return mx;
}


int readBMPfile(char *fname,ImageData **img)
{	//Ū�i�@�i�Ϥ��s�bimg�o���и�
	int i,c;
	int errcode=0;
	BITMAPFILEHEADER BMPFile;
	int	fsize;
	BITMAPINFOHEADER BMPInfo;
	BITMAPCOREHEADER BMPCore;
	int	colors;
	int	colorTableSize;
	int	bitsSize;
	int	BISize;
	int x,y;
	int mx,my,depth;
	int pad;
	int mxb,myb;
	int isPM=FALSE;	
	FILE *fp;
	//DIB�����c
    WORD    HEAD_bfType;
    DWORD   HEAD_bfSize;
    WORD    HEAD_bfReserved1;
    WORD    HEAD_bfReserved2;
    DWORD   HEAD_bfOffBits;
    DWORD   INFO_bfSize;
    Pixel palet[MAXCOLORS];
    Pixel setcolor;

	if((fp=fopen(fname,"rb"))==NULL) {
		return -1;
	}


	if(!freadWORD(&HEAD_bfType,fp)) {//�}�lŪ�ɮ������p�G�SŪ��N��$ABORT
		errcode=-2;
		goto $ABORT;
	}
	if (HEAD_bfType != 0x4d42) {//0x4D42="BM" �OBITMAP��TYPE
		errcode=-10;				//�ΨӽT�w�O���O.bmp��
		goto $ABORT;
	}

	if(!freadDWORD(&HEAD_bfSize,fp)) {//Ū�ɮפj�p
		errcode=-10;
		goto $ABORT;
	}

	if(!freadWORD(&HEAD_bfReserved1,fp)) {
		errcode=-10;
		goto $ABORT;
	}

	if(!freadWORD(&HEAD_bfReserved2,fp)) {
		errcode=-10;
		goto $ABORT;
	}

	if(!freadDWORD(&HEAD_bfOffBits,fp)) {
		errcode=-10;
		goto $ABORT;
	}

	if(!freadDWORD(&INFO_bfSize,fp)) {
		errcode=-10;
		goto $ABORT;
	}


	if (INFO_bfSize == 40/*sizeof(BITMAPINFOHEADER)*/ || INFO_bfSize == 12/*sizeof(BITMAPCOREHEADER)*/) {
		BMPInfo.biSize =	INFO_bfSize;

		if(INFO_bfSize == sizeof(BITMAPCOREHEADER)) {
			WORD tmp;
			isPM =	TRUE;

			if(!freadWORD(&tmp,fp)) {
				errcode=-10;
				goto $ABORT;
			}
			BMPInfo.biWidth=tmp;

			if(!freadWORD(&tmp,fp)) {
				errcode=-10;
				goto $ABORT;
			}
			BMPInfo.biHeight=tmp;

			if(!freadWORD(&(BMPInfo.biPlanes),fp)) {
				errcode=-10;
				goto $ABORT;
			}

			if(!freadWORD(&(BMPInfo.biBitCount),fp)) {
				errcode=-10;
				goto $ABORT;
			}
		}
		else {		

			if(!freadDWORD(&(BMPInfo.biWidth),fp)) {
				errcode=-10;
				goto $ABORT;
			}

			if(!freadDWORD(&(BMPInfo.biHeight),fp)) {
				errcode=-10;
				goto $ABORT;
			}

			if(!freadWORD(&(BMPInfo.biPlanes),fp)) {
				errcode=-10;
				goto $ABORT;
			}

			if(!freadWORD(&(BMPInfo.biBitCount),fp)) {
				errcode=-10;
				goto $ABORT;
			}
		}

		if(!isPM) {

			if(!freadDWORD(&(BMPInfo.biCompression),fp)) {
				errcode=-10;
				goto $ABORT;
			}

			if(!freadDWORD(&(BMPInfo.biSizeImage),fp)) {
				errcode=-10;
				goto $ABORT;
			}

			if(!freadDWORD(&(BMPInfo.biXPelsPerMeter),fp)) {
				errcode=-10;
				goto $ABORT;
			}

			if(!freadDWORD(&(BMPInfo.biYPelsPerMeter),fp)) {
				errcode=-10;
				goto $ABORT;
			}

			if(!freadDWORD(&(BMPInfo.biClrUsed),fp)) {
				errcode=-10;
				goto $ABORT;
			}

			if(!freadDWORD(&(BMPInfo.biClrImportant),fp)) {
				errcode=-10;
				goto $ABORT;
			}
		}
	}
	else {
		errcode=-10;
		goto $ABORT;
	}
	mx=BMPInfo.biWidth;
	my=BMPInfo.biHeight;
	depth=BMPInfo.biBitCount;

	if(depth!=8 && depth!=24) {
		errcode=-3;
		goto $ABORT;
	}

	if(BMPInfo.biCompression!=BI_RGB) {
		errcode=-20;
		goto $ABORT;
	}

	if(BMPInfo.biClrUsed==0) {
		colors	= countOfDIBColorEntries(BMPInfo.biBitCount);
	}
	else {
		colors	= BMPInfo.biClrUsed;
	}


	if (!isPM)	{
		for(i=0;i<colors;i++) {

			c=fgetc(fp);
			if(c==EOF) {
				errcode=-10;
				goto $ABORT;
			}
			palet[i].b=c;

			c=fgetc(fp);
			if(c==EOF) {
				errcode=-10;
				goto $ABORT;
			}
			palet[i].g=c;

			c=fgetc(fp);
			if(c==EOF) {
				errcode=-10;
				goto $ABORT;
			}
			palet[i].r=c;

			c=fgetc(fp);
			if(c==EOF) {
				errcode=-10;
				goto $ABORT;
			}
		}
	} else {
		for(i=0;i<colors;i++) {
			c=fgetc(fp);
			if(c==EOF) {
				errcode=-10;
				goto $ABORT;
			}
			palet[i].b=c;
			c=fgetc(fp);
			if(c==EOF) {
				errcode=-10;
				goto $ABORT;
			}
			palet[i].g=c;
			c=fgetc(fp);
			if(c==EOF) {
				errcode=-10;
				goto $ABORT;
			}
			palet[i].r=c;
		}
	}

	*img=createImage(mx,my,24);
	mxb=getDIBxmax(mx,depth);
	pad=mxb-mx*depth/8;

	for(y=my-1;y>=0;y--) {
		for(x=0;x<mx;x++) {
			if(depth==8) {	
				c=fgetc(fp);
				if(c==EOF) {
					errcode=-20;
					goto $ABORT;
				}
				setcolor.r=palet[c].r;
				setcolor.g=palet[c].g;
				setcolor.b=palet[c].b;
			}
			else if(depth==24) {
				c=fgetc(fp);
				if(c==EOF) {
					errcode=-20;
					goto $ABORT;
				}
				setcolor.b=c;
				c=fgetc(fp);
				if(c==EOF) {
					errcode=-20;
					goto $ABORT;
				}
				setcolor.g=c;
				c=fgetc(fp);
				if(c==EOF) {
					errcode=-20;
					goto $ABORT;
				}
				setcolor.r=c;
			}
			setPixel(*img,x,y,&setcolor);
		}

		for(i=0;i<pad;i++) {
			c=fgetc(fp);
			if(c==EOF) {
				errcode=-20;
				goto $ABORT;
			}
		}
	}
$ABORT:	
	fclose(fp);
	return errcode;
}


int writeBMPfile(char *fname,ImageData *img)
{	//�N�@�i�w�g��ҥHpixel�s�n��ImageData���c�u����X���@�i�Ϥ�
	FILE *fp;
	BITMAPFILEHEADER bfn;
	int w,h,rw;
	int mxb,pad;
	int depth;	
	int pbyte;	
	int palsize;	
	int x,y,i;
	int saveloop,saverest;
	int	iBytes;
	unsigned int wsize;
	Pixel pix;

	w=img->width;	//���o�Ϥ��e
	h=img->height;	//���o�Ϥ���
	depth=img->depth; //�����Ϥ���depeth(���m�Ƕ�...


	if(depth!=24) {//
		//errcode=-3;
		goto $abort1;
	}

	if(depth==24) {
		pbyte=1;
	}
	else {
		pbyte=depth/8;
	}
	if(depth>=24) {
		palsize=0;
	}
	else {
		palsize=256;
	}
	rw=getDIBxmax(w,depth);

	bfn.bfType=0x4d42;	//'BM'
	bfn.bfSize=14+40+/*sizeof(BITMAPFILEHEADER) +sizeof(BITMAPINFOHEADER) +*/
			   palsize*4/*sizeof(RGBQUAD)*/ +
			   rw*h*pbyte;
	bfn.bfReserved1=0;
	bfn.bfReserved2=0;
	bfn.bfOffBits=14+40/*sizeof(BITMAPFILEHEADER) +sizeof(BITMAPINFOHEADER)*/ +
			      palsize*4/*sizeof(RGBQUAD)*/;

	if((fp=fopen(fname,"wb"))==NULL) {
		goto $abort1;
	}
	//�g�JDIB�����c
	fwriteWORD(bfn.bfType,fp);
	fwriteDWORD(bfn.bfSize,fp);
	fwriteWORD(bfn.bfReserved1,fp);
	fwriteWORD(bfn.bfReserved2,fp);
	fwriteDWORD(bfn.bfOffBits,fp);
	fwriteDWORD(40/*sizeof(BITMAPINFOHEADER)*/,fp);	//biSize
	fwriteDWORD(w,fp);		// biWidth
	fwriteDWORD(h,fp);		// biHeight
	fwriteWORD(1,fp);		// biPlanes
	fwriteWORD(depth,fp);	// biBitCount
	fwriteDWORD(BI_RGB,fp);	// biCompression
	fwriteDWORD(0,fp);	// biSizeImage
	fwriteDWORD(300,fp);	// biXPelsPerMeter
	fwriteDWORD(300,fp);	// biYPelsPerMeter
	fwriteDWORD(0,fp);		// biClrUsed
	fwriteDWORD(0,fp);		// biClrImportant


	pad=rw-w*depth/8;

	for(y=h-1;y>=0;y--) {
		for(x=0;x<w;x++) {
			getPixel(img,x,y,&pix);
			fputc(pix.b,fp);
			fputc(pix.g,fp);
			fputc(pix.r,fp);
		}

		for(i=0;i<pad;i++) {
			fputc(0,fp);
		}
	}
	return 0;
$abort1:
	return 1;
$abort2:
	fclose(fp);
	return 1;
}

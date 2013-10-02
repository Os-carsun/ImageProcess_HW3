#include "image.h"
#include <math.h>
#include <stdlib.h>



#define PAI   3.1415926

#define WSIZE 27


void setMaskValue(double *G,int msize,int i,int j,double val)
{
	G[i*msize+j]=val;
}

void makeGauss(double scale,int msize,double *Gx,double *Gy,double *Gxy,double *Gxx,double *Gyy)
{
    int    i,j,k,m,skip;
    int    x0,y0,value ;
    int    vmin,vmax,range,val;
    double vx,vy,vxy,vxx,vyy,r0,exponent,dvalue;

    x0 = y0 = WSIZE/2;
    for(i=0; i<msize; i++){
    	for(j=0; j < msize; j++){
    	    r0 = (double)((i-y0)*(i-y0)+(j-x0)*(j-x0));

		    exponent = 1.0/(2.0*PAI*scale) * exp(-r0/(2.0*scale));
		    vx  = (-(double)(j-x0)/scale) * exponent;
		    vy  = (-(double)(i-y0)/scale) * exponent;
		    vxy = ((double)((i-y0)*(j-x0))/(scale*scale)) * exponent;
		    vxx = ((double)(j-x0)*(j-x0)/(scale*scale)-1.0/scale) * exponent;
		    vyy = ((double)(i-y0)*(i-y0)/(scale*scale)-1.0/scale) * exponent;

		    setMaskValue(Gx,msize,i,j,vx);
		    setMaskValue(Gy,msize,i,j,vy);
		    setMaskValue(Gxy,msize,i,j,vxy);
		    setMaskValue(Gxx,msize,i,j,vxx);
		    setMaskValue(Gyy,msize,i,j,vyy);
    	}
    }
}

int absi(int i)
{
	if(i<0) return -i;
	return i;
}

int* getFilter(int i)
{
	int fil1[9]={
		 1, 1, 1,
		 1,-2, 1,
		-1,-1,-1};
	int fil2[9]={
		 1, 1, 1,
		 1,-2,-1,
		 1,-1,-1};
	int fil3[9]={
		 1, 1,-1,
		 1,-2,-1,
		 1, 1,-1};
	int fil4[9]={
		 1,-1,-1,
		 1,-2,-1,
		 1, 1, 1};
	int fil5[9]={
		-1,-1,-1,
		 1,-2, 1,
		 1, 1, 1};
	int fil6[9]={
		-1,-1, 1,
		-1,-2, 1,
		 1, 1, 1};
	int fil7[9]={
		-1, 1, 1,
		-1,-2, 1,
		-1, 1, 1};
	int fil8[9]={
		 1, 1, 1,
		-1,-2, 1,
		-1,-1, 1};

	switch(i) {
	case 0:
		return fil1;
	case 1:
		return fil2;
	case 2:
		return fil3;
	case 3:
		return fil4;
	case 4:
		return fil5;
	case 5:
		return fil6;
	case 6:
		return fil7;
	case 7:
		return fil8;
	}
}


int effect_1(ImageData *img,ImageData *outimg)
{

	int fil1[9]={
	 0, 0, 0,
	 1, 0,-1,
	 0, 0, 0};
	int fil2[9]={
	 0, 1, 0,
	 0, 0, 0,
	 0,-1, 0};

	int x,y;
	int i;
	int val;
	int xx,yy;
	int rrx,ggx,bbx;
	int rry,ggy,bby;
	Pixel col;
	int sadr;
	int x1,y1,x2,y2;

	x1=0;
	y1=0;
	x2=img->width-1;
	y2=img->height-1;

	for(y=y1;y<=y2;y++) {
		for(x=x1;x<=x2;x++) {
			rrx=ggx=bbx=0;
			rry=ggy=bby=0;
			sadr=0;
			for(yy=0;yy<3;yy++) {
				for(xx=0;xx<3;xx++) {
					val = getPixel(img,x+xx-1,y+yy-1,&col);

					rrx+= col.r*fil1[sadr];
					ggx+= col.g*fil1[sadr];
					bbx+= col.b*fil1[sadr];
					rry+= col.r*fil2[sadr];
					ggy+= col.g*fil2[sadr];
					bby+= col.b*fil2[sadr];
					sadr++;
				}
			}
			col.r=(int)(sqrt((double)(rrx*rrx+rry*rry))/2.0);
			col.g=(int)(sqrt((double)(ggx*ggx+ggy*ggy))/2.0);
			col.b=(int)(sqrt((double)(bbx*bbx+bby*bby))/2.0);
			setPixel(outimg,x,y,&col);	
		}
	}
	return TRUE;
}

int effect_2(ImageData *img,ImageData *outimg)
{

	int fil1[9]={
	 1, 0,-1,
	 2, 0,-2,
	 1, 0,-1};
	int fil2[9]={
	 1, 2, 1,
	 0, 0, 0,
	-1,-2,-1};

	int x,y;
	int i;
	int val;
	int xx,yy;
	int rrx,ggx,bbx;
	int rry,ggy,bby;
	Pixel col;
	int sadr;
	int x1,y1,x2,y2;

	x1=0;
	y1=0;
	x2=img->width-1;
	y2=img->height-1;

	for(y=y1;y<=y2;y++) {
		for(x=x1;x<=x2;x++) {
			rrx=ggx=bbx=0;
			rry=ggy=bby=0;
			sadr=0;
			for(yy=0;yy<3;yy++) {
				for(xx=0;xx<3;xx++) {
					val = getPixel(img,x+xx-1,y+yy-1,&col);

					rrx+= col.r*fil1[sadr];
					ggx+= col.g*fil1[sadr];
					bbx+= col.b*fil1[sadr];
					rry+= col.r*fil2[sadr];
					ggy+= col.g*fil2[sadr];
					bby+= col.b*fil2[sadr];
					sadr++;
				}
			}
			col.r=(int)(sqrt((double)(rrx*rrx+rry*rry))/8.0);
			col.g=(int)(sqrt((double)(ggx*ggx+ggy*ggy))/8.0);
			col.b=(int)(sqrt((double)(bbx*bbx+bby*bby))/8.0);
			setPixel(outimg,x,y,&col);	
		}
	}
	return TRUE;
}

int effect_3(ImageData *img,ImageData *outimg)
{
	int robert1[9]={
	 0, 0, 0,
	 0, 1, 0,
	 0, 0,-1};
	int robert2[9]={
	 0, 0, 0,
	 0, 0, 1,
	 0,-1, 0};

	int x,y;
	int i;
	int val;
	int xx,yy;
	int rrx,ggx,bbx;
	int rry,ggy,bby;
	int rr,gg,bb;
	Pixel col;
	int sadr;
	int x1,y1,x2,y2;

	x1=0;
	y1=0;
	x2=img->width-1;
	y2=img->height-1;

	for(y=y1;y<=y2;y++) {
		for(x=x1;x<=x2;x++) {
			rrx=ggx=bbx=0;
			rry=ggy=bby=0;
			sadr=0;
			for(yy=0;yy<3;yy++) {
				for(xx=0;xx<3;xx++) {
					val = getPixel(img,x+xx-1,y+yy-1,&col);

					rrx+= col.r*robert1[sadr];
					ggx+= col.g*robert1[sadr];
					bbx+= col.b*robert1[sadr];
					rry+= col.r*robert2[sadr];
					ggy+= col.g*robert2[sadr];
					bby+= col.b*robert2[sadr];
					sadr++;
				}
			}
			rrx=absi(rrx);
			ggx=absi(ggx);
			bbx=absi(bbx);
			rry=absi(rry);
			ggy=absi(ggy);
			bby=absi(bby);
			if(rrx>rry) rr=rrx;
			else rr=rry;
			if(ggx>ggy) gg=ggx;
			else gg=ggy;
			if(bbx>bby) bb=bbx;
			else bb=bby;
			col.r = rr;
			col.g = gg;
			col.b = bb;
			setPixel(outimg,x,y,&col);	
		}
	}
	return TRUE;
}

int effect_4(ImageData *img,ImageData *outimg,int som)
{

	int fil_x[11*11]={
	-14, 	-16, 	-17, 	-19, 	-20, 	0, 	20, 	19, 	17, 	16,	14,
	-16, 	-18, 	-20, 	-22, 	-24, 	0, 	24, 	22, 	20, 	18,	16,
	-17, 	-20, 	-24, 	-28, 	-32, 	0, 	32, 	28, 	24, 	20,	17,
	-19, 	-22, 	-28, 	-35, 	-45, 	0, 	45, 	35, 	28, 	22,	19,
	-20, 	-24, 	-32, 	-45, 	-71, 	0, 	71, 	45, 	32, 	24,	20,
	-20, 	-25, 	-33, 	-50, 	-100, 	0, 	100,	50, 	33, 	25,	20,
	-20, 	-24, 	-32, 	-45, 	-71, 	0, 	71, 	45, 	32, 	24,	20,
	-19, 	-22, 	-28, 	-35, 	-45, 	0, 	45, 	35, 	28, 	22,	19,
	-17, 	-20, 	-24, 	-28, 	-32, 	0, 	32, 	28, 	24, 	20,	17,
	-16, 	-18, 	-20, 	-22, 	-24, 	0, 	24, 	22, 	20, 	18,	16,
	-14, 	-16, 	-17, 	-19, 	-20, 	0, 	20, 	19, 	17, 	16,	14
	};

	int fil_y[11*11]={
	-14, 	-16, 	-17, 	-19, 	-20,  -20, 	-20, 	-19, 	-17, 	-16,	-14,
	-16, 	-18, 	-20, 	-22, 	-24,  -25, 	-24, 	-22, 	-20, 	-18,	-16,
	-17, 	-20, 	-24, 	-28, 	-32,  -33, 	-32, 	-28, 	-24, 	-20,	-17,
	-19, 	-22, 	-28, 	-35, 	-45,  -50, 	-45, 	-35, 	-28, 	-22,	-19,
	-20, 	-24, 	-32, 	-45, 	-71, -100, 	-71, 	-45, 	-32, 	-24,	-20,
	  0, 	  0, 	  0, 	  0, 	  0, 	0, 	0,		  0, 	  0, 	  0,	  0,
	20, 	24, 	32, 	45, 	71,   100, 	71, 	45, 	 32, 	 24,	20,
	19, 	22, 	28, 	35, 	45,    50, 	45, 	35, 	 28, 	 22,	19,
	17, 	20, 	24, 	28, 	32,    33, 	32, 	28, 	 24, 	 20,	17,
	16, 	18, 	20, 	22, 	24,    25, 	24, 	22, 	 20, 	 18,	16,
	14, 	16, 	17, 	19, 	20,    20, 	20, 	19, 	 17, 	 16,	14
	};

	int x,y;
	int i;
	int val;
	int xx,yy;
	int rrx,ggx,bbx;
	int rry,ggy,bby;
	int rr,gg,bb;
	Pixel col;
	int sadr;
	int x1,y1,x2,y2;
	int area,sum;

	x1=0;
	y1=0;
	x2=img->width-1;
	y2=img->height-1;

	area=som*2+1;
	sum=0;
	for(yy=0;yy<area;yy++) {
		sadr=(yy+(5-som))*11+(5-som);
		for(xx=0;xx<area;xx++) {
			sum+=absi(fil_x[sadr]);
			sadr++;
		}
	}
	for(y=y1;y<=y2;y++) {
		for(x=x1;x<=x2;x++) {
			rrx=ggx=bbx=0;
			rry=ggy=bby=0;
			for(yy=0;yy<area;yy++) {
				sadr=(yy+(5-som))*11+(5-som);
				for(xx=0;xx<area;xx++) {
					val = getPixel(img,x+xx-1,y+yy-1,&col);

					rrx+= (int)col.r*fil_x[sadr];
					ggx+= (int)col.g*fil_x[sadr];
					bbx+= (int)col.b*fil_x[sadr];
					rry+= (int)col.r*fil_y[sadr];
					ggy+= (int)col.g*fil_y[sadr];
					bby+= (int)col.b*fil_y[sadr];
					sadr++;
				}
			}
			rr=(int)(sqrt(rrx*rrx+rry*rry)/(double)sum);
			gg=(int)(sqrt(ggx*ggx+ggy*ggy)/(double)sum);
			bb=(int)(sqrt(bbx*bbx+bby*bby)/(double)sum);
			col.r = rr;
			col.g = gg;
			col.b = bb;
			setPixel(outimg,x,y,&col);	
		}
	}
	return TRUE;
}

int effect_5(ImageData *img,ImageData *outimg)
{



	int i;
	int val;
	int x,y;
	int xx,yy;
	int hh;
	int du,dd,dl,dr;
	int endn;
	int c1,c2;
	int rrx,ggx,bbx;
	int rr,gg,bb;
	Pixel col;
	int *fil;
	int sadr;
	int rate,vh;
	int maxR,maxG,maxB;
	int x1,y1,x2,y2;

	x1=0;
	y1=0;
	x2=img->width-1;
	y2=img->height-1;

	for(y=y1;y<=y2;y++) {
		for(x=x1;x<=x2;x++) {
			maxR=maxG=maxB= -10000;
			for(i=0;i<8;i++) {
				fil=getFilter(i);
				rrx=ggx=bbx=0;
				sadr=0;
				for(yy=0;yy<3;yy++) {
					for(xx=0;xx<3;xx++) {
						val = getPixel(img,x+xx-1,y+yy-1,&col);	
						rrx+= col.r*fil[sadr];
						ggx+= col.g*fil[sadr];
						bbx+= col.b*fil[sadr];
						sadr++;
					}
				}
				if(maxR<rrx) maxR=rrx;
				if(maxG<ggx) maxG=ggx;
				if(maxB<bbx) maxB=bbx;
			}
			maxR/=10;
			maxG/=10;
			maxB/=10;
			col.r = maxR;
			col.g = maxG;
			col.b = maxB;
			setPixel(outimg,x,y,&col);	
		}
	}
	return TRUE;
}

int effect_6(ImageData *img,ImageData *outimg,double scale,double th,int msize)
{
	int x,y;
	int i;
	int val;
	int xx,yy;
	int wt,ht;
	int xo,yo;
	int rrx,ggx,bbx;
	int rry,ggy,bby;
	int rr,gg,bb,gray;
	Pixel col;
	int sadr,adr;
	int x1,y1,x2,y2;

	double Lx,Ly,Lxy,Lxx,Lyy;
	double *lap,*edge;
	double *Gx,*Gy,*Gxy,*Gxx,*Gyy;

	Gx=(double*)malloc(sizeof(double)*msize*msize);
	Gy=(double*)malloc(sizeof(double)*msize*msize);
	Gxy=(double*)malloc(sizeof(double)*msize*msize);
	Gxx=(double*)malloc(sizeof(double)*msize*msize);
	Gyy=(double*)malloc(sizeof(double)*msize*msize);

	makeGauss(scale,msize,Gx,Gy,Gxy,Gxx,Gyy);

	x1=0;
	y1=0;
	x2=img->width-1;
	y2=img->height-1;


	xo = yo = WSIZE/2;

	wt=x2-x1+1;
	ht=y2-y1+1;
	lap = (double *)malloc(wt*ht*sizeof(double));
	edge = (double *)malloc(wt*ht*sizeof(double));

	for(y=y1;y<=y2;y++) {
		for(x=x1;x<=x2;x++) {
			Lx=Ly=Lxy=Lxx=Lyy=0.0;
			for(yy=0;yy<msize;yy++) {
				for(xx=0;xx<msize;xx++) {
					val = getPixel(img,x+xx-xo,y+yy-yo,&col);

					rr=col.r;
					gg=col.g;
					bb=col.b;
					gray=(bb*28+77*rr+gg*151)/256;

					Lx  += (double)gray * Gx[yy*msize+xx];
					Ly  += (double)gray * Gy[yy*msize+xx];
					Lxy += (double)gray * Gxy[yy*msize+xx];
					Lxx += (double)gray * Gxx[yy*msize+xx];
					Lyy += (double)gray * Gyy[yy*msize+xx];

				}
			}
			lap[x+y*wt] =(Lx*Lx*Lxx+2.0*Lxy*Lx*Ly+Ly*Ly*Lyy);
			edge[x+y*wt] =sqrt(Lx*Lx+Ly*Ly);
		}
	}
	for(y=y1;y<=y2-1;y++) {
		for(x=x1;x<=x2-1;x++) {
			adr=x+y*wt;
			col.r=0;
			col.g=0;
			col.b=0;

			if(edge[adr]>=th && (lap[adr]*lap[adr+1]<0.0 || lap[adr]*lap[adr+wt]<0.0) ) {
				col.r=255;
				col.g=255;
				col.b=255;
			}
			setPixel(outimg,x,y,&col);	
		}
	}
	free(lap);
	free(edge);

	return TRUE;

}

main(int ac,char *av[])
{
	ImageData *img,*outimg;
	int res;
	int p1;
	double p2,p3;
	int selector;
	char fileName[200];
	char outFileName[200];

	printf("input file name:\t");
	scanf("%s",fileName);

	printf("output file name:\t");
	scanf("%s",outFileName);
	res=readBMPfile(fileName,&img);

	if(res<0) {
		printf("not this file \n");
		return -1;
	}

	outimg=createImage(img->width,img->height,24);
	
	printf("select effect\n");
	scanf("%d",&selector);
	
	switch(selector){
		case 1:
			effect_1(img,outimg);
			break;
		case 2:
			effect_2(img,outimg);
			break;
		case 3:
			effect_3(img,outimg);
			break;
		case 4:
			printf("som?\n");
			scanf("%d",&p1);
			effect_4(img,outimg,p1);
			break;
		case 5:
			effect_5(img,outimg);
			break;
		case 6:
			printf("scale\n");
			scanf("%lf",&p2);
			printf("th\n");
			scanf("%lf",&p3);
			printf("msize\n");
			scanf("%d",&p1);
			effect_6(img,outimg,p2,p3,p1);
			break;

	}

	writeBMPfile(outFileName,outimg);
	disposeImage(img);
	disposeImage(outimg);

}




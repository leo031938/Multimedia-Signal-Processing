#include <stdio.h>
#include <stdlib.h>
#include <math.h>
float PI=3.1415927;

/*construct a structure of BMP header*/
typedef struct Bmpheader{
    unsigned short identifier; // 0x0000
    unsigned int filesize; // 0x0002
    unsigned short reserved; // 0x0006
    unsigned short reserved2;
    unsigned int bitmap_dataoffset; // 0x000A
    unsigned int bitmap_headersize; // 0x000E
    unsigned int width; // 0x0012
    unsigned int height; // 0x0016
    unsigned short planes; // 0x001A
    unsigned short bits_perpixel; // 0x001C
    unsigned int compression; // 0x001E
    unsigned int bitmap_datasize; // 0x0022
    unsigned int hresolution; // 0x0026
    unsigned int vresolution; // 0x002A
    unsigned int usedcolors; // 0x002E
    unsigned int importantcolors; // 0x0032
    unsigned int palette; // 0x0036
} Bitmap;

/*construct a structure of RGB*/
typedef struct RGB{
    int R;
    int G;
    int B;
} ImgRGB;
typedef struct rl{
    int R[64];
    int L[64];
}  ImgRL;

ImgRL RL1[25715], RL2[25715], RL3[25715] ;
Bitmap readheader(FILE* fp);
ImgRGB** malloc_2D(int row, int col);
void InputData(FILE* fp,ImgRGB **array,int H,int W);
void output_bmp(ImgRGB **RGB,FILE* outfile,Bitmap bmpheader);
void RGBtoYCbCr(float *R,float *G,float *B,float *Y,float *Cb,float *Cr);
void FDCT(float *in_data,float *out_data);
void IDCT(float *in_data,float *out_data);
void quantization(float *DCT, float *quantized, int n);
void dequantization(float*quantized, float *DCT, int n);
void YCbCrtoRGB(float *Y,float*Cb,float *Cr,float *R,float *G,float *B);
void zigzag(float *quantized, int *zigzaged);
void dezigzag(int *zigzaged, float *dezigzag);
void RLE(int *zigzaged, ImgRL *RL, int picnum);
void deRLE(ImgRL *RL, int *deRL, int picnum);
int main(int argc,char *argv[]){
//    FILE *fp=fopen(argv[1],"rb");
//    FILE *fw=fopen(argv[2],"wb");
    FILE *fp=fopen("input.bmp","rb");
    FILE *fw=fopen("output.bmp","wb");
    //read header
    Bitmap bmpheader=readheader(fp);
    ImgRGB **Data_RGB=malloc_2D(bmpheader.height, bmpheader.width);
    
    //read data
    InputData(fp,Data_RGB,bmpheader.height,bmpheader.width);
    printf("%d\n%d\n",bmpheader.height,bmpheader.width);
    /*put other else function here*/
    
    
    int countu=bmpheader.width>>3,countv=bmpheader.height>>3,h=8,w=8,u=0,v=0,x=0,y=0;
    int picnum=0;
    float Y=0, Cb=0, Cr=0;
    float tempY[8][8]={0},tempCb[8][8]={0}, tempCr[8][8]={0};
    float DCTY[8][8]={0}, DCTCb[8][8]={0}, DCTCr[8][8]={0};
    float Q1[8][8]={0},Q2[8][8]={0},Q3[8][8]={0};
    float DQ1[8][8]={0},DQ2[8][8]={0},DQ3[8][8]={0};
    float IDCTY[8][8]={0},IDCTCb[8][8]={0}, IDCTCr[8][8]={0};
    float buf1[8][8]={0}, buf2[8][8]={0}, buf3[8][8]={0};
    int zigzag1[8][8]={0}, zigzag2[8][8]={0}, zigzag3[8][8]={0};
    float dezigzag1[8][8]={0}, dezigzag2[8][8]={0}, dezigzag3[8][8]={0};
    int deRLE1[8][8]={0}, deRLE2[8][8]={0}, deRLE3[8][8]={0};
    int frontDC1=0, frontDC2=0, frontDC3=0;
    int temp1=0, temp2=0, temp3=0;
    int tempp1=0, tempp2=0, tempp3=0;
    
    //8*8 blocking
    printf("%d\t%d\n",countu,countv);
    
    for(u=0;u<=countu;u++){
        if(u==countu){
            w=bmpheader.width%8;
        }
        else w = 8;
        for(v=0;v<=countv;v++){
            if(v==countv){
                h=bmpheader.height%8;
            }
            else h = 8;
            for( x=0;x<w;x++){
                for( y=0;y<h;y++){
                    //get data
                    float R=Data_RGB[x+8*u][y+8*v].R;
                    float G=Data_RGB[x+8*u][y+8*v].G;
                    float B=Data_RGB[x+8*u][y+8*v].B;
                    //printf("%f\t%f\t%f\n",R, G, B);
                    RGBtoYCbCr(&R,&G,&B,&Y,&Cb,&Cr);// RGB->YCbCr
                    tempY[x][y]=Y-128;
                    tempCb[x][y]=Cb-128;
                    tempCr[x][y]=Cr-128;
                    //printf("%f\t%f\t%f\t%f\t%f\t%f\n", Y,Cb,Cr , tempY[x][y], tempCb[x][y],tempCr[x][y]);
                }
            }
            
            //DCT
            FDCT(&tempY[0][0], &DCTY[0][0]);
            FDCT(&tempCb[0][0], &DCTCb[0][0]);
            FDCT(&tempCr[0][0], &DCTCr[0][0]);
            
            /*int m=0,n=0;
            for(m=0;m<8;m++){
                for(n=0;n<8;n++){
                    printf("%f\t%f\t%f\n", DCTY[m][n], DCTCb[m][n], DCTCr[m][n]);
                }
            }
    printf("\n");*/
    
            //Quantization
            quantization(&DCTY[0][0], &Q1[0][0],0);
            quantization(&DCTCb[0][0], &Q2[0][0],1);
            quantization(&DCTCr[0][0], &Q3[0][0],1);
            
            /*int m=0,n=0;
             for(m=0;m<8;m++){
                 for(n=0;n<8;n++){
                     printf("%f\t%f\t%f\n", Q1[m][n], Q2[m][n], Q3[m][n]);
                }
             }
             printf("\n");*/
            
            //zigzag
            zigzag(&Q1[0][0], &zigzag1[0][0]);
            zigzag(&Q2[0][0], &zigzag2[0][0]);
            zigzag(&Q3[0][0], &zigzag3[0][0]);
            
            /*int m=0,n=0;
            for(m=0;m<8;m++){
                for(n=0;n<8;n++){
                    printf("%.3f\t", Q1[m][n]);
                }
                printf("\n");
            }*/
        
           /*int m=0,n=0;
            for(m=0;m<8;m++){
                for(n=0;n<8;n++){
                    printf("%d\t", zigzag1[m][n]);
                }
                printf("\n");
            }
            printf("eeeeee\n")*/
            
            
            //RLE
            RLE(&zigzag1[0][0], &RL1[0], picnum);
            RLE(&zigzag2[0][0], &RL2[0], picnum);
            RLE(&zigzag3[0][0], &RL3[0], picnum);
            
            /*for(m=0;m<64;m++){
                    printf("%d\t%d\n", RL1[picnum].R[m],  RL1[picnum].L[m]);
                }
                printf("\n");
            printf("eeeeee\n");
             */
            //printf("%d\t", zigzag1[0][0]);
            
            //DC
            //recode DC value
            temp1=zigzag1[0][0];
            temp2=zigzag2[0][0];
            temp3=zigzag3[0][0];
            
            //Get the first value
            zigzag1[0][0]=zigzag1[0][0]-frontDC1;
            zigzag2[0][0]=zigzag2[0][0]-frontDC2;
            zigzag3[0][0]=zigzag3[0][0]-frontDC3;
            
            //printf("%d\t",zigzag1[0][0]);
            
            //Record front  DC value
            frontDC1=temp1;
            frontDC2=temp2;
            frontDC3=temp3;

            //deDC
            //record diff value
            zigzag1[0][0]=zigzag1[0][0]+tempp1;
            zigzag2[0][0]=zigzag2[0][0]+tempp2;
            zigzag3[0][0]=zigzag3[0][0]+tempp3;
            //printf("%d\t%d\n", tempp1, zigzag1[0][0]);
            tempp1=zigzag1[0][0];
            tempp2=zigzag2[0][0];
            tempp3=zigzag3[0][0];
            
            
            deRLE1[0][0]=zigzag1[0][0];
            deRLE2[0][0]=zigzag2[0][0];
            deRLE3[0][0]=zigzag3[0][0];
            //deRLE
            deRLE(&RL1[0], &deRLE1[0][0], picnum);
            deRLE(&RL2[0], &deRLE2[0][0], picnum);
            deRLE(&RL3[0], &deRLE3[0][0], picnum);
            /*int m=0;
            for(m=0;m<8;m++){
                for(n=0;n<8;n++){
                    printf("%d\t", deRLE1[m][n]);
                }
                printf("\n");
            }
            printf("ffffffffff\n");*/
            
            picnum++;
            //printf("%d\n",picnum);
            //dezigzag
            dezigzag(&deRLE1[0][0], &dezigzag1[0][0]);
            dezigzag(&deRLE2[0][0], &dezigzag2[0][0]);
            dezigzag(&deRLE3[0][0], &dezigzag3[0][0]);
            //Inverse Quantization
            dequantization(&dezigzag1[0][0], &DQ1[0][0],0);
            dequantization(&dezigzag2[0][0], &DQ2[0][0],1);
            dequantization(&dezigzag3[0][0], &DQ3[0][0],1);
            /*
            for(m=0;m<8;m++){
                for(n=0;n<8;n++){
                    printf("%f\t%f\t%f\n", DQ1[m][n], DQ2[m][n], DQ3[m][n]);
                }
            }*/
            //Inverse DCT
            IDCT(&DQ1[0][0], &IDCTY[0][0]);
            IDCT(&DQ2[0][0], &IDCTCb[0][0]);
            IDCT(&DQ3[0][0], &IDCTCr[0][0]);
            /*
            for(m=0;m<8;m++){
                for(n=0;n<8;n++){
                    printf("%f\t%f\t%f\n", IDCTY[m][n], IDCTCb[m][n], IDCTCr[m][n]);
                }
            }*/
            // YCbCr->RGB
            YCbCrtoRGB(&IDCTY[0][0], &IDCTCb[0][0], &IDCTCr[0][0], &buf1[0][0],&buf2[0][0],&buf3[0][0]);
            //put data
            for( x=0;x<w;++x){
                for( y=0;y<h;++y){
                    Data_RGB[x+8*u][y+8*v].R=buf1[x][y];
                    Data_RGB[x+8*u][y+8*v].G=buf2[x][y];
                    Data_RGB[x+8*u][y+8*v].B=buf3[x][y];
                }
            }
            //printf("%d\t%d\n",u,v);
        }
    }

    
    //output bmp
    output_bmp(Data_RGB,fw,bmpheader);
    
    return 0;
}


/*read header*/
Bitmap readheader(FILE* fp){
    Bitmap x;
    fread(&x.identifier,sizeof(x.identifier),1,fp);
    fread(&x.filesize,sizeof(x.filesize),1,fp);
    fread(&x.reserved,sizeof(x.reserved),1,fp);
    fread(&x.reserved2,sizeof(x.reserved2),1,fp);
    fread(&x.bitmap_dataoffset,sizeof(x.bitmap_dataoffset),1,fp);
    fread(&x.bitmap_headersize,sizeof(x.bitmap_headersize),1,fp);
    fread(&x.width,sizeof(x.width),1,fp);
    fread(&x.height,sizeof(x.height),1,fp);
    fread(&x.planes,sizeof(x.planes),1,fp);
    fread(&x.bits_perpixel,sizeof(x.bits_perpixel),1,fp);
    fread(&x.compression,sizeof(x.compression),1,fp);
    fread(&x.bitmap_datasize,sizeof(x.bitmap_datasize),1,fp);
    fread(&x.hresolution,sizeof(x.hresolution),1,fp);
    fread(&x.vresolution,sizeof(x.vresolution),1,fp);
    fread(&x.usedcolors,sizeof(x.usedcolors),1,fp);
    fread(&x.importantcolors,sizeof(x.importantcolors),1,fp);
    return x;
}

/*input data without header into RGB structure*/
void InputData(FILE* fp,ImgRGB **array,int H, int W){
    int temp,i,j;
    for( i=0;i<H;i++){
        for( j=0;j<W;j++){
            temp=fgetc(fp);
            array[i][j].B=temp;
            temp=fgetc(fp);
            array[i][j].G=temp;
            temp=fgetc(fp);
            array[i][j].R=temp;
        }
    }
}

/* A function of making two dimensions memory locate array*/
ImgRGB** malloc_2D(int row, int col){
    ImgRGB **Array, *Data;
    int i;
    Array=(ImgRGB**)malloc(row*sizeof(ImgRGB *));
    Data=(ImgRGB*)malloc(row*col*sizeof(ImgRGB));
    for(i=0; i<row; i++,Data+=col){
        Array[i] = Data;
    }
    return Array;
}

/*output header and data*/
void output_bmp(ImgRGB **RGB,FILE* outfile,Bitmap bmpheader){
    int x,y;
    fwrite(&bmpheader.identifier, sizeof(short), 1 , outfile);
    fwrite(&bmpheader.filesize, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.reserved, sizeof(short), 1 , outfile);
    fwrite(&bmpheader.reserved2, sizeof(short), 1 , outfile);
    fwrite(&bmpheader.bitmap_dataoffset, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.bitmap_headersize, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.width, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.height, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.planes, sizeof(short), 1 , outfile);
    fwrite(&bmpheader.bits_perpixel, sizeof(short), 1 , outfile);
    fwrite(&bmpheader.compression, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.bitmap_datasize, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.hresolution, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.vresolution, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.usedcolors, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.importantcolors, sizeof(int), 1 , outfile);
    
    for ( x=0; x<bmpheader.height; x++){
        for( y=0; y<bmpheader.width; y++){
            fwrite(&RGB[x][y].B, sizeof(char),1,outfile);
            fwrite(&RGB[x][y].G, sizeof(char),1,outfile);
            fwrite(&RGB[x][y].R, sizeof(char),1,outfile);
        }
    }
}
void RGBtoYCbCr(float *R,float *G,float *B,float *Y,float *Cb,float *Cr){
    
    *Y = (0.299 * *R + 0.587 * *G + 0.114 * *B);
    *Cb = (-0.16874 * *R - 0.3313 * *G + 0.50000 * *B ) + 128;
    *Cr = (0.50000 * *R - 0.4187* *G - 0.0813* *B ) + 128;
    
}
void FDCT(float *in_data,float *out_data)
{
    float data, sqrt_2;
    int i, j, u, v, n=8;
    
    sqrt_2 = sqrt(2.0);
    for (u=0; u<n; u++)
        for (v=0; v<n; v++)
        {
            data = 0;
            for (i=0; i<n; i++){
                for (j=0; j<n; j++){
                    data += in_data[i*n+j] * cos((2*i+1)*u*PI/(2*n)) *
                    cos((2*j+1)*v*PI/(2*n));
                }
            }
            if (u == 0) data /= sqrt_2;
            if (v == 0) data /= sqrt_2;
            data = data*2/n;
            out_data[u*n+v] = data;
        }
}
void IDCT(float *in_data,float *out_data)
{
    
    float data, sub_data, sqrt_2;
    int i, j, u, v, n=8;
    
    sqrt_2 = sqrt(2.0);
    for (i=0; i<n; i++){
        for (j=0; j<n; j++)
        {
            data = 0;
            for (u=0; u<n; u++)
            {
                for (v=0; v<n; v++)
                {
                    sub_data = in_data[u*n+v] * cos((2*i+1)*u*PI/(2*n)) *
                    cos((2*j+1)*v*PI/(2*n));
                    if (u == 0) sub_data /= sqrt_2;
                    if (v == 0) sub_data /= sqrt_2;
                    
                    data += sub_data;
                }
            }
            data = data*0.25;
            out_data[i*n+j] = data;
        }
    }
}
int light[8][8]={
    16, 11, 10, 16, 24, 40, 51, 61,
    12, 12, 14, 19, 26, 58, 60, 55,
    14, 13, 16, 24, 40, 67, 69, 56,
    14, 17, 22, 29, 51, 87, 80, 62,
    18, 22, 37, 56, 68, 109, 103, 77,
    24, 35, 55, 64, 81, 104, 113, 92,
    49, 64, 78, 87, 103, 121, 120, 101,
    72, 92, 95, 98, 112, 100, 103, 99
};
int color[8][8]={
    17, 18, 24, 47, 99, 99, 99, 99,
    18, 21, 26, 66, 99, 99, 99, 99,
    24, 26, 56, 99, 99, 99, 99, 99,
    47, 66, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99
};
void quantization(float *DCT, float *quantized, int n){
    int u, v;
    for (u=0; u<8; u++){
        for (v=0; v<8; v++)
            if(n==0){
                quantized[u*8+v]=DCT[u*8+v]/light[u][v];
            }
            else if(n==1){
                quantized[u*8+v]=DCT[u*8+v]/color[u][v];
            }
    }
}
void dequantization(float *quantized, float *DCT, int n){
    int u, v;
    for (u=0; u<8; u++){
        for (v=0; v<8; v++){
            if(n==0){
                DCT[u*8+v]=quantized[u*8+v]*light[u][v];
            }
            else if(n==1){
                DCT[u*8+v]=quantized[u*8+v]*color[u][v];
            }
        }
    }
}
void YCbCrtoRGB(float *Y,float *Cb,float *Cr,float *R,float *G,float *B){
    int u, v;
    for (u=0; u<8; u++){
        for (v=0; v<8; v++){
            R[u+8*v]  = Y[u+8*v]  + 1.402*Cr[u+8*v] +128 ;
            G[u+8*v]  = Y[u+8*v]  - 0.34414*Cb[u+8*v]  - 0.71414*Cr[u+8*v] +128 ;
            B[u+8*v]  = Y[u+8*v]  + 1.772*Cb[u+8*v] +128 ;
        }
    }
}
int scan_order[64] = {
    0,  1,  8, 16,  9,  2,  3, 10,
    17, 24, 32, 25, 18, 11,  4,  5,
    12, 19, 26, 33, 40, 48, 41, 34,
    27, 20, 13,  6,  7, 14, 21, 28,
    35, 42, 49, 56, 57, 50, 43, 36,
    29, 22, 15, 23, 30, 37, 44, 51,
    58, 59, 52, 45, 38, 31, 39, 46,
    53, 60, 61, 54, 47, 55, 62, 63
};
int descan_order[64] = {
    0, 1, 5, 6, 14, 15, 27, 28,
    2, 4, 7, 13, 16, 26, 29, 42,
    3, 8, 12, 17, 25, 30, 41, 43,
    9, 11, 18, 24, 31, 40, 44, 53,
    10, 19, 23, 32, 39, 45, 52, 54,
    20, 22, 33, 38, 46, 41, 55, 60,
    21, 34, 37, 47, 50, 56, 59, 61,
    35, 36, 48, 49, 57, 58, 62, 63
};
void zigzag(float *quantized, int *zigzaged){
    int i=0;
    for(i=0; i<64; i++){
        zigzaged[i]=(int)((quantized[scan_order[i]]));
    }
}
void dezigzag(int *zigzaged, float *dezigzag){
    int i=0;
    for(i=0; i<64; i++){
        dezigzag[i]=(float)(zigzaged[descan_order[i]]);
    }
}
void RLE(int *zigzaged, ImgRL *RL, int picnum){
    int i=0, j=0, zerocount=0;
    
    for(i=1; i<64; i++){
        if(zigzaged[i]!=0){
            RL[picnum].R[j]=zerocount;
            RL[picnum].L[j]=zigzaged[i];
            zerocount=0;
            j++;
        }
        else if(zigzaged[i]==0){
            zerocount++;
            }
        }
    }


void deRLE(ImgRL *RL, int *deRL, int picnum){
    int i=0, j=0;
    for(i=1; i<64; i++){
        
        if(RL[picnum].R[j]==0){
            deRL[i]=RL[picnum].L[j];
            j++;
        }
        
        else if(RL[picnum].R[j]!=0){
            for(int k=0; k<RL[picnum].R[j]; k++){
                deRL[i]=0;
                i++;
            }
            deRL[i]=RL[picnum].L[j];
            j++;
        }
    }
}

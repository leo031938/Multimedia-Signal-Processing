#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define double_pi 6.283158//define the value which DFT transform will use it

struct Wave{
     char RIFF[4];
     int  all_size;
     char wave[4];
     char fmt_subchunk[4];
     int fmt_size;
     short fmt_type;
     short channels;
     int sample_rate;
     int byte_rate;
     short block_align;
     short bits_per_sample;
     char data[4];
     int data_size;
     short data8;
    };

int main(int argc,char*argv[])
{

    struct Wave wav;//build .wav Header file format
    FILE *input;
    FILE *output;
    input = fopen(argv[1],"rb");
    output = fopen(argv[2],"wb");

    //read Header file
    fread(&wav.RIFF,4,1,input);
    printf("%.4s\n",wav.RIFF);

    fread(&wav.all_size,4,1,input);
    printf("%d\n",wav.all_size);

    fread(&wav.wave,4,1,input);
    printf("%.4s\n",wav.wave);

    fread(&wav.fmt_subchunk,4,1,input);
    printf("%.4s\n",wav.fmt_subchunk);

    fread(&wav.fmt_size,4,1,input);
    printf("%d\n",wav.fmt_size);

    fread(&wav.fmt_type,2,1,input);
    printf("%d\n",wav.fmt_type);

    fread(&wav.channels,2,1,input);
    printf("%d\n",wav.channels);

    fread(&wav.sample_rate,4,1,input);
    printf("%d\n",wav.sample_rate);

    fread(&wav.byte_rate,4,1,input);
    printf("%d\n",wav.byte_rate);

    fread(&wav.block_align,2,1,input);
    printf("%d\n",wav.block_align);

    fread(&wav.bits_per_sample,2,1,input);
    printf("%d\n",wav.bits_per_sample);

    fread(&wav.data,4,1,input);
    printf("%.4s\n",wav.data);

    fread(&wav.data_size,4,1,input);
    printf("%d\n",wav.data_size);

    //read data
    float x[28000];//2 bytes per bytes so 56000/2 samples
    int i;
    for(i=0;i<28000;i++){
        fread(&wav.data8,2,1,input);
        //printf("%d\n",wav.data8);
        float xdata=wav.data8;
        /*long complement=0;

        if(xdata>32765){
            complement=xdata-65536;//the number bigger than 32765 is negative,so we need to solve that
        }
        else{
            complement=xdata;//keeping intact
        }*/
        float float_data = (float)(xdata/32767);//transform it between -1~1
        x[i]=float_data;
        //printf("%f",x[i]);
        }
    fclose(input);
    printf("  *  \n * * \n*****\n");

    //DFT tranform with LPF
    float Xre[28000];
    float Xim[28000];
    int k;
    for(i=0;i<28000;i++){
        Xre[i]=0;
        Xim[i]=0;
        for(k=0;k<28000;k++){
            Xre[i]+=x[k]*cos(i*k*double_pi/28000);
            Xim[i]-=x[k]*sin(i*k*double_pi/28000);
        }
            if(i>12250 && i<15750){
                Xre[i]=0;
                Xim[i]=0;
            }
            else{
                Xre[i]=Xre[i]/2;
                Xim[i]=Xim[i]/2;
            }
    }
    printf("IDFT\n");
    //IDFT tranform
    float Yim[28000];
    float Yre[28000];
    short A[28000];
    for(i=0;i<28000;i++){
        Yim[i]=0;
        Yre[i]=0;
        for(k=0;k<28000;k++){
            Yre[i] += Xre[k]*cos(i*k*double_pi/28000)-Xim[k]*sin(double_pi*i*k/28000);
            Yim[i] += Xim[k]*cos(i*k*double_pi/28000)+Xre[k]*sin(double_pi*i*k/28000);
            }
    }
    for(i=0;i<28000;i++){
        Yre[i] = Yre[i]/28000;
        Yim[i] = Yim[i]/28000;
        A[i]=(short)(sqrt(Yre[i]*Yre[i]+Yim[i]*Yim[i])*32768);
        if(x[i]<0){
            A[i] = -A[i];//translate the signal
        }
    }
    fwrite(&wav,44,1,output);//write in the header
    for(i=0;i<28000;i++){
            printf("%d  =  %f  %d\n",i,x[i],A[i]);
            fwrite(&A[i],2,1,output);
        }
    fclose(output);
    printf("0987759579\n");

    return 0;
}



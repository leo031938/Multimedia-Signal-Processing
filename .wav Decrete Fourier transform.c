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
     int data8;
    };

int main(int argc,char*argv[])
{

    struct Wave wav;//build .wav Header file format
    FILE *input;
    FILE *output_magnitude;
    FILE *output_phase;
    input = fopen(argv[1],"rb");
    output_magnitude= fopen(argv[2],"wb");
    output_phase= fopen(argv[3],"wb");

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
        unsigned long xdata=wav.data8;
        long complement=0;

        if(xdata>32765){
            complement=xdata-65536;//the number bigger than 32765 is negative,so we need to solve that
        }
        else{
            complement=xdata;//keeping intact
        }
        double float_data = (double)(complement/(double)32768);//transform it between -1~1
            x[i]=float_data;
        }
    fclose(input);

    //DFT tranform
    float Xre[28000];
    float Xim[28000];
    float P[28000];
    float D[28000];
    float A[28000];

    int k;

    for(i=0;i<28000;i++){
        Xre[i]=0;
        Xim[i]=0;
        for(k=0;k<28000;k++){
            Xre[i]+=x[k]*cos(i*k*double_pi/28000);
            Xim[i]-=x[k]*sin(i*k*double_pi/28000);
        }
        P[i]=Xre[i]*Xre[i]+Xim[i]*Xim[i];
        A[i]=20*log10(P[i]/20);//compute the amplitude
        D[i]=atan(Xim[i]/Xre[i]);//compute the phase
        fwrite(&A[i],sizeof(A[i]),1,output_magnitude);
        fwrite(&D[i],sizeof(D[i]),1,output_phase);
        }
    fclose(output_magnitude);
    fclose(output_phase);
    return 0;
}



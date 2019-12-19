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
    short datahaha;
};
int main()
{
    struct Wave in1,in2;
    FILE *input1,*input2;
    input1 = fopen("vowel-16k.wav","rb");
    input2 = fopen("vowel-8k.wav","rb");
    fread(&in1,sizeof(in1),1,input1);
    fread(&in2,sizeof(in2),1,input2);
    printf("%d\n",in1.data_size);
    printf("%d\n",in2.data_size);

    //read Header file
     /*fread(&wav.RIFF,4,1,input1);
     fread(&wav.all_size,4,1,input1);
     fread(&wav.wave,4,1,input);
     fread(&wav.fmt_subchunk,4,1,input);
     fread(&wav.fmt_size,4,1,input);
     fread(&wav.fmt_type,2,1,input);
     fread(&wav.channels,2,1,input);
     fread(&wav.sample_rate,4,1,input);
     fread(&wav.byte_rate,4,1,input);
     fread(&wav.block_align,2,1,input);
     fread(&wav.bits_per_sample,2,1,input);
     fread(&wav.data,4,1,input);
     fread(&wav.data_size,4,1,input);*/

    struct Wave wav1;
    wav1.RIFF[0]='R';
    wav1.RIFF[1]='I';
    wav1.RIFF[2]='F';
    wav1.RIFF[3]='F';
    wav1.all_size=32036;
    wav1.wave[0]='W';
    wav1.wave[1]='A';
    wav1.wave[2]='V';
    wav1.wave[3]='E';
    wav1.fmt_subchunk[0]='f';
    wav1.fmt_subchunk[1]='m';
    wav1.fmt_subchunk[2]='t';
    wav1.fmt_subchunk[3]=' ';
    wav1.fmt_size=16;
    wav1.fmt_type=1;
    wav1.channels=1;
    wav1.sample_rate=16000;
    wav1.byte_rate=32000;//SampleRate * NumChannels * BitsPerSample/8
    wav1.block_align=2;
    wav1.bits_per_sample=16;
    wav1.data[0]='d';
    wav1.data[1]='a';
    wav1.data[2]='t';
    wav1.data[3]='a';
    wav1.data_size=32000;

    FILE *fp0;
    FILE *fp1;
    FILE *fp2;
    FILE *fp3;

    FILE *fp4;
    FILE *fp5;
    FILE *fp6;
    FILE *fp7;

    fp0=fopen("cos_050Hz_16k.wav","wb");
    fp1=fopen("cos_200Hz_16k.wav","wb");
    fp2=fopen("cos_055Hz_16k.wav","wb");
    fp3=fopen("cos_220Hz_16k.wav","wb");

    fp4=fopen("cos_050Hz_8k.wav","wb");
    fp5=fopen("cos_200Hz_8k.wav","wb");
    fp6=fopen("cos_055Hz_8k.wav","wb");
    fp7=fopen("cos_220Hz_8k.wav","wb");


    fwrite(&wav1.RIFF,1,4,fp0);
    printf("%.4s\n",wav1.RIFF);
    fwrite(&wav1.all_size,4,1,fp0);
    printf("%d\n",wav1.all_size);
    fwrite(&wav1.wave,1,4,fp0);
    printf("%.4s\n",wav1.wave);

    fwrite(&wav1.fmt_subchunk,1,4,fp0);
    printf("%.4s\n",wav1.fmt_subchunk);
    fwrite(&wav1.fmt_size,4,1,fp0);
    printf("%d\n",wav1.fmt_size);
    fwrite(&wav1.fmt_type,2,1,fp0);
    printf("%d\n",wav1.fmt_type);
    fwrite(&wav1.channels,2,1,fp0);
    printf("%d\n",wav1.channels);
    fwrite(&wav1.sample_rate,4,1,fp0);
    printf("%d\n",wav1.sample_rate);
    fwrite(&wav1.byte_rate,4,1,fp0);
    printf("%d\n",wav1.byte_rate);
    fwrite(&wav1.block_align,2,1,fp0);
    printf("%d\n",wav1.block_align);
    fwrite(&wav1.bits_per_sample,2,1,fp0);
    printf("%d\n",wav1.bits_per_sample);
    fwrite(wav1.data,1,4,fp0);
    printf("%.4s\n",wav1.data);
    fwrite(&wav1.data_size,4,1,fp0);
    printf("%d\n",wav1.data_size);

    fwrite(&wav1.RIFF,1,4,fp1);
    printf("%.4s\n",wav1.RIFF);
    fwrite(&wav1.all_size,4,1,fp1);
    printf("%d\n",wav1.all_size);
    fwrite(&wav1.wave,1,4,fp1);
    printf("%.4s\n",wav1.wave);

    fwrite(&wav1.fmt_subchunk,1,4,fp1);
    printf("%.4s\n",wav1.fmt_subchunk);
    fwrite(&wav1.fmt_size,4,1,fp1);
    printf("%d\n",wav1.fmt_size);
    fwrite(&wav1.fmt_type,2,1,fp1);
    printf("%d\n",wav1.fmt_type);
    fwrite(&wav1.channels,2,1,fp1);
    printf("%d\n",wav1.channels);
    fwrite(&wav1.sample_rate,4,1,fp1);
    printf("%d\n",wav1.sample_rate);
    fwrite(&wav1.byte_rate,4,1,fp1);
    printf("%d\n",wav1.byte_rate);
    fwrite(&wav1.block_align,2,1,fp1);
    printf("%d\n",wav1.block_align);
    fwrite(&wav1.bits_per_sample,2,1,fp1);
    printf("%d\n",wav1.bits_per_sample);
    fwrite(wav1.data,1,4,fp1);
    printf("%.4s\n",wav1.data);
    fwrite(&wav1.data_size,4,1,fp1);
    printf("%d\n",wav1.data_size);

    fwrite(&wav1.RIFF,1,4,fp2);
    printf("%.4s\n",wav1.RIFF);
    fwrite(&wav1.all_size,4,1,fp2);
    printf("%d\n",wav1.all_size);
    fwrite(&wav1.wave,1,4,fp2);
    printf("%.4s\n",wav1.wave);

    fwrite(&wav1.fmt_subchunk,1,4,fp2);
    printf("%.4s\n",wav1.fmt_subchunk);
    fwrite(&wav1.fmt_size,4,1,fp2);
    printf("%d\n",wav1.fmt_size);
    fwrite(&wav1.fmt_type,2,1,fp2);
    printf("%d\n",wav1.fmt_type);
    fwrite(&wav1.channels,2,1,fp2);
    printf("%d\n",wav1.channels);
    fwrite(&wav1.sample_rate,4,1,fp2);
    printf("%d\n",wav1.sample_rate);
    fwrite(&wav1.byte_rate,4,1,fp2);
    printf("%d\n",wav1.byte_rate);
    fwrite(&wav1.block_align,2,1,fp2);
    printf("%d\n",wav1.block_align);
    fwrite(&wav1.bits_per_sample,2,1,fp2);
    printf("%d\n",wav1.bits_per_sample);
    fwrite(wav1.data,1,4,fp2);
    printf("%.4s\n",wav1.data);
    fwrite(&wav1.data_size,4,1,fp2);
    printf("%d\n",wav1.data_size);

    fwrite(&wav1.RIFF,1,4,fp3);
    printf("%.4s\n",wav1.RIFF);
    fwrite(&wav1.all_size,4,1,fp3);
    printf("%d\n",wav1.all_size);
    fwrite(&wav1.wave,1,4,fp3);
    printf("%.4s\n",wav1.wave);

    fwrite(&wav1.fmt_subchunk,1,4,fp3);
    printf("%.4s\n",wav1.fmt_subchunk);
    fwrite(&wav1.fmt_size,4,1,fp3);
    printf("%d\n",wav1.fmt_size);
    fwrite(&wav1.fmt_type,2,1,fp3);
    printf("%d\n",wav1.fmt_type);
    fwrite(&wav1.channels,2,1,fp3);
    printf("%d\n",wav1.channels);
    fwrite(&wav1.sample_rate,4,1,fp3);
    printf("%d\n",wav1.sample_rate);
    fwrite(&wav1.byte_rate,4,1,fp3);
    printf("%d\n",wav1.byte_rate);
    fwrite(&wav1.block_align,2,1,fp3);
    printf("%d\n",wav1.block_align);
    fwrite(&wav1.bits_per_sample,2,1,fp3);
    printf("%d\n",wav1.bits_per_sample);
    fwrite(wav1.data,1,4,fp3);
    printf("%.4s\n",wav1.data);
    fwrite(&wav1.data_size,4,1,fp3);
    printf("%d\n",wav1.data_size);


    wav1.all_size=16036;
    wav1.sample_rate=8000;
    wav1.byte_rate=16000;
    wav1.data_size=16000;


    fwrite(&wav1.RIFF,1,4,fp4);
    printf("%.4s\n",wav1.RIFF);
    fwrite(&wav1.all_size,4,1,fp4);
    printf("%d\n",wav1.all_size);
    fwrite(&wav1.wave,1,4,fp4);
    printf("%.4s\n",wav1.wave);

    fwrite(&wav1.fmt_subchunk,1,4,fp4);
    printf("%.4s\n",wav1.fmt_subchunk);
    fwrite(&wav1.fmt_size,4,1,fp4);
    printf("%d\n",wav1.fmt_size);
    fwrite(&wav1.fmt_type,2,1,fp4);
    printf("%d\n",wav1.fmt_type);
    fwrite(&wav1.channels,2,1,fp4);
    printf("%d\n",wav1.channels);
    fwrite(&wav1.sample_rate,4,1,fp4);
    printf("%d\n",wav1.sample_rate);
    fwrite(&wav1.byte_rate,4,1,fp4);
    printf("%d\n",wav1.byte_rate);
    fwrite(&wav1.block_align,2,1,fp4);
    printf("%d\n",wav1.block_align);
    fwrite(&wav1.bits_per_sample,2,1,fp4);
    printf("%d\n",wav1.bits_per_sample);
    fwrite(wav1.data,1,4,fp4);
    printf("%.4s\n",wav1.data);
    fwrite(&wav1.data_size,4,1,fp4);
    printf("%d\n",wav1.data_size);


    fwrite(&wav1.RIFF,1,4,fp5);
    printf("%.4s\n",wav1.RIFF);
    fwrite(&wav1.all_size,4,1,fp5);
    printf("%d\n",wav1.all_size);
    fwrite(&wav1.wave,1,4,fp5);
    printf("%.4s\n",wav1.wave);

    fwrite(&wav1.fmt_subchunk,1,4,fp5);
    printf("%.4s\n",wav1.fmt_subchunk);
    fwrite(&wav1.fmt_size,4,1,fp5);
    printf("%d\n",wav1.fmt_size);
    fwrite(&wav1.fmt_type,2,1,fp5);
    printf("%d\n",wav1.fmt_type);
    fwrite(&wav1.channels,2,1,fp5);
    printf("%d\n",wav1.channels);
    fwrite(&wav1.sample_rate,4,1,fp5);
    printf("%d\n",wav1.sample_rate);
    fwrite(&wav1.byte_rate,4,1,fp5);
    printf("%d\n",wav1.byte_rate);
    fwrite(&wav1.block_align,2,1,fp5);
    printf("%d\n",wav1.block_align);
    fwrite(&wav1.bits_per_sample,2,1,fp5);
    printf("%d\n",wav1.bits_per_sample);
    fwrite(wav1.data,1,4,fp5);
    printf("%.4s\n",wav1.data);
    fwrite(&wav1.data_size,4,1,fp5);
    printf("%d\n",wav1.data_size);

    fwrite(&wav1.RIFF,1,4,fp6);
    printf("%.4s\n",wav1.RIFF);
    fwrite(&wav1.all_size,4,1,fp6);
    printf("%d\n",wav1.all_size);
    fwrite(&wav1.wave,1,4,fp6);
    printf("%.4s\n",wav1.wave);

    fwrite(&wav1.fmt_subchunk,1,4,fp6);
    printf("%.4s\n",wav1.fmt_subchunk);
    fwrite(&wav1.fmt_size,4,1,fp6);
    printf("%d\n",wav1.fmt_size);
    fwrite(&wav1.fmt_type,2,1,fp6);
    printf("%d\n",wav1.fmt_type);
    fwrite(&wav1.channels,2,1,fp6);
    printf("%d\n",wav1.channels);
    fwrite(&wav1.sample_rate,4,1,fp6);
    printf("%d\n",wav1.sample_rate);
    fwrite(&wav1.byte_rate,4,1,fp6);
    printf("%d\n",wav1.byte_rate);
    fwrite(&wav1.block_align,2,1,fp6);
    printf("%d\n",wav1.block_align);
    fwrite(&wav1.bits_per_sample,2,1,fp6);
    printf("%d\n",wav1.bits_per_sample);
    fwrite(wav1.data,1,4,fp6);
    printf("%.4s\n",wav1.data);
    fwrite(&wav1.data_size,4,1,fp6);
    printf("%d\n",wav1.data_size);

    fwrite(&wav1.RIFF,1,4,fp7);
    printf("%.4s\n",wav1.RIFF);
    fwrite(&wav1.all_size,4,1,fp7);
    printf("%d\n",wav1.all_size);
    fwrite(&wav1.wave,1,4,fp7);
    printf("%.4s\n",wav1.wave);

    fwrite(&wav1.fmt_subchunk,1,4,fp7);
    printf("%.4s\n",wav1.fmt_subchunk);
    fwrite(&wav1.fmt_size,4,1,fp7);
    printf("%d\n",wav1.fmt_size);
    fwrite(&wav1.fmt_type,2,1,fp7);
    printf("%d\n",wav1.fmt_type);
    fwrite(&wav1.channels,2,1,fp7);
    printf("%d\n",wav1.channels);
    fwrite(&wav1.sample_rate,4,1,fp7);
    printf("%d\n",wav1.sample_rate);
    fwrite(&wav1.byte_rate,4,1,fp7);
    printf("%d\n",wav1.byte_rate);
    fwrite(&wav1.block_align,2,1,fp7);
    printf("%d\n",wav1.block_align);
    fwrite(&wav1.bits_per_sample,2,1,fp7);
    printf("%d\n",wav1.bits_per_sample);
    fwrite(wav1.data,1,4,fp7);
    printf("%.4s\n",wav1.data);
    fwrite(&wav1.data_size,4,1,fp7);
    printf("%d\n",wav1.data_size);


    //x(t)=10000cos(2pift)w(t) f=50,200,55,220
    short s0[16000],s1[16000],s2[16000],s3[16000],s4[8000],s5[8000],s6[8000],s7[8000];
    int N=16000;
    int i;
    for(i=0;i<N;i++){
        s0[i]=10000*cos(double_pi*50*i/N);
        s1[i]=10000*cos(double_pi*200*i/N);
        s2[i]=10000*cos(double_pi*55*i/N);
        s3[i]=10000*cos(double_pi*220*i/N);
        fwrite(&s0[i],2,1,fp0);
        fwrite(&s1[i],2,1,fp1);
        fwrite(&s2[i],2,1,fp2);
        fwrite(&s3[i],2,1,fp3);
    }

    int M=8000;
    for(i=0;i<M;i++){
        s4[i]=10000*cos(double_pi*50*i/M);
        s5[i]=10000*cos(double_pi*200*i/M);
        s6[i]=10000*cos(double_pi*55*i/M);
        s7[i]=10000*cos(double_pi*220*i/M);
        fwrite(&s4[i],2,1,fp4);
        fwrite(&s5[i],2,1,fp5);
        fwrite(&s6[i],2,1,fp6);
        fwrite(&s7[i],2,1,fp7);
    }


    FILE *fp01;
    FILE *fp02;
    FILE *fp03;
    FILE *fp04;
    FILE *fp11;
    FILE *fp12;
    FILE *fp13;
    FILE *fp14;
    FILE *fp21;
    FILE *fp22;
    FILE *fp23;
    FILE *fp24;
    FILE *fp31;
    FILE *fp32;
    FILE *fp33;
    FILE *fp34;
    FILE *fp41;
    FILE *fp42;
    FILE *fp43;
    FILE *fp44;
    FILE *fp51;
    FILE *fp52;
    FILE *fp53;
    FILE *fp54;
    FILE *fp61;
    FILE *fp62;
    FILE *fp63;
    FILE *fp64;
    FILE *fp71;
    FILE *fp72;
    FILE *fp73;
    FILE *fp74;
    FILE *fp81;
    FILE *fp82;
    FILE *fp83;
    FILE *fp84;
    FILE *fp91;
    FILE *fp92;
    FILE *fp93;
    FILE *fp94;


    fp01=fopen("cos_050Hz_16k.{Set1}.txt","wb");
    fp02=fopen("cos_050Hz_16k.{Set2}.txt","wb");
    fp03=fopen("cos_050Hz_16k.{Set3}.txt","wb");
    fp04=fopen("cos_050Hz_16k.{Set4}.txt","wb");

    fp11=fopen("cos_200Hz_16k.{Set1}.txt","wb");
    fp12=fopen("cos_200Hz_16k.{Set2}.txt","wb");
    fp13=fopen("cos_200Hz_16k.{Set3}.txt","wb");
    fp14=fopen("cos_200Hz_16k.{Set4}.txt","wb");

    fp21=fopen("cos_055Hz-16k.{Set1}.txt","wb");
    fp22=fopen("cos_055Hz-16k.{Set2}.txt","wb");
    fp23=fopen("cos_055Hz-16k.{Set3}.txt","wb");
    fp24=fopen("cos_055Hz-16k.{Set4}.txt","wb");

    fp31=fopen("cos_220Hz_16k.{Set1}.txt","wb");
    fp32=fopen("cos_220Hz_16k.{Set2}.txt","wb");
    fp33=fopen("cos_220Hz_16k.{Set3}.txt","wb");
    fp34=fopen("cos_220Hz_16k.{Set4}.txt","wb");

    fp41=fopen("cos_050Hz_8k.{Set1}.txt","wb");
    fp42=fopen("cos_050Hz_8k.{Set2}.txt","wb");
    fp43=fopen("cos_050Hz_8k.{Set3}.txt","wb");
    fp44=fopen("cos_050Hz_8k.{Set4}.txt","wb");

    fp51=fopen("cos_200Hz_8k.{Set1}.txt","wb");
    fp52=fopen("cos_200Hz_8k.{Set2}.txt","wb");
    fp53=fopen("cos_200Hz_8k.{Set3}.txt","wb");
    fp54=fopen("cos_200Hz_8k.{Set4}.txt","wb");

    fp61=fopen("cos_055Hz_8k.{Set1}.txt","wb");
    fp62=fopen("cos_055Hz_8k.{Set2}.txt","wb");
    fp63=fopen("cos_055Hz_8k.{Set3}.txt","wb");
    fp64=fopen("cos_055Hz_8k.{Set4}.txt","wb");

    fp71=fopen("cos_220Hz_8k.{Set1}.txt","wb");
    fp72=fopen("cos_220Hz_8k.{Set2}.txt","wb");
    fp73=fopen("cos_220Hz_8k.{Set3}.txt","wb");
    fp74=fopen("cos_220Hz_8k.{Set4}.txt","wb");

    fp81=fopen("vowel-16k.{Set1}.txt","wb");
    fp82=fopen("vowel-16k.{Set2}.txt","wb");
    fp83=fopen("vowel-16k.{Set3}.txt","wb");
    fp84=fopen("vowel-16k.{Set4}.txt","wb");

    fp91=fopen("vowel-8k.{Set1}.txt","wb");
    fp92=fopen("vowel-8k.{Set2}.txt","wb");
    fp93=fopen("vowel-8k.{Set3}.txt","wb");
    fp94=fopen("vowel-8k.{Set4}.txt","wb");







//setcos_050Hz-16k.txt
    /*fread(&wav1.RIFF,4,1,fp0);
     fread(&wav1.all_size,4,1,fp0);
     fread(&wav1.wave,4,1,fp0);
     fread(&wav1.fmt_subchunk,4,1,fp0);
     fread(&wav1.fmt_size,4,1,fp0);
     fread(&wav1.fmt_type,2,1,fp0);
     fread(&wav1.channels,2,1,fp0);
     fread(&wav1.sample_rate,4,1,fp0);
     fread(&wav1.byte_rate,4,1,fp0);
     fread(&wav1.block_align,2,1,fp0);
     fread(&wav1.bits_per_sample,2,1,fp0);
     fread(&wav1.data,4,1,fp0);
     fread(&wav1.data_size,4,1,fp0);
     for(i=0;i<N;i++){
     fread(&fp0,2,1,data[i]);
     }*/

//set1
    int D,F,j,W,k,y;
    short x[512];
    double Xre[512],Xim[512],A[512];
    F=200;//frame
    W=16000*0.005;//window size
    D=16000*0.008;//DFT window size
    M=16000*0.005;//frame size
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)
            x[j]=s0[i*M+j]*1;
            else x[j]=0;
        }
        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp01,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp01,"\n");
    }

//set2
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s0[i*M+j]*(0.54-0.46*cos(double_pi*j/(M-1)));
            else x[j]=0;
        }
        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp02,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp02,"\n");
    }

//set3
    F=100;//frame
    W=16000*0.02;//window size
    D=16000*0.032;//DFT window size
    M=16000*0.01;//frame size

    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s0[i*M+j]*1;
            else x[j]=0;
        }

        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp03,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp03,"\n");
    }

//set4
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s0[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else x[j]=0;
        }

        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp04,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp04,"\n");
    }
//set cos_200Hz-16k.txt
//set1
    F=200;//frame
    W=16000*0.005;//window size
    D=16000*0.008;//DFT window size
    M=16000*0.005;//frame size
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(M>j)x[j]=s1[i*M+j]*1;
            else x[j]=0;
        }

    for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp11,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp11,"\n");
    }

//set2
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s1[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else x[j]=0;
        }
    for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp12,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp12,"\n");
    }

//set3
    F=100;//frame
    W=16000*0.02;//window size
    D=16000*0.032;//DFT window size
    M=16000*0.01;//frame size

    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s1[i*M+j]*1;
            else x[j]=0;
        }

    for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp13,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp13,"\n");
    }
//set4

    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s1[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else x[j]=0;
        }

        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp14,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp14,"\n");
    }
//setcos_55Hz-16k.txt
//set1

    F=200;//frame
    W=16000*0.005;//window size
    D=16000*0.008;//DFT window size
    M=16000*0.005;//frame size
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s2[i*M+j]*1;
            else x[j]=0;
        }

        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp21,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp21,"\n");
    }
//set2

    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s2[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else x[j]=0;
        }

        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp22,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp22,"\n");
    }

//set3
    F=100;//frame
    W=16000*0.02;//window size
    D=16000*0.032;//DFT window size
    M=16000*0.01;//frame size

    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s2[i*M+j]*1;
            else x[j]=0;
        }

        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp23,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp23,"\n");
    }

//set4
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s2[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else x[j]=0;
        }

        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp24,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp24,"\n");
    }
//setcos_220Hz-16k.txt
//set1
    F=200;//frame
    W=16000*0.005;//window size
    D=16000*0.008;//DFT window size
    M=16000*0.005;//frame size
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s3[i*M+j]*1;
            else x[j]=0;
        }

        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp31,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp31,"\n");
    }

//set2
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s3[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else x[j]=0;
        }

    for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp32,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp32,"\n");
    }

//set3
    F=100;//frame
    W=16000*0.02;//window size
    D=16000*0.032;//DFT window size
    M=16000*0.01;//frame size

    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s3[i*M+j]*1;
            else x[j]=0;
        }
    for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp33,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp33,"\n");
    }
//set4

    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s3[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else x[j]=0;
        }

    for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp34,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp34,"\n");
    }

//cos_050Hz-8k.txt
//set1

    int P=8000;
    F=200;//frame
    W=P*0.005;//window size
    D=P*0.008;//DFT window size
    M=P*0.005;//frame size
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s4[i*M+j]*1;
            else x[j]=0;
        }
        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp41,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp41,"\n");
    }
    //set2
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s4[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else x[j]=0;
        }

        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp42,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp42,"\n");
    }

//set3
    F=100;//frame
    W=P*0.02;//window size
    D=P*0.032;//DFT window size
    M=P*0.01;//frame size

    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(M>j)x[j]=s4[i*M+j]*1;
            else x[j]=0;
        }

        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp43,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp43,"\n");
    }
//set4

    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s4[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else x[j]=0;
        }
        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp44,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp44,"\n");
    }
//setcos_200Hz-8k.txt
//set1

    F=200;//frame
    W=P*0.005;//window size
    D=P*0.008;//DFT window size
    M=P*0.005;//frame size
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s5[i*M+j]*1;
            else x[j]=0;
        }

        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp51,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp51,"\n");
    }

//set2
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s5[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else x[j]=0;
        }

   for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp52,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp52,"\n");
    }

//set3
    F=100;//frame
    W=P*0.02;//window size
    D=P*0.032;//DFT window size
    M=P*0.01;//frame size

    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(M>j)x[j]=s5[i*M+j]*1;
            else x[j]=0;
        }

        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp53,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp53,"\n");
    }

//set4
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s5[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else x[j]=0;
        }

        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp54,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp54,"\n");
    }

//setcos_55Hz-8k.txt
//set1
    F=200;//frame
    W=P*0.005;//window size
    D=P*0.008;//DFT window size
    M=P*0.005;//frame size
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s6[i*M+j]*1;
            else x[j]=0;
        }

    for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp61,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp61,"\n");
    }
//set2

    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s6[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else x[j]=0;
        }
        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp62,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp62,"\n");
    }

//set3
    F=100;//frame
    W=P*0.02;//window size
    D=P*0.032;//DFT window size
    M=P*0.01;//frame size

    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(M>j)x[j]=s6[i*M+j]*1;
            else x[j]=0;
        }
for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp63,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp63,"\n");
    }

//set4
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s6[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else x[j]=0;
        }
        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp64,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp64,"\n");
    }
//setcos_220Hz-8k.txt
//set1

    F=200;//frame
    W=P*0.005;//window size
    D=P*0.008;//DFT window size
    M=P*0.005;//frame size
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s7[i*M+j]*1;
            else x[j]=0;
        }

    for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp71,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp71,"\n");
    }

//set2

    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s7[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else x[j]=0;
        }

        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp72,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp72,"\n");
    }

//set3
    F=100;//frame
    W=P*0.02;//window size
    D=P*0.032;//DFT window size
    M=P*0.01;//frame size

    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(M>j)x[j]=s7[i*M+j]*1;
            else x[j]=0;
        }

        for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp73,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp73,"\n");
    }

//set4
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)x[j]=s7[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else x[j]=0;
        }

     for(k=0;k<D;k++){
            Xre[k]=0;
            Xim[k]=0;
            for(y=0;y<D;y++){
                Xre[k]+=x[y]*cos(k*y*double_pi/D);
                Xim[k]-=x[y]*sin(k*y*double_pi/D);
            }
        A[k]=20*log10(sqrt(Xre[k]*Xre[k]+Xim[k]*Xim[k]));
        fprintf(fp74,"%f ",A[k]);
        printf("%f\n",A[k]);
        }
    fprintf(fp74,"\n");
    }


//vowel-16k.{Set1}.txt
//set1

    short z[1824];
    double Xrez[1824],Ximz[1824],Az[1824];
    short data[57000];
    for(i=0;i<57000;i++){
        fread(&in1.datahaha,2,1,input1);
        data[i]=in1.datahaha;
    }

    F=712;//frame
    W=16000*0.005;//window size
    D=16000*0.008;//DFT window size
    M=16000*0.005;//frame size
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)
            z[j]=data[i*M+j]*1;
            else z[j]=0;

        }
        for(k=0;k<D;k++){
            Xrez[k]=0;
            Ximz[k]=0;
            for(y=0;y<D;y++){
                Xrez[k]+=z[y]*cos(k*y*double_pi/D);
                Ximz[k]-=z[y]*sin(k*y*double_pi/D);
            }
        Az[k]=20*log10(sqrt(Xrez[k]*Xrez[k]+Ximz[k]*Ximz[k]));
        fprintf(fp81,"%f ",Az[k]);
        }
    fprintf(fp81,"\n");
    }
//set2
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)
            z[j]=data[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else z[j]=0;
        }

        for(k=0;k<D;k++){
            Xrez[k]=0;
            Ximz[k]=0;
            for(y=0;y<D;y++){
                Xrez[k]+=z[y]*cos(k*y*double_pi/D);
                Ximz[k]-=z[y]*sin(k*y*double_pi/D);
            }
        Az[k]=20*log10(sqrt(Xrez[k]*Xrez[k]+Ximz[k]*Ximz[k]));
        fprintf(fp82,"%f ",Az[k]);
        }
    fprintf(fp82,"\n");
    }

//set3

    F=356;//frame
    W=16000*0.02;//window size
    D=16000*0.032;//DFT window size
    M=16000*0.01;//frame size
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)
            z[j]=data[i*M+j]*1;
            else z[j]=0;

        }
        for(k=0;k<D;k++){
            Xrez[k]=0;
            Ximz[k]=0;
            for(y=0;y<D;y++){
                Xrez[k]+=z[y]*cos(k*y*double_pi/D);
                Ximz[k]-=z[y]*sin(k*y*double_pi/D);
            }
        Az[k]=20*log10(sqrt(Xrez[k]*Xrez[k]+Ximz[k]*Ximz[k]));
        fprintf(fp83,"%f ",Az[k]);
        }
    fprintf(fp83,"\n");
    }
//set4

    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)
            z[j]=data[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else z[j]=0;
        }
        for(k=0;k<D;k++){
            Xrez[k]=0;
            Ximz[k]=0;
            for(y=0;y<D;y++){
                Xrez[k]+=z[y]*cos(k*y*double_pi/D);
                Ximz[k]-=z[y]*sin(k*y*double_pi/D);
            }
        Az[k]=20*log10(sqrt(Xrez[k]*Xrez[k]+Ximz[k]*Ximz[k]));
        fprintf(fp84,"%f ",Az[k]);
        }
    fprintf(fp84,"\n");
    }

//vowel-8k.{Set1}.txt
//set1

    for(i=0;i<57000;i++){
        fread(&in2.datahaha,2,1,input2);
        data[i]=in2.datahaha;
    }
    F=712;//frame
    W=8000*0.005;//window size
    D=8000*0.008;//DFT window size
    M=8000*0.005;//frame size
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)
            z[j]=data[i*M+j]*1;
            else z[j]=0;
        }
        for(k=0;k<D;k++){
            Xrez[k]=0;
            Ximz[k]=0;
            for(y=0;y<D;y++){
                Xrez[k]+=z[y]*cos(k*y*double_pi/D);
                Ximz[k]-=z[y]*sin(k*y*double_pi/D);
            }
        Az[k]=20*log10(sqrt(Xrez[k]*Xrez[k]+Ximz[k]*Ximz[k]));
        fprintf(fp91,"%f ",Az[k]);
        }
    fprintf(fp91,"\n");
    }

//set2
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)
            z[j]=data[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));
            else z[j]=0;
        }
        for(k=0;k<D;k++){
            Xrez[k]=0;
            Ximz[k]=0;
            for(y=0;y<D;y++){
                Xrez[k]+=z[y]*cos(k*y*double_pi/D);
                Ximz[k]-=z[y]*sin(k*y*double_pi/D);
            }
        Az[k]=20*log10(sqrt(Xrez[k]*Xrez[k]+Ximz[k]*Ximz[k]));
        fprintf(fp92,"%f ",Az[k]);
        }
    fprintf(fp92,"\n");
    }
//set3

    F=356;//frame
    W=8000*0.02;//window size
    D=8000*0.032;//DFT window size
    M=8000*0.01;//frame size
    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)
            z[j]=data[i*M+j]*1;
            else z[j]=0;

        }
        for(k=0;k<D;k++){
            Xrez[k]=0;
            Ximz[k]=0;
            for(y=0;y<D;y++){
                Xrez[k]+=z[y]*cos(k*y*double_pi/D);
                Ximz[k]-=z[y]*sin(k*y*double_pi/D);
            }
        Az[k]=20*log10(sqrt(Xrez[k]*Xrez[k]+Ximz[k]*Ximz[k]));
        fprintf(fp93,"%f ",Az[k]);
        }
    fprintf(fp93,"\n");
    }
//set4

    for(i=0;i<F;i++){
        for(j=0;j<D;j++){
            if(W>j)
            z[j]=data[i*M+j]*(0.54-0.46*cos(double_pi*j/(W-1)));

            else z[j]=0;
        }
        for(k=0;k<D;k++){
            Xrez[k]=0;
            Ximz[k]=0;
            for(y=0;y<D;y++){
                Xrez[k]+=z[y]*cos(k*y*double_pi/D);
                Ximz[k]-=z[y]*sin(k*y*double_pi/D);
            }
        Az[k]=20*log10(sqrt(Xrez[k]*Xrez[k]+Ximz[k]*Ximz[k]));
        fprintf(fp94,"%f ",Az[k]);
        }
    fprintf(fp94,"\n");
    }
}


package com.fftdemo;

/**
 * Created by Administrator on 2018/7/10.
 */
import android.util.Log;

import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

public class fft {
    public String TAG="fft";
    public static double PAI =3.14159265358979;
    public static int N_FFT=1024;
    public static int  M_of_N_FFT = 10;
    public static int Npart2_of_N_FFT = 512;
    public static int Npart4_of_N_FFT =256;
    public static double [] SIN_TABLE_of_N_FFT=new double[256 + 1];
    public complex_of_N_FFT [] data_of_N_FFT=new complex_of_N_FFT[1024];


    public void CREATE_SIN_TABLE()
    {
        int i = 0;
        for (i = 0; i <= Npart4_of_N_FFT; i++)
        {
            SIN_TABLE_of_N_FFT[i] = sin(PAI * i / Npart2_of_N_FFT);//SIN_TABLE[i] = sin(PI2*i/N);
        }
    }
    public double Sin_find(double x)
    {
        int i = (int)(N_FFT * x);
        i = i >> 1;
        if (i > Npart4_of_N_FFT)//注意：i已经转化为0~N之间的整数了！
        {
            //不会超过N/2
            i = Npart2_of_N_FFT - i;//i = i - 2*(i-Npart4);
        }
        return SIN_TABLE_of_N_FFT[i];
    }
    public double Cos_find(double x)
    {
        int i = (int)(N_FFT * x);
        i = i >> 1;
        if (i < Npart4_of_N_FFT)//注意：i已经转化为0~N之间的整数了！
        {
            //不会超过N/2
            //i = Npart4 - i;
            return SIN_TABLE_of_N_FFT[Npart4_of_N_FFT - i];
        }
        else//i>Npart4 && i<N/2
        {
            //i = i - Npart4;
            return -SIN_TABLE_of_N_FFT[i - Npart4_of_N_FFT];
        }
    }
    public void ChangeSeat(complex_of_N_FFT DataInput[])
    {
        int nextValue, nextM, i, k, j = 0;
        complex_of_N_FFT temp;
        nextValue = N_FFT / 2;                  //变址运算，即把自然顺序变成倒位序，采用雷德算法
        nextM = N_FFT - 1;
        for (i = 0; i < nextM; i++)
        {
            if (i < j)                                        //如果i<j,即进行变址
            {
                temp = DataInput[j];
                DataInput[j] = DataInput[i];
                DataInput[i] = temp;
            }
            k = nextValue;                //求j的下一个倒位序
            while (k <= j)                            //如果k<=j,表示j的最高位为1
            {
                j = j - k;                                      //把最高位变成0
                k = k / 2;                                      //k/2，比较次高位，依次类推，逐个比较，直到某个位为0
            }
            j = j + k;                                          //把0改为1
        }
    }

    short freadint16_live(byte data0,byte data1)
    {

        byte []buffer=new byte[2];
        buffer[0]=data0;
        buffer[1]=data1;
        int []buffer1=new int[2];
        for (int i = 0;i < 2;i++)
        {
            buffer1[i] = buffer[i];
            if (buffer1[i] < 0)
            {
                buffer1[i] = 256 + buffer1[i];
            }
        }
        int j;
        j = buffer1[1];
        j = (j << 8) + buffer1[0];
        return (short) j;
    }

    void FFT()
    {
        int L = 0, B = 0, J = 0, K = 0;
        int step = 0, KB = 0;
        //ElemType P=0;
        double angle;
        complex_of_N_FFT W=new complex_of_N_FFT();
        complex_of_N_FFT Temp_XX=new complex_of_N_FFT();
        ChangeSeat(data_of_N_FFT);//变址
        //CREATE_SIN_TABLE();
        for (L = 1; L <= M_of_N_FFT; L++)
        {
            step = 1 << L;//2^L
            B = step >> 1;//B=2^(L-1)
            for (J = 0; J < B; J++)
            {
                //P = (1<<(M-L))*J;//P=2^(M-L) *J
                angle = (double)J / B;                        //这里还可以优化
                W.imag = -Sin_find(angle);         //用C++该函数课声明为inline
                W.real = Cos_find(angle);         //用C++该函数课声明为inline
                //W.real =  cos(angle*PI);
                //W.imag = -sin(angle*PI);
                for (K = J; K < N_FFT; K = K + step)
                {
                    KB = K + B;
                    //Temp_XX = XX_complex(data[KB],W);
                    //用下面两行直接计算复数乘法，省去函数调用开销
                    Temp_XX.real = data_of_N_FFT[KB].real * W.real - data_of_N_FFT[KB].imag * W.imag;
                    Temp_XX.imag = W.imag * data_of_N_FFT[KB].real + data_of_N_FFT[KB].imag * W.real;
                    data_of_N_FFT[KB].real = data_of_N_FFT[K].real - Temp_XX.real;
                    data_of_N_FFT[KB].imag = data_of_N_FFT[K].imag - Temp_XX.imag;
                    data_of_N_FFT[K].real = data_of_N_FFT[K].real + Temp_XX.real;
                    data_of_N_FFT[K].imag = data_of_N_FFT[K].imag + Temp_XX.imag;
                }
            }
        }
    }

    public double[] livefft(byte pcm[])
    {

        int len=2048;
        len=len/2;
        int k=0;
        short [] w = new short[len];
        double [] pcmdata = new double[len];
        CREATE_SIN_TABLE();

        for(int i = 0; i < 1024; i++)
        {
            data_of_N_FFT[i] = new complex_of_N_FFT();
        }
        for (int i = 0; i < len; i++)
        {
            w[i] = freadint16_live(pcm[i*2],pcm[i*2+1]);
          // Log.e(TAG, "livefft: "+pcm[i*2]+"===="+pcm[i*2]);
        }
        for (int i = 0; i < N_FFT; i++)//制造输入序列
        {
            data_of_N_FFT[i].real = w[i];
            data_of_N_FFT[i].real /= 32768.0;///=====ba
            data_of_N_FFT[i].imag = 0;
        }
        FFT();
        double add= 0;
        for (int i = 0; i < N_FFT; i++)
        {
            double aa= data_of_N_FFT[i].real * data_of_N_FFT[i].real + data_of_N_FFT[i].imag * data_of_N_FFT[i].imag;
            pcmdata[i] = sqrt(aa);
            add=add+pcmdata[i];
        }
        Log.e(TAG, "livefft:====== "+add/N_FFT);
        Log.e(TAG, "length is:====== "+pcmdata.length);
        return pcmdata;
    }




}

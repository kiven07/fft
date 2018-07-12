using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.ComponentModel;
using System.Drawing;

using System.IO;
using System.Collections;
using UnityEngine;
public class fft : MonoBehaviour
{




    public struct complex_of_N_FFT                         //定义复数结构体
    {
        public double real, imag;
    };
    public char[] head1 = { 'd', 'a', 't', 'a' };
    public const int N_FFT = 1024;                            //傅里叶变换的点数
    public const int M_of_N_FFT = 10;                       //蝶形运算的级数，N = 2^M
    public const int Npart2_of_N_FFT = N_FFT / 2;          //创建正弦函数表时取PI的1/2
    public const int Npart4_of_N_FFT = N_FFT / 4;          //创建正弦函数表时取PI的1/4
    public double[] SIN_TABLE_of_N_FFT = new double[Npart4_of_N_FFT + 1];
    public complex_of_N_FFT[] data_of_N_FFT = new complex_of_N_FFT[N_FFT];


    public int freadint32(BinaryReader br)
    {
        int i;
        byte[] b = new byte[4];
        b = br.ReadBytes(4);
        Debug.Log(b[3]);
        Debug.Log(b[2]);
        Debug.Log(b[1]);
        Debug.Log(b[0]);
        i = b[3];
        i = (i << 8) + b[2];
        i = (i << 8) + b[1];
        i = (i << 8) + b[0];
        return i;
    }
    public Int16 freadint16(BinaryReader br)
    {
        int i;
        byte[] b = new byte[2];
        b = br.ReadBytes(2);
        i = b[1];
        i = (i << 8) + b[0];
        return (Int16)i;
    }

    public bool RunFFT(string FileName, ArrayList wData)
    {
        FileStream fs;
        BinaryReader br;

       
        int i, j;
        int k = 0;
        wData.Clear();
        CREATE_SIN_TABLE();
        fs = new FileStream(FileName, FileMode.Open);
        br = new BinaryReader(fs);
        fs.Seek(38, SeekOrigin.Begin);
        foreach (char c in head1)
        {
            if (c != br.ReadByte())
            {
                Debug.Log("错误");
                fs.Close();
            }
        }
        int len = freadint32(br);
        len = len / 2;
        Int16[] w = new Int16[len];
        for (i = 0; i < len; i++)
        {
            w[i] = freadint16(br);
            if (i == 2001)
            {
                Int16 a2 = w[i];
                Debug.Log(a2);
            }
        }
        fs.Close();
        j = 0;
        while (len > N_FFT)
        {
            for (i = 0; i < N_FFT; i++)//制造输入序列
            {
                data_of_N_FFT[i].real = w[j];
                data_of_N_FFT[i].real /= 32768.0;
                data_of_N_FFT[i].imag = 0;
                j++;
                //printf("%f ", data_of_N_FFT[i].real);
            }
            FFT();
            len -= N_FFT;
          
            
            for (i = 0; i < N_FFT; i++)
            {
                k++;
                wData.Add(Math.Sqrt(data_of_N_FFT[i].real * data_of_N_FFT[i].real + data_of_N_FFT[i].imag * data_of_N_FFT[i].imag));
               

                if (k == 2001)
                {
                    double aa = Math.Sqrt(data_of_N_FFT[i].real * data_of_N_FFT[i].real + data_of_N_FFT[i].imag * data_of_N_FFT[i].imag);

                    double bb = (double)wData[2000];
                    Debug.Log(Math.Sqrt(data_of_N_FFT[i].real * data_of_N_FFT[i].real + data_of_N_FFT[i].imag * data_of_N_FFT[i].imag));
                }
            }
           

        }
        return true;
    }


    public void CREATE_SIN_TABLE()
    {
        int i = 0;
        for (i = 0; i <= Npart4_of_N_FFT; i++)
        {
            SIN_TABLE_of_N_FFT[i] = Math.Sin(Math.PI * i / Npart2_of_N_FFT);//SIN_TABLE[i] = sin(PI2*i/N);
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

    public void ChangeSeat(complex_of_N_FFT[] DataInput)
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

    public void FFT()
    {
        int L = 0, B = 0, J = 0, K = 0;
        int step = 0, KB = 0;
        //ElemType P=0;
        double angle;
        complex_of_N_FFT W, Temp_XX;

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

}


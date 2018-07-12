#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"wav_fft.h"


const int N_FFT = 1024;                            //傅里叶变换的点数
const int M_of_N_FFT = 10;                       //蝶形运算的级数，N = 2^M
const int Npart2_of_N_FFT = 512;          //创建正弦函数表时取PI的1/2
const int Npart4_of_N_FFT =256;          //创建正弦函数表时取PI的1/4
double SIN_TABLE_of_N_FFT[256 + 1];
complex_of_N_FFT  data_of_N_FFT[1024];
void CREATE_SIN_TABLE()
{
	int i = 0;
	for (i = 0; i <= Npart4_of_N_FFT; i++)
	{
		SIN_TABLE_of_N_FFT[i] = sin(PAI * i / Npart2_of_N_FFT);//SIN_TABLE[i] = sin(PI2*i/N);
	}
}
double Sin_find(double x)
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
double Cos_find(double x)
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


void ChangeSeat(complex_of_N_FFT DataInput[])
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

int freadint32(FILE * br)
{
	int i;
	char * buffer;
	size_t result;
	buffer = (char*)malloc(sizeof(char) * 4);
	result = fread(buffer, 1, 4, br);
	int * buffer1;
	buffer1 = (int*)malloc(sizeof(int) * 4);
	for (int i = 0;i < 4;i++)
	{
		buffer1[i] = buffer[i];
		if (buffer1[i] < 0)
		{
			buffer1[i] = 256 + buffer1[i];
		}
	}
	i = buffer1[3];
	i = (i << 8) + buffer1[2];
	i = (i << 8) + buffer1[1];
	i = (i << 8) + buffer1[0];
	free(buffer1);
	free(buffer);
	return i;
}

short freadint16(FILE * br)
{
	int i;
	char * buffer;
	size_t result;
	buffer = (char*)malloc(sizeof(char) * 2);
	result = fread(buffer, 1, 2, br);
	int * buffer1;
	buffer1 = (int*)malloc(sizeof(int) * 2);
	for (int i = 0;i < 2;i++)
	{
		buffer1[i] = buffer[i];
		if (buffer1[i] < 0)
		{
			buffer1[i] = 256 + buffer1[i];
		}
	}
	i = buffer1[1];
	i = (i << 8) + buffer1[0];

	free(buffer);
	free(buffer1);
	return i;
}
void FFT()
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
int RunFFt(char *wavpath) {
	FILE * pFile;
	long lSize;
	int k = 0;
	int j;
	char * buffer;
	size_t result;
	CREATE_SIN_TABLE();
	pFile = fopen(wavpath, "rb");
	//获取文件大小
	fseek(pFile, 0, SEEK_END);
	lSize = ftell(pFile);
	if (pFile == NULL)
	{
		fputs("File error", stderr); exit(1);
	}
	//设置位置
	fseek(pFile, 0, SEEK_SET);
	buffer = (char*)malloc(sizeof(char)*lSize);
	int temp=0;
	int bb=0;
	while(temp==0)
	{
		result = fread(buffer, 1, 4, pFile);
		if(buffer[0]=='d'&&buffer[1]=='a'&&buffer[2]=='t'&&buffer[3]=='a')
		{
			temp=1;
		}
		else
		{
			fseek(pFile, -3, SEEK_CUR);
		}
		bb++;
	}
	printf("add is %d\n",bb);
	int len = freadint32(pFile);  //获取data数据长度
	printf("==len== is %d\n",len);
	len = len / 2;
	short * w;
	w = (short*)malloc(sizeof(short) * len);
	//动态数组接收数据
	for (int i = 0; i < len; i++)
	{

		w[i] = freadint16(pFile);

	}
	fclose(pFile);
	j = 0;
	while (len > N_FFT)
	{
		for (int i = 0; i < N_FFT; i++)//制造输入序列
		{
			data_of_N_FFT[i].real = w[j];
			data_of_N_FFT[i].real /= 32768.0;
			data_of_N_FFT[i].imag = 0;
			j++;
			//printf("%f ", data_of_N_FFT[i].real);
		}
		FFT();
		len -= N_FFT;
		double add=0;
		for (int i = 0; i < N_FFT; i++)
		{
			wavdata[k] = sqrt(data_of_N_FFT[i].real * data_of_N_FFT[i].real + data_of_N_FFT[i].imag * data_of_N_FFT[i].imag);
			add=add+wavdata[k];
			//printf("%f\n",add);
			k++;
		}
		printf("%f\n",add/1024);
	}
	free(w);
	free(buffer);
	return k;
}




#define PAI 3.14159265358979
typedef struct complex_of_N_FFT                         //定义复数结构体
{
	double real, imag;
}complex_of_N_FFT;

double wavdata[50000000];
double a[10000000];
//声明方法
int RunFFt(char *wavpath);

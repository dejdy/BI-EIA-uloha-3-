/*  Short job 2
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <algorithm>
#include <vector>
#include <parallel/algorithm>
#define MAXIMAL_M 20

using namespace std;


void multiply_XY(const unsigned int * array_x,const unsigned int *array_y,const float *array_a,const float *vec_x,float *vec_y,const unsigned int n,const unsigned int nnz)
{
  unsigned int i,x,y;

  for(i=0;i<n;i++) vec_y[i]=0.0f;

  for(i=0;i<nnz;i++)
  {
    x=array_x[i];
    y=array_y[i];
    vec_y[y]+=vec_x[x]*array_a[i];
  }
}


void multiplyT_XY(const unsigned int * array_x, const unsigned int *array_y,const float *array_a,const float *vec_x,float *vec_y,const unsigned int n,const unsigned int nnz)
{
  unsigned int i,x,y;
  for(i=0;i<n;i++) vec_y[i]=0.0;

  for(i=0;i<nnz;i++)
  {
    x=array_x[i];
    y=array_y[i];
    vec_y[x]+=vec_x[y]*array_a[i];
    //printf("MULT %i,%i = %f *%f\n",x,y,array_a[i],vec_x[y]);
  }
}



unsigned int check_sol(float *a,float *b, int n)
{
  int i,j=0,k;
  for(i=0;i<n;i++)
  {
    k=0;
    if (fabs(b[i])>1.0e-3)
    {
      if (fabs(1.0-a[i]/b[i])>1e-3)
        k=1;
    }
    else
    {
      if (fabs( a[i]-b[i])>1e-3)
        k=1;
    }
    if (k)
    {
      if (j<5) printf("Error %i poz %i %g,%g\n",j,i,a[i],b[i]);
      //fflush(stdout);
      j++;
    }

  }
  return j;
}

float num_points(int it,float t)
{
  if (it<=20) return 0.0;
  if ((20<it)&&(it<31)) return ((float)(it-20))/10.0f * 4.0f;
  if ((it==31)) return (11.2f/t) * 9.0f;

  printf("Strange ????");
  return 0.0;
}

void multiply(const unsigned int * cols, const unsigned int *rowPtr,const float *array_a,const float *vec_x,float *vec_y,const unsigned int n)
{
    #pragma omp parallel for
    for (unsigned int i = 0; i < n; i++)
    {
       vec_y[i] = 0.0;
       for (unsigned int k = rowPtr[i]; k < rowPtr[i+1]; k++)
       {
          vec_y[i] += array_a[k]*vec_x[cols[k]];
       }
    }
}

typedef struct matrix
{
    int arr_x;
    int arr_y;
    float val;

} Matrix;



int main(void) {
unsigned long int i,j,k,ok,x,y;
unsigned int n,nnz,rep,kolik,kolik2;
int p1,p2;
double t,start,end,pend,pstart;

float * vec_x;
float * vec_y;
float * vec_y2;
float * vec_t;
float * vec_t2;

  unsigned int * array_x;
  unsigned int * array_y;
  float * array_a;



  n=2000000; //400000;//
//const unsigned int m=2;

    vec_x=new float[n];
    vec_y=new float[n];
    vec_y2=new float[n];
    vec_t=new float[n];
    vec_t2=new float[n];

  // main (parallel) computation
  array_x=new unsigned int[28*n];
  array_y=new unsigned int[28*n];
  array_a=new float[28*n];
  unsigned long int pomo[100];

  nnz=0;
  for(i=0;i<n;i++)
  {
    y=(i*2791+11)%n;
    j=2+std::min(100*y/n,100-100*y/n); // j=3; //
    kolik=0;
    for(k=0;k<j;k++)
    {
      x=(k*102931+y*80147)%n;
      kolik2=1;
      for(ok=0;ok<kolik;ok++)
      {
        if (pomo[ok]==x) { kolik2=0; break;}
      }
      if (kolik2==0) continue;
      pomo[kolik++]=x;
      array_x[nnz]=x;
      array_y[nnz]=y;
      p1=x%29;
      p2=y%31;
      array_a[nnz]=0.009-0.03f*p1+p2*0.01f; //((float)(x%97-49))*0.01+((float)(y%19-9))*0.07;
      //printf("%i (%i,%i) = %f\n",nnz,x,y,array_a[nnz]);
      nnz++;
    }
  }
  printf("n=%i nnz=%i\n",n,nnz);
  fflush(stdout);

  start=omp_get_wtime();

  /*
    you can add some precomputation here
    zde mozno pridat preproccesing
  */

    vector<Matrix> arr(nnz);
    unsigned int * rowPtrCSR = new unsigned int[nnz];
    unsigned int * rowPtrCSC = new unsigned int[nnz];
    float * valsCSR = new float[nnz];
    float * valsCSC = new float[nnz];
    unsigned int * colsCSR = new unsigned int[nnz];
    unsigned int * colsCSC = new unsigned int[nnz];


   for(int i = 0; i < nnz; i++)
   {
        arr[i].arr_x = array_x[i];
        arr[i].arr_y = array_y[i];
        arr[i].val = array_a[i];
   }

    __gnu_parallel::sort(arr.begin(), arr.end(), [](const Matrix & a, const Matrix & b) {return a.arr_y < b.arr_y;});

    rowPtrCSR[0] = arr[0].arr_y;
    int rowPtrSize = 1;
    int temp = arr[0].arr_y;
    for(int i = 0; i < nnz; i++)
    {
        colsCSR[i] = arr[i].arr_x;
        valsCSR[i] = arr[i].val;
        if(arr[i].arr_y != temp)
        {
            rowPtrCSR[rowPtrSize] = i;
            rowPtrSize++;
            temp = arr[i].arr_y;
        }
    }
    rowPtrCSR[rowPtrSize] = nnz;

    __gnu_parallel::sort(arr.begin(), arr.end(), [](const Matrix & a, const Matrix & b) {return a.arr_x < b.arr_x;});

    rowPtrCSC[0] = arr[0].arr_x;
    rowPtrSize = 1;
    temp = arr[0].arr_x;
    for(int i = 0; i < nnz; i++)
    {
        colsCSC[i] = arr[i].arr_y;
        valsCSC[i] = arr[i].val;
        if(arr[i].arr_x != temp)
        {
            rowPtrCSC[rowPtrSize] = i;
            rowPtrSize++;
            temp = arr[i].arr_x;
        }
    }
    rowPtrCSC[rowPtrSize] = nnz;


  rep=0;
  while(rep<30)
  {

    pstart=omp_get_wtime();
    for(i=0;i<n;i++)
    {
      p1=i%23-j%19;
      vec_x[i]=rep*0.002+p1*0.007; //((rep+j)%7-5)*0.1+((rep+i)%17-8)*0.03;
    }

    multiply_XY(array_x,array_y,array_a,vec_x,vec_t,n,nnz);
    multiplyT_XY(array_x,array_y,array_a,vec_x,vec_t2,n,nnz);
    pend=omp_get_wtime();

    start+=pend-pstart;
    // part for modification (beginning)


    multiply(colsCSR, rowPtrCSR, valsCSR, vec_x, vec_y, n);
    multiply(colsCSC, rowPtrCSC, valsCSC, vec_x, vec_y2, n);

    // part for modification (end)
    end=omp_get_wtime();

  i=0;
  i+=check_sol(vec_y,vec_t,n);
  i+=check_sol(vec_y2,vec_t2,n);
  if (i)
  {
  printf("errors =%lu !!!!\n",i);
  fflush(stdout);
  }
  ok=i;
  if ((end-start)>=15.0) ok=1;
  //printf("\n");
  //t=end-start;
  //printf("timeXY=%g\n",t);
    if (ok)  break;
    rep++;
  } // of rep
  t=end-start;
  printf("it=%i      time=%g\n",rep+1,t);
  t=num_points(rep+1,t);
  printf("points=%g\n",t);
    delete[] vec_x;
    delete[] vec_y;
    delete[] vec_y2;
    delete[] vec_t;
    delete[] vec_t2;

  delete[] array_a;
  delete[] array_x;
  delete[] array_y;
  delete[] rowPtrCSC;
  delete[] rowPtrCSR;
  delete[] valsCSC;
  delete[] valsCSR;
  delete[] colsCSR;
  delete[] colsCSC;


  return 0;
}


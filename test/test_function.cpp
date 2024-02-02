#include <R.h>
#include <stdio.h>
#include <cblas.h>
//#include <lapacke.h>


//TEST functions --------------------------------
//tests the performance on an array
//malloc and realloc testings
extern "C" {
  
  /*extern void dgels_(char* trans, int* m, int* n,
                     int* nrhs, double* a, int* lda,
                     double* b, int* ldb,
                     double* work, int* lwork, int* info);
  */
  

  
  
  void double_array(int *y, int *x, int *xl){
    
    //devuelve el puntero a un entero
    int* r = (int *) malloc(15 * sizeof(int *));
    for(int i = 0; i<15; i++){
      r[i] = 1;
      //printf("%d [%d]\n", *(r+i), (int) sizeof(int * ));
    }
    r = (int *) realloc(r, 30 * sizeof(int *));
    char *s = (char *) "max";
    for(int i = 0; i<30; i++){
      printf("%d [%d] [%d] %c\n", *(r+i), (int) sizeof(int * ), (int) strlen(s),
             s[1] );
    }
    
    free(r);
    
    //nótese aquí que el array comienza en posición 0;
    for(int i = 0; i< *xl; i++ )
      y[i] = x[i] * 2;
    
  }
}




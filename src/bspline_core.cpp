#include <R.h>
#include <stdio.h>
#include <stdlib.h> // pulls in declaration of malloc, free, realloc
#include <string.h> // pulls in declaration for strlen.
#include <cblas.h>
#include <math.h>
//#include <stdexcept> //for C++ exception handling


extern "C" {
  
  /*extern void dgtsv_(const long *Np, const long *NRHSp, double *DL,
                     double *D, double *DU, double *B, const long *LDBp,
                     long *INFOp);*/
  
  //LAPACK extern matrix multiplication routine
  extern void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, 
                     double *ALPHA, double * A, int * LDA, double * B,
                     int *LDB, double *BETA, double * C, int *LDC);
  
  // solves minimum square problems by QR decomposition
  extern void dgels_(char* trans, int* m, int* n,
                     int* nrhs, double* a, int* lda,
                     double* b, int* ldb,
                     double* work, int* lwork, int* info);
  
  //it makes the A %*% t(A)
  extern void dsyrk_(char *uplo, char *trans, int *n, int *k, 
                     double *alpha, double *a, int *lda, 
                     double *beta, double *c, int *ldc);
  
  void bs_eval(double *x, int *degree, 
               double *support, int *n_support, 
               double *B, int *nx_B, int **ny_B);
  
  void bspline_eval(double *x, int *degree, 
                    double *support, int *n_support, 
                    double *B, int *nx_B, int *ny_B);
  
  void interpolate_bspline (double *x, int *nx,
                            double *cp, int *ny_cp, int *nx_cp, 
                            int *degree, double *knots, int *n_knots, 
                            int *derivate, double * R, int *is_periodic,
                            int *errcode);
  
  /**
   * prompt.mode == distance for the footpoint computation (implementation)
   */
  void fu_distance( double * prompt_begin, double * prompt_end, int * pr_length, double * cp, 
                    int * cp_x, int * cp_y, int * deg, double * knots, int * knots_n, int * is_periodic, 
                    double * x, int * xcols, int * xrows, double * u, int * un){
    
    double * v = (double *) malloc( *pr_length * *cp_y * sizeof(double *));
    double * param = (double *) malloc( *pr_length * sizeof(double *));
    
    double inc = (*prompt_end - *prompt_begin) / (double)  (*pr_length-1);
    *(param + 0) = * prompt_begin;
    for(int i=1; i < * pr_length; i++) *(param+i) = *(param+i-1) + inc;
    
    int errcode=0, derivate=0;
    interpolate_bspline(param, pr_length,cp, cp_x, cp_y, deg, knots, knots_n, &derivate, v, is_periodic, &errcode);
    
    
    for(int j=0; j< *xrows; j++ ){
      double m_dis = INFINITY;
      int m_ind = -1;  
      for(int i = 0; i< *pr_length; i++){
        double dis = 0;
        
        for(int k=0; k< *xcols; k++){ //*xcols or *cp_y, they must be the same anyway
          double el = v[i + *pr_length * k] - x[j + *xrows * k];
          dis = dis + el * el;
        }
        
        if(dis < m_dis){
          m_dis = dis;
          m_ind = i;
        }
      }
      u[j] = param[m_ind];
    }

    free(param);
    free(v);
    
  }


  /**
   * footpoint implementation for the directional algorithm in the footpoint computation
   */
  void alg_directional(double * cp, 
                    int * cp_x, int * cp_y, int * deg, double * knots, int * knots_n, int * is_periodic, 
                    double * x, int * xcols, int * xrows, double * u, int * un, double * tolerance, double * alpha, 
                    double * vz, int * max_iterations){
    
    double * x0 = (double *) malloc(*cp_y * sizeof(double *));
    double * xd = (double *) malloc(*cp_y * sizeof(double *));
    double * xdd = (double *) malloc(*cp_y * sizeof(double *));
    double * t = (double *) malloc(1 * sizeof(double *));
    int * tn = (int *) malloc(1 * sizeof(int *) );
    double * p = (double *) malloc (*cp_y * sizeof(double *));
    
    *tn = 1;

    double mink=INFINITY, maxk=-INFINITY;
    for(int i = 0; i< *knots_n; i++){
      if(knots[i] < mink) mink = knots[i];
      if(knots[i] > maxk) maxk = knots[i];
    }
    
    for(int i=0; i< *xrows; i++){
      
      *t = u[i];
      
      //RECURRENCE ...  
      //2.-) performs Newton-Gauss minimization on each point
      double tol = INFINITY;
      int count = 0, b = 1;

      //printf("%f %i\n", *t, *tn);
      while(tol > * tolerance){

        int errcode=0, derivate=0;
        interpolate_bspline(t, tn,cp, cp_x, cp_y, deg, knots, knots_n, &derivate, x0, is_periodic, &errcode);
        derivate = 1;
        interpolate_bspline(t, tn,cp, cp_x, cp_y, deg, knots, knots_n, &derivate, xd, is_periodic, &errcode);
        derivate = 2;
        interpolate_bspline(t, tn,cp, cp_x, cp_y, deg, knots, knots_n, &derivate, xdd, is_periodic, &errcode);
        
       
        //printf("%f %f %f\n", x0[0], x0[1], x0[2]);
        //printf("%f %f %f\n", x0[0], x0[2], x0[4]);
        

        double sum_pxd = 0, sum_xdxd = 0, sum_pxdd = 0; //sum_pvz = 0, sum_xdvz = 0, sum_xddvz = 0;
        for(int k = 0; k < *cp_y; k++){ 
          p[k] = x0[k] - x[i + *xrows * k];
          sum_pxd += p[k] * xd[k];
          //sum_pvz += p[k] * vz[k];
          //sum_xdvz += xd[k] * vz[k];
          sum_xdxd += xd[k] * xd[k];
          sum_pxdd += p[k] * xdd[k];
          //sum_xddvz += xdd[k] * vz[k];
        }
        

        double fd = (1 - *alpha) * sum_pxd; //+ *alpha * sum_pvz * sum_xdvz;
        double fdd= (1 - *alpha) * ( sum_xdxd + sum_pxdd ); //+ *alpha * (sum_xdvz * sum_xdvz + sum_pvz * sum_xddvz);
        
        double inc_t = -fd /fdd;
        count ++;
        tol = sum_pxd * sum_pxd;
        *t +=  inc_t;

        //printf("%f %f %f\n", xd[0], xd[1], xd[2]);
        
        if(count > *max_iterations || *t < mink || *t > maxk){
          b = 0;
          break;
        }
        
        //break;
      }
      
      if(b == 1) { u[i] = *t; }
      
      
    }
    
    free(x0);
    free(xd);
    free(xdd);
    free(t);
    free(tn);
    free(p);
    
  }
  
  /**
   * footpoint implementation for the rogers algorithm in the footpoint computation
   */
  void alg_rogers(double * cp, 
                       int * cp_x, int * cp_y, int * deg, double * knots, int * knots_n, int * is_periodic, 
                       double * x, int * xcols, int * xrows, double * u, int * un, double * tolerance, double * alpha, 
                       double * vz, int * max_iterations, double * step){
    
    double * x0 = (double *) malloc(*cp_y * sizeof(double *));
    double * xd = (double *) malloc(*cp_y * sizeof(double *));
    
    double * t = (double *) malloc(1 * sizeof(double *));
    int * tn = (int *) malloc(1 * sizeof(int *) );
    //double * p = (double *) malloc (*cp_y * sizeof(double *));
    double * qs = (double *) malloc ( *cp_y * sizeof(double *) );
    
    *tn = 1;
    
    double mink=INFINITY, maxk=-INFINITY;
    for(int i = 0; i< *knots_n; i++){
      if(knots[i] < mink) mink = knots[i];
      if(knots[i] > maxk) maxk = knots[i];
    }
    
    for(int i=0; i< *xrows; i++){
      
      *t = u[i];
      
      //RECURRENCE ...  
      //2.-) performs Newton-Gauss minimization on each point
      double tol = INFINITY;
      int count = 0;
      
      while(tol > * tolerance){
        
        int errcode=0, derivate=0;
        interpolate_bspline(t, tn,cp, cp_x, cp_y, deg, knots, knots_n, &derivate, x0, is_periodic, &errcode);
        derivate = 1;
        interpolate_bspline(t, tn,cp, cp_x, cp_y, deg, knots, knots_n, &derivate, xd, is_periodic, &errcode);

        double sum_xdxc = 0, sum_xdxd = 0;
        for(int k = 0; k < *cp_y; k++){
          sum_xdxc += *(xd + k) * ( *(x + i + *xrows * k) - *(x0 + k) );
          sum_xdxd += *(xd + k) * *(xd + k);
        }

        double dt=0, tolarg=0;
        double f = sum_xdxc / sum_xdxd;
        for(int k = 0; k< *cp_y; k++){
          qs[k] = x0[k] + f * xd[k];
          tolarg += (qs[k] - x0[k]) * (qs[k] - x0[k]);
        }
        dt = *step * f;
        
        //printf("%f %f %f %f %f ***\n",dt, f, tolarg, sum_xdxc, sum_xdxd);
        
        if(*t + dt < mink || *t + dt > maxk) // if the point is projected to the tip
          break;
        
        
        *t += dt;
        tol = sqrt(tolarg); 
        count ++;
        
        if(count > *max_iterations)
          break;

      }
      
      u[i] = *t;
      
      
    }
    
    free(x0);
    free(xd);
    free(t);
    free(tn);
    free(qs);
    
  }
  
  
  /**
   * 
   */
  void bspline_footpoint(double * u, int * un,double * x, int * xcols, int *xrows, double * cp, int * cp_x, int * cp_y,
                         int *deg, double * knots, int * knots_n, int * pm, 
                         double * prompt_begin, double * prompt_end,
                         int * init_scan_points, int * algorithm, int * is_periodic, double * tolerance, double * alpha, 
                         double * vz, int * max_iterations, double * step){
    
    //selectes the prompt.mode first
    //if(*is_u_null == 1){
    switch(*pm){
      case 1 :
      case 2 :
      default:
        //fprintf(stderr,"prompt.mode orthogonal or directional not allowed in native mode... changing to prompt mode distance\n");
        fu_distance(prompt_begin, prompt_end, init_scan_points, cp, cp_x, cp_y, deg, knots, knots_n, is_periodic,
                    x, xcols, xrows, u, un);
        
    }
    //}
    
    
    // THE parameter correction algorithm is here ...
    // algorithm none = 0, Rogers = 1, orthogonal = 2, directional = 3.
    switch(*algorithm){
      case 0:
        //fprintf(stderr,"no algorithm selected. Prompted values returned instead ...\n");
        break;
      case 1:
        //fprintf(stderr,"Rogers algorithm selected.\n");
        alg_rogers(cp, cp_x, cp_y, deg, knots, knots_n, is_periodic, x, xcols, xrows, 
                   u, un, tolerance, alpha, vz, max_iterations, step);
        break;
      case 2:
        //fprintf(stderr,"Algorithm orthogonal not implemented in native mode. Changing to directional algorithm.\n");
      case 3:
      default : 
        //fprintf(stderr,"Directional algorithm selected.\n");
        alg_directional(cp, cp_x, cp_y, deg, knots, knots_n, is_periodic,
                     x, xcols, xrows, u, un, tolerance, alpha, vz, max_iterations);
    }
  
  
  }
  
  
  
  /**
   * It finds the optimal solution of ||B - A*X||. This is an internal function.
   * The arguments have the same meaning as in the lapack entry for dgels database.
   * 
   * Use this function like this:     bspline_dgels(&m, &n, &nrhs, A, b);
   * assuming of course that all the arguments come from the R environment an here are
   * handled by pointers.
   */
  double * dgels(int m, int n, int nrhs, double *A, double *b){
    
    char transa ='N';
    int info, lwork=-1, lda=m, ldb=m;
    double* work;
    
    /* this block makes a query in order to find the working space for the problem */
    double query;
    dgels_(&transa, &m, &n, &nrhs, A, &lda, b, &ldb, &query, &lwork, &info);
    lwork=(int)query;
    
    /* the working space for the problem is reserved and the solution to the problem is found */
    work= (double *) malloc(lwork*sizeof(double));
    if(work==NULL){fprintf(stderr,"memory allocation failed\n");exit(1);}
    
    dgels_(&transa, &m, &n, &nrhs, A, &lda, b, &ldb, work, &lwork, &info);
    
    //copies the right part of the result in other allocation
    double * result = (double *) malloc(n * nrhs * sizeof(double));
    
    
    if(n > m){fprintf(stderr,"number of columns in dgels is lower than number of rows"); exit(1);}
    
    int i,j;
    for(j=0; j<nrhs; j++ )
      for(i=0; i<n; i++)
        result[i+j * n] = b[i + j * m];

    free(work);
    
    return(result);
  }

  
  /**
   * u, nu - the sequence of parameter values corresponding to the best footpoints guesses (and its length)
   * x, xn - the data points in the form of a matrix of row-vectors (the linearization in R is made by as.vector(x) or as.numeric(x) )
   * degree - the degree of the polynomial spline
   * support, n_support - the current support of nodes, and its length
   * R - the result matrix in the dgels minimization problem. 
   * regularization_term - the regularization term as displayed in the R function
   */
  void bspline_solve(double *u, int *nu, int *degree, 
                     double *support, int *n_support, double *x,int *xn, double * R, double * regularization_term){

    // This block obtains the spline matrix evaluation for the current problem
     int ny_B = *n_support - *degree - 1;
     int nx_B = *  nu;
     double * B = (double *) malloc( nx_B * ny_B * sizeof(double *));
     bspline_eval(u, degree, support, n_support, B, &ny_B, &nx_B);
     
     
     // This block generates the Nprod in the R implementation, N %*% t(N) only using lapack functions
     double * BB = (double *) malloc (ny_B * ny_B * sizeof(double *));
     //dsyrk(nx_B, ny_B, B, BB);
     char uplo = 'L', trans = 'T';
     double alpha = 1, beta = 0;
     int lda = nx_B, ldc = ny_B;
     dsyrk_(&uplo, &trans, &ny_B, &nx_B, &alpha, B, &lda, &beta, BB, &ldc);
     uplo = 'U';
     dsyrk_(&uplo, &trans, &ny_B, &nx_B, &alpha, B, &lda, &beta, BB, &ldc);
     for(int i = 0; i<ny_B;i++) //temporal regularization term
       BB[i+i*ny_B] = BB[i+i*ny_B] * *regularization_term * *regularization_term;
     
     
     // This block generates the product N %*% x in the R implementation, using lapack functions
     char transa = 'T', transb = 'N';
     int m = ny_B, n = *xn, k= nx_B;
     lda = k;
     int ldb = k;
     ldc = m;
     double * Bx = (double *) malloc( m * n *  sizeof(double *) );
     dgemm_(&transa, &transb, &m,&n,&k, &alpha, B, &lda, x, &ldb, &beta, Bx, &ldc);

     
     // This block solves the minimum square problem
     int m2 = ny_B, n2 = ny_B, k2 = *xn;
     double * result = dgels(m2,n2,k2,BB,Bx);
     
     // insterts the data into the result matrix
     for(int i = 0; i< n2 * k2; i++) R[i] = result[i];

     free(B);
     free(BB);
     free(Bx);
     free(result);
  }
  
 /*
  * this function represents the native part for the function
  * bspline_support
  * 
  * *degree is not required !!! it is only used for the multiplicites, which
  * have been alreade precalculated. WE HAVE TO ELIMINATE *degree param
  */
 void bspline_support(int *degree, double *knots, int *n_knots, 
                      int *multiplicities, double *support, int *errcode) {
   
   
   //tests if the knots are unique and increasing
   for(int i=0; i < *n_knots - 1; i++ )
     if(knots[i] >= knots[i+1]){
       *errcode = 1;
       return;
       //alternative take:
       //throw handling does not work apparently; here's the line
       //throw std::invalid_argument( "knots must be unique"); 
     }
    
   
   int c = 0;
   for(int i = 0; i < *n_knots; i++)
     for(int j = 0; j < multiplicities[i]; j++){
       support[c] = knots[i];
       c++;
     }
   
 }
  
/*
 * this function represents the native part for the function
 * bspline_find_interval_index.
 */
 void bspline_find_interval_index(double* x, double * support, 
                                  int *n_support, int *r, int *errcode){
     for(int i = 0; i < *n_support -1; i++)
       if(*x >= support[i] && *x < support[i+1]){ 
         *r=i+1;
         break;
       }
       
     int i=0;
     if(*r == 0 && *x >= support[*n_support-1] )
        while(support[*n_support-i-1] == support[*n_support-i-2]){
           i++;
           *r = *r - 1;
         }
 }
 /*
  * bspline_find_interval_index implementation for inner purposes.
  */
 void bs_fii(double **x, double *support, 
                  int **n_support, int **k){
  
    for(int i = 0; i < **n_support -1; i++)
      if(**x >= support[i] && **x < support[i+1]){ 
        **k=i+1;
        break;
      }
      
    int i=0;
    if(*k == 0 && **x >= support[**n_support-1] )
      while(support[**n_support-i-1] == support[**n_support-i-2]){
        i++;
        **k = **k - 1;
      }
     
 }
  
/*
 * this function represents the native part for the function
 * bspline.
 */
 void bspline(double *x, int *k, int *degree, 
              double *support, int *n_support, 
              double *B, int * n_B, int *errcode){
   
   //ALL CODE HERE
   //makes a default spline if the last point has been reached
   //i.e. the user tries to interpolate the last point
   if(*x == support[*n_support-1]){
     *errcode=2;
     B[*n_support-*degree-2]=1;
     return;
   }
   
   //finds the interval index in case it is not specified
   if(*k==-1)
     bs_fii(&x, support, &n_support, &k);
   
   //1.) builds the 0-degree B-spline
   B[*k-1] = 1;
   
   //2.) evaluates using the recursion theorem for B-splines (de Boor)
   int p = 0;   //the current degree in the routine
   double t1,t2;
   while(p < *degree ){
       p = p+1;
       
       // it evaluates the middle b-splines. Both terms may be not-null
       for(int i=(*k-p-1); i < *k; i++) {
         t1 = 0; t2 = 0;
         
         if( support[i+p] > support[i] )
           t1 = B[i] * (*x - support[i]) / ( support[i+p] - support[i] ) ;

         if( support[i+p+1] > support[i+1] )
           t2 = B[ (i % *n_B)+1 ] * ( support[i+p+1] - *x ) / ( support[i+p+1] - support[i+1] ) ;
         
         B[i] = t1 + t2;
       }
    }
   
     
  }  
  
 /*
  * bspline function implementation for inner purposes
  */
 void bs(double *x, int *k, int **degree, 
         double *support, int **n_support, 
         double *B, int **n_B){
   //ALL CODE HERE
   //makes a default spline if the last point has been reached
   //i.e. the user tries to interpolate the last point
   if(*x == support[**n_support-1]){
     B[**n_support-**degree-2]=1;
     return;
   }
   
   //finds the interval index in case it is not specified
   if(*k==-1)
     bs_fii(&x, support, &*n_support, &k);
   
   //1.) builds the 0-degree B-spline
   B[*k-1] = 1;
   
   //2.) evaluates using the recursion theorem for B-splines (de Boor)
   int p = 0;   //the current degree in the routine
   double t1,t2;
   while(p < **degree ){
     p = p+1;
     
     // it evaluates the middle b-splines. Both terms may be not-null
     for(int i=(*k-p-1); i < *k; i++) {
       t1 = 0; t2 = 0;
       
       if( support[i+p] > support[i] )
         t1 = B[i] * (*x - support[i]) / ( support[i+p] - support[i] ) ;
       
       if( support[i+p+1] > support[i+1] )
         t2 = B[ (i % **n_B)+1 ] * ( support[i+p+1] - *x ) / ( support[i+p+1] - support[i+1] ) ;
       
       B[i] = t1 + t2;
     }
   }
 }
     /*
  * this function represents the native part for the function
  * bspline_eval.
  */  
 void bspline_eval(double *x, int *degree, 
                   double *support, int *n_support, 
                   double *B, int *nx_B, int *ny_B){
   
   double xi;
   int k,i,j;
   
   double * Brow = (double *) malloc(*nx_B * sizeof(double));
   //sets to zero
   for(i = 0; i<*nx_B; i++)
     Brow[i] = 0;
   
   for(j = 0; j<*ny_B; j++){
    
    xi = x[j];
    k = -1;
    bs(&xi, &k, &degree, support, &n_support, Brow, &nx_B);
    
    //assigns the calculated values
    for(i = 0; i<*nx_B; i++){
        B[j + *ny_B*i] = Brow[i];
        Brow[i] = 0;
    }
   
   }
   free(Brow); 
   
 }
  
 /*
  * this function represents the native part for the function
  * recompute_control_points
  * v : array of control points
  * ny_v : number of rows (= number of control points)
  * nx_v : number of columns (= number of coordinates of the control points)
  * U : the knots support
  * degree : the spline degree
  */
 void recompute_control_points(double *v,
                               int *ny_v,
                               int *nx_v,
                               double *U,
                               int *degree){
   
   int i,j, nyn_v, co;
   double degd = (double) *degree;
   //nyn_v = *ny_v - 1; 
   double * new_cpoints = (double *) malloc(*ny_v * *nx_v * sizeof(double));
   
   //makes the proper calculation
   for(j = 0; j <  *ny_v; j++)    //iterates over the rows (coordinates)
   for(i = 0; i <  *nx_v; i++)    //iterates over the cols 
       new_cpoints[j* *nx_v + i] = ( degd/( U[j + *degree + 1] - U[j + 1] ) ) * (v[j + 1 + i * *ny_v] - v[j + i * *ny_v]);
   
   
   //reasigns to old variable
   for(j = 0; j <  *ny_v; j++)    
   for(i = 0; i <  *nx_v; i++) {
       co = j * *nx_v + i;
       v[co] = new_cpoints[co];
   }

   free(new_cpoints);
 }
     
     
  
  /*
   * this function represents the native part of the function
   * transform_support_periodic
   * degree : an integer expressing the degree of the spline
   * support : a double array, expressing a clamped support
   * n_support: the length of the support (number of knots) 
   */
  void transform_support_periodic(int *degree, 
                                  double *support,
                                  int *n_support){
    
    const int a  = *n_support-*degree-1, b=*n_support-(2 * *degree), aminusb=a-b+1;
    double* lint = (double *) malloc( aminusb * sizeof(double)); 
    double sum = 0;
    
    for(int i=a; i>=b; i--) lint[aminusb-i+b-1]=support[i] - support[i-1];  
    
    for(int i=0; i<aminusb; i++) {
      sum = 0;
      for(int j=i; j<aminusb; j++) sum += lint[j];
      lint[i] = sum;
    }
    
    double* sup_per = (double *) malloc( *n_support * sizeof(double)); //allocation for the periodic support
    for(int i = 0; i<*n_support; i++) sup_per[i] = support[i];         //makes a copy of the support
    
    for(int i=0; i< *degree;i++) sup_per[i] = support[*degree - 1] - lint[i]; //makes the lower bounding of the periodic spline

    const int a2 = 1 + *degree,b2 = 2* *degree;
    for(int i = a2; i<=b2; i++) lint[i-a2] = support[i] - support[i-1];
    double * lint2 = (double *) malloc ( sizeof(double) * aminusb);
    for(int i=0; i<aminusb; i++) {
      sum = 0;
      for(int j=0; j<=i; j++) sum += lint[j];
      lint2[i] = sum;
    }
    
    for(int i = a+1; i<= *n_support-1; i++ ) sup_per[i] = support[*n_support - *degree-1] + lint2[i-a-1];
    for(int i = 0; i < *n_support; i++) support[i] = sup_per[i];
        
    free(sup_per);
    free(lint);
    free(lint2);
  }
  
  
  // INTERPOLATE_BSPLINE implementation here ----------------------------------------------
  /*
   * transform_support_periodic for inner purposes
   */
  void trans_sup_per(int *degree, 
                                  double *support,
                                  int *n_support){
    
    const int a  = *n_support-*degree-1, b=*n_support-(2 * *degree), aminusb=a-b+1;
    double* lint = (double *) malloc( aminusb * sizeof(double)); 
    double sum = 0;
    
    for(int i=a; i>=b; i--) lint[aminusb-i+b-1]=support[i] - support[i-1];  
    
    for(int i=0; i<aminusb; i++) {
      sum = 0;
      for(int j=i; j<aminusb; j++) sum += lint[j];
      lint[i] = sum;
    }
    
    double* sup_per = (double *) malloc( *n_support * sizeof(double)); //allocation for the periodic support
    for(int i = 0; i<*n_support; i++) sup_per[i] = support[i];         //makes a copy of the support
    
    for(int i=0; i< *degree;i++) sup_per[i] = support[*degree - 1] - lint[i]; //makes the lower bounding of the periodic spline
    
    const int a2 = 1 + *degree,b2 = 2* *degree;
    for(int i = a2; i<=b2; i++) lint[i-a2] = support[i] - support[i-1];
    double * lint2 = (double *) malloc ( sizeof(double) * aminusb);
    for(int i=0; i<aminusb; i++) {
      sum = 0;
      for(int j=0; j<=i; j++) sum += lint[j];
      lint2[i] = sum;
    }
    
    for(int i = a+1; i<= *n_support-1; i++ ) sup_per[i] = support[*n_support - *degree-1] + lint2[i-a-1];
    for(int i = 0; i < *n_support; i++) support[i] = sup_per[i];
    
    free(sup_per);
    free(lint);
    free(lint2);
  }
 
  
  /*
   * bspline_eval function for inner purposes 
   */
  void bs_eval(double *x, int *degree, 
               double *support, int *n_support, 
               double *B, int *nx_B, int **ny_B){
    double xi;
    int k,i,j;
    
    double * Brow = (double *) malloc(*nx_B * sizeof(double));
    //sets to zero
    for(i = 0; i<*nx_B; i++)
      Brow[i] = 0;
    
    for(j = 0; j<**ny_B; j++){
      
      xi = x[j];
      k = -1;
      bs(&xi, &k, &degree, support, &n_support, Brow, &nx_B);
      
      //assigns the calculated values
      for(i = 0; i<*nx_B; i++){
        B[j + **ny_B*i] = Brow[i];
        Brow[i] = 0;
      }
      
    }
    free(Brow); 
  }
     
  /*
   * recompute_control_points implementation for inner purposes 
   * v : array of control points
   * ny_v : number of rows (= number of control points)
   * nx_v : number of columns (= number of coordinates of the control points)
   * U : the knots support
   * degree : the spline degree
   */
  void recomp_cp           (double *v,
                                int *ny_v,
                                int *nx_v,
                                double *U,
                                int *degree){
    
    int i,j, n_cp = *ny_v * *nx_v;
    double degd = (double) *degree;
    double * new_cpoints = (double *) malloc(n_cp * sizeof(double));

    //printf("%d : %d : %d\n", *ny_v, **nx_v, *degree);
      
    //makes the proper calculation
    for(i = 0; i <  *ny_v; i++)    //iterates over the rows (coordinates)
      for(j = 0; j <  *nx_v; j++)    //iterates over the cols 
        new_cpoints[j* *ny_v + i] = ( degd/( U[i + *degree + 1] - U[i + 1] ) ) * 
                                    (v[i + 1 + j * (*ny_v + 1)] - v[i + j * (*ny_v + 1)]);
      
    
    v = (double *) realloc(v, sizeof(double) * n_cp);
    for(i = 0; i <  *ny_v; i++)    //iterates over the rows (coordinates)
      for(j = 0; j <  *nx_v; j++)    //iterates over the cols 
        v[j * *ny_v + i] = new_cpoints[j * *ny_v + i];
    
    free(new_cpoints);
    
  }
  
  /*
   * bspline_support implementation for inner purposes
   */
  void bs_sup(double * support, double *knots, int **n_knots, 
              int *multiplicities, int *sum_multiplicities, int **errcode) {
    
    support = (double *) realloc(support,*sum_multiplicities * sizeof(double));
    
    //tests if the knots are unique and increasing
    for(int i=0; i < **n_knots - 1; i++ )
      if(knots[i] >= knots[i+1]){
        **errcode = 1;
        return;        
      }
      
    int c = 0;
    for(int i = 0; i < **n_knots; i++)
      for(int j = 0; j < multiplicities[i]; j++){
        support[c] = knots[i];
        c++;
      }
    
  }
  
  /*
   * compute multiplicities for inner purposes only
   */
  int * compute_multiplicites( int *degree,
                              int **n_knots,
                              int *sum_multiplicities
                              ){
    int * multiplicities = (int *) malloc(**n_knots * sizeof(int));
    for (int i = 0; i < **n_knots; i++){
      if(i == 0 || i == **n_knots - 1)
        *(multiplicities + i) = *degree + 1;
      else
        *(multiplicities + i) = 1;
      *sum_multiplicities += multiplicities[i];
    }
    
    return(multiplicities) ;
  }
  
  
  
  /*
   * the first thing we need to ensure is optimal matrix multiplication, using the known examples
   * 
   * NOTE: we need to include the header in the R function that prepares for the periodic splines
   */
  void interpolate_bspline (double *x, int *nx,
                            double *cp, int *ny_cp, int *nx_cp, 
                            int *degree, double *knots, int *n_knots, 
                            int *derivate, double * R, int *is_periodic,
                            int *errcode){
    
    if(*derivate < 0 || *derivate > *degree){
      fprintf(stderr,"[Max] the derivate cannot be negative or higher than the degree.\n");
      exit(1);
    }
    if(*n_knots <= 1){
      fprintf(stderr,"[Max] the minimum length for the knot sequence is 2. Either increase the number of control points or decrease the current degree.\n");
      exit(1);
    }
    if(*ny_cp <= 1){
      fprintf(stderr,"[Max] the minimum number of control points is 2, don't be silly!\n");
      exit(1);
    }
    
    //orders the knots
    int flag;
    do{
      for(int i=0; i< *n_knots -1; i++){
        int flag = 0;
        if(knots[i] > knots[i+1]){
          double tr = knots[i+1];
          knots[i+1] = knots[i];
          knots[i] = tr;
          flag = 1;
        }
      }
    }while(flag == 1);
    
    // ensures the equivalence bewteen knots limits and the control points
    double mink=INFINITY, maxk=-INFINITY;
    for(int i = 0; i< *n_knots; i++){
      if(knots[i] < mink) mink = knots[i];
      if(knots[i] > maxk) maxk = knots[i];
    }
      
    if(*is_periodic == 1)
      for(int i = 0; i< *nx; i++) if(x[i] == maxk) x[i] = mink;
      
      
    
    

    //assign the multiplicities values
    int sum_multiplicities = 0, cur_deg = *degree, nny_cp = *ny_cp, nnx_cp = *nx_cp, 
        *multiplicities;
    double * support, * ccp, * B, alpha = 1, beta=0;
    double * kknots;
    char transa = 'n', transb = 'n';
    int n_cp;
    int * nn_knots = (int *) malloc(sizeof(int*));
    //makes a copy of the control points
    if(*is_periodic == 1){
      nny_cp = *ny_cp + *degree;
      n_cp = nny_cp * nnx_cp;
      
      ccp = (double *) malloc(n_cp * sizeof(double));
      for(int i=0; i<*ny_cp; i++) for(int j=0; j<*nx_cp; j++) 
      *(ccp + i + nny_cp*j) = *(cp + i + *ny_cp*j );
      
      for(int i=*ny_cp; i<nny_cp; i++) for(int j=0; j<*nx_cp; j++) {
        *(ccp + i + nny_cp*j) = *(cp + (i-*ny_cp) + *ny_cp*j);
      }
      
      // changes the knots sequence too
      // knots <- seq(from=knots[1], to=knots[length(knots)], length.out=nrow(vc)-degree+1)
      *nn_knots = nny_cp - *degree+1;
      kknots = (double *) malloc (*nn_knots * sizeof(double *));
      double inc = (maxk - mink) / (double)  (*nn_knots-1);
      *(kknots + 0) = mink;
      for(int i=1; i < *nn_knots; i++) *(kknots+i) = *(kknots+i-1) + inc;
      //fprintf(stderr,"[Max] knots have been evenly distributed over the interval to fit the periodic spline operation. The original knot sequence is lost.\n");
    
      
    }else{
      *nn_knots = *n_knots;
      n_cp = nny_cp * nnx_cp;
      ccp = (double *) malloc(n_cp * sizeof(double));
      kknots = (double *) malloc (*nn_knots * sizeof(double *));
      for(int i=0; i<n_cp; i++) *(ccp + i) = *(cp + i);
      for(int i=0; i < *nn_knots; i++) *(kknots+i) = *(knots+i);
    }
    
    //generates the support before computation
    cur_deg = *degree;
    sum_multiplicities = 0;
    multiplicities = compute_multiplicites(&cur_deg, &nn_knots, &sum_multiplicities);
    support = (double *) malloc(sum_multiplicities * sizeof(double));

    if(*derivate > 0)
       for(int i = 1; i <= *derivate; i++){
         
         cur_deg = *degree - i + 1;
         sum_multiplicities = 0;
         multiplicities = compute_multiplicites(&cur_deg, &nn_knots, &sum_multiplicities);
         bs_sup(support, kknots, &nn_knots, multiplicities, &sum_multiplicities, &errcode );
         if(*is_periodic == 1)
           trans_sup_per(&cur_deg, support, &sum_multiplicities);
         
         nny_cp -= 1;
         //double * test;
         recomp_cp(ccp,&nny_cp,&nnx_cp , support ,&cur_deg); 
         
    }
    
    cur_deg = *degree - *derivate;
    sum_multiplicities = 0;
    multiplicities = compute_multiplicites(&cur_deg, &nn_knots, &sum_multiplicities);
    bs_sup(support, kknots, &nn_knots, multiplicities, &sum_multiplicities, &errcode );
    if(*is_periodic == 1)
      trans_sup_per(&cur_deg, support, &sum_multiplicities);
    
    //evaluates the spline on the parameters values (x vector)
    int Bnx = sum_multiplicities - cur_deg - 1;
    B = (double *) malloc(*nx * Bnx * sizeof(double));
    bs_eval(x, &cur_deg, support, &sum_multiplicities, B, &Bnx, &nx);

    //applies products of two matrices at last
    
    dgemm_(&transa, &transb, nx, nx_cp, &Bnx, &alpha, B, 
           nx, ccp, &Bnx, &beta, R, nx);
    
    // this block corrects the last evaluation. The result at the right end of the parameter space (knots) 
    // must match the position of the last 
    // control point
    for(int i =0; i< *nx; i++){
        if( x[i] >=  maxk)
        for(int j = 0; j< *nx_cp; j++) R[i + j* *nx] = ccp[(j+1) *nny_cp-1]; //
    }
    
    
    //free memory space 
    free(multiplicities);
    free(support);
    free(ccp);
    free(B);
    free(nn_knots);
    free(kknots);
    
  }
  
  
  

  
}




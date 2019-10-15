/*
 * COMP 633:  sequential implementation for n bodies (pa1a)
 *            V3:  half-pair interactions
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define  NMAX    50000


/*
 * attributes of a body in 2D
 */
struct Body {
  double x, y;
  double m;
  double vx,vy;  // velocity
  double fx,fy;  // total force (instead of acceleration)
};

/*
 * n-body system
 */
static struct Body     P[NMAX];

double G      = 6.673E-11;
double DeltaT = 0.001;
int    nts    = 4; 
int    n;


int main(int argc, char * argv[])
{

  int h, i, j, k;
  double t1,t2;

  // how many bodies
  if (argc > 1) {
    n = atoi(argv[1]);
    n = (n > 0 && n <= NMAX) ? n : NMAX;
  }

  /* 
   * initial state: place n bodies along diagonal of a unit square, 
   * each with unit mass and zero initial velocity and acceleration
   * and total force
   */
  for (i = 0; i < n; i++) {
    P[i].x = P[i].y = (double) i / (double) n;
    P[i].m = 1.0;
    P[i].vx = P[i].vy = 0.0;
    P[i].fx = P[i].fy = 0.0;
  }

  /*
   *  n-body simulation for nts steps
   */
  t1 = omp_get_wtime();

  for (k = 1; k <= nts; k++) {
 #pragma omp parallel shared(P) private(i, j)
   {
   #pragma omp for schedule(guided)
    for (i = 0; i < n; i++) {

      double Fx = 0.0, Fy = 0.0;

      /* 
       * half-pairs:  each pairwise interaction is evaluated once 
       */
      for (j = i + 1; j < n;  j++) {

        double rx, ry, d2, d, c;
        
        /* interaction caluclation
        */
        rx = P[j].x - P[i].x;   
        ry = P[j].y - P[i].y;   
        d2 = (rx * rx) + (ry * ry);
        c  = P[j].m * P[i].m / (d2 * sqrt(d2));           
        Fx += c * rx;              
        Fy += c * ry;  
<<<<<<< HEAD
        #pragma omp atomic
        P[j].fx -= c * rx;  // reciprocal interaction
        #pragma omp atomic
=======
        #pragma omp atomic 
        P[j].fx -= c * rx;  // reciprocal interaction
        #pragma omp atomic 
>>>>>>> 78740f9a34447983ecc4e6f9fe937f7ffe251f06
        P[j].fy -= c * ry;             

        } /* j */
	//Lock P[i]

      #pragma omp atomic 
      P[i].fx +=  Fx;
      
      #pragma omp atomic 
      P[i].fy +=  Fy;
	//Unlock P[i]
    } /* i */
   }/* End Parallel */
    /* advance velocities and positions 
     */
 #pragma omp parallel shared(P) private(i, j)
   {
   #pragma omp for schedule(guided)
    for (i = 0; i < n; i++) {
      P[i].x  += P[i].vx * DeltaT;
      P[i].y  += P[i].vy * DeltaT;
      P[i].vx += ((G * P[i].fx) / P[i].m) * DeltaT; 
      P[i].vy += ((G * P[i].fy) / P[i].m) * DeltaT; 
      P[i].fx = P[i].fy = 0.0;  // reset forces for next time step
    }
   }/* End Parallel */
  } /* k */

  t2 = omp_get_wtime();
  
  /* report result in units of millions of interactions per second 
   */
  {
    double interactions = nts * (double) n * (double) n;
    double dt = t2 - t1;
    printf("%.5g \n",1E-6 * interactions / dt); 
  }

} 




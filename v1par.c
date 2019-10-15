/*
 * COMP 633:  sequential implementation for n bodies (pa1a)
 *            V1:  all-pair interactions
 *                 direct transcription of spec
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define  NMAX    10000

/*
 * attributes of a body in 2D
 */
struct Body {
  double x,y;    // position
  double m;      // mass
  double vx,vy;  // velocity
  double ax,ay;  // acceleration
};

/*
 * n-body system
 */
static struct Body P[NMAX];

double G      = 6.673E-11;
double DeltaT = 0.001;
int    nts    = 6;     // number of time steps
int    n;              // number of bodies



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
   * initial state: place n bodies along diagonal of unit square, 
   * each with unit mass and zero initial velocity and acceleration
   */
 #pragma omp parallel shared(P) private(i)
	{
	   #pragma omp for schedule(guided)
	  for (i = 0; i < n; i++) {
	    P[i].x = P[i].y = (double) i / (double) n;
	    P[i].m = 1.0;
	    P[i].vx = P[i].vy = 0.0;
	    P[i].ax = P[i].ay = 0.0;
	  }
	}


  /*
   *  n-body simulation for nts steps
   */
t1 = omp_get_wtime();
for (k = 1; k <= nts; k++) {
 #pragma omp parallel shared(P) private(k,i,j)
	{

	   #pragma omp for schedule(guided)
	    for (i = 0; i < n; i++) {

	      /*  force accumulator for body i
	       */
	      double fx = 0.0;
	      double fy = 0.0;

	      /*  interact with all bodies (except self)
	       */
	      for (j = 0 ; j < n;  j++) {

		if (j != i) { 
		  double rx, ry, d, d3, c;

		  /* direct transcription of f_ij 
		   */
		  rx = P[j].x - P[i].x;   
		  ry = P[j].y - P[i].y;
		  d  = sqrt((rx * rx) + (ry * ry));
		  d3 = pow(d,3.0);
		  fx += G * P[i].m * P[j].m * rx / d3;
		  fy += G * P[i].m * P[j].m * ry / d3;
		} /* if */

	      }/* for j */

	      P[i].ax = fx / P[i].m;
	      P[i].ay = fy / P[i].m;

	    }/* for i */
	 }/*End parallel*/
    /* update velocities and positions 
     */
 #pragma omp parallel shared(P) private(k,i,j)
	{
	   #pragma omp for schedule(guided)
    for (i = 0; i < n; i++) {
      P[i].x  += P[i].vx * DeltaT;
      P[i].y  += P[i].vy * DeltaT;
      P[i].vx += P[i].ax * DeltaT; 
      P[i].vy += P[i].ay * DeltaT; 
    }
   }/*End parallel*/
 } /* for k */
  t2 = omp_get_wtime();

  /* report result in units of millions of interactions per second 
   */
  {
    double interactions = nts*n*n;
    double dt = t2 - t1;
    printf("%.5g \n", 1E-6 * interactions / dt); 
  }
}

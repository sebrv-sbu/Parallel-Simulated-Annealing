#ifdef ICC
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/***********************************************************************
 * fit_log.c created 6-02 Lorraine Greenwald
 * USE LOG BASE 10 rather than natural log  
 * need to measure the vertical drop (y-intercept) between the curves 
 * to show improvement in the algorithm.  
 ***********************************************************************
 * When working with power law curves, we can fit a straight
 * line by using log log coordinates.  
 * originally used natural log, need to use log base 10.  6-14-02  
 * This is hard coded using mode=5.
 * The input is expected in  X  Y format so no switching needed.
 * Do NOT normalize the X value because this screws up y-intercept values. 
 * Before this code is run, 
 * need to get the two column input file using awk on the .average file 
 * TO COMPILE:  gcc -lm fit_log.c then mv a.out fit_log
 * TO RUN:  ./fit_log dsj1000_gen2.4021.45
 * The output is  X  with fitted Y coordinates.  First line of outfile has
 * a (y-intercept) and b (slope).  outfile gets overwritten so remember to:
 * mv outfile dsj1000_gen2.4021.log_ab 
 *************************************************************************/

      /********************************************************************
       * a is the y intercept and b is the slope for the straight line    *
       *             y = a +(b * x)                                       *
       *             log10(y) = a + [ b *log10(x)]  this is log base 10   *
       * returns fitted y: y=pow( 10, a+b*log(x) )  10 raised to the power*
       ********************************************************************
       * see the modeling of Data Chapter 15 in Numerical Recipes      *
       * the Chi Square fitting assuuming sigma i is 1.                *
       * The x1 is Sx; x2 is Sxx; y1 is Sy; y2 is Sxy   page 662       * 
       *              x1+=log10(x[i]);                                 *
       *              x2+=log10(x[i])*log10(x[i]);                     *
       *              y1+=log10(y[i]);                                 *
       *              y2+=log(10x[i])*log10(y[i]);                     *
       * d is delta, d1 and d2 are used for formula 15.2.6 page 662    * 
       *****************************************************************/

double det3(double x11, double x12, double x13,
	    double x21, double x22, double x23,
	    double x31, double x32, double x33);
double det4(double x11, double x12, double x13, double x14,
	    double x21, double x22, double x23, double x24,
	    double x31, double x32, double x33, double x34,
	    double x41, double x42, double x43, double x44);
int main(argc, argv)
int argc;
char **argv;
{
char            s[80];
FILE            *infile,*outptr;  
 int             i, j,m,n1,n2;
int             norm=0;	    /*LG 06-02 default do NOT normalize */

int             switch_col = 0; /*LG 03-02 input X Y default do NOT switch_col */
int             mode=5;    /*LG 03-02 default mode =5 */    
                            /* 1-parabolic */
                              /* 2-exponential */
double           x[1000], y[1000],min;
double          x1,x2,x3,x4,x5,x6,x7,x8,
		y1,y2,y3,y4,y5,d,d1,d2,d3,d4,d5;
double          a,a0,a1,a2,a3,a4,b,beta;
char            fname[80];
char            str[80];

            x1=0.0; x2=0.0; x3=0.0; x4=0.0;
	    x5=0.0; x6=0.0; x7=0.0; x8=0.0;
            y1=0.0; y2=0.0; y3=0.0; y4=0.0; y5=0.0;
            m = 0;
/********************************************
 * use argv to get parameters    
 ********************************************/
	if (argc != 2)
	{
	printf ("Usage: fit <input file> <mode 5=log10log10 hardcoded> \n");
 	exit (1);
	}

	infile = fopen (argv[1], "r");
	if (infile == NULL)
	{
		printf ("FIT: Error in opening input file\n");
		exit (1);
	}

	outptr = fopen ("outfile", "w");
	if (outptr == NULL)
		{
		printf ("FIT: Error in  opening output file\n");
		exit (1);
	}

     
	/*****************************************************
         commented out for now...if ever need these options 
         *****************************************************
            printf("\nEnter mode: ");
            printf("\n\t1-parabolic");
            printf("\n\t2-exponential");
            printf("\n\t3-logarithmic");
            printf("\n\t4-quadratic log-log");
            printf("\n\t  log(y) = a0 + a1*log(x) + a2*[log(x)]^2");
            printf("\n\t5-linear log10-log10");
            printf("\n\t  log(y) = a0 + a1*log(x)");
            printf("\n\t6-linear");
            printf("\n\t7-quadratic w/o linear term");
            printf("\n\t8-quadratic a0=0, a1=1");
            printf("\n\t  y = x - a2*x^2");
	    printf("\n\t9-y = x - a(x^b)");
	    printf("\n\t10- y-x = a(x^b)");
	    printf("\n\t11-cubic");
	    printf("\n\t12-4th deg poly");
	    printf("\n\t13- xy=c");
	    printf("\n\t14- y = ax + b/x");
	    printf("\n\t");
            scanf("%d",&mode);
         *****************************************************
	 ****force mode =5 **********
         ********************************************************/
	    if(mode==9 || mode==10) {
	       printf("\n\tInput beta: ");
	       scanf("%lf",&beta);
	    }
	    /*****************************************************
	    ****default norm =0 *********
            printf("\nNormalize x?(1-yes, 0-no): ");
            scanf("%d",&norm);
            printf("\nSwitch x, y column?(1-yes, 0-no): ");
            scanf("%d",&switch_col);
	    ****force switch_col =0 *********
            ******************************************************/
            printf("\nThe average is in \"outfile\"\n");
 
            while(NULL != fgets(s,80,infile)) {
                  i = m++;
                   if(switch_col != 1) {
                   sscanf(s,"%lf %lf\n",&x[i],&y[i]);
                   }
                   if(switch_col == 1) {
                     sscanf(s,"%lf %lf\n",&y[i],&x[i]);
                   }
		   /*
                 printf("%d %g %g\n",m,x[i],y[i]);
		 */
                   if(norm == 1) {
                     if(i==0) min=x[0];
                     x[i]=x[i]/min;
                   }
                   if(mode==1 || mode==11 || mode==12) {
                     x1+=x[i];
                     x2+=x[i]*x[i];
                     x3+=x[i]*x[i]*x[i];
                     x4+=x[i]*x[i]*x[i]*x[i];  
                     y1+=y[i];
                     y2+=y[i]*x[i];
                     y3+=y[i]*x[i]*x[i];
                   }
                   if(mode==2) {
                     x1+=x[i];
                     x2+=x[i]*x[i];
                     y1+=log(y[i]);
                     y2+=x[i]*log(y[i]);
                   }
                   if(mode==3) {
                     x1+=log(x[i]);
                     x2+=log(x[i])*log(x[i]);
                     y1+=y[i];
                     y2+=log(x[i])*y[i];
                   }
                   if(mode==4) {
                     x1+=log(x[i]);
                     x2+=log(x[i])*log(x[i]);
                     x3+=log(x[i])*log(x[i])*log(x[i]);
                     x4+=log(x[i])*log(x[i])*log(x[i])*log(x[i]);
                     y1+=log(y[i]);
                     y2+=log(x[i])*log(y[i]);
                     y3+=log(x[i])*log(x[i])*log(y[i]);
                   }
                   if(mode==5) {
                     x1+=log10(x[i]);
                     x2+=log10(x[i])*log10(x[i]);
                     y1+=log10(y[i]);
                     y2+=log10(x[i])*log10(y[i]);
                   }
                   if(mode==6) {
                     x1+=x[i];
                     x2+=x[i]*x[i];
                     y1+=y[i];
                     y2+=y[i]*x[i];
                   }
		   if(mode==7) {
		     x1+=x[i]*x[i];
		     x2+=x[i]*x[i]*x[i]*x[i];
		     y1+=y[i];
		     y2+=y[i]*x[i]*x[i];
		   }
		   if(mode==8) {
		     x1+=x[i];
		     x2+=x[i]*x[i];
		     y1+=y[i];
		   }
		   if(mode==9) {
		     x1+=x[i];
		     x2+=pow(x[i],beta);
		     y1+=y[i];
		   }
		   if(mode==10) {
		     x2+=pow(x[i],beta);
		     y1+=(y[i]-m);
		     printf("%d %lf\n",m,y[i]-m);
		   }
		   if(mode==11 || mode==12) {
		     x5+=x[i]*x[i]*x[i]*x[i]*x[i];
		     x6+=x[i]*x[i]*x[i]*x[i]*x[i]*x[i];
		     y4+=x[i]*x[i]*x[i]*y[i];
		   }
		   if(mode==12) {
		     x7+=x[i]*x[i]*x[i]*x[i]*x[i]*x[i]*x[i];
		     x8+=x[i]*x[i]*x[i]*x[i]*x[i]*x[i]*x[i]*x[i];
		     y5+=x[i]*x[i]*x[i]*x[i]*y[i];
		   }
		   if(mode==13) {
		     x1+=x[i];
		     y1+=y[i];
		   }
		   if(mode==14) {
		     x2+=x[i]*x[i];
		     x3+= 1.0 / (x[i]*x[i]);
		     y2+= x[i]*y[i];
		     y3+= y[i] / x[i];
		   }

                }

      if( mode == 1 || mode == 4) {
      d= m* (x2*x4-x3*x3)-x1*(x1*x4-x2*x3)+x2*(x1*x3-x2*x2);
      d1=y1*(x2*x4-x3*x3)-x1*(y2*x4-y3*x3)+x2*(y2*x3-y3*x2);
      d2=m* (y2*x4-y3*x3)-y1*(x1*x4-x2*x3)+x2*(x1*y3-x2*y2);
      d3=m* (x2*y3-x3*y2)-x1*(x1*y3-x2*y2)+y1*(x1*x3-x2*x2);
      a0=d1/d;
      a1=d2/d;
      a2=d3/d;
printf("a0=%g a1=%g a2=%g\n",a0,a1,a2);
fprintf(outptr,"TitleText: a0=%g a1=%g a2=%g\n",a0,a1,a2);
      }
      if ( mode==2 || mode==3 || mode==5 || mode==6 || mode==7) {
      printf("4\n");
      /*****************************************************************
       * see the modeling of Data Chapter 15 in Numerical Recipes      *
       * d is delta, d1 and d2 are used for formula 15.2.6 page 662    * 
       * m is the number of samples and sigma i is assumed to be 1     *
       *****************************************************************/
         d = m*x2 - x1*x1;
         d1 = y1*x2 - y2*x1;
         d2 = m*y2 - x1*y1;
         if(mode==2) a = exp(d1/d);
         if(mode==3 || mode==5 || mode==6 || mode==7) a = d1 / d;
         b = d2/d;
       /*****************************************************************
       * a is the y intercept and b is the slope for                    *
       * the best fit straight line: y = a +(b * x)                     *
       *****************************************************************/

printf("x1=%lf x2=%lf y1=%lf y2=%lf\n",x1,x2,y1,y2);
printf("d=%lf d1=%lf d2=%lf\n",d,d1,d2);
printf("a=%g b=%lf\n",a,b);
fprintf(outptr,"fit coordinates: a   %lf    b   %lf\n",a,b);
      }
      if (mode == 8 || mode == 9) {
	 b=(x1-y1)/x2;
	 fprintf(outptr,"TitleText: b=%lf\n",b);
      }
      if (mode == 10) {
	 b=y1/x2;
	 fprintf(outptr,"TitleText: b=%lf\n",b);
      }
      if (mode == 11) {
d = m * det3(x2,x3,x4,x3,x4,x5,x4,x5,x6) -
    x1 * det3(x1,x3,x4,x2,x4,x5,x3,x5,x6) +
    x2 * det3(x1,x2,x4,x2,x3,x5,x3,x4,x6) -
    x3 * det3(x1,x2,x3,x2,x3,x4,x3,x4,x5);
d1 = y1 * det3(x2,x3,x4,x3,x4,x5,x4,x5,x6) -
     x1 * det3(y2,x3,x4,y3,x4,x5,y4,x5,x6) +
     x2 * det3(y2,x2,x4,y3,x3,x5,y4,x4,x6) -
     x3 * det3(y2,x2,x3,y3,x3,x4,y4,x4,x5);
d2 = m * det3(y2,x3,x4,y3,x4,x5,y4,x5,x6) -
     y1 * det3(x1,x3,x4,x2,x4,x5,x3,x5,x6) +
     x2 * det3(x1,y2,x4,x2,y3,x5,x3,y4,x6) -
     x3 * det3(x1,y2,x3,x2,y3,x4,x3,y4,x5);
d3 = m * det3(x2,y2,x4,x3,y3,x5,x4,y4,x6) -
     x1 * det3(x1,y2,x4,x2,y3,x5,x3,y4,x6) +
     y1 * det3(x1,x2,x4,x2,x3,x5,x3,x4,x6) -
     x3 * det3(x1,x2,y2,x2,x3,y3,x3,x4,y4);
d4 = m * det3(x2,x3,y2,x3,x4,y3,x4,x5,y4) -
     x1 * det3(x1,x3,y2,x2,x4,y3,x3,x5,y4) +
     x2 * det3(x1,x2,y2,x2,x3,y3,x3,x4,y4) -
     y1 * det3(x1,x2,x3,x2,x3,x4,x3,x4,x5);
     a0 = d1/d;
     a1 = d2/d;
     a2 = d3/d;
     a3 = d4/d;
     a4 = d5/d;
     fprintf(outptr,"TitleText: a=%lf b=%lf c=%lf d=%lf\n",
	     a0,a1,a2,a3);
}
      if (mode == 12) {
d =  m * det4(x2,x3,x4,x5, x3,x4,x5,x6, x4,x5,x6,x7, x5,x6,x7,x8) -
    x1 * det4(x1,x3,x4,x5, x2,x4,x5,x6, x3,x5,x6,x7, x4,x6,x7,x8) +
    x2 * det4(x1,x2,x4,x5, x2,x3,x5,x6, x3,x4,x6,x7, x4,x5,x7,x8) -
    x3 * det4(x1,x2,x3,x5, x2,x3,x4,x6, x3,x4,x5,x7, x4,x5,x6,x8) +
    x4 * det4(x1,x2,x3,x4, x2,x3,x4,x5, x3,x4,x5,x6, x4,x5,x6,x7);
d1 = y1 *det4(x2,x3,x4,x5, x3,x4,x5,x6, x4,x5,x6,x7, x5,x6,x7,x8) -
     x1 *det4(y2,x3,x4,x5, y3,x4,x5,x6, y4,x5,x6,x7, y5,x6,x7,x8) +
     x2 *det4(y2,x2,x4,x5, y3,x3,x5,x6, y4,x4,x6,x7, y5,x5,x7,x8) -
     x3 *det4(y2,x2,x3,x5, y3,x3,x4,x6, y4,x4,x5,x7, y5,x5,x6,x8) +
     x4 *det4(y2,x2,x3,x4, y3,x3,x4,x5, y4,x4,x5,x6, y5,x5,x6,x7);
d2 = m * det4(y2,x3,x4,x5, y3,x4,x5,x6, y4,x5,x6,x7, y5,x6,x7,x8) -
     y1 *det4(x1,x3,x4,x5, x2,x4,x5,x6, x3,x5,x6,x7, x4,x6,x7,x8) +
     x2 *det4(x1,y2,x4,x5, x2,y3,x5,x6, x3,y4,x6,x7, x4,y5,x7,x8) -
     x3 *det4(x1,y2,x3,x5, x2,y3,x4,x6, x3,y4,x5,x7, x4,y5,x6,x8) +
     x4 *det4(x1,y2,x3,x4, x2,y3,x4,x5, x3,y4,x5,x6, x4,y5,x6,x7);
d3 = m * det4(x2,y2,x4,x5, x3,y3,x5,x6, x4,y4,x6,x7, x5,y5,x7,x8) -
     x1 *det4(x1,y2,x4,x5, x2,y3,x5,x6, x3,y4,x6,x7, x4,y5,x7,x8) +
     y1 *det4(x1,x2,x4,x5, x2,x3,x5,x6, x3,x4,x6,x7, x4,x5,x7,x8) -
     x3 *det4(x1,x2,y2,x5, x2,x3,y3,x6, x3,x4,y4,x7, x4,x5,y5,x8) +
     x4 *det4(x1,x2,y2,x4, x2,x3,y3,x5, x3,x4,y4,x6, x4,x5,y5,x7);
d4 = m * det4(x2,x3,y2,x5, x3,x4,y3,x6, x4,x5,y4,x7, x5,x6,y5,x8) -
     x1 *det4(x1,x3,y2,x5, x2,x4,y3,x6, x3,x5,y4,x7, x4,x6,y5,x8) +
     x2 *det4(x1,x2,y2,x5, x2,x3,y3,x6, x3,x4,y4,x7, x4,x5,y5,x8) -
     y1 *det4(x1,x2,x3,x5, x2,x3,x4,x6, x3,x4,x5,x7, x4,x5,x6,x8) +
     x4 *det4(x1,x2,x3,y2, x2,x3,x4,y3, x3,x4,x5,y4, x4,x5,x6,y5);
d5 = m * det4(x2,x3,x4,y2, x3,x4,x5,y3, x4,x5,x6,y4, x5,x6,x7,y5) -
     x1 *det4(x1,x3,x4,y2, x2,x4,x5,y3, x3,x5,x6,y4, x4,x6,x7,y5) +
     x2 *det4(x1,x2,x4,y2, x2,x3,x5,y3, x3,x4,x6,y4, x4,x5,x7,y5) -
     x3 *det4(x1,x2,x3,y2, x2,x3,x4,y3, x3,x4,x5,y4, x4,x5,x6,y5) +
     y1 *det4(x1,x2,x3,x4, x2,x3,x4,x5, x3,x4,x5,x6, x4,x5,x6,x7);
     a0 = d1/d;
     a1 = d2/d;
     a2 = d3/d;
     a3 = d4/d;
     a4 = d5/d;
     fprintf(outptr,"TitleText: a=%g b=%g c=%g d=%g e=%g\n",
             a0,a1,a2,a3,a4);
}
      if (mode == 13) {
        a = y1/x1;
	fprintf(outptr,"TitleText: a=%lf\n",a);
      }
      if (mode == 14) {
	d = x2*x3 - m*m;
	d1 = y2*x3 - m*y3;
	d2 = x2*y3 - m*y2;
	a0 = d1/d;
	a1 = d2/d;
	fprintf(outptr,"TitleText: a=%lf b=%lf\n",a0,a1);
      }
      
      for ( j=0; j<m; j++ ) {
          if(mode==1) y[j]=a0+a1*x[j]+a2*x[j]*x[j];
          if(mode==2) y[j]=a*exp(b*x[j]);
          if(mode==3) y[j]=a + b * log(x[j]);
          if(mode==4) y[j]=exp( a0+a1*log(x[j])+a2*log(x[j])*log(x[j]) );
          if(mode==5) y[j]=pow( 10, a+b*log10(x[j]) );
          if(mode==6) y[j]=a + b * x[j];
          if(mode==7) y[j]=a + b * x[j]*x[j];
          if(mode==8) y[j]=x[j] - b * x[j]*x[j];
          if(mode==9) y[j]=x[j] - b * pow(x[j],beta);
          if(mode==10) y[j]= x[j] + b * pow(x[j],beta);
	  if(mode==11) y[j]= a0+a1*x[j]+a2*x[j]*x[j]+
			     a3*x[j]*x[j]*x[j];
          if(mode==12) y[j]= a0+a1*x[j]+a2*x[j]*x[j]+
			     a3*x[j]*x[j]*x[j]+
			     a4*x[j]*x[j]*x[j]*x[j];
          if(mode==13) y[j] = a / x[j];
	  if(mode==14) y[j] = a0*x[j] + a1/x[j];
          if(norm==1) x[j]*=min;

	  /*************************************************
	   * need this in output file
	   *************************************************/
          if(switch_col != 1) {
            fprintf(outptr, "%g %g\n",x[j],y[j]);
          }
          if(switch_col == 1) {
            fprintf(outptr, "%g %g\n",y[j],x[j]);
          }
	   /*************************************************
	   * need this in output file
	   *************************************************/
          

      } /* end for j loop */
             fclose(infile);
             fclose(outptr);

}

double det3(double x11, double x12, double x13,
	   double x21, double x22, double x23,
           double x31, double x32, double x33)

{
   double y;

   y = x11*(x22*x33-x32*x23) -
       x12*(x21*x33-x31*x23) +
       x13*(x21*x32-x31*x22);
   return y;
}

double det4(double x11, double x12, double x13, double x14,
	    double x21, double x22, double x23, double x24,
            double x31, double x32, double x33, double x34,
	    double x41, double x42, double x43, double x44)
{
   double y;

   y = x11*det3( x22,x23,x24, x32,x33,x34, x42,x43,x44) -
       x12*det3( x21,x23,x24, x31,x33,x34, x41,x43,x44) +
       x13*det3( x21,x22,x24, x31,x32,x34, x41,x42,x44) -
       x14*det3( x21,x22,x23, x31,x32,x33, x41,x42,x43);
   return y;
}










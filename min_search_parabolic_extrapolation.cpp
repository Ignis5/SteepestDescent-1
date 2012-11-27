////////////////////////////////////////////////////////////////////////////////
// File: min_search_parabolic_extrapolation.c                                 //
// Routine(s):                                                                //
//    Min_Search_Parabolic_Extrapolation                                      //
////////////////////////////////////////////////////////////////////////////////

#include <float.h>                        // required for DBL_MAX, DBL_EPSILON
#include <math.h>                         // required for fabs(), sqrt()

#define SQRT5 2.236067977499789696

static const double alpha = 0.5 * (3.0 - SQRT5);
static double min_displacement;

//                        Internally Defined Routines 

static double Parabolic_Minimum( double x[], double y[],
                                                double cut_off_scale_factor ); 
static double New_Test_Point( double x[], double y[], double a, double b, 
                                                double cut_off_scale_factor );

////////////////////////////////////////////////////////////////////////////////
//  int Min_Search_Parabolic_Extrapolation(double (*f)(double), double x[],   //
//            double fx[], double a, double b,  double cut_off_scale_factor)  //
//                                                                            //
//  Description:                                                              //
//     This routine can be used to bracket the minimum of a unimodal function //
//     or a local minimum of a multimodal function within a preassigned       //
//     closed interval.                                                       //
//                                                                            //
//     Given points (x[0], f(x[0])) and (x[2], f(x[2])) with                  //
//     a <= x[0] < x[2] <= b, the parabolic extrapolation method brackets a   //
//     local minimum of f(x).                                                 //
//                                                                            //
//  Arguments:                                                                //
//     double (*f)(double)                                                    //
//        Pointer to a user-defined function of a real variable (type double) //
//        returning a real number (type double).                              //
//     double x[]                                                             //
//        On input, x[0] and x[2] are defined as two initial starting points  //
//        where a < x[0], x[2] < b.  On output, x[0] is the lower endpoint of //
//        the final subinterval which contains the local minimum, x[2] is the //
//        upper endpoint and x[1] is a point between x[0] and x[2].           //
//     double fx[]                                                            //
//        On input, fx[0] is the value of f(x) evaluated at x[0] and fx[2] is //
//        the value of f(x) evaluated at x[2].  On output, fx[i] is the value //
//        of f(x) evaluated at the final value of x[i], i = 0,1,2.            //
//     double a                                                               //
//        The lower bound of the interval on which the function f is defined. //
//        If a is -infinity, then set a = -DBL_MAX where DBL_MAX is defined   //
//        in <float.h>.                                                       //
//     double b                                                               //
//        The upper bound of the interval on which the function f is defined. //
//        If b is infinity, then set b = DBL_MAX where DBL_MAX is defined in  //
//        <float.h>.                                                          //
//     double cutoff_scale_factor                                             //
//        A parameter which limits the displacement in any single step.       //
//                                                                            //
//  Return Values:                                                            //
//      0  Success                                                            //
//     -1  Process failed to find a value x such that a local minimum is      //
//         found in the interval [a, b], i.e. the minimum is probably located //
//         on the boundary.                                                   //
//     -2  Process failed, the initial values of x[0] and x[2] are equal.     //
//                                                                            //
//  Example:                                                                  //
//     #include <float.h>                                                     //
//     extern double f(double);                                               //
//     double a = -DBL_MAX;                                                   //
//     double b = DBL_MAX;                                                    //
//     double x[3];                                                           //
//     double fx[3];                                                          //
//     double cut_off_scale_factor = 10.0;                                    //
//     int err;                                                               //
//                                                                            //
//     (set x[0], x[2] )                                                      //
//                                                                            //
//     fx[0] = f(x[0]);                                                       //
//     fx[2] = f(x[2]);                                                       //
//                                                                            //
//     err = Min_Search_Parabolic_Extrapolation(f, x, fx, a, b,               //
//                                                     cut_off_scale_factor)  //
//     if (err < 0 ) printf("No local minimum found\n");                      //
//     else printf("local minimum is between (%12.6e,%12.6e) and              //
//                                (%12.6e,%12.6e)\n", x[0],fx[0],x[2],fx[2]); //
////////////////////////////////////////////////////////////////////////////////

int Min_Search_Parabolic_Extrapolation(double (*f)(double), double x[],
            double fx[], double a, double b,  double cut_off_scale_factor)
{
   double xmin, ymin;
   int return_code = 0;

   min_displacement = sqrt(DBL_EPSILON);

                          // Verify that x[0] < x[2]

   if ( x[0] == x[2] ) return -2;
   if ( x[0] > x[2] ) {
      x[1] = x[0]; fx[1] = fx[0]; 
      x[0] = x[2]; fx[0] = fx[2];
      x[2] = x[1]; fx[2] = fx[1];
   }

                   // Set x[1] such that x[0] < x[1] < x[2]

   if ( fx[0] < fx[2] ) x[1] = x[0] + alpha * (x[2] - x[0]);
      x[1] = x[2] - alpha * (x[2] - x[0]);
   fx[1] = f(x[1]);

              // Iterate until either f(x[0]) > f(x[1]) < f(x[2])
              // or x[0] == a or x[2] == b.

   while ( ( fx[0] <= fx[1] ) || ( fx[2] <= fx[1] ) ) {
      xmin = New_Test_Point( x, fx, a, b, cut_off_scale_factor );
      if ( (xmin == a) || (xmin == b) ) { return_code = -1; break; }
      ymin = f(xmin);
      if ( xmin < x[0] ) {
         x[2] = x[1]; fx[2] = fx[1];
         x[1] = x[0]; fx[1] = fx[0];
         x[0] = xmin; fx[0] = ymin;
      } 
      else if ( xmin > x[2] ) {
         x[0] = x[1]; fx[0] = fx[1];
         x[1] = x[2]; fx[1] = fx[2];
         x[2] = xmin; fx[2] = ymin;
      }
      else if ( xmin < x[1] )
         if  (( ymin >= fx[0] ) && ( ymin <= fx[1] )) {
            x[2] = x[1]; fx[2] = fx[1];
            x[1] = xmin; fx[1] = ymin;
         }
         else if ( ymin > fx[1] ) { x[0] = xmin; fx[0] = ymin; break; }
         else { x[2] = x[1]; fx[2] = fx[1]; x[1] = xmin; fx[1] = ymin; break; }
      else 
         if  (( ymin >= fx[2] ) && ( ymin <= fx[1] )) {
            x[0] = x[1]; fx[0] = fx[1];
            x[1] = xmin; fx[1] = ymin;
         }
         else if ( ymin > fx[1] ) { x[2] = xmin; fx[2] = ymin; break; }
         else if ( ymin < fx[2] ) { x[1] = xmin; fx[1] = ymin; break; }
   }

   return return_code;
}


////////////////////////////////////////////////////////////////////////////////
//  static double New_Test_Point( double x[], double y[], double a, double b, //
//                     double min_displacement, double cut_off_scale_factor ) //
//                                                                            //
//  Description:                                                              //
//     Given three points (x[0],y[0]), (x[1],y[1]), and (x[2],y[2]) such that //
//     x[0] < x[1] < x[2] find the location of the next test point, normally  //
//     the minimum, xmin, of the parabola through these points.               //
//                                                                            //
//     For C compilers which support the 'inline' type qualifier, this routine//
//     could be declared with the 'inline' qualifier.                         //
//                                                                            //
//  Return Values:                                                            //
//     The next test point.                                                   //
////////////////////////////////////////////////////////////////////////////////
static double New_Test_Point( double x[], double y[], double a, double b, 
                                                 double cut_off_scale_factor ) 
{
   double xmin;
   double eta;

   xmin = Parabolic_Minimum( x, y, cut_off_scale_factor );

                 // Compare the location of the minimum of the
                 // parabolic fit to the left-most node.

   eta = ( x[0] == 0.0 ) ? 1.0 : fabs(x[0]);
   if ( fabs(xmin - x[0]) < (min_displacement * eta) ) 
      xmin -= 100.0 * min_displacement * eta;
   else {
                 // Compare the location of the minimum of the
                 // parabolic fit to the center node.

      eta = ( x[1] == 0.0 ) ? 1.0 : fabs(x[1]);
      if ( fabs(xmin - x[1]) < (min_displacement * eta) ) 
         if ( y[2] < y[0] )  
            xmin = x[1] + alpha * (x[2] - x[1]);
         else 
            xmin = x[1] - alpha * (x[1] - x[0]);
      else {
                 // Compare the location of the minimum of the
                 // parabolic fit to the right-most node.

         eta = ( x[2] == 0.0 ) ? 1.0 : fabs(x[2]);
         if ( fabs(x[2] - xmin) < (min_displacement * eta) ) 
            xmin += 100.0 * min_displacement * eta;
      }
   }

   if (xmin < a) xmin = a;
   if (xmin > b) xmin = b;

   return xmin;
}


////////////////////////////////////////////////////////////////////////////////
//  static double Parabolic_Minimum( double x[], double y[],                  //
//                                            double cut_off_scale_factor )   //
//                                                                            //
//  Description:                                                              //
//     Given three points (x[0],y[0]), (x[1],y[1]), and (x[2],y[2]) such that //
//     x[0] < x[1] < x[2] and y[1] > the line through (x[0],y[0]) and         //
//     (x[2],y[2]) at x[1] find the location of the minimum, xmin, of the     //
//      parabola through these points.                                        //
//                                                                            //
//     For C compilers which support the 'inline' type qualifier, this routine//
//     could be declared with the 'inline' qualifier.                         //
//                                                                            //
//  Return Values:                                                            //
//     The location of the minimum.                                           //
////////////////////////////////////////////////////////////////////////////////
static double Parabolic_Minimum( double x[], double y[],
                                                  double cut_off_scale_factor ) 
{
   double d1 = (x[2] - x[1]) * (y[1] - y[0]);
   double d2 = (x[1] - x[0]) * (y[2] - y[1]);
   double denominator = 2.0 * (d1 - d2);
   double numerator = (x[2] + x[1]) * d1 - (x[1] + x[0]) * d2;
   double x_lb = x[0] - cut_off_scale_factor * (x[1] - x[0]);
   double x_ub = x[2] + cut_off_scale_factor * (x[2] - x[1]);
   double ymax =( y[2] - y[0]) * (x[1] - x[0]) / (x[2] - x[0]) + y[0];
   double xmin;

   if (y[1] >= ymax) xmin = ( y[0] < y[2] ) ? x_lb : x_ub;
   else {
      if ( x_lb * denominator <= numerator ) xmin = x_lb;
      else if ( x_ub * denominator >= numerator ) xmin = x_ub;
      else  xmin = numerator / denominator;
   }
   return xmin;
}

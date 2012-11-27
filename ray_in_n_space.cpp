////////////////////////////////////////////////////////////////////////////////
// File: ray_in_n_space.c                                                     //
// Routines:                                                                  //
//    Ray_in_n_Space                                                          //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Ray_in_n_Space(double v[], double initial_pt[], double direction[],  //
//                                                    double distance, int n) //
//                                                                            //
//  Description:                                                              //
//     Calculates the vector v = initial_pt + distance * direction.           //
//     If direction has unit norm, the v lies at a distance |distance| from   //
//     the initial_pt, i.e. || v - initial_pt || = | distance |.              //
//     In general, || v - initial_pt || = |distance| || direction ||.         //
//     The argument distance can be positive or negative.  The argument       //
//     direction need not have unit norm.                                     //
//                                                                            //
//  Arguments:                                                                //
//     double v[]                                                             //
//            Pointer to the first element of the resultant vector v[n].      //
//     double initial_pt[]                                                    //
//            Pointer to the first element of the vector initial_pt[n].       //
//     double direction[]                                                     //
//            Pointer to the first element of the vector direction[n].        //
//     double distance                                                        //
//            Distance to travel from the initial point, initial_pt, along    //
//            a line with direction 'direction', assuming that 'direction'    //
//            has unit norm.                                                  //
//     int    n                                                               //
//            The dimension of the vectors, v, initial_pt, and direction.     //
//                                                                            //
//  Return Values:                                                            //
//     address to v.                                                          //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double v[N], init_pt[N], d[N], distance;                               //
//     int i;                                                                 //
//                                                                            //
//     (your code to intialize the vectors init_pt and d, and the scalar)     //
//     (distance.                                                       )     //
//                                                                            //
//     Ray_in_n_Space(v, init_pt, d, distance, N);                            //
//     for (i = 0; i < N; i++)                                                //
//        printf(" v[%3d] = %12.6f\n", i,v[i]);                               //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //

void Ray_in_n_Space(double v[], double initial_pt[], double direction[],
                                                      double distance, int n )
{
   for (--n; n >= 0; n--) 
     v[n]  = initial_pt[n] + distance * direction[n];
}

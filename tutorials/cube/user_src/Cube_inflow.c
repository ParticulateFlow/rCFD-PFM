#include "udf.h"

/***********************************************************************
  vprofile.c
  UDF for specifying a steady-state inlet profile boundary condition for cubic dispersion
  power-law; k; epsilon.
  (C) Yaxing Du (2020)
************************************************************************/


DEFINE_PROFILE(inlet_x_v,thread,index)
{
 
 real x[ND_ND];
 real z;
 face_t f;

 begin_f_loop(f,thread)

  {
   
   F_CENTROID(x,f,thread);
   z = x[2];
   F_PROFILE(f,thread,index) = (double) 5.83* pow(z,0.19);

  }

 end_f_loop(f,thread)

}

DEFINE_PROFILE(inlet_x_k,thread,index)
{
 
 real x[ND_ND];
 real z;

 face_t f;

 begin_f_loop(f,thread)

  {
   
   F_CENTROID(x,f,thread);
   z = x[2];
   
   F_PROFILE(f,thread,index) = (double) 0.202* pow(z,0.0004);

  }

 end_f_loop(f,thread)

}

DEFINE_PROFILE(inlet_x_dissp,thread,index)
{
 
 real x[ND_ND];
 real z;

 face_t f;

 begin_f_loop(f,thread)

  {
   
   F_CENTROID(x,f,thread);
   z = x[2];
   
   F_PROFILE(f,thread,index) = (double) 2.585* pow(z,0.00006);

  }

 end_f_loop(f,thread)

}

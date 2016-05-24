//
//  KaspyCycler.hpp
//  kaspy_cuda
//
//  Created by Andrei Koulikov on 24.05.16.
//
//

#ifndef KaspyCycler_hpp
#define KaspyCycler_hpp

#include <stdio.h>

#include "fortran_vars.h"





class KaspyCycler
{
public:
    
    KaspyCycler(fortran_common_vars * for_vars, fortran_common_arrays * for_arrays,
                fortran_ffloats * for_floats, fortran_wind_data * for_wind_data,
                float * for_press, float * for_press0,
                float * for_uwd, float * for_uwd0, float * for_vwd, float * for_vwd0) :
    m_width(F_DATA_WIDTH),
    m_height(F_DATA_HEIGHT),
    m_fVars(for_vars),
    m_fArrays(for_arrays),
    m_fFloats(for_floats),
    m_fWindData(for_wind_data),
    itime6(0),
    itime6_old(0),
    m_press(for_press),
    m_press0(for_press0),
    m_uwd(for_uwd),
    m_uwd0(for_uwd0),
    m_vwd(for_vwd),
    m_vwd0(for_vwd0)
    {
        //ersetb();
        printf("cycler alloc with dht %f, array marker %32.32f \n", m_fVars->dht, m_fArrays->marker);
        printf("end marker %32.32f \n", m_fArrays->end_marker);
        printf("ff end %32.32f \n", m_fFloats->end_marker);
        printf("wind data: \n \t %d \t %d \t %d \n \t %d \t %d \t %d \n \t %d \t %d \t %d \n",
               m_fWindData->kx, m_fWindData->ky, m_fWindData->kt,
               m_fWindData->kxu, m_fWindData->kyu, m_fWindData->ktu,
               m_fWindData->kxv, m_fWindData->kyv, m_fWindData->ktv);
    }
    
 
    void loadData();

    
    void findElves();
    
    void makeWsurf(float ro_ratio);
    
    
    
    int m_width;
    int m_height;
    
    fortran_common_vars * m_fVars;
    fortran_common_arrays * m_fArrays;
    fortran_ffloats * m_fFloats;
    fortran_wind_data * m_fWindData;
    
    float * m_press;
    float * m_press0;
    float * m_uwd;
    float * m_uwd0;
    float * m_vwd;
    float * m_vwd0;
    
    
private:
    float ftim;
    float btim;
    int itime6;
    int itime6_old;
    
    
    float * g_fbu;
    float * g_fbv;
    float * g_ffu;
    float * g_ffv;
    
    float * g_fxb;
    float * g_fxf;
    float * g_fyb;
    float * g_fyf;
    
    float * g_fb;
    float * g_ff;
    
    float * g_wusurf;
    float * g_wvsurf;
    
    float * g_dum;
    float * g_dvm;
    
    float * g_d;
};





#endif /* KaspyCycler_hpp */

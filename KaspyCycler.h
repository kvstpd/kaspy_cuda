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
    
    KaspyCycler(fortran_common_vars * for_vars, fortran_common_arrays * for_arrays, fortran_ffloats * for_floats) :
    m_width(F_DATA_WIDTH),
    m_height(F_DATA_HEIGHT),
    m_fVars(for_vars),
    m_fArrays(for_arrays),
    m_fFloats(for_floats)
    {
        //ersetb();
        printf("cycler alloc with var marker %f, array marker %32.32f \n", m_fVars->marker, m_fArrays->marker);
        printf("end marker %32.32f \n", m_fArrays->end_marker);
        printf("ff end %32.32f \n", m_fFloats->end_marker);
    }
    
    
    void findElves();
    
    void makeWsurf(float ro_ratio);
    
    
    
    int m_width;
    int m_height;
    
    fortran_common_vars * m_fVars;
    fortran_common_arrays * m_fArrays;
    fortran_ffloats * m_fFloats;
};





#endif /* KaspyCycler_hpp */


#include <stdio.h>
// #include <stdlib.h>

#ifdef _WIN64
    #include <windows.h>
#else
    #include <sys/time.h>
#endif

#include <math.h>

#include "fortran_vars.h"


void ersetb(void)
{
    setbuf(stdout,NULL); /* set output to unbuffered */
}


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
        //printf("cycler alloc with var marker %f, array marker %32.32f \n", m_fVars->marker, m_fArrays->marker);
        //printf("end marker %32.32f \n", m_fArrays->end_marker);
        //printf("ff end %32.32f \n", m_fFloats->end_marker);
    }
    
    
    void findElves()
    {
        //printf("arrays is set to %llxd \n", (long long)m_fArrays );

        
        float * elves = &(m_fArrays->elf[0][0]);
        
        float elf_min = elves[0];
        float elf_max = elves[0];
         
        for (int i=1; i<F_DATA_SIZE; i++)
        {
            if (elves[i] > elf_max)
            {
                elf_max = elves[i];
            }
            
            if (elves[i] < elf_min)
            {
                elf_min = elves[i];
            }
        }
        
        //printf("C SAYS: time is %f, elf min is %f, elf max is %f \n",m_fVars->timeh, elf_min, elf_max);
    }

   
/* 
 DO J=2,JMM1
 DO I=2,IMM1
	uw=(btim*fbu(i,j)+ftim*ffu(i,j))
	vw=(btim*fbv(i,j)+ftim*ffv(i,j))
	speed=sqrt(uw**2+vw**2) !******************************************************
 !      speed=0
	windc=1.0e-3*(0.8+speed*0.065)*ro_ratio*speed
 WUSURF(I,J)=-windc*uw
 1 	*.25*(DUM(I,J+1)+DUM(I+1,J)+DUM(I-1,J)+DUM(I,J-1))+
 2  0.5*(d(i,j)+d(i-1,j))*(btim*FxB(i,j)+ftim*FxF(i,j))
 WVSURF(I,J)=-windc*vw
 1 	*.25*(DVM(I,J+1)+DVM(I+1,J)+DVM(I-1,J)+DVM(I,J-1))+
 2  0.5*(d(i,j)+d(i,j-1))*(btim*FyB(i,j)+ftim*FyF(i,j))
 end do
 end do
 */
    void makeWsurf(float ro_ratio)
    {
        float ftim = fmodf((float)m_fVars->timeh6, 1.0f);
        float btim = 1.0f - ftim;
        float uw, vw, speed, windc;
        
        int ji, jp1i, jip1, jim1, jm1i;
        
        float * g_fbu = &m_fFloats->fbu[0][0];
        float * g_fbv = &m_fFloats->fbv[0][0];
        float * g_ffu = &m_fFloats->ffu[0][0];
        float * g_ffv = &m_fFloats->ffv[0][0];

        float * g_wusurf = &m_fArrays->wusurf[0][0];
        float * g_wvsurf = &m_fArrays->wvsurf[0][0];

        float * g_dum = &m_fArrays->dum[0][0];
        float * g_dvm = &m_fArrays->dvm[0][0];

        float * g_d = &m_fArrays->d[0][0];

        
        for (int j=1; j<(m_height-1); j++ )
        {
            for (int i=1; i<(m_width-1); i++ )
            {
                ji = j * m_width + i;
                jp1i = ji + m_width;
                jip1 = ji + 1;
                jim1 = ji - 1;
                jm1i = ji - m_width;
                
                uw = btim * (g_fbu[ji]) + ftim * (g_ffu[ji]);
                vw = btim * (g_fbv[ji]) + ftim * (g_ffv[ji]);
                
                speed = sqrtf(uw*uw + vw*vw);
                windc = 0.001f * (0.8f + speed * 0.065f) * ro_ratio * speed;
                
                g_wusurf[ji] = -windc * uw *
                    0.25f * (g_dum[jp1i]+g_dum[jip1]+g_dum[jim1]+g_dum[jm1i])
                + 0.5f * (g_d[ji] + g_d[jim1]) * (btim * m_fFloats->fxb[j][i] + ftim * m_fFloats->fxf[j][i]);
                
                g_wvsurf[ji] = -windc * vw *
                0.25f * (g_dvm[jp1i]+g_dvm[jip1]+g_dvm[jim1]+g_dvm[jm1i])
                + 0.5f * (g_d[ji] + g_d[jm1i]) * (btim * m_fFloats->fyb[j][i] + ftim * m_fFloats->fyf[j][i]);
            }
        }
    }
    
    
    int m_width;
    int m_height;
    
    fortran_common_vars * m_fVars;
    fortran_common_arrays * m_fArrays;
    fortran_ffloats * m_fFloats;
};


double last_measured_time = 0.0;

KaspyCycler * cycler = 0;


void _i_cycler_time()
{
    double new_time;
    
#ifdef _WIN64
    new_time = GetTickCount()/1000.0;
#else
    struct timeval time;
    gettimeofday(&time, NULL); 

    new_time =  time.tv_sec*1.0 + time.tv_usec/1000000.0;
#endif
    
    ersetb();
    
    printf("time now is %10.10f\n", new_time);
    
    if (last_measured_time > 0.0)
    {
        printf("time elapsed after last measurement is %10.10f\n", new_time - last_measured_time);
    }
    
    last_measured_time = new_time;
}


// Intel Fortran WIN naming
extern "C" void CYCLER_TIME()
{
    _i_cycler_time();
}

// GFortran Unix naming
extern "C" void cycler_time_()
{
    _i_cycler_time();
}



int _i_cycler_init(float * vars_marker, double * arrays_marker, double * ffloats_marker)
{
    if ((vars_marker == 0) || (arrays_marker == 0)|| (ffloats_marker == 0))
    {
        return -1;
    }
    
    
    _i_cycler_time();
    
    cycler = new KaspyCycler((fortran_common_vars *)vars_marker,
                             (fortran_common_arrays *)arrays_marker, (fortran_ffloats *)ffloats_marker);
     
    
    return 0;
}


// Intel Fortran WIN naming
extern "C" int CYCLER_INIT(float * vars_marker, double * arrays_marker, double * ffloats_marker)
{
    return _i_cycler_init(vars_marker, arrays_marker, ffloats_marker);
}

// GFortran Unix naming
extern "C" int cycler_init_(float * vars_marker, double * arrays_marker, double * ffloats_marker)
{
    return _i_cycler_init(vars_marker, arrays_marker, ffloats_marker);
}


// GFortran Unix naming
extern "C" void cycler_wsurf_(int * icycler, float * ro_ratio)
{
    if (cycler)
    {
        cycler->makeWsurf(*ro_ratio);
    }
}

// Intel Fortran WIN naming
extern "C" void CYCLER_WSURF(int * icycler, float * ro_ratio)
{
    if (cycler)
    {
        cycler->makeWsurf(*ro_ratio);
    }
}


void _i_cycler_destroy(int * icycler)
{
    _i_cycler_time();
    
    if (cycler)
    {
        delete cycler;
    }
}

// GFortran Unix naming
extern "C" void cycler_destroy_(int * icycler)
{
    _i_cycler_destroy(icycler);
}

// Intel Fortran WIN naming
extern "C" void CYCLER_DESTROY(int * icycler)
{
    _i_cycler_destroy(icycler);
}



//
//  KaspyCycler.cpp
//  kaspy_cuda
//
//  Created by Andrei Koulikov on 24.05.16.
//
//

#include "KaspyCycler.h"







void KaspyCycler::findElves()
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
void KaspyCycler::makeWsurf(float ro_ratio)
{
    float ftim = fmodf((float)m_fVars->timeh6, 1.0f);
    float btim = 1.0f - ftim;
    float uw, vw, speed, windc;
    
    int ji, jp1i, jip1, jim1, jm1i;
    
    float * g_fbu = &m_fFloats->fbu[0][0];
    float * g_fbv = &m_fFloats->fbv[0][0];
    float * g_ffu = &m_fFloats->ffu[0][0];
    float * g_ffv = &m_fFloats->ffv[0][0];

    float * g_fxb = &m_fFloats->fxb[0][0];
    float * g_fxf = &m_fFloats->fxf[0][0];
    float * g_fyb = &m_fFloats->fyb[0][0];
    float * g_fyf = &m_fFloats->fyf[0][0];
   
    
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
            + 0.5f * (g_d[ji] + g_d[jim1]) * (btim * g_fxb[ji] + ftim * g_fxf[ji]);
            
            g_wvsurf[ji] = -windc * vw *
            0.25f * (g_dvm[jp1i]+g_dvm[jip1]+g_dvm[jim1]+g_dvm[jm1i])
            + 0.5f * (g_d[ji] + g_d[jm1i]) * (btim * g_fyb[ji] + ftim * g_fyf[ji]);
        }
    }
}

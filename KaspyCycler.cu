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



void KaspyCycler::loadData()
{
    g_fbu = &m_fFloats->fbu[0][0];
    g_fbv = &m_fFloats->fbv[0][0];
    g_ffu = &m_fFloats->ffu[0][0];
    g_ffv = &m_fFloats->ffv[0][0];
    
    g_fxb = &m_fFloats->fxb[0][0];
    g_fxf = &m_fFloats->fxf[0][0];
    g_fyb = &m_fFloats->fyb[0][0];
    g_fyf = &m_fFloats->fyf[0][0];

    
    g_fb = &m_fFloats->fb[0][0];
    g_ff = &m_fFloats->ff[0][0];
   
    
    g_wusurf = &m_fArrays->wusurf[0][0];
    g_wvsurf = &m_fArrays->wvsurf[0][0];
    
    g_dum = &m_fArrays->dum[0][0];
    g_dvm = &m_fArrays->dvm[0][0];
    
    g_d = &m_fArrays->d[0][0];
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
    m_fVars->timeh6 = (m_fVars->timeh / m_fVars->dht) + 1.0f;

    float timeh6 = m_fVars->timeh6;
    
    int pressSize = m_fWindData->kx * m_fWindData->ky;
    int windUSize = m_fWindData->kxu * m_fWindData->kyu;
    int windVSize = m_fWindData->kxv * m_fWindData->kyv;
    
    itime6 = (int)timeh6;

    ftim = (timeh6 - itime6);
    btim = 1.0f - ftim;
    
    if (itime6 > itime6_old)
    {
        itime6_old = itime6;
        
        memcpy(g_fxb, g_fxf, F_DATA_SIZE * sizeof(float));
        memcpy(g_fyb, g_fyf, F_DATA_SIZE * sizeof(float));
        memcpy(g_fb, g_ff, F_DATA_SIZE * sizeof(float));
        memcpy(g_fbu, g_ffu, F_DATA_SIZE * sizeof(float));
        memcpy(g_fbv, g_ffv, F_DATA_SIZE * sizeof(float));

        setbuf(stdout,NULL);
        
        //printf("press size is %d\n", pressSize );
        
        //printf("press 000 is %f press0 00 is %f\n", m_press[0], m_press0[0]);

        //printf("copy pressure from %#018llx to %#018llx\n", m_press, m_press0);

        memcpy(m_press0, m_press + (itime6 - 1) * pressSize, pressSize * sizeof(float));
        memcpy(m_uwd0, m_uwd + (itime6 - 1) * windUSize, windUSize * sizeof(float));
        memcpy(m_vwd0, m_vwd + (itime6 - 1) * windVSize, windVSize * sizeof(float));
        
    }
    
        /*press0(:,:)=press(:,:,itime6)
        call getnewpressureVAR(kx,ky,XKI,XKA,YKI,YKA,PRESS0,
                               1 FF,fxf,fyf)
        uwd0(:,:)=uwd(:,:,itime6)
        call getnewwindVAR(kxu,kyu,XKUI,XKUA,YKUI,YKUA,uwd0,ffu)
        vwd0(:,:)=vwd(:,:,itime6)
        call getnewwindVAR(kxv,kyv,XKVI,XKVA,YKVI,YKVA,vwd0,ffv)*/


            
    float uw, vw, speed, windc;
    int ji, jp1i, jip1, jim1, jm1i;

    
    
    
    ftim = fmodf((float)m_fVars->timeh6, 1.0f);
    btim = 1.0f - ftim;
    
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

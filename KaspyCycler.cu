//
//  KaspyCycler.cpp
//  kaspy_cuda
//
//  Created by Andrei Koulikov on 24.05.16.
//
//

#include "KaspyCycler.h"

// getnewpressureVAR(kx,ky,XKI,XKA,YKI,YKA,PRESS0,FF,fxf,fyf)

#ifdef _WIN64
extern "C"  void GETNEWPRESSUREVAR(int * kx, int * ky, float * xki, float * xka, float * yki, float * yka,
								   float * press0, float * ff, float * fxf, float * fyf);

extern "C"  void GETNEWWINDVAR(int * kxu, int * kyu, float * xkui, float * xkua,
							   float * ykui, float * ykua, float * uwd0, float * ffu);

#else

extern "C"  void getnewpressurevar_(int * kx, int * ky, float * xki, float * xka, float * yki, float * yka,
                              float * press0, float * ff, float * fxf, float * fyf);

extern "C"  void getnewwindvar_(int * kxu, int * kyu, float * xkui, float * xkua,
                               float * ykui, float * ykua, float * uwd0, float * ffu);

#endif





// call getnewwindVAR(kxu,kyu,XKUI,XKUA,YKUI,YKUA,uwd0,ffu)






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
float * g_dx;
float * g_dy;

float * g_fluxua;
float * g_fluxva;

float * g_ua;
float * g_va;

float * g_el;
float * g_elf;
float * g_elb;





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

void KaspyCycler::sendDataToGPU()
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
    g_dx = &m_fArrays->dx[0];
    g_dy = &m_fArrays->dy[0];


    g_fluxua = &m_fArrays->fluxua[0][0];
    g_fluxva = &m_fArrays->fluxva[0][0];
    
    g_ua = &m_fArrays->ua[0][0];
    g_va = &m_fArrays->va[0][0];
    
    g_el = &m_fArrays->el[0][0];
    g_elf = &m_fArrays->elf[0][0];
    g_elb = &m_fArrays->elb[0][0];
}

void KaspyCycler::getDataToCPU()
{
    
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

        //setbuf(stdout,NULL);
        
        //printf("press size is %d\n", pressSize );
        
        //printf("press 000 is %f press0 00 is %f\n", m_press[0], m_press0[0]);

        //printf("copy pressure from %#018llx to %#018llx\n", m_press, m_press0);

        memcpy(m_press0, m_press + (itime6 - 1) * pressSize, pressSize * sizeof(float));
		
#ifdef _WIN64
		GETNEWPRESSUREVAR(&m_fWindData->kx, &m_fWindData->ky, &m_fWindData->xki, &m_fWindData->xka,
						  &m_fWindData->yki, &m_fWindData->yka, m_press0, g_ff, g_fxf, g_fyf);
#else
		getnewpressurevar_(&m_fWindData->kx, &m_fWindData->ky, &m_fWindData->xki, &m_fWindData->xka,
						   &m_fWindData->yki, &m_fWindData->yka, m_press0, g_ff, g_fxf, g_fyf);
#endif
		
		
        memcpy(m_uwd0, m_uwd + (itime6 - 1) * windUSize, windUSize * sizeof(float));
		
		
#ifdef _WIN64
		GETNEWWINDVAR(&m_fWindData->kxu, &m_fWindData->kyu, &m_fWindData->xkui, &m_fWindData->xkua,
					  &m_fWindData->ykui, &m_fWindData->ykua, m_uwd0, g_ffu);
#else
		getnewwindvar_(&m_fWindData->kxu, &m_fWindData->kyu, &m_fWindData->xkui, &m_fWindData->xkua,
		              &m_fWindData->ykui, &m_fWindData->ykua, m_uwd0, g_ffu);
#endif

		
        memcpy(m_vwd0, m_vwd + (itime6 - 1) * windVSize, windVSize * sizeof(float));
		
#ifdef _WIN64
		GETNEWWINDVAR(&m_fWindData->kxv, &m_fWindData->kyv, &m_fWindData->xkvi, &m_fWindData->xkva,
					  &m_fWindData->ykvi, &m_fWindData->ykva, m_vwd0, g_ffv);
#else
		getnewwindvar_(&m_fWindData->kxv, &m_fWindData->kyv, &m_fWindData->xkvi, &m_fWindData->xkva,
		              &m_fWindData->ykvi, &m_fWindData->ykva, m_vwd0, g_ffv);		
#endif
        
		
		
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
    
    for (int j=1; j<m_height; j++ )
    {
        for (int i=1; i<m_width; i++ )
        {
            if ((j<(m_height-1)) && i<(m_width-1))
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
            
            //DO 405 J=2,JM
            //DO 405 I=2,IM
            //FLUXUA(I,J)=.25E0*(D(I,J)+D(I-1,J))*(DY(j)+DY(j))*UA(I,J)
            //405  FLUXVA(I,J)=.25E0*(D(I,J)+D(I,J-1))*(DX(j)+DX(j-1))*VA(I,J)

            g_fluxua[ji] = 0.25f * (g_d[ji] + g_d[jim1]) * (g_dy[j] + g_dy[j] /*???*/) * g_ua[ji];
            g_fluxva[ji] = 0.25f * (g_d[ji] + g_d[jm1i]) * (g_dx[j] + g_dx[j-1] ) * g_va[ji];
            
        }
    }
    
    
    /// HERE SHOULD START A NEW CUDA CALL TO KEEP fluxua fluxva synced
   
    /*DO 410 J=2,JMM1
    DO 410 I=2,IMM1
    410 ELF(I,J)=ELB(I,J)
    1    -DTE2*(FLUXUA(I+1,J)-FLUXUA(I,J)+FLUXVA(I,J+1)-FLUXVA(I,J))
    2                    / ART(J) */
    
    float dte2 = m_fVars->dte * 2.0f;
    
    for (int j=1; j<(m_height-1); j++ )
    {
        float artj = m_fArrays->art[j];
        
        for (int i=1; i<(m_width-1); i++ )
        {
            ji = j * m_width + i;
            jp1i = ji + m_width;
            jip1 = ji + 1;
            
            g_elf[ji] = g_elb[ji] - dte2 *
                (g_fluxua[jip1] - g_fluxua[ji] + g_fluxva[jp1i] - g_fluxva[ji]) /  artj;
            
        }
    }
 }








/*
 C     SURROUNDING
 DO J=2,KY+1
 DO I=2,KX+1
 PKK(I,J)=PK(I-1,J-1)
 END DO
 END DO
 DO J=2,KY+1
 PKK(1,J)=2*PKK(2,J)-PKK(3,J)
 PKK(KX+2,J)=2*PKK(KX+1,J)-PKK(KX,J)
 END DO
 DO I=1,KX+2
 PKK(I,1)=2*PKK(I,2)-PKK(I,3)
 PKK(I,KY+2)=2*PKK(I,KY+1)-PKK(I,KY)
 END DO
 CALL GETBICUBIC(KX+2,KY+2,50,PKK,C)
 
 
 
 do j=1,Ny
 y=ymi+(j-1)*dy
 j0=(y-yki)/dky+1
 if (j0<1) j0=1
 if (j0>ky-1) j0=ky-1
 u=(y-(yki+(j0-1)*dky))/dky
 
 do i=1,Nx
 x=xmi+(i-1)*dx
 i0=(x-xki)/dkx+1
 if (i0<1) I0=1
 if (i0>kx-1) i0=kx-1
 t=(x-(xki+(i0-1)*dkx))/dkx
 ay=0.
 a2=0.
 a1=0.
 DO K=4,1,-1
 ay=t*ay+((c(K,4,i0,j0)*u+c(k,3,i0,j0))*u+c(K,2,i0,j0))*u+
 1		  c(K,1,i0,j0)
 a2=t*a2+(3.*c(K,4,i0,j0)*u+2.*c(K,3,i0,j0))*u+c(K,2,i0,j0)
 a1=u*a1+(3.*c(4,K,i0,j0)*t+2.*c(3,K,i0,j0))*t+c(2,K,i0,j0)
 END DO
 a1=a1/dkx/c2/cos(c1*y)
 a2=a2/dky/c2
 
 p(i,j)=ay
 px(i,j)=a1
 py(i,j)=a2
 
 end do
 END DO
 
 CALL GETPRESScube(KX,KY,KX,PRESS0,XKI,XKA,YKI,YKA,
 1                     IM,JM,IM,FF,FXF,FYF,XMI,XMA,YMI,YMA)
 subroutine getpressCUBE(kx,ky,kd,pk,xki,xka,yki,yka,
 1                     nx,ny,nd,P,px,py,xmi,xma,ymi,yma)
 */


void KaspyCycler::getNewPressure()
{
	int kx = m_fWindData->kx;
	int ky = m_fWindData->ky;
	float kd = kx;
	float * pk = m_press0;
	float xki = m_fWindData->xki;
	float xka = m_fWindData->xka;
	float yki = m_fWindData->yki;
	float yka = m_fWindData->yka;
	int nx = F_DATA_WIDTH;
	int ny = F_DATA_HEIGHT;
	int nd = F_DATA_WIDTH;
	
	float * p = g_ff;
	float * px = g_fxf;
	float * py = g_fyf;
	
	float xmi = m_fVars->xmi;
	float xma = m_fVars->xma;
	float ymi = m_fVars->xmi;
	float yma = m_fVars->xma;
	
	float pkk[50][50];
	float c[50][50][4][4];
	
	float c1=3.1415926/180.0;
	float c2=111111.0f;
	
	float dky=(yka-yki)/(ky-1.0f);
	float  dkx=(xka-xki)/(kx-1.0f);
 
 	float dy=(yma-ymi)/(ny-1.0f);
 	float dx=(xma-xmi)/(nx-1.0f);
	
/*DO J=2,KY+1
 DO I=2,KX+1
 PKK(I,J)=PK(I-1,J-1)
 END DO
 END DO
 DO J=2,KY+1
 PKK(1,J)=2*PKK(2,J)-PKK(3,J)
 PKK(KX+2,J)=2*PKK(KX+1,J)-PKK(KX,J)
 END DO
 DO I=1,KX+2
 PKK(I,1)=2*PKK(I,2)-PKK(I,3)
 PKK(I,KY+2)=2*PKK(I,KY+1)-PKK(I,KY)
 END DO*/
	
	for (int j=1; j<ky; j++ )
	{
		for (int i=1; i<kx; i++ )
		{
			//int ji = j * kx + i;
			//int jm1i = ji - kx;
			//int jim1 = ji -  1;
			//int jm1im1 = jim1 - 1;
			
			pkk[j][i] = pk[j * (kx - 1) + i - 1];
		}
	}
	
	for (int j=1; j<ky; j++ )
	{
		pkk[j][0] = 2*pkk[j][1] - pkk[j][2];
		pkk[j][kx+1] = 2*pkk[j][kx] - pkk[j][kx-1];
	}
	
	for (int i=1; i<(kx+1); i++ )
	{
		pkk[0][i] = 2*pkk[1][i] - pkk[2][i];
		pkk[ky+1][i] = 2*pkk[ky][i] - pkk[ky-1][i];
	}
	
}







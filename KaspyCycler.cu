//
//  KaspyCycler.cpp
//  kaspy_cuda
//
//  Created by Andrei Koulikov on 24.05.16.
//
//

#include "KaspyCycler.h"

void getbicubic(int nx, int ny, int nd, float * z, float * c);
void bcucof(float * y,float * y1,float * y2, float * y12,float d1,float d2,float * cc);




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

float * g_uab;
float * g_vab;

float * g_uaf;
float * g_vaf;


float * g_el;
float * g_elf;
float * g_elb;

float * g_fsm;

float * g_tps;


float * g_advua;
float * g_advva;

float * g_aru;
float * g_arv;

float * g_wubot;
float * g_wvbot;
float * g_cbc;

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
	
	g_advua = &m_fArrays->advua[0][0];
	g_advva = &m_fArrays->advua[0][0];
	
    g_ua = &m_fArrays->ua[0][0];
    g_va = &m_fArrays->va[0][0];

	g_uab = &m_fArrays->uab[0][0];
	g_vab = &m_fArrays->vab[0][0];

	g_uaf = &m_fArrays->uaf[0][0];
	g_vaf = &m_fArrays->vaf[0][0];
	
	
    g_el = &m_fArrays->el[0][0];
    g_elf = &m_fArrays->elf[0][0];
    g_elb = &m_fArrays->elb[0][0];
	
	g_fsm = &m_fArrays->fsm[0][0];
	
	g_tps = &m_fArrays->tps[0][0];
	
	g_aru = &m_fArrays->aru[0];
	g_arv = &m_fArrays->arv[0];

	
	g_wubot = &m_fArrays->wubot[0][0];
	g_wvbot = &m_fArrays->wvbot[0][0];

	g_cbc = &m_fArrays->cbc[0][0];
}

void KaspyCycler::getDataToCPU()
{
    
}




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



        memcpy(m_press0, m_press + (itime6 - 1) * pressSize, pressSize * sizeof(float));
		
		getWindPressure('p');

		
        memcpy(m_uwd0, m_uwd + (itime6 - 1) * windUSize, windUSize * sizeof(float));
		
		
		getWindPressure('u');

        memcpy(m_vwd0, m_vwd + (itime6 - 1) * windVSize, windVSize * sizeof(float));
		
		getWindPressure('v');

    }
	
	
    float uw, vw, speed, windc;
    int ji, jp1i, jip1, jim1, jm1i, jp1ip1, jm1im1, jp1im1, jm1ip1;


    
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
            


            g_fluxua[ji] = 0.25f * (g_d[ji] + g_d[jim1]) * (g_dy[j] + g_dy[j] ) * g_ua[ji];
            g_fluxva[ji] = 0.25f * (g_d[ji] + g_d[jm1i]) * (g_dx[j] + g_dx[j-1] ) * g_va[ji];
            
        }
    }
    
    
    /// HERE SHOULD START A NEW CUDA CALL TO KEEP fluxua fluxva synced
	
    
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
 


	/// BCOND 1
	float tide_l = m_fVars->tide_l;
	
	for (int j=1; j<m_height; j++ )
	{
		g_elf[j * m_width + 1] = tide_l;
		g_elf[j * m_width + m_width - 2] = tide_l;
		
		g_elf[j * m_width] = tide_l;
		g_elf[j * m_width + m_width - 1] = tide_l;
	}
	
	for (int i=1; i<m_width; i++ )
	{
		g_elf[i] =  g_elf[i + m_width];
		
		g_elf[i + m_width * (m_height - 1)  ] =  g_elf[i + m_width * (m_height - 2)];
	}
	
	for (int j=1; j<m_height; j++ )
	{
		for (int i=1; i<m_width; i++ )
		{
			ji = j * m_width + i;
			
			g_elf[ji] *= g_fsm[ji];
		}
	}

	
	if (m_fVars->iint % 10 == 0)
	{
		//ADVAVE()
		
		//       ADVUA=0
		//		FLUXUA=0
		
		//memset(g_advua, 0, F_DATA_SIZE * sizeof(float));
		//memset(g_fluxua, 0, F_DATA_SIZE * sizeof(float));
		
		
		float aam2d = m_fArrays->aam2d;
		
		for (int j=1; j<m_height; j++ )
		{
			for (int i=1; i<(m_width-1); i++ )
			{
				ji = j * m_width + i;
				jip1 = ji + 1;
				jim1 = ji - 1;
				
				/*g_fluxua[ji] = g_dy[j] * (.125f * ((g_d[ji + 1]+g_d[ji])*g_ua[ji + 1]
						+(g_d[ji]+g_d[ji - 1])*g_ua[ji])
										  *(g_ua[ji + 1]+g_ua[ji])
										  - g_d[ji]*2.0f*aam2d*(g_uab[ji + 1]-g_uab[ji])/g_dx[j]);*/
				g_fluxua[ji]=g_dy[j]*(.125e0*((g_d[jip1]+g_d[ji])*g_ua[jip1]
											  +(g_d[ji]+g_d[jim1])*g_ua[ji])
									  *(g_ua[jip1]+g_ua[ji])
									  -g_d[ji]*2.e0*aam2d*(g_uab[jip1]-g_uab[ji])/g_dx[j]);
				
				
			}
		}
		
		
		for (int j=1; j<m_height; j++ )
		{
			for (int i=1; i<m_width; i++ )
			{
				ji = j * m_width + i;
				jp1i = ji + m_width;
				jip1 = ji + 1;
				jim1 = ji - 1;
				jm1i = ji - m_width;
				jm1im1 = jm1i  - 1;
				
				/*g_tps[ji] =(g_d[ji]+g_d[jim1]+g_d[jm1i]+g_d[jm1im1]) *aam2d
				*((g_uab[ji]-g_uab[jm1i]) /(4.0f*g_dy[j])+(g_vab[ji]-g_vab[jim1]) /(4.0f*g_dx[j]) );
				
				g_fluxva[ji]=(.125f*((g_d[ji]+g_d[jm1i])*g_va[ji]
									 +(g_d[jim1]+g_d[jm1im1])*g_va[jim1])
							  *(g_ua[ji]+g_va[jm1i]) - g_tps[ji])*g_dx[j];*/
				
				g_tps[ji]=(g_d[ji]+g_d[jim1]+g_d[jm1i]+g_d[jm1im1])
				*aam2d
				*((g_uab[ji]-g_uab[jm1i])
				  /(4*g_dy[j])
				  +(g_vab[ji]-g_vab[jim1])
				  /(4*g_dx[j]) );
				
				g_fluxva[ji]=(.125e0*((g_d[ji]+g_d[jm1i])*g_va[ji]
									  +(g_d[jim1]+g_d[jm1im1])*g_va[jim1])
							  *(g_ua[ji]+g_ua[jm1i])
							  -g_tps[ji])*g_dx[j];
				
			}
		}

		
		for (int j=1; j<(m_height-1); j++ )
		{
			for (int i=1; i<(m_width-1); i++ )
			{
				ji = j * m_width + i;
				jim1 = ji - 1;
				jp1i = ji + m_width;
				
				g_advua[ji]=(g_fluxua[ji]-g_fluxua[jim1]
						   +g_fluxva[jp1i]-g_fluxva[ji])/g_aru[j];
			}
			
		}
		
		//memset(g_advva, 0, F_DATA_SIZE * sizeof(float));
		//memset(g_fluxva, 0, F_DATA_SIZE * sizeof(float));
		
		
		for (int j=1; j<(m_height-1); j++ )
		{
			for (int i=1; i<m_width; i++ )
			{
				ji = j * m_width + i;
				jp1i = ji + m_width;
				jip1 = ji + 1;
				jim1 = ji - 1;
				jm1i = ji - m_width;
				jm1im1 = jm1i  - 1;
				
				
			 	g_fluxva[ji]=g_dx[j]*(.125e0*((g_d[jp1i]+g_d[ji])
									       *g_va[jp1i]+(g_d[ji]+g_d[jm1i])*g_va[ji])
									      *(g_va[jp1i]+g_va[ji])
								         -g_d[ji]*2.e0*aam2d*(g_vab[jp1i]-g_vab[ji])/g_dy[j]);
				
			}
		}
		
		
		for (int j=1; j<m_height; j++ )
		{
			for (int i=1; i<m_width; i++ )
			{
				ji = j * m_width + i;
				jp1i = ji + m_width;
				jip1 = ji + 1;
				jim1 = ji - 1;
				jm1i = ji - m_width;
				jm1im1 = jm1i  - 1;
				
				
				g_fluxua[ji]=(.125e0*((g_d[ji]+g_d[jim1])*g_ua[ji]
									         +(g_d[jm1i]+g_d[jm1im1])*g_ua[jm1i])*
							                        (g_va[jim1]+g_va[ji])
							  -g_tps[ji])*g_dy[j];
			}
		}
		
		for (int j=1; j<(m_height-1); j++ )
		{
			for (int i=1; i<(m_width-1); i++ )
			{
				ji = j * m_width + i;
				jp1i = ji + m_width;
				jip1 = ji + 1;
				jim1 = ji - 1;
				jm1i = ji - m_width;
				jm1im1 = jm1i  - 1;
				
				g_advva[ji]=(g_fluxua[jip1]-g_fluxua[ji]
							         +g_fluxva[ji]-g_fluxva[jm1i])/g_arv[j];
			}
		}
	
		
		for (int j=1; j<(m_height-1); j++ )
		{
			for (int i=1; i<(m_width-1); i++ )
			{
				ji = j * m_width + i;
				jp1i = ji + m_width;
				jip1 = ji + 1;
				jim1 = ji - 1;
				jm1i = ji - m_width;
				jm1im1 = jm1i  - 1;
				
				jp1im1 = jp1i - 1;
				jm1ip1 = jm1i + 1;

				g_wubot[ji]=-0.5e0*(g_cbc[ji]+g_cbc[jim1])
				     *sqrtf(g_uab[ji]*g_uab[ji]+powf(.25e0*(g_vab[ji]
											  +g_vab[jp1i]+g_vab[jim1]+g_vab[jp1im1]), 2) )*g_uab[ji];
				
				g_wvbot[ji]=-0.5e0*(g_cbc[ji]+g_cbc[jm1i])
				    *sqrtf((.25e0*(g_uab[ji]+g_uab[jip1]
								  +g_uab[jm1i]+g_uab[jm1ip1]))**2+g_vab[ji]*g_vab[ji])*g_vab[ji];
				
			}
		}
		
		
		
	}
	
	
	
	
	
}




void KaspyCycler::getWindPressure(char uv)
{
	int kx, ky, kd, nx, ny, nd;
	float * p;
	float * px;
	float * py;
	float * pk;
	float xki, xka, yki, yka, xmi, xma, ymi, yma;
	float pkkd[50][50];
	float cd[50][50][4][4];
	
	float * pkk = &pkkd[0][0];
	float * c = &cd[0][0][0][0];
	
	if (uv == 'u')
	{
		kx = m_fWindData->kxu;
		ky = m_fWindData->kyu;
		pk = m_uwd0;
		
		xki = m_fWindData->xkui;
		xka = m_fWindData->xkua;
		yki = m_fWindData->ykui;
		yka = m_fWindData->ykua;
		
		p = g_ffu;
	}
	else if (uv == 'v')
	{
		kx = m_fWindData->kxv;
		ky = m_fWindData->kyv;
		pk = m_vwd0;
		
		xki = m_fWindData->xkvi;
		xka = m_fWindData->xkva;
		yki = m_fWindData->ykvi;
		yka = m_fWindData->ykva;
		
		p = g_ffv;
	}
	else if (uv == 'p')
	{
		kx = m_fWindData->kx;
		ky = m_fWindData->ky;
		//float kd = kx;
		pk = m_press0;
		xki = m_fWindData->xki;
		xka = m_fWindData->xka;
		yki = m_fWindData->yki;
		yka = m_fWindData->yka;

		
		p = g_ff;
		px = g_fxf;
		py = g_fyf;
	}
	else
	{
		// don't know what to do
		return;
	}
	
	kd = kx;
	
	nx = F_DATA_WIDTH;
	ny = F_DATA_HEIGHT;
	nd = F_DATA_WIDTH;
	
	xmi = m_fVars->xmi;
	xma = m_fVars->xma;
	ymi = m_fVars->ymi;
	yma = m_fVars->yma;
	
	float c1=3.1415926/180.0;
	float c2=111111.0f;
	
	
	float dky=(yka-yki)/(ky-1.0f);
	float  dkx=(xka-xki)/(kx-1.0f);
 
	float dy=(yma-ymi)/(ny-1.0f);
	float dx=(xma-xmi)/(nx-1.0f);
	
	
	for (int j=1; j<=ky; j++ )
	{
		for (int i=1; i<=kx; i++ )
		{
			pkk[j * 50 + i] = pk[(j - 1) * kd + i - 1];
		}
	}

	
	for (int j=1; j<=ky; j++ )
	{
		pkk[j*50+0] = 2.0f*pkk[j*50+1] - pkk[j*50+2];
		pkk[j*50+kx+1] = 2.0f*pkk[j*50+kx] - pkk[j*50+kx-1];
	}
	
	
	for (int i=0; i<=(kx+1); i++ )
	{
		pkk[0*50+i] = 2.0f*pkk[1*50+i] - pkk[2*50+i];
		pkk[(ky+1)*50+i] = 2.0f*pkk[ky*50+i] - pkk[(ky-1)*50+i];
	}
	
	getbicubic(kx + 2,ky + 2, 50, pkk,c);
	
	for (int j=0; j<ny; j++ )
	{
		float y = ymi + j*dy;
		int j0 = (int)((y - yki)/dky);
		
		if (j0 < 0)
		{
			j0 = 0;
		}
		
		if (j0 > ky-2)
		{
			j0 = ky-2;
		}
		
		float u = (y - (yki + j0*dky))/dky;
		
		for (int i=0; i<nx; i++ )
		{
			float x = xmi + i * dx;
			int i0 = (int)((x - xki)/dkx);
			
			if (i0 < 0) i0 = 0;
			
			if (i0 > kx-2) i0 = kx-2;
			
			float t = ( x - (xki + i0*dkx) )/dkx;
			
			float ay = 0.0f;
			float a2 = 0.0f;
			float a1 = 0.0f;
			
			int ji = j * nx + i;
			
			for (int k=3; k>=0; k-- )
			{
				ay = t*ay+((c[j0 * 800 + i0 * 16 + 3 * 4 + k] * u + c[j0 * 800 + i0 * 16 + 2 * 4 + k])*u
						   + c[j0 * 800 + i0 * 16 + 1 * 4 + k])*u + c[j0 * 800 + i0 * 16 + 0 * 4 + k];
			}
			
			if (uv == 'p')
			{
				for (int k=3; k>=0; k-- )
				{
					a2 = t*a2 + (3.0f*c[j0 * 800 + i0 * 16 + 3 * 4 + k]*u
								 + 2.0f*c[j0 * 800 + i0 * 16 + 2 * 4 + k])*u+c[j0 * 800 + i0 * 16 + 1 * 4 + k];
					
					a1 = u*a1 + (3.0f*c[j0 * 800 + i0 * 16 + k * 4 + 3]*t +
								 2.0f*c[j0 * 800 + i0 * 16 + k * 4 + 2])*t+c[j0 * 800 + i0 * 16 + k * 4 + 1];
					
				}
				
				a1 = a1/dkx/c2/cosf(c1*y);
				a2 = a2/dky/c2;
				
				px[ji] = a1;
				py[ji] = a2;
			}
			
			p[ji] = ay;
			
		}
		
	}
	

}







void getbicubic(int nx, int ny, int nd, float * z, float * c)
{
	float d1 = 1.0f;
	float d2 = 1.0f;
	
	float y[4];
	float y1[4];
	float y2[4];
	float y12[4];
	float cc[4][4];
	
	
	for (int j=1; j<ny-2; j++ )
	{
		for (int i=1; i<nx-2; i++ )
		{
			/*
			 Y(1)=Z(I,J)
			 Y(2)=Z(I+1,J)
			 Y(3)=Z(I+1,J+1)
			 Y(4)=Z(I,J+1)
			 */
			y[0] = z[j * nd + i];
			y[1] = z[j * nd + i + 1];
			y[2] = z[(j+1) * nd + i + 1];
			y[3] = z[(j+1) * nd + i];
			
			/*
			 Y1(1)=0.5*(Z(I+1,J)-Z(I-1,J))
			 Y1(4)=0.5*(Z(I+1,J+1)-Z(I-1,J+1))
			 Y1(2)=0.5*(Z(I+2,J)  -Z(I,J))
			 Y1(3)=0.5*(Z(I+2,J+1)-Z(I,J+1))
			 */
			y1[0] = 0.5f * (z[j * nd + i + 1] - z[j * nd + i - 1]);
			y1[3] = 0.5f * (z[(j+1) * nd + i + 1] - z[(j+1) * nd + i - 1]);
			y1[1] = 0.5f * (z[j * nd + i + 2] - z[j * nd + i]);
			y1[2] = 0.5f * (z[(j+1) * nd + i + 2] - z[(j+1) * nd + i]);

			
			/*
			 Y2(1)=0.5*(Z(I,J+1)  -Z(I,J-1))
			 Y2(2)=0.5*(Z(I+1,J+1)-Z(I+1,J-1))
			 Y2(3)=0.5*(Z(I+1,J+2)-Z(I+1,J))
			 Y2(4)=0.5*(Z(I,J+2)-Z(I,J))
			 */
			y2[0] = 0.5f * (z[(j+1) * nd + i] - z[(j-1) * nd + i]);
			y2[1] = 0.5f * (z[(j+1) * nd + i + 1] - z[(j-1) * nd + i + 1]);
			y2[2] = 0.5f * (z[(j+2) * nd + i + 1] - z[(j) * nd + i + 1]);
			y2[3] = 0.5f * (z[(j+2) * nd + i] - z[j * nd + i]);
			
			
			/*
			 Y12(1)=0.25*(Z(I+1,J+1)-Z(I+1,J-1)-Z(I-1,J+1)+Z(I-1,J-1))
			 Y12(2)=0.25*(Z(I+2,J+1)-Z(I+2,J-1)-Z(I,J+1)+Z(I,J-1))
			 Y12(3)=0.25*(Z(I+2,J+2)-Z(I+2,J)-Z(I,J+2)+Z(I,J))
			 Y12(4)=0.25*(Z(I+1,J+2)-Z(I+1,J)-Z(I-1,J+2)+Z(I-1,J))
			 */
			y12[0] = 0.25f * (z[(j+1) * nd + i + 1] - z[(j-1) * nd + i + 1]
							  - z[(j+1) * nd + i - 1] + z[(j-1) * nd + i - 1]);
			y12[1] = 0.25f * (z[(j+1) * nd + i + 2] - z[(j-1) * nd + i + 2]
							  - z[(j+1) * nd + i] + z[(j-1) * nd + i]);
			y12[2] = 0.25f * (z[(j+2) * nd + i + 2] - z[(j) * nd + i + 2]
							  - z[(j+2) * nd + i] + z[j * nd + i]);
			y12[3] = 0.25f * (z[(j+2) * nd + i + 1] - z[(j) * nd + i + 1]
							  - z[(j+2) * nd + i -1] + z[(j) * nd + i -1]);
	
			
			bcucof(&y[0],&y1[0],&y2[0],&y12[0],d1,d2,&cc[0][0]);
			
			for (int k=0; k<4; k++ )
			{
				for (int l=0; l<4; l++ )
				{
					//printf("\nk is %d l is %d\n", k, l);
					c[(j-1)* 800 + (i-1) * 16 + l * 4 + k ] = cc[l][k];
				}
			}
			
			
		}
	 }
	
}




void bcucof(float * y,float * y1,float * y2, float * y12,float d1,float d2,float * cc)
{
	float xx;
	float cl[16];
	
	float x[16];
	
	float wt[] = {
		1,0,-3,2,0,0,0,0,-3,0,9,-6,2,0,-6,4,
		0,0,0,0,0,0,0,0,3,0,-9,6,-2,0,6,-4,
		0,0,0,0,0,0,0,0,0,0,9,-6,0,0,-6,4,
		0,0,3,-2,0,0,0,0,0,0,-9,6,0,0,6,-4,
		0,0,0,0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,
		0,0,0,0,0,0,0,0,-1,0,3,-2,1,0,-3,2,
		0,0,0,0,0,0,0,0,0,0,-3,2,0,0,3,-2,
		0,0,0,0,0,0,3,-2,0,0,-6,4,0,0,3,-2,
		0,1,-2,1,0,0,0,0,0,-3,6,-3,0,2,-4,2,
		0,0,0,0,0,0,0,0,0,3,-6,3,0,-2,4,-2,
		0,0,0,0,0,0,0,0,0,0,-3,3,0,0,2,-2,
		0,0,-1,1,0,0,0,0,0,0,3,-3,0,0,-2,2,
		0,0,0,0,0,1,-2,1,0,-2,4,-2,0,1,-2,1,
		0,0,0,0,0,0,0,0,0,-1,2,-1,0,1,-2,1,
		0,0,0,0,0,0,0,0,0,0,1,-1,0,0,-1,1,
		0,0,0,0,0,0,-1,1,0,0,2,-2,0,0,-1,1
 	};
	
	//float d1 = *pd1;
	//float d2 = *pd2;
	
	
	float d1d2 = d1 * d2;

	for (int i=0; i<4; i++ )
	{
		x[i] = y[i];
		x[i + 4] = y1[i] * d1;
		x[i + 8] = y2[i] * d2;
		x[i + 12] = y12[i] * d1d2;
	}
	
	for (int i=0; i<16; i++ )
	{
		xx = 0.0f;
		
		for (int k=0; k<16; k++ )
		{
			xx += wt[i + k*16] * x[k];
		}
		
		cl[i] = xx;
	}
	
	int l = 0;
	
	for (int i=0; i<4; i++ )
	{
		for (int j=0; j<4; j++ )
		{
			cc[j*4 + i] = cl[l++];
		}
	}
	
}







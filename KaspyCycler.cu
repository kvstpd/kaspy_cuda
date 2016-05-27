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
		
		getNewWind('p');
		
		
		
        memcpy(m_uwd0, m_uwd + (itime6 - 1) * windUSize, windUSize * sizeof(float));
		
		
		getNewWind('u');
		

		
        memcpy(m_vwd0, m_vwd + (itime6 - 1) * windVSize, windVSize * sizeof(float));
		
		getNewWind('v');
		

        
		
		
    }
	


            
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
            


            g_fluxua[ji] = 0.25f * (g_d[ji] + g_d[jim1]) * (g_dy[j] + g_dy[j] ) * g_ua[ji];
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




void KaspyCycler::getNewWind(char uv)
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




void KaspyCycler::getNewPressure()
{
	int kx = m_fWindData->kx;
	int ky = m_fWindData->ky;
	//float kd = kx;
	float * pk = m_press0;
	float xki = m_fWindData->xki;
	float xka = m_fWindData->xka;
	float yki = m_fWindData->yki;
	float yka = m_fWindData->yka;
	int nx = F_DATA_WIDTH;
	int ny = F_DATA_HEIGHT;
	//int nd = F_DATA_WIDTH;
	
	float * p = g_ff;
	float * px = g_fxf;
	float * py = g_fyf;
	
	float xmi = m_fVars->xmi;
	float xma = m_fVars->xma;
	float ymi = m_fVars->ymi;
	float yma = m_fVars->yma;
	
	float pkkd[50][50];
	float cd[50][50][4][4];
	
	float * pkk = &pkkd[0][0];
	float * c = &cd[0][0][0][0];
	
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
			pkk[j * 50 + i] = pk[(j - 1) * kx + i - 1];
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
			
			for (int k=3; k>=0; k-- )
			{
				ay = t*ay+((c[j0 * 800 + i0 * 16 + 3 * 4 + k] * u + c[j0 * 800 + i0 * 16 + 2 * 4 + k])*u
						   + c[j0 * 800 + i0 * 16 + 1 * 4 + k])*u + c[j0 * 800 + i0 * 16 + 0 * 4 + k];
				
				a2 = t*a2 + (3.0f*c[j0 * 800 + i0 * 16 + 3 * 4 + k]*u
							 + 2.0f*c[j0 * 800 + i0 * 16 + 2 * 4 + k])*u+c[j0 * 800 + i0 * 16 + 1 * 4 + k];
				
				a1 = u*a1 + (3.0f*c[j0 * 800 + i0 * 16 + k * 4 + 3]*t +
							 2.0f*c[j0 * 800 + i0 * 16 + k * 4 + 2])*t+c[j0 * 800 + i0 * 16 + k * 4 + 1];
				
			}
			
			a1 = a1/dkx/c2/cosf(c1*y);
			a2 = a2/dky/c2;
			
			int ji = j * nx + i;
			
			p[ji] = ay;
			px[ji] = a1;
			py[ji] = a2;
			
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







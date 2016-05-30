//
//  KaspyCycler.cpp
//  kaspy_cuda
//
//  Created by Andrei Koulikov on 24.05.16.
//
//

#include "KaspyCycler.h"


#include "cuda_runtime.h"
#include "device_launch_parameters.h"



void getbicubic(int nx, int ny, int nd, float * z, float * c);
void bcucof(float * y,float * y1,float * y2, float * y12,float d1,float d2,float * cc);


float grav = 9.806;


float * g_fbu = 0;
float * g_fbv = 0;
float * g_ffu = 0;
float * g_ffv = 0;

float * g_fxb = 0;
float * g_fxf = 0;
float * g_fyb = 0;
float * g_fyf = 0;

float * g_fb = 0;
float * g_ff = 0;

float * g_wusurf = 0;
float * g_wvsurf = 0;

float * g_dum = 0;
float * g_dvm = 0;

float * g_d = 0;
float * g_dx = 0;
float * g_dy = 0;

float * g_fluxua = 0;
float * g_fluxva = 0;

float * g_ua = 0;
float * g_va = 0;

float * g_uab = 0;
float * g_vab = 0;

float * g_uaf = 0;
float * g_vaf = 0;


float * g_el = 0;
float * g_elf = 0;
float * g_elb = 0;

float * g_fsm = 0;

float * g_tps = 0;


float * g_advua = 0;
float * g_advva = 0;

float * g_aru = 0;
float * g_arv = 0;

float * g_wubot = 0;
float * g_wvbot = 0;
float * g_cbc = 0;

float * g_cor = 0;

float * g_h = 0;

float * g_press0 = 0;
float * g_uwd0 = 0;
float * g_vwd0 = 0;


__constant__ __device__  float g_grav = 9.806;

__constant__ __device__  int  g_width;
__constant__ __device__  int  g_height;
__constant__ __device__  int  g_widthm1;
__constant__ __device__  int  g_heightm1;

__constant__ __device__ int g_ewidth;


/*, float * g_fbu, float * g_ffu, float * g_fbv, float * g_ffv, float * g_dum, float * g_dvm, float * g_d, float * g_wusurf, float * g_wvsurf, float * g_fluxua, float * g_fluxva, float * g_dx, float * g_dy, float * g_ua, float * g_va, float * g_fxf, float * g_fyf, float * g_fxb, float *  g_fyb*/

__global__ void surf_and_flux_1(float ftim, float ro_ratio)
{
	/*int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	int ji = j * g_width + i;
 	int jp1i = ji + g_width;
 	int jip1 = ji + 1;
 	int jim1 = ji - 1;
 	int jm1i = ji - g_width;*/
	
	float btim = 1.0f - ftim;
	
	
	/*if (i < g_widthm1 && j < g_heightm1)
	{
		float uw = btim * (g_fbu[ji]) + ftim * (g_ffu[ji]);
		float vw = btim * (g_fbv[ji]) + ftim * (g_ffv[ji]);
		
		float speed = sqrtf(uw*uw + vw*vw);
		float windc = 0.001f * (0.8f + speed * 0.065f) * ro_ratio * speed;
		
		g_wusurf[ji] = -windc * uw *
		0.25f * (g_dum[jp1i]+g_dum[jip1]+g_dum[jim1]+g_dum[jm1i])
		+ 0.5f * (g_d[ji] + g_d[jim1]) * (btim * g_fxb[ji] + ftim * g_fxf[ji]);
		
		g_wvsurf[ji] = -windc * vw *
		0.25f * (g_dvm[jp1i]+g_dvm[jip1]+g_dvm[jim1]+g_dvm[jm1i])
		+ 0.5f * (g_d[ji] + g_d[jm1i]) * (btim * g_fyb[ji] + ftim * g_fyf[ji]);
	}
	
	if (i < g_width && j < g_height)
	{
		g_fluxua[ji] = 0.25f * (g_d[ji] + g_d[jim1]) * (g_dy[j] + g_dy[j] ) * g_ua[ji];
		g_fluxva[ji] = 0.25f * (g_d[ji] + g_d[jm1i]) * (g_dx[j] + g_dx[j-1] ) * g_va[ji];
	}*/
}


/*
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
 }*/





void KaspyCycler::findElves()
{
	/// DO CUDA REDUCTION instead of copying back to host mem
	
	float * h_elf =  &m_fArrays->elf[0][0];
	
	cudaError_t err = cudaMemcpy(h_elf, g_elf,  m_height * m_width * sizeof(float), cudaMemcpyDeviceToHost);
	
	if (err != cudaSuccess)
	{
		fprintf(stderr, "Failed to update host array ELF  (error code %s)!\n", cudaGetErrorString(err));
	}
	

	
	float elf_min = h_elf[0];
    float elf_max = h_elf[0];
    
    for (int i=1; i<F_DATA_SIZE; i++)
    {
        if (h_elf[i] > elf_max)
        {
            elf_max = h_elf[i];
        }
        
        if (h_elf[i] < elf_min)
        {
            elf_min = h_elf[i];
        }
    }
	
	m_fVars->elfmin =  elf_min;
	m_fVars->elfmax =  elf_max;
}


void KaspyCycler::sendDataToGPU()
{
	//int ewidth = ((int)m_pitch) / sizeof(float);
	int wm1 = m_width - 1 ;
	int hm1 = m_height - 1 ;
	
	if ( (cudaMemcpyToSymbol(g_width, &m_width, sizeof(int))  == cudaSuccess)
		&& (cudaMemcpyToSymbol(g_height, &m_height, sizeof(int))  == cudaSuccess)
		&&(cudaMemcpyToSymbol(g_widthm1, &wm1, sizeof(int))  == cudaSuccess)
		&& (cudaMemcpyToSymbol(g_heightm1, &hm1, sizeof(int))  == cudaSuccess)
		//&& (cudaMemcpyToSymbol(g_ewidth, &ewidth,  sizeof(int))  == cudaSuccess)
		)
	{
		printf("GPU constant memory filled\n");
		
		//int test_ewidth = 0;
		
		//cudaMemcpyFromSymbol(&test_ewidth, g_ewidth, sizeof(int));
		
		
		//printf("pitched width is now %d\n", test_ewidth);
	}
	else
	{
		printf("GPU memory copy error (error code %s)!\n", cudaGetErrorString(cudaGetLastError()));
	}
	
	
	size_t s_data_size =  m_height * m_width *  sizeof(float);
	
	if ( (cudaMemcpy(g_fbu,&m_fFloats->fbu[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_fbv,&m_fFloats->fbv[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_ffu,&m_fFloats->ffu[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_ffv,&m_fFloats->ffv[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_fb,&m_fFloats->fb[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_ff,&m_fFloats->ff[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_fxb,&m_fFloats->fxb[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_fxf,&m_fFloats->fxf[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_fyb,&m_fFloats->fyb[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_fyf,&m_fFloats->fyf[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_wusurf,&m_fArrays->wusurf[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_wvsurf,&m_fArrays->wvsurf[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_dum,&m_fArrays->dum[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_dvm,&m_fArrays->dvm[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_d, &m_fArrays->d[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		
		&& (cudaMemcpy(g_fluxua,&m_fArrays->fluxua[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_fluxva,&m_fArrays->fluxva[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		
		&& (cudaMemcpy(g_ua,&m_fArrays->ua[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_va,&m_fArrays->va[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_uab,&m_fArrays->uab[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_vab,&m_fArrays->vab[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_uaf,&m_fArrays->uaf[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_vaf,&m_fArrays->vaf[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_el,&m_fArrays->el[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_elb,&m_fArrays->elb[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_elf,&m_fArrays->elf[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_fsm,&m_fArrays->fsm[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_tps,&m_fArrays->tps[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_advua,&m_fArrays->advua[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_advva,&m_fArrays->advva[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_wubot,&m_fArrays->wubot[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_wvbot,&m_fArrays->wvbot[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_cbc,&m_fArrays->cbc[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_h,&m_fArrays->h[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		
		
		&& (cudaMemcpy(g_cor, &m_fArrays->cor[0], m_height * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_aru, &m_fArrays->aru[0], m_height * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_arv, &m_fArrays->arv[0],  m_height * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_dx, &m_fArrays->dx[0], m_height * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_dy, &m_fArrays->dy[0], m_height * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess))
		
	{
		printf("GPU memory filled\n");
	}
	else
	{
		printf("GPU memory copy error!\n");

	}
	
	/*g_fbu = &m_fFloats->fbu[0][0];
    g_fbv = &m_fFloats->fbv[0][0];
    g_ffu = &m_fFloats->ffu[0][0];
    g_ffv = &m_fFloats->ffv[0][0];
    
    g_fxb = &m_fFloats->fxb[0][0];
    g_fxf = &m_fFloats->fxf[0][0];
    g_fyb = &m_fFloats->fyb[0][0];
    g_fyf = &m_fFloats->fyf[0][0];
    
    
    g_fb = &m_fFloats->fb[0][0];
    g_ff = &m_fFloats->ff[0][0];
	
	g_h = &m_fArrays->h[0][0];
    
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
	g_advva = &m_fArrays->advva[0][0];
	
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
	
	g_cor = &m_fArrays->cor[0];*/
}

void KaspyCycler::getDataToCPU()
{
    
}




void KaspyCycler::makeWsurf(float ro_ratio)
{
    m_fVars->timeh6 = (m_fVars->timeh / m_fVars->dht) + 1.0f;

    float timeh6 = (float)m_fVars->timeh6;
	
    int pressSize = m_fWindData->kx * m_fWindData->ky;
    int windUSize = m_fWindData->kxu * m_fWindData->kyu;
    int windVSize = m_fWindData->kxv * m_fWindData->kyv;
	
	size_t s_width =  m_width *  sizeof(float);
    
    itime6 = (int)timeh6;

    //ftim = (timeh6 - itime6);
    //btim = 1.0f - ftim;
    
    if (itime6 > itime6_old)
    {
        itime6_old = itime6;
        
        /*memcpy(g_fxb, g_fxf, F_DATA_SIZE * sizeof(float));
        memcpy(g_fyb, g_fyf, F_DATA_SIZE * sizeof(float));
        memcpy(g_fb, g_ff, F_DATA_SIZE * sizeof(float));
        memcpy(g_fbu, g_ffu, F_DATA_SIZE * sizeof(float));
        memcpy(g_fbv, g_ffv, F_DATA_SIZE * sizeof(float));*/
		
		if ( (cudaMemcpy(g_fxb,g_fxf, F_DATA_SIZE * sizeof(float), cudaMemcpyDeviceToDevice) == cudaSuccess)
			&& (cudaMemcpy(g_fyb,g_fyf, F_DATA_SIZE * sizeof(float), cudaMemcpyDeviceToDevice) == cudaSuccess)
			&& (cudaMemcpy(g_fb,g_ff, F_DATA_SIZE * sizeof(float), cudaMemcpyDeviceToDevice) == cudaSuccess)
			&& (cudaMemcpy(g_fbu,g_ffu, F_DATA_SIZE * sizeof(float), cudaMemcpyDeviceToDevice) == cudaSuccess)
			&& (cudaMemcpy(g_fbv,g_ffv, F_DATA_SIZE * sizeof(float), cudaMemcpyDeviceToDevice) == cudaSuccess)
			)
		{
			printf("ff arrays reset\n");
		}
		else
		{
			printf("GPU memory copy error!\n");
		}
			
			

		size_t s_p_width = m_fWindData->kx * sizeof(float);

       // memcpy(g_press0, m_press + (itime6 - 1) * pressSize, pressSize * sizeof(float));
		
		if ( (cudaMemcpy(g_press0, m_press + (itime6 - 1) * pressSize, pressSize * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess) )
		{
			printf("pressure data copied \n");
		}
		else
		{
			printf("GPU memory copy error!\n");
		}
		
		
		
		getWindPressure('p');

		

		size_t s_wu_width = m_fWindData->kxu * sizeof(float);
		
		//memcpy(g_uwd0, m_uwd + (itime6 - 1) * windUSize, windUSize * sizeof(float));
		
		if ( (cudaMemcpy(g_uwd0, m_uwd + (itime6 - 1) * windUSize, windUSize * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess) )
		{
			printf("wind U data copied \n");
		}
		else
		{
			printf("GPU memory copy error!\n");
		}
		
		
		getWindPressure('u');

        //memcpy(g_vwd0, m_vwd + (itime6 - 1) * windVSize, windVSize * sizeof(float));
		size_t s_wv_width = m_fWindData->kxv * sizeof(float);
		
		
		if ( (cudaMemcpy(g_vwd0, m_vwd + (itime6 - 1) * windVSize, windVSize * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess) )
		{
			printf("wind V data copied \n");
		}
		else
		{
			printf("GPU memory copy error!\n");
		}
		
		
		getWindPressure('v');

    }
	
	
    float uw, vw, speed, windc;
    int ji, jp1i, jip1, jim1, jm1i, jp1ip1, jm1im1, jp1im1, jm1ip1, jl, jlm1, j1, j2, j3, jli, jlm1i, j1i, j2i, j3i;
	
	
	float ftim = fmodf(timeh6, 1.0f);
	
	//int threadsPerBlock = 64;
	
	//int blocksPerGridJ = (m_height + threadsPerBlock - 1) / threadsPerBlock;
	//int blocksPerGridI = (m_width + threadsPerBlock - 1) / threadsPerBlock;
	
	/*dim3 threadsPerSquareBlock(4, 4);
	
	dim3 numSquareBlocks((m_width + threadsPerSquareBlock.x - 1) / threadsPerSquareBlock.x, (m_height + threadsPerSquareBlock.y - 1) / threadsPerSquareBlock.y);
	
	cudaError_t err = cudaSuccess;*/
	
	
	//surf_and_flux_1<<<numSquareBlocks, threadsPerSquareBlock>>>(ftim, ro_ratio);
	setbuf(stdout,NULL);
	printf("before call\n");
	setbuf(stdout,NULL);
	
	//surf_and_flux_1<<<4, 8>>>(ftim, ro_ratio);
	
	setbuf(stdout,NULL);
	printf("after call\n");
	setbuf(stdout,NULL);
	
	/*,  g_fbu,  g_ffu,  g_fbv,  g_ffv,  g_dum,  g_dvm,  g_d,  g_wusurf,  g_wvsurf,  g_fluxua,  g_fluxva,  g_dx,  g_dy,  g_ua,  g_va,  g_fxf,  g_fyf, g_fxb,  g_fyb);*/

    
    /*
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
    }*/
    
    
    /// HERE SHOULD START A NEW CUDA CALL TO KEEP fluxua fluxva synced
	
	
	setbuf(stdout,NULL);
	printf("b call 2\n");
	setbuf(stdout,NULL);
	
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
	
	
	setbuf(stdout,NULL);
	printf("a call 2\n");
	setbuf(stdout,NULL);
 


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
	{//ADVAVE()
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
				g_fluxua[ji]=g_dy[j]*(.125f*((g_d[jip1]+g_d[ji])*g_ua[jip1]
											  +(g_d[ji]+g_d[jim1])*g_ua[ji])
									  *(g_ua[jip1]+g_ua[ji])
									  -g_d[ji]*2.0f*aam2d*(g_uab[jip1]-g_uab[ji])/g_dx[j]);
				
				
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
				
				g_fluxva[ji]=(.125f*((g_d[ji]+g_d[jm1i])*g_va[ji]
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
				
				
			 	g_fluxva[ji]=g_dx[j]*(.125f*((g_d[jp1i]+g_d[ji])
									       *g_va[jp1i]+(g_d[ji]+g_d[jm1i])*g_va[ji])
									      *(g_va[jp1i]+g_va[ji])
								         -g_d[ji]*2.0f*aam2d*(g_vab[jp1i]-g_vab[ji])/g_dy[j]);
				
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
				
				
				g_fluxua[ji]=(.125f*((g_d[ji]+g_d[jim1])*g_ua[ji]
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
				    *sqrtf(powf(.25e0*(g_uab[ji]+g_uab[jip1]
								  +g_uab[jm1i]+g_uab[jm1ip1]), 2)+g_vab[ji]*g_vab[ji])*g_vab[ji];
				
			}
		}
		
		
		// END ADVAVE();
	}
	
	float alpha =  0.225f;
	float dte = m_fVars->dte;
	
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
			jp1im1 = jp1i - 1;
			jm1ip1 = jm1i + 1;
			
			float uaf1= g_advua[ji]   -0.25f*(g_cor[j]*g_d[ji]*(g_va[jp1i]+g_va[ji])
					                 +g_cor[j]*g_d[jim1]*(g_va[jp1im1]+g_va[jim1]) )
			         +0.5f*grav*g_dy[j]/g_aru[j]*(g_d[ji]+g_d[jim1])
			             *( (1.0f-2.0f*alpha)*(g_el[ji]-g_el[jim1])
							            +alpha*(g_elb[ji]-g_elb[jim1]+g_elf[ji]-g_elf[jim1]) )
			+g_wusurf[ji]-g_wubot[ji];
			
			g_uaf[ji]=
			         ((g_h[ji]+g_elb[ji]+g_h[jim1]+g_elb[jim1])*g_uab[ji]
					                   -4.e0*dte*uaf1)  /(g_h[ji]+g_elf[ji]+g_h[jim1]+g_elf[jim1]);
			
		}
	}
	
	
	for (int j=1; j<m_height; j++ )
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
			
			float vaf1=g_advva[ji]
			+.25f*(  g_cor[j]*g_d[ji]*(g_ua[jip1]+g_ua[ji])
				  +g_cor[j-1]*g_d[jm1i]*(g_ua[jm1ip1]+g_ua[jm1i]) )
			+0.5f*grav*g_dx[j]/g_arv[j]*(g_d[ji]+g_d[jm1i])
			*( (1.0f-2.0f*alpha)*(g_el[ji]-g_el[jm1i])
			  +alpha*(g_elb[ji]-g_elb[jm1i]+g_elf[ji]-g_elf[jm1i]) )
			+ g_wvsurf[ji]-g_wvbot[ji];
			
			g_vaf[ji]= ((g_h[ji]+g_elb[ji]+g_h[jm1i]+g_elb[jm1i])*g_vab[ji]
						-4.0f*dte*vaf1) /(g_h[ji]+g_elf[ji]+g_h[jm1i]+g_elf[jm1i]);
			
		}
		
	}
	
	
	
		/// BCOND 2
	float gae;
	
	for (int j=1; j<(m_height-1); j++ )
	{
		j1 =  j * m_width;
		j2 =  j1 + 1;
		j3 =  j1 + 2;
		
		jl = j1 + m_width -1;
		jlm1 = jl - 1;
		
		if(g_dum[jl] > 0.5f)
		{
			gae = dte*sqrtf(grav*g_h[jl])/g_dx[j];
			
			g_uaf[jl] = gae*g_ua[jlm1]+(1.0f-gae)*g_ua[jl];
		}
		else
		{
			g_uaf[jl] = 0.0f;
		}

		g_vaf[jl]=0.0f;
		
		if(g_dum[j2] > 0.5f)
		{
			gae = dte*sqrtf(grav*g_h[j2])/g_dx[j];
			g_uaf[j2]=gae*g_ua[j3]+(1.-gae)*g_ua[j2];
		}
		else
		{
			g_uaf[j2]=0.0f;
		}
		
		g_uaf[j1]=g_uaf[j2];
		g_vaf[j1]=0.0;
		
	}
	
	
	

	for (int i=1; i<(m_width-1); i++ )
	{
		jli = m_width * (m_height - 1) + i;
		jlm1i = jli - m_width;
		
		j1i = i;
		
		j2i = m_width + j1i;
		
		j3i = m_width + j2i;
		
		
		if (g_dvm[jli] > 0.5f)
		{
			gae = dte*sqrtf(grav*g_h[jli])/g_dy[m_height-1];
			
			g_vaf[jli] = gae*g_va[jlm1i]+(1.0f-gae)*g_va[jli];
		}
		else
		{
			g_vaf[jli]=0.0f;
		}

		g_uaf[jli]=0.0;

		if (g_dvm[j2i] > 0.5f)
		{
			gae=dte*sqrtf(grav*g_h[j2i])/g_dy[0];
			
			g_vaf[j2i]=gae*g_va[j3i]+(1.-gae)*g_va[j2i];
		}
		else
		{
			g_vaf[j2i]=0.0f;
		}
		

		g_vaf[j1i]=g_vaf[j1i];
		g_uaf[j1i]=0.0f;
	}
	
	/// must separate cuda calls here
	
	for (int j=1; j<m_height; j++ )
	{
		for (int i=1; i<m_width; i++ )
		{
			ji = j * m_width + i;
			
			g_uaf[ji] = g_uaf[ji] * g_dum[ji];
			g_vaf[ji] = g_vaf[ji] * g_dvm[ji];
		}
	}
	// END BCOND 2
	
	
	
	float vmaxl = 100.0f;
	
	
	float tpsmax = 0.0f;
	
	int imax = 0;
	int jmax = 0;
	
	for (int j=1; j<m_height; j++ )
	{
		for (int i=1; i<m_width; i++ )
		{
			ji = j * m_width + i;
			
			g_tps[ji] = sqrtf(g_uaf[ji]*g_uaf[ji] + g_vaf[ji]*g_vaf[ji]);
			
			if (g_tps[ji] > tpsmax)
			{
				tpsmax = g_tps[ji];
				imax = i;
				jmax = j;
			}
		}
	}
	
	
	if (tpsmax > vmaxl)
	{
		setbuf(stdout,NULL);
		
		printf("vamax>vmax!!! at i=%d, j=%d \n", imax,jmax);
		
		exit(-1);
	}

	
	
	float smoth = 0.10f;
	
	
	for (int j=1; j<m_height; j++ )
	{
		for (int i=1; i<m_width; i++ )
		{
			ji = j * m_width + i;
			
			g_ua[ji]=g_ua[ji]+0.5f*smoth*(g_uab[ji]-2.0f*g_ua[ji]+g_uaf[ji]);
			g_va[ji]=g_va[ji]+0.5f*smoth*(g_vab[ji]-2.0f*g_va[ji]+g_vaf[ji]);
			g_el[ji]=g_el[ji]+0.5f*smoth*(g_elb[ji]-2.0f*g_el[ji]+g_elf[ji]);
			g_elb[ji]=g_el[ji];  // OP
			g_el[ji]=g_elf[ji];  // OP
			g_d[ji]=g_h[ji]+g_el[ji];
			g_uab[ji]=g_ua[ji];  // OP
			g_ua[ji]=g_uaf[ji];  // OP
			g_vab[ji]=g_va[ji];  // OP
			g_va[ji]=g_vaf[ji];  // OP
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
		pk = g_uwd0;
		
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
		pk = g_vwd0;
		
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
		pk = g_press0;
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






int KaspyCycler::init_device()
{
	int device_count = 0;
	
	if (m_gpu_device >= 0)
	{
		// already initialized
		printf("CUDA device is already initiaized\n");
		
		return m_gpu_device;
	}
	
	
	cudaGetDeviceCount(&device_count);
	
	for (int i = 0 ; i < device_count ; ++i)
	{
		cudaDeviceProp properties;
		cudaGetDeviceProperties(&properties, i);
		
		if (properties.major > 1 || (properties.major == 1 && properties.minor >= 1))
		{
			m_gpu_device = i;
			
			printf("Running on GPU %d (%s) \n",i ,properties.name);
			break;
		}
		else
		{
			printf("GPU %d (%s) does not support CUDA Dynamic Parallelism\n", i ,properties.name);
		}
	}
	
	
	if (m_gpu_device == -1)
	{
		printf("No suitable device found!\n");
		return m_gpu_device;
	}
	
	if (cudaSetDevice(m_gpu_device) == cudaSuccess)
	{
		printf("device set OK\n");
	}
	else
	{
		printf("unable to set device!\n");
		m_gpu_device = -1;
	}
	
	
	
	// Allocate GPU memory.
	if ( (cudaMallocManaged((void **)&g_fbu, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_fbv, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_ffu, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_ffv, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_fb, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_ff, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_fxb, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_fxf, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_fyb, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_fyf, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_wusurf, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_wvsurf, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_dum, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_dvm, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_d, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		
		&& (cudaMallocManaged((void **)&g_fluxua, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_fluxva, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		
		&& (cudaMallocManaged((void **)&g_ua, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_va, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_uab, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_vab, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_uaf, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_vaf, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_el, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_elb, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_elf, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_fsm, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_tps, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_advua, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_advva, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_wubot, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_wvbot, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_cbc, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_h, m_height*m_width * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		
		&& (cudaMallocManaged((void **)&g_press0,  m_fWindData->ky *  m_fWindData->kx * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		
		&& (cudaMallocManaged((void **)&g_uwd0, m_fWindData->kyu * m_fWindData->kxu * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_vwd0, m_fWindData->kyv * m_fWindData->kxv * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)

		
		&& (cudaMallocManaged((void **)&g_cor, m_height * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_aru, m_height * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_arv, m_height * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_dx, m_height * sizeof(float), cudaMemAttachGlobal) == cudaSuccess)
		&& (cudaMallocManaged((void **)&g_dy, m_height * sizeof(float), cudaMemAttachGlobal) == cudaSuccess))
	{
		printf("GPU memory allocated\n");
		
	}
	else
	{
		printf("GPU memory allocation error!\n");
		deinit_device();
		return m_gpu_device;
	}
	
	
	
	
	return m_gpu_device;
}


void KaspyCycler::deinit_device()
{
	if (m_gpu_device >=0)
	{
		
		if (g_fbu)
		{
			cudaFree(g_fbu);
		}
		
		if (g_fbv)
		{
			cudaFree(g_fbv);
		}
		
		if (g_ffu)
		{
			cudaFree(g_ffu);
		}
		
		if (g_ffv)
		{
			cudaFree(g_ffv);
		}
		
		
		if (g_fb)
		{
			cudaFree(g_fb);
		}
		
		if (g_ff)
		{
			cudaFree(g_ff);
		}
		
		
		if (g_fxb)
		{
			cudaFree(g_fxb);
		}
		
		if (g_fxf)
		{
			cudaFree(g_fxf);
		}
		
		if (g_fyb)
		{
			cudaFree(g_fyb);
		}
		
		if (g_fyf)
		{
			cudaFree(g_fyf);
		}
		
		if (g_wusurf)
		{
			cudaFree(g_wusurf);
		}
		
		if (g_wvsurf)
		{
			cudaFree(g_wvsurf);
		}
		
		if (g_dum)
		{
			cudaFree(g_dum);
		}
		
		if (g_dvm)
		{
			cudaFree(g_dvm);
		}
		
		if (g_d)
		{
			cudaFree(g_d);
		}
		
		if (g_fluxua)
		{
			cudaFree(g_fluxua);
		}
		
		if (g_fluxva)
		{
			cudaFree(g_fluxva);
		}
		
		if (g_ua)
		{
			cudaFree(g_ua);
		}
		
		if (g_va)
		{
			cudaFree(g_va);
		}
		
		if (g_uaf)
		{
			cudaFree(g_uaf);
		}
		
		if (g_vaf)
		{
			cudaFree(g_vaf);
		}
		
		if (g_uab)
		{
			cudaFree(g_uab);
		}
		
		if (g_vab)
		{
			cudaFree(g_vab);
		}
		
		if (g_el)
		{
			cudaFree(g_el);
		}
		
		if (g_elb)
		{
			cudaFree(g_elb);
		}
		
		if (g_elf)
		{
			cudaFree(g_elf);
		}
		
		if (g_fsm)
		{
			cudaFree(g_fsm);
		}
		
		if (g_tps)
		{
			cudaFree(g_tps);
		}
		
		if (g_advua)
		{
			cudaFree(g_advua);
		}
		
		if (g_advva)
		{
			cudaFree(g_advva);
		}
		
		if (g_aru)
		{
			cudaFree(g_aru);
		}
		
		if (g_arv)
		{
			cudaFree(g_arv);
		}
		
		if (g_wubot)
		{
			cudaFree(g_wubot);
		}
		
		if (g_wvbot)
		{
			cudaFree(g_wvbot);
		}
		
		if (g_cbc)
		{
			cudaFree(g_cbc);
		}
		
		if (g_h)
		{
			cudaFree(g_h);
		}
		
		if (g_cor)
		{
			cudaFree(g_cor);
		}
		
		if (g_press0)
		{
			cudaFree(g_press0);
		}
		
		
		if (g_uwd0)
		{
			cudaFree(g_uwd0);
		}
		
		if (g_vwd0)
		{
			cudaFree(g_vwd0);
		}
		
		if (cudaDeviceReset() == cudaSuccess)
		{
			printf("GPU device reset ok\n");
		}
		
		m_gpu_device = -1;
	}
}





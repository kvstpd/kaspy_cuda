//
//  KaspyCycler.cpp
//  kaspy_cuda
//
//  Created by Andrei Koulikov on 24.05.16.
//
//

#define CUB_STDERR

#include "InitValues.h"

#include "KaspyCycler.h"


#include "cuda_runtime.h"
#include "device_launch_parameters.h"


#include "cub/cub.cuh"


using namespace cub;


extern InitValues * initValues;

extern "C" void WRITEGRD(int * NX, int * NY, int * NDX, float * Z, float * XMI, float * XMA, float * YMI, float * YMA,const char * NAME);





void getbicubic(int nx, int ny, int nd, float * z, float * c);
void bcucof(float * y,float * y1,float * y2, float * y12,float d1,float d2,float * cc);


float grav = 9.806;

CachingDeviceAllocator  g_allocator(true);



int * g_stations_x = 0;
int * g_stations_y = 0;

float * g_station_elves = 0;

__device__ int * dev_stations_x = 0;
__device__ int * dev_stations_y = 0;

__device__ int dev_nstations = 0;

__device__ float * dev_station_elves = 0;



__device__ float * dev_fbu = 0;
__device__ float * dev_fbv = 0;
__device__ float * dev_ffu = 0;
__device__ float * dev_ffv = 0;

__device__ float * dev_fxb = 0;
__device__ float * dev_fxf = 0;
__device__ float * dev_fyb = 0;
__device__ float * dev_fyf = 0;

__device__ float * dev_fb = 0;
__device__ float * dev_ff = 0;

__device__ float * dev_wusurf = 0;
__device__ float * dev_wvsurf = 0;

__device__ float * dev_dum = 0;
__device__ float * dev_dvm = 0;

__device__ float * dev_d = 0;
__device__ float * dev_dx = 0;
__device__ float * dev_dy = 0;

__device__ float * dev_fluxua = 0;
__device__ float * dev_fluxva = 0;

__device__ float * dev_ua = 0;
__device__ float * dev_va = 0;

__device__ float * dev_uab = 0;
__device__ float * dev_vab = 0;

__device__ float * dev_uaf = 0;
__device__ float * dev_vaf = 0;


__device__ float * dev_el = 0;
__device__ float * dev_elf = 0;
__device__ float * dev_elb = 0;

__device__ float * dev_elf_r = 0;

__device__ float * dev_fsm = 0;

__device__ float * dev_tps = 0;


__device__ float * dev_advua = 0;
__device__ float * dev_advva = 0;

__device__ float * dev_aru = 0;
__device__ float * dev_arv = 0;

__device__ float * dev_wubot = 0;
__device__ float * dev_wvbot = 0;
__device__ float * dev_cbc = 0;

__device__ float * dev_cor = 0;

__device__ float * dev_h = 0;

__device__ float * dev_press = 0;
__device__ float * dev_uwd = 0;
__device__ float * dev_vwd = 0;

//__device__ float * dev_press0 = 0;
//__device__ float * dev_uwd0 = 0;
//__device__ float * dev_vwd0 = 0;

__device__ float * dev_art = 0;


__device__ float * dev_p = 0;
__device__ float * dev_pk = 0;
__device__ float * dev_px = 0;
__device__ float * dev_py = 0;


__device__ float * dev_temp = 0;




__device__ float * dev_sel = 0;
__device__ float * dev_ssel = 0;
__device__ float * dev_sfel = 0;
__device__ float * dev_sfa = 0;
__device__ float * dev_ssfa = 0;
__device__ float * dev_sfar = 0;
__device__ float * dev_ssfar = 0;
__device__ float * dev_sfelr = 0;
__device__ float * dev_su = 0;
__device__ float * dev_sv = 0;
__device__ float * dev_ssu = 0;
__device__ float * dev_ssv = 0;
__device__ float * dev_ssuv = 0;
__device__ float * dev_ssue = 0;
__device__ float * dev_ssve = 0;



float * g_sel = 0;
float * g_ssel = 0;
float * g_sfel = 0;
float * g_sfa = 0;
float * g_ssfa = 0;
float * g_sfar = 0;
float * g_ssfar = 0;
float * g_sfelr = 0;
float * g_su = 0;
float * g_sv = 0;
float * g_ssu = 0;
float * g_ssv = 0;
float * g_ssuv = 0;
float * g_ssue = 0;
float * g_ssve = 0;


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

float * g_elf_r = 0;

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

float * g_press = 0;
float * g_uwd = 0;
float * g_vwd = 0;

float * g_press0 = 0;
float * g_uwd0 = 0;
float * g_vwd0 = 0;

float * g_art = 0;


__constant__ __device__  float dev_grav = 9.806f;
__constant__ __device__  float dev_ro_ratio = 1.29f/1020.0f;

__constant__ __device__  int  dev_width;
__constant__ __device__  int  dev_height;
__constant__ __device__  int  dev_widthm1;
__constant__ __device__  int  dev_heightm1;

__constant__ __device__ int dev_ewidth;

__constant__ __device__ float dev_dte;
__constant__ __device__ float dev_dte2;
__constant__ __device__ float dev_aam2d;

__constant__ __device__ float dev_tide_l = 0.0f;

__constant__ __device__ float dev_alpha = 0.225f;

__constant__ __device__ float dev_vmaxl = 100.0f;;

__device__ int dev_should_stop = 0;

__constant__ __device__ float dev_smoth = 0.10f;

__constant__ __device__ float dev_xmi;// = m_fVars->xmi;
__constant__ __device__ float dev_xma;// = m_fVars->xma;
__constant__ __device__ float dev_ymi;// = m_fVars->ymi;
__constant__ __device__ float dev_yma;// = m_fVars->yma;



__constant__ __device__ float dev_c1 = 3.1415926/180.0;
__constant__ __device__ float dev_c2 = 111111.0f;


__constant__ __device__ float dev_wt[] = {
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


__device__ float dev_pkk[50 * 50];
__device__ float dev_c[50 * 50 * 4 * 4];



__device__ void dev_bcucof(float * y,float * y1,float * y2, float * y12,float d1,float d2,float * cc)
{
	float xx;
	float cl[16];
	
	float x[16];
	
	
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
			xx += dev_wt[i + k*16] * x[k];
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


__global__ void dev_make_p(int nx, int ny, int kx, int ky, float dx, float dy, float dkx, float dky, float xki, float yki)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	
	
	if (i < nx && j < ny)
	{
		float y = dev_ymi + j*dy;
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
		
		
		float x = dev_xmi + i * dx;
		
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
			ay = t*ay+((dev_c[j0 * 800 + i0 * 16 + 3 * 4 + k] * u + dev_c[j0 * 800 + i0 * 16 + 2 * 4 + k])*u
					   + dev_c[j0 * 800 + i0 * 16 + 1 * 4 + k])*u + dev_c[j0 * 800 + i0 * 16 + 0 * 4 + k];
		}
		
		
		if (dev_px != 0)
		{
			for (int k=3; k>=0; k-- )
			{
				a2 = t*a2 + (3.0f*dev_c[j0 * 800 + i0 * 16 + 3 * 4 + k]*u
							 + 2.0f*dev_c[j0 * 800 + i0 * 16 + 2 * 4 + k])*u+dev_c[j0 * 800 + i0 * 16 + 1 * 4 + k];
				
				a1 = u*a1 + (3.0f*dev_c[j0 * 800 + i0 * 16 + k * 4 + 3]*t +
							 2.0f*dev_c[j0 * 800 + i0 * 16 + k * 4 + 2])*t+dev_c[j0 * 800 + i0 * 16 + k * 4 + 1];
				
			}
			
			a1 = a1/dkx/dev_c2/cosf(dev_c1*y);
			a2 = a2/dky/dev_c2;
			
			dev_px[ji] = a1;
			dev_py[ji] = a2;
		}
		
		dev_p[ji] = ay;
		
	}
}




__global__ void dev_bicubic(int nx, int ny, int nd)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	
	float d1 = 1.0f;
	float d2 = 1.0f;
	
	float y[4];
	float y1[4];
	float y2[4];
	float y12[4];
	float cc[4][4];
	
	
	
	
	if (i > 0 && j > 0 && i < (nx - 2) && j < (ny - 2))
	{
		
		y[0] = dev_pkk[j * nd + i];
		y[1] = dev_pkk[j * nd + i + 1];
		y[2] = dev_pkk[(j+1) * nd + i + 1];
		y[3] = dev_pkk[(j+1) * nd + i];
		
		y1[0] = 0.5f * (dev_pkk[j * nd + i + 1] - dev_pkk[j * nd + i - 1]);
		y1[3] = 0.5f * (dev_pkk[(j+1) * nd + i + 1] - dev_pkk[(j+1) * nd + i - 1]);
		y1[1] = 0.5f * (dev_pkk[j * nd + i + 2] - dev_pkk[j * nd + i]);
		y1[2] = 0.5f * (dev_pkk[(j+1) * nd + i + 2] - dev_pkk[(j+1) * nd + i]);
		
		
		y2[0] = 0.5f * (dev_pkk[(j+1) * nd + i] - dev_pkk[(j-1) * nd + i]);
		y2[1] = 0.5f * (dev_pkk[(j+1) * nd + i + 1] - dev_pkk[(j-1) * nd + i + 1]);
		y2[2] = 0.5f * (dev_pkk[(j+2) * nd + i + 1] - dev_pkk[(j) * nd + i + 1]);
		y2[3] = 0.5f * (dev_pkk[(j+2) * nd + i] - dev_pkk[j * nd + i]);
		
		
		y12[0] = 0.25f * (dev_pkk[(j+1) * nd + i + 1] - dev_pkk[(j-1) * nd + i + 1]
						  - dev_pkk[(j+1) * nd + i - 1] + dev_pkk[(j-1) * nd + i - 1]);
		y12[1] = 0.25f * (dev_pkk[(j+1) * nd + i + 2] - dev_pkk[(j-1) * nd + i + 2]
						  - dev_pkk[(j+1) * nd + i] + dev_pkk[(j-1) * nd + i]);
		y12[2] = 0.25f * (dev_pkk[(j+2) * nd + i + 2] - dev_pkk[(j) * nd + i + 2]
						  - dev_pkk[(j+2) * nd + i] + dev_pkk[j * nd + i]);
		y12[3] = 0.25f * (dev_pkk[(j+2) * nd + i + 1] - dev_pkk[(j) * nd + i + 1]
						  - dev_pkk[(j+2) * nd + i -1] + dev_pkk[(j) * nd + i -1]);
		
		
		dev_bcucof(&y[0],&y1[0],&y2[0],&y12[0],d1,d2,&cc[0][0]);
		
		for (int k=0; k<4; k++ )
		{
			for (int l=0; l<4; l++ )
			{
				//printf("\nk is %d l is %d\n", k, l);
				dev_c[(j-1)* 800 + (i-1) * 16 + l * 4 + k ] = cc[l][k];
			}
		}
	}
}




__global__ void dev_pkk_ij(int kx, int ky, int kd)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	if (i > 0 && j > 0 && i <= kx && j <= ky)
	{
		dev_pkk[j * 50 + i] = dev_pk[(j - 1) * kd + i - 1];
	}
}


__global__ void dev_pkk_j(int kx, int ky)
{
	int j = blockDim.x * blockIdx.x + threadIdx.x;
	
	if (j > 0 && j <= ky)
	{
		dev_pkk[j*50+0] = 2.0f*dev_pkk[j*50+1] - dev_pkk[j*50+2];
		dev_pkk[j*50+kx+1] = 2.0f*dev_pkk[j*50+kx] - dev_pkk[j*50+kx-1];
	}
}


__global__ void dev_pkk_i(int kx, int ky)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	
	if (i > 0 && i <= kx+1)
	{
		dev_pkk[0*50+i] = 2.0f*dev_pkk[1*50+i] - dev_pkk[2*50+i];
		dev_pkk[(ky+1)*50+i] = 2.0f*dev_pkk[ky*50+i] - dev_pkk[(ky-1)*50+i];
	}
}


/**/

__global__ void surf_and_flux_1(float ftim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	int ji = j * dev_width + i;
	int jp1i = ji + dev_width;
	int jip1 = ji + 1;
	int jim1 = ji - 1;
	int jm1i = ji - dev_width;
	
	float btim = 1.0f - ftim;
	
	
	if (i < dev_widthm1 && j < dev_heightm1)
	{
		float uw = btim * (dev_fbu[ji]) + ftim * (dev_ffu[ji]);
		float vw = btim * (dev_fbv[ji]) + ftim * (dev_ffv[ji]);
		
		float speed =  hypotf(uw, vw); //sqrtf(uw*uw + vw*vw);
		float windc = 0.001f * (0.8f + speed * 0.065f) * dev_ro_ratio * speed;
		
		dev_wusurf[ji] = -windc * uw *
		0.25f * (dev_dum[jp1i]+dev_dum[jip1]+dev_dum[jim1]+dev_dum[jm1i])
		+ 0.5f * (dev_d[ji] + dev_d[jim1]) * (btim * dev_fxb[ji] + ftim * dev_fxf[ji]);
		
		dev_wvsurf[ji] = -windc * vw *
		0.25f * (dev_dvm[jp1i]+dev_dvm[jip1]+dev_dvm[jim1]+dev_dvm[jm1i])
		+ 0.5f * (dev_d[ji] + dev_d[jm1i]) * (btim * dev_fyb[ji] + ftim * dev_fyf[ji]);
	}
	
	if (i < dev_width && j < dev_height)
	{
		dev_fluxua[ji] = 0.25f * (dev_d[ji] + dev_d[jim1]) * (dev_dy[j] + dev_dy[j] ) * dev_ua[ji];
		dev_fluxva[ji] = 0.25f * (dev_d[ji] + dev_d[jm1i]) * (dev_dx[j] + dev_dx[j-1] ) * dev_va[ji];
	}
}


__global__ void elf_and_flux_2()
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	int ji = j * dev_width + i;
	int jp1i = ji + dev_width;
	int jip1 = ji + 1;
	
	if (i > 0 && j > 0 && i < dev_widthm1 && j < dev_heightm1)
	{
		dev_elf[ji] = dev_elb[ji] - dev_dte2 *
		(dev_fluxua[jip1] - dev_fluxua[ji] + dev_fluxva[jp1i] - dev_fluxva[ji]) / dev_art[j];
	}

}


__global__ void bcond_1_j()
{
	int j = blockDim.x * blockIdx.x + threadIdx.x;
	
	if (j > 0 && j < dev_height)
	{
		dev_elf[j * dev_width + 1] = dev_tide_l;
		dev_elf[j * dev_width + dev_width - 2] = dev_tide_l;
		
		dev_elf[j * dev_width] = dev_tide_l;
		dev_elf[j * dev_width + dev_width - 1] = dev_tide_l;
	}
}

__global__ void bcond_1_i()
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	
	if (i > 0 && i< dev_width)
	{
		dev_elf[i] =  dev_elf[i + dev_width];
		
		dev_elf[i + dev_width * (dev_height - 1)  ] =  dev_elf[i + dev_width * (dev_height - 2)];
	}
}


__global__ void bcond_1_ji()
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	int ji = j * dev_width + i;
	
	
	if (i > 0 && j > 0 && i < dev_width && j < dev_height)
	{
		dev_elf[ji] *= dev_fsm[ji];
	}
	
}

__global__ void uaf_and_vaf_3()
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	int ji = j * dev_width + i;
	
 	int jp1i = ji + dev_width;
 	int jip1 = ji + 1;
 	int jim1 = ji - 1;
 	int jm1i = ji - dev_width;
 	//int jm1im1 = jm1i  - 1;
 	int jp1im1 = jp1i - 1;
 	int jm1ip1 = jm1i + 1;
	
	
	

	if (i > 0 && j > 0)
	{
		if (i < dev_width && j < dev_heightm1)
		{

			
			
			float uaf1= dev_advua[ji]   -0.25f*(dev_cor[j]*dev_d[ji]*(dev_va[jp1i]+dev_va[ji])
												+dev_cor[j]*dev_d[jim1]*(dev_va[jp1im1]+dev_va[jim1]) )
			+0.5f*dev_grav*dev_dy[j]/dev_aru[j]*(dev_d[ji]+dev_d[jim1])
			*( (1.0f-2.0f*dev_alpha)*(dev_el[ji]-dev_el[jim1])
			  +dev_alpha*(dev_elb[ji]-dev_elb[jim1]+dev_elf[ji]-dev_elf[jim1]) )
			+dev_wusurf[ji]-dev_wubot[ji];
			
			dev_uaf[ji]=
			((dev_h[ji]+dev_elb[ji]+dev_h[jim1]+dev_elb[jim1])*dev_uab[ji]
			 -4.e0*dev_dte*uaf1)  /(dev_h[ji]+dev_elf[ji]+dev_h[jim1]+dev_elf[jim1]);
		}
		
		if (i < dev_widthm1 && j < dev_height)
		{
			float vaf1=dev_advva[ji]
			+.25f*(  dev_cor[j]*dev_d[ji]*(dev_ua[jip1]+dev_ua[ji])
				   +dev_cor[j-1]*dev_d[jm1i]*(dev_ua[jm1ip1]+dev_ua[jm1i]) )
			+0.5f*dev_grav*dev_dx[j]/dev_arv[j]*(dev_d[ji]+dev_d[jm1i])
			*( (1.0f-2.0f*dev_alpha)*(dev_el[ji]-dev_el[jm1i])
			  +dev_alpha*(dev_elb[ji]-dev_elb[jm1i]+dev_elf[ji]-dev_elf[jm1i]) )
			+ dev_wvsurf[ji]-dev_wvbot[ji];
			
			dev_vaf[ji]= ((dev_h[ji]+dev_elb[ji]+dev_h[jm1i]+dev_elb[jm1i])*dev_vab[ji]
						  -4.0f*dev_dte*vaf1) /(dev_h[ji]+dev_elf[ji]+dev_h[jm1i]+dev_elf[jm1i]);
			
		}
		
	}
	
}


__global__ void bcond_2_j()
{
	int j = blockDim.x * blockIdx.x + threadIdx.x;
	
	int j1 =  j * dev_width;
	int j2 =  j1 + 1;
	int j3 =  j1 + 2;
	int jl = j1 + dev_widthm1;
	int jlm1 = jl - 1;
	
	float gae;
	
	if (j > 0 && j < dev_heightm1)
	{
		if(dev_dum[jl] > 0.5f)
		{
			gae = dev_dte*sqrtf(dev_grav*dev_h[jl])/dev_dx[j];
			
			dev_uaf[jl] = gae*dev_ua[jlm1]+(1.0f-gae)*dev_ua[jl];
		}
		else
		{
			dev_uaf[jl] = 0.0f;
		}
		
		dev_vaf[jl]=0.0f;
		
		if(dev_dum[j2] > 0.5f)
		{
			gae = dev_dte*sqrtf(dev_grav*dev_h[j2])/dev_dx[j];
			dev_uaf[j2]=gae*dev_ua[j3]+(1.0f-gae)*dev_ua[j2];
		}
		else
		{
			dev_uaf[j2]=0.0f;
		}
		
		dev_uaf[j1]=dev_uaf[j2];
		dev_vaf[j1]=0.0f;

	}
}

__global__ void bcond_2_i()
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	
	int jli = dev_width * (dev_heightm1) + i;
	int jlm1i = jli - dev_width;
	
	int j1i = i;
	
	int j2i = dev_width + j1i;
	
	int j3i = dev_width + j2i;
	
	float gae;
	
	if (i > 0 && i< dev_widthm1)
	{
		if (dev_dvm[jli] > 0.5f)
		{
			gae = dev_dte * sqrtf(dev_grav * dev_h[jli]) / dev_dy[dev_heightm1];
			
			dev_vaf[jli] = gae * dev_va[jlm1i]+(1.0f-gae)*dev_va[jli];
		}
		else
		{
			dev_vaf[jli]=0.0f;
		}
		
		dev_uaf[jli]=0.0;
		
		if (dev_dvm[j2i] > 0.5f)
		{
			gae=dev_dte*sqrtf(dev_grav*dev_h[j2i])/dev_dy[0];
			
			dev_vaf[j2i]=gae*dev_va[j3i]+(1.-gae)*dev_va[j2i];
		}
		else
		{
			dev_vaf[j2i]=0.0f;
		}
		
		
		dev_vaf[j1i]=dev_vaf[j1i];
		dev_uaf[j1i]=0.0f;
	}
}


__global__ void bcond_2_ji()
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	int ji = j * dev_width + i;
	
	
	if (i > 0 && j > 0 && i < dev_width && j < dev_height)
	{
		dev_uaf[ji] = dev_uaf[ji] * dev_dum[ji];
		dev_vaf[ji] = dev_vaf[ji] * dev_dvm[ji];
	}
 

}

__global__ void tps_and_other_arrays_4()
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	int ji = j * dev_width + i;
	
	
	if (i > 0 && j > 0 && i < dev_width && j < dev_height)
	{
		dev_tps[ji] = hypotf(dev_uaf[ji], dev_vaf[ji]);  ///sqrtf(dev_uaf[ji]*dev_uaf[ji] + dev_vaf[ji]*dev_vaf[ji]);
		
		/*if (dev_tps[ji] > dev_vmaxl)
		{
			dev_should_stop = 1;
		}*/
		
		dev_ua[ji]=dev_ua[ji]+0.5f*dev_smoth*(dev_uab[ji]-2.0f*dev_ua[ji]+dev_uaf[ji]);
		dev_va[ji]=dev_va[ji]+0.5f*dev_smoth*(dev_vab[ji]-2.0f*dev_va[ji]+dev_vaf[ji]);
		dev_el[ji]=dev_el[ji]+0.5f*dev_smoth*(dev_elb[ji]-2.0f*dev_el[ji]+dev_elf[ji]);
		//dev_elb[ji]=dev_el[ji];  // OP
		//dev_el[ji]=dev_elf[ji];  // OP
		dev_d[ji]=dev_h[ji]+dev_elf[ji];
		//dev_uab[ji]=dev_ua[ji];  // OP
		//dev_ua[ji]=dev_uaf[ji];  // OP
		//dev_vab[ji]=dev_va[ji];  // OP
		//dev_va[ji]=dev_vaf[ji];  // OP
		
	}
	
	
}

__global__ void swap_arrays_5()
{
	float * t;
	
	t = dev_elb;
	dev_elb = dev_el;
	dev_el = dev_elf;
	dev_elf = t;
	
	
	t = dev_uab;
	dev_uab = dev_ua;
	dev_ua = dev_uaf;
	dev_uaf = t;
	
	t = dev_vab;
	dev_vab = dev_va;
	dev_va = dev_vaf;
	dev_vaf = t;
}


__global__ void adv_fluxes_1()
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	int ji = j * dev_width + i;
	int jip1 = ji + 1;
	int jim1 = ji - 1;
	int jm1i = ji - dev_width;
	int jm1im1 = jm1i  - 1;
 
	
	
	if (i > 0 && j > 0)
	{
		if (i < dev_widthm1 && j < dev_height)
		{
			dev_fluxua[ji]=dev_dy[j]*(.125f*((dev_d[jip1]+dev_d[ji])*dev_ua[jip1]
										 +(dev_d[ji]+dev_d[jim1])*dev_ua[ji])
								  *(dev_ua[jip1]+dev_ua[ji])
								  -dev_d[ji]*2.0f*dev_aam2d*(dev_uab[jip1]-dev_uab[ji])/dev_dx[j]);
		}
		
		
		if (i < dev_width && j < dev_height)
		{
			dev_tps[ji]=(dev_d[ji]+dev_d[jim1]+dev_d[jm1i]+dev_d[jm1im1])
			*dev_aam2d
			*((dev_uab[ji]-dev_uab[jm1i])
			  /(4.0f*dev_dy[j])
			  +(dev_vab[ji]-dev_vab[jim1])
			  /(4.0f*dev_dx[j]) );
			
			dev_fluxva[ji]=(.125f*((dev_d[ji]+dev_d[jm1i])*dev_va[ji]
								 +(dev_d[jim1]+dev_d[jm1im1])*dev_va[jim1])
						  *(dev_ua[ji]+dev_ua[jm1i])
						  -dev_tps[ji])*dev_dx[j];
			
		}

	}
}




__global__ void adv_advua_1()
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	int ji = j * dev_width + i;
	int jim1 = ji - 1;
 	int jp1i = ji + dev_width;
	
	if (i > 0 && j > 0 && i < dev_widthm1 && j < dev_heightm1)
	{
		dev_advua[ji]=(dev_fluxua[ji]-dev_fluxua[jim1]
					 +dev_fluxva[jp1i]-dev_fluxva[ji])/dev_aru[j];
		
	}
}

__global__ void adv_fluxes_2()
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	int ji = j * dev_width + i;
	int jim1 = ji - 1;
	int jp1i = ji + dev_width;
	int jm1i = ji - dev_width;
	int jm1im1 = jm1i  - 1;
 
	
	
	if (i > 0 && j > 0)
	{
		if (i < dev_width && j < dev_heightm1)
		{
			dev_fluxva[ji]=dev_dx[j]*(.125f*((dev_d[jp1i]+dev_d[ji])
										 *dev_va[jp1i]+(dev_d[ji]+dev_d[jm1i])*dev_va[ji])
								  *(dev_va[jp1i]+dev_va[ji])
								  -dev_d[ji]*2.0f*dev_aam2d*(dev_vab[jp1i]-dev_vab[ji])/dev_dy[j]);
		}
		
		
		if (i < dev_width && j < dev_height)
		{
			dev_fluxua[ji]=(.125f*((dev_d[ji]+dev_d[jim1])*dev_ua[ji]
								 +(dev_d[jm1i]+dev_d[jm1im1])*dev_ua[jm1i])*
						  (dev_va[jim1]+dev_va[ji])
						  -dev_tps[ji])*dev_dy[j];
			
		}
		
	}
}


__global__ void adv_advva_2()
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	int ji = j * dev_width + i;
	int jip1 = ji + 1;
	int jm1i = ji - dev_width;
	
	if (i > 0 && j > 0 && i < dev_widthm1 && j < dev_heightm1)
	{
		dev_advva[ji]=(dev_fluxua[jip1]-dev_fluxua[ji]
					 +dev_fluxva[ji]-dev_fluxva[jm1i])/dev_arv[j];
		
	}
}

__global__ void adv_bot_3()
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	int ji = j * dev_width + i;
 	int jp1i = ji + dev_width;
 	int jip1 = ji + 1;
 	int jim1 = ji - 1;
 	int jm1i = ji - dev_width;
 
 	int jp1im1 = jp1i - 1;
 	int jm1ip1 = jm1i + 1;
	
	if (i > 0 && j > 0 && i < dev_widthm1 && j < dev_heightm1)
	{
		dev_wubot[ji]=-0.5f*(dev_cbc[ji]+dev_cbc[jim1])
		* hypotf(dev_uab[ji], 0.25f*(dev_vab[ji] +dev_vab[jp1i]+dev_vab[jim1]+dev_vab[jp1im1]))
		*dev_uab[ji];
		
		dev_wvbot[ji]=-0.5f*(dev_cbc[ji]+dev_cbc[jm1i])
		* hypotf(.25f*(dev_uab[ji]+dev_uab[jip1]+dev_uab[jm1i]+dev_uab[jm1ip1]), dev_vab[ji])
		* dev_vab[ji];
	}
}



__global__ void dev_fill_station_data(int khour)
{
	int n = blockDim.x * blockIdx.x + threadIdx.x;
	
	int ji;
	
	if (n < dev_nstations)
	{
		ji = (dev_stations_y[n] - 1) * dev_width + dev_stations_x[n] - 1;
		
		dev_station_elves[(khour-1) * dev_nstations + n] = dev_el[ji];
		
		//printf("set st to %f", dev_el[ji]);
		
		/*if(n==0)
		{
			printf("khour is %d elf is %f\n", khour, dev_station_elves[khour * dev_nstations + n]);
		}*/
	}
	

}


__global__ void dev_statistics_1(float ftim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	int ji = j * dev_width + i;
	
	if (i < dev_width && j < dev_height)
	{
		dev_sel[ji] += dev_el[ji];
		dev_ssel[ji] += dev_el[ji]*dev_el[ji];
	}
	
	//float btim = 1.0f - ftim;
	
	//float fa = (btim * dev_fb[ji] + ftim * dev_ff[ji] - 100.0f)/10.0f;
	
	
	
	//dev_sfa[ji] += fa;
	//dev_ssfa[ji] += fa*fa;
	

	//dev_sfel[ji] += dev_el[ji] * fa;
	
	
}

__global__ void dev_statistics_finalize(float nstat)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	int ji = j * dev_width + i;
	
	
	
	if (i < dev_width && j < dev_height)
	{
		//dev_sfa[ji] /= nstat;
		dev_sel[ji] /= nstat;
		
		//dev_sfar[ji] /= nstat;
		
		float prev_ssel = dev_ssel[ji];
		
		//dev_ssfa[ji] = (dev_ssfa[ji]/nstat - dev_sfa[ji]*dev_sfa[ji]) * 10000.0f;
		dev_ssel[ji] = (dev_ssel[ji]/nstat - dev_sel[ji]*dev_sel[ji]) * 10000.0f;
		
		
		if (dev_ssel[ji] >= 8570.0f)
		{
			printf("something wrong at i=%d, j=%d, with nstat=%f, sel=%f, prev=%f\n", i, j, nstat, dev_sel[ji], prev_ssel);
		}

	}
	
}



float * KaspyCycler::getElves()
{
	return g_el;
}

float * KaspyCycler::getSurface()
{
	return g_fsm;
}


void KaspyCycler::findElves()
{
	/// DO CUDA REDUCTION instead of copying back to host mem



	DeviceReduce::Min(d_temp_storage, temp_storage_bytes, g_elf, g_elf_r, F_DATA_SIZE);
	cudaMemcpy(&m_fVars->elfmin, g_elf_r,  sizeof(float), cudaMemcpyDeviceToHost);
	
	DeviceReduce::Max(d_temp_storage, temp_storage_bytes, g_elf, g_elf_r, F_DATA_SIZE);
	cudaMemcpy(&m_fVars->elfmax, g_elf_r,  sizeof(float), cudaMemcpyDeviceToHost);
	
	
}


void KaspyCycler::sendDataToGPU()
{
	//int ewidth = ((int)m_pitch) / sizeof(float);
	int wm1 = m_width - 1 ;
	int hm1 = m_height - 1 ;
	float dte = (float)m_fVars->dte;
	float dte2 = (float)m_fVars->dte * 2.0f;
	float tide_l = (float)m_fVars->tide_l;
	
	float aam2d = m_fArrays->aam2d;
	
	float xmi = m_fVars->xmi;
	float xma = m_fVars->xma;
	float ymi = m_fVars->ymi;
	float yma = m_fVars->yma;
	
	//float xki = m_fWindData->xki;
	//float xka = m_fWindData->xka;
	//float yki = m_fWindData->yki;
	//float yka = m_fWindData->xka;
	
	
	if ( (cudaMemcpyToSymbol(dev_width, &m_width, sizeof(int))  == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_height, &m_height, sizeof(int))  == cudaSuccess)
		&&(cudaMemcpyToSymbol(dev_widthm1, &wm1, sizeof(int))  == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_heightm1, &hm1, sizeof(int))  == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_dte, &dte, sizeof(float))  == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_dte2, &dte2, sizeof(float))  == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_tide_l, &tide_l, sizeof(float))  == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_aam2d, &aam2d, sizeof(float))  == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_xmi, &xmi, sizeof(float))  == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_xma, &xma, sizeof(float))  == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_ymi, &ymi, sizeof(float))  == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_yma, &yma, sizeof(float))  == cudaSuccess)
		//&& (cudaMemcpyToSymbol(dev_xki, &xki, sizeof(float))  == cudaSuccess)
		//&& (cudaMemcpyToSymbol(dev_xka, &xka, sizeof(float))  == cudaSuccess)
		//&& (cudaMemcpyToSymbol(dev_yki, &yki, sizeof(float))  == cudaSuccess)
		//&& (cudaMemcpyToSymbol(dev_yka, &yka, sizeof(float))  == cudaSuccess)
		//&& (cudaMemcpyToSymbol(dev_ewidth, &ewidth,  sizeof(int))  == cudaSuccess)
		)
	{
		printf("GPU constant memory filled\n");
		
		
		//int test_i = 0;
		//int test_f = 0;
		
		
	}
	else
	{
		printf("GPU memory copy error (error code %s)!\n", cudaGetErrorString(cudaGetLastError()));
		
		deinit_device();
		
		exit(-1);
	}
	
	
	size_t s_data_size =  m_height * m_width *  sizeof(float);

	size_t press_data_size =  m_fWindData->ky *  m_fWindData->kx * m_fWindData->kt * sizeof(float);
	size_t uwd_data_size =  m_fWindData->kyu *  m_fWindData->kxu * m_fWindData->ktu * sizeof(float);
	size_t vwd_data_size =  m_fWindData->kyv *  m_fWindData->kxv * m_fWindData->ktv * sizeof(float);
	
	//printf("have m_h ADDR as %llx\n", (unsigned long long)m_h);

	
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
		&& (cudaMemcpy(g_h, &m_fArrays->h[0][0], s_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		
		&& (cudaMemcpy(g_press,&m_press[0],press_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_uwd,&m_uwd[0],uwd_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_vwd,&m_vwd[0],vwd_data_size, cudaMemcpyHostToDevice) == cudaSuccess)
		
		
		
		&& (cudaMemcpy(g_cor, &m_fArrays->cor[0], m_height * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_aru, &m_fArrays->aru[0], m_height * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_arv, &m_fArrays->arv[0],  m_height * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_art, &m_fArrays->art[0],  m_height * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_dx, &m_fArrays->dx[0], m_height * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_dy, &m_fArrays->dy[0], m_height * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess)
		
		&& (cudaMemcpy(g_stations_x, m_stations_x, m_stations * sizeof(int), cudaMemcpyHostToDevice) == cudaSuccess)
		&& (cudaMemcpy(g_stations_y, m_stations_y, m_stations * sizeof(int), cudaMemcpyHostToDevice) == cudaSuccess)
		)
		
	{
		printf("GPU memory filled\n");
	}
	else
	{
		printf("GPU memory copy error (error code %s)!\n", cudaGetErrorString(cudaGetLastError()));

		deinit_device();
		
		exit(-1);
	}
	
	
}

void KaspyCycler::writeStatistics()
{
	size_t s_data_size =  m_height * m_width *  sizeof(float);
	
	dim3 threadsPerSquareBlock(initValues->m_cuda_threads_2d_x, initValues->m_cuda_threads_2d_y);
	
	dim3 numSquareBlocks((m_width + threadsPerSquareBlock.x - 1) / threadsPerSquareBlock.x, (m_height + threadsPerSquareBlock.y - 1) / threadsPerSquareBlock.y);
	
	
	float * host_buf =  (float *) malloc(s_data_size );
	
	float * gpu_buf =  g_sel;
	
	const char * stat_filename = "sel.grd\0                           ";
	
	
	if (host_buf)
	{
		printf("have nstat %d!\n", m_nstat);
		dev_statistics_finalize<<< numSquareBlocks, threadsPerSquareBlock>>>((float)m_nstat);
		
		
		cudaError_t err = cudaMemcpy(host_buf, gpu_buf,  s_data_size, cudaMemcpyDeviceToHost);
		
		if (err == cudaSuccess)
		{
			// CALL WRITEGRD(IM,JM,IM,SSEL,xmi,xma,ymi,yma,NAME)
			
			WRITEGRD(&m_width, &m_height, &m_width, host_buf, &m_fVars->xmi, &m_fVars->xma,&m_fVars->ymi, &m_fVars->yma,stat_filename);

		}
		else
		{
			fprintf(stderr, "Failed to update statistics data  (error code %s)!\n", cudaGetErrorString(err));
		}
		
		
		
		free(host_buf);
	}
	else
	{
		printf("memory allocation failed!\n");
		
		deinit_device();
		
		exit(-1);
	}
	
}

void KaspyCycler::getDataToCPU()
{
	float * h_el =  &m_fArrays->el[0][0];
	

	cudaMemcpyFromSymbol (&g_el, dev_el, sizeof(float *), 0,cudaMemcpyDeviceToHost);
	
	
	cudaError_t err = cudaMemcpy(h_el, g_el,  m_height * m_width * sizeof(float), cudaMemcpyDeviceToHost);
	
	if (err != cudaSuccess)
	{
		fprintf(stderr, "Failed to update host array EL  (error code %s)!\n", cudaGetErrorString(err));
		
		deinit_device();
		
		exit(-1);
	}
	
	
	err = cudaMemcpy(m_station_elves, g_station_elves,  (m_duration-1) * m_stations * sizeof(float), cudaMemcpyDeviceToHost);
	
	if (err != cudaSuccess)
	{
		fprintf(stderr, "Failed to update station data  (error code %s)!\n", cudaGetErrorString(err));
		
		deinit_device();
		
		exit(-1);
	}
}




void KaspyCycler::makeWsurf()
{
	//cudaError_t err;
	
	float hours =  (m_duration - 1 )*m_fVars->dht;
	
	
	int n_iterations = int(hours * 3600.0/m_fVars->dte) + 2;
	
	
	int iold = 0;
	
	int ihour_s = int(600.0 / m_fVars->dte);
	
	float timeh;
	float timeh6;
	float ftim;
	int itimeh;
	int itime6;
	int itime6_old  = 0;
	
    int pressSize = m_fWindData->kx * m_fWindData->ky;
    int windUSize = m_fWindData->kxu * m_fWindData->kyu;
    int windVSize = m_fWindData->kxv * m_fWindData->kyv;
	
	
	int threadsPerBlock = initValues->m_cuda_threads_1d;
	
	int blocksPerGridJ = (m_height + threadsPerBlock - 1) / threadsPerBlock;
	int blocksPerGridI = (m_width + threadsPerBlock - 1) / threadsPerBlock;
	int blocksPerStations = (m_stations + threadsPerBlock - 1) / threadsPerBlock;
	
	dim3 threadsPerSquareBlock(initValues->m_cuda_threads_2d_x, initValues->m_cuda_threads_2d_y);
	
	dim3 numSquareBlocks((m_width + threadsPerSquareBlock.x - 1) / threadsPerSquareBlock.x, (m_height + threadsPerSquareBlock.y - 1) / threadsPerSquareBlock.y);
	
	
	
	for (int i=0; i<n_iterations; i++)
	{
		m_fVars->timeh = i * m_fVars->dte / 3600.0;
		
		m_fVars->timeh6 = (m_fVars->timeh / m_fVars->dht) + 1.0;
		
		timeh = (float)m_fVars->timeh;
		timeh6 = (float)m_fVars->timeh6;
		
		itimeh=int(timeh);
		
		ftim = fmodf(timeh6, 1.0f);
		
		itime6 = (int)timeh6;
		
		
		if ( (i % ihour_s) == 1)
		{
			findElves();
			printf("elves t=%f level=%f,%f \n", timeh, m_fVars->elfmin, m_fVars->elfmax);
		}
		
		
		if (itimeh > iold)
		{
			iold=itimeh;
		
			
			dev_fill_station_data<<< blocksPerStations, threadsPerBlock>>>(itimeh);
			
			/// STATISTICS HERE

			//
			dev_statistics_1<<< numSquareBlocks, threadsPerSquareBlock>>>(ftim);
			
			m_nstat++;
			
			//printf("STTT %d\n", m_nstat);
		}
		
		
		
		
		
		
		
		if (itime6 > itime6_old)
		{
			itime6_old = itime6;
			
			float * p_temp;
			
			p_temp = g_fxb;
			g_fxb = g_fxf;
			g_fxf = p_temp;
			
			p_temp = g_fyb;
			g_fyb = g_fyf;
			g_fyf = p_temp;
			
			p_temp = g_fb;
			g_fb = g_ff;
			g_ff = p_temp;
			
			p_temp = g_fbu;
			g_fbu = g_ffu;
			g_ffu = p_temp;
			
			p_temp = g_fbv;
			g_fbv = g_ffv;
			g_ffv = p_temp;
			
			cudaMemcpyToSymbol(dev_fxf, &g_fxf, sizeof(float *));
			cudaMemcpyToSymbol(dev_fxb, &g_fxb, sizeof(float *));
			cudaMemcpyToSymbol(dev_fyf, &g_fyf, sizeof(float *));
			cudaMemcpyToSymbol(dev_fyb, &g_fyb, sizeof(float *));
			cudaMemcpyToSymbol(dev_ff, &g_ff, sizeof(float *));
			cudaMemcpyToSymbol(dev_fb, &g_fb, sizeof(float *));
			cudaMemcpyToSymbol(dev_ffu, &g_ffu, sizeof(float *));
			cudaMemcpyToSymbol(dev_fbu, &g_fbu, sizeof(float *));
			cudaMemcpyToSymbol(dev_ffv, &g_ffv, sizeof(float *));
			cudaMemcpyToSymbol(dev_fbv, &g_fbv, sizeof(float *));
			
			
			
			g_press0 = g_press + (itime6 - 1) * pressSize;
			getWindPressure('p');
			
			g_uwd0 =  g_uwd + (itime6 - 1) * windUSize;
			getWindPressure('u');
			
			g_vwd0 = g_vwd + (itime6 - 1) * windVSize;
			getWindPressure('v');
			
			//printf("PPP\n");
		}
		
		
		surf_and_flux_1<<<numSquareBlocks, threadsPerSquareBlock>>>(ftim);
		
		/*err = cudaGetLastError();
		
		if (err != cudaSuccess)
		{
			printf("error calling surf_and_flux_1 kernel!  (error code %s)!\n", cudaGetErrorString(err));
			return;
		}*/
		
		elf_and_flux_2<<<numSquareBlocks, threadsPerSquareBlock>>>();
		
		/*if (cudaGetLastError() != cudaSuccess)
		{
			printf("error calling elf_and_flux_2 kernel! \n");
		}*/
		
		/// BCOND 1
		
		
		bcond_1_j<<< blocksPerGridJ, threadsPerBlock>>>();
		
		/*if (cudaGetLastError() != cudaSuccess)
		{
			printf("error calling bcond_1_j kernel! \n");
		}*/
		
		bcond_1_i<<< blocksPerGridI, threadsPerBlock>>>();
		
		/*if (cudaGetLastError() != cudaSuccess)
		{
			printf("error calling bcond_1_i kernel! \n");
		}*/
		
		bcond_1_ji<<< numSquareBlocks, threadsPerSquareBlock>>>();
		
		/*if (cudaGetLastError() != cudaSuccess)
		{
			printf("error calling bcond_1_ji kernel! \n");
		}*/
		
		
		
		if (i % 10 == 9)
		{
			//ADVAVE()
			//       ADVUA=0 ?
			//		FLUXUA=0 ?
			
			
			adv_fluxes_1<<< numSquareBlocks, threadsPerSquareBlock>>>();
			
			/*if (cudaGetLastError() != cudaSuccess)
			{
				printf("error calling adv_fluxes_1 kernel! \n");
			}*/
			
			adv_advua_1<<< numSquareBlocks, threadsPerSquareBlock>>>();
			
			/*if (cudaGetLastError() != cudaSuccess)
			{
				printf("error calling adv_advua_1 kernel! \n");
			}*/
			
			adv_fluxes_2<<< numSquareBlocks, threadsPerSquareBlock>>>();
			
			if (cudaGetLastError() != cudaSuccess)
			{
				printf("error calling adv_fluxes_2 kernel! \n");
			}
			
			adv_advva_2<<< numSquareBlocks, threadsPerSquareBlock>>>();
			
			/*if (cudaGetLastError() != cudaSuccess)
			{
				printf("error calling adv_advva_2 kernel! \n");
			}*/
			
			
		 adv_bot_3<<< numSquareBlocks, threadsPerSquareBlock>>>();
			
			/*if (cudaGetLastError() != cudaSuccess)
			{
				printf("error calling adv_bot_3 kernel! \n");
			}*/
			
			
			// END ADVAVE();
		}
		
		
		uaf_and_vaf_3<<<numSquareBlocks, threadsPerSquareBlock>>>();
		
		/*if (cudaGetLastError() != cudaSuccess)
		{
			printf("error calling uaf_and_vaf_3 kernel! \n");
		}*/
		
		
	 bcond_2_j<<< blocksPerGridJ, threadsPerBlock>>>();
		
		/*if (cudaGetLastError() != cudaSuccess)
		{
			printf("error calling bcond_2_j kernel! \n");
		}*/
		
	 bcond_2_i<<< blocksPerGridI, threadsPerBlock>>>();
		
		/*if (cudaGetLastError() != cudaSuccess)
		{
			printf("error calling bcond_2_i kernel! \n");
		}*/
		
		bcond_2_ji<<< numSquareBlocks, threadsPerSquareBlock>>>();
		
		/*if (cudaGetLastError() != cudaSuccess)
		{
			printf("error calling bcond_2_ji kernel! \n");
		}*/
		
		tps_and_other_arrays_4<<<numSquareBlocks, threadsPerSquareBlock>>>();
		
		/*if (cudaGetLastError() != cudaSuccess)
		{
			printf("error calling tps_and_other_arrays_4 kernel! \n");
		}*/
		
		
		swap_arrays_5<<<1, 1>>>();
		
		/*if (cudaGetLastError() != cudaSuccess)
		{
			printf("error calling swap_arrays_5 kernel! \n");
		}*/

	}
	
	
	getDataToCPU();
	
	
	FILE * hnd = fopen(initValues->m_output_stations, "w");
	
	if (hnd!= NULL)
	{
		for (int i=1; i<m_duration; i++)
		{
			fprintf(hnd, "elves t=%d ", i);
			
			for (int k=0; k<m_stations; k++)
			{
				fprintf(hnd, " %9.3f ", m_station_elves[(i-1) * m_stations + k]  );
			}
						
			fprintf(hnd, " \n");
		}
		
		fclose(hnd);
	}
	
	

	
	

}


void KaspyCycler::getWindPressure(char uv)
{
	int kx, ky, kd, nx, ny;//, nd;

	float xki, xka, yki, yka, xmi, xma, ymi, yma;

	
	float * zero = 0;
	
	if (uv == 'u')
	{
		kx = m_fWindData->kxu;
		ky = m_fWindData->kyu;
		//pk = g_uwd0;
		
		xki = m_fWindData->xkui;
		xka = m_fWindData->xkua;
		yki = m_fWindData->ykui;
		yka = m_fWindData->ykua;
		
		cudaMemcpyToSymbol(dev_p, &g_ffu, sizeof(float *));
		cudaMemcpyToSymbol(dev_pk, &g_uwd0, sizeof(float *));
		cudaMemcpyToSymbol(dev_px, &zero, sizeof(float *));
		cudaMemcpyToSymbol(dev_py, &zero, sizeof(float *));
		
		//p = g_ffu;
	}
	else if (uv == 'v')
	{
		kx = m_fWindData->kxv;
		ky = m_fWindData->kyv;
		//pk = g_vwd0;
		
		xki = m_fWindData->xkvi;
		xka = m_fWindData->xkva;
		yki = m_fWindData->ykvi;
		yka = m_fWindData->ykva;
		
		//p = g_ffv;
		cudaMemcpyToSymbol(dev_p, &g_ffv, sizeof(float *));
		cudaMemcpyToSymbol(dev_pk, &g_vwd0, sizeof(float *));
		cudaMemcpyToSymbol(dev_px, &zero, sizeof(float *));
		cudaMemcpyToSymbol(dev_py, &zero, sizeof(float *));
		
	}
	else if (uv == 'p')
	{
		kx = m_fWindData->kx;
		ky = m_fWindData->ky;
		//float kd = kx;
		//pk = g_press0;
		xki = m_fWindData->xki;
		xka = m_fWindData->xka;
		yki = m_fWindData->yki;
		yka = m_fWindData->yka;

		cudaMemcpyToSymbol(dev_p, &g_ff, sizeof(float *));
		cudaMemcpyToSymbol(dev_pk, &g_press0, sizeof(float *));
		cudaMemcpyToSymbol(dev_px, &g_fxf, sizeof(float *));
		cudaMemcpyToSymbol(dev_py, &g_fyf, sizeof(float *));
		
		
		//p = g_ff;
		//px = g_fxf;
		//py = g_fyf;
	}
	else
	{
		// don't know what to do
		return;
	}
	
	kd = kx;
	
	nx = F_DATA_WIDTH;
	ny = F_DATA_HEIGHT;
	//nd = F_DATA_WIDTH;
	
	xmi = m_fVars->xmi;
	xma = m_fVars->xma;
	ymi = m_fVars->ymi;
	yma = m_fVars->yma;
	
	//float c1=3.1415926/180.0;
	//float c2=111111.0f;
	
	
	float dky=(yka-yki)/(ky-1.0f);
	float dkx=(xka-xki)/(kx-1.0f);
 
	float dy=(yma-ymi)/(ny-1.0f);
	float dx=(xma-xmi)/(nx-1.0f);
	

	
	int threadsPerBlock = initValues->m_cuda_threads_1d;

	
	dim3 threadsPerSquareBlock(initValues->m_cuda_threads_2d_x, initValues->m_cuda_threads_2d_y);
	
	
	dim3 numSquareBlocks((kx  + threadsPerSquareBlock.x ) / threadsPerSquareBlock.x, (ky  + threadsPerSquareBlock.y ) / threadsPerSquareBlock.y);
	
	
	dev_pkk_ij<<<numSquareBlocks, threadsPerSquareBlock>>>(kx, ky, kd);
	
	int blocksPerGridJ = (ky + threadsPerBlock) / threadsPerBlock;
	int blocksPerGridI = (kx + 1 + threadsPerBlock) / threadsPerBlock;
	
	dev_pkk_j<<<threadsPerBlock, blocksPerGridJ>>>(kx, ky);
	dev_pkk_i<<<threadsPerBlock, blocksPerGridI>>>(kx, ky);
	


	
	dim3 numNBlocks(((nx) + threadsPerSquareBlock.x - 1) / threadsPerSquareBlock.x, ((ny) + threadsPerSquareBlock.y - 1) / threadsPerSquareBlock.y);
	
	
	dev_bicubic<<<numNBlocks, threadsPerSquareBlock>>>(kx + 2, ky + 2, 50);

	
	dev_make_p<<<numNBlocks, threadsPerSquareBlock>>>(nx, ny, kx, ky, dx, dy, dkx, dky, xki, yki);
	


}








int KaspyCycler::init_device()
{
	int device_count = 0;
	
	size_t square_size = m_height*m_width * sizeof(float);
	
	//setbuf(stdout,NULL);
	//printf("\n\nIII device is %d\n\n", m_gpu_device);
	//setbuf(stdout,NULL);
	
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
	if ( (cudaMalloc((void **)&g_fbu, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_fbv, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_ffu, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_ffv, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_fb, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_ff, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_fxb, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_fxf, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_fyb, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_fyf, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_wusurf, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_wvsurf, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_dum, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_dvm, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_d, square_size) == cudaSuccess)
		
		&& (cudaMalloc((void **)&g_fluxua, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_fluxva, square_size) == cudaSuccess)
		
		&& (cudaMalloc((void **)&g_ua, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_va, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_uab, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_vab, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_uaf, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_vaf, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_el, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_elb, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_elf, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_elf_r, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_fsm, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_tps, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_advua, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_advva, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_wubot, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_wvbot, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_cbc, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_h, square_size) == cudaSuccess)
		
		&& (cudaMalloc((void **)&g_press,  m_fWindData->ky *  m_fWindData->kx * m_fWindData->kt * sizeof(float)) == cudaSuccess)
		
		&& (cudaMalloc((void **)&g_uwd, m_fWindData->kyu * m_fWindData->kxu * m_fWindData->ktu * sizeof(float)) == cudaSuccess)
		&& (cudaMalloc((void **)&g_vwd, m_fWindData->kyv * m_fWindData->kxv * m_fWindData->ktv * sizeof(float)) == cudaSuccess)

		
		&& (cudaMalloc((void **)&g_cor, m_height * sizeof(float)) == cudaSuccess)
		&& (cudaMalloc((void **)&g_aru, m_height * sizeof(float)) == cudaSuccess)
		&& (cudaMalloc((void **)&g_arv, m_height * sizeof(float)) == cudaSuccess)
		&& (cudaMalloc((void **)&g_art, m_height * sizeof(float)) == cudaSuccess)
		
		&& (cudaMalloc((void **)&g_dx, m_height * sizeof(float)) == cudaSuccess)
		&& (cudaMalloc((void **)&g_dy, m_height * sizeof(float)) == cudaSuccess)
		
		&& (cudaMalloc((void **)&g_stations_x, m_stations * sizeof(int)) == cudaSuccess)
		&& (cudaMalloc((void **)&g_stations_y, m_stations * sizeof(int)) == cudaSuccess)
		&& (cudaMalloc((void **)&g_station_elves, (m_duration -1) *  m_stations * sizeof(float)) == cudaSuccess)
		
		&& (cudaMalloc((void **)&g_sel, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_ssel, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_sfel, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_sfa, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_ssfa, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_sfar, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_ssfar, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_sfelr, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_su, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_sv, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_ssu, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_ssv, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_ssuv, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_ssue, square_size) == cudaSuccess)
		&& (cudaMalloc((void **)&g_ssve, square_size) == cudaSuccess)
		
		
		)
	{
		printf("GPU memory allocated\n");
		
	}
	else
	{
		printf("GPU memory allocation error!\n");
		deinit_device();
		return m_gpu_device;
	}
	
	
	if ( (cudaMemset (g_sel, 0, square_size) == cudaSuccess)
		&& (cudaMemset (g_ssel, 0, square_size) == cudaSuccess)
		&& (cudaMemset (g_sfel, 0, square_size) == cudaSuccess)
		&& (cudaMemset (g_sfa, 0, square_size) == cudaSuccess)
		&& (cudaMemset (g_ssfa, 0, square_size) == cudaSuccess)
		&& (cudaMemset (g_sfar, 0, square_size) == cudaSuccess)
		&& (cudaMemset (g_ssfar, 0, square_size) == cudaSuccess)
		&& (cudaMemset (g_sfelr, 0, square_size) == cudaSuccess)
		&& (cudaMemset (g_su, 0, square_size) == cudaSuccess)
		&& (cudaMemset (g_sv, 0, square_size) == cudaSuccess)
		&& (cudaMemset (g_ssu, 0, square_size) == cudaSuccess)
		&& (cudaMemset (g_ssv, 0, square_size) == cudaSuccess)
		&& (cudaMemset (g_ssuv, 0, square_size) == cudaSuccess)
		&& (cudaMemset (g_ssue, 0, square_size) == cudaSuccess)
		&& (cudaMemset (g_ssve, 0, square_size) == cudaSuccess)
		
		
		)
	{
		printf("GPU memory zeroed\n");
		
	}
	else
	{
		printf("GPU memory fill error!\n");
		deinit_device();
		return m_gpu_device;
	}
	
	
	
	d_temp_storage = 0;
	temp_storage_bytes = 0;
	
	DeviceReduce::Min(d_temp_storage, temp_storage_bytes, g_elf, g_elf_r, F_DATA_SIZE);
	g_allocator.DeviceAllocate(&d_temp_storage, temp_storage_bytes);
	
	
	if ( (cudaMemcpyToSymbol(dev_fbu, &g_fbu, sizeof(g_fbu)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_fbv, &g_fbv, sizeof(g_fbv)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_ffu, &g_ffu, sizeof(g_ffu)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_ffv, &g_ffv, sizeof(g_ffv)) == cudaSuccess)
		
		&& (cudaMemcpyToSymbol(dev_fxb, &g_fxb, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_fxf, &g_fxf, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_fyb, &g_fyb, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_fyf, &g_fyf, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_fb, &g_fb, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_ff, &g_ff, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_wusurf, &g_wusurf, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_wvsurf, &g_wvsurf, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_dum, &g_dum, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_dvm, &g_dvm, sizeof(float *)) == cudaSuccess)

		&& (cudaMemcpyToSymbol(dev_d, &g_d, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_dx, &g_dx, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_dy, &g_dy, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_fluxua, &g_fluxua, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_fluxva, &g_fluxva, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_ua, &g_ua, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_va, &g_va, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_uab, &g_uab, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_vab, &g_vab, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_uaf, &g_uaf, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_vaf, &g_vaf, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_el, &g_el, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_elf, &g_elf, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_elf_r, &g_elf_r, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_elb, &g_elb, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_fsm, &g_fsm, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_tps, &g_tps, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_advua, &g_advua, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_advva, &g_advva, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_aru, &g_aru, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_arv, &g_arv, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_art, &g_art, sizeof(float *)) == cudaSuccess)

		&& (cudaMemcpyToSymbol(dev_wubot, &g_wubot, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_wvbot, &g_wvbot, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_cbc, &g_cbc, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_cor, &g_cor, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_h, &g_h, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_press, &g_press, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_uwd, &g_uwd, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_vwd, &g_vwd, sizeof(float *)) == cudaSuccess)
		
		&& (cudaMemcpyToSymbol(dev_stations_x, &g_stations_x, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_stations_y, &g_stations_y, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_station_elves, &g_station_elves, sizeof(float *)) == cudaSuccess)
		
		&& (cudaMemcpyToSymbol(dev_nstations, &m_stations, sizeof(int)) == cudaSuccess)
		
		&& (cudaMemcpyToSymbol(dev_sel, &g_sel, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_ssel, &g_ssel, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_sfel, &g_sfel, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_sfa, &g_sfa, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_ssfa, &g_ssfa, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_sfar, &g_sfar, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_ssfar, &g_ssfar, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_sfelr, &g_sfelr, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_su, &g_su, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_sv, &g_sv, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_ssu, &g_ssu, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_ssv, &g_ssv, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_ssuv, &g_ssuv, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_ssue, &g_ssue, sizeof(float *)) == cudaSuccess)
		&& (cudaMemcpyToSymbol(dev_ssve, &g_ssve, sizeof(float *)) == cudaSuccess)

		
		
		)
	{
		printf("Device pointers initialized\n");
	}
	else
	{
		printf("Pointer init error!\n");
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
		
		if (g_elf_r)
		{
			cudaFree(g_elf_r);
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
		
		if (g_art)
		{
			cudaFree(g_art);
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
		
		if (g_press)
		{
			cudaFree(g_press);
		}
		
		
		if (g_uwd)
		{
			cudaFree(g_uwd);
		}
		
		if (g_vwd)
		{
			cudaFree(g_vwd);
		}
		
		if (g_sel)
		{
			cudaFree(g_sel);
		}
		
		
		if (g_ssel)
		{
			cudaFree(g_ssel);
		}
		
		if (g_sfel)
		{
			cudaFree(g_sfel);
		}
		
		
		if (g_sfa)
		{
			cudaFree(g_sfa);
		}
		if (g_ssfa)
		{
			cudaFree(g_ssfa);
		}
		if (g_sfar)
		{
			cudaFree(g_sfar);
		}
		if (g_ssfar)
		{
			cudaFree(g_ssfar);
		}
		if (g_sfelr)
		{
			cudaFree(g_sfelr);
		}
		
		if (g_su)
		{
			cudaFree(g_su);
		}
		if (g_sv)
		{
			cudaFree(g_sv);
		}
		if (g_ssu)
		{
			cudaFree(g_ssu);
		}
		if (g_ssv)
		{
			cudaFree(g_ssv);
		}
		if (g_ssuv)
		{
			cudaFree(g_ssuv);
		}
		if (g_ssue)
		{
			cudaFree(g_ssue);
		}
		if (g_ssve)
		{
			cudaFree(g_ssve);
		}
		
		if (d_temp_storage)
		{
			g_allocator.DeviceFree(d_temp_storage);
		}
		
		if (cudaDeviceReset() == cudaSuccess)
		{
			printf("GPU device reset ok\n");
		}
		
		m_gpu_device = -1;
	}
}





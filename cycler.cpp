
#include <stdio.h>
// #include <stdlib.h>

#ifdef _WIN64
    #include <windows.h>
#else
    #include <sys/time.h>
#endif

#include <math.h>


#include "multithreading.h"


#include "fortran_vars.h"

#include "InitValues.h"

#include "KaspyCycler.h"

#include "DrawArrayWindow.h"




double last_measured_time = 0.0;

KaspyCycler * cycler = 0;

DrawArrayWindow  * defaultGlWindow = 0;

InitValues * initValues = 0;


float * c_press = 0;
float * c_uwd = 0;
float * c_vwd = 0;



int * c_stations_x = 0;
int * c_stations_y = 0;

float * c_station_elves = 0;

CUTThread cycler_thread;


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
    
    //ersetb();
    setbuf(stdout,NULL);
    
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

// Intel Fortran WIN naming
extern "C" void CYCLER_READ_INI()
{
	initValues = new InitValues();
	
	initValues->read();
}

// GFortran Unix naming
extern "C" void cycler_read_ini_()
{
	initValues = new InitValues();
	
	initValues->read();
}

// Intel Fortran WIN naming
extern "C" void CYCLER_GET_INT_PARAM(char * parName, int * param)
{
	if (initValues)
	{
		initValues->get_int_parameter(parName, param);
	}
}

// GFortran Unix naming
extern "C" void cycler_get_int_param_(char * parName, int * param)
{
	if (initValues)
	{
		initValues->get_int_parameter(parName, param);
	}
}

// Intel Fortran WIN naming
extern "C" void CYCLER_GET_STRING_PARAM(char * parName, char * param)
{
	if (initValues)
	{
		initValues->get_string_parameter(parName, param);
	}
}

// GFortran Unix naming
extern "C" void cycler_get_string_param_(char * parName, char * param)
{ 
	if (initValues)
	{
		initValues->get_string_parameter(parName, param);
	}
}


void _i_cycler_init(int * icycler, float * vars_marker, double * arrays_marker, double * ffloats_marker,
                int * wind_marker)
{
	cycler_read_ini_();
	
    if ((initValues == 0) || (vars_marker == 0) || (arrays_marker == 0)
        || (ffloats_marker == 0) || (wind_marker == 0))
    {
		*icycler = -1;
		
        return ; 
    } 
	
	
     
    //_i_cycler_time();

	 
	if (initValues->m_show_window)
	{
		defaultGlWindow = new DrawArrayWindow();
	}
	
	
	/// read_grd(char * name, int * nx, int * ny, int * nz,
	//	float * xmi, float * xma, float * ymi, float * yma,
	//float * zmi, float * zma, float ** z = 0, float multiplier)
	
	//ALLOCATE (PRESS(KX,KY,KT))
	//c ,PRESS0(KX,KY))
	//CALL READGR3(KX,KY,KT,XKI,XKA,YKI,YKA,TKI,TKA,namep,PRESS)
	
	//CALL READDIMGR3(KXU,KYU,KTU,nameu)
	//ALLOCATE (UWD(KXU,KYU,KTU))
	//c ,UWD0(KXU,KYU))
	//CALL READGR3
	//1 (KXU,KYU,KTU,XKUI,XKUA,YKUI,YKUA,TKUI,TKUA,nameu,UWD)
	
	
	fortran_wind_data * w_data = (fortran_wind_data *)wind_marker;
	fortran_common_vars * common_vars = (fortran_common_vars *)vars_marker;
	
	
	initValues->read_grd(initValues->m_pressure_grd, &w_data->kx, &w_data->ky, &w_data->kt, &w_data->xki, &w_data->xka, &w_data->yki, &w_data->yka, &w_data->tki, &w_data->tka, &c_press, 0.001f );

	initValues->read_grd(initValues->m_u_wind_grd, &w_data->kxu, &w_data->kyu, &w_data->ktu, &w_data->xkui, &w_data->xkua, &w_data->ykui, &w_data->ykua, &w_data->tkui, &w_data->tkua, &c_uwd );
	
	initValues->read_grd(initValues->m_v_wind_grd, &w_data->kxv, &w_data->kyv, &w_data->ktv, &w_data->xkvi, &w_data->xkva, &w_data->ykvi, &w_data->ykva, &w_data->tkvi, &w_data->tkva, &c_vwd );
	
	
	int nstations = 0;
	
	initValues->read_stations(&nstations, &c_stations_x, &c_stations_y, &c_station_elves);
	
	printf("read %d stations\n", nstations);
	
	//int size = w_data->kxv * w_data->kyv * w_data->ktv;
	
	//initValues->save_z("c_vwd.ttt", c_vwd, size, w_data->kxv);
	

	///
	
	common_vars->dht = (w_data->tka - w_data->tki) / (float)(w_data->kt - 1);
	
	//printf("tka %f tki %f kt %d dht %f\n",w_data->tka, w_data->tki, w_data->kt,  common_vars->dht );
	
	
	
	_i_cycler_time();
	
    cycler = new KaspyCycler(common_vars, (fortran_common_arrays *)arrays_marker,
							 (fortran_ffloats *)ffloats_marker, w_data,
                             c_press, c_uwd, c_vwd, nstations, c_stations_x, c_stations_y, c_station_elves, initValues->m_duration);
	
	
	*icycler = 0;
     
    
    return;
}


// Intel Fortran WIN naming
extern "C" void CYCLER_CREATE(int * icycler, float * vars_marker, double * arrays_marker, double * ffloats_marker,
                             int * wind_marker)
{
    _i_cycler_init(icycler, vars_marker, arrays_marker, ffloats_marker, wind_marker);
}

// GFortran Unix naming
extern "C" void cycler_create_(int * icycler, float * vars_marker, double * arrays_marker, double * ffloats_marker,
                              int * wind_marker)
{
    _i_cycler_init(icycler, vars_marker, arrays_marker, ffloats_marker, wind_marker);
}



// GFortran Unix naming
extern "C" void cycler_load_(int * icycler)
{
    if (cycler)
    {
		int device = cycler->init_device();
		
		if (device >= 0)
		{
			setbuf(stdout,NULL);
			printf("before init GL\n");
			
			if (initValues && initValues->m_show_window)
			{
				if (defaultGlWindow->gl_init(device) < 0)
				{
					printf("unable to init GL window!\n");
					*icycler = -1;
					
					return;
				}
			}

				
			cycler->sendDataToGPU();
			
			if (initValues && initValues->m_show_window)
			{
				defaultGlWindow->set_data_to_display(cycler->getElves(), cycler->getSurface(), F_DATA_WIDTH, F_DATA_HEIGHT, F_DATA_WIDTH);
			}

		}
		else
		{
			printf("unable to init CUDA device!\n");
			
			*icycler = -1;
		}
		
    }
	else
	{
		*icycler = -1;
	}
	
	return;
}

// Intel Fortran WIN naming
extern "C" void CYCLER_LOAD(int * icycler)
{
	if (cycler)
	{
		int device = cycler->init_device();
		
		if (device >= 0)
		{
			setbuf(stdout,NULL);
			printf("before init GL\n");
			
			if (initValues && initValues->m_show_window)
			{
				if (defaultGlWindow->gl_init(device) < 0)
				{
					printf("unable to init GL window!\n");
					*icycler = -1;
					
					return;
				}
			}
			
			
			cycler->sendDataToGPU();
			
			if (initValues && initValues->m_show_window)
			{
				defaultGlWindow->set_data_to_display(cycler->getElves(), cycler->getSurface(), F_DATA_WIDTH, F_DATA_HEIGHT, F_DATA_WIDTH);
			}
			
		}
		else
		{
			printf("unable to init CUDA device!\n");
			
			*icycler = -1;
		}
		
	}
	else
	{
		*icycler = -1;
	}
	
	return;
}




// GFortran Unix naming
extern "C" void cycler_get_data_back_(int * icycler)
{
    if (cycler)
    {
        cycler->getDataToCPU();
    }
}

// Intel Fortran WIN naming
extern "C" void CYCLER_GET_DATA_BACK(int * icycler)
{
    if (cycler)
    {
        cycler->getDataToCPU();
    }
}


CUT_THREADROUTINE cycler_work(void * data)
{
	if (cycler)
	{
		cycler->makeWsurf();
	}
}





// Intel Fortran WIN naming
extern "C" void CYCLER_WSURF(int * icycler)
{
	cycler_thread = cutStartThread(cycler_work, icycler);
	
	cutEndThread(cycler_thread);
	
	
	cutDestroyThread(cycler_thread);
}

// GFortran Unix naming
extern "C" void cycler_wsurf_(int * icycler)
{
	CYCLER_WSURF(icycler);
}


// GFortran Unix naming
extern "C" void cycler_find_elves_(int * icycler)
{
	if (cycler)
	{
		cycler->findElves();
	}
}

// Intel Fortran WIN naming
extern "C" void CYCLER_FIND_ELVES(int * icycler)
{
	if (cycler)
	{ 
		cycler->findElves();
	}
}




void _i_cycler_destroy(int * icycler)
{
    _i_cycler_time();
	
	if (c_press)
	{
		free(c_press);
	}
	
	if (c_uwd)
	{
		free(c_uwd);
	}
	
	if (c_vwd)
	{
		free(c_vwd);
	}
	
	if (c_stations_x)
	{
		free(c_stations_x);
	}
	
	if (c_stations_y)
	{
		free(c_stations_y);
	}
	
	if (c_station_elves)
	{
		free(c_station_elves);
	}
	
	if (defaultGlWindow)
	{
		defaultGlWindow->gl_deinit();
		
		delete defaultGlWindow;
		
		defaultGlWindow = 0;
	}
	
	if (cycler)
    {
		cycler->deinit_device();
		
        delete cycler;
    }
	
	if (initValues)
	{
		delete initValues;
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


 
/*extern "C" void READDIMGR3(int * nx,int * ny,int * nz,char * name)
{ 
	if (initValues)
	{
		initValues->read_grd(name, nx,ny,nz);
	}
}

extern "C" void readdimgr3_(int * nx,int * ny,int * nz,char * name)
{
	if (initValues)
	{
		initValues->read_grd(name, nx,ny,nz);
	}
}

extern "C" void READGR3(int * nx,int * ny,int * nz,float * xmi,float * xma,float * ymi,float * yma,float * zmi,float * zma,
						char * name,float * z)
{  
	if (initValues)
	{
		initValues->read_grd(name, nx,ny,nz,xmi,xma,ymi,yma,zmi,zma,z);
	} 
}
 
extern "C" void readgr3_(int * nx,int * ny,int * nz,float * xmi,float * xma,float * ymi,float * yma,float * zmi,float * zma,
						 char * name,float * z)
{
	//printf("here\n"); 
	//printf("z is %f \n", z[0]); 
	  
	if (initValues) 
	{
		initValues->read_grd(name, nx,ny,nz,xmi,xma,ymi,yma,zmi,zma,z);
	}
}*/


extern "C" void SAVE_Z(int * nx,int * nsize,float * z, const char * name)
{
	//printf("here\n");
	//printf("z is %f \n", z[0]);
	
	if (initValues)
	{
		initValues->save_z(name, z, *nsize, *nx);
	}
}



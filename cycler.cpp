
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

fortran_wind_data * w_data = 0;
fortran_common_vars * common_vars = 0;
fortran_common_arrays * common_arrays = 0;

fortran_ffloats * common_floats = 0;



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



void _i_cycler_init(int * icycler)
{
	initValues = new InitValues();
	
	initValues->read();
	
	if (initValues->m_show_window)
	{
		defaultGlWindow = new DrawArrayWindow();
	}
	

	
	w_data = (fortran_wind_data *) malloc(sizeof(fortran_wind_data) );
	
	common_vars = (fortran_common_vars *) malloc(sizeof(fortran_common_vars) );
	
	common_arrays = (fortran_common_arrays *) malloc(sizeof(fortran_common_arrays) );
	
	common_floats = (fortran_ffloats *) malloc(sizeof(fortran_ffloats) );
	
	
	initValues->read_grd(initValues->m_pressure_grd, &w_data->kx, &w_data->ky, &w_data->kt, &w_data->xki, &w_data->xka, &w_data->yki, &w_data->yka, &w_data->tki, &w_data->tka, &c_press, 0.001f );

	initValues->read_grd(initValues->m_u_wind_grd, &w_data->kxu, &w_data->kyu, &w_data->ktu, &w_data->xkui, &w_data->xkua, &w_data->ykui, &w_data->ykua, &w_data->tkui, &w_data->tkua, &c_uwd );
	
	initValues->read_grd(initValues->m_v_wind_grd, &w_data->kxv, &w_data->kyv, &w_data->ktv, &w_data->xkvi, &w_data->xkva, &w_data->ykvi, &w_data->ykva, &w_data->tkvi, &w_data->tkva, &c_vwd );
	
	
	int nstations = 0;
	
	initValues->read_stations(&nstations, &c_stations_x, &c_stations_y, &c_station_elves);
	
	printf("read %d stations\n", nstations);
	

	common_vars->dht = (w_data->tka - w_data->tki) / (float)(w_data->kt - 1);
	

	
	_i_cycler_time();
	
    cycler = new KaspyCycler(common_vars, common_arrays,
							 common_floats, w_data,
                             c_press, c_uwd, c_vwd, nstations, c_stations_x, c_stations_y, c_station_elves, initValues->m_duration);
	
	
	*icycler = 0;
     
    
    return;
}


void _i_cycler_load(int * icycler)
{
	if (cycler && initValues)
	{
		initValues->tide_from_grd(common_arrays, common_vars, initValues->m_initial_tide_grd);
		
		int device = cycler->init_device();
		
		if (device >= 0)
		{
			setbuf(stdout,NULL);
			printf("before init GL\n");
			
			if (initValues && initValues->m_show_window)
			{
				if (defaultGlWindow->gl_init(device, initValues->m_frames_per_second) < 0)
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





CUT_THREADPROC cycler_work(void * data)
{
	if (cycler)
	{
		cycler->makeWsurf();
		
		cycler->writeStatistics("ssel");
		
		cycler->writeStatistics("ssfar");
		cycler->writeStatistics("sfelr");
		
		cycler->writeStatistics("ssu");
		cycler->writeStatistics("ssv");
		cycler->writeStatistics("ssuv");
		cycler->writeStatistics("ssue");
		cycler->writeStatistics("ssve");
		
	}
	
	CUT_THREADEND;
}





void _i_cycler_wsurf(int * icycler)
{
	cycler_thread = cutStartThread((CUT_THREADROUTINE)cycler_work, icycler);
	
	
	if (initValues->m_show_window)
	{
		defaultGlWindow->gl_show();
	}
	
	cutEndThread(cycler_thread);
	
	
	cutDestroyThread(cycler_thread);
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




int main(int argc, const char * argv[])
{

	
	int icycler = -1;
	
	_i_cycler_init(&icycler);
	
	if (icycler < 0)
	{
		printf("cycler creation failed!\n");
		exit(-1);
	}
	
	_i_cycler_load(&icycler);
	
	if (icycler < 0)
	{
		printf("data loading failed!\n");
		exit(-1);
	}
	
	_i_cycler_wsurf(&icycler);
	
	_i_cycler_destroy(&icycler);
}



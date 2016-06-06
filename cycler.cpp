
#include <stdio.h>
// #include <stdlib.h>

#ifdef _WIN64
    #include <windows.h>
#else
    #include <sys/time.h>
#endif

#include <math.h>

#include "fortran_vars.h"

#include "KaspyCycler.h"

#include "DrawArrayWindow.h"




double last_measured_time = 0.0;

KaspyCycler * cycler = 0;

DrawArrayWindow  * defaultGlWindow = 0;


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



int _i_cycler_init(float * vars_marker, double * arrays_marker, double * ffloats_marker,
                int * wind_marker,
                   float * press,  float * uwd, float * vwd)
{
    if ((vars_marker == 0) || (arrays_marker == 0)
        || (ffloats_marker == 0) || (press == 0)
        || (wind_marker == 0))
    {
        return -1;
    }
    
    
    _i_cycler_time();
	
	defaultGlWindow = new DrawArrayWindow();
	
    cycler = new KaspyCycler((fortran_common_vars *)vars_marker,
                             (fortran_common_arrays *)arrays_marker, (fortran_ffloats *)ffloats_marker,
                             (fortran_wind_data *)wind_marker,
                             press, uwd, vwd);
     
    
    return 0;
}


// Intel Fortran WIN naming
extern "C" int CYCLER_CREATE(float * vars_marker, double * arrays_marker, double * ffloats_marker,
                             int * wind_marker,
                             float * press, float * uwd, float * vwd)
{
    return _i_cycler_init(vars_marker, arrays_marker, ffloats_marker, wind_marker,
                          press, uwd, vwd);
}

// GFortran Unix naming
extern "C" int cycler_create_(float * vars_marker, double * arrays_marker, double * ffloats_marker,
                              int * wind_marker,
                              float * press,  float * uwd,  float * vwd)
{
    return _i_cycler_init(vars_marker, arrays_marker, ffloats_marker, wind_marker,
                          press, uwd,  vwd);
}



// GFortran Unix naming
extern "C" int cycler_load_(int * icycler)
{
    if (cycler)
    {
		int device = cycler->init_device();
		
		if (device >= 0)
		{
			if (defaultGlWindow->gl_init(device) >= 0)
			{
				
				cycler->sendDataToGPU();
				
				defaultGlWindow->set_data_to_display(cycler->getElves(), cycler->getSurface(), F_DATA_WIDTH, F_DATA_HEIGHT, F_DATA_WIDTH);
				
				return 1;
			}
			else
			{
				printf("unable to init GL window!\n");
				return -1;
			}

		}
		else
		{
			printf("unable to init CUDA device!\n");
			return -1;
		}
		
    }
	
	return -1;
}

// Intel Fortran WIN naming
extern "C" int CYCLER_LOAD(int * icycler)
{
    if (cycler)
    {
		if (cycler->init_device() >= 0)
		{
			cycler->sendDataToGPU();
			return 1;
		}
		else
		{
			printf("unable to init CUDA device!\n");
			return -1;
		}
    }
	
	return -1;
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



// GFortran Unix naming
extern "C" void cycler_wsurf_(int * icycler)
{
    if (cycler)
    {
        cycler->makeWsurf();
    }
}

// Intel Fortran WIN naming
extern "C" void CYCLER_WSURF(int * icycler)
{
    if (cycler)
    {
        cycler->makeWsurf();
    }
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



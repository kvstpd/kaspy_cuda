
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



void _i_cycler_init(int * icycler, float * vars_marker, double * arrays_marker, double * ffloats_marker,
                int * wind_marker,
                   float * press,  float * uwd, float * vwd)
{
    if ((vars_marker == 0) || (arrays_marker == 0)
        || (ffloats_marker == 0) || (press == 0)
        || (wind_marker == 0))
    {
		*icycler = -1;
		
        return ;
    }
    
    
    _i_cycler_time();
	
	defaultGlWindow = new DrawArrayWindow();
	
    cycler = new KaspyCycler((fortran_common_vars *)vars_marker,
                             (fortran_common_arrays *)arrays_marker, (fortran_ffloats *)ffloats_marker,
                             (fortran_wind_data *)wind_marker,
                             press, uwd, vwd);
	
	
	*icycler = 0;
     
    
    return;
}


// Intel Fortran WIN naming
extern "C" void CYCLER_CREATE(int * icycler, float * vars_marker, double * arrays_marker, double * ffloats_marker,
                             int * wind_marker,
                             float * press, float * uwd, float * vwd)
{
    _i_cycler_init(icycler, vars_marker, arrays_marker, ffloats_marker, wind_marker,
                          press, uwd, vwd);
}

// GFortran Unix naming
extern "C" void cycler_create_(int * icycler, float * vars_marker, double * arrays_marker, double * ffloats_marker,
                              int * wind_marker,
                              float * press,  float * uwd,  float * vwd)
{
    _i_cycler_init(icycler, vars_marker, arrays_marker, ffloats_marker, wind_marker,
                          press, uwd,  vwd);
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
			if (defaultGlWindow->gl_init(device) >= 0)
			{
				
				cycler->sendDataToGPU();
				
				defaultGlWindow->set_data_to_display(cycler->getElves(), cycler->getSurface(), F_DATA_WIDTH, F_DATA_HEIGHT, F_DATA_WIDTH);
				
				//*icycler = 1;
			}
			else
			{
				printf("unable to init GL window!\n");
				*icycler = -1;
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
			if (defaultGlWindow->gl_init(device) >= 0)
			{
				
				cycler->sendDataToGPU();
				
				defaultGlWindow->set_data_to_display(cycler->getElves(), cycler->getSurface(), F_DATA_WIDTH, F_DATA_HEIGHT, F_DATA_WIDTH);
				
				*icycler = 1;
			}
			else
			{
				printf("unable to init GL window!\n");
				*icycler = -1;
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
	
	if (defaultGlWindow && (cycler->m_fVars->iint % 50 == 0 ) )
	{
		defaultGlWindow->gl_draw_frame();
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



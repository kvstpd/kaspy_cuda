//
//  KaspyCycler.hpp
//  kaspy_cuda
//
//  Created by Andrei Koulikov on 24.05.16.
//
//

#ifndef KaspyCycler_hpp
#define KaspyCycler_hpp

#include <stdlib.h>
#include <stdio.h>

#include "fortran_vars.h"


class KaspyCycler
{
public:
    
    KaspyCycler(fortran_common_vars * for_vars, fortran_common_arrays * for_arrays,
                fortran_ffloats * for_floats, fortran_wind_data * for_wind_data,
                float * for_press, float * for_uwd, float * for_vwd, int n_stations, int * c_stations_x, int * c_stations_y, float * c_station_elves, int n_duration) :
    m_width(F_DATA_WIDTH),
    m_height(F_DATA_HEIGHT),
    m_fVars(for_vars),
    m_fArrays(for_arrays),
    m_fFloats(for_floats),
    m_fWindData(for_wind_data),
    m_press(for_press),
//    m_press0(for_press0),
    m_uwd(for_uwd),
//    m_uwd0(for_uwd0),
    m_vwd(for_vwd),
	m_stations(n_stations),
	m_stations_x(c_stations_x),
	m_stations_y(c_stations_y),
	m_station_elves(c_station_elves),
	m_duration(n_duration),
	m_nstat(0),
//    m_vwd0(for_vwd0),
	m_gpu_device(-1),
	d_temp_storage(0),
	temp_storage_bytes(0)
    {
        setbuf(stdout,NULL);
        
        printf("cycler alloc with dht %f, array marker %32.32f \n", m_fVars->dht, m_fArrays->marker);
        printf("end marker %32.32f \n", m_fArrays->end_marker);
        printf("ff end %32.32f \n", m_fFloats->end_marker);
        printf("wind data: \n \t %d \t %d \t %d \n \t %d \t %d \t %d \n \t %d \t %d \t %d \n",
               m_fWindData->kx, m_fWindData->ky, m_fWindData->kt,
               m_fWindData->kxu, m_fWindData->kyu, m_fWindData->ktu,
               m_fWindData->kxv, m_fWindData->kyv, m_fWindData->ktv);
		
		
		//setbuf(stdout,NULL);
		//printf("\n\nCCC device is %d\n\n", m_gpu_device);
		//setbuf(stdout,NULL);
    }
    
 
	int init_device();
	void deinit_device();
	
    void sendDataToGPU();
    void getDataToCPU();
	
	
	void writeStatistics();


	void getWindPressure(char uv);
	void getWindPressureG(char uv);

    void findElves();
	
	float * getElves();
	float * getSurface();
	
    void makeWsurf();
    
    
    
    int m_width;
    int m_height;
 
	int m_stations;
	int m_duration;
	
	int m_nstat;

    fortran_common_vars * m_fVars;
    fortran_common_arrays * m_fArrays;
    fortran_ffloats * m_fFloats;
    fortran_wind_data * m_fWindData;
	
	float * m_h;
	
    float * m_press;
//    float * m_press0;
    float * m_uwd;
//    float * m_uwd0;
    float * m_vwd;
//    float * m_vwd0;
	int * m_stations_x;
	int * m_stations_y;
	float * m_station_elves;
	
private:
    //float ftim;
    //float btim;

    
	int m_gpu_device;

	size_t m_pitch;
	
	size_t m_press_pitch;
	size_t m_wu_pitch;
	size_t m_wv_pitch;
	
	void            * d_temp_storage;
	size_t          temp_storage_bytes;
};





#endif /* KaspyCycler_hpp */

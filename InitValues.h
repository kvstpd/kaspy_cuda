

#ifndef InitValues_h
#define InitValues_h




class InitValues
{
public:
	InitValues(char * init_file = 0) :
	m_iterations(4),
	m_cuda_threads_1d(256),
	m_cuda_threads_2d_x(8),
	m_cuda_threads_2d_y(8),
	m_iters_per_frame(10),
	m_show_window(false)
	{
		if (init_file == 0)
		{
			m_file_name = "Kaspy.ini";
		}
		else
		{
			m_file_name = init_file;
		}
	}
	
	void read();
	
	void scan_parameter(char * parName, char * parValue);
	
	void get_int_parameter(char * parName, int * param);
	
	
	void get_string_parameter(char * parName, char * param);
	
	void read_grd(char * name, int * nx, int * ny, int * nz,
							 float * xmi, float * xma, float * ymi, float * yma,
				  float * zmi, float * zma, float * z);

	
	void save_z(const char * name, float * z, int zSize, int nCols);

	
	unsigned int m_iterations;
	unsigned int m_cuda_threads_1d;
	unsigned int m_cuda_threads_2d_x;
	unsigned int m_cuda_threads_2d_y;
	unsigned int m_show_window;
	unsigned int m_iters_per_frame;
	
	char m_pressure_grd[128];
	char m_u_wind_grd[128];
	char m_v_wind_grd[128];
	char m_initial_tide_grd[128];
	char m_stations_file[128];
	
	
private:
	const char * m_file_name;
	
};


#endif /* InitValues_h */

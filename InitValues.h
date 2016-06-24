

#ifndef InitValues_h
#define InitValues_h


class InitValues
{
public:
	InitValues(char * init_file)
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
	
	void read()
	{
		
	}
	
	
	unsigned int m_iterations;
	unsigned int m_cuda_threads_1d;
	unsigned int m_cuda_threads_2d_x;
	unsigned int m_cuda_threads_2d_y;
	
private:
	const char * m_file_name;
	
};


#endif /* InitValues_h */

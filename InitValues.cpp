//
//  InitValues.cpp
//  kaspy_cuda
//
//  Created by pprd on 29.06.16.
//
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "InitValues.h"





void InitValues::read()
{
	FILE * hnd = fopen(m_file_name, "r");
	
	if (hnd!= NULL)
	{
		char lineChars[128];
		
		int iValue;
		
		char * line = &lineChars[0];
		
		char * parName;
		char * parValue;
		char * equals;
		
		while ((line = fgets(line, sizeof(lineChars), hnd)))
		{
			switch (line[0])
			{
				case ';':
				case '[':
				case ' ':
				case '\n':
					continue;
					break;
					
				
				default:
					equals = strchr(line, '=');
					
					if (equals != NULL)
					{
						equals[0] = 0;
						parName = line;
						parValue = equals + 1;
						
						scan_parameter(parName,parValue);
					}
					
					break;
			};
			//
		}
		
		
		fclose(hnd);
	}
	
	
}

void InitValues::scan_parameter(char * parName, char * parValue)
{
	if (strcmp("iterations", parName) == 0)
	{
		sscanf(parValue, "%d", &m_iterations);
		
		printf("have to iterate %d times\n", m_iterations);
	}
	else if (strcmp("show_window", parName) == 0)
	{
		sscanf(parValue, "%d", &m_show_window);
		
		printf("show_window: %d\n", m_show_window);
	}
	else if (strcmp("cuda_threads_1d", parName) == 0)
	{
		sscanf(parValue, "%d", &m_cuda_threads_1d);
		
		printf("Will have %d 1-D CUDA threads\n", m_cuda_threads_1d);
	}
	else if (strcmp("cuda_threads_2d_x", parName) == 0)
	{
		sscanf(parValue, "%d", &m_cuda_threads_2d_x);
		
		printf("Will have %d 2-D CUDA threads on X axis\n", m_cuda_threads_2d_x);
	}
	else if (strcmp("cuda_threads_2d_y", parName) == 0)
	{
		sscanf(parValue, "%d", &m_cuda_threads_2d_y);
		
		printf("Will have %d 2-D CUDA threads on Y axis\n", m_cuda_threads_2d_y);
	}
	else if (strcmp("iters_per_frame", parName) == 0)
	{
		sscanf(parValue, "%d", &m_iters_per_frame);
		
		printf("Will make %d iterations before drawing new frame\n",m_iters_per_frame);
	}
	else if (strcmp("pressure_grd", parName) == 0)
	{
		sscanf(parValue, "%s", m_pressure_grd);
		
		printf("Pressure GRD file is %s\n", m_pressure_grd);
	}
	else if (strcmp("u_wind_grd", parName) == 0)
	{
		sscanf(parValue, "%s", m_u_wind_grd);
		
		printf("U wind GRD file is %s\n", m_u_wind_grd);
	}
	else if (strcmp("v_wind_grd", parName) == 0)
	{
		sscanf(parValue, "%s", m_v_wind_grd);
		
		printf("V wind GRD file is %s\n", m_v_wind_grd);
	}
	
	
}

void InitValues::get_int_parameter(char * parName, int * param)
{
	if (strcmp("iterations", parName) == 0)
	{
		*param = m_iterations;
	}
	else if (strcmp("show_window", parName) == 0)
	{
		*param = m_show_window;
	}
	else if (strcmp("cuda_threads_1d", parName) == 0)
	{
		*param = m_cuda_threads_1d;
	}
	else if (strcmp("cuda_threads_2d_x", parName) == 0)
	{
		*param = m_cuda_threads_2d_x;
	}
	else if (strcmp("cuda_threads_2d_y", parName) == 0)
	{
		*param = m_cuda_threads_2d_y;
	}
	else if (strcmp("iters_per_frame", parName) == 0)
	{
		*param = m_iters_per_frame;
	}
}


void InitValues::get_string_parameter(char * parName, char * param)
{
	if (strcmp("pressure_grd", parName) == 0)
	{
		strncpy(param, m_pressure_grd, strlen(m_pressure_grd)+1 );
	}
}


void InitValues::read_grd(char * name, int * nx, int * ny, int * nz,
			  float * xmi, float * xma, float * ymi, float * yma,
			  float * zmi, float * zma, float ** z, float multiplier)
{
	char lineChars[1024];
	char * line = &lineChars[0];
	
	int readValues = 0;
	

	
	FILE * hnd = fopen(name, "r");
	
	if (hnd!= NULL)
	{
		fgets(line, sizeof(lineChars), hnd);
		
		if ((line = fgets(line, sizeof(lineChars), hnd)))
		{
			sscanf(line, "%d %d %d", nx, ny, nz);
		}
		else
		{
			printf("GR3 file read error!\n");
			return;
		}
		
		
		if (z && xmi && xma && ymi && yma && zmi && zma)
		{
			if ((line = fgets(line, sizeof(lineChars), hnd)))
			{
				sscanf(line, "%f %f", xmi, xma);
			}
			
			if ((line = fgets(line, sizeof(lineChars), hnd)))
			{
				sscanf(line, "%f %f", ymi, yma);
			}
			
			if ((line = fgets(line, sizeof(lineChars), hnd)))
			{
				sscanf(line, "%f %f", zmi, zma);
			}
			
			
			int maxValues =  (*nx) * (*ny) * (*nz);
			
			char * val;
			
			bool mustEnd = false;
			
			int nchars;
			
			float * zz = (float *) malloc(maxValues * sizeof(float) );
			
			
	
			
			if (!zz)
			{
				printf("memory allocation error!\n");
				return;
			}
			
			float * zz0 = zz;
			
			*z = zz0;
			
			
			while (!mustEnd && ((line = fgets(line, sizeof(lineChars), hnd))))
			{
				val = line;
				
				while ((val = strchr(val, ' ')))
				{
					if (++readValues > maxValues)
					{
						//printf("reached ix %d iy %d iz %d\n", ix, iy, iz);
						mustEnd = true;
						break;
					}
					
					
					sscanf(val, "%f%n", zz++, &nchars );
					val += nchars;
				}
				
			}
			
			//printf("read %d values\n", readValues - 1);
			
			
			
			if (multiplier != 1.0f)
			{
				while (zz0 < zz)
				{
					*zz0 *= multiplier;
					zz0++;
				}
			}	
		}
		
		
		fclose(hnd);
		

	}
}


void InitValues::save_z(const char * name, float * z, int zSize, int nCols)
{
	FILE * hnd = fopen(name, "w");
	
	if (hnd!= NULL)
	{
		for (int i = 0; i < zSize; i++)
		{
			if (i % nCols == 0)
			{
				fprintf(hnd, "\n");
			}
			
			fprintf(hnd, " %9.1f", z[i]);
		}
		
		fclose(hnd);
	}
}

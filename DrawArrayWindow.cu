
#include <GL/glew.h>

#include <GL/freeglut.h>


// CUDA utilities and system includes
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

// Includes
#include <stdlib.h>
#include <stdio.h>



#include "DrawArrayWindow.h"


void display(void);
void keyboard(unsigned char key, int /*x*/, int /*y*/);
void reshape(int x, int y);
void cleanup(void);


static DrawArrayWindow * currentWindow = 0;



__global__ void gpu_gen_buffer(int width, int height, unsigned char * buf, float kw, float kh, int zw, const float * z, const float * h, float brighten)//,  float m_b2y)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;

	int zi  = int( kw * i);
	int zj  =  int( kh * j);

	int ind;
	float zval, hval;

	if ((i < width) && (j < height))
	{
		ind = (i + (height - j - 1) * width) * 4;

		zval = z[zi + zw * zj] * brighten * 255.0f;

		hval = h[zi + zw * zj];

		if (hval <= 0.f)
		{
			if (hval <= -2500.f)
			{
				buf[ind] =  140;//int(fminf(zval, 255.0f)); //r
				buf[ind + 1] = 100; //g
				buf[ind + 2] = 0;
			}
			else if (hval <= -1000.f)
			{
				buf[ind] =  150;//int(fminf(zval, 255.0f)); //r
				buf[ind + 1] = 130; //g
				buf[ind + 2] = 0;
			}
			else if (hval <= -250.f)
			{
				buf[ind] =  40;//int(fminf(zval, 255.0f)); //r
				buf[ind + 1] = 120; //g
				buf[ind + 2] = 30;
			}
			else
			{
				buf[ind] =  20;//int(fminf(zval, 255.0f)); //r
				buf[ind + 1] = 100; //g
				buf[ind + 2] = 20;
			}

		}
		else if (zval >=0.0f)
		{
			int zvv = int(fminf(zval, 255.0f));

			buf[ind] =  0xff;//int(fminf(zval, 255.0f)); //r
			buf[ind + 1] = 0xff - zvv; //g
			buf[ind + 2] = 0xff - zvv; //b
		}
		else
		{
			int zvv = int(fminf(-zval, 255.0f));

			buf[ind] = 0xff - zvv;  //r
			buf[ind + 1] = 0xff - zvv; //g
			buf[ind + 2] = 0xff; //b
		}



		//buf[ind] += i % 0xff; //r
		//buf[ind + 1] = 0x00; //g
		//buf[ind + 2] += ind % 0xff; //b
		buf[ind + 3] = 0xff; //a
	}
	

}



DrawArrayWindow::DrawArrayWindow()
{
	m_cuda_gResource = 0;
	//m_texArray = 0;
	m_pixData = 0;

	m_gl_texId = 0;
	m_gl_bufId = 0;

	m_width = DAW_START_WINDOW_WIDTH;
	m_height = DAW_START_WINDOW_HEIGHT;


	m_gpu_original = 0;
	m_original_width  = 0;
	m_original_height = 0;

	m_brighten = 1.0f;

}

void DrawArrayWindow::set_data_to_display(float * gpu_zdata, float * gpu_hdata, int w, int h, int p_w)
{
	m_gpu_original = gpu_zdata;
	m_gpu_h = gpu_hdata;

	m_original_width  = w;
	m_original_height = h;

	m_original_pitched_width = p_w;
}


int DrawArrayWindow::gl_init(int device)
{
	int argc = 0;

	setbuf(stdout,NULL);
    printf("starting to init GL\n");

    glutInit(&argc, 0);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowSize(m_width, m_height);
    glutCreateWindow("Ocean");

    glewInit();

    if (!glewIsSupported("GL_VERSION_1_5 GL_ARB_vertex_buffer_object GL_ARB_pixel_buffer_object"))
    {
        fprintf(stderr, "Error: no GL extensions found!\n");

		return -1;
    }
	
	cudaError_t err = cudaGetLastError();

	cudaGLSetGLDevice(device);

	err = cudaGetLastError();


    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutReshapeFunc(reshape);


	 glutCloseFunc(cleanup);


    glDisable(GL_DEPTH_TEST);
	glEnable(GL_TEXTURE_2D);  


	 if (this->gl_rebuild_texture(m_width, m_height) < 0 )
	 {
        fprintf(stderr, "Error: texture creation error!\n");

		return -1;
	 }
	 
 

	currentWindow = this;

	setbuf(stdout,NULL);
    printf("GL init finished\n");
    
	return 0;
}


void DrawArrayWindow::clear_gl_objects()
{
	glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);

		

	if (m_cuda_gResource)
	{
		cudaGraphicsUnregisterResource(m_cuda_gResource);

		m_cuda_gResource = 0;
	}


	if (m_gl_texId)
	{
		glDeleteTextures(1, &m_gl_texId);
		m_gl_texId = 0;
	}


	if (m_gl_bufId)
	{
		glDeleteBuffers(1, &m_gl_bufId);
		m_gl_bufId = 0;	
	}


	if (m_pixData)
	{
		delete m_pixData;
		m_pixData = 0;
	}
}

int DrawArrayWindow::gl_rebuild_texture(int width, int height)
{
	this->clear_gl_objects();


	m_width = width;
	m_height = height;


	// 4 bytes per pixel
	int dataSize = sizeof(unsigned char) * m_width * m_height * 4;
	 
	m_pixData = (unsigned char *) malloc(dataSize);
 
  

	glGenTextures(1, &m_gl_texId);        

	glBindTexture(GL_TEXTURE_2D, m_gl_texId);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, m_width, m_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, m_pixData);


	glGenBuffers(1, &m_gl_bufId);
	glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB,m_gl_bufId);
	glBufferData(GL_PIXEL_UNPACK_BUFFER_ARB, dataSize, m_pixData, GL_STREAM_COPY);
	  

	if (cudaGraphicsGLRegisterBuffer(&m_cuda_gResource, m_gl_bufId, cudaGraphicsMapFlagsWriteDiscard)  == cudaSuccess)
	{
		fprintf(stderr, "GL texture linked to CUDA\n");
	}
	else
	{
		fprintf(stderr, "Error: unable to link GL texture to CUDA!\n");

		return -1;
	}

	return 0;
}


void DrawArrayWindow::gl_draw_frame()
{
	display();
	glutMainLoopEvent();
}



void DrawArrayWindow::gl_deinit()
{
	this->clear_gl_objects();

	if (currentWindow)
	{
		currentWindow = 0;
	}
}



void DrawArrayWindow::fill_texture()
{
	if (m_cuda_gResource && m_gpu_original)
	{
		size_t num_bytes;
		unsigned char * devPixData;

		cudaError_t err = cudaSuccess;


		dim3 threadsPerSquareBlock(16, 16);

		dim3 numSquareBlocks((m_width + threadsPerSquareBlock.x - 1) / threadsPerSquareBlock.x, (m_height + threadsPerSquareBlock.y - 1) / threadsPerSquareBlock.y);


		cudaGraphicsMapResources(1, &m_cuda_gResource, 0);
  
		cudaGraphicsResourceGetMappedPointer((void**)&devPixData, &num_bytes, m_cuda_gResource);
 

		gpu_gen_buffer<<<numSquareBlocks, threadsPerSquareBlock>>>(m_width, m_height, devPixData, (0.0f + m_original_width) / m_width, (0.0f + m_original_height) / m_height, m_original_pitched_width, m_gpu_original, m_gpu_h, m_brighten);//,  m_gpu_b2x);


		err = cudaGetLastError();

		if (err != cudaSuccess)
		{
			fprintf(stderr, "DRAW: Failed to launch GPU algorithm  (error code %s)!\n", cudaGetErrorString(err));
			//exit(EXIT_FAILURE);
		}
   
		//cudaDeviceSynchronize();

		cudaGraphicsUnmapResources(1, &m_cuda_gResource, 0);


	}

}



DrawArrayWindow::~DrawArrayWindow()
{
	this->gl_deinit();
}


void display(void)
{
	if (currentWindow)
	{
		currentWindow->fill_texture();
	
		glClear(GL_COLOR_BUFFER_BIT);

		glColor3f(1.0f, 1.0f, 1.0f);
		glBindTexture(GL_TEXTURE_2D, currentWindow->m_gl_texId);
		glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB,currentWindow->m_gl_bufId);
 
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, currentWindow->m_width, currentWindow->m_height, GL_RGBA, GL_UNSIGNED_BYTE, 0);



		glBegin(GL_QUADS);
		glVertex2f(0, 0);
		glTexCoord2f(0, 0);
		glVertex2f(0, 1);
		glTexCoord2f(1, 0);
		glVertex2f(1, 1);
		glTexCoord2f(1, 1);
		glVertex2f(1, 0);
		glTexCoord2f(0, 1);
		glEnd();
		glBindTexture(GL_TEXTURE_2D, 0);
		glutSwapBuffers();

	}
}



void keyboard(unsigned char key, int /*x*/, int /*y*/)
{

    switch (key)
    {
        case 27:
        case 'q':
        case 'Q':
            printf("Shutting down...\n");

            glutDestroyWindow(glutGetWindow());
            break;

        default:
            break;
    }
}

void reshape(int x, int y)
{
    glViewport(0, 0, x, y);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 1, 0, 1, 0, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	if (currentWindow)
	{
		currentWindow->gl_rebuild_texture(x, y);
	}	
}

void cleanup(void)
{
	if (currentWindow)
	{
		currentWindow->clear_gl_objects();
	}
	
}
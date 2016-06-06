#define DAW_START_WINDOW_WIDTH 1024
#define DAW_START_WINDOW_HEIGHT 1024


struct cudaArray;
struct cudaGraphicsResource;


class DrawArrayWindow
{

public:

	DrawArrayWindow();



	~DrawArrayWindow();

	void set_data_to_display(float * gpu_zdata, float * gpu_hdata, int w, int h, int p_w);
	int gl_init(int device);

	int gl_rebuild_texture(int width, int height);

	void gl_deinit();

	void gl_draw_frame();

	void clear_gl_objects();


	void fill_texture();


	//unsigned int m_gl_bufferId;
	unsigned int m_gl_texId;
	unsigned int m_gl_bufId;


    cudaGraphicsResource * m_cuda_gResource;

	int m_width;
	int m_height;

	float m_brighten;

private:

	//cudaArray * m_texArray;

	unsigned char * m_pixData;

	float * m_gpu_original;

	float * m_gpu_h;

	int m_original_width;
	int m_original_height;

	int m_original_pitched_width;


};




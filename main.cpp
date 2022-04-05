//
#include <cstdio>
#include <vector>
//
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"

class Image
{
public:
	Image()
	{
	}
	void load(const char* filename)
	{
		int comp;
		float* pix = stbi_loadf(filename, &width_, &height_, &comp, 0);
		pixels_.resize(width_* height_);
		for(int i=0;i<pixels_.size();++i)
		{
			pixels_[i] = {pix[i*3+0], pix[i*3+1], pix[i*3+2]};
		}
		STBI_FREE(pix);
	}
	void savePFM(const char* filename)
	{
		FILE* file = fopen(filename,"wb");
		fprintf(file, "%s",  "PF\n");
		fprintf(file,"%ld %ld\n", width_, height_);
		fprintf(file, "%s",  "-1.000000\n");
		fwrite(pixels_.data(),  sizeof(Pixel), pixels_.size(), file);
		fclose(file);
	}
private:
	struct Pixel
	{
		double r;
		double g;
		double b;
	};
	std::vector<Pixel> pixels_;
	int32_t width_  = 0;
	int32_t height_ = 0;
};

void main()
{
	// TODO: hdrÇÉçÅ[ÉhÇµÇƒpmfÇ≈ï€ë∂Ç∑ÇÈ
	Image image;
	image.load("../pt.hdr");
	image.savePFM("../pt.pmf");
}
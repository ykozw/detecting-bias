#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <vector>
//
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"

#define MAX_P 1.f
#define M_PI 3.141592653589f


// computes area left of t for t-distribution
// with v degrees of freedom and t<=0
// in double
// from stack overflow:
static inline double tcdf64(const double t, const int v)
{
	// Algorithm 3: The Integral of Student's t-distribution (Applied statistics 1968 vol 17 p 189, B.E. Cooper)
	const double x = v / (v + t * t);
	double c = 1.0f, s = 1.0f;
	int ioe = v % 2, k = 2 + ioe;
	if (v < 1) return 0.0f;
	if (v >= 4)
	{
		while (k <= v - 2)
		{
			c *= x - x / k;
			s += c;
			k += 2;
		}
	}
	c = t / sqrtf(v);
	if (ioe != 1)
		return 0.5f + 0.5f * sqrtf(x) * c * s;
	return 0.5f + ((v == 1 ? 0 : x * c * s) + atanf(c)) / M_PI;
}

// computes area left of t for t-distribution
// with v degrees of freedom and t<=0
// in float
// from stack overflow:
static inline float tcdf(float t, int v)
{
	// Algorithm 3: The Integral of Student's t-distribution (Applied statistics 1968 vol 17 p 189, B.E. Cooper)
	const float b = v / (v + t * t);
	float c = 1.0f, s = 1.0f;
	int ioe = v % 2, k = 2 + ioe;
	if (v < 1) return 0.0f;
	if (v >= 4)
	{
		while (k <= v - 2)
		{
			c *= b - b / k;
			s += c;
			k += 2;
		}
	}
	c = t / sqrtf(v);
	if (ioe != 1)
		return 0.5f + 0.5f * sqrtf(b) * c * s;
	return 0.5f + ((v == 1 ? 0 : b * c * s) + atanf(c)) / M_PI;
}

// color coding
static inline void viridis_quintic(
	float x,        // input float between 0 and 1
	float* col)     // output colour ramp
{
	x = fminf(1.0f, x);
	const float x2 = x * x;
	const float x3 = x2 * x;
	const float x4 = x2 * x2;
	const float x5 = x3 * x2;
	col[0] = +0.280268003 - 0.143510503 * x + 2.225793877 * x2 - 14.815088879 * x3 + 25.212752309 * x4 - 11.772589584 * x5;
	col[1] = -0.002117546 + 1.617109353 * x - 1.909305070 * x2 + 2.701152864 * x3 - 1.685288385 * x4 + 0.178738871 * x5;
	col[2] = +0.300805501 + 2.614650302 * x - 12.019139090 * x2 + 28.933559110 * x3 - 33.491294770 * x4 + 13.762053843 * x5;
}

// approximate erfi Abramowitz , M., and Stegun , I. A. 1964. Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables

/* Indirect computation for debugging:
   float erf_approx3(float t)
   {
   const float x = 1.0f / (1.0f + 0.4747f * t);
   const float a1 = 0.3480242;
   const float    a2 = -0.0958798;
   const float    a3 = 0.7478556;
   printf("t = %f, ", t);
   printf("x = %f, ", x);
   printf("e^(-x*x) = %f, ", exp(-x*x));
   printf("a1 + a2 + a3 = %f ", a1 + a2 + a3);
   const float val= 1.0 - (a1 * x + a2 * x * x + a3 * x * x * x) * exp(- t * t);
   printf("erf(%f)=%f, ", t, val);
   return val;
   }

   float phi(float x)
   {
   return 0.5 * (1.f + erf_approx3(x / sqrt(2.f)));
   }
   */

static inline float gauss_cdf(float t)
{
	float tt = t > 0 ? t : -t;
	tt = tt / sqrt(2.0f);

	const float x = 1.f / (1.f + 0.47047 * tt);
	const float erf = 1.0 - x * (0.3480242
		+ x * (-0.0958798
			+ 0.7478556 * x))
		* exp(-tt * tt);
	return 1.f - erf;
}

//
struct Pixel
{
	double r;
	double g;
	double b;
public:
	Pixel()
		:r(0.0), g(0.0), b(0.0)
	{
	}
	Pixel(float v)
		:r(v), g(v), b(v)
	{
	}
	Pixel(double _r, double _g, double _b)
		:r(_r), g(_g), b(_b)
	{
	}
	double& operator[](int32_t idx)
	{
		return *(&r + idx);
	}
	const double& operator[](int32_t idx)const
	{
		return *(&r + idx);
	}
	Pixel operator + (const Pixel& rhs) const
	{
		return Pixel(
			r + rhs.r,
			g + rhs.g,
			b + rhs.b
		);
	}
	Pixel operator * (const Pixel& rhs) const
	{
		return Pixel(
			r * rhs.r,
			g * rhs.g,
			b * rhs.b
		);
	}
	Pixel& operator += (const Pixel& rhs)
	{
		*this = *this + rhs;
		return *this;
	}
};

//
template<typename PixelType>
class Image
{
public:
public:
	Image()
	{
	}
	Image(int32_t width, int32_t height)
	{
		width_ = width;
		height_ = height;
		pixels_.resize(width * height);
	}
	void load(const char* filename)
	{
		int comp;
		float* pix = stbi_loadf(filename, &width_, &height_, &comp, 0);
		pixels_.resize(width_ * height_);
		for (int i = 0; i < pixels_.size(); ++i)
		{
			pixels_[i] = { pix[i * 3 + 0], pix[i * 3 + 1], pix[i * 3 + 2] };
		}
		STBI_FREE(pix);
	}
	void save(const char* filename)
	{
		std::vector<float> arr(pixels_.size() * 3);
		for (int pi = 0; pi < pixels_.size(); ++pi)
		{
			const auto& p = pixels_[pi];
			arr[pi * 3 + 0] = p.r;
			arr[pi * 3 + 1] = p.g;
			arr[pi * 3 + 2] = p.b;
		}
		stbi_write_hdr(filename, width_, height_, 3, arr.data());
	}


	int32_t width() const
	{
		return width_;
	}
	int32_t height() const
	{
		return height_;
	}
	PixelType& operator()(int32_t x, int32_t y)
	{
		return pixels_[x + y * width_];
	}
	const PixelType& operator()(int32_t x, int32_t y) const
	{
		return pixels_[x + y * width_];
	}
private:

	std::vector<PixelType> pixels_;
	int32_t width_ = 0;
	int32_t height_ = 0;
};
//
Image<Pixel> sum_of_square_welch(const Image<Pixel>& src, const int block_size, float scale)
{
	const int tw = src.width() / block_size;
	const int th = src.height() / block_size;
	Image<Pixel> ret(tw, th);
	for (int tx = 0; tx < tw; ++tx)
	{
		for (int ty = 0; ty < th; ++ty)
		{
			const int bx = tx * block_size;
			const int by = ty * block_size;
			Pixel sum;
			for (int ox = 0; ox < block_size; ++ox)
			{
				for (int oy = 0; oy < block_size; ++oy)
				{
					auto t = src(bx + ox, by + oy);
					sum += t * t;
				}
			}
			ret(tx, ty) = sum * Pixel(scale);
		}
	}
	return ret;
}
//
Image<Pixel> sum_welch(const Image<Pixel>& src, const int block_size, float scale)
{
	const int tw = src.width() / block_size;
	const int th = src.height() / block_size;
	Image<Pixel> ret(tw, th);
	for (int tx = 0; tx < tw; ++tx)
	{
		for (int ty = 0; ty < th; ++ty)
		{
			const int bx = tx * block_size;
			const int by = ty * block_size;
			Pixel sum;
			for (int ox = 0; ox < block_size; ++ox)
			{
				for (int oy = 0; oy < block_size; ++oy)
				{
					sum += src(bx + ox, by + oy);
				}
			}
			ret(tx, ty) = sum * Pixel(scale);
		}
	}
	return ret;
}

/*
 * === Usage:
 *
 * Pass first and second image as arguments
 *
 * If the image is called image_fb00.pfm there
 * should also be a file called image_welch (no extension).
 *
 * Also works with double frame buffer image_fbd00.pfm.
 *
 * The welch file doesn't need a header, the width
 * and height are extracted from the .pfm image.
 *
 *
 * === Useful resources:
 *
 * Sample Variance:
 * http://www.statisticshowto.com/sample-variance/
 *
 * Welch test intro and example:
 * https://www.youtube.com/watch?v=2-ecXltt2vI
 * https://www.youtube.com/watch?v=gzrmHpA54Sc
 *
 * p-value table
 * https://math.stackexchange.com/questions/808474/p-value-for-lower-upper-tailed-t-test
 */
int main(int argc, char* argv[])
{
	//
	float p_scale = 1.f;
	const char* outname_prefix = "";
	int replace_with_gauss = 0;
	int gauss_nu_threshold = 0;
	int visualization_mode = 0; //0 = regular or scaled color scale, 1 = square root
	int block_size = 32;
	int vis_three = 0;
	//
	printf("pscale %f\n", p_scale);
	printf("blocksize %d\n", block_size);
	printf("replace gauss? %d, threshold %d \n", replace_with_gauss, gauss_nu_threshold);
	printf("Write to prefix %s\n", outname_prefix);
	printf("Vis mode: %d\n", visualization_mode);

	// read img1 and img2 (only needed to extract width, height and debugging info)
	Image<Pixel> pixels1;
	pixels1.load("../bin/pt128spp.hdr");
	const int32_t width = pixels1.width();
	const int32_t height = pixels1.height();
	Image<Pixel> pixels2;
	pixels2.load("../bin/pt1024spp.hdr");
	Image<Pixel> output(pixels1.width(), pixels1.height());
	const uint64_t w_wd = pixels1.width() / block_size;
	const uint64_t w_ht = pixels1.height() / block_size;

	//
	const double welchsamples1 = 1024.0f;
	const double welchsamples2 = 1024.0f;
	//
	const Image<Pixel> welch1_2 = sum_of_square_welch(pixels1, block_size, 1.0f);
	const Image<Pixel> welch2_2 = sum_of_square_welch(pixels2, block_size, 1.0f);
	const Image<Pixel> welch1_1 = sum_welch(pixels1, block_size, 1.0f);
	const Image<Pixel> welch2_1 = sum_welch(pixels2, block_size, 1.0f);
	//
	std::vector<double> pvals(w_wd * w_ht * 3);
	int pvalcnt = 0; // only fill p-values where they were actually computed. This avoids distortion of the p-value histogram if the images are partially black.
	//
	std::vector<int> nuvals(w_wd * w_ht * 3);

	

	const double n1 = welchsamples1;
	const double n2 = welchsamples2;

	printf("N1 = %f, N2 = %f\n", n1, n2);
	// this works if the last field of the squared welch files additionally contains the spp count:
	//printf("spp1 = %f, spp2 = %f\n", welch1_2[w_wd*w_ht*3+1], welch2_2[w_wd*w_ht*3+1]);

	// avoid division by zero
	const double n1inv = n1 > 0.0 ? 1.0 / n1 : 0.0;
	const double n2inv = n2 > 0.0 ? 1.0 / n2 : 0.0;

	// iterate w_ht * w_wd pixel blocks
	for (int j = 0; j < w_ht; j++)
	{
		for (int i = 0; i < w_wd; i++)
		{
			// gather sample sum and sum of squares per color channel
			double sumX1[3] = { 0.0 };
			double sumX2[3] = { 0.0 };
			for (int k = 0; k < 3; k++)
			{
				sumX1[k] = welch1_1(i, j)[k];
				sumX2[k] = welch2_1(i, j)[k];
			}

			double s1_2[3], s2_2[3], tmp[3], t[3];
			int nu[3];
			float cdf = 1.0;
			float p_values[3];

			// compute statistics per color channel
			for (int k = 0; k < 3; k++)
			{
				// unbiased sample variance
				s1_2[k] = (1.0 / (n1 - 1.0)) * (welch1_2(i, j)[k] - (n1inv * sumX1[k]) * sumX1[k]);
				s2_2[k] = (1.0 / (n2 - 1.0)) * (welch2_2(i, j)[k] - (n2inv * sumX2[k]) * sumX2[k]);
				//
				assert(s1_2[k] >= 0.0);
				assert(s2_2[k] >= 0.0);

				tmp[k] = s1_2[k] * n1inv + s2_2[k] * n2inv;
				if (tmp[k] == 0 || !(tmp[k] < FLT_MAX))
				{
					cdf = 1.0f;
					continue;
				}
				assert(tmp[k] >= 0.0);
				assert(tmp[k] == tmp[k]);

				// test lower bound of cdf interval:
				// t = (X_1 - X_2) / sqrt (s_1²/N_1 + s_2²/N_2)

				// compute t-statistic
				t[k] = -fabs(sumX1[k] * n1inv - sumX2[k] * n2inv) / sqrtf(tmp[k]);
				if (!(t[k] == t[k]) // take care of NAN
					|| t[k] == 0 || !(t[k] < FLT_MAX))
				{
					cdf = 1.0f;
					continue;
				}

				// approximate degrees of freedom with Welch-Satterthwaite equation
				nu[k] = roundf(tmp[k] * tmp[k] / (
					s1_2[k] * s1_2[k] * n1inv * n1inv / (n1 - 1.0f) +
					s2_2[k] * s2_2[k] * n2inv * n2inv / (n2 - 1.0f)
					));
				if (i == 0 && j == 0) // logging for first pixel block
				{
					printf("tmp[%d] = %f\n", k, tmp[k]);
					printf("s1_2[%d] = %f\n", k, s1_2[k]);
					printf("s2_2[%d] = %f\n", k, s2_2[k]);
					printf("n1inv = %f\n", n1inv);
					printf("n2inv = %f\n", n2inv);
					printf("nu[%d] = %d\n", k, nu[k]);
				}

				// now test the null-hypothesis that the two means are equal, using
				// a two-tailed test on the t-distribution:
				const float p_value = (replace_with_gauss && nu[k] >= gauss_nu_threshold) ? gauss_cdf(t[k]) : 2.f * tcdf64(t[k], nu[k]);
				cdf = fminf(cdf, p_value); // this cdf is the p-value. The smallest p-value is at the "worst" color channel

				assert(p_value >= -0.00001f && p_value <= 1.f);
				p_values[k] = p_value;


				// if we do this we also get a p-value of 0 for every index where no p-value was actually computed, e.g. if t[k] or tmp[k] = 0:
				// This has the advantage of knowing which p-value corresponds to each pixel, but will skew the p-value histogram:
				// pvals[3*(j*w_wd+i)+k] = p_value;

				// Instead, we keep track of how many pvalues were actually computed.
				// These values are useful for histograms, but we can't reconstruct
				// which value corresponds to which pixel:
				pvals[pvalcnt] = p_value;
				pvalcnt++;

				nuvals[3 * (j * w_wd + i) + k] = nu[k];
			}

			// visualize p-values in color map
			if (vis_three)
			{
				for (int kk = 0; kk < 3; kk++)
				{
					float col[3];
					if (visualization_mode == 1)
					{
						p_values[kk] = sqrt(p_values[kk]);
					}
					const float confidence = fmaxf(0.f, 1.f - (p_values[kk] > 0.f ? (p_scale * p_values[kk]) : 0.f));
					viridis_quintic(confidence / MAX_P, col);

					for (int jj = 0; jj < block_size; jj++)
					{
						for (int ii = (kk * block_size) / 3; ii < ((kk + 1) * block_size) / 3; ii++)
						{
							for (int k = 0; k < 3; k++)
							{
								output(block_size * i + ii, block_size * j + jj)[k] = col[k];
							}
						}
					}
				}
			}
			else
			{
				float col[3];
				//the color method is yellow for large values and pink for values close to zero. The p-value is "good" when it is close to one, so use 1-cdf as argument here
				if (visualization_mode == 1)
				{
					cdf = sqrtf(cdf);
				}
				const float confidence = fmaxf(0.f, 1.f - (cdf > 0.f ? (p_scale * cdf) : 0.f));
				viridis_quintic(confidence / MAX_P, col);
				for (int jj = 0; jj < block_size; jj++)
				{
					for (int ii = 0; ii < block_size; ii++)
					{
						for (int k = 0; k < 3; k++)
						{
							output(block_size * i + ii, block_size * j + jj)[k] = col[k];
						}
					}
				}
			}
		}
	}

	// draw color scale:
	if (false)
	{
		for (int j = 0; j < height; j++)
		{
			for (int i = width - 16; i < width - 8; i++)
			{
				float f = 1.0f - j / (height - 1.0f);
				for (int k = 0; k < 3; k++)
				{
					auto& p = output(i, j);
					p[k] = f;
				}
			}
		}

		for (int j = 0; j < height; j++) for (int i = width - 8; i < width; i++)
		{
			float f = 1.0f - j / (height - 1.0f), col[3];
			viridis_quintic(f / MAX_P, col);
			for (int k = 0; k < 3; k++)
			{
				output(i, j)[k] = col[k];
			}
		}
	}

	output.save("../bin/result.hdr");

	if (true)
	{
		// write out p-values
		FILE* f1 = fopen("pvalues.txt", "wb");
		if (f1)
		{
			// print number of p-values that were actually computed
			fprintf(f1, "%d ", pvalcnt);
			for (int k = 0; k < pvalcnt; k++)
			{
				fprintf(f1, "%f ", pvals[k]);
			}
			fclose(f1);
		}
	}

	if (true)
	{
		// write out nu values
		FILE* f2 = fopen("nuvalues.txt", "wb");
		if (f2)
		{
			for (int k = 0; k < w_wd * w_ht * 3; k++)
			{
				fprintf(f2, "%d ", nuvals[k]);
			}
			fclose(f2);
		}
}

#if 0 // debug print the cdf function
	for (int k = 0; k < 100; k++)
	{
		float t = -4.0 + 8.0 * k / 100.0f;
		fprintf(stderr, "%g %g\n", t, tcdf(t, 7));
	}
#endif
}


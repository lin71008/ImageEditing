///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <random>
#include <map>

using namespace std;

// constants
#define RED   0  // red channel
#define GREEN 1  // green channel
#define BLUE  2  // blue channel
#define ALPHA 3  // alpha channel
    
const unsigned char RGBA32_BLACK[] = {   0,   0,   0, 255 };
const unsigned char RGBA32_WHITE[] = { 255, 255, 255, 255 };
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color

// snippet
#define POS_XY(SIZE, X_MAX, Y, X) \
    ((SIZE) * ((Y) * (X_MAX) + (X)))

#define GET_RGBA32(COLOR_MAP, POS, COLORS) \
    (COLORS)[ RED   ] = (COLOR_MAP)[ (POS) + RED   ];\
    (COLORS)[ GREEN ] = (COLOR_MAP)[ (POS) + GREEN ];\
    (COLORS)[ BLUE  ] = (COLOR_MAP)[ (POS) + BLUE  ];\
    (COLORS)[ ALPHA ] = (COLOR_MAP)[ (POS) + ALPHA ]

#define SET_RGBA32(COLOR_MAP, POS, COLORS) \
    (COLOR_MAP)[ (POS) + RED   ] = (COLORS)[ RED   ];\
    (COLOR_MAP)[ (POS) + GREEN ] = (COLORS)[ GREEN ];\
    (COLOR_MAP)[ (POS) + BLUE  ] = (COLORS)[ BLUE  ];\
    (COLOR_MAP)[ (POS) + ALPHA ] = (COLORS)[ ALPHA ]

#define DIFF_RGBA32_NAIVE(COLOR_MAP, POS, COLORS) \
    (COLOR_MAP)[ (POS) + RED   ] -= (COLORS)[ RED   ];\
    (COLOR_MAP)[ (POS) + GREEN ] -= (COLORS)[ GREEN ];\
    (COLOR_MAP)[ (POS) + BLUE  ] -= (COLORS)[ BLUE  ];\
    (COLOR_MAP)[ (POS) + ALPHA ] -= (COLORS)[ ALPHA ]

#define ADD_RGBA32_NAIVE(COLOR_MAP, POS, COLORS, SCALE) \
    (COLOR_MAP)[ (POS) + RED   ] += ((COLORS)[ RED   ]) * (SCALE);\
    (COLOR_MAP)[ (POS) + GREEN ] += ((COLORS)[ GREEN ]) * (SCALE);\
    (COLOR_MAP)[ (POS) + BLUE  ] += ((COLORS)[ BLUE  ]) * (SCALE);\
    (COLOR_MAP)[ (POS) + ALPHA ] += ((COLORS)[ ALPHA ]) * (SCALE)

#define GET_R8G8B8A8(COLOR_MAP, POS, R, G, B, A) \
    (R) = (COLOR_MAP)[ (POS) + RED   ];\
    (G) = (COLOR_MAP)[ (POS) + GREEN ];\
    (B) = (COLOR_MAP)[ (POS) + BLUE  ];\
    (A) = (COLOR_MAP)[ (POS) + ALPHA ]

#define SET_R8G8B8A8(COLOR_MAP, POS, R, G, B, A) \
    (COLOR_MAP)[ (POS) + RED   ] = R;\
    (COLOR_MAP)[ (POS) + GREEN ] = G;\
    (COLOR_MAP)[ (POS) + BLUE  ] = B;\
    (COLOR_MAP)[ (POS) + ALPHA ] = A

// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

    return res;
}// Binomial


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char *d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
	    data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char   *rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
        {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
	    return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
	    cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
	    return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
	    width = height = 0;
	    return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
    for (int i = 0; i < 4 *(width * height); i += 4)
    {
        int r, g, b, a, gray;
        GET_R8G8B8A8(data, i, r, g, b, a);
        gray = 0.299 * r + 0.587 * g + 0.114 * b;
        SET_R8G8B8A8(data, i, gray, gray, gray, a);
    }
    return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
    for (int i = 0; i < 4 *(width * height); i += 4)
    {
        unsigned char r, g, b, a;
        GET_R8G8B8A8(data, i, r, g, b, a);
        SET_R8G8B8A8(data, i, r & 224, g & 224, b & 192, a);  // R3G3B2
    }
    return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
    map <tuple<int, int, int>, int> reduce_color_counter;
    
    // reduce channel size to 5
    for (int i = 0; i < 4 * width * height; i += 4)
    {
        int r, g, b, a;
        GET_R8G8B8A8(data, i, r, g, b, a);
        tuple<int, int, int> key = make_tuple(r & 248, g & 248, b & 248);
        if (reduce_color_counter.find(key) != reduce_color_counter.end())
        {
            reduce_color_counter[key]++;
        }
        else
        {
            reduce_color_counter[key] = 1;
        }
    }
    
    // reduce color table
    vector<tuple<int, int, int>> reduce_color;
    vector<pair<tuple<int, int, int>, int>> pre_reduce_color(reduce_color_counter.begin(), reduce_color_counter.end());
    sort(pre_reduce_color.begin(), pre_reduce_color.end(), [](const auto& x, const auto& y) {return x.second > y.second; });
    for (int i = 0; i < 256 && i < reduce_color_counter.size(); i++)
    {
        reduce_color.push_back(pre_reduce_color[i].first);
    }

    for (int i = 0; i < 4 * width * height; i += 4)
    {
        int r, g, b, a;
        GET_R8G8B8A8(data, i, r, g, b, a);
        tuple<int, int, int> key = make_tuple(r & 248, g & 248, b & 248);
        auto it = find(reduce_color.begin(), reduce_color.end(), key);
        if (!(it != reduce_color.end()))
        {
            // find closest chosen color
            int k, closest_factor = 1e6, closest_color_index = 0;
            for (int j = 0; j < 256; ++j)
            {
                k = (r - get<0>(reduce_color[j])) * (r - get<0>(reduce_color[j])) + \
                    (g - get<1>(reduce_color[j])) * (g - get<1>(reduce_color[j])) + \
                    (b - get<2>(reduce_color[j])) * (b - get<2>(reduce_color[j]));
                
                if (k < closest_factor)
                {
                    closest_factor = k;
                    closest_color_index = j;
                } 
            }
            SET_R8G8B8A8(data, i, get<0>(reduce_color[closest_color_index]), \
                                  get<1>(reduce_color[closest_color_index]), \
                                  get<2>(reduce_color[closest_color_index]), a);
        }
        else
        {
            SET_R8G8B8A8(data, i, r & 248, g & 248, b & 248, a);
        }
    }
    return true;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
    To_Grayscale();

	for (int i = 0; i < 4 *(width * height); i += 4)
    {
        if (data[i] > 128)
        {
            SET_RGBA32(data, i, RGBA32_WHITE);
        }
        else
        {
            SET_RGBA32(data, i, RGBA32_BLACK);
        }
    }
    return true;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
    To_Grayscale();

    default_random_engine generator;
    uniform_real_distribution<float> distribution(-0.2, 0.2);

    for (int i = 0; i < 4 * (width * height); i+=4)
    {
        if ((data[i] / 255.0) + distribution(generator) > 0.5)
        {
            SET_RGBA32(data, i, RGBA32_WHITE);
        }
        else
        {
            SET_RGBA32(data, i, RGBA32_BLACK);
        }
    }
    return true;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
    To_Grayscale();

    // expend channel size
    float* Image;
    Image = new float [width * height];
    memset(Image, 0, sizeof(float) * ((size_t) width) * ((size_t) height));
    for (int i = 0; i < width * height; i++)
    {
        Image[i] = data[4 * i] / 255.0;
    }

    float old_pixel, new_pixel, quant_error;
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            old_pixel = Image[POS_XY(1, width, i, j)];
            new_pixel = (old_pixel > 0.5) ? 1.0 : 0.0;
            quant_error = old_pixel - new_pixel;

            if (!(j == width - 1))                   Image[POS_XY(1, width, i, j + 1)]     += quant_error * 7 / 16;
            if (!(i == height - 1 || j == 0))         Image[POS_XY(1, width, i + 1, j - 1)] += quant_error * 3 / 16;
            if (!(i == height - 1))                   Image[POS_XY(1, width, i + 1, j)]     += quant_error * 5 / 16;
            if (!(i == height - 1 || j == width - 1)) Image[POS_XY(1, width, i + 1, j + 1)] += quant_error * 1 / 16;
        }
    }

    // reduce channel size
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (Image[POS_XY(1, width, i, j)] > 0.5)
            {
                SET_RGBA32(data, POS_XY(4, width, i, j), RGBA32_WHITE);
            }
            else
            {
                SET_RGBA32(data, POS_XY(4, width, i, j), RGBA32_BLACK);
            }
        }
    }
    delete[] Image;

    return true;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
    To_Grayscale();

    // initial
    int avg = 0;
    vector<int> bright;

    // get avg
    for (int i = 0; i < 4 * width * height; i += 4)
    {
        bright.push_back(data[i]);
        avg += data[i];
    }

    // pivot: N % brightness
    int pivot = (width * height)  - (avg / 255);

    sort(bright.begin(), bright.end(), [](const auto& x, const auto& y) {return x < y; });
    int threshold = bright[pivot];

    for (int i = 0; i < 4 * width * height; i += 4)
    {
        if (data[i] > threshold)
        {
            SET_RGBA32(data, i, RGBA32_WHITE);
        }
        else
        {
            SET_RGBA32(data, i, RGBA32_BLACK);
        }
    }
    return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
    To_Grayscale();

    // constants
    const float mask[4][4] = {{ 0.7059, 0.3529, 0.5882, 0.2353},
                              { 0.0588, 0.9412, 0.8235, 0.4118},
                              { 0.4706, 0.7647, 0.8824, 0.1176},
                              { 0.1765, 0.5294, 0.2941, 0.6471}};

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (data[POS_XY(4, width, i, j)] / 255.0 > mask[i % 4][j % 4])
            {
                SET_RGBA32(data, POS_XY(4, width, i, j), RGBA32_WHITE);
            }
            else
            {
                SET_RGBA32(data, POS_XY(4, width, i, j), RGBA32_BLACK);
            }
        }
    }
    return true;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
    // constants
    const int avaliable_red_channel_vlaue[]   = { 0, 36, 73, 109, 146, 182, 219, 255 };  // 3 bits  // => total 8 bits
    const int avaliable_green_channel_vlaue[] = { 0, 36, 73, 109, 146, 182, 219, 255 };  // 3 bits  // => total 8 bits
    const int avaliable_blue_channel_vlaue[]  = { 0, 85, 170, 255 };       // 2 bits                // => total 8 bits
    
    auto find_closest_value = [](const int Arr[], const int Arr_size, int n) -> int\
    { \
        int result_index = 0; \
        for (int i = 0; i < Arr_size; ++i) \
        { \
            if (n - Arr[i] > 0) \
            { \
                result_index = i; \
            } \
        } \
        return Arr[result_index]; \
    };

    auto find_closest_color = [=](int color[], int result[]) -> void \
    { \
        result[ RED   ] = (int) find_closest_value(avaliable_red_channel_vlaue,   8, color[ RED   ]); \
        result[ GREEN ] = (int) find_closest_value(avaliable_green_channel_vlaue, 8, color[ GREEN ]); \
        result[ BLUE  ] = (int) find_closest_value(avaliable_blue_channel_vlaue,  4, color[ BLUE  ]); \
        result[ ALPHA ] = color[ ALPHA ]; \
    };

    int old_pixel[4], new_pixel[4], quant_error[4];

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            GET_RGBA32(data, POS_XY(4, width, i, j), old_pixel);
            find_closest_color(old_pixel, new_pixel);
            SET_RGBA32(data, POS_XY(4, width, i, j), new_pixel);
            SET_RGBA32(quant_error, 0, old_pixel);
            DIFF_RGBA32_NAIVE(quant_error, 0, new_pixel);

            if (!(j == width - 1))                   { ADD_RGBA32_NAIVE(data, POS_XY(4, width, i, j + 1),     quant_error, (7.0 / 16)); }
            if (!(i == height - 1 || j == 0))        { ADD_RGBA32_NAIVE(data, POS_XY(4, width, i + 1, j - 1), quant_error, (3.0 / 16)); }
            if (!(i == height -1))                   { ADD_RGBA32_NAIVE(data, POS_XY(4, width, i + 1, j),     quant_error, (5.0 / 16)); }
            if (!(i == height -1 || j == width - 1)) { ADD_RGBA32_NAIVE(data, POS_XY(4, width, i + 1, j + 1), quant_error, (1.0 / 16)); }
        }
    }

    return true;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
    }

    return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
    ClearToBlack();
    return false;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    ClearToBlack();
    return false;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    ClearToBlack();
    return false;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N( unsigned int N )
{
    ClearToBlack();
   return false;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    ClearToBlack();
    return false;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    ClearToBlack();
    return false;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
    ClearToBlack();
    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    ClearToBlack();
    return false;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    ClearToBlack();
    return false;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    ClearToBlack();
    return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
               unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}


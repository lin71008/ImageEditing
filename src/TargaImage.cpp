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
    
const unsigned char RGBA32_BLANK[] = {   0,   0,   0,   0 };
const unsigned char RGBA32_BLACK[] = {   0,   0,   0, 255 };
const unsigned char RGBA32_WHITE[] = { 255, 255, 255, 255 };
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color

// snippet
#define POS_XY(SIZE, X_MAX, Y, X) \
    ((SIZE) * ((Y) * (X_MAX) + (X)))

#define GET_RGBA32(COLOR_MAP, POS, COLORS) \
    (COLORS)[ RED   ] = (COLOR_MAP)[ (POS) + RED   ]; \
    (COLORS)[ GREEN ] = (COLOR_MAP)[ (POS) + GREEN ]; \
    (COLORS)[ BLUE  ] = (COLOR_MAP)[ (POS) + BLUE  ]; \
    (COLORS)[ ALPHA ] = (COLOR_MAP)[ (POS) + ALPHA ]; \

#define SET_RGBA32(COLOR_MAP, POS, COLORS) \
    (COLOR_MAP)[ (POS) + RED   ] = (COLORS)[ RED   ]; \
    (COLOR_MAP)[ (POS) + GREEN ] = (COLORS)[ GREEN ]; \
    (COLOR_MAP)[ (POS) + BLUE  ] = (COLORS)[ BLUE  ]; \
    (COLOR_MAP)[ (POS) + ALPHA ] = (COLORS)[ ALPHA ]; \

#define DIFF_RGBA32_NAIVE(COLOR_MAP, POS, COLORS) \
    (COLOR_MAP)[ (POS) + RED   ] -= (COLORS)[ RED   ]; \
    (COLOR_MAP)[ (POS) + GREEN ] -= (COLORS)[ GREEN ]; \
    (COLOR_MAP)[ (POS) + BLUE  ] -= (COLORS)[ BLUE  ]; \
    (COLOR_MAP)[ (POS) + ALPHA ] -= (COLORS)[ ALPHA ]; \

#define ADD_RGBA32_NAIVE(COLOR_MAP, POS, COLORS, SCALE) \
    (COLOR_MAP)[ (POS) + RED   ] += (((COLORS)[ RED   ]) * (SCALE)); \
    (COLOR_MAP)[ (POS) + GREEN ] += (((COLORS)[ GREEN ]) * (SCALE)); \
    (COLOR_MAP)[ (POS) + BLUE  ] += (((COLORS)[ BLUE  ]) * (SCALE)); \

#define GET_R8G8B8A8(COLOR_MAP, POS, R, G, B, A) \
    (R) = (COLOR_MAP)[ (POS) + RED   ]; \
    (G) = (COLOR_MAP)[ (POS) + GREEN ]; \
    (B) = (COLOR_MAP)[ (POS) + BLUE  ]; \
    (A) = (COLOR_MAP)[ (POS) + ALPHA ]; \

#define SET_R8G8B8A8(COLOR_MAP, POS, R, G, B, A) \
    (COLOR_MAP)[ (POS) + RED   ] = R; \
    (COLOR_MAP)[ (POS) + GREEN ] = G; \
    (COLOR_MAP)[ (POS) + BLUE  ] = B; \
    (COLOR_MAP)[ (POS) + ALPHA ] = A; \

static unsigned int __GET_RGBA32_SAFE = 0;
#define GET_RGBA32_SAFE(COLOR_MAP, SIZE, Y_MAX, X_MAX, Y, X, COLORS) \
    if ( (X) >= 0 && (X) < (X_MAX) && (Y) >= 0 && (Y) < (Y_MAX)) \
    { \
        GET_RGBA32((COLOR_MAP), POS_XY((SIZE), (X_MAX), (Y), (X)), (COLORS));\
        __GET_RGBA32_SAFE = 1; \
    } \
    else \
    { \
        SET_RGBA32((COLORS), 0, RGBA32_BLANK); \
        __GET_RGBA32_SAFE = 0; \
    }

#define GET_2D_CONVOLUTION(F, F_SIZE, F_Y_MAX, F_X_MAX, G, G_Y_MAX, G_X_MAX, RESULT, Y, X, SCALE) \
    { \
        bool __not_blank = false; \
        double __pixel_buffer_0[4], __pixel_buffer_1[4]; \
        SET_RGBA32(__pixel_buffer_0, 0, RGBA32_BLANK);\
        for (int __i = 0; __i < (G_Y_MAX); __i++) \
        { \
            for (int __j = 0; __j < (G_X_MAX); __j++) \
            { \
                GET_RGBA32_SAFE((F), (F_SIZE), (F_Y_MAX), (F_X_MAX), (Y) + __i, (X) + __j, __pixel_buffer_1);\
                if ( __GET_RGBA32_SAFE ) \
                { \
                    __not_blank = true; \
                } \
                ADD_RGBA32_NAIVE(__pixel_buffer_0, 0, __pixel_buffer_1, (G)[__i][__j]); \
            } \
        } \
        if ( __not_blank ) \
        { \
            SET_RGBA32((RESULT), 0, RGBA32_BLACK); \
            ADD_RGBA32_NAIVE((RESULT), 0, __pixel_buffer_0, (SCALE)); \
        } \
        else \
        { \
            SET_RGBA32((RESULT), 0, RGBA32_BLANK); \
        } \
    }

#define SET_RGBA32_CUT(COLOR_MAP, POS, COLORS) \
    (COLOR_MAP)[ (POS) + RED   ] = ((COLORS)[ RED   ]) > 0 ? (((COLORS)[ RED   ]) < 255 ? ((COLORS)[ RED   ]) : 255) : 0; \
    (COLOR_MAP)[ (POS) + GREEN ] = ((COLORS)[ GREEN ]) > 0 ? (((COLORS)[ GREEN ]) < 255 ? ((COLORS)[ GREEN ]) : 255) : 0; \
    (COLOR_MAP)[ (POS) + BLUE  ] = ((COLORS)[ BLUE  ]) > 0 ? (((COLORS)[ BLUE  ]) < 255 ? ((COLORS)[ BLUE  ]) : 255) : 0; \
    (COLOR_MAP)[ (POS) + ALPHA ] = ((COLORS)[ ALPHA ]) > 0 ? (((COLORS)[ ALPHA ]) < 255 ? ((COLORS)[ ALPHA ]) : 255) : 0; \


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
        int factor = abs(n - Arr[0]), result_index = 0; \
        for (int i = 0; i < Arr_size; ++i) \
        { \
            if (abs(n - Arr[i]) < factor) \
            { \
                factor = abs(n - Arr[i]);\
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

    // expend channel size
    float* Image = new float[4 * (width * height)];
    for (int i = 0; i < 4 * (width * height); i += 4)
    {
        float pixel_buffer[4];
        GET_RGBA32(data, i, pixel_buffer);
        SET_RGBA32(Image, i, pixel_buffer);
    }

    int old_pixel[4], new_pixel[4], quant_error[4];

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            GET_RGBA32(Image, POS_XY(4, width, i, j), old_pixel);
            find_closest_color(old_pixel, new_pixel);
            SET_RGBA32(data, POS_XY(4, width, i, j), new_pixel);
            SET_RGBA32(quant_error, 0, RGBA32_BLACK);
            SET_RGBA32(quant_error, 0, old_pixel);
            DIFF_RGBA32_NAIVE(quant_error, 0, new_pixel);

            if (!(j == width - 1))                   { ADD_RGBA32_NAIVE(Image, POS_XY(4, width, i, j + 1),     quant_error, (7.0 / 16)); }
            if (!(i == height - 1 || j == 0))        { ADD_RGBA32_NAIVE(Image, POS_XY(4, width, i + 1, j - 1), quant_error, (3.0 / 16)); }
            if (!(i == height -1))                   { ADD_RGBA32_NAIVE(Image, POS_XY(4, width, i + 1, j),     quant_error, (5.0 / 16)); }
            if (!(i == height -1 || j == width - 1)) { ADD_RGBA32_NAIVE(Image, POS_XY(4, width, i + 1, j + 1), quant_error, (1.0 / 16)); }
        }
    }

    delete[] Image;

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
    // constants
    const float scale = 1.0 / 25;
    const int kernel[][5] = {{ 1, 1, 1, 1, 1 },
                             { 1, 1, 1, 1, 1 },
                             { 1, 1, 1, 1, 1 },
                             { 1, 1, 1, 1, 1 },
                             { 1, 1, 1, 1, 1 }};

    // new Image
    float* Image = new float[4 * (width * height)];

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            float pixel_buffer[4];
            GET_2D_CONVOLUTION(data, 4, height, width, kernel, 5, 5, pixel_buffer, i-2, j-2, scale);
            SET_RGBA32(Image, POS_XY(4, width, i, j), pixel_buffer);
        }
    }

    // update
    for (int i = 0; i < 4 * (width * height); i += 4)
    {
        float pixel_buffer[4];
        GET_RGBA32(Image, i, pixel_buffer);
        SET_RGBA32(data, i, pixel_buffer);
    }
    delete[] Image;

    return true;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    // constants
    const float scale = 1.0 / 81;
    const int kernel[][5] = {{ 1, 2, 3, 2, 1 },
                             { 2, 4, 6, 4, 2 },
                             { 3, 6, 9, 6, 3 },
                             { 2, 4, 6, 4, 2 },
                             { 1, 2, 3, 2, 1 }};

    // new Image
    float* Image = new float[4 * (width * height)];

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            float pixel_buffer[4];
            GET_2D_CONVOLUTION(data, 4, height, width, kernel, 5, 5, pixel_buffer, i-2, j-2, scale);
            SET_RGBA32(Image, POS_XY(4, width, i, j), pixel_buffer);
        }
    }

    // update
    for (int i = 0; i < 4 * (width * height); i += 4)
    {
        float pixel_buffer[4];
        GET_RGBA32(Image, i, pixel_buffer);
        SET_RGBA32(data, i, pixel_buffer);
    }
    delete[] Image;

    return true;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    return Filter_Gaussian_N(5);
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N( unsigned int N )
{
    // generate kernel
    double scale = 1.0;
    double **kernel = new double*[N];
    for (int i = 0; i < N; ++i)
    {
        kernel[i] = new double[N];
        for (int j = 0; j < N; ++j)
        {
            kernel[i][j] = Binomial(N - 1, i) / (1 << (N - 1)) * \
                           Binomial(N - 1, j) / (1 << (N - 1));
        }
    }

    // new Image
    int* Image = new int[4 * (width * height)];

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            double pixel_buffer[4];
            GET_2D_CONVOLUTION(data, 4, height, width, kernel, N, N, pixel_buffer, i - ((int) (N / 2)), j - ((int) (N / 2)), scale);
            SET_RGBA32(Image, POS_XY(4, width, i, j), pixel_buffer);
        }
    }

    // update
    for (int i = 0; i < 4 * (width * height); i += 4)
    {
        double pixel_buffer[4];
        GET_RGBA32(Image, i, pixel_buffer);
        SET_RGBA32(data, i, pixel_buffer);
    }
    delete[] Image;

    for (int i = 0; i < N; ++i)
    {
        delete[] kernel[i];
    }
    delete[] kernel;

    return true;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    // constants
    const float scale = 1.0 / 81;
    const int kernel[][5] = {{ -1, -2, -3, -2, -1 },
                             { -2, -4, -6, -4, -2 },
                             { -3, -6, 72, -6, -3 },
                             { -2, -4, -6, -4, -2 },
                             { -1, -2, -3, -2, -1 }};

    // new Image
    double* Image = new double[4 * (width * height)];

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            double pixel_buffer[4];
            GET_2D_CONVOLUTION(data, 4, height, width, kernel, 5, 5, pixel_buffer, i-2, j-2, scale);
            SET_RGBA32(Image, POS_XY(4, width, i, j), pixel_buffer);
        }
    }

    // update
    for (int i = 0; i < 4 * (width * height); i += 4)
    {
        double pixel_buffer[4];
        GET_RGBA32(Image, i, pixel_buffer);
        SET_RGBA32_CUT(data, i, pixel_buffer);
    }
    delete[] Image;

    return true;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    // constants
    const float scale = 1.0 / 81;
    const int kernel[][5] = {{ -1, -2,  -3, -2, -1 },
                             { -2, -4,  -6, -4, -2 },
                             { -3, -6, 153, -6, -3 },
                             { -2, -4,  -6, -4, -2 },
                             { -1, -2,  -3, -2, -1 }};

    // new Image
    double* Image = new double[4 * (width * height)];

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            double pixel_buffer[4];
            GET_2D_CONVOLUTION(data, 4, height, width, kernel, 5, 5, pixel_buffer, i-2, j-2, scale);
            SET_RGBA32(Image, POS_XY(4, width, i, j), pixel_buffer);
        }
    }

    // update
    for (int i = 0; i < 4 * (width * height); i += 4)
    {
        double pixel_buffer[4];
        GET_RGBA32(Image, i, pixel_buffer);
        SET_RGBA32_CUT(data, i, pixel_buffer);
    }
    delete[] Image;

    return true;
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
    // constants
    const int Brush_Size[] = { 7, 3, 1 };

    // copied of original image
    int* Original = new int[4 * (width * height)];
    for (int i = 0; i < 4 * (width * height); i += 4)
    {
        int pixel_buffer[4];
        GET_RGBA32(data, i, pixel_buffer);
        SET_RGBA32(Original, i, pixel_buffer);
    }

    // new canvas
    ClearToBlack();

    // brush size = ...
    for (int b = 0; b < 3; ++b)
    {
        // using box filter to get avg color in [-b, -b] to [b, b]
        int box_filter_size = Brush_Size[b];
        double box_filter_scale = 1 / box_filter_size;
        double** box_filter_kernal = new double* [box_filter_size];
        for (int i = 0; i < box_filter_size; ++i)
        {
            box_filter_kernal[i] = new double [box_filter_size];
            for (int j = 0; j < box_filter_size; ++j)
            {
                box_filter_kernal[i][j] = 1 / box_filter_size;
            }
        }

        // using gauss filter to get brush color in [-r, -r] to [r, r]
        int gaussian_filter_size = (2 * Brush_Size[b] + 1);
        double gaussian_filter_scale = 1.0;
        double** gaussian_filter_kernal = new double* [gaussian_filter_size];
        for (int i = 0; i < gaussian_filter_size; ++i)
        {
            gaussian_filter_kernal[i] = new double [gaussian_filter_size];
            for (int j = 0; j < gaussian_filter_size; ++j)
            {
                gaussian_filter_kernal[i][j] = (Binomial((gaussian_filter_size - 1), i) / (1 << (gaussian_filter_size - 1))) * \
                                               (Binomial((gaussian_filter_size - 1), j) / (1 << (gaussian_filter_size - 1)));
            }
        }

        vector<Stroke> stroke_pool;

        for (int i = 0; i < height; i += Brush_Size[b])
        {
            for (int j = 0; j < width; j += Brush_Size[b])
            {
                double current_color[4], brush_color[4];
                GET_2D_CONVOLUTION(Original, 4, height, width, \
                                   gaussian_filter_kernal, gaussian_filter_size, gaussian_filter_size, \
                                   brush_color, ((int) i - Brush_Size[b]), ((int) j - Brush_Size[b]), 1.0);

                GET_2D_CONVOLUTION(data, 4, height, width, \
                                   box_filter_kernal, box_filter_size, box_filter_size, \
                                   current_color, ((int) i - ((int) Brush_Size[b] / 2)), ((int) j - ((int) Brush_Size[b] / 2)), box_filter_scale);

                // first layer or color diff > 100
                if (!b || abs(brush_color[0] + brush_color[1] + brush_color[2] - current_color[0] - current_color[1] - current_color[2]) > 125)
                {
                    int x = j, y = i, large_error = 0;
                    for (int di = -((int) Brush_Size[b] / 2); b && di < ((int) Brush_Size[b] / 2); ++di)
                    {
                        for (int dj = -((int) Brush_Size[b] / 2); b && dj < ((int) Brush_Size[b] / 2); ++dj)
                        {
                            int pixel_buffer[4];
                            GET_RGBA32_SAFE(Original, 4, height, width, i + di, j + dj, pixel_buffer);
                            if (__GET_RGBA32_SAFE)
                            {
                                int pixel_buffer_error = abs((pixel_buffer[ RED ] + pixel_buffer[ GREEN ] + pixel_buffer[ BLUE ]) - \
                                                             (current_color[ RED ] + current_color[ GREEN ] + current_color[ BLUE ]));
                                if (pixel_buffer_error > large_error)
                                {
                                    large_error = pixel_buffer_error;
                                    y = i + di;
                                    x = j + dj;
                                }
                            }
                        }
                    }
                    stroke_pool.push_back(Stroke(Brush_Size[b], x, y, brush_color[ RED ], brush_color[ GREEN ], brush_color[ BLUE ], brush_color[ ALPHA ]));
                }
            }
        }

        random_shuffle(begin(stroke_pool), stroke_pool.end());

        for (auto i = stroke_pool.begin(); i != stroke_pool.end(); i++)
        {
            Paint_Stroke(*i);
        }

        for (int i = 0; i < box_filter_size; ++i)
        {
            delete[] box_filter_kernal[i];
        }
        delete[] box_filter_kernal;

        for (int i = 0; i < gaussian_filter_size; ++i)
        {
            delete[] gaussian_filter_kernal[i];
        }
        delete[] gaussian_filter_kernal;
    }

    delete[] Original;

    return true;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    return Resize(0.5);
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    return Resize(2.0);
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    // constants
    const float scale_00 = 1.0 / 16;
    const int kernel_00[][3] = {{ 1, 2, 1 },
                                { 2, 4, 2 },
                                { 1, 2, 1 }};
    const float scale_01 = 1.0 / 32;
    const int kernel_01[][3] = {{ 1, 2, 1 },
                                { 3, 6, 3 },
                                { 3, 6, 3 },
                                { 1, 2, 1 }};
    const float scale_10 = 1.0 / 32;
    const int kernel_10[][4] = {{ 1, 3, 3, 1 },
                                { 2, 6, 6, 2 },
                                { 1, 3, 3, 1 }};
    const float scale_11 = 1.0 / 64;
    const int kernel_11[][4] = {{ 1, 3, 3, 1 },
                                { 3, 9, 9, 3 },
                                { 3, 9, 9, 3 },
                                { 1, 3, 3, 1 }};

    // new Image
    int* Image = new int[4 * (((int)(width * scale)) * ((int)(height * scale)))];

    for (int i = 0; i < ((int)(height * scale)); ++i)
    {
        for (int j = 0; j < ((int)(width * scale)); ++j)
        {
            int pixel_buffer[4];
            if (!(i % 2) && !(j % 2))
            {
                GET_2D_CONVOLUTION(data, 4, height, width, kernel_00, 3, 3, pixel_buffer, ((int)(i / scale))-1, ((int)(j / scale))-1, scale_00);
            }
            else if ((i % 2) && !(j % 2))
            {
                GET_2D_CONVOLUTION(data, 4, height, width, kernel_01, 4, 3, pixel_buffer, ((int)(i / scale))-1, ((int)(j / scale))-1, scale_01);
            }
            else if (!(i % 2) && (j % 2))
            {
                GET_2D_CONVOLUTION(data, 4, height, width, kernel_10, 3, 4, pixel_buffer, ((int)(i / scale))-1, ((int)(j / scale))-1, scale_10);
            }
            else
            {
                GET_2D_CONVOLUTION(data, 4, height, width, kernel_11, 4, 4, pixel_buffer, ((int)(i / scale))-1, ((int)(j / scale))-1, scale_11);
            }
            SET_RGBA32(Image, POS_XY(4, ((int)(width * scale)), i, j), pixel_buffer);
        }
    }

    // update
    free(this->data);
    this->data = new unsigned char[4 * (((int)(width * scale)) * ((int)(height * scale)))];
    for (int i = 0; i < 4 * (((int)(width * scale)) * ((int)(height * scale))); i += 4)
    {
        int pixel_buffer[4];
        GET_RGBA32(Image, i, pixel_buffer);
        SET_RGBA32_CUT(data, i, pixel_buffer);
    }
    width = ((int)(width * scale));
    height = ((int)(height * scale));
    delete[] Image;

    return true;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    // constants
    const float radians = DegreesToRadians(angleDegrees);

    // new Image
    float* Image = new float[4 * (width * height)];

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            int x = (j - (width / 2)) *  cos(radians) + (i - (height / 2)) * sin(radians) + (width / 2);
            int y = (j - (width / 2)) * -sin(radians) + (i - (height / 2)) * cos(radians) + (height / 2);
            float pixel_buffer[4];
            GET_RGBA32_SAFE(data, 4, height, width, y, x, pixel_buffer);
            SET_RGBA32(Image, POS_XY(4, width, i, j), pixel_buffer);
        }
    }

    // update
    for (int i = 0; i < 4 * (width * height); i += 4)
    {
        float pixel_buffer[4];
        GET_RGBA32(Image, i, pixel_buffer);
        SET_RGBA32_CUT(data, i, pixel_buffer);
    }
    delete[] Image;

    return true;
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


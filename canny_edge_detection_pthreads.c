#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "time.h"

#define MAX_PIXEL_VALUE 255
#define WEAK_PIXEL_VALUE 25
#define MIN_PIXEL_VALUE 0

#define GAUSS_KERNEL_SIZE 5
#define GAUSS_KERNEL_SIGMA 1.4
#define SOBEL_KERNEL_SIZE 3
#define HALF_CIRCLE_ANGLE 180.0
#define LOW_THRESHOLD_RATIO 0.05
#define HIGH_THRESHOLD_RATIO 0.09

#define ERROR_CODE -1

#pragma pack(1)

typedef struct
{
    unsigned char fileMarker1;    /* 'B' */
    unsigned char fileMarker2;    /* 'M' */
    unsigned int bfSize;          /* File's size */
    unsigned short unused1;       /* Aplication specific */
    unsigned short unused2;       /* Aplication specific */
    unsigned int imageDataOffset; /* Offset to the start of image data */
} bmp_fileheader;

typedef struct
{
    unsigned int biSize; /* Size of the info header - 40 bytes */
    signed int width;    /* Width of the image */
    signed int height;   /* Height of the image */
    unsigned short planes;
    unsigned short bitPix;      /* Number of bits per pixel = 3 * 8 */
    unsigned int biCompression; /* Type of compression */
    unsigned int biSizeImage;   /* Size of the image data */
    int biXPelsPerMeter;
    int biYPelsPerMeter;
    unsigned int biClrUsed;
    unsigned int biClrImportant;
    unsigned int redMask;
    unsigned int greenMask;
    unsigned int blueMask;
    unsigned int alphaMask;
    unsigned int CSType;
    int redX;
    int redY;
    int redZ;
    int greenX;
    int greenY;
    int greenZ;
    int blueX;
    int blueY;
    int blueZ;
    unsigned int gammaRed;
    unsigned int gammaGreen;
    unsigned int gammaBlue;
} bmp_infoheader;

typedef struct
{
    double Green;
    double Red;
    double Blue;
} bmp_image;

typedef struct
{
    int num_threads;
    int thread_id;
} thread_arg;

pthread_barrier_t barrier;
pthread_mutex_t mutex;

#pragma pack()

bmp_fileheader *fileHeader;
bmp_infoheader *infoHeader;
bmp_image **physicalImage;
bmp_image **theta;
bmp_image max_pixel1;
bmp_image max_pixel2;
double** gauss_kernel;

// Function to read the file header and the info header.
void read_data(FILE *image, bmp_fileheader *fileHeader, bmp_infoheader *infoHeader)
{
    fread(&fileHeader->fileMarker1, sizeof(unsigned char), 1, image);
    fread(&fileHeader->fileMarker2, sizeof(unsigned char), 1, image);
    fread(&fileHeader->bfSize, sizeof(unsigned int), 1, image);
    fread(&fileHeader->unused1, sizeof(unsigned short), 1, image);
    fread(&fileHeader->unused2, sizeof(unsigned short), 1, image);
    fread(&fileHeader->imageDataOffset, sizeof(unsigned int), 1, image);

    fread(&infoHeader->biSize, sizeof(unsigned int), 1, image);
    fread(&infoHeader->width, sizeof(signed int), 1, image);
    fread(&infoHeader->height, sizeof(signed int), 1, image);
    fread(&infoHeader->planes, sizeof(unsigned short), 1, image);
    fread(&infoHeader->bitPix, sizeof(unsigned short), 1, image);
    fread(&infoHeader->biCompression, sizeof(unsigned int), 1, image);
    fread(&infoHeader->biSizeImage, sizeof(unsigned int), 1, image);
    fread(&infoHeader->biXPelsPerMeter, sizeof(int), 1, image);
    fread(&infoHeader->biYPelsPerMeter, sizeof(int), 1, image);
    fread(&infoHeader->biClrUsed, sizeof(unsigned int), 1, image);
    fread(&infoHeader->biClrImportant, sizeof(unsigned int), 1, image);

    fread(&infoHeader->redMask, sizeof(unsigned int), 1, image);
    fread(&infoHeader->greenMask, sizeof(unsigned int), 1, image);
    fread(&infoHeader->blueMask, sizeof(unsigned int), 1, image);
    fread(&infoHeader->alphaMask, sizeof(unsigned int), 1, image);
    fread(&infoHeader->CSType, sizeof(unsigned int), 1, image);
    fread(&infoHeader->redX, sizeof(int), 1, image);
    fread(&infoHeader->redY, sizeof(int), 1, image);
    fread(&infoHeader->redZ, sizeof(int), 1, image);
    fread(&infoHeader->greenX, sizeof(int), 1, image);
    fread(&infoHeader->greenY, sizeof(int), 1, image);
    fread(&infoHeader->greenZ, sizeof(int), 1, image);
    fread(&infoHeader->blueX, sizeof(int), 1, image);
    fread(&infoHeader->blueY, sizeof(int), 1, image);
    fread(&infoHeader->blueZ, sizeof(int), 1, image);
    fread(&infoHeader->gammaRed, sizeof(unsigned int), 1, image);
    fread(&infoHeader->gammaGreen, sizeof(unsigned int), 1, image);
    fread(&infoHeader->gammaBlue, sizeof(unsigned int), 1, image);
}

// Function to read the bmp data (the image).
void read_physicalImage(bmp_image **physicalImage, FILE *image, bmp_infoheader *infoHeader,
                        bmp_fileheader *fileHeader)
{
    unsigned char blue;
    unsigned char green;
    unsigned char red;

    fseek(image, fileHeader->imageDataOffset, SEEK_SET);
    for (int i = infoHeader->height - 1; i >= 0; i--)
    {
        for (int j = 0; j < infoHeader->width; j++)
        {
            fread(&blue, sizeof(unsigned char), 1, image);
            fread(&green, sizeof(unsigned char), 1, image);
            fread(&red, sizeof(unsigned char), 1, image);

            physicalImage[i][j].Blue = blue;
            physicalImage[i][j].Green = green;
            physicalImage[i][j].Red = red;
        }

        if (infoHeader->width % 4 != 0)
        {
            fseek(image, infoHeader->width % 4, SEEK_CUR);
        }
    }
}

// Function to write the bmp file header and the bmp info header.
void print(FILE *out, bmp_fileheader *fileHeader, bmp_infoheader *infoHeader)
{
    fwrite(&fileHeader->fileMarker1, sizeof(unsigned char), 1, out);
    fwrite(&fileHeader->fileMarker2, sizeof(unsigned char), 1, out);
    fwrite(&fileHeader->bfSize, sizeof(unsigned int), 1, out);
    fwrite(&fileHeader->unused1, sizeof(unsigned short), 1, out);
    fwrite(&fileHeader->unused2, sizeof(unsigned short), 1, out);
    fwrite(&fileHeader->imageDataOffset, sizeof(unsigned int), 1, out);

    fwrite(&infoHeader->biSize, sizeof(unsigned int), 1, out);
    fwrite(&infoHeader->width, sizeof(signed int), 1, out);
    fwrite(&infoHeader->height, sizeof(signed int), 1, out);
    fwrite(&infoHeader->planes, sizeof(unsigned short), 1, out);
    fwrite(&infoHeader->bitPix, sizeof(unsigned short), 1, out);
    fwrite(&infoHeader->biCompression, sizeof(unsigned int), 1, out);
    fwrite(&infoHeader->biSizeImage, sizeof(unsigned int), 1, out);
    fwrite(&infoHeader->biXPelsPerMeter, sizeof(int), 1, out);
    fwrite(&infoHeader->biYPelsPerMeter, sizeof(int), 1, out);
    fwrite(&infoHeader->biClrUsed, sizeof(unsigned int), 1, out);
    fwrite(&infoHeader->biClrImportant, sizeof(unsigned int), 1, out);

    fwrite(&infoHeader->redMask, sizeof(unsigned int), 1, out);
    fwrite(&infoHeader->greenMask, sizeof(unsigned int), 1, out);
    fwrite(&infoHeader->blueMask, sizeof(unsigned int), 1, out);
    fwrite(&infoHeader->alphaMask, sizeof(unsigned int), 1, out);
    fwrite(&infoHeader->CSType, sizeof(unsigned int), 1, out);
    fwrite(&infoHeader->redX, sizeof(int), 1, out);
    fwrite(&infoHeader->redY, sizeof(int), 1, out);
    fwrite(&infoHeader->redZ, sizeof(int), 1, out);
    fwrite(&infoHeader->greenX, sizeof(int), 1, out);
    fwrite(&infoHeader->greenY, sizeof(int), 1, out);
    fwrite(&infoHeader->greenZ, sizeof(int), 1, out);
    fwrite(&infoHeader->blueX, sizeof(int), 1, out);
    fwrite(&infoHeader->blueY, sizeof(int), 1, out);
    fwrite(&infoHeader->blueZ, sizeof(int), 1, out);
    fwrite(&infoHeader->gammaRed, sizeof(unsigned int), 1, out);
    fwrite(&infoHeader->gammaGreen, sizeof(unsigned int), 1, out);
    fwrite(&infoHeader->gammaBlue, sizeof(unsigned int), 1, out);
}

// Function to write the bmp data (the image).
void print_physicalImage(bmp_image **physicalImage, FILE *out, bmp_fileheader *fileHeader,
    bmp_infoheader *infoHeader) {
    int padding = 0;
    unsigned char blue, green, red;

    fseek(out, fileHeader->imageDataOffset, SEEK_SET);
    for(int i = infoHeader->height - 1; i >= 0; i--){
        for(int j = 0; j < infoHeader->width; j++){
            blue = physicalImage[i][j].Blue;
            green = physicalImage[i][j].Green;
            red = physicalImage[i][j].Red;

            if (i <= 3 || i >= infoHeader->height - 4 || j <= 2 || j >= infoHeader->width - 2) {
                unsigned char zero = 0;
                fwrite(&zero, sizeof(unsigned char), 1, out);
                fwrite(&zero, sizeof(unsigned char), 1, out);
                fwrite(&zero, sizeof(unsigned char), 1, out);
            } else {
                fwrite(&blue, sizeof(unsigned char), 1, out);
                fwrite(&green, sizeof(unsigned char), 1, out);
                fwrite(&red, sizeof(unsigned char), 1, out);
            }
        }

        fwrite(&padding, sizeof(unsigned char), infoHeader->width % 4, out);
    }
}

// Function to convert a color photo to black and white.
void rgb2gray(bmp_image ***physicalImage, bmp_infoheader *infoHeader, thread_arg arg)
{
    double noColor;

    int start, end;
    start = arg.thread_id * (double)infoHeader->height / arg.num_threads;
    end =
        (arg.thread_id + 1) * (double)infoHeader->height / arg.num_threads < infoHeader->height
            ? (arg.thread_id + 1) * (double)infoHeader->height / arg.num_threads
            : infoHeader->height;

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < infoHeader->width; j++)
        {
            noColor = ((*physicalImage)[i][j].Blue * 0.114021 +
                       (*physicalImage)[i][j].Green * 0.587043 +
                       (*physicalImage)[i][j].Red * 0.298936) /
                      3;

            (*physicalImage)[i][j].Blue = noColor;
            (*physicalImage)[i][j].Green = noColor;
            (*physicalImage)[i][j].Red = noColor;
        }
    }
}

// Function to compute the gaussian kernel.
double **compute_gaussian_kernel(int size, double sigma)
{
    int k = size / 2;
    double normal = 1.0 / (2 * M_PI * sigma * sigma);

    double **x = (double **)calloc(size, sizeof(double *));
    for (int i = 0; i < size; i++)
    {
        x[i] = (double *)calloc(size, sizeof(double));
    }
    double **y = (double **)calloc(size, sizeof(double *));
    for (int i = 0; i < size; i++)
    {
        y[i] = (double *)calloc(size, sizeof(double));
    }
    double **result = (double **)calloc(size, sizeof(double *));
    for (int i = 0; i < size; i++)
    {
        result[i] = (double *)calloc(size, sizeof(double));
    }

    for (int i = -k; i < k + 1; i++)
    {
        for (int j = 0; j < size; j++)
        {
            x[i + k][j] = i;
            y[j][i + k] = i;
        }
    }

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            result[i][j] = exp(-(x[i][j] * x[i][j] + y[i][j] * y[i][j]) / (2.0 * sigma * sigma)) * normal;
        }
    }

    return result;
}

void free_gaussian_kernel(double*** gauss_kernel, int height) {
    for (int i = 0; i < height; i++) {
        free((*gauss_kernel)[i]);
    }
    free(*gauss_kernel);
    *gauss_kernel = NULL;
}

void free_image(bmp_image ***image, int height) {
    for (int i = 0; i < height; i++) {
        free((*image)[i]);
    }
    free(*image);
    *image = NULL;
}

// Function to blur the image by applying the gaussian kernel on the image.
void image_noise_reduction(bmp_image ***physicalImage, bmp_infoheader *infoHeader,
                           double **gauss_kernel, int k, thread_arg arg)
{
    bmp_image **img = *physicalImage;

    int start, end;
    start = arg.thread_id * (double)infoHeader->height / arg.num_threads;
    end =
        (arg.thread_id + 1) * (double)infoHeader->height / arg.num_threads < infoHeader->height
            ? (arg.thread_id + 1) * (double)infoHeader->height / arg.num_threads
            : infoHeader->height;

    bmp_image **result = (bmp_image **)calloc(end - start, sizeof(bmp_image *));
    for (int i = 0; i < end - start; i++)
    {
        result[i] = (bmp_image *)calloc(infoHeader->width, sizeof(bmp_image));
    }

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < infoHeader->width; j++)
        {
            double new_fade = 0;
            for (int m = -k; m < k + 1; m++)
            {
                for (int n = -k; n < k + 1; n++)
                {
                    if (i + m >= 0 && j + n >= 0 && i + m < infoHeader->height && j + n < infoHeader->width)
                    {
                        new_fade += img[i + m][j + n].Blue * gauss_kernel[m + k][n + k];
                    }
                }
            }

            result[i - start][j].Blue = new_fade;
            result[i - start][j].Green = new_fade;
            result[i - start][j].Red = new_fade;
        }
    }

    pthread_barrier_wait(&barrier);

    for (int i = start; i < end; i++)
    {
        (*physicalImage)[i] = result[i - start];
    }
}

// Function to detects the edge intensity and direction.
void sobel_filters(bmp_image ***physicalImage, bmp_infoheader *infoHeader, thread_arg arg)
{
    int Kx[SOBEL_KERNEL_SIZE][SOBEL_KERNEL_SIZE] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    int Ky[SOBEL_KERNEL_SIZE][SOBEL_KERNEL_SIZE] = {{1, 2, 1}, {0, 0, 0}, {-1, -2, -1}};

    int start, end;
    start = arg.thread_id * (double)infoHeader->height / arg.num_threads;
    end =
        (arg.thread_id + 1) * (double)infoHeader->height / arg.num_threads < infoHeader->height
            ? (arg.thread_id + 1) * (double)infoHeader->height / arg.num_threads
            : infoHeader->height;

    bmp_image **img = *physicalImage;
    bmp_image **result = (bmp_image **)calloc(end - start, sizeof(bmp_image *));

    for (int i = 0; i < end - start; i++)
    {
        result[i] = (bmp_image *)calloc(infoHeader->width, sizeof(bmp_image));
    }

    int k = SOBEL_KERNEL_SIZE / 2;
    for (int i = start; i < end; i++) {
        for(int j = 0; j < infoHeader->width; j++) {
            double new_fade_x = 0;
            double new_fade_y = 0;
            for (int m = -k; m < k + 1; m++) {
                for (int n = -k; n < k + 1; n++) {
                    if (j + n >= 0 && j + n < infoHeader->width && i + m >= 0 && i + m < infoHeader->height) {
                        new_fade_x += (double)img[i + m][j + n].Blue * Kx[m + k][n + k];
                        new_fade_y += (double)img[i + m][j + n].Blue * Ky[m + k][n + k];
                    }
                }
            }

            double new_val = sqrt(new_fade_x * new_fade_x + new_fade_y * new_fade_y);
            result[i - start][j].Blue = new_val;
            result[i - start][j].Green = new_val;
            result[i - start][j].Red = new_val;

            theta[i][j].Blue = atan2(new_fade_x, new_fade_y);
            theta[i][j].Green = atan2(new_fade_x, new_fade_y);
            theta[i][j].Red = atan2(new_fade_x, new_fade_y);
        }
    }

    pthread_barrier_wait(&barrier);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < infoHeader->width; j++) {
            pthread_mutex_lock(&mutex);
            if (img[i][j].Red > max_pixel1.Red && img[i][j].Green > max_pixel1.Green 
                && img[i][j].Blue > max_pixel1.Blue) {
                max_pixel1.Red = img[i][j].Red;
                max_pixel1.Blue = img[i][j].Blue;
                max_pixel1.Green = img[i][j].Green;
            }
            pthread_mutex_unlock(&mutex);
        }
    }

    pthread_barrier_wait(&barrier);

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < infoHeader->width; j++)
        {
            result[i - start][j].Blue = (double)result[i - start][j].Blue / max_pixel1.Blue * 255;
            result[i - start][j].Green = (double)result[i - start][j].Green / max_pixel1.Green * 255;
            result[i - start][j].Red = (double)result[i - start][j].Red / max_pixel1.Red * 255;
        }
    }

    pthread_barrier_wait(&barrier);

    for (int i = start; i < end; i++)
    {
        (*physicalImage)[i] = result[i - start];
    }
}

// Function to perform non-maximum suppression to thin out the edges.
void non_max_suppression(bmp_image ***physicalImage, bmp_image **theta_img,
                         bmp_infoheader *infoHeader, thread_arg arg)
{
    bmp_image **img = *physicalImage;

    int start, end;
    start = arg.thread_id * (double)infoHeader->height / arg.num_threads;
    end =
        (arg.thread_id + 1) * (double)infoHeader->height / arg.num_threads < infoHeader->height
            ? (arg.thread_id + 1) * (double)infoHeader->height / arg.num_threads
            : infoHeader->height;

    bmp_image **result = (bmp_image **)calloc(end - start, sizeof(bmp_image *));
    for (int i = 0; i < end - start; i++)
    {
        result[i] = (bmp_image *)calloc(infoHeader->width, sizeof(bmp_image));
    }

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < infoHeader->width; j++)
        {
            theta_img[i][j].Blue = theta_img[i][j].Blue * HALF_CIRCLE_ANGLE / M_PI;
            if (theta_img[i][j].Blue < 0)
            {
                theta_img[i][j].Blue += HALF_CIRCLE_ANGLE;
            }
        }
    }

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < infoHeader->width; j++)
        {
            double q = 255;
            double r = 255;

            if (((0 <= theta_img[i][j].Blue < 22.5) || (157.5 <= theta_img[i][j].Blue <= 180)) 
                && (j + 1 < infoHeader->width) && (j-1>=0))
            {
                q = img[i][j + 1].Blue;
                r = img[i][j - 1].Blue;
            }
            else if (22.5 <= theta_img[i][j].Blue < 67.5 && (i+1<infoHeader->height) 
                && (i-1>=0) && (j-1>=0) && (j+1<infoHeader->width)) 
            {
                q = img[i + 1][j - 1].Blue;
                r = img[i - 1][j + 1].Blue;
            }
            else if (67.5 <= theta_img[i][j].Blue < 112.5 && (i+1<infoHeader->height) 
                && (i-1>=0)) 
            {
                q = img[i + 1][j].Blue;
                r = img[i - 1][j].Blue;
            }
            else if (112.5 <= theta_img[i][j].Blue < 157.5 && (i+1<infoHeader->height) 
                && (i-1>=0) && (j-1>=0) && (j+1<infoHeader->width))
            {
                q = img[i - 1][j - 1].Blue;
                r = img[i + 1][j + 1].Blue;
            }

            if ((img[i][j].Blue >= q) && (img[i][j].Blue >= r))
            {
                result[i - start][j].Blue = img[i][j].Blue;
                result[i - start][j].Green = img[i][j].Blue;
                result[i - start][j].Red = img[i][j].Blue;
            }
            else
            {
                result[i - start][j].Blue = 0;
                result[i - start][j].Green = 0;
                result[i - start][j].Red = 0;
            }
        }
    }

    for (int i = start; i < end; i++)
    {
        (*physicalImage)[i] = result[i - start];
    }
}

// Function to identifying 3 kinds of pixels: strong, weak, and non-relevant.
void double_threshold(bmp_image ***physicalImage, bmp_infoheader *infoHeader,
                      double low_threshold_ratio, double high_threshold_ratio, thread_arg arg)
{
    bmp_image **img = *physicalImage;

    int start, end;
    start = arg.thread_id * (double)infoHeader->height / arg.num_threads;
    end =
        (arg.thread_id + 1) * (double)infoHeader->height / arg.num_threads < infoHeader->height
            ? (arg.thread_id + 1) * (double)infoHeader->height / arg.num_threads
            : infoHeader->height;

    bmp_image **result = (bmp_image **)calloc(end - start, sizeof(bmp_image *));
    for (int i = 0; i < end - start; i++)
    {
        result[i] = (bmp_image *)calloc(infoHeader->width, sizeof(bmp_image));
    }

    pthread_barrier_wait(&barrier);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < infoHeader->width; j++) {
            pthread_mutex_lock(&mutex);
            if (img[i][j].Red > max_pixel2.Red && img[i][j].Green > max_pixel2.Green 
                && img[i][j].Blue > max_pixel2.Blue) {
                max_pixel2.Red = img[i][j].Red;
                max_pixel2.Blue = img[i][j].Blue;
                max_pixel2.Green = img[i][j].Green;
            }
            pthread_mutex_unlock(&mutex);
        }
    }

    pthread_barrier_wait(&barrier);

    bmp_image high_threshold;
    high_threshold.Red = max_pixel2.Red * high_threshold_ratio;
    high_threshold.Blue = max_pixel2.Blue * high_threshold_ratio;
    high_threshold.Green = max_pixel2.Green * high_threshold_ratio;

    bmp_image low_threshold;
    low_threshold.Red = high_threshold.Red * low_threshold_ratio;
    low_threshold.Blue = high_threshold.Blue * low_threshold_ratio;
    low_threshold.Green = high_threshold.Green * low_threshold_ratio;

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < infoHeader->width; j++)
        {
            if (img[i][j].Red >= high_threshold.Red && img[i][j].Green >= high_threshold.Green && img[i][j].Blue >= high_threshold.Blue)
            {
                result[i - start][j].Red = MAX_PIXEL_VALUE;
                result[i - start][j].Green = MAX_PIXEL_VALUE;
                result[i - start][j].Blue = MAX_PIXEL_VALUE;
            }
            else if (img[i][j].Red >= low_threshold.Red && img[i][j].Green >= low_threshold.Green && img[i][j].Blue >= low_threshold.Blue)
            {
                result[i - start][j].Red = WEAK_PIXEL_VALUE;
                result[i - start][j].Green = WEAK_PIXEL_VALUE;
                result[i - start][j].Blue = WEAK_PIXEL_VALUE;
            }
            else
            {
                result[i - start][j].Red = MIN_PIXEL_VALUE;
                result[i - start][j].Green = MIN_PIXEL_VALUE;
                result[i - start][j].Blue = MIN_PIXEL_VALUE;
            }
        }
    }

    pthread_barrier_wait(&barrier);

    for (int i = start; i < end; i++)
    {
        (*physicalImage)[i] = result[i - start];
    }
}

int is_strong_pixel(bmp_image pixel)
{
    if (pixel.Red == MAX_PIXEL_VALUE && pixel.Green == MAX_PIXEL_VALUE && pixel.Blue == MAX_PIXEL_VALUE)
    {
        return 1;
    }

    return 0;
}

// Function to transform weak pixels into strong ones.
void hysteresis(bmp_image ***physicalImage, bmp_infoheader *infoHeader, thread_arg arg)
{
    bmp_image **img = *physicalImage;

    int start, end;
    start = arg.thread_id * (double)infoHeader->height / arg.num_threads;
    end =
        (arg.thread_id + 1) * (double)infoHeader->height / arg.num_threads < infoHeader->height
            ? (arg.thread_id + 1) * (double)infoHeader->height / arg.num_threads
            : infoHeader->height;

    bmp_image **result = (bmp_image **)calloc(end - start, sizeof(bmp_image *));
    for (int i = 0; i < end - start; i++)
    {
        result[i] = (bmp_image *)calloc(infoHeader->width, sizeof(bmp_image));
    }

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < infoHeader->width; j++)
        {
            if (img[i][j].Red == WEAK_PIXEL_VALUE && img[i][j].Green == WEAK_PIXEL_VALUE && img[i][j].Blue == WEAK_PIXEL_VALUE)
            {
                if ((i - 1 >= 0 && j - 1 >= 0 && is_strong_pixel(img[i - 1][j - 1])) ||
                    (i - 1 >= 0 && is_strong_pixel(img[i - 1][j])) ||
                    (j - 1 >= 0 && is_strong_pixel(img[i][j - 1])) ||
                    (i + 1 < infoHeader->height && j + 1 < infoHeader->width && is_strong_pixel(img[i + 1][j + 1])) ||
                    (i + 1 < infoHeader->height && is_strong_pixel(img[i + 1][j])) ||
                    (j + 1 < infoHeader->width && is_strong_pixel(img[i][j + 1])) ||
                    (i + 1 < infoHeader->height && j - 1 >= 0 && is_strong_pixel(img[i + 1][j - 1])) ||
                    (i - 1 >= 0 && j + 1 < infoHeader->width && is_strong_pixel(img[i - 1][j + 1])))
                {
                    result[i - start][j].Red = MAX_PIXEL_VALUE;
                    result[i - start][j].Blue = MAX_PIXEL_VALUE;
                    result[i - start][j].Green = MAX_PIXEL_VALUE;
                }
                else
                {
                    result[i - start][j].Red = MIN_PIXEL_VALUE;
                    result[i - start][j].Blue = MIN_PIXEL_VALUE;
                    result[i - start][j].Green = MIN_PIXEL_VALUE;
                }
            }
            else
            {
                result[i - start][j].Red = img[i][j].Red;
                result[i - start][j].Blue = img[i][j].Blue;
                result[i - start][j].Green = img[i][j].Green;
            }
        }
    }

    pthread_barrier_wait(&barrier);

    for (int i = start; i < end; i++)
    {
        (*physicalImage)[i] = result[i - start];
    }
}

void *thread_func(void *args)
{
    thread_arg arg = *(thread_arg *)args;

    // Turn image into black and white.
    rgb2gray(&physicalImage, infoHeader, arg);
    pthread_barrier_wait(&barrier);

    // Blur image.
    image_noise_reduction(&physicalImage, infoHeader, gauss_kernel, GAUSS_KERNEL_SIZE / 2, arg);
    pthread_barrier_wait(&barrier);

    // Detect the edge intensity and direction.
    sobel_filters(&physicalImage, infoHeader, arg);
    pthread_barrier_wait(&barrier);

    // Thin out the edges.
    non_max_suppression(&physicalImage, theta, infoHeader, arg);
    pthread_barrier_wait(&barrier);

    // Identifying strong, weak, and non-relevant pixels.
    double_threshold(&physicalImage, infoHeader, LOW_THRESHOLD_RATIO, HIGH_THRESHOLD_RATIO, arg);
    pthread_barrier_wait(&barrier);

    // Edge tracking by hysteresis.
    hysteresis(&physicalImage, infoHeader, arg);

    return NULL;
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        printf("The program should be ran ./<program_name> <number_of_threads>\n");
        exit(-1);
    }

    int num_threads = atoi(argv[1]);
    FILE *inputImage = fopen("wp.bmp", "rb");

    // Read the image file header and info header.
    fileHeader = calloc(1, sizeof(bmp_fileheader));
    if (fileHeader == NULL) {
        exit(ERROR_CODE);
    }
    infoHeader = calloc(1, sizeof(bmp_infoheader));
    if (infoHeader == NULL) {
        exit(ERROR_CODE);
    }
    read_data(inputImage, fileHeader, infoHeader);

    // Read the image data.
    physicalImage = calloc(infoHeader->height, sizeof(bmp_image *));
    if (physicalImage == NULL) {
        exit(ERROR_CODE);
    }
    for (int i = 0; i < infoHeader->height; i++) {
        physicalImage[i] = calloc(infoHeader->width, sizeof(bmp_image));
        if (physicalImage[i] == NULL) {
            exit(ERROR_CODE);
        }
    }
    read_physicalImage(physicalImage, inputImage, infoHeader, fileHeader);

    fclose(inputImage);

    time_t start = clock();

    int r;
    long id;
    void *status;
    thread_arg arguments[num_threads];
    pthread_t threads[num_threads];

    // Alloc memory for theta.
    theta = calloc(infoHeader->height, sizeof(bmp_image *));
    if (theta == NULL){
        exit(ERROR_CODE);
    }
    for (int i = 0; i < infoHeader->height; i++) {
        theta[i] = calloc(infoHeader->width, sizeof(bmp_image));
    }

    // Initialize maximum pixels at 0.
    max_pixel1.Red = MIN_PIXEL_VALUE;
    max_pixel1.Green = MIN_PIXEL_VALUE;
    max_pixel1.Blue = MIN_PIXEL_VALUE;

    max_pixel2.Red = MIN_PIXEL_VALUE;
    max_pixel2.Green = MIN_PIXEL_VALUE;
    max_pixel2.Blue = MIN_PIXEL_VALUE;

    // Compute the gaussian kernel.
    gauss_kernel = compute_gaussian_kernel(GAUSS_KERNEL_SIZE, GAUSS_KERNEL_SIGMA);

    pthread_barrier_init(&barrier, NULL, num_threads);
    pthread_mutex_init(&mutex, NULL);

    // Create threads.
    for (int i = 0; i < num_threads; i++) {
        arguments[i].num_threads = num_threads;
        arguments[i].thread_id = i;
        r = pthread_create(&threads[i], NULL, thread_func, &arguments[i]);
    }

    // Close threads.
    for (int i = 0; i < num_threads; i++) {
        r = pthread_join(threads[i], NULL);
    }

    time_t end = clock();

    printf("Total time taken: %f\n", ((double)(end - start)  / CLOCKS_PER_SEC));
    
    // Write image.
    FILE *out = fopen("image_pthreads.bmp", "wb");
    print(out, fileHeader, infoHeader);
    print_physicalImage(physicalImage, out, fileHeader, infoHeader);
    fclose(out);

    free_image(&physicalImage, infoHeader->height);
    free_gaussian_kernel(&gauss_kernel, GAUSS_KERNEL_SIZE);
    free_image(&theta, infoHeader->height);
    free(fileHeader);
    free(infoHeader);

    return 0;
}

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

#define ROOT                    0

#define MAX_PIXEL_VALUE 		255
#define WEAK_PIXEL_VALUE 		25
#define MIN_PIXEL_VALUE 		0

#define PADDING_SIZE            2
#define GAUSS_KERNEL_SIZE 		5
#define GAUSS_KERNEL_SIGMA 		1.4
#define SOBEL_KERNEL_SIZE 		3
#define HALF_CIRCLE_ANGLE 		180.0
#define LOW_THRESHOLD_RATIO 	0.05
#define HIGH_THRESHOLD_RATIO 	0.09

#define ERROR_CODE 				-1

#pragma pack(1)

typedef struct {
    unsigned char  fileMarker1;		/* 'B' */
    unsigned char  fileMarker2; 	/* 'M' */
    unsigned int   bfSize; 			/* File's size */
    unsigned short unused1; 		/* Aplication specific */
    unsigned short unused2; 		/* Aplication specific */
    unsigned int   imageDataOffset; /* Offset to the start of image data */
} bmp_fileheader;

typedef struct {
    unsigned int   biSize; 			/* Size of the info header - 40 bytes */
    signed int     width; 			/* Width of the image */
    signed int     height; 			/* Height of the image */
    unsigned short planes;
    unsigned short bitPix;  		/* Number of bits per pixel = 3 * 8 */
    unsigned int   biCompression; 	/* Type of compression */
    unsigned int   biSizeImage; 	/* Size of the image data */
    int            biXPelsPerMeter;
    int            biYPelsPerMeter;
    unsigned int   biClrUsed;
    unsigned int   biClrImportant;
    unsigned int   redMask;
    unsigned int   greenMask;
    unsigned int   blueMask;
    unsigned int   alphaMask;
    unsigned int   CSType;
    int            redX;
    int            redY;
    int            redZ;
    int            greenX;
    int            greenY;
    int            greenZ;
    int            blueX;
    int            blueY;
    int            blueZ;
    unsigned int   gammaRed;
    unsigned int   gammaGreen;
    unsigned int   gammaBlue;
} bmp_infoheader;

typedef struct {
    double Green;
    double Red;
    double Blue;
} bmp_image;

#pragma pack()

const int num_items = 3;
int block_lengths[3] = {1, 1, 1};
MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
MPI_Datatype MPI_PIXEL_TYPE;
MPI_Aint offsets[3];


// Function to read the file header and the info header.
void read_data(FILE *image, bmp_fileheader *fileHeader, bmp_infoheader *infoHeader) {
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
	bmp_fileheader * fileHeader) {
	unsigned char blue;
    unsigned char green;
    unsigned char red;

    fseek(image, fileHeader->imageDataOffset, SEEK_SET);
    for(int i = infoHeader->height - 1; i >= 0; i--){
        for(int j = 0; j < infoHeader->width; j++){
            fread(&blue, sizeof(unsigned char), 1, image);
            fread(&green, sizeof(unsigned char), 1, image);
            fread(&red, sizeof(unsigned char), 1, image);

            physicalImage[i][j].Blue = blue;
            physicalImage[i][j].Green= green;
            physicalImage[i][j].Red = red;
		}

        if(infoHeader->width % 4 != 0) {
            fseek(image, infoHeader->width % 4, SEEK_CUR);
        }
    }
}

// Function to write the bmp file header and the bmp info header.
void print(FILE *out, bmp_fileheader * fileHeader, bmp_infoheader *infoHeader) {
    fwrite(&fileHeader->fileMarker1, sizeof(unsigned char), 1, out);
    fwrite(&fileHeader->fileMarker2, sizeof(unsigned char), 1, out);
    fwrite(&fileHeader->bfSize, sizeof(unsigned int), 1, out);
    fwrite(&fileHeader->unused1, sizeof(unsigned short), 1, out);
    fwrite(&fileHeader->unused2, sizeof(unsigned short), 1, out);
    fwrite(&fileHeader->imageDataOffset, sizeof(unsigned int), 1,out);

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
void rgb2gray(bmp_image ***physicalImage, int height, int width){
    double noColor;

    for(int i = 2; i < height - 2; i++) {
        for(int j = 0; j < width; j++) {
            noColor = ((*physicalImage)[i][j].Blue *  0.114021 + 
                (*physicalImage)[i][j].Green * 0.587043 +
                (*physicalImage)[i][j].Red * 0.298936) / 3;

            (*physicalImage)[i][j].Blue = noColor;
            (*physicalImage)[i][j].Green = noColor;
            (*physicalImage)[i][j].Red = noColor;
        }
    }
}

// Function to compute the gaussian kernel.
double** compute_gaussian_kernel(int size, double sigma) {
    int k = size / 2;
    double normal = 1.0 / (2 * M_PI * sigma * sigma);

    double** x = (double**)calloc(size, sizeof(double *));
    for (int i = 0; i < size; i++) {
        x[i] = (double*)calloc(size, sizeof(double));
    }
    double** y = (double**)calloc(size, sizeof(double *));
    for (int i = 0; i < size; i++) {
        y[i] = (double*)calloc(size, sizeof(double));
    }
    double** result = (double**)calloc(size, sizeof(double *));
    for (int i = 0; i < size; i++) {
        result[i] = (double*)calloc(size, sizeof(double));
    }

    for (int i = -k; i < k + 1; i++) {
        for (int j = 0; j < size; j++) {
            x[i + k][j] = i;
            y[j][i + k] = i;
        }
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            result[i][j] = exp(-(x[i][j] * x[i][j] + y[i][j] * y[i][j]) 
            	/ (2.0 * sigma * sigma)) * normal;
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
void image_noise_reduction(bmp_image ***physicalImage, int height, int width, 
	double** gauss_kernel, int k) {
    bmp_image **img = *physicalImage;

    bmp_image **result = (bmp_image **)calloc(height, sizeof(bmp_image *));
    for (int i = 0; i < height; i++) {
        result[i] = (bmp_image *)calloc(width, sizeof(bmp_image));
    }

    for (int i = 2; i < height - 2; i++) {
        for(int j = 0; j < width; j++) {
            double new_fade = 0;
            for (int m = -k; m < k + 1; m++) {
                for (int n = -k; n < k + 1; n++) {
                    if (j + n >= 0 && j + n < width)
                        new_fade += img[i + m][j + n].Blue * gauss_kernel[m + k][n + k];
                }
            }

            result[i][j].Blue = new_fade;
            result[i][j].Green = new_fade;
            result[i][j].Red = new_fade;
        }
    }

    free_image(physicalImage, height);
    *physicalImage = result;
}

// Function to detects the edge intensity and direction.
bmp_image** sobel_filters(bmp_image ***physicalImage, int height, int width) {
    int Kx[SOBEL_KERNEL_SIZE][SOBEL_KERNEL_SIZE] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    int Ky[SOBEL_KERNEL_SIZE][SOBEL_KERNEL_SIZE] = {{1, 2, 1}, {0, 0, 0}, {-1, -2, -1}};

    bmp_image **img = *physicalImage;
    bmp_image **result = (bmp_image **)calloc(height, sizeof(bmp_image *));
    bmp_image **theta_img = (bmp_image **)calloc(height, sizeof(bmp_image *));

    for (int i = 0; i < height; i++) {
        result[i] = (bmp_image *)calloc(width, sizeof(bmp_image));
        theta_img[i] = (bmp_image *)calloc(width, sizeof(bmp_image));
    }

    int k = SOBEL_KERNEL_SIZE / 2;

    for (int i = 2; i < height - 2; i++) {
        for(int j = 0; j < width; j++) {
            double new_fade_x = 0;
            double new_fade_y = 0;
            for (int m = -k; m < k + 1; m++) {
                for (int n = -k; n < k + 1; n++) {
                    if (j + n >= 0 && j + n < width) {
                        new_fade_x += (double)img[i + m][j + n].Blue * Kx[m + k][n + k];
                        new_fade_y += (double)img[i + m][j + n].Blue * Ky[m + k][n + k];
                    }
                }
            }

            double new_val = sqrt(new_fade_x * new_fade_x + new_fade_y * new_fade_y);
            result[i][j].Blue = new_val;
            result[i][j].Green = new_val;
            result[i][j].Red = new_val;

            theta_img[i][j].Blue = atan2(new_fade_x, new_fade_y);
            theta_img[i][j].Green = atan2(new_fade_x, new_fade_y);
            theta_img[i][j].Red = atan2(new_fade_x, new_fade_y);
        }
    }

    bmp_image max_pixel;
    max_pixel.Red = MIN_PIXEL_VALUE;
    max_pixel.Green = MIN_PIXEL_VALUE;
    max_pixel.Blue = MIN_PIXEL_VALUE;

    for (int i = 2; i < height - 2; i++) {
        for (int j = 0; j < width; j++) {
            if (result[i][j].Red > max_pixel.Red && result[i][j].Green > max_pixel.Green 
                && result[i][j].Blue > max_pixel.Blue) {
                max_pixel.Red = result[i][j].Red;
                max_pixel.Blue = result[i][j].Blue;
                max_pixel.Green = result[i][j].Green;
            }
        }
    }

    MPI_Send(&max_pixel, 1, MPI_PIXEL_TYPE, ROOT, 0, MPI_COMM_WORLD);
    MPI_Recv(&max_pixel, 1, MPI_PIXEL_TYPE, ROOT, 0, MPI_COMM_WORLD, NULL);

    for (int i = 2; i < height-2; i++) {
        for(int j = 0; j < width; j++) {
            result[i][j].Blue = (double)result[i][j].Blue / max_pixel.Blue * 255;
            result[i][j].Green = (double)result[i][j].Green / max_pixel.Green * 255;
            result[i][j].Red = (double)result[i][j].Red / max_pixel.Red * 255;
        }
    }

    free_image(physicalImage, height);
    *physicalImage = result;

    return theta_img;
}

// Function to perform non-maximum suppression to thin out the edges.
void non_max_suppression(bmp_image ***physicalImage, bmp_image **theta_img,
    int height, int width) {
    bmp_image **img = *physicalImage;
    
    bmp_image **result = (bmp_image **)calloc(height, sizeof(bmp_image *));
    for (int i = 0; i < height; i++) {
        result[i] = (bmp_image *)calloc(width, sizeof(bmp_image));
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            theta_img[i][j].Blue = theta_img[i][j].Blue * HALF_CIRCLE_ANGLE / M_PI;
            if (theta_img[i][j].Blue < 0) {
                theta_img[i][j].Blue += HALF_CIRCLE_ANGLE;
            }
        }
    }

    for (int i = 2; i < height - 2; i++) {
        for (int j = 0; j < width; j++) {
            double q = 255;
            double r = 255;

            if (((0 <= theta_img[i][j].Blue < 22.5) || (157.5 <= theta_img[i][j].Blue <= 180))
                && (j + 1 < width) && (j-1>=0)) {
                q = img[i][j+1].Blue;
                r = img[i][j-1].Blue;
            } else if ((22.5 <= theta_img[i][j].Blue < 67.5) && (i+1<height) && (i-1>=0) && (j-1>=0) && (j+1<width)) {
                q = img[i+1][j-1].Blue;
                r = img[i-1][j+1].Blue;
            } else if (67.5 <= theta_img[i][j].Blue < 112.5&& (i+1<height) && (i-1>=0)) {
                q = img[i+1][j].Blue;
                r = img[i-1][j].Blue;
            } else if (112.5 <= theta_img[i][j].Blue < 157.5 && (i+1<height) && (i-1>=0) && (j-1>=0) && (j+1<width)){
                q = img[i-1][j-1].Blue;
                r = img[i+1][j+1].Blue;
            }

            if ((img[i][j].Blue >= q) && (img[i][j].Blue >= r)) {
                result[i][j].Blue = img[i][j].Blue;
                result[i][j].Green = img[i][j].Blue;
                result[i][j].Red = img[i][j].Blue;
            }
            else {
                result[i][j].Blue = 0;
                result[i][j].Green = 0;
                result[i][j].Red = 0;
            }
        }
    }

    free_image(physicalImage, height);
    *physicalImage = result;
}

// Function to identifying 3 kinds of pixels: strong, weak, and non-relevant.
void double_threshold(bmp_image ***physicalImage, int height, int width,
	double low_threshold_ratio, double high_threshold_ratio) {
	bmp_image **img = *physicalImage;

	bmp_image **result = (bmp_image **)calloc(height, sizeof(bmp_image *));
    for (int i = 0; i < height; i++) {
        result[i] = (bmp_image *)calloc(width, sizeof(bmp_image));
    }

	bmp_image max_pixel;
	max_pixel.Red = MIN_PIXEL_VALUE;
	max_pixel.Green = MIN_PIXEL_VALUE;
	max_pixel.Blue = MIN_PIXEL_VALUE;

	for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
        	if (img[i][j].Red > max_pixel.Red && img[i][j].Green > max_pixel.Green 
        		&& img[i][j].Blue > max_pixel.Blue) {
        		max_pixel.Red = img[i][j].Red;
        		max_pixel.Blue = img[i][j].Blue;
        		max_pixel.Green = img[i][j].Green;
        	}
        }
    }

    MPI_Send(&max_pixel, 1, MPI_PIXEL_TYPE, ROOT, 0, MPI_COMM_WORLD);
    MPI_Recv(&max_pixel, 1, MPI_PIXEL_TYPE, ROOT, 0, MPI_COMM_WORLD, NULL);

    bmp_image high_threshold;
    high_threshold.Red = max_pixel.Red * high_threshold_ratio;
    high_threshold.Blue = max_pixel.Blue * high_threshold_ratio;
    high_threshold.Green = max_pixel.Green * high_threshold_ratio;

    bmp_image low_threshold;
    low_threshold.Red = high_threshold.Red * low_threshold_ratio;
    low_threshold.Blue = high_threshold.Blue * low_threshold_ratio;
    low_threshold.Green = high_threshold.Green * low_threshold_ratio;

    for (int i = 2; i < height - 2; i++) { 
        for (int j = 0; j < width; j++) {
        	if (img[i][j].Red >= high_threshold.Red && img[i][j].Green >= high_threshold.Green 
        		&& img[i][j].Blue >= high_threshold.Blue) {
        		result[i][j].Red = MAX_PIXEL_VALUE;
        		result[i][j].Green = MAX_PIXEL_VALUE;
        		result[i][j].Blue = MAX_PIXEL_VALUE;
        	} else if (img[i][j].Red >= low_threshold.Red && img[i][j].Green >= low_threshold.Green 
        		&& img[i][j].Blue >= low_threshold.Blue) {
        		result[i][j].Red = WEAK_PIXEL_VALUE;
        		result[i][j].Green = WEAK_PIXEL_VALUE;
        		result[i][j].Blue = WEAK_PIXEL_VALUE;
        	} else {
        		result[i][j].Red = MIN_PIXEL_VALUE;
        		result[i][j].Green = MIN_PIXEL_VALUE;
        		result[i][j].Blue = MIN_PIXEL_VALUE;
        	}
        }
    }

    free_image(physicalImage, height);
    *physicalImage = result;
}

int is_strong_pixel(bmp_image pixel) {
	if (pixel.Red == MAX_PIXEL_VALUE && pixel.Green == MAX_PIXEL_VALUE 
		&& pixel.Blue == MAX_PIXEL_VALUE) {
		return 1;
	}

	return 0;
}

// Function to transform weak pixels into strong ones.
void hysteresis(bmp_image ***physicalImage, int height, int width) {
	bmp_image **img = *physicalImage;

	bmp_image **result = (bmp_image **)calloc(height, sizeof(bmp_image *));
    for (int i = 0; i < height; i++) {
        result[i] = (bmp_image *)calloc(width, sizeof(bmp_image));
    }

	for (int i = 2; i < height - 2; i++) {
        for (int j = 0; j < width ; j++) {
        	if (img[i][j].Red == WEAK_PIXEL_VALUE && img[i][j].Green == WEAK_PIXEL_VALUE 
        		&& img[i][j].Blue == WEAK_PIXEL_VALUE) {
        		if ((i - 1 >= 0 && j - 1 >= 0 && is_strong_pixel(img[i-1][j-1])) ||
        			(i - 1 >= 0 && is_strong_pixel(img[i-1][j])) ||
        			(j - 1 >= 0 && is_strong_pixel(img[i][j-1])) ||
        			(i + 1 < height && j + 1 < width && is_strong_pixel(img[i+1][j+1])) ||
        			(i + 1 < height && is_strong_pixel(img[i+1][j])) ||
        			(j + 1 < width && is_strong_pixel(img[i][j+1])) ||
        			(i + 1 < height && j - 1 >= 0 && is_strong_pixel(img[i+1][j-1])) ||
        			(i - 1 >= 0 && j + 1 < width && is_strong_pixel(img[i-1][j+1]))) {
        			result[i][j].Red = MAX_PIXEL_VALUE;
        			result[i][j].Blue = MAX_PIXEL_VALUE;
        			result[i][j].Green = MAX_PIXEL_VALUE;
        		} else {
        			result[i][j].Red = MIN_PIXEL_VALUE;
        			result[i][j].Blue = MIN_PIXEL_VALUE;
        			result[i][j].Green = MIN_PIXEL_VALUE;
        		}
        	} else {
        		result[i][j].Red = img[i][j].Red;
        		result[i][j].Blue = img[i][j].Blue;
        		result[i][j].Green = img[i][j].Green;
        	}
    	}
    }

    free_image(physicalImage, height);
    *physicalImage = result;	
}

// Function to update the neighbors of each chunk of pixels. 
void update_padding(bmp_image ***physicalImage, int height, int width, int rank, int num_proc) {
    bmp_image **image = *physicalImage;

    if (rank != num_proc - 1) {
        MPI_Send(image[height - 4], width, MPI_PIXEL_TYPE, rank + 1, 0, MPI_COMM_WORLD);
        MPI_Send(image[height - 3], width, MPI_PIXEL_TYPE, rank + 1, 0, MPI_COMM_WORLD);

        MPI_Recv(image[height - 2], width, MPI_PIXEL_TYPE, rank + 1, 0, MPI_COMM_WORLD, NULL);
        MPI_Recv(image[height - 1], width, MPI_PIXEL_TYPE, rank + 1, 0, MPI_COMM_WORLD, NULL);
    }
    
    if (rank != ROOT + 1) {
        MPI_Recv(image[0], width, MPI_PIXEL_TYPE, rank - 1, 0, MPI_COMM_WORLD, NULL);
        MPI_Recv(image[1], width, MPI_PIXEL_TYPE, rank - 1, 0, MPI_COMM_WORLD, NULL);
            
        MPI_Send(image[2], width, MPI_PIXEL_TYPE, rank - 1, 0, MPI_COMM_WORLD);
        MPI_Send(image[3], width, MPI_PIXEL_TYPE, rank - 1, 0, MPI_COMM_WORLD);
    }

    *physicalImage = image;
}

void send_gaussian_kernel(double** gauss_kernel, int gauss_kernel_size, int dest_proc) {
    gauss_kernel_size = GAUSS_KERNEL_SIZE;
    MPI_Send(&gauss_kernel_size, 1, MPI_INT, dest_proc, 0, MPI_COMM_WORLD);
    for (int j = 0; j < gauss_kernel_size; j++) {
        MPI_Send(gauss_kernel[j], gauss_kernel_size, MPI_DOUBLE, dest_proc, 0, MPI_COMM_WORLD);
    }
}

double** recv_gaussian_kernel(int *gauss_kernel_size) {
    MPI_Recv(gauss_kernel_size, 1, MPI_INT, ROOT, 0, MPI_COMM_WORLD, NULL);

    double** gauss_kernel = (double **) calloc(*gauss_kernel_size, sizeof(double *));
    if(gauss_kernel == NULL){
        exit(ERROR_CODE);
    }

    for (int i = 0; i < *gauss_kernel_size; i++) {
        gauss_kernel[i] = (double *) calloc(*gauss_kernel_size, sizeof(double));
        if(gauss_kernel[i] == NULL){
            exit(ERROR_CODE);
        }
            
        MPI_Recv(gauss_kernel[i], *gauss_kernel_size, MPI_DOUBLE, ROOT, 0, MPI_COMM_WORLD, NULL);
    }

    return gauss_kernel;
}

void compute_max_pixel(int num_proc) {
    bmp_image max_pixel;
    bmp_image partial_max_pixel;

    max_pixel.Red = MIN_PIXEL_VALUE;
    max_pixel.Green = MIN_PIXEL_VALUE;
    max_pixel.Blue = MIN_PIXEL_VALUE;

    // Gether the maximum pixels of each chunk.
    for(int i = 1; i < num_proc; i++) {
        MPI_Recv(&partial_max_pixel, 1, MPI_PIXEL_TYPE, i, 0, MPI_COMM_WORLD, NULL);
        if (partial_max_pixel.Red > max_pixel.Red && partial_max_pixel.Green > max_pixel.Green 
            && partial_max_pixel.Blue > max_pixel.Blue) {
            max_pixel.Red = partial_max_pixel.Red;
            max_pixel.Green = partial_max_pixel.Green;
            max_pixel.Blue = partial_max_pixel.Blue;
        }
    }
    
    // Send the miximum pixel.
    for(int i = 1; i < num_proc; i++) {
        MPI_Send(&max_pixel, 1, MPI_PIXEL_TYPE, i, 0, MPI_COMM_WORLD);
    }
}

int main(int argc, char *argv[]) {
    bmp_image **image;
    int num_proc, rank;
    int gauss_kernel_size;
    double** gauss_kernel;
    int num_lines, num_columns;

    offsets[0] = offsetof(bmp_image, Green);
    offsets[1] = offsetof(bmp_image, Red);
    offsets[2] = offsetof(bmp_image, Blue);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Type_create_struct(num_items, block_lengths, offsets, types, &MPI_PIXEL_TYPE);
    MPI_Type_commit(&MPI_PIXEL_TYPE);

    double start_time, end_time;
    if (rank == ROOT) {
        FILE *inputImage = fopen("wp.bmp", "rb");

        // Read the image file header and info header.
        bmp_fileheader *fileHeader = calloc(1,sizeof(bmp_fileheader));
        if(fileHeader == NULL){
            exit(ERROR_CODE);
        }
        bmp_infoheader *infoHeader = calloc(1,sizeof(bmp_infoheader));
        if(infoHeader == NULL){
            exit(ERROR_CODE);
        }
        read_data(inputImage, fileHeader, infoHeader);

        // Read the image data.
        image = calloc(infoHeader -> height, sizeof(bmp_image *));
        if(image == NULL){
            exit(ERROR_CODE);
        }
        for(int i = 0; i < infoHeader -> height; i++) {
            image[i] = calloc(infoHeader -> width, sizeof(bmp_image));
            if(image[i] == NULL) {
                exit(ERROR_CODE);
            }
        }
        read_physicalImage(image, inputImage, infoHeader, fileHeader);

        fclose(inputImage);

        start_time = MPI_Wtime();

        // Compute the gaussian kernel.
        gauss_kernel = compute_gaussian_kernel(GAUSS_KERNEL_SIZE, GAUSS_KERNEL_SIGMA);

        // Devide the image into chuncks.
        int start[num_proc], end[num_proc];
        for (int i = 1; i < num_proc; i++) {
            start[i] = (i - 1) * (double)infoHeader->height / (num_proc - 1);
            end[i] = i * (double)infoHeader->height / (num_proc - 1) < infoHeader->height
                ? i * (double)infoHeader->height / (num_proc - 1)
                : infoHeader->height;

            num_lines = end[i] - start[i];
            num_columns = infoHeader->width;

            // Send the pixel chuncks to the worker processes.
            MPI_Send(&num_lines, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&num_columns, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            for (int j = 0; j < num_lines; j++) {
                MPI_Send(image[j + start[i]], num_columns, MPI_PIXEL_TYPE, i, j, MPI_COMM_WORLD);
            }

            // Send the gaussian kernel to the worker processes.
            send_gaussian_kernel(gauss_kernel, gauss_kernel_size, i);
        }
        num_lines = infoHeader->height;

        // Compute the maximum pixel for sobel_filters.
        compute_max_pixel(num_proc);

        // Compute the maximum pixel for double_threshold.
        compute_max_pixel(num_proc);

        // Reconstruct the image.
        for (int i = 1; i < num_proc; i++) {
            num_lines = end[i] - start[i];
            for (int j = 0; j < num_lines; j++) {
                MPI_Recv(image[j + start[i]], num_columns, MPI_PIXEL_TYPE, i, j, MPI_COMM_WORLD, NULL);
            }
        }

        end_time = MPI_Wtime();

        printf("Total time: %f\n", (double)(end_time - start_time));

        // Write the image.
        FILE *out = fopen("image_mpi.bmp", "wb");
        print(out, fileHeader, infoHeader);
        print_physicalImage(image, out, fileHeader, infoHeader);
        fclose(out);

        free(fileHeader);
        free(infoHeader);
    } else {
        // Receive the image size from proc 0.
        MPI_Recv(&num_lines, 1, MPI_INT, ROOT, 0, MPI_COMM_WORLD, NULL);
        MPI_Recv(&num_columns, 1, MPI_INT, ROOT, 0, MPI_COMM_WORLD, NULL);

        // Receive the image chunck from proc 0 (add 2 lines with zeros above and 
        // below as padding).
        num_lines += 2 * PADDING_SIZE;
        image = calloc(num_lines, sizeof(bmp_image *));
        if(image == NULL){
            exit(ERROR_CODE);
        }
        for(int i = 0; i < num_lines; i++) {
            image[i] = calloc(num_columns, sizeof(bmp_image));
            if(image[i] == NULL) {
                exit(ERROR_CODE);
            }
        }
        for(int i = PADDING_SIZE; i < num_lines - PADDING_SIZE; i++) {
            MPI_Recv(image[i], num_columns, MPI_PIXEL_TYPE, ROOT, i - PADDING_SIZE, MPI_COMM_WORLD, NULL);
        }

        // Receive the gaussian kernel from proc 0.
        gauss_kernel = recv_gaussian_kernel(&gauss_kernel_size);

        // Turn image into black and white.
        rgb2gray(&image, num_lines, num_columns);

        // Blur image.
        update_padding(&image, num_lines, num_columns, rank,num_proc);
        image_noise_reduction(&image, num_lines, num_columns, gauss_kernel, gauss_kernel_size / 2);

        // Detect the edge intensity and direction.
        update_padding(&image, num_lines, num_columns, rank,num_proc);
        bmp_image **theta_image = sobel_filters(&image, num_lines, num_columns);

        // Thin out the edges.
        update_padding(&image, num_lines, num_columns, rank,num_proc);
        non_max_suppression(&image, theta_image, num_lines, num_columns);

        // Identifying strong, weak, and non-relevant pixels.
        update_padding(&image, num_lines, num_columns, rank,num_proc);
        double_threshold(&image, num_lines, num_columns, LOW_THRESHOLD_RATIO, HIGH_THRESHOLD_RATIO);

        // Edge tracking by hysteresis.
        update_padding(&image, num_lines, num_columns, rank,num_proc);
        hysteresis(&image, num_lines, num_columns);

        // Transmit to proc 0 the image chuncks (ignore padding).
        for (int i = PADDING_SIZE; i < num_lines - PADDING_SIZE; i++) {
            MPI_Send(image[i], num_columns, MPI_PIXEL_TYPE, ROOT, i - PADDING_SIZE, MPI_COMM_WORLD);
        }
    }

    free_image(&image, num_lines);
    free_gaussian_kernel(&gauss_kernel, gauss_kernel_size);

    MPI_Finalize();

    return 0;
}

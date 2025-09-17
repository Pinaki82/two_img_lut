// Last Change: 2025-09-17  Wednesday: 03:21:50 PM
#ifndef LUTGEN_H
#define LUTGEN_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

/** A 3-component float vector (RGB) **/
typedef struct {
  double r;
  double g;
  double b;
} Vec3;

/** A cell in the LUT: accumulated sum + weight **/
typedef struct {
  double r_sum;
  double g_sum;
  double b_sum;
  double weight;
} LUTCell;

/**
   Generate a 3D LUT mapping from source and target images.

   Parameters:
     src_data, tgt_data: pointers to source & target image pixel data,
       assumed to be 8-bit RGB, stride = width * 3, no alpha.
     width, height: dimensions of images (must match).
     lut_size: the size N of the LUT (N×N×N).
     out_lut: preallocated array of LUTCell of size N*N*N, to be filled.

   Return: 0 on success, non-zero on error.
*/
int generate_lut(
        const unsigned char *src_data,
        const unsigned char *tgt_data,
        int width,
        int height,
        int lut_size,
        LUTCell *out_lut
);

/**
   Normalise a LUT: for each cell, divide accumulated sums by weight.
   Cells with zero weight are set to some marker (e.g. negative).
*/
void normalize_lut(
        LUTCell *lut,
        int lut_size
);

/**
   Fill missing cells in LUT (cells with no data) via nearest neighbour or other method.
*/
void fill_lut_missing_cells(
        LUTCell *lut,
        int lut_size
);

/**
   Write LUT to .cube file.

   Parameters:
     lut: array of LUTCell (normalised)
     lut_size: size N
     filename: output path

   Return: 0 on success, non-zero on error.
*/
int write_cube_lut(
        const LUTCell *lut,
        int lut_size,
        const char *filename
);

#ifdef __cplusplus
}
#endif

#endif // LUTGEN_H



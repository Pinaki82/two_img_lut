// lutgen.c
// gcc -g -std=c11 -Wall -Wextra -pedantic two_img_lut.c -o two_img_lut -lm
// ./two_img_lut source.png target.png lut_size output.cube
// Example: ./two_img_lut mpv-shot0001.png mpv-shot0001-corrected.png 33 00010.cube
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Replace with your image loading library, e.g. stb_image
#define STB_IMAGE_IMPLEMENTATION
#include "stb-master/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb-master/stb_image_write.h"
#include "lutgen.h"

typedef struct {
    double r, g, b;
    double w;
} Cell;

static int clamp_int(int x, int lo, int hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s source.png target.png lut_size output.cube\n", argv[0]);
        return 1;
    }
    const char *pathA = argv[1];
    const char *pathB = argv[2];
    int N = atoi(argv[3]);
    const char *outcube = argv[4];
    if (N < 2) {
        fprintf(stderr, "lut_size must be >= 2\n");
        return 1;
    }

    int wA, hA, chA;
    unsigned char *dataA = stbi_load(pathA, &wA, &hA, &chA, 3);
    if (!dataA) {
        fprintf(stderr, "Failed to load %s\n", pathA);
        return 1;
    }
    int wB, hB, chB;
    unsigned char *dataB = stbi_load(pathB, &wB, &hB, &chB, 3);
    if (!dataB) {
        fprintf(stderr, "Failed to load %s\n", pathB);
        stbi_image_free(dataA);
        return 1;
    }
    if (wA != wB || hA != hB) {
        fprintf(stderr, "Source and target image sizes do not match\n");
        stbi_image_free(dataA);
        stbi_image_free(dataB);
        return 1;
    }

    // Allocate LUT cells
    size_t total = (size_t)N * N * N;
    Cell *lut = (Cell*) malloc(total * sizeof(Cell));
    if (!lut) {
        fprintf(stderr, "Out of memory allocating LUT\n");
        stbi_image_free(dataA);
        stbi_image_free(dataB);
        return 1;
    }
    for (size_t i = 0; i < total; i++) {
        lut[i].r = lut[i].g = lut[i].b = 0.0;
        lut[i].w = 0.0;
    }

    // Accumulate over all pixels
    for (int y = 0; y < hA; y++) {
        for (int x = 0; x < wA; x++) {
            int idx = (y * wA + x) * 3;
            double rA = dataA[idx] / 255.0;
            double gA = dataA[idx+1] / 255.0;
            double bA = dataA[idx+2] / 255.0;
            double rB = dataB[idx] / 255.0;
            double gB = dataB[idx+1] / 255.0;
            double bB = dataB[idx+2] / 255.0;

            // map A colour to LUT space
            double fi = rA * (N - 1);
            double fj = gA * (N - 1);
            double fk = bA * (N - 1);
            int i0 = (int)floor(fi);
            int j0 = (int)floor(fj);
            int k0 = (int)floor(fk);
            double di = fi - i0;
            double dj = fj - j0;
            double dk = fk - k0;

            for (int dii = 0; dii <= 1; dii++) {
                int ii = i0 + dii;
                if (ii < 0 || ii >= N) continue;
                double wi = (dii == 0) ? (1.0 - di) : di;
                for (int djj = 0; djj <= 1; djj++) {
                    int jj = j0 + djj;
                    if (jj < 0 || jj >= N) continue;
                    double wj = (djj == 0) ? (1.0 - dj) : dj;
                    for (int dkk = 0; dkk <= 1; dkk++) {
                        int kk = k0 + dkk;
                        if (kk < 0 || kk >= N) continue;
                        double wk = (dkk == 0) ? (1.0 - dk) : dk;
                        double w = wi * wj * wk;
                        size_t cell_index = ( (ii * N + jj) * N + kk );
                        lut[cell_index].r += rB * w;
                        lut[cell_index].g += gB * w;
                        lut[cell_index].b += bB * w;
                        lut[cell_index].w += w;
                    }
                }
            }
        }
    }

    // Normalize (average) cells
    for (size_t i = 0; i < total; i++) {
        if (lut[i].w > 0.0) {
            lut[i].r /= lut[i].w;
            lut[i].g /= lut[i].w;
            lut[i].b /= lut[i].w;
        } else {
            // mark missing by setting to negative
            lut[i].r = lut[i].g = lut[i].b = -1.0;
        }
    }

    // Fill missing cells by nearest neighbor
    // For each missing cell, search outward in 3D for a filled one
    // (Simple brute force search)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                size_t ci = ( (i * N + j) * N + k );
                if (lut[ci].r >= 0.0) continue;  // already filled

                // search radius
                int maxr = N;  // up to whole cube
                int best_i=-1, best_j=-1, best_k=-1;
                double best_dist = 1e30;
                for (int dii = -maxr; dii <= maxr; dii++) {
                    int ii = i + dii;
                    if (ii < 0 || ii >= N) continue;
                    for (int djj = -maxr; djj <= maxr; djj++) {
                        int jj = j + djj;
                        if (jj < 0 || jj >= N) continue;
                        for (int dkk = -maxr; dkk <= maxr; dkk++) {
                            int kk = k + dkk;
                            if (kk < 0 || kk >= N) continue;
                            size_t ci2 = ( (ii * N + jj) * N + kk );
                            if (lut[ci2].r < 0.0) continue;  // still missing
                            // distance in cell space
                            double dist = dii*dii + djj*djj + dkk*dkk;
                            if (dist < best_dist) {
                                best_dist = dist;
                                best_i = ii; best_j = jj; best_k = kk;
                            }
                        }
                    }
                }
                if (best_i >= 0) {
                    size_t ci2 = ( (best_i * N + best_j) * N + best_k );
                    lut[ci].r = lut[ci2].r;
                    lut[ci].g = lut[ci2].g;
                    lut[ci].b = lut[ci2].b;
                } else {
                    // fallback: set to identity
                    lut[ci].r = i / (double)(N - 1);
                    lut[ci].g = j / (double)(N - 1);
                    lut[ci].b = k / (double)(N - 1);
                }
            }
        }
    }

    // Write .cube file
    FILE *f = fopen(outcube, "w");
    if (!f) {
        fprintf(stderr, "Failed to open %s for writing\n", outcube);
        free(lut);
        stbi_image_free(dataA);
        stbi_image_free(dataB);
        return 1;
    }

    fprintf(f, "# Generated by lutgen\n");
    fprintf(f, "TITLE \"Generated LUT from %s -> %s\"\n", pathA, pathB);
    fprintf(f, "LUT_3D_SIZE %d\n", N);
    fprintf(f, "DOMAIN_MIN 0.0 0.0 0.0\n");
    fprintf(f, "DOMAIN_MAX 1.0 1.0 1.0\n");

    // Order: r fastest, then g, then b (common ordering)
    for (int b = 0; b < N; b++) {
        for (int g = 0; g < N; g++) {
            for (int r = 0; r < N; r++) {
                size_t ci = ( (r * N + g) * N + b );
                double orr = lut[ci].r;
                double ogg = lut[ci].g;
                double obb = lut[ci].b;
                // clamp
                if (orr < 0.0) orr = r / (double)(N - 1);
                if (ogg < 0.0) ogg = g / (double)(N - 1);
                if (obb < 0.0) obb = b / (double)(N - 1);
                fprintf(f, "%.6f %.6f %.6f\n", orr, ogg, obb);
            }
        }
    }

    fclose(f);

    // Cleanup
    free(lut);
    stbi_image_free(dataA);
    stbi_image_free(dataB);

    printf("Wrote LUT to %s\n", outcube);
    return 0;
}


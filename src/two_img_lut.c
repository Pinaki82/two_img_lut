// Last Change: 2025-09-17  Wednesday: 07:46:59 PM
// two_img_lut.c

/*
   two_img_lut — generate 3D LUT from two images
   Author: ChatGPT
   Date: 9 September 2025

   License: MIT-0 License

   Copyright (c) 2025 Pinaki Gupta

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including, without limitation, the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/


// gcc -g -std=c11 -Wall -Wextra -pedantic two_img_lut.c -o two_img_lut -lm
// ./two_img_lut source.png target.png lut_size[optional] output.cube

/*

  Example:
  ./two_img_lut mpv-shot0001.png mpv-shot0001-corrected.png 00010.cube
  Or,
  ./two_img_lut mpv-shot0001.png mpv-shot0001-corrected.png 33 00010.cube

*/

#define _POSIX_C_SOURCE 199309L

#ifndef __linux__
  #error "This program currently supports only Linux"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

// stb_image - single-translation-unit usage
#define STB_IMAGE_IMPLEMENTATION
#include "external/stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "external/stb/stb_image_write.h"

#define DEFAULT_LUT_SIZE 129
#define MAX_FILL_RADIUS 3 /* kept for config, not used in BFS */

// Cell stored as floats for memory efficiency and speed
typedef struct {
  float r, g, b; /* accumulated (or averaged) color */
  double w;      /* accumulated weight */
} Cell;

/* coordinate triple */
typedef struct {
  int i, j, k;
} Coord;

/* circular queue of Coord */
typedef struct {
  Coord *data;
  size_t head, tail, capacity;
} Queue;

static Queue *queue_create(size_t capacity) {
  Queue *q = malloc(sizeof(Queue));

  if(!q) {
    return NULL;
  }

  q->data = malloc(sizeof(Coord) * capacity);

  if(!q->data) {
    free(q);
    return NULL;
  }

  q->head = q->tail = 0;
  q->capacity = capacity;
  return q;
}

static void queue_free(Queue *q) {
  if(!q) {
    return;
  }

  free(q->data);
  free(q);
}

static int queue_is_empty(const Queue *q) {
  return q->head == q->tail;
}

/* returns 1 on success, 0 if full */
static int queue_push(Queue *q, Coord c) {
  size_t next = (q->tail + 1) % q->capacity;

  if(next == q->head) {
    return 0;  /* full */
  }

  q->data[q->tail] = c;
  q->tail = next;
  return 1;
}

/* returns 1 on success (pop), 0 if empty */
static int queue_pop(Queue *q, Coord *out) {
  if(queue_is_empty(q)) {
    return 0;
  }

  *out = q->data[q->head];
  q->head = (q->head + 1) % q->capacity;
  return 1;
}

/* multi-source BFS fill that averages already-filled neighbours */
static void bfs_fill(Cell *lut, int N) {
  size_t total = (size_t)(N * N * N);
  unsigned char *visited = calloc(total, 1);

  if(!visited) {
    return;
  }

  Queue *q = queue_create(total + 16);

  if(!q) {
    free(visited);
    return;
  }

  /* enqueue every already-filled cell (marker: r >= 0) */
  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) {
      for(int k = 0; k < N; ++k) {
        size_t idx = (size_t)((size_t)i * (size_t)N + (size_t)j) * (size_t)N + (size_t)k;

        if(lut[idx].r >= 0.0f) {
          visited[idx] = 1;
          queue_push(q, (Coord) {
            i, j, k
          });
        }
      }
    }
  }

  const int di[6] = {+1, -1, 0, 0, 0, 0};
  const int dj[6] = { 0, 0, +1, -1, 0, 0};
  const int dk[6] = { 0, 0, 0, 0, +1, -1};
  Coord cur;

  while(queue_pop(q, &cur)) {
    size_t cur_idx = ((size_t)cur.i * (size_t)N + (size_t)cur.j) * (size_t)N + (size_t)cur.k;

    for(int n = 0; n < 6; ++n) {
      int ni = cur.i + di[n];
      int nj = cur.j + dj[n];
      int nk = cur.k + dk[n];

      if(ni < 0 || ni >= N || nj < 0 || nj >= N || nk < 0 || nk >= N) {
        continue;
      }

      size_t nidx = ((size_t)ni * (size_t)N + (size_t)nj) * (size_t)N + (size_t)nk;

      if(visited[nidx]) {
        continue;
      }

      /* compute average of already-filled neighbours of (ni,nj,nk) */
      double sr = 0.0, sg = 0.0, sb = 0.0;
      int count = 0;

      for(int m = 0; m < 6; ++m) {
        int mi = ni + di[m], mj = nj + dj[m], mk = nk + dk[m];

        if(mi < 0 || mi >= N || mj < 0 || mj >= N || mk < 0 || mk >= N) {
          continue;
        }

        size_t midx = ((size_t)mi * (size_t)N + (size_t)mj) * (size_t)N + (size_t)mk;

        if(lut[midx].r >= 0.0f) {
          sr += lut[midx].r;
          sg += lut[midx].g;
          sb += lut[midx].b;
          ++count;
        }
      }

      if(count > 0) {
        lut[nidx].r = (float)(sr / count);
        lut[nidx].g = (float)(sg / count);
        lut[nidx].b = (float)(sb / count);
      }

      else {
        /* fallback: copy from current (rare) */
        lut[nidx].r = lut[cur_idx].r;
        lut[nidx].g = lut[cur_idx].g;
        lut[nidx].b = lut[cur_idx].b;
      }

      visited[nidx] = 1;
      queue_push(q, (Coord) {
        ni, nj, nk
      });
    }
  }

  queue_free(q);
  free(visited);
}

int main(int argc, char *argv[]) {
  const char *pathA = NULL, *pathB = NULL, *outcube = NULL;
  int lut_size = DEFAULT_LUT_SIZE;

  if(argc == 4) {
    pathA = argv[1];
    pathB = argv[2];
    outcube = argv[3];
  }

  else if(argc == 5) {
    pathA = argv[1];
    pathB = argv[2];
    lut_size = atoi(argv[3]);
    outcube = argv[4];
  }

  else {
    fprintf(stderr, "Usage: %s source.png target.png [lut_size] output.cube\n", argv[0]);
    return 1;
  }

  if(lut_size < 2 || lut_size > 2000) {
    fprintf(stderr, "Error: lut_size must be >=2 and reasonably small (you passed %d)\n", lut_size);
    return 1;
  }

  struct timespec t_start, t_end;

  if(clock_gettime(CLOCK_MONOTONIC, &t_start) != 0) {
    perror("clock_gettime (start)");
  }

  int wA, hA, chA;
  unsigned char *dataA = stbi_load(pathA, &wA, &hA, &chA, 3);

  if(!dataA) {
    fprintf(stderr, "Failed to load %s\n", pathA);
    return 1;
  }

  int wB, hB, chB;
  unsigned char *dataB = stbi_load(pathB, &wB, &hB, &chB, 3);

  if(!dataB) {
    fprintf(stderr, "Failed to load %s\n", pathB);
    stbi_image_free(dataA);
    return 1;
  }

  if(wA != wB || hA != hB) {
    fprintf(stderr, "Source and target image sizes do not match\n");
    stbi_image_free(dataA); stbi_image_free(dataB);
    return 1;
  }

  size_t total = (size_t)lut_size * lut_size * lut_size;
  Cell *lut = malloc(total * sizeof(Cell));

  if(!lut) {
    fprintf(stderr, "Out of memory allocating LUT (%zu cells)\n", total);
    stbi_image_free(dataA);
    stbi_image_free(dataB);
    return 1;
  }

  /* INITIALIZE accumulators TO ZERO (important!) */
  for(size_t i = 0; i < total; ++i) {
    lut[i].r = 0.0f;
    lut[i].g = 0.0f;
    lut[i].b = 0.0f;
    lut[i].w = 0.0;
  }

  /* Accumulate: distribute each pixel's target colour into up to 8 surrounding LUT cells.
     Parallelized with OpenMP if available; uses atomic updates on cell accumulators. */
#if defined(_OPENMP)
  #pragma omp parallel for schedule(dynamic)
#endif

  for(int y = 0; y < hA; ++y) {
    for(int x = 0; x < wA; ++x) {
      int idx_px = (y * wA + x) * 3;
      double rA = dataA[idx_px + 0] / 255.0;
      double gA = dataA[idx_px + 1] / 255.0;
      double bA = dataA[idx_px + 2] / 255.0;
      double rB = dataB[idx_px + 0] / 255.0;
      double gB = dataB[idx_px + 1] / 255.0;
      double bB = dataB[idx_px + 2] / 255.0;
      double fi = rA * (lut_size - 1);
      double fj = gA * (lut_size - 1);
      double fk = bA * (lut_size - 1);
      int i0 = (int)floor(fi);
      int j0 = (int)floor(fj);
      int k0 = (int)floor(fk);
      double di = fi - i0;
      double dj = fj - j0;
      double dk = fk - k0;

      for(int dii = 0; dii <= 1; ++dii) {
        int ii = i0 + dii;

        if(ii < 0 || ii >= lut_size) {
          continue;
        }

        double wi = (dii == 0) ? (1.0 - di) : di;

        for(int djj = 0; djj <= 1; ++djj) {
          int jj = j0 + djj;

          if(jj < 0 || jj >= lut_size) {
            continue;
          }

          double wj = (djj == 0) ? (1.0 - dj) : dj;

          for(int dkk = 0; dkk <= 1; ++dkk) {
            int kk = k0 + dkk;

            if(kk < 0 || kk >= lut_size) {
              continue;
            }

            double wk = (dkk == 0) ? (1.0 - dk) : dk;
            double w = wi * wj * wk;
            size_t cell_index = ((size_t)ii * lut_size + (size_t)jj) * lut_size + (size_t)kk;
            /* atomic adds to preserve correctness when parallelized */
#if defined(_OPENMP)
            #pragma omp atomic
#endif
            lut[cell_index].w += w;
#if defined(_OPENMP)
            #pragma omp atomic
#endif
            lut[cell_index].r += (float)(rB * w);
#if defined(_OPENMP)
            #pragma omp atomic
#endif
            lut[cell_index].g += (float)(gB * w);
#if defined(_OPENMP)
            #pragma omp atomic
#endif
            lut[cell_index].b += (float)(bB * w);
          }
        }
      }
    }
  } /* end accumulation */

  /* Normalize: convert accumulated sums into average RGB and mark missing cells with -1 */
  for(size_t i = 0; i < total; ++i) {
    if(lut[i].w > 0.0) {
      lut[i].r = (float)(lut[i].r / lut[i].w);
      lut[i].g = (float)(lut[i].g / lut[i].w);
      lut[i].b = (float)(lut[i].b / lut[i].w);
    }

    else {
      /* missing marker */
      lut[i].r = lut[i].g = lut[i].b = -1.0f;
    }
  }

  /* Fill holes using multi-source BFS averaging */
  bfs_fill(lut, lut_size);
  /* Write .cube file (ASCII) */
  FILE *f = fopen(outcube, "w");

  if(!f) {
    fprintf(stderr, "Failed to open %s for writing\n", outcube);
    free(lut);
    stbi_image_free(dataA);
    stbi_image_free(dataB);
    return 1;
  }

  fprintf(f, "# Generated by two_img_lut\n");
  fprintf(f, "TITLE \"Generated LUT from %s -> %s\"\n", pathA, pathB);
  fprintf(f, "LUT_3D_SIZE %d\n", lut_size);
  fprintf(f, "DOMAIN_MIN 0.0 0.0 0.0\n");
  fprintf(f, "DOMAIN_MAX 1.0 1.0 1.0\n");

  /* cube ordering: r fastest, then g, then b (common) */
  for(int b = 0; b < lut_size; ++b) {
    for(int g = 0; g < lut_size; ++g) {
      for(int r = 0; r < lut_size; ++r) {
        size_t ci = ((size_t)r * lut_size + (size_t)g) * lut_size + (size_t)b;
        double orr = lut[ci].r;
        double ogg = lut[ci].g;
        double obb = lut[ci].b;

        if(orr < 0.0) {
          orr = r / (double)(lut_size - 1);
        }

        if(ogg < 0.0) {
          ogg = g / (double)(lut_size - 1);
        }

        if(obb < 0.0) {
          obb = b / (double)(lut_size - 1);
        }

        /* clamp (defensive) */
        if(orr < 0.0) {
          orr = 0.0;
        }

        if(ogg < 0.0) {
          ogg = 0.0;
        }

        if(obb < 0.0) {
          obb = 0.0;
        }

        if(orr > 1.0) {
          orr = 1.0;
        }

        if(ogg > 1.0) {
          ogg = 1.0;
        }

        if(obb > 1.0) {
          obb = 1.0;
        }

        fprintf(f, "%.6f %.6f %.6f\n", orr, ogg, obb);
      }
    }
  }

  fclose(f);
  free(lut);
  stbi_image_free(dataA);
  stbi_image_free(dataB);

  if(clock_gettime(CLOCK_MONOTONIC, &t_end) != 0) {
    perror("clock_gettime (end)");
  }

  else {
    long sec = t_end.tv_sec - t_start.tv_sec;
    long nsec = t_end.tv_nsec - t_start.tv_nsec;

    if(nsec < 0) {
      sec -= 1;
      nsec += 1000000000L;
    }

    double elapsed = (double)sec + (double)nsec / 1e9;
    printf("Wrote %s (LUT size %d). Elapsed time: %.3f s\n", outcube, lut_size, elapsed);
  }

  return 0;
}



/*
	Copyright (c) 2023 ByteBit/xtreme8000

	This file is part of DXTgopt.

	DXTgopt is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	DXTgopt is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with DXTgopt.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "lodepng.h"

#define blue(c) ((c)&31)
#define green(c) (((c) >> 5) & 63)
#define red(c) (((c) >> 11) & 31)
#define rgb32(r, g, b) (((int)(b) << 16) | ((int)(g) << 8) | (int)(r))

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

struct color_result {
	uint32_t error;
	uint16_t col1, col2;
};

struct color_result color_res_min(struct color_result a,
								  struct color_result b) {
	return a.error < b.error ? a : b;
}

struct color_result color_res_create() {
	return (struct color_result) {.error = UINT32_MAX};
}

#pragma omp declare reduction(min : struct color_result : omp_out              \
								  = color_res_min(omp_out, omp_in))            \
	initializer(omp_priv = color_res_create())

uint64_t time_get() {
	struct timespec now;
	timespec_get(&now, TIME_UTC);
	return ((uint64_t)now.tv_sec) * 1000 + ((uint64_t)now.tv_nsec) / 1000000;
}

void complete_palette(uint32_t col1, uint32_t col2, int* a, int* b, int* c,
					  int* d, int* e) {
	// https://github.com/svn2github/libsquish/blob/master/colourblock.cpp#L140
	// https://github.com/GNOME/gimp/blob/master/plug-ins/file-dds/dxt.c#L119
	// https://github.com/dolphin-emu/dolphin/blob/master/Source/Core/VideoCommon/TextureDecoder_Generic.cpp#L144
	a[0] = (red(col1) << 3) | (red(col1) >> 2);
	a[1] = (green(col1) << 2) | (green(col1) >> 4);
	a[2] = (blue(col1) << 3) | (blue(col1) >> 2);

	b[0] = (red(col2) << 3) | (red(col2) >> 2);
	b[1] = (green(col2) << 2) | (green(col2) >> 4);
	b[2] = (blue(col2) << 3) | (blue(col2) >> 2);

	c[0] = (a[0] * 2 + b[0]) / 3;
	c[1] = (a[1] * 2 + b[1]) / 3;
	c[2] = (a[2] * 2 + b[2]) / 3;

	d[0] = (a[0] + b[0] * 2) / 3;
	d[1] = (a[1] + b[1] * 2) / 3;
	d[2] = (a[2] + b[2] * 2) / 3;

	e[0] = (a[0] + b[0] + 1) / 2;
	e[1] = (a[1] + b[1] + 1) / 2;
	e[2] = (a[2] + b[2] + 1) / 2;
}

uint32_t col_pal_diff(uint32_t imgcol, int* pal) {
	int r = imgcol & 0xFF;
	int g = (imgcol >> 8) & 0xFF;
	int b = imgcol >> 16;

	int x = (r - pal[0]) * (r - pal[0]);
	int y = (g - pal[1]) * (g - pal[1]);
	int z = (b - pal[2]) * (b - pal[2]);
	return x + y + z;
}

void palette_error(uint32_t* image, int* a, int* b, int* c, int* d, int* e,
				   uint32_t* total_error) {
	total_error[0] = 0;
	total_error[1] = 0;

	for(size_t k = 0; k < 16; k++) {
		uint32_t err1 = col_pal_diff(*image, a);
		uint32_t err2 = col_pal_diff(*image, b);
		uint32_t err3 = col_pal_diff(*image, c);
		uint32_t err4 = col_pal_diff(*image, d);
		uint32_t err5 = col_pal_diff(*image, e);
		image++;

		total_error[0] += min(min(min(err1, err2), err3), err4);
		total_error[1] += min(min(err1, err2), err5);
	}
}

int best_fit(uint32_t image, bool mode4, int* a, int* b, int* c, int* d,
			 int* e) {
	int best = 0;
	uint32_t err = col_pal_diff(image, a);
	uint32_t err2 = col_pal_diff(image, b);
	uint32_t err3 = col_pal_diff(image, mode4 ? c : e);

	if(err2 < err) {
		err = err2;
		best = 1;
	}

	if(err3 < err) {
		err = err3;
		best = 2;
	}

	if(mode4) {
		uint32_t err4 = col_pal_diff(image, d);

		if(err4 < err) {
			err = err4;
			best = 3;
		}
	}

	return best;
}

void print_block_info(uint32_t* image, uint16_t* colors) {
	assert(image && colors);

	int a[3], b[3], c[3], d[3], e[3];
	complete_palette(colors[0], colors[1], a, b, c, d, e);

	uint32_t total_error[2];
	palette_error(image, a, b, c, d, e, total_error);

	if(colors[1] < colors[0]) {
		printf("4 color mode\n");
		printf("block error: %u\n", total_error[0]);
		printf("color0: (%u, %u, %u)\n", a[0], a[1], a[2]);
		printf("color1: (%u, %u, %u)\n", b[0], b[1], b[2]);
		printf("color2: (%u, %u, %u)\n", c[0], c[1], c[2]);
		printf("color3: (%u, %u, %u)\n", d[0], d[1], d[2]);
	} else {
		printf("3 color mode\n");
		printf("block error: %u\n", total_error[1]);
		printf("color0: (%u, %u, %u)\n", a[0], a[1], a[2]);
		printf("color1: (%u, %u, %u)\n", b[0], b[1], b[2]);
		printf("color2: (%u, %u, %u)\n", e[0], e[1], e[2]);
	}
}

void encode_16block(uint32_t* image, uint16_t* colors, uint8_t* output) {
	assert(image && colors && output);

	int a[3], b[3], c[3], d[3], e[3];
	complete_palette(colors[0], colors[1], a, b, c, d, e);

	output[0] = colors[0] & 0xFF;
	output[1] = colors[0] >> 8;

	output[2] = colors[1] & 0xFF;
	output[3] = colors[1] >> 8;

	bool mode4 = colors[1] < colors[0];

	for(size_t k = 0; k < 16; k += 4) {
		uint8_t res = best_fit(image[k + 0], mode4, a, b, c, d, e);
		res |= (best_fit(image[k + 1], mode4, a, b, c, d, e) << 2);
		res |= (best_fit(image[k + 2], mode4, a, b, c, d, e) << 4);
		res |= (best_fit(image[k + 3], mode4, a, b, c, d, e) << 6);
		output[4 + k / 4] = res;
	}
}

void compute_16block_fast(uint32_t* image, uint16_t* colors) {
	assert(image && colors);

	int col_min[3] = {255, 255, 255};
	int col_max[3] = {0, 0, 0};

	for(size_t k = 0; k < 16; k++) {
		int col[3] = {image[k] & 0xFF, (image[k] >> 8) & 0xFF, image[k] >> 16};

		for(size_t i = 0; i < 3; i++) {
			col_min[i] = min(col_min[i], col[i]);
			col_max[i] = max(col_max[i], col[i]);
		}
	}

	colors[0] = ((col_min[0] >> 3) << 11) | ((col_min[1] >> 2) << 5)
		| (col_min[2] >> 3);
	colors[1] = ((col_max[0] >> 3) << 11) | ((col_max[1] >> 2) << 5)
		| (col_max[2] >> 3);

	if(colors[0] < colors[1]) {
		uint16_t tmp = colors[0];
		colors[0] = colors[1];
		colors[1] = tmp;
	}
}

void compute_16block(uint32_t* image, uint16_t* colors) {
	assert(image && colors);

	struct color_result result = color_res_create();

#pragma omp parallel for schedule(dynamic) reduction(min : result)
	for(size_t col1 = 0; col1 <= 0xFFFF; col1++) {
		for(size_t col2 = 0; col2 < col1; col2++) {
			int a[3], b[3], c[3], d[3], e[3];
			complete_palette(col1, col2, a, b, c, d, e);

			uint32_t total_error[2];
			palette_error(image, a, b, c, d, e, total_error);

			if(total_error[0] < total_error[1]) {
				result = color_res_min(result,
									   (struct color_result) {
										   .error = total_error[0],
										   .col1 = col1,
										   .col2 = col2,
									   });
			} else {
				result = color_res_min(result,
									   (struct color_result) {
										   .error = total_error[1],
										   .col1 = col2,
										   .col2 = col1,
									   });
			}
		}
	}

	colors[0] = result.col1;
	colors[1] = result.col2;
}

// endianness tolerant
void fwrite_32(FILE* f, uint32_t* data, size_t length) {
	for(size_t k = 0; k < length; k++)
		fwrite((uint8_t[]) {data[k] & 0xFF, (data[k] >> 8) & 0xFF,
							(data[k] >> 16) & 0xFF, data[k] >> 24},
			   sizeof(uint8_t), sizeof(uint32_t), f);
}

int main(int argc, char** argv) {
	if(argc < 3) {
		printf("Usage: ./dxtgopt input.png output.dds\n");
		return 1;
	}

	uint8_t* image;
	unsigned width, height;
	if(lodepng_decode32_file(&image, &width, &height, argv[1])) {
		printf("Input image \"%s\" is not a valid image file\n", argv[1]);
		return 1;
	}

	if(width % 4 || height % 4) {
		printf("Input image dimensions (%u,%u) are not a multiple of 4\n",
			   width, height);
		free(image);
		return 1;
	}

	FILE* f = fopen(argv[2], "wb");

	if(!f) {
		printf("Output file could not be created\n");
		free(image);
		return 1;
	}

	size_t total_blocks = width * height / 16;

	fwrite_32(
		f,
		(uint32_t[]) {
			0x20534444, 124, 0x00081007, width,		 height, total_blocks * 8,
			0,			1,	 0,			 0,			 0,		 0,
			0,			0,	 0,			 0,			 0,		 0,
			0,			32,	 0x00000004, 0x31545844, 0,		 0,
			0,			0,	 0,			 0x00001000, 0,		 0,
			0,			0,
		},
		32);

	size_t blocks = 0;

	for(size_t y = 0; y < height; y += 4) {
		for(size_t x = 0; x < width; x += 4) {
			uint32_t tmp[16];
			for(size_t py = 0; py < 4; py++) {
				for(size_t px = 0; px < 4; px++) {
					uint8_t* imgcol = image + (x + px + (y + py) * width) * 4;
					tmp[px + py * 4] = rgb32(imgcol[0], imgcol[1], imgcol[2]);
				}
			}

			uint16_t colors[2];
			uint64_t time_start = time_get();
			compute_16block(tmp, colors);
			uint64_t time_end = time_get();

			uint8_t result[8];
			encode_16block(tmp, colors, result);
			print_block_info(tmp, colors);

			fwrite(result, sizeof(uint8_t), sizeof(result), f);

			printf("took %lu ms\n", time_end - time_start);
			printf("progress: %0.1f%%\n",
				   (float)(++blocks) / (float)total_blocks * 100.0F);
		}
	}

	fclose(f);
	free(image);

	return 0;
}

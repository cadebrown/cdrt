
#include "cdrt.h"

// for the image loading
#define STB_IMAGE_IMPLEMENTATION
#include "./stb_image.h"

img_t cdrt_img_load(char *fname) {
    int width, height, channels;
    pix_t* data = (pix_t*)stbi_load(fname, &width, &height, &channels, STBI_rgb_alpha);
    if (data == NULL) {
        fprintf(stderr, "Failed to load img '%s'\n", fname);
        return NULL;
    }

    img_t self = malloc(sizeof(*self));
    self->data = data;
    self->w = width;
    self->h = height;

    return self;
}

// create an image for data
img_t cdrt_img_create(int w, int h, pix_t* data) {
    img_t self = malloc(sizeof(*self));
    self->w = w;
    self->h = h;
    self->data = malloc(sizeof(*self->data) * w * h);
    if (data != NULL) {
        // copy our data
        memcpy(self->data, data, sizeof(pix_t) * w * h);
    } else {
        // zero out
        int i;
        for (i = 0; i < w * h; ++i) {
            self->data[i] = PIX_RGBA(0, 0, 0, 0);
        }
    }
    return self;
}

// create an image of a single color
img_t cdrt_img_color(pix_t color) {
    return cdrt_img_create(1, 1, (pix_t[]){color});
}

void cdrt_img_free(img_t img) {
    // just free the data
    free(img->data);
    // reset
    img->data = NULL;
    img->w = img->h = 0;
}

void cdrt_img_out(img_t img, char* fname) {
    char* ext = strrchr(fname, '.');

    if (ext == NULL) {
        fprintf(stderr, "Failed to output img to '%s' (no extension)\n", fname);
    } else if (strcmp(ext, ".bmp") == 0) {
        // output bitmap

        FILE* fp = fopen(fname, "wb");
        if (fp == NULL) {
            fprintf(stderr, "Failed to open file '%s' for 'wb' mode\n", fname);
            return;
        }

        // calculate file size
        int filesize = 54 + 3 * img->w * img->h;

        // file header
        unsigned char file_H[14] = {
            'B','M', 
            00, 00, 00, 00, 
            00, 00, 
            00, 00, 
            54, 00, 00, 00
        };

        // info header
        unsigned char info_H[40] = {
            40, 00, 00, 00,
            00, 00, 00, 00,
            00, 00, 00, 00,
            01, 00,
            24, 00,
        };

        // for padding 
        unsigned char pad[4] = {
            00, 00, 00, 00
        };

        // set ints, 1 byte at a time
        file_H[ 2] = (unsigned char)(filesize);
        file_H[ 3] = (unsigned char)(filesize >>  8);
        file_H[ 4] = (unsigned char)(filesize >> 16);
        file_H[ 5] = (unsigned char)(filesize >> 24);

        // code to expand integers into byte format
        info_H[ 4] = (unsigned char)(img->w);
        info_H[ 5] = (unsigned char)(img->w >>  8);
        info_H[ 6] = (unsigned char)(img->w >> 16);
        info_H[ 7] = (unsigned char)(img->w >> 24);

        info_H[ 8] = (unsigned char)(img->h);
        info_H[ 9] = (unsigned char)(img->h >>  8);
        info_H[10] = (unsigned char)(img->h >> 16);
        info_H[11] = (unsigned char)(img->h >> 24);

        // write them both to the file
        fwrite(file_H, 1, sizeof(file_H), fp);
        fwrite(info_H, 1, sizeof(info_H), fp);

        // write in reverse order
        int i, j, k;
        for (j = img->h - 1; j >= 0; j--) {
            for (i = 0; i < img->w; i++) {
                k = i + j * img->w;

                //k = j * img.h + i;
                pix_t cc = img->data[k];

                // instead of RGB, bmp files take BGR
                fwrite(&cc.b, 1, 1, fp);
                fwrite(&cc.g, 1, 1, fp);
                fwrite(&cc.r, 1, 1, fp);
            }
            // pad the file to DWORD len
            fwrite(pad, 1, (4-(img->w*3)%4)%4, fp);
        }
        // close the file pointer
        fclose(fp);

    } else {
        fprintf(stderr, "Failed to output img to '%s' (unknown extension: '%s')\n", fname, ext);
    }

}



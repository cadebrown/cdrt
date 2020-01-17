/* mat.c - constructing materials */

#include "cdrt.h"

mat_t cdrt_mat_new_diffuse(img_t tex_diff) {
    mat_t self = malloc(sizeof(*self));

    self->tex_diff = tex_diff;
    self->tex_spec = NULL;

    return self;
}



mat_t cdrt_mat_new_specular(img_t tex_spec) {
    mat_t self = malloc(sizeof(*self));

    self->tex_diff = NULL;
    self->tex_spec = tex_spec;

    return self;
}



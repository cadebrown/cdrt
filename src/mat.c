/* mat.c - constructing materials */

#include "cdrt.h"

mat_t cdrt_mat_new_diffuse(img_t tex_diff) {
    mat_t self = malloc(sizeof(*self));

    self->tex_diff = tex_diff;
    self->tex_spec = NULL;
    self->tex_refr = NULL;

    self->IOR = IOR_DEFAULT;

    return self;
}

mat_t cdrt_mat_new_specular(img_t tex_spec) {
    mat_t self = malloc(sizeof(*self));

    self->tex_diff = NULL;
    self->tex_spec = tex_spec;
    self->tex_refr = NULL;

    self->IOR = IOR_DEFAULT;

    return self;
}

mat_t cdrt_mat_new_refractive(img_t tex_refr, float IOR) {
    mat_t self = malloc(sizeof(*self));

    self->tex_diff = NULL;
    self->tex_spec = NULL;
    self->tex_refr = tex_refr;

    self->IOR = IOR;

    return self;
}


mat_t cdrt_mat_new_SR(img_t tex_spec, img_t tex_refr, float IOR) {
    mat_t self = malloc(sizeof(*self));

    self->tex_diff = NULL;
    self->tex_spec = tex_spec;
    self->tex_refr = tex_refr;

    self->IOR = IOR;

    return self;
}




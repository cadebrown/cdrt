/* obj.c - for creating objects */

#include "cdrt.h"


obj_t cdrt_obj_new_empty() {
    obj_t self = malloc(sizeof(*self));
    self->type = OBJT_NONE;
    self->T = M4X4_I;
    self->parent = NULL;
    self->n_sub = 0;
    self->sub = NULL;
    return self;
}


obj_t cdrt_obj_new_plane() {
    obj_t self = malloc(sizeof(*self));
    self->type = OBJT_MODEL;
    self->T = M4X4_I;
    self->parent = NULL;
    self->n_sub = 0;
    self->sub = NULL;
    self->model.type = MODEL_PLANE;
    self->model.mat = NULL;
    return self;
}


obj_t cdrt_obj_new_sphere() {
    obj_t self = malloc(sizeof(*self));
    self->type = OBJT_MODEL;
    self->T = M4X4_I;
    self->parent = NULL;
    self->n_sub = 0;
    self->sub = NULL;
    self->model.type = MODEL_SPHERE;
    self->model.mat = NULL;
    return self;
}

obj_t cdrt_obj_new_camera(float FOV) {
    obj_t self = malloc(sizeof(*self));
    self->type = OBJT_CAMERA;
    self->T = M4X4_I;
    self->parent = NULL;
    self->n_sub = 0;
    self->sub = NULL;
    self->camera.FOV = FOV;
    return self;
}

// creates a new directional light
obj_t cdrt_obj_new_dirlight(v3 dir, v3 color) {
    obj_t self = malloc(sizeof(*self));
    self->type = OBJT_DIRLIGHT;
    self->T = M4X4_I;
    self->parent = NULL;
    self->n_sub = 0;
    self->sub = NULL;
    self->dirlight.dir = v3_unit(dir);
    self->dirlight.color = color;
    return self;
}

obj_t cdrt_obj_new_pointlight(v3 color) {
    obj_t self = malloc(sizeof(*self));
    self->type = OBJT_POINTLIGHT;
    self->T = M4X4_I;
    self->parent = NULL;
    self->n_sub = 0;
    self->sub = NULL;
    self->pointlight.color = color;
    return self;
}


void cdrt_obj_addsub(obj_t self, obj_t new_sub) {
    assert(new_sub->parent == NULL);

    int idx = self->n_sub++;
    self->sub = realloc(self->sub, sizeof(*self->sub) * self->n_sub);

    // record the reference
    new_sub->parent = self;
    self->sub[idx] = new_sub;

}

void cdrt_obj_free(obj_t self) {
    int i;
    for (i = 0; i < self->n_sub; ++i) {
        cdrt_obj_free(self->sub[i]);
    }

    free(self->sub);
    free(self);
}




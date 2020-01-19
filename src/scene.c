/* scene.c - scene management */

#include "cdrt.h"

scene_t cdrt_scene_new(obj_t root) {
    scene_t self = malloc(sizeof(*self));
    self->root = root;

    self->cache.cam = NULL;

    self->cache.n_planes = 0;
    self->cache.planes = NULL;
    self->cache.ext_planes = NULL;

    self->cache.n_spheres = 0;
    self->cache.spheres = NULL;
    self->cache.ext_spheres = NULL;

    self->cache.n_dirlights = 0;
    self->cache.dirlights = NULL;

    self->cache.n_pointlights = 0;
    self->cache.pointlights = NULL;

    return self;
}


static void __recache(scene_t self, obj_t cur, int depth, m4x4 cur_T) {
    // record local space change
    cur_T = m4x4_mul(cur_T, cur->T);
    cur->cache.baked_T = cur_T;

    // recurse on all sub objects
    int i;
    for (i = 0; i < cur->n_sub; ++i) {
        __recache(self, cur->sub[i], depth+1, cur_T);
    }

    // now, set special properties
    if (cur->type == OBJT_CAMERA) {
        // set the camera
        assert(self->cache.cam == NULL);
        self->cache.cam = cur;
    } else if (cur->type == OBJT_MODEL) {
        if (cur->model.type == MODEL_SPHERE) {
            int idx = self->cache.n_spheres++;
            self->cache.spheres = realloc(self->cache.spheres, sizeof(*self->cache.spheres) * self->cache.n_spheres);
            self->cache.ext_spheres = realloc(self->cache.ext_spheres, sizeof(*self->cache.ext_spheres) * self->cache.n_spheres);

            // set the values
            self->cache.spheres[idx].xyz = TRANSFORM_POS(cur->cache.baked_T);
            self->cache.spheres[idx].w = 0.6f;
            
            // and set external data
            self->cache.ext_spheres[idx].src = cur;

            // and calculate the correct matrix
            // essentially, create an orthonormal basis of our rows (or close to it)
            // then inverse it
            m3x3 UVR;
            UVR.rows[0] = v3_unit(cur_T.rows[0].xyz);
            UVR.rows[1] = v3_unit(cur_T.rows[1].xyz);
            UVR.rows[2] = v3_unit(cur_T.rows[2].xyz);
            UVR = m3x3_inv(UVR);
            self->cache.ext_spheres[idx].UVR = UVR;


        } else if (cur->model.type == MODEL_PLANE) {
            // add a plane to the cache
            int idx = self->cache.n_planes++;
            self->cache.planes = realloc(self->cache.planes, sizeof(*self->cache.planes) * self->cache.n_planes);
            self->cache.ext_planes = realloc(self->cache.ext_planes, sizeof(*self->cache.ext_planes) * self->cache.n_planes);

            // calculate the normal (get the 'up' direction)
            self->cache.planes[idx].xyz = cur->cache.baked_T.Y.xyz;

            // calculate distance from origin along the normal
            self->cache.planes[idx].w = v3_dot(TRANSFORM_POS(cur->cache.baked_T), self->cache.planes[idx].xyz);

            // create a UV mapping object that can be used to get the tex coords
            m3x3 UVT;
            UVT.X = cur_T.X.xyz;
            UVT.Y = cur_T.Y.xyz;
            UVT.Z = cur_T.Z.xyz;
            UVT = m3x3_inv(UVT);

            // and set external data
            self->cache.ext_planes[idx].src = cur;
            self->cache.ext_planes[idx].UVT = UVT;



        }
    } else if (cur->type == OBJT_DIRLIGHT) {
        // add dirlight to the cache
        int idx = self->cache.n_dirlights++;
        self->cache.dirlights = realloc(self->cache.dirlights, sizeof(*self->cache.dirlights) * self->cache.n_dirlights);

        // set it
        self->cache.dirlights[i] = cur->dirlight;
    } else if (cur->type == OBJT_POINTLIGHT) {
        // add pointlight to the cache
        int idx = self->cache.n_pointlights++;
        self->cache.pointlights = realloc(self->cache.pointlights, sizeof(*self->cache.pointlights) * self->cache.n_pointlights);

        // set it
        self->cache.pointlights[i] = (struct pointlight_cache_s){cur->pointlight, TRANSFORM_POS(cur_T)};
    }
}

// recalculate the cache
void cdrt_scene_recache(scene_t self) {

    self->cache.cam = NULL;
    self->cache.n_pointlights = self->cache.n_dirlights = self->cache.n_planes = self->cache.n_spheres = 0;

    // start the cache recursively, with the identity matrix
    __recache(self, self->root, 0, M4X4_I);
}

/* cdrt.h - header file for the Chemical Development RayTracing library */

#pragma once
#ifndef CDRT_H__
#define CDRT_H__

// standard headers
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

// -lm flag
#include <math.h>

// for timing
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

// include our vector types
#include <cd_vector.h>

/* represents a ray that can be casted */
typedef struct {

    // origin in 3D space
    v3 origin;

    // the direction (which should be a unit directory)
    v3 dir;

} ray_t;

// evaluates a point on a ray, given a distance along it
static inline v3 ray_eval(ray_t ray, float dist) {
    return v3_add(ray.origin, v3_scale(ray.dir, dist));
}


/* a single physical object */
typedef struct obj_s* obj_t;


/* raycasthit - represents a 'hit' on an object */
typedef struct {

    // the object that was hit
    obj_t hit;

    // distance to the hit
    float dist;
    
    // UV, texture coordinates
    v2 uv;

    // the point at which the hit occured
    v3 point;

    // the surface normal at this point which was hit
    v3 normal;

} rch_t;

/* pixel values */
typedef union {
    
    uint8_t _[4];

    struct { uint8_t r, g, b, a; };

} pix_t;


// construct a pixel
#define PIX_RGBA(_r, _g, _b, _a) ((pix_t){_r, _g, _b, _a})

// add 2 pixels
#define PIX_ADD(_a, _b) PIX_RGBA((_a).r+(_b).r, (_a).g+(_b).g, (_a).b+(_b).b, (_a).a+(_b).a)

/* a 2D image */
typedef struct img_t {

    // width and height of the image
    int w, h;

    // the row-major data of the image
    // 01234
    // 56789
    // ...
    pix_t* data;

}* img_t;

// no image
#define IMG_NONE ((struct img_t){-1, -1, NULL})

// load an image from a given file name.
img_t cdrt_img_load(char *fname);

// create an image for data
img_t cdrt_img_create(int w, int h, pix_t* data);

// create an image of a single color
img_t cdrt_img_color(pix_t color);


// output image to a given filename
void cdrt_img_out(img_t img, char* fname);

// free an image's resources
void cdrt_img_free(img_t img);


/* material structure */
typedef struct mat_t {


    // texture for the diffuse coloring
    // you can check whether the data==NULL to see if diffuse is enabled
    img_t tex_diff;

    // texture for the specular coloring
    // you can check whether the data==NULL to see if reflection is enabled
    img_t tex_spec;

    
}* mat_t;


// empty material
#define MAT_EMPTY ((struct mat_t){NULL, NULL})


// create a new diffuse material
mat_t cdrt_mat_new_diffuse(img_t tex_diff);

// create a new specular material
mat_t cdrt_mat_new_specular(img_t tex_spec);



enum {
    // none-type, i.e. is not rendered, and may be used as a spacer
    OBJT_NONE,

    // camera-type, i.e. is the camera that the scene will be rendered with
    OBJT_CAMERA,

    // point light, from a single location
    OBJT_POINTLIGHT,

    // directional light, i.e. like a sun that is infinite distance away
    OBJT_DIRLIGHT,

    // model-type, i.e. is a rendered geometry within a scene
    OBJT_MODEL

};

struct obj_s {

    // current local transform space
    m4x4 T;

    // the type (OBJT_*) of the object
    int type;

    union {

        // if type==OBJT_CAMERA, this describes it
        struct obj_camera {
            
            // field of view, in radians
            float FOV;

        } camera;

        // if type==OBJT_MODEL
        struct obj_model {

            // the model type
            enum {

                // means this represents a plane.
                // in that case, its transform's `UP` represents the normal vector
                // and it is infinite
                MODEL_PLANE = 0,

                // means this represents a sphere
                // in that case, the transform's position is the center of the spheres,
                // and the scales/etc are the stretch factors of it
                MODEL_SPHERE,

            } type;

            // the model's material
            mat_t mat;

        } model;

        // if type==OBJT_POINTLIGHT
        struct pointlight_s {

            v3 color;
        
        } pointlight;

        // if type==OBJT_DIRLIGHT
        struct dirlight_s {

            // the (unit vector) direction of the light
            v3 dir;

            // the (rgb) color of the directional light
            v3 color;

        } dirlight;
    };

    // the pointer to the parent object
    obj_t parent;

    // number of sub objects 
    int n_sub;

    // the sub object
    obj_t* sub;


    // internal cache computed by scene_recache
    struct {
        // the baked transform from all the items
        m4x4 baked_T;
    } cache;

};

// get teh position from a transform
#define TRANSFORM_POS(_tra) V3((_tra).r0c3, (_tra).r1c3, (_tra).r2c3)

// creates a new empty object
obj_t cdrt_obj_new_empty();

// creates a new plane object
obj_t cdrt_obj_new_plane();

// creates a new sphere object
obj_t cdrt_obj_new_sphere();

// creates a new camera object
obj_t cdrt_obj_new_camera(float FOV);

// creates a new directional light
obj_t cdrt_obj_new_dirlight(v3 dir, v3 color);

// creates a new point light, given a color
obj_t cdrt_obj_new_pointlight(v3 color);


// add a new sub object, which must not have a parent
void cdrt_obj_addsub(obj_t self, obj_t new_sub);

// frees an object and its resources
void cdrt_obj_free(obj_t self);


/* a rendering scene */

typedef struct scene_t {

    // the root object of the scene
    obj_t root;


    // cached values to be used internally by the renderer
    struct {

        // the camera to use
        obj_t cam;


        // number of planes
        int n_planes;

        // packed result of plane
        // xyz: normal, w: dot(centerxyz, normal) aka distance off of starting point
        v4* planes;

        // extra data about the planes
        struct {

            // source object
            obj_t src;

            // translation matrix to map UV coordinates
            m3x3 UVT;

        }* ext_planes;


        // number of spheres
        int n_spheres;

        // packed result of spheres:
        // xyz -> position, w -> radius^2
        v4* spheres;

        // extra data about the sphere
        struct {
            // the source object
            obj_t src;

            // rotation matrix to map UV coordinates
            m3x3 UVR;
            
        }* ext_spheres;


        // number of dirlights in the scene
        int n_dirlights;

        // packed array of dirlights
        struct dirlight_s* dirlights;

        // number of pointlights in the scene
        int n_pointlights;

        // packed array of pointlights
        struct pointlight_cache_s {
            // the light
            struct pointlight_s ptl;

            // the position
            v3 pos;
        }* pointlights;

    } cache;


}* scene_t;


// create a new scene object 
scene_t cdrt_scene_new(obj_t root);

// recalculate cache results for a scene
void cdrt_scene_recache(scene_t self);


// main rendering function, to output to an image
void cdrt_render(img_t dest, scene_t scene);

#endif


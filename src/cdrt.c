
#include "cdrt.h"


int main(int argc, char** argv) {

    /* load images/textures */

    img_t t_map = cdrt_img_load("./worldmap.png");
    if (!t_map) return 1;

    img_t t_chk = cdrt_img_load("./chk.png");
    if (!t_chk) return 1;

    img_t t_red = cdrt_img_color(PIX_RGBA(255, 0, 0, 0));
    if (!t_red) return 1;

    /* construct our materials we will use */

    // material of the entire world texture, mercator projection
    mat_t m_map = cdrt_mat_new_specular(t_map);

    // material of a checkerboard pattern
    mat_t m_chk = cdrt_mat_new_diffuse(t_chk);

    // material of red
    //mat_t m_red = cdrt_mat_new_SR(t_red, t_red, 2.0f);
    mat_t m_red = cdrt_mat_new_specular(t_red);


    /* construct object hierarcy */

    // the (empty) root of the entire scene
    obj_t root = cdrt_obj_new_empty();

    // create a camera with an FOV of 70deg
    obj_t cam = cdrt_obj_new_camera(70.0f);

    // the ground, which is a flat, infinite plane
    obj_t ground = cdrt_obj_new_plane();
    ground->model.mat = m_chk;
    ground->T = m4T_xyz(0.0f, -1.0f, 0.0f);
    //ground->T = m4x4_mul(m4T_rot_z(3.1415926535*.5), m4T_xyz(0.0f, 4.0f, 0.0f));
    
    printf("MAT: " M4X4_FMT "\n", M4X4__(ground->T));

    // create 2 spheres

    obj_t s0 = cdrt_obj_new_sphere();
    s0->model.mat = m_map;
    s0->T = m4T_xyz(100.0f, 0.5f, 2.0f);

    obj_t s1 = cdrt_obj_new_sphere();
    s1->model.mat = m_red;
    s1->T = m4T_xyz(0.0f, 0.1f, 1.5f);
    
    // create the sun light
    obj_t sun = cdrt_obj_new_dirlight(V3(1, -1, 1), V3_1);

    // create a point light
    obj_t ptl = cdrt_obj_new_pointlight(V3(4.4, 4.2, 4.3));
    ptl->T = m4T_xyz(0.0f, -.7f, 1.1f);

    // construct hierarchy
    cdrt_obj_addsub(root, cam);

    cdrt_obj_addsub(root, ground);

    cdrt_obj_addsub(root, s0);
    cdrt_obj_addsub(root, s1);

    cdrt_obj_addsub(root, sun);

    cdrt_obj_addsub(root, ptl);


    // create a scene containing the root of the hierarchy
    scene_t scn = cdrt_scene_new(root);

    // create blank output
    img_t output = cdrt_img_create(1024, 1024, NULL);

    // render to image
    cdrt_render(output, scn);

    // output the image
    cdrt_img_out(output, "out.bmp");

    return 0;
}


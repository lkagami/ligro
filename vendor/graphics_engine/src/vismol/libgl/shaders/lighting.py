#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Shared lighting block for every lit fragment shader.
#
#  [EN] 2026-10-13: the ambient/diffuse/specular math used to be copied
#  into each fragment shader (spheres, sticks, impostors, surface, cartoon),
#  all reading one `uniform Light my_light`. It now lives here, once. A
#  fragment shader opts in by writing the line
#
#      #include <vm_lighting>
#
#  among its declarations and calling vm_lighting() from main();
#  VismolGLCore.create_shader() swaps that line for LIGHTING_GLSL before
#  compiling (GLSL itself has no #include, so a shader that is compiled
#  some other way fails loudly instead of silently losing its lighting).
#  The uniforms are uploaded by VismolGLCore.load_lights().
#
#  Frame convention (unchanged): normals are transformed by model_mat only
#  and the view matrix carries no rotation, so light directions are fixed
#  relative to the CAMERA -- the lights do not turn with the molecule.
#
#  With a single light and vm_shading_mode == 0 the result is identical to
#  the previous per-shader code.

LIGHTING_INCLUDE = "#include <vm_lighting>"

MAX_LIGHTS = 4

SHADING_MODES = ("standard", "flat", "toon", "matte", "glossy", "rim")
# gl_parameters "shading_mode" -> vm_shading_mode uniform (3 is taken by
# NORMAL_MAP_MODE below)
SHADING_MODE_IDS = {"standard": 0, "flat": 1, "toon": 2, "matte": 4, "glossy": 5, "rim": 6}
# Internal vm_shading_mode value: the shaders output normals instead of
# colours (VismolGLCore._read_solid_depth()). Not selectable by the user.
NORMAL_MAP_MODE = 3

LIGHTING_GLSL = """
#define VM_MAX_LIGHTS 4
uniform int   vm_n_lights;
uniform vec3  vm_light_dir[VM_MAX_LIGHTS];  // unit vectors pointing TOWARDS each light
uniform vec3  vm_light_int[VM_MAX_LIGHTS];
uniform vec3  vm_ambient;                   // ambient_coef * key-light intensity
uniform float vm_shininess;
uniform int   vm_shading_mode;              // 0 standard, 1 flat, 2 toon, 3 normal map,
                                            // 4 matte, 5 glossy, 6 rim
uniform float vm_flat_softness;
uniform float vm_toon_bands;

// normal      -- unit surface normal
// cam_to_frag -- unit vector from the camera to the fragment
// col         -- base colour
// wrap        -- 0.0 = Lambert, 1.0 = half-Lambert (surfaces)
vec3 vm_lighting(vec3 normal, vec3 cam_to_frag, vec3 col, float wrap){
    if (vm_shading_mode == 3) {
        // Not a look: writes the surface normal as a colour, for the
        // outline pass of the image export (NORMAL_MAP_MODE below).
        return normal * 0.5 + 0.5;
    }
    if (vm_shading_mode == 1) {
        // Flat: the base colour, only slightly darkened towards the
        // silhouette so shapes still read as round. Ignores the lights.
        float facing = max(0.0, dot(normal, -cam_to_frag));
        return col * (1.0 - vm_flat_softness * (1.0 - facing));
    }
    // Matte: soft wrapped diffuse and no highlight at all.
    if (vm_shading_mode == 4) wrap = max(wrap, 0.6);
    vec3 lit = vm_ambient * col;
    for (int i = 0; i < vm_n_lights; i++) {
        float diffuse_coef = max(0.0, dot(normal, vm_light_dir[i]));
        diffuse_coef = mix(diffuse_coef, (diffuse_coef + 1.0) * 0.5, wrap);
        float specular_coef = 0.0;
        if (diffuse_coef > 0.0 && vm_shading_mode != 4)
            specular_coef = pow(max(0.0, dot(cam_to_frag, reflect(vm_light_dir[i], normal))), vm_shininess);
        if (vm_shading_mode == 5) {
            // Glossy: a small, bright, white highlight (polished plastic)
            // on top of the normal diffuse term.
            lit += diffuse_coef * col * vm_light_int[i]
                 + 1.6 * pow(specular_coef, 5.0) * vm_light_int[i];
            continue;
        }
        if (vm_shading_mode == 2) {
            // Toon: a few hard diffuse steps and a hard-edged highlight.
            diffuse_coef = ceil(diffuse_coef * vm_toon_bands) / vm_toon_bands;
            specular_coef = step(0.5, specular_coef);
        }
        vec3 diffuse = diffuse_coef * col * vm_light_int[i];
        lit += diffuse + specular_coef * vm_light_int[i] * (vec3(1.0) - diffuse);
    }
    if (vm_shading_mode == 6) {
        // Rim: brightens surfaces seen edge-on, so silhouettes stand out
        // (most useful on a dark background).
        float rim = pow(1.0 - max(0.0, dot(normal, -cam_to_frag)), 3.0);
        lit += 0.7 * rim * (0.5 * col + vec3(0.5));
    }
    return lit;
}
"""

# Extra lights added to the key light (gl_parameters "light_position" /
# "light_intensity") by each "light_preset". Directions point towards the
# light, in the camera-fixed frame: +x right, +y up, +z towards the viewer.
LIGHT_PRESETS = {
    "headlight"  : [],
    "three_point": [{"direction": [ 0.8, 0.1, 0.6], "intensity": [0.25, 0.25, 0.25]},
                    {"direction": [ 0.0, 0.6,-0.8], "intensity": [0.30, 0.30, 0.30]}],
}
# Key-light direction used by presets that move it off the camera axis
# ("headlight" keeps gl_parameters["light_position"]).
PRESET_KEY_DIRECTION = {
    "three_point": [-0.45, 0.50, 0.74],
}


def inject(shader_source):
    """ Returns shader_source with the include line replaced by the GLSL
        lighting block (unchanged if the shader does not use it). """
    return shader_source.replace(LIGHTING_INCLUDE, LIGHTING_GLSL)

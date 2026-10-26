#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Independent Cartoon representation -- shader pair.
#
#  [EN] Freshly written for the new ribbon renderer (see ribbon_geometry.py/
#  ribbon_backbone.py in vismol/utils/ and CartoonRepresentation in
#  libgl/representations.py) -- replaces the shader that shipped with the
#  cartoon implementation being removed. `CameraMatrices` (the UBO block:
#  view_mat/proj_mat), the `model_mat` uniform, and the shared lighting
#  block (shaders/lighting.py, `#include <vm_lighting>`) are NOT a
#  cartoon-specific convention to reinvent -- they're this whole engine's
#  own fixed per-frame binding contract (VismolGLCore.update_camera_ubo()/
#  load_lights(), shared by every representation type's shader), so this
#  new shader still declares them exactly as every other shader in this
#  package does; independence from the old cartoon code is about the
#  RIBBON algorithm and this file's own structure, not about redeclaring
#  the graphics engine's shared interface under different names.
# ============================================================================

v_shader_triangles = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
uniform mat4 model_mat;

in vec3 vert_coord;
in vec3 vert_color;
in vec3 vert_norm;

out vec3 frag_color;
out vec3 frag_coord;
out vec3 frag_norm;

void main () {
    frag_color = vert_color;
    // [EN] Matches every other shader in this package's own convention
    // (see e.g. shaders/sticks.py): normals are transformed by model_mat
    // alone, not the full model-view matrix -- the light directions
    // (vm_light_dir, see shaders/lighting.py) are likewise fixed,
    // un-transformed directions (VismolGLCore.load_lights(), no per-frame
    // view transform), so keeping normals in that SAME space is what makes
    // the dot-product lighting math correct here, not an inconsistency.
    frag_norm = mat3(transpose(inverse(model_mat))) * vert_norm;
    frag_coord = (view_mat * model_mat * vec4(vert_coord, 1.0)).xyz;
    gl_Position = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
}
"""

f_shader_triangles = """
#version 330

#include <vm_lighting>
// [EN] 2026-09-26 user request ("alterar a transparencia da representacao
// de cartoon"): plain alpha uniform, same role as SurfaceRepresentation's
// own "surf_alpha" -- see CartoonRepresentation.draw_representation() for
// where this gets set + the matching glBlendFunc/glDepthMask dance.
uniform float cartoon_alpha;

// [EN] 2026-09-28 user request ("a representacao deve ser sensivel ao
// fog... melhora muito a nocao de profundidade"): same linear depth-fog
// uniforms every other lit shader in this package already declares (see
// shaders/sticks.py) -- VismolGLCore.load_fog() was ALREADY being called
// for this program by CartoonRepresentation.draw_representation(), but
// silently did nothing (_uniform_changed() skips a uniform the program
// doesn't declare), so the ribbon was the one thing on screen that
// ignored the fog.
uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;
uniform float fog_intensity;

in vec3 frag_coord;
in vec3 frag_color;
in vec3 frag_norm;

out vec4 final_color;

void main () {
    vec3 normal = normalize(frag_norm);

    // [EN] Two-sided lighting: a ribbon is a THIN, open surface (helix/
    // strand cross-sections are only ~0.5 A tall) -- draw_representation()
    // disables face culling for it, so both the "front" and "back" of the
    // SAME physical ribbon surface are visible at once from many camera
    // angles. Without this flip, whichever side's stored normal happens
    // to face away from the camera/light renders dark while the
    // neighbouring, front-facing side of that same continuous surface
    // stays lit -- alternating light/dark blotches on one ribbon.
    // gl_FrontFacing is true for whichever side the current triangle's
    // winding order faces the camera; flipping the normal on the back
    // side keeps the same physical surface point lit consistently no
    // matter which side the camera currently sees.
    if (!gl_FrontFacing) {
        normal = -normal;
    }

    vec3 lit_color = vm_lighting(normal, normalize(frag_coord), frag_color, 0.0);

    // [EN] frag_coord is already in view space (see v_shader_triangles),
    // so |z| is the eye distance -- same formula as shaders/sticks.py,
    // except the factor is clamped to [0,1] (beyond fog_end the unclamped
    // mix() would extrapolate PAST the fog colour). Only RGB is fogged;
    // alpha stays cartoon_alpha so transparency is unaffected.
    float dist = abs(frag_coord.z);
    if (dist >= fog_start) {
        float fog_factor = 1.0 - fog_intensity * (1.0 - clamp((fog_end - dist) / (fog_end - fog_start), 0.0, 1.0));
        lit_color = mix(fog_color.rgb, lit_color, fog_factor);
    }

    final_color = vec4(lit_color, cartoon_alpha);
}
"""

# [EN] 2026-09-26 user request ("tornar a representacao clicavel... plot o
# c-alpha trace como referencia para encontrar os residuos"): a REAL
# picking shader for Cartoon -- the "cartoon_sel" program used to be
# compiled from the SAME v_shader_triangles/f_shader_triangles pair as the
# visible pass (see _compile_shader_cartoon() in vismol_glcore.py before
# this fix), which is wrong for picking: that fragment shader applies real
# lighting, so the picking-ID colours it would output are NOT the exact
# byte-for-byte IDs _pick() decodes (GL.glReadPixels + r+g*256+b*256*256).
# This pair is flat/unlit -- outputs vert_color completely unchanged, same
# convention as every other representation's own "_sel" shader (compare
# shaders/lines.py's sel_vertex_shader_lines/sel_fragment_shader_lines).
# Draws CartoonRepresentation's own separate Ca-trace pick geometry (see
# _build_pick_trace_geometry() in representations.py), NOT the visible
# ribbon mesh -- a ribbon has no 1:1 vertex-to-atom correspondence to give
# real per-vertex pick colours to, but its OWN backbone Ca trace does.
sel_vertex_shader_cartoon = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
uniform mat4 model_mat;

in vec3 vert_coord;
in vec3 vert_color;

out vec3 frag_color;

void main () {
    frag_color = vert_color;
    gl_Position = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
}
"""

sel_fragment_shader_cartoon = """
#version 330

in vec3 frag_color;
out vec4 final_color;

void main () {
    final_color = vec4(frag_color, 1.0);
}
"""

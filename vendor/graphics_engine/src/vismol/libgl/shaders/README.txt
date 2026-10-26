EasyHybrid3 — superfície de orbitais como triângulos (flat shading)
=====================================================================

Contexto: no branch dev (src/graphics_engine como submódulo), a
representação de superfície (SurfaceRepresentation) desenhava a malha
como GL_LINES e o shader "surface" estava carregando o par de shaders
de linhas em vez do par de triângulos com Phong que já existe em
shaders/surface.py. O geometry shader de triângulos também dependia
de um atributo "vert_normal" que nunca é alimentado pela VAO (não há
buffer de normais criado em nenhum lugar do código).

O que este patch faz:
  1. shaders/surface.py — geometry_shader_surface passa a usar a normal
     de face (norm_vec, já calculada via produto vetorial das arestas
     do triângulo) em vez do atributo vert_normal inexistente. Resultado:
     flat shading (uma normal constante por triângulo).
  2. vismol_glcore.py — _compile_shader_surface() agora carrega
     vertex_shader_surface / fragment_shader_surface / geometry_shader_surface
     em vez do trio "_lines".
  3. representations.py — draw_representation() da SurfaceRepresentation
     desenha com GL_TRIANGLES em vez de GL_LINES.

Isso NÃO implementa smooth shading (normal por vértice) — para isso
seria preciso calcular normais no momento da geração da malha
(surface_analysis_window.py / pDynamo) e criar um buffer de normais
novo na VAO. Fica como próximo passo, se o flat shading não for
suficiente visualmente.

Como aplicar
------------
Opção A — patch (recomendado, preserva o resto do arquivo intacto):

    cd EasyHybrid3/src/graphics_engine
    git apply /caminho/para/surface_triangles.patch
    git diff --stat     # conferir que só os 3 arquivos mudaram

Se "git apply" reclamar de contexto (por ter mudado desde o clone
usado aqui), use:

    patch -p1 < /caminho/para/surface_triangles.patch

Opção B — substituir os arquivos inteiros manualmente:
  Copie os 3 arquivos de full_files/ para:
    src/graphics_engine/src/vismol/libgl/shaders/surface.py
    src/graphics_engine/src/vismol/libgl/vismol_glcore.py
    src/graphics_engine/src/vismol/libgl/representations.py

Importante: este ambiente é headless (sem display/GPU), então o patch
foi validado apenas por análise estática (sintaxe Python OK, lógica do
GLSL consistente com o resto do pipeline). O teste visual real só pode
ser feito na sua máquina, rodando o EasyHybrid de fato.

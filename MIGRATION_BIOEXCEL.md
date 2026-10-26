# LiGRO — migração inicial da dinâmica molecular para BioExcel

## Alterações nesta versão

- Criado `ligro/core/pipelines/dynamics.py`, que coordena NVT, NPT e MD de produção por meio da interface `Runner`.
- O pipeline chama `BiobbRunner`, que delega `grompp` e `mdrun` aos BioExcel Building Blocks em `ligro/biobb_gmx.py`.
- `save_tprfile()` prepara NVT/NPT e gera `md.tpr` sem executar a produção.
- `run_simulation()` executa o pipeline, valida os resultados e copia os arquivos existentes para `<projeto>_results`.
- O resultado `md.xtc` é gerado como saída de produção; não se tenta renomear `md_2.xtc` inexistente.
- Adicionado `pack_topology()` para criar um arquivo de topologia zipado atualizado com o `trp.top` e arquivos `.itp` locais.
- Adicionados testes unitários sem execução de GROMACS real.

## Validação

A suíte unitária foi executada no ambiente de trabalho: 40 testes passaram. A compilação de todos os módulos Python também passou. Isso não substitui uma simulação científica real; os wrappers BioExcel precisam ser verificados no ambiente do usuário com a versão instalada de `biobb-gromacs`.

## Limite importante desta migração

Esta é uma migração inicial do estágio NVT/NPT/MD e da minimização da proteína que já estava em `ligro/core`. Os caminhos legados de preparação proteína-ligante/cofator ainda contêm chamadas diretas ao GROMACS para solvatação/ionização/minimização e geração de restraints; os fluxos LIE também ainda usam `queue.sh`/`queue2.sh`. Portanto, não se deve considerar a migração integral de todos os workflows concluída ainda.

O próximo passo necessário é portar a preparação específica de ligantes/cofatores e os workflows LIE para Building Blocks equivalentes, validar as interfaces reais de `genrestr`/`trjconv`, e só então remover as filas legadas restantes sem alterar a semântica científica.

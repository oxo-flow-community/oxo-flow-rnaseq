# STAR index (auto-built)

The STAR index is an upstream input artifact (PREPARE_GENOME is not
ported), so the repo ships no binary index. The engine's `star_index`
reference builder (see `[[references]]` in main.oxoflow) builds it from
the shipped fixture genome when the `SAindex` sentinel is missing — a
fresh clone runs end to end.

For real data, point `config.star_index` at your own index directory
(containing `SA`, `SAindex`, `genomeParameters.txt`, ...); an existing
index wins over the auto-build.

# Reproducir los candidatos de fósforo

Desde la raíz del repositorio:

```bash
Rscript --vanilla tests/test_phosphorus_indicators.R
Rscript --vanilla analysis/run_phosphorus_indicators.R
(cd data/derived && sha256sum -c SHA256SUMS)
(cd data/phosphorus_indicators && sha256sum -c SHA256SUMS)
```

Se reutilizan `dplyr`, `ggplot2`, `patchwork` y `scales`, ya usados por el proyecto.
PNG/PDF/SVG usan dispositivos de R; se necesita Cairo para PDF/SVG.
CI analiza el snapshot sin `arrow` ni acceso a las grandes fuentes locales.

Para regenerar la extracción, con `arrow` instalado:

```bash
Rscript --vanilla analysis/prepare_phosphorus_inputs.R /ruta/a/la/fuente/data
```

Los ocho archivos originales necesarios y sus hashes están en
`data/phosphorus_indicators/sources.csv`. Se leen TPM anotados por biblioteca,
EggNOG completo, taxonomía y etiquetas de muestras. No se modifican las fuentes.
El snapshot conserva 717 IDs (505 procariotas y 212 eucariotas) en 6.741
observaciones positivas. El número de IDs no equivale a organismos o copias
génicas únicas: puede incluir redundancia entre ensamblajes/transcritos.

## Entradas y trazabilidad

| Archivo en `data/phosphorus_indicators/` | Contenido |
| --- | --- |
| `transcripts.csv` | ID original, compartimento, tratamiento, tiempo, biblioteca técnica, TPM, linaje, KO procesados/completos, PFAM y clasificación conservadora. |
| `original_annotation_hits.csv` | Hit EggNOG original con referencia, e-value, score, descripción, nivel de anotación y todos los KO. |
| `original_taxonomy.csv` | Asignación original del clasificador, sin reinterpretar sus encabezados eucarióticos. |
| `taxon_totals.csv` | RNA de todos los genes de los taxones incluidos, antes del filtro funcional. |
| `library_qc.csv`, `sample_metadata.csv` | Totales de bibliotecas y etiquetas WP2 originales. R1/R2/R3 no implican replicación biológica. |
| `sources.csv`, `SHA256SUMS` | Identidad de fuentes y snapshot. `rows_wp2` cuenta filas WP2 de expresión/metadatos y filas seleccionadas de anotación, según `role`. |

Los predictores provienen de `data/derived`: se conservan las cifras exactas de
Figura 6 para NP, PO4, DIN y PP; Figura 1 añade otros ocho predictores disponibles.
No se reconstruye N:P dividiendo medias ni se imputan límites de detección.

## Métodos y parámetros

Selección por tokens KO exactos de la lista **completa** del TSV EggNOG:
`pstA=K02038`, `pstB=K02036`, `pstC=K02037`, `pstS=K02040`,
`pit=K03306/K14640`, `phoA=K01077`. K14640 mantiene la inclusión del catálogo
histórico y lleva una bandera. K01077 es una familia amplia de fosfatasas
alcalinas: `phoA` es una etiqueta operativa, especialmente en eucariotas.

`strict` excluye IDs con más de un KO distinto; es una sensibilidad conservadora,
no un umbral de confianza validado. `inclusive` admite esas asignaciones una
vez por gen compatible: no deben sumarse genes diferentes como un total sin
solapamientos. Las subunidades `pst` se mantienen separadas. La razón auxiliar
`pst_pit` suma explícitamente sus TPM medios y divide por `pit`; un denominador
cero/ausente produce NA y no se añaden pseudocuentas.

Los TPM ya suman un millón por biblioteca y compartimento. No se renormaliza
el subconjunto de genes de P. Las omisiones de la tabla positiva dispersa se
completan como no detecciones operativas en las bibliotecas existentes, no
como ausencia biológica. No se crea la biblioteca eucariota R3 de R+P 0 h.
Se suman ORFs por taxón/gen/biblioteca y se promedian bibliotecas técnicas por
estado incluyendo no detecciones; se conservan mínimo, máximo y detección.

La contribución usa el TPM del gen de **todos** los taxones del mismo gen,
compartimento y biblioteca, incluidos grupos ajenos al objetivo y sin resolver.
La fracción dentro del taxón divide por todo su RNA, incluidos otros genes.
Se promedian las fracciones por biblioteca; no son cocientes de medias ni
actividad por célula. Denominador sin RNA produce NA. Sin DNA o abundancia
celular específica no se separa regulación de cambios poblacionales.

Modelos nuevos: `respuesta ~ predictor`, sin transformar, para TPM y fracción
dentro del taxón, con seis estados. Se exportan pendiente, intercepto, R²,
rho de Spearman y pendientes al retirar cada estado. Estabilidad significa
conservación del signo, no de magnitud. Los p suponen independencia y se
etiquetan como inadecuados para inferencia de tratamiento. BH cubre todas las
combinaciones taxón–gen–predictor–normalización por modo. No se selecciona por
el mejor p ni se ajustan modelos multivariantes a seis estados.

La regla exploratoria por gen pide detección en ≥80% de bibliotecas, pendiente
negativa de la fracción con PO4, signo conservado al retirar cualquier estado
y R>R+P en ambos tiempos. No prueba causalidad, especificidad o significación.
Familias/órdenes/clases/grupos no se presentan como organismos resueltos.
El ranking por taxón prefiere `pstS` (luego `pstB/C/A`, finalmente `pit`) en
bacterias y `phoA` en eucariotas. Prioriza coherencia en subunidades, detección
y finalmente contribución. Los nombres GTDB provisionales conservan esa condición.

## Salidas

Todas se escriben en `results/phosphorus_indicators/`:

- `candidate_taxa_*.csv`: ranking por taxón y grupo.
- `candidate_ranking_*.csv`: evidencia por taxón/gen e IDs originales.
- `taxon_gene_{libraries,states}_*.csv`: métricas antes/después de agregar técnicos.
- `associations_*.csv`: todas las asociaciones, p/q exploratorios y sensibilidad.
- `contrasts_*.csv`, `pst_pit_by_state_*.csv`: contrastes descriptivos y razones.
- `annotation_sensitivity.csv`: comparación de modos sin ocultar ambigüedad.
- `historical_regression_audit.csv`: métricas LM históricas y sus escalas.
- `figures/`: 17 figuras, cada una en PNG de 300 dpi, PDF y SVG.

Se dibujan hasta cuatro taxones por grupo del ranking estricto para legibilidad;
las tablas conservan todos. Las rectas no muestran intervalos inferenciales.
Los rangos de condiciones son mínimo–máximo técnico, no confianza biológica.
Los históricos permanecen intactos; véase la [auditoría](phosphorus_audit.md).

## Datos adicionales necesarios

Estos archivos se incorporarían en `data/phosphorus_indicators/` como nuevas
entradas revisadas; el módulo actual no afirma leerlos automáticamente:

| Archivo propuesto | Columnas mínimas | Uso |
| --- | --- | --- |
| `taxonomy_verified.csv` | `compartment`, `feature_id`, `taxon_id`, `accepted_name`, `accepted_rank`, `lineage`, `reference_accession`, `database_version`, `method`, `identity`, `coverage`, `support`, `review_status` | Resolver rangos eucarióticos y género dinoflagelado con confianza evaluable. |
| `sample_design_verified.csv` | `compartment`, `library_id`, `source_sample_id`, `bag_id`, `treatment`, `time_h`, `technical_replicate`, `biological_replicate`, `batch`, `collection_datetime` | Resolver etiquetas contradictorias sin inventar replicación. |
| `functional_validation.csv` | `feature_id`, `sequence_accession`, `amino_acid_sequence`, `domain_hits`, `active_site_support`, `phylogeny_support`, `localization`, `assay`, `validation_status` | Validar función y organismo de origen; descartar contaminación. |
| `population_controls.csv` | `source_sample_id`, `taxon_id`, `dna_gene_copies` o `cells_per_ml`, `method`, `units`, `detection_limit`, `uncertainty` | Separar población e inversión de RNA con muestras emparejadas. |
| `independent_validation.csv` | `experiment_id`, `bag_or_culture_id`, `taxon_id`, `P_treatment`, `time_h`, `PO4`, `DIN`, `NP_definition`, `marker_measurement`, `growth_or_P_response`, `units`, `detection_limits` | Validar respuesta/especificidad en unidades independientes y frente a otros estresores. |

El último requisito exige nuevos experimentos/muestras. Las secuencias pueden
recuperarse de ensamblajes originales o reprocesando ENA PRJEB94162; no se
inventaron ni se reensamblaron en esta ampliación.

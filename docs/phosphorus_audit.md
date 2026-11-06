# Auditoría de INTERES/WP2 para candidatos indicadores de fósforo

Fecha: 6 de noviembre de 2026. La auditoría se realizó antes de implementar el
módulo nuevo. Se conservaron sin cambios los datos derivados, las figuras y
`analysis/paper_figures.R` de la rama principal.

## Fuentes y alcance de la reproducción

El proyecto correcto es
[erickdelgadillo/dealing-with-p-metat](https://github.com/erickdelgadillo/dealing-with-p-metat),
asociado al [artículo de Delgadillo-Nuño y colaboradores](https://doi.org/10.1093/ismeco/ycag035).
El punto de partida auditado es el commit
[`6a3d0c927d16515f684e27d675cad58bc5dff60d`](https://github.com/erickdelgadillo/dealing-with-p-metat/tree/6a3d0c927d16515f684e27d675cad58bc5dff60d).
La documentación de `data/README.md` identifica el origen de las tablas compactas
y sus diferencias conocidas. El cuaderno histórico se conserva en
[`ccffdc1eb624ae15052bd612df6069e236712b30`](https://github.com/erickdelgadillo/dealing-with-p-metat/blob/ccffdc1eb624ae15052bd612df6069e236712b30/legacy/original/Script_WP2_Plots_V12.Rmd).

La rama principal conserva entradas suficientes para reconstruir gráficamente
las figuras 1–6 y los suplementos 1–2. Esa reproducción empieza en tablas
derivadas: no repite ensamblaje, cuantificación, anotación ni modelos DGE
originales. Las tablas de contribución de la Figura 3 agrupan categorías
funcionales; no permiten atribuir por sí solas `pst`, `pit` o `phoA` a un género.
La Figura 6 conserva seis medias por serie, sin identificadores de transcritos
ni bibliotecas individuales. Tampoco permite desagregar organismos por sí sola.

Para la ampliación se recuperó una fuente local de TPM anotados por biblioteca.
El snapshot `data/phosphorus_indicators/` conserva identificadores originales,
las tres columnas KO procesadas y todos los KO del EggNOG original,
la taxonomía original, las etiquetas de muestras, los
totales por taxón y los totales de biblioteca. `sources.csv` registra nombres
relativos y SHA-256 de las fuentes; `SHA256SUMS` identifica los archivos
compactos. Su regeneración requiere los archivos originales documentados; el
análisis del snapshot no requiere esas rutas locales ni `arrow`. La selección
y la bandera de ambigüedad consultan la lista completa de KO original.
Los hits funcionales y las filas originales de taxonomía se guardan aparte.

La copia de trabajo local partía de `7f38fcd254937a01db9ec5f197abe1b80b5fa3cb`
(merge del PR 15). El conector devolvió la base remota indicada arriba.
Se compararon los identificadores de blob: el script histórico, las doce
entradas derivadas y las ocho figuras comunes son idénticos en ambas bases.
La publicación añade archivos nuevos conservando la base remota.

## Cómo se hicieron las regresiones históricas

El cuaderno histórico lee hojas de `Correlations.xlsx`, en lugar de reconstruir
esas hojas desde los transcritos. Las seis observaciones son C, R y R+P en 0 y
72 h. La variable `pst_pit` representa una razón agregada, no la expresión
individual de `pst` o `pit`.

| Panel/hoja | Respuesta y predictor | Curva histórica | Modelo del que se extrae p |
| --- | --- | --- | --- |
| A / `High_affinity` | `pst:pit` frente a PO4, por orden bacteriano | NLS: `y = a*x^b` | `lm(log(y) ~ log(x))` |
| B / `Low_affinity` | `nptA` de diatomeas frente a PO4 | Potencia retransformada del ajuste log-log | `lm(log(y) ~ log(x))` |
| C / `P_mobilization` | `ppA` frente a DIN | NLS: `y = a*exp(b*x)` | `lm(log(y) ~ x)` |
| D / `Organic_P` | `phoA` de diatomeas/Pseudomonadales frente a N:P | Lineal | `lm(y ~ x)` |
| E / `Organic_p_and_Phosphonate` | Genes de dinoflagelados frente a producción primaria | NLS exponencial | `lm(log(y) ~ x)` |
| F / `Phosphonate` | `phnC` de Rhodobacterales frente a producción primaria | Lineal | `lm(y ~ x)` |

En particular, la asociación original de `phoA` de dinoflagelados es con
producción primaria; no es una regresión original frente a N:P. Las funciones
del cuaderno aparecen en las líneas 2931–2945, 3003–3015, 3043–3070,
3136–3148, 3206–3249 y 3321–3333 del commit histórico enlazado.

El código histórico aplica `p * sqrt(6/100)`, denominado «Good» en sus
comentarios. El factor es 0,244949: reduce el valor convencional aproximadamente
cuatro veces. No es una corrección por comparaciones múltiples ni resuelve la
dependencia experimental. Se conserva en una columna separada para explicar
los resultados históricos; no se usa para clasificar candidatos.

Además, casi todos los R² y p dibujados en la Figura 6 son texto fijo. En los
paneles no lineales se mezclan el ajuste NLS de la curva en escala original y
el modelo lineal de los logaritmos usado para p. Un R² de una de esas escalas
no debe atribuirse a la otra.

Los seis estados cubren PO4=0,057819–1,557318 µmol/L,
N:P=1,904998–44,020919, DIN=0,336318–3,076897 µmol/L y
producción primaria=4,046045–10,685464 µg C/L/h. Son rangos de estados
observados, no rangos de validez de un bioindicador. La distribución de PO4
está muy concentrada en valores bajos, con un estado R+P de valor alto;
las pendientes deben examinarse con sensibilidad a estados influyentes.

## Tabla recalculada del snapshot histórico

`analysis/audit_historical_regressions.R` recalcula las doce series presentes
en `figure_6_correlations.csv`. El R² de la tabla siguiente corresponde
exclusivamente al modelo lineal indicado arriba: escala logarítmica de la
respuesta en A/B/C/E y escala original en D/F. La serie Bacteroidota/`ppA`
está en los datos aunque no se dibuja en el panel final.

| Serie | Gen/razón | R² del LM | p convencional* | p con factor histórico | q BH de las 12 series* |
| --- | --- | ---: | ---: | ---: | ---: |
| Rhodobacterales / PO4 | `pst:pit` | 0,892907 | 0,004464 | 0,001094 | 0,053572 |
| Pseudomonadales / PO4 | `pst:pit` | 0,749662 | 0,025795 | 0,006318 | 0,061564 |
| Diatomeas / PO4 | `nptA` | 0,702263 | 0,037235 | 0,009121 | 0,061564 |
| Diatomeas / DIN | `ppA` | 0,718498 | 0,033051 | 0,008096 | 0,061564 |
| Bacteroidota / DIN | `ppA` | 0,758570 | 0,023903 | 0,005855 | 0,061564 |
| Diatomeas / N:P | `phoA` | 0,559808 | 0,087120 | 0,021340 | 0,104545 |
| Pseudomonadales / N:P | `phoA` | 0,483437 | 0,125121 | 0,030648 | 0,125121 |
| Dinoflagelados / PP | `acp5` | 0,804654 | 0,015360 | 0,003762 | 0,061439 |
| Dinoflagelados / PP | `pho_3` | 0,827392 | 0,011886 | 0,002911 | 0,061439 |
| Dinoflagelados / PP | `phoA` | 0,485229 | 0,124126 | 0,030404 | 0,125121 |
| Dinoflagelados / PP | `pepM` | 0,688370 | 0,041043 | 0,010053 | 0,061564 |
| Rhodobacterales / PP | `phnC` | 0,619730 | 0,063091 | 0,015454 | 0,084121 |

\* Cada serie tiene seis medias. Los valores p y q son diagnósticos bajo el
supuesto de independencia, que este diseño no satisface para inferencia de
tratamiento. BH se aplica a las doce series almacenadas, no a la familia
desconocida de todas las asociaciones exploradas originalmente. Ninguna de las
doce q es menor de 0,05. Esto no sustituye una evaluación con unidades
experimentales independientes ni demuestra ausencia de respuesta biológica.

Discrepancias comprobadas con las anotaciones fijas:

- En A, las etiquetas R²=0,893/0,750 coinciden con los modelos log-log.
  Los R² en escala original de los NLS de Rhodobacterales/Pseudomonadales
  son 0,699962/0,709396.
- En C, diatomeas/`ppA`: el R² log-lineal es 0,718498 y el del NLS en
  escala original es 0,826670; la etiqueta fija dice 0,720.
- En E, dinoflagelados/`phoA`: la etiqueta R²=0,775 y p=0,005 no se
  reproduce con el snapshot y las funciones conservadas. Se obtiene
  R² log-lineal=0,485229, R² NLS=0,968845, p convencional=0,124126 y
  p con factor histórico=0,030404.

Los NLS anteriores se comprobaron en memoria durante la auditoría, con las
fórmulas e inicializaciones del cuaderno. La exportación automática histórica
contiene las métricas LM y su escala; no presenta esos R² de NLS como si
procedieran del mismo ajuste.

## Resolución taxonómica y experimental

Las asignaciones bacterianas permiten distinguir numerosos géneros y también
familias, órdenes y nombres provisionales GTDB. Una familia o un código
provisional no equivale a una especie identificada. Se conservan los rangos
y nombres originales para revisar cada asignación; no se proponen especies.

Los encabezados de rango eucariótico no se pueden interpretar literalmente:
la columna `Species` puede contener un género como `Leptocylindrus`, mientras
que otras columnas contienen clados en niveles desplazados. El extractor
acepta únicamente las correspondencias revisadas de este conjunto. Para
`phoA`, `Leptocylindrus` es una asignación a género; los dinoflagelados
quedan como Dinophyceae o Suessiales sin género resuelto. No se convierten
automáticamente valores de `Genus`/`Species` en nombres de organismos.

El artículo describe una bolsa por tratamiento y bibliotecas técnicas:
18 procarióticas y 17 eucarióticas, con dos bibliotecas eucarióticas R+P en
0 h. No se inventa una tercera biblioteca. Algunas leyendas emplean el término
«biological replicates», en tensión con la descripción metodológica; para
esta ampliación se adopta la unidad experimental descrita en métodos. La
bolsa es la unidad del tratamiento y sus tiempos son medidas dependientes.
No puede estimarse un efecto generalizable de tratamiento a partir de una
bolsa por tratamiento. Véanse los métodos y la discusión del
[artículo](https://academic.oup.com/ismecommun/article/6/1/ycag035/8490866).

## Interpretación de las métricas nuevas

- Los TPM están normalizados por compartimento. No se comparan directamente
  magnitudes procarióticas con eucarióticas ni se reescalan únicamente los genes
  seleccionados a un millón.
- La contribución de un taxón usa como denominador el total del gen de todos
  los taxones aceptados dentro del mismo compartimento y biblioteca. En modo
  estricto se excluyen anotaciones con varios KO; el modo inclusivo muestra la
  sensibilidad a esas anotaciones.
- La fracción dentro del taxón divide por todo el RNA asignado a ese taxón,
  incluidos genes ajenos al metabolismo de P. Es una composición del RNA,
  no una medición de transcripción por célula. Faltan abundancia poblacional
  específica o DNA pareado para separar población y regulación.
- Las fuentes son dispersas y contienen TPM positivos. Una ausencia de fila
  se trata como no detección operativa al completar las bibliotecas observadas,
  no como prueba de ausencia biológica. Primero se incluyen esos ceros en las
  medias técnicas. Una fracción con denominador cero queda como NA.
- Las regresiones nuevas son descriptivas sobre seis estados. Se informa la
  sensibilidad al retirar un estado y BH para la familia explícita de taxón,
  gen, predictor y métrica de cada modo; no se utiliza p como puerta de
  clasificación. Retirar un estado no elimina toda la dependencia entre tiempos.
- Las diferencias R frente a R+P y sus rangos técnicos sirven para priorizar
  validaciones. No son ensayos replicados de efecto de fosfato. La prioridad
  por coherencia de respuesta es una regla exploratoria, no una puntuación
  validada de bioindicador.

`pst`, `pit` y la ortología `phoA` apoyan hipótesis funcionales, pero su presencia
o expresión no demuestra estrés por P en cualquier taxón. Las asignaciones
KO, incluida la etiqueta histórica K14640 para `pit`, requieren revisar
secuencias, dominios y función antes de diseñar ensayos específicos. Validar
candidatos requiere manipulación de P con bolsas o cultivos independientes,
otras causas de estrés y medidas funcionales complementarias.

La inspección de referencias de algunas anotaciones eucarióticas `phoA`
encontró correspondencias con referencias cyanobacterianas o metazoarias.
Un alineamiento con e-value muy bajo y puntuación alta sustenta homología,
pero no resuelve por sí solo identidad funcional, organismo de origen,
localización subcelular ni actividad fosfatasa. Deben revisarse dominios,
filogenia y contexto del transcrito, y descartarse contaminación o asignación
taxonómica incorrecta, antes de considerar validado un marcador eucariótico.
La taxonomía del transcrito y el linaje de la referencia de anotación se
mantienen como evidencias distintas.

## Reproducción y comprobaciones

Desde la raíz del repositorio:

```bash
Rscript --vanilla tests/test_phosphorus_indicators.R
Rscript --vanilla analysis/audit_historical_regressions.R
Rscript --vanilla analysis/run_phosphorus_indicators.R
```

La auditoría exporta `results/phosphorus_indicators/historical_regression_audit.csv`.
El módulo exporta tablas y figuras PNG/PDF/SVG en el mismo directorio. Los
rangos de las comparaciones son mínimo–máximo de bibliotecas técnicas; no
son intervalos de confianza biológicos. Los scripts no escriben en las
carpetas históricas.

Las pruebas verifican anotación exacta, rangos eucarióticos, duplicados,
denominadores, ceros operativos, número real de bibliotecas y regresiones
conocidas. Los únicos datos sintéticos se encuentran señalados dentro de las
pruebas. CI repite las pruebas y el análisis del snapshot, comprueba sus
checksums y exige que los archivos históricos permanezcan idénticos al commit
auditado. Se regeneraron las ocho figuras originales en una carpeta temporal
con R exit 0; hubo avisos de puntos DGE fuera de los límites de escala del
script histórico. No se reemplazaron los archivos comprometidos.
La regeneración de las figuras originales y la comparación visual
son comprobaciones distintas de esos controles de integridad; no se atribuye
reproducción estadística completa a una comprobación de hashes.

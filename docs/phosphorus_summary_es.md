# Resultado científico para la propuesta de limitación por fósforo

La ampliación permite priorizar tres asignaciones para validación:
**Amylibacter**, el **linaje GTDB provisional HTCC2207** y **Leptocylindrus**.
No permite identificar un dinoflagelado concreto ni proponer especies. Ninguno
es todavía un bioindicador validado de limitación por P.

| Grupo | Resultado | Evidencia descriptiva y límite principal |
| --- | --- | --- |
| Rhodobacterales | **Amylibacter**; alternativa Planktomarina | Cuatro subunidades `pst` y `pit` detectados. Cinco IDs de `pstS` aportan señal en 18/18 bibliotecas. R/R+P de TPM: 118,7 a 0 h y 37,8 a 72 h; fracción dentro del taxón: 130,2 y 54,5. No identifica especie, aislado ni regulación por célula. |
| Pseudomonadales | **HTCC2207**, género provisional de Porticoccaceae | Cuatro subunidades `pst` con patrón coherente según la regla exploratoria. `pstS` detectado en 18/18 bibliotecas: contribución media de 25,8% al TPM procariótico estricto del gen. R/R+P de TPM: 26,6 y 2,09; fracción dentro del taxón: 30,7 y 2,60. `pit` cambia en sentido distinto. Es un linaje de referencia, no una especie identificada. |
| Diatomeas | **Leptocylindrus**, género | Un ID K01077, `k141_2734405.p1`, detectado en 16/17 bibliotecas; contribución media de 44,9% al TPM eucariótico estricto de K01077. R/R+P de TPM: 18,8 y 91,4; fracción dentro del taxón: 18,9 y 82,2. La pendiente con N:P cambia de signo al retirar un estado. La referencia EggNOG es cyanobacteriana: exige verificar origen y función antes de llamarlo `phoA` algal validado. |
| Dinoflagelados | **No hay organismo concreto defendible** | En modo estricto quedan tres IDs K01077 asignados a Dinophyceae sin género; detección agregada 12/17 bibliotecas y R no supera R+P en ambos tiempos. El modo inclusivo añade anotaciones ambiguas y Suessiales, también sin género. |

Las contribuciones son medias de proporciones por biblioteca dentro de cada
compartimento y gen. Las razones describen dos bolsas concretas: no estiman
un efecto replicado del fosfato ni garantizan especificidad. Los IDs pueden
incluir redundancia entre ensamblajes; no equivalen a copias génicas únicas.

## Solidez de la evidencia

Los modelos de fracción dentro del taxón frente a N:P dan R²=0,536 para
Amylibacter/`pstS`, 0,451 para HTCC2207/`pstS` y 0,406 para
Leptocylindrus/K01077. Los p convencionales son 0,098, 0,144 y 0,173,
respectivamente, bajo un supuesto de independencia que el diseño no satisface.
Ninguna de las 2.424 asociaciones estrictas tiene q BH <0,05; mínimo 0,713.
Son cifras de transparencia, no pruebas confirmatorias.

Hay una bolsa por tratamiento. Técnicos y tiempos no aportan réplicas
biológicas independientes. La pendiente negativa con PO4 conserva su signo
al retirar estados en los tres candidatos; esto no demuestra estabilidad de
magnitud ni elimina confusión por tratamiento/tiempo. Cambios de RNA pueden
reflejar población o inversión transcripcional. No hay control poblacional
específico emparejado en las entradas usadas.

La auditoría histórica detectó que el código multiplicaba p por √(6/100) y
que las etiquetas de Figura 6 son mayoritariamente fijas. Para `phoA` de
diatomeas frente a N:P se recupera p convencional=0,087, frente al 0,021
obtenido tras ese factor. La etiqueta de `phoA` dinoflagelado frente a
producción primaria no se reproduce con el snapshot y modelos conservados.
Se documenta sin cambiar código, datos o figuras históricos.

## Qué falta

1. Recuperar secuencias y verificar dominios, función, filogenia y localización;
   separar taxonomía del transcrito y referencia funcional, especialmente para
   K01077 eucariótico. Un e-value bajo no valida por sí solo función u origen.
2. Resolver el linaje dinoflagelado con rangos y evidencia por ID correctos.
3. Medir DNA, copias génicas o células del taxón en muestras emparejadas.
4. Validar respuesta y alivio por P en bolsas/cultivos independientes, con
   controles de N, crecimiento y otros estresores y actividad enzimática.

La cautela funcional se apoya en trabajos sobre
[PhoA atípica de dinoflagelados](https://pubmed.ncbi.nlm.nih.gov/26379645/),
[regulación y actividad en cinco dinoflagelados](https://doi.org/10.3390/microorganisms7080232)
y [fosfatasa de diatomeas con P suficiente](https://pubmed.ncbi.nlm.nih.gov/33608874/).
Para controlar abundancia, véase la
[evaluación metodológica de DNA/RNA pareado](https://pmc.ncbi.nlm.nih.gov/articles/PMC8275336/).

## Texto breve reutilizable para responder a la directora

He desagregado los datos de INTERES/WP2 conservando trazabilidad por transcrito.
Los resultados priorizan Amylibacter y el linaje provisional HTCC2207 para
transporte de fosfato, y una asignación al género Leptocylindrus con ortología
de fosfatasa alcalina. Los patrones persisten parcialmente al expresar la
señal respecto al RNA total del taxón, pero siguen siendo exploratorios.
Para dinoflagelados no se puede identificar un género concreto. No presentaría
ninguno como bioindicador validado: faltan validación funcional/taxonómica,
control poblacional y replicación independiente. Además, la auditoría detectó
un factor que reducía los p históricos y etiquetas no siempre recalculadas;
lo he documentado y conservado intactos los resultados originales.

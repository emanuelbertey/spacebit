# Proyecto: auditor (Au) + plus + app Rust

## Estructura

```
auditor/    fuentes C del Au con el parche imagen-cero y sus tests
plus/       fuentes C del Plus (BG_0 = 37162531*N/50) y tests del serial
create/     app Rust que embebe al Au por FFI (con sus propios native/)
create_plus/ app Rust que embebe al Plus por FFI (shim + tests)
security/   auditoria de parametros + estimacion MSIS (scripts y README)
reportes/   informes
test/       build.sh + RNG portable + tests de falsificacion
.github/    workflows: build (C/Rust) + seguridad
```

Dos variantes nomas: **auditor** y **plus**. Un `.exe` por variante.

## Seguridad (resumen)

Los esquemas vienen de **MatRiCT** (IEEE S&P 2022), basados en **retículos
(module-SIS / module-LWE)**, post-cuánticos. El secreto no es de 256 bits:
es un vector de coeficientes ternarios `{-1,0,1}` (~3855 bits en Au, ~3246
en Plus); los 256 bits son el seed AES-256 y el serial/hash.

| | Au | Plus |
|---|---|---|
| `d` | 64 | 256 |
| `q` | 2147221513 | 3489661057 |
| clave | 38×64 coef. | 8×256 coef. |
| tope escalado | `TG × N/100` | `BG_0 = 37162531 × N/50` |

Parches: serial cero rechazado en `verify.c`, texto firmado dentro del hash,
topes que escalan solos, y RNG portable para Windows.

Detalle completo en `security/README.md`.

## Workflows

- `build-binarios.yml`: compila `auditor.exe`, `plus.exe` y los tests de
  1000 llaves / falsificación (MSYS2/MinGW). No ejecuta tests.
- `build-rust.yml`: `cargo build --release` en `create/` y `create_plus/`.
- `test-rust.yml`: manual; compila los `.exe` de test de ambos crates.
- `security.yml`: auditoría de parámetros + estimación MSIS con los scripts
  de `main/`. Publica resultados como artifact y en el resumen del job.
  Cada job tiene `timeout-minutes: 60`.

## Tests

- `auditor/test_todos.c`: 4 pruebas del Au en un exe (honesta/auditoría,
  firma falsa, falsificador, observador).
- `auditor/test_falsificador.c`, `auditor/test_firma_falsa.c`,
  `auditor/test_auditor.c`: fuentes de esas pruebas.
- `test/falsificacion_rango_alto.c`: valores ×10 a ×1000 del recomendado.
- `test/falsificacion_agresiva.c`: bordes de norma, serial ajeno, anillo
  sustituido y prueba armada desde cero.
- `plus/`: `test.c`, `tp1.c`, `keygen_zero.c`.
- `create/tests/v4.rs` (6 tests) y `create_plus/tests/plus.rs` (3 tests).

## RAG

Indice consultable de todo el repo (C/Rust/Python/docs), sin APIs ni modelos
descargados: búsqueda por tokens con puntaje, más búsqueda exacta de
`#define` y firmas de funciones. El índice se guarda en `.indice.json`
(ignorado por git) y se reconstruye con `--rebuild`.

- `python test-d128/rag.py` construye/actualiza el índice.
- `python test-d128/rag.py "texto a buscar"` consulta por similitud (opción `-k N`).
- `python test-d128/rag.py --def NAME` / `--fn NAME` búsqueda exacta.

Ver `test/README.md` y `security/README.md`.

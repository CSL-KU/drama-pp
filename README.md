# DRAMA++: A fast and robust DRAM address-map reverse-engineering tool

This code is based on [DRAMA](https://github.com/isec-tugraz/drama) and includes the following improvements:

- Added support for the ARM64 architecture.
- Implemented a faster GF(2) solver with polynomial-time complexity (the original version has exponential-time complexity).
- Fixed a logical bug that prevented high-order physical address bits from being considered.
- Fixed a logical bug that left the `base` address in the address pool when it should have been added to the set.
- Improved timing measurements.
- Additional changes for improved functionality, usability, reliability, and portability.

To see all changes, run:

```
git diff c5c83471...HEAD re/measure.cpp
```

## Usage

### Prerequisites

- Linux (x86-64 or ARM64) with a recent kernel.
- `g++` and `make` installed.
- Permission to read `/proc/self/pagemap` (often requires `sudo`).
- The tool attempts huge pages first and falls back to regular pages if unavailable.

### Build

```
cd re
make
```

This produces the `measure` binary in `re/`.

The region-aware offline solver lives in `re/gf2_bank_solver.py`.

### Run

Measure DRAM bank functions and save them to `map.txt`:

```
./measure [-b <start bit>] [-e <end bit>] [-c <cpu core>] [-r <scale factor>] [-m <memory size>] [-g <memory size in GB>] [-i <outer loops>] [-j <inner loops>] [-s <expected sets>] [-q <sets for early quit>] [-t <threshold cycles>] [-v <verbosity>] [-f <output file>]

```

Notes:
- `-s`: expected sets = DIMMs × ranks × banks (e.g., 1×2×8 = 16).
- `-m`/`-g`: memory to map. `-m` accepts MB by default and also supports `M`/`G` suffixes (e.g., `-m 1024`, `-m 1G`).
- `-c`: pin to a CPU core (you can also use `taskset`).
- `-r`: timing scale factor (advanced tuning).
- `-i`/`-j`: outer/inner loop counts; ARM64 may benefit from a higher `-j`.
- `-t`: timing threshold (cycles) to override auto gap detection.
- `-b`/`-e`: search bit window (defaults: 5..40).
- `-q`: stop after N sets are found.
- `-v`: verbosity level.
- `-f`: output file for discovered functions (default `map.txt`).

### Outputs

- `setN.txt`: physical addresses of each discovered same-bank set.
- `map.txt`: one line per XOR function with the physical address bit indices.
- `recovered_bank_mapping.txt`: output of `gf2_bank_solver.py`.
    On symmetric systems this is the legacy flat format.
    On asymmetric systems it can be a piecewise format with region blocks and selector terms.

The repository also includes a validated asymmetric example mapping and matching manual region file:

- `re/recovered_bank_mapping-hynix-2x4GB+elpida-1x8GB-correct.txt`
- `re/regions-hynix-2x4GB+elpida-1x8GB-correct.txt`

### Offline Solver

The Python solver can still recover a single flat mapping from `setN.txt` files:

```
python3 gf2_bank_solver.py --files set*.txt
```

For asymmetric systems, it also supports region-aware recovery:

```
python3 gf2_bank_solver.py --files set*.txt --auto-regions
```

Useful tuning flags for asymmetric layouts:

- `--auto-top-splits <N>`: limit how many promising selector terms are explored per unresolved branch.
- `--auto-max-selector-bits <N>`: cap the XOR width of one selector term.
- `--auto-max-depth <N>`: cap recursive region splitting depth.
- `--auto-diagnostics-limit <N>`: bound how many unresolved branches and suggested follow-up splits are reported on failure.
- `--auto-trace-depth <N>`: keep `--verbose` auto-region tracing readable by limiting how deep recursive branch logs go.

To save the recovered selector blocks as a starter manual regions file:

```
python3 gf2_bank_solver.py --files set*.txt --auto-regions --emit-regions-file regions.txt
```

For example:

```
python3 gf2_bank_solver.py --files set*.txt --auto-regions --auto-top-splits 12 --auto-diagnostics-limit 6 --auto-trace-depth 1
```

If the heuristic search is not sufficient, provide explicit region definitions:

```
python3 gf2_bank_solver.py --files set*.txt --regions-file regions.txt
```

Manual region files support two syntaxes:

```
region 0x40 0x0
region 0x40 0x40
```

or the more general selector-block form:

```
region
selector 0 6
selector 1 14 17

region
selector 1 6
selector 0 16 19
```

Each `selector` line means that the XOR of the listed physical-address bits must equal the given `0` or `1` value.

### Example

1 DIMM, 1 channel, 2 ranks, 8 banks (16 sets), mapping 1 GB:

```
sudo ./measure -s 16 -g 1
```

## DRAM Bank Map Database
See: [Found-DRAM-BankMap.md](./Found-DRAM-BankMap.md) for examples of discovered DRAM bank-mapping functions.

## Speed Comparison

| Platform | DRAMA | DRAMA++ |
|-----------------------|------------|----------|
| Xeon E3-1220 v5 (64 banks) | 54.5s<sup>1</sup> | 3.4s<sup>2</sup> |
| Raspberry Pi 4 (8 banks) | N/A | 0.6s |

- <sup>1</sup> DRAMA option used: `-s 64 -n 10` (the default `n=5000` took more than 10 minutes and did not recover the map).
- <sup>2</sup> DRAMA++ option used: `-s 64` (manually setting a threshold such as `-t 300` can make it even faster and more reliable).

## Limitations
DRAMA++ no longer silently accepts a bogus single global mapping on asymmetric configurations. The offline solver can now:

- fail loudly when the discovered sets cannot be labeled by one unique global XOR mapping,
- recover piecewise local mappings from explicit region definitions, and
- attempt recursive heuristic region discovery using raw-bit and XOR-based selector terms.

The automatic region search is best-effort. Some highly interleaved asymmetric layouts may still require a manual `--regions-file` to fully separate the local mappings.

When auto-region search cannot finish, the solver now reports the unresolved selector branch, the affected `setN.txt` files, and a short list of promising next selector terms to try manually.

## Citation

If you use this tool, please cite:

    @inproceedings{sullivan2026rtas,
        title = {{Per-Bank Memory Bandwidth Regulation for Predictable and Performant Real-Time Systems}},
        author = {Connor Sullivan and Amin Mamandipoor and Cole Strickler and Heechul Yun},
        booktitle = {IEEE Real-Time and Embedded Technology and Applications Symposium (RTAS)},
        year = {2026},
        month = {May}
    }

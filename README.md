# Custom XSPEC Cross-ARF Packed Model

This directory contains the custom XSPEC `xcpkg` package for packed cross-ARF fitting.
It is intended for the packed `xcrsarfrun` workflow only.

The package provides:

- `xcrsarf`: a local XSPEC mixing model derived from XSPEC's `crossarf` model.
- `xcrsarfrun`: an XSPEC command that reads one config file, loads spectra and responses, creates the packed model, and ties parameters/constants.
- packed source mixing through a single `xpack` model, with optional source-only (no spectra) point-source/AGN models added as normal XSPEC source slots.

## Files

- `XCrossArf.cxx`, `XCrossArf.h`: packed cross-ARF mixing model.
- `XCrossArfSetup.cxx`: `xcrsarfrun` setup command and config parser.
- `xcrsarfWrapper.cxx`, `xcrsarfWrapper.h`: XSPEC wrapper entry point.
- `lmodel_xcrsarf.dat`: XSPEC local-model definition used by the build.
- `xcrsarf_packed_config.example`: packed-mode config example.
- `build_xcrsarf.sh`: builds the `xcpkg` package into `build/xcpkg`.

## Build

From this directory:

```bash
cd /path/to/Custom_Xspec_CrossARF_Model
./build_xcrsarf.sh
```

The build output is stable and local to this folder:

```text
/path/to/Custom_Xspec_CrossARF_Model/build/xcpkg/libxcpkg.dylib
```

## Load In XSPEC

Start XSPEC from the same HEASoft environment, then load the package:

```text
XSPEC12> lmod xcpkg /path/to/Custom_Xspec_CrossARF_Model/build/xcpkg
```

## Run A Packed Config

Use a packed config file:

```text
XSPEC12> xcrsarfrun /path/to/config.txt
```

The config should begin with:

```text
mode packed
sources N
```

In packed mode, observed region sources are evaluated inside one model named `xpack`.
For large A/B runs this avoids loading the full XSPEC source-slot response matrix.

## Model Blocks

The source model blocks may use the convenient XSPEC-like source syntax:

```text
model 1:s1 apec*const
10.1 0.1 1.0 1.0 25.0 25.0
0.4 -0.01 0.1 0.1 1.0 1.0
0.0546 -0.01 0.01 0.01 0.30 0.30
0.002 1e-6 1e-6 1e-6 1e-1 1e-1
1.0 -1 0.0 0.0 2.0 2.0
```

`xcrsarfrun` strips the `model 1:s1` convenience wrapper and loads the source into `xpack`.
Source 1 is the default template for any source without an explicit model block.

If the config mixes `apec*const` and `po*const` sources, `xcrsarfrun` expands them to the shared packed expression `(apec+po)*const` and freezes/zeros the inactive component rows.
The constant must remain the last parameter.

## Spectrum Lines

Observed sources need a `spectrum` row:

```text
spectrum SPEC REGION OBS_LABEL PHA RMF DIAGONAL_ARF [MIX_GROUP]
```

Example:

```text
spectrum 1 1 A 70201001002/arfRun18RegA1.pha 70201001002/arfRun18RegA1.rmf 70201001002/arfRun18RegA1_1.arf
spectrum 15 1 B 70201001002/arfRun18RegB1.pha 70201001002/arfRun18RegB1.rmf 70201001002/arfRun18RegB1_1.arf
```

`SPEC` should be the XSPEC data group number. Keeping these contiguous is simplest.
`REGION` is the physical source/region number.
`OBS_LABEL` is the constant label, such as `A` or `B`.

For multiple observations that reuse A/B labels but should mix only within their own observation, add `MIX_GROUP`:

```text
spectrum 1 1 A obs001/regA1.pha obs001/regA1.rmf obs001/regA1_1.arf obs001_A
spectrum 2 1 B obs001/regB1.pha obs001/regB1.rmf obs001/regB1_1.arf obs001_B
spectrum 3 1 A obs002/regA1.pha obs002/regA1.rmf obs002/regA1_1.arf obs002_A
spectrum 4 1 B obs002/regB1.pha obs002/regB1.rmf obs002/regB1_1.arf obs002_B
```

## ARF Lines

Every source-to-target contribution is declared with:

```text
arf SOURCE_REGION TARGET_SPEC ARF_PATH
```

Example:

```text
arf 1 1 obsA/reg01_1.arf
arf 2 1 obsA/reg02_1.arf
arf 1 15 obsB/reg01_1.arf
arf 2 15 obsB/reg02_1.arf
```

`TARGET_SPEC` is the `SPEC` number from a `spectrum` row.

## Observation Constants

Set the reference observation label with:

```text
reference_obs A
```

The last parameter of each source model is treated as the observation constant.
The reference label is fixed to 1. Other labels are tied across all spectra and sources.
For example, all B constants tie together, including source-only point-source models.

## Source-Only Point Sources/AGN

A point source can have a model and cross-ARFs without an observed spectrum.
To do that:

1. Include it in `sources N`.
2. Add a normal `model SOURCE:sSOURCE ...` block.
3. Do not add a `spectrum` row for that source number.
4. Add `arf SOURCE TARGET_SPEC ARF_PATH` rows for every target spectrum it contributes to.

Example:

```text
sources 16

model 15:s15 po*const
1.176 -0.1 0.1 0.1 3.0 3.0
6.1e-5 1e-6 1e-6 1e-6 1e-1 1e-1
1.0 -1 0.0 0.0 2.0 2.0

arf 15 1 70201001002/arfRun18RegA15_1.arf
arf 15 15 70201001002/arfRun18RegB15_1.arf
```

In this case `xcrsarfrun` creates `xpack` for the observed region sources, then creates a normal XSPEC source model such as `s15` in source slot 15.
The point-source constants are tied directly to the corresponding `xpack` constants.

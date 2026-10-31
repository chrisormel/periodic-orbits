!NOTE: This readme file has been automatically generated!

# Periodic-orbit cubes and fitted planetary systems

This repository provides periodic-orbit (fPO) solution cubes and fitted orbital
solutions for the PDS 70, GJ 876, and HR 8799 planetary systems, together with
Python examples for loading the solutions and integrating them with REBOUND.

You are welcome to use these data and utilities. Please cite the accompanying
preprint/paper by **Ormel, Yi, and Huang** when using them in your work.
Full bibliographic details and a preprint link will be added here.

## Requirements

Use Python with NumPy, pandas, Astropy, and REBOUND:

```bash
python -m pip install numpy pandas astropy rebound
```

The current helpers use the REBOUND 5.2 interface. Matplotlib is also needed
to run the plotting example in `POfit/Solutions/runsolutions.py`:

```bash
python -m pip install matplotlib
```

## PO: loading the fPO solution cubes

`PO/Cubes/` contains the solution tables as pickled pandas DataFrames.
`PO/Helper/` contains the supporting functions and the example `runpo.py`.

Run the example from its directory:

```bash
cd PO/Helper
python runpo.py
```

At the bottom of `runpo.py`, select the table (`ftable`), cube-grid identifier
(`cid`), and solution key (`key`). A grid node can contain several solutions,
identified by keys such as `A0`; some nodes contain no solutions.
Set `npo` and `npt` to choose the integration duration and number of output times.

The example loads the grid and solution vectors, prints the node properties,
constructs a REBOUND simulation, and prints a sample of orbital quantities.
The two functions can also be used from Python while working in `PO/Helper/`:

```python
import numpy as np
from runpo import load_solution, run_solution

xv, xs, tpo, chain = load_solution('../Cubes/321low.pkl', 25000, 'A0')
tarr = np.linspace(0, tpo, 21)
run_solution(chain, xv, xs, tarr, tpo)
```

Here `xv` contains the grid parameters (period-ratio offset and planet masses),
`xs` is the PO solution vector, and `tpo` is the PO period. Planets are ordered
from the innermost outward. The simulation uses Jacobi orbital elements,
`G = 1`, central mass `= 1`, and an initial inner-planet period of `2*pi`.
Masses are relative to the central mass, and printed angles are in radians.

To construct a simulation for your own analysis:

```python
from runpo import POsol

posim = POsol(chain, *xv)
sim = posim.setup_sim(xs)
sim.integrate(tpo)
```

## POfit: fitted solutions for observed systems

`POfit/Solutions/` contains our best-fit orbital solutions to the astrometry
and radial-velocity (RV) data of PDS 70, GJ 876, and HR 8799, and the independent
REBOUND loader `runsolutions.py`. It does not require the original fitting code.

Files are named by system, solution number, and optimization flag, for example:

```text
gj876-1.pkl
gj876-1-optim.pkl
pds70-1.pkl
hr8799-1.pkl
```

The `-optim` suffix selects the locally optimized version of a solution.
Check the directory for the available system names and solution numbers.

From `POfit/Solutions/`, load a fitted system with:

```python
import runsolutions as runsol

sim, metadata = runsol.run_named_solution('gj876', 1, optim=True)
tpo = metadata['tpo']
sim.integrate(sim.t + tpo)
```

The return values are the REBOUND simulation and a metadata dictionary
containing the PO period (`tpo`) and output units (`units`). By default,
masses are in stellar-mass units, distances are in AU for RV solutions or mas
for astrometric solutions, and times are in days.

Each saved file contains the planet-input dictionaries (`planL`), shared
parameters and source units (`dother`), the initial epoch (`tstart_mjd`), and
the dimensionless PO period (`tpo`). The loader converts that stored period
to the simulation's time unit in `metadata['tpo']`.

### Choosing output units

Use `outunits='mass-length-time'`, for example:

```python
sim, metadata = runsol.run_named_solution(
    'gj876', 1, optim=True, outunits='mstar-au-year'
)
```

`mstar` keeps masses normalized to the central star. For physical masses,
use `msun`, or the shortcuts `SI` (`kg-m-s`) and `cgs` (`g-cm-s`).

Converting an angular solution to physical distances requires a stellar mass
or parallax, unless already supplied in the saved metadata:

```python
sim, metadata = runsol.run_named_solution(
    'hr8799', 1, optim=True,
    outunits='msun-au-year', dadd={'mstar': 1.5}
)
```

This stellar mass is an illustrative choice. Alternatively, supply
`dadd={'parallax': ...}` with parallax in mas. A manual `mstar` is in the
requested physical mass unit; for normalized `mstar` output it is in solar
masses. Saved stellar masses are always in solar masses.

The fitted `gN` remains authoritative. For AU input it determines the stellar
mass. For mas input, a stellar mass determines parallax, or vice versa.
When both are supplied for mas input, stellar mass takes priority and the
loader warns that the redundant parallax is ignored.

`sim.t` labels the epoch of the initial PO configuration. With default output
it is in MJD days; with year output it is in Julian years since the MJD origin,
rather than calendar decimal years. It generally precedes the earliest
observation. Use `sim.t + duration` as an integration target when starting
from the loaded state.

For a complete integration and plotting example:

```bash
cd POfit/Solutions
python runsolutions.py
```

The example writes `plot.png` in the working directory.

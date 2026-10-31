"""Load an exported orbital solution and initialize REBOUND."""

from pathlib import Path
import pickle
import warnings

import rebound
import numpy as np
from astropy import units as u
from astropy.constants import G


SOLUTION_DIR = Path(__file__).resolve().parent


def solution_path(system, number, optim=False):
    """Return the conventional path for an exported solution."""
    suffix = "-optim" if optim else ""
    return SOLUTION_DIR / f"{system}-{number}{suffix}.pkl"


def load_solution(path):
    """Load and validate saved REBOUND inputs and the PO period."""
    path = Path(path)
    if not path.is_absolute():
        path = SOLUTION_DIR / path

    with path.open("rb") as stream:
        solution = pickle.load(stream)

    required = {"planL", "dother", "tstart_mjd", "tpo"}
    missing = required.difference(solution)
    if missing:
        names = ", ".join(sorted(missing))
        raise ValueError(f"{path} is missing required fields: {names}")
    return solution

def parse_outunits(outunits):
    """Return REBOUND names and Astropy units for a physical unit system."""
    shortcut = {"cgs": "g-cm-s", "si": "kg-m-s"}
    spec = shortcut.get(outunits.lower(), outunits)
    fields = spec.split("-")
    if len(fields) != 3:
        raise ValueError(
            "outunits must be 'default', 'cgs', 'SI', or mass-length-time"
        )

    mass_name, length_name, time_name = (field.lower() for field in fields)
    aliases = {
        "msun": ("msun", u.Msun),
        "msolar": ("msun", u.Msun),
        "solmass": ("msun", u.Msun),
        "year": ("yr", u.yr),
        "years": ("yr", u.yr),
    }

    def get_unit(name):
        if name == "mstar":
            return "mstar", None
        if name in aliases:
            return aliases[name]
        try:
            return name, u.Unit(name)
        except ValueError as exc:
            raise ValueError(f"unsupported output unit: {name}") from exc

    mass_name, mass_unit = get_unit(mass_name)
    length_name, length_unit = get_unit(length_name)
    time_name, time_unit = get_unit(time_name)
    if mass_unit is not None and not mass_unit.is_equivalent(u.kg):
        raise ValueError(f"{mass_name} is not a mass unit")
    if length_name == "mas" and mass_name == "mstar":
        pass
    elif not length_unit.is_equivalent(u.m):
        raise ValueError(f"{length_name} is not a length unit")
    if not time_unit.is_equivalent(u.s):
        raise ValueError(f"{time_name} is not a time unit")
    return (
        (mass_name, length_name, time_name),
        (mass_unit, length_unit, time_unit),
    )


def physical_scales(dother, dadd, mass_unit, length_unit):
    """Return stellar mass and one input length in requested output units."""
    inunits = dother.get("inunits")
    if not isinstance(inunits, dict):
        raise ValueError("physical conversion requires dother['inunits']")
    if inunits.get("mass") != "mstar" or inunits.get("time") != "day":
        raise ValueError(f"unsupported input unit system: {inunits}")

    source_length = inunits.get("length")
    if source_length not in ("au", "mas"):
        raise ValueError(f"unsupported input length unit: {source_length}")

    gmstar_au_day = float(dother["gN"])
    G_au_day = G.to_value(u.au**3 / (u.Msun * u.day**2))

    # Saved mstar values are always in Msun.  A dadd override follows the
    # requested output mass unit, as specified by the public API.
    if "mstar" in dadd:
        mstar = dadd["mstar"] * mass_unit
    elif "mstar" in dother:
        mstar = dother["mstar"] * u.Msun
    else:
        mstar = None

    parallax = dadd.get("parallax", dother.get("parallax"))
    if parallax is not None:
        parallax = float(parallax)  # always expressed in mas

    if source_length == "au":
        inferred_mass = (gmstar_au_day / G_au_day) * u.Msun
        if "mstar" in dadd:
            warnings.warn(
                "[runsolutions]Ignoring manual mstar because gN determines the stellar mass for an au solution.",
                UserWarning, stacklevel=2,)
        mstar = inferred_mass
        length_scale = (1.0 * u.au).to_value(length_unit)
    else:
        if mstar is None and parallax is None:
            raise ValueError(
                "a mas solution requires mstar or parallax for physical units"
            )
        if mstar is None:
            mstar = (
                gmstar_au_day / (G_au_day * parallax**3)
            ) * u.Msun
        elif parallax is None:
            mstar_msun = mstar.to_value(u.Msun)
            parallax = (gmstar_au_day / (G_au_day * mstar_msun))**(1/3)
        else:
            warnings.warn(
                "[runsolutions]:Ignoring manual parallax because mstar and the fitted gN determine the distance scale.",
                UserWarning, stacklevel=2,)
            mstar_msun = mstar.to_value(u.Msun)
            parallax = (gmstar_au_day / (G_au_day * mstar_msun))**(1/3)
        length_scale = ((1.0 / parallax) * u.au).to_value(length_unit)

    return mstar.to_value(mass_unit), length_scale


def normalized_length_scale(dother, dadd, length_name, length_unit):
    """Convert one input length while retaining stellar-mass units."""
    inunits = dother.get("inunits")
    if not isinstance(inunits, dict):
        raise ValueError("unit conversion requires dother['inunits']")
    if inunits.get("mass") != "mstar" or inunits.get("time") != "day":
        raise ValueError(f"unsupported input unit system: {inunits}")

    source_length = inunits.get("length")
    if source_length not in ("au", "mas"):
        raise ValueError(f"unsupported input length unit: {source_length}")
    if source_length == length_name:
        return 1.0

    parallax = dadd.get("parallax", dother.get("parallax"))
    mstar = dadd.get("mstar", dother.get("mstar"))

    if source_length == "mas" and mstar is not None:
        if parallax is not None:
            warnings.warn(
                "Ignoring manual parallax because mstar and the fitted gN "
                "determine the distance scale.",
                UserWarning,
                stacklevel=2,
            )
        G_au_day = G.to_value(u.au**3 / (u.Msun * u.day**2))
        parallax = (
            float(dother["gN"]) / (G_au_day * float(mstar))
        )**(1/3)
    elif source_length == "au" and mstar is not None and parallax is not None:
        origin = "manual" if "mstar" in dadd else "stored"
        warnings.warn(
            f"Ignoring {origin} mstar because parallax determines the "
            "angular scale of an au solution.",
            UserWarning,
            stacklevel=2,
        )

    if parallax is None:
        raise ValueError(
            f"converting {source_length} to {length_name} requires parallax"
        )
    parallax = float(parallax)  # always expressed in mas

    if source_length == "mas":
        one_input_length = (1.0 / parallax) * u.au
        return one_input_length.to_value(length_unit)
    if length_name == "mas":
        return parallax
    return (1.0 * u.au).to_value(length_unit)


def setup_simulation(solution, outunits="default", dadd=None):
    """Construct a REBOUND simulation and return ``(sim, metadata)``.

    outunits:
    - default unit system as described in Paper 1
        * mass in stellar mass units
        * distances in au (RV) OR mas (astrometry)
        * time in units of days

    otherwise, provide as {munit}-{dunit}-{tunit} or cgs/SI. Examples:

        msun-au-year, cgs, SI

    when providing such physical units, you may need to provide additional
    information through "dadd", e.g., when "msun" in outunits:

        dadd = {'mstar': 1.0} # --> assumes 1 solar mass star.

    If already present in the "dother" solution object, "dadd" overwrites 
    these.
    """

    planL = solution["planL"]
    dother = solution["dother"].copy()
    dadd = {} if dadd is None else dict(dadd)

    sim = rebound.Simulation()
    if outunits.lower() == "default":
        stellar_mass = 1.0
        length_scale = 1.0
        time_scale = 1.0
        sim.G = dother["gN"]
        output_units = dict(dother["inunits"])
    else:
        names, units = parse_outunits(outunits)
        mass_name, length_name, time_name = names
        mass_unit, length_unit, time_unit = units
        time_scale = (1.0 * u.day).to_value(time_unit)
        if mass_name == "mstar":
            stellar_mass = 1.0
            length_scale = normalized_length_scale(
                dother, dadd, length_name, length_unit
            )
            sim.G = dother["gN"] * length_scale**3 / time_scale**2
            output_units = {
                "mass": "mstar",
                "length": length_name,
                "time": time_name,
            }
        else:
            stellar_mass, length_scale = physical_scales(
                dother, dadd, mass_unit, length_unit
            )
            sim.units = (length_name, time_name, mass_name)
            sim.G = G.to_value(length_unit**3 / (mass_unit * time_unit**2))
            output_units = {
                "mass": mass_name,
                "length": length_name,
                "time": time_name,
            }

    sim.add(m=stellar_mass)

    for index, saved_args in enumerate(planL):
        particle_args = saved_args.copy()
        particle_args["m"] *= stellar_mass
        if "a" in particle_args:
            particle_args["a"] *= length_scale
        if "P" in particle_args:
            particle_args["P"] *= time_scale
        if "prat" in particle_args:
            if index == 0:
                raise ValueError("the innermost planet cannot use 'prat'")
            particle_args["P"] = (
                sim.particles[index].P * particle_args.pop("prat")
            )
        sim.add(
            inc=dother["inc"],
            Omega=dother["Omega"],
            **particle_args,
        )

    sim.move_to_com()
    sim.t = solution["tstart_mjd"] * time_scale

    metadata = {
        "tpo": solution["tpo"] * sim.particles[1].P / (2*np.pi),
        "units": output_units,
    }
    return sim, metadata


def run_solution(path, outunits="default", dadd=None):
    """Load *path* and return ``(simulation, metadata)``."""
    return setup_simulation(load_solution(path), outunits=outunits, dadd=dadd)


def run_named_solution(system, number, optim=False, outunits="default", 
                       dadd=None):
    """Load a named solution and return ``(simulation, metadata)``."""
    return run_solution(solution_path(system, number, optim=optim), 
                        outunits=outunits, 
                        dadd=dadd,)


if __name__=="__main__":
    #run an example
    import numpy as np

    sim, metadata = run_named_solution('pds70', 1, optim=False)
    ps = sim.particles
    #tstart = time where the fPO cycle that contained the first 
    #data point began
    tstart = sim.t
    tpo = metadata["tpo"]
    nplan = len(ps)-1
    nres = nplan-1
    jres = 1

    npo = 5
    npt = npo *200
    tarr = np.linspace(tstart, tstart+tpo*npo, npo*npt)
    ntime = len(tarr)

    #calculate the 2BR angles
    phiarr = np.zeros((ntime,nres,2))
    eccarr = np.zeros((ntime,nplan))
    delarr = np.zeros((ntime,nplan-1)) #resonance offset
    for i,time in enumerate(tarr):
        sim.integrate(time)
        for k in range(nres):
            phiarr[i,k,0] = (2*ps[k+2].l -ps[k+1].l -ps[k+1].pomega)%(2*np.pi) #inner
            phiarr[i,k,1] = (2*ps[k+2].l -ps[k+1].l -ps[k+2].pomega)%(2*np.pi) #outer
            delarr[i,k] = ps[k+2].P/ps[k+1].P *jres/(jres+1) -1 #our definition
        for k in range(nplan):
            eccarr[i,k] = ps[k+1].e

    #plot stuff
    import pylab as pl
    fg, axL = pl.subplots(1,4, figsize=(9,4), constrained_layout=True)

    for k in range(nres):
        axL[0].plot(tarr,phiarr[:,k,0], '.', ms=1) #inner resonance
        axL[1].plot(tarr,phiarr[:,k,1], '.', ms=1) #inner resonance
        axL[3].plot(tarr,delarr[:,k], lw=0.5) #resonance offset
    for k in range(nplan):
        axL[2].plot(tarr,eccarr[:,k]) #eccentricity

    for ax in axL:
        ax.set_xlabel('time [mjd]')

    for ax in axL[:2]:
        ax.set_ylim(0,2*np.pi)

    axL[0].set_ylabel('inner 2B angle')
    axL[1].set_ylabel('outer 2B angle')
    axL[2].set_ylabel('eccentricity')
    axL[3].set_ylabel('resonant offset')

    fg.savefig('plot.png')
    pl.close(fg)

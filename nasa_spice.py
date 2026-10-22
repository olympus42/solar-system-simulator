"""
nasa_spice.py
=============
Real Solar System data from NASA's SPICE Toolkit.

SPICE is the system NASA's Jet Propulsion Laboratory (JPL) uses to answer one
question for every mission it flies: *where is everything, and when?* Cassini,
Juno, New Horizons and the Mars rovers all plan with it.

SPICE keeps its knowledge in files called **kernels**. We need three:

    naif0012.tls   leap seconds   -> turns "2026-10-22" into SPICE time
    de440s.bsp     ephemeris      -> where every planet is, 1849-2150
    gm_de440.tpc   constants      -> G x mass of the Sun and each planet

They are downloaded once from NASA (about 32 MB) into the `kernels/` folder.

We talk to SPICE from Python through SpiceyPy (`pip install spiceypy`).

Units match solar_system.py: AU, years, solar masses.
"""

import os
import ssl
import sys
import urllib.request
from pathlib import Path

import numpy as np
import spiceypy as spice

from solar_system import BODIES

# ----------------------------------------------------------------------
# 1. Getting the kernels from NASA
# ----------------------------------------------------------------------

NAIF = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels"
KERNEL_URLS = {
    "naif0012.tls": f"{NAIF}/lsk/naif0012.tls",
    "de440s.bsp":   f"{NAIF}/spk/planets/de440s.bsp",
    "gm_de440.tpc": f"{NAIF}/pck/gm_de440.tpc",
}

# Kernels live next to this file, unless you point SPICE_KERNELS somewhere else.
KERNEL_DIR = Path(os.environ.get("SPICE_KERNELS", Path(__file__).parent / "kernels"))

AU_KM = 149_597_870.700          # 1 AU in km (IAU 2012 definition, exact)
YEAR_S = 365.25 * 86_400.0       # 1 Julian year in seconds


def download_kernels():
    """Download any kernel that isn't in kernels/ yet. Safe to run again."""
    KERNEL_DIR.mkdir(parents=True, exist_ok=True)
    for name, url in KERNEL_URLS.items():
        target = KERNEL_DIR / name
        if target.exists() or _alternative_exists(name):
            continue
        print(f"Downloading {name} from NASA ...", flush=True)
        try:
            _download(url, target)
        except Exception as err:          # network or certificate trouble
            sys.exit(_download_help(name, url, err))


def _alternative_exists(name):
    # The full DE440 (de440.bsp, 114 MB, years 1550-2650) works just as well.
    return name == "de440s.bsp" and (KERNEL_DIR / "de440.bsp").exists()


def _download(url, target):
    tmp = target.with_suffix(target.suffix + ".part")
    with urllib.request.urlopen(url, timeout=60) as resp, open(tmp, "wb") as out:
        total = int(resp.headers.get("Content-Length") or 0)
        done = 0
        while chunk := resp.read(1 << 20):
            out.write(chunk)
            done += len(chunk)
            if total:
                print(f"\r  {done / 1e6:5.1f} / {total / 1e6:.1f} MB", end="", flush=True)
    tmp.rename(target)
    print()


def _download_help(name, url, err):
    hint = ""
    if isinstance(err, ssl.SSLError) or "CERTIFICATE" in str(err).upper():
        hint = ("\nOn a Mac with Python from python.org, run the\n"
                "'Install Certificates.command' file in Applications/Python 3.x, then try again.")
    return (f"\nCould not download {name}: {err}{hint}\n"
            f"You can also download it in your browser:\n  {url}\n"
            f"and put it in: {KERNEL_DIR}\n")


# ----------------------------------------------------------------------
# 2. Loading them into SPICE
# ----------------------------------------------------------------------

_loaded = False


def load_kernels():
    """Make sure the kernels are here, then hand them to SPICE (once)."""
    global _loaded
    if _loaded:
        return
    download_kernels()
    # Load every kernel file in the folder: .tls, .bsp, .tpc
    for path in sorted(KERNEL_DIR.iterdir()):
        if path.suffix in (".tls", ".bsp", ".tpc"):
            spice.furnsh(str(path))
    _loaded = True


def to_et(date):
    """'2026-10-22' or '2026-10-22 18:30' (UTC) -> SPICE ephemeris time.

    Ephemeris time (ET) = seconds since noon on 1 January 2000 (TDB).
    The leap-seconds kernel is what makes this conversion exact.
    """
    load_kernels()
    return spice.str2et(str(date))


# ----------------------------------------------------------------------
# 3. Asking NASA where the planets are
# ----------------------------------------------------------------------

# SPICE numbers every object. Planets with moons are tracked by the
# *barycenter* of the planet + its moons (that's what gravity "sees" from far
# away). For Earth that means the Earth-Moon barycenter.
SPICE_ID = {
    "Sun": 10,
    "Mercury": 1, "Venus": 2, "Earth": 3, "Mars": 4,
    "Jupiter": 5, "Saturn": 6, "Uranus": 7, "Neptune": 8,
}

# Positions are measured from the Solar System Barycenter (SPICE id 0) -- the
# true centre of mass, which the Sun itself wobbles around -- in the plane of
# Earth's orbit (the ecliptic), so x-y is the familiar top-down view.
FRAME = "ECLIPJ2000"
ORIGIN = "0"


def gm_au_yr(name):
    """G x mass of a body, from NASA's constants kernel, in AU^3 / yr^2."""
    load_kernels()
    gm_km3_s2 = spice.bodvrd(str(SPICE_ID[name]), "GM", 1)[1][0]
    return gm_km3_s2 * YEAR_S**2 / AU_KM**3


def real_state(date, names=None):
    """The real Solar System at `date`, ready for solar_system.simulate().

    Returns (state, G):
        state = (names, masses [Msun], pos [AU], vel [AU/yr]) in 3D
        G     = gravitational constant in AU^3 / (Msun * yr^2)

    G comes from NASA's measured G*M_sun, so our Newton engine and NASA's
    ephemeris use exactly the same physics constants.
    """
    load_kernels()
    names = list(names or SPICE_ID.keys())
    et = to_et(date)

    G = gm_au_yr("Sun")                                   # because M_sun = 1
    masses = np.array([gm_au_yr(n) / G for n in names])

    pos = np.zeros((len(names), 3))
    vel = np.zeros((len(names), 3))
    for i, name in enumerate(names):
        s, _ = spice.spkezr(str(SPICE_ID[name]), et, FRAME, "NONE", ORIGIN)
        pos[i] = np.array(s[:3]) / AU_KM                  # km   -> AU
        vel[i] = np.array(s[3:]) / AU_KM * YEAR_S         # km/s -> AU/yr
    return (names, masses, pos, vel), G


def real_positions(names, et_start, times_yr):
    """Where NASA says each body really is at each time (AU).

    times_yr: years after et_start.  Returns array (len(times), n_bodies, 3).
    """
    load_kernels()
    out = np.zeros((len(times_yr), len(names), 3))
    for k, t in enumerate(times_yr):
        et = et_start + t * YEAR_S
        for i, name in enumerate(names):
            p, _ = spice.spkpos(str(SPICE_ID[name]), et, FRAME, "NONE", ORIGIN)
            out[k, i] = np.array(p) / AU_KM
    return out


def sky_today(date):
    """Distances right now: from the Sun, from Earth, and light travel time."""
    load_kernels()
    et = to_et(date)
    rows = []
    for name in SPICE_ID:
        if name in ("Sun", "Earth"):
            continue
        target = str(SPICE_ID[name])
        from_sun, _ = spice.spkpos(target, et, FRAME, "NONE", "SUN")
        # From Earth itself (399), with light-time correction: we see each
        # planet where it was when the light left it.
        from_earth, lt = spice.spkpos(target, et, FRAME, "LT+S", "399")
        rows.append((name,
                     np.linalg.norm(from_sun) / AU_KM,
                     np.linalg.norm(from_earth) / AU_KM,
                     lt / 60.0))
    return rows


if __name__ == "__main__":
    # Quick check: `python3 nasa_spice.py` downloads the kernels and prints today's sky.
    from datetime import datetime, timezone
    today = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M")
    print(f"The Solar System on {today} UTC, according to NASA JPL (DE440):\n")
    print(f"{'Planet':<9}{'from Sun':>11}{'from Earth':>13}{'light time':>13}")
    for name, sun, earth, minutes in sky_today(today):
        print(f"{name:<9}{sun:>8.3f} AU{earth:>10.3f} AU{minutes:>9.1f} min")
    assert list(SPICE_ID) == list(BODIES), "SPICE_ID and BODIES should list the same bodies"

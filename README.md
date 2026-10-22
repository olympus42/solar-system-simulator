# Solar System Simulator

An N-body gravity simulator that builds the entire Solar System from a single
equation — Newton's law of universal gravitation:

```
F = G * m1 * m2 / r^2
```

No orbits are hard-coded. You give the Sun and eight planets their real masses
and starting distances, and gravity produces the orbits on its own. The
simulator then **verifies itself** by reproducing Kepler's Third Law.

![orbits](orbits.png)

## What it does

- Simulates the Sun + 8 planets as a gravitational N-body system
- Integrates the motion with **velocity-Verlet** (an energy-conserving method)
- Plots the emergent orbits (`orbits.png`)
- Checks energy conservation over 165 years (`energy_conservation.png`)
- Verifies the physics against **Kepler's Third Law** (`T² = a³`)
- Starts from NASA's **real** planet positions (SPICE + JPL DE440) and checks
  itself against NASA, year by year

## The result that proves it works

Running the simulation and measuring each planet's orbital period straight from
the trajectory gives:

```
Planet      a [AU]   T sim [yr]     T^2/a^3
----------------------------------------------------------
Mercury      0.387        0.241      1.0038
Venus        0.723        0.616      1.0030
Earth        1.000        1.002      1.0033
Mars         1.524        1.885      1.0038
Jupiter      5.203       11.884      1.0027
Saturn       9.537       29.377      0.9949
Uranus      19.191       83.749      0.9924
Neptune     30.070      164.352      0.9935
```

Every planet lands on `T²/a³ ≈ 1.0` — Kepler's Third Law, falling out of the
simulation on its own. Earth's year comes out to 1.002 years. That agreement is
the whole point: the code isn't told the orbits, it *discovers* them, and the
laws of astronomy confirm it got them right.

## Run it

```bash
pip install -r requirements.txt
python3 solar_system.py
```

Then open `orbits.png` and `energy_conservation.png`.

## Real data from NASA (SPICE)

The circles above are a model. This part uses the **real** Solar System.

[SPICE](https://naif.jpl.nasa.gov/naif/) is the toolkit NASA's Jet Propulsion
Laboratory uses to know where everything is: Cassini, Juno, New Horizons and the
Mars rovers all navigate with it. `nasa_spice.py` asks SPICE (through
[SpiceyPy](https://github.com/AndrewAnnex/SpiceyPy)) where each planet actually
is, using JPL's **DE440** ephemeris.

```bash
pip install -r requirements.txt
python3 real_solar_system.py          # start today, run 20 years
python3 animate.py --real             # animation starting from today's sky
python3 nasa_spice.py                 # just print today's distances
```

The first run downloads three NASA files ("kernels", about 32 MB) into `kernels/`:

| Kernel | What it holds |
|---|---|
| `naif0012.tls` | Leap seconds, to turn a calendar date into SPICE time |
| `de440s.bsp` | Positions and velocities of the Sun and planets, 1849–2150 |
| `gm_de440.tpc` | G × mass of the Sun and every planet |

### What it does

1. Prints where every planet is **today**: distance from the Sun, from Earth, and
   how long its light takes to reach you.
2. Draws the real sky (`real_sky.png`).
3. Starts **your own Newton engine** from NASA's exact positions and velocities,
   runs it forward 20 years, and measures how far each planet ends up from where
   NASA says it really is (`newton_vs_nasa.png`).

### What we found

Example run starting 22 October 2026, time step 0.37 days:

```
Planet       after 1 yr    after 5 yr   after 10 yr   after 20 yr
Mercury      287,133 km     1.97 M km     3.37 M km     6.61 M km
Venus         42,956 km    186,521 km    374,301 km    763,857 km
Earth         12,507 km     62,504 km    125,594 km    250,466 km
Mars           2,309 km     12,826 km     25,804 km     49,426 km
Jupiter           17 km        129 km        423 km        746 km
Saturn             2 km         11 km         30 km        111 km
Uranus             0 km          1 km          3 km          9 km
Neptune            0 km          2 km          7 km         30 km
```

A few hundred lines of Python, one equation from 1687, and Neptune is within
30 km of NASA's position after 20 years.

**Why is Mercury the worst?** Halve the time step and every error drops by
almost exactly **4×** (Mercury after 5 years: 8.0 → 2.0 → 0.50 → 0.12 million km
for steps of 0.73, 0.37, 0.18 and 0.09 days). Error ∝ dt² is the fingerprint of
velocity-Verlet, a *second-order* method. So the physics is right; the step is
too coarse for Mercury, which circles the Sun every 88 days. The outer planets
move so slowly that the same step is plenty.

Shrink the step far enough and the remaining gap is physics Newton didn't
include: the asteroid belt and Pluto pulling on Mars, and **general relativity**
bending Mercury's orbit (43 arcseconds per century, the effect Einstein
explained in 1915).

Note: "Earth" here is the Earth–Moon barycenter, the point the two orbit
around, which is what the Sun's gravity actually moves. Same for planets with
moons.

## Make it yours (start here)

The best way to learn is to change things and see what happens:

1. **Change the time step.** In `solar_system.py`, set `DT = 0.05` and rerun.
   Watch the energy-conservation plot get worse — you'll *see* why numerical
   method choice matters.
2. **Add a moon.** Add an entry to `BODIES` and give it a starting position near
   a planet. Can you get it to orbit the planet instead of the Sun?
3. **Break Jupiter.** Multiply Jupiter's mass by 1000 and watch it disturb the
   inner planets. This is literally how astronomers hunt for unseen planets.
4. **Add a comet** on a highly elliptical path (start it slow and far out).
5. **Animate it.** Done: `python3 animate.py`. Try `--real` to start from today.
6. **Go back in time.** `python3 real_solar_system.py --date 1977-09-05` starts
   the day Voyager 1 launched.
7. **Beat Mercury.** Find the largest `--dt` that keeps Mercury within one
   Earth diameter (12,742 km) of NASA after 1 year.
8. **Add Pluto.** It's SPICE id `9`. Add it to `BODIES` and `SPICE_ID`.
9. **Add Einstein.** Add the general-relativity correction to Mercury's
   acceleration and watch the remaining error shrink.

Each of these teaches a real concept — numerical stability, orbital mechanics,
perturbation. That's Phase 1 in action.

## Why this project

It's the first tool of a longer plan: learn to code by pointing it at physics
you already understand, ship something real and public, and build the habit of
verifying your results against known laws. Newton's gravity → Kepler's laws is
the perfect first proof that your code does real science.

## License

MIT — see `LICENSE`.

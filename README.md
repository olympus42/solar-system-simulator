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
5. **Animate it.** Use `matplotlib.animation` to watch the planets move in real
   time.

Each of these teaches a real concept — numerical stability, orbital mechanics,
perturbation. That's Phase 1 in action.

## Why this project

It's the first tool of a longer plan: learn to code by pointing it at physics
you already understand, ship something real and public, and build the habit of
verifying your results against known laws. Newton's gravity → Kepler's laws is
the perfect first proof that your code does real science.

## License

MIT — see `LICENSE`.

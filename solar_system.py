"""
solar_system.py
================
A real N-body gravity simulator for the Solar System.

Everything here comes from ONE equation you already know -- Newton's law of
universal gravitation:

        F = G * m1 * m2 / r^2

We give the Sun and the eight planets their real masses and starting distances,
then let gravity do the rest. No orbits are hard-coded: the ellipses you see are
an *emergent* result of integrating Newton's law step by step.

We use a clean unit system standard in celestial mechanics:
    - distance  in AU              (1 AU = Earth-Sun distance)
    - time      in years
    - mass      in solar masses    (Sun = 1)
In these units the gravitational constant is exactly G = 4 * pi^2, which is a
lovely consequence of Earth taking 1 year to travel 1 AU-radius orbit.

Author: Starman
License: MIT
"""

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# 1. Physical constants and initial conditions (real astronomical data)
# ----------------------------------------------------------------------

G = 4.0 * np.pi**2          # gravitational constant in AU^3 / (Msun * yr^2)

# name: (mass [Msun], semi-major axis a [AU], colour, marker size)
BODIES = {
    "Sun":     (1.0,        0.000,  "#FDB813", 10),
    "Mercury": (1.651e-7,   0.387,  "#8c8c8c", 3),
    "Venus":   (2.447e-6,   0.723,  "#e6c07b", 5),
    "Earth":   (3.003e-6,   1.000,  "#3a6ea5", 5),
    "Mars":    (3.213e-7,   1.524,  "#c1440e", 4),
    "Jupiter": (9.543e-4,   5.203,  "#d8a47f", 8),
    "Saturn":  (2.857e-4,   9.537,  "#e3d9b0", 7),
    "Uranus":  (4.365e-5,  19.191,  "#a7d3d4", 6),
    "Neptune": (5.150e-5,  30.070,  "#4b70dd", 6),
}


def initial_state():
    """Build starting positions and velocities.

    Each planet starts on the +x axis at its semi-major axis distance, moving
    in the +y direction at the circular orbital speed v = sqrt(G * M_sun / a).
    The Sun is given a small velocity so the total momentum of the system is
    zero (otherwise the whole solar system would slowly drift across the plot).
    """
    names = list(BODIES.keys())
    masses = np.array([BODIES[n][0] for n in names])
    n = len(names)

    pos = np.zeros((n, 2))   # (x, y) for each body, in AU
    vel = np.zeros((n, 2))   # (vx, vy) for each body, in AU/yr

    for i, name in enumerate(names):
        a = BODIES[name][1]
        if a == 0.0:
            continue                       # the Sun starts at the origin
        pos[i] = [a, 0.0]
        speed = np.sqrt(G * masses[0] / a)  # circular speed around the Sun
        vel[i] = [0.0, speed]

    # Zero the total momentum by giving the Sun a recoil velocity.
    total_momentum = (masses[:, None] * vel).sum(axis=0)
    vel[0] = -total_momentum / masses[0]

    return names, masses, pos, vel


# ----------------------------------------------------------------------
# 2. The physics engine
# ----------------------------------------------------------------------

def accelerations(pos, masses, G=G):
    """Return the gravitational acceleration on every body from every other.

    This is the heart of the simulator. For each pair of bodies we apply
    Newton's law of gravitation. The acceleration on body i is:

        a_i = sum over j != i of  G * m_j * (r_j - r_i) / |r_j - r_i|^3

    Works in 2D or 3D: `pos` has one row per body and one column per axis.
    """
    n = len(masses)
    acc = np.zeros_like(pos)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = pos[j] - pos[i]                 # vector from i to j
            r = np.sqrt(d @ d) + 1e-12          # distance (+tiny value avoids /0)
            acc[i] += G * masses[j] * d / r**3
    return acc


def simulate(dt, n_steps, state=None, G=G):
    """Integrate the whole system forward using velocity-Verlet.

    Velocity-Verlet is used instead of naive Euler because it conserves energy
    over long runs -- essential for a physically believable orbit. We record the
    full trajectory so we can plot it and check the physics afterwards.

    By default the planets start on circles (see `initial_state`). Pass
    `state=(names, masses, pos, vel)` to start from anywhere else -- for example
    the real positions NASA measures, loaded by `nasa_spice.py`.
    """
    names, masses, pos, vel = state if state is not None else initial_state()
    n, dims = pos.shape

    trajectory = np.zeros((n_steps, n, dims))
    energy = np.zeros(n_steps)

    acc = accelerations(pos, masses, G)
    for step in range(n_steps):
        # velocity-Verlet update
        pos = pos + vel * dt + 0.5 * acc * dt**2
        new_acc = accelerations(pos, masses, G)
        vel = vel + 0.5 * (acc + new_acc) * dt
        acc = new_acc

        trajectory[step] = pos
        energy[step] = total_energy(pos, vel, masses, G)

    return names, masses, trajectory, energy


def total_energy(pos, vel, masses, G=G):
    """Total mechanical energy (kinetic + potential). Should stay ~constant."""
    kinetic = 0.5 * np.sum(masses * np.sum(vel**2, axis=1))
    potential = 0.0
    n = len(masses)
    for i in range(n):
        for j in range(i + 1, n):
            r = np.sqrt(np.sum((pos[j] - pos[i])**2)) + 1e-12
            potential -= G * masses[i] * masses[j] / r
    return kinetic + potential


# ----------------------------------------------------------------------
# 3. Verification -- does the simulation obey Kepler's Third Law?
# ----------------------------------------------------------------------

def measure_periods(names, trajectory, dt):
    """Measure each planet's orbital period by timing a full revolution.

    We track the angle of each planet around the Sun and detect when it wraps
    past 2*pi. This gives a period measured *from the simulation itself*, which
    we then compare against Kepler's Third Law.
    """
    periods = {}
    for i, name in enumerate(names):
        if name == "Sun":
            continue
        x = trajectory[:, i, 0]
        y = trajectory[:, i, 1]
        angle = np.unwrap(np.arctan2(y, x))       # continuous angle in radians
        total_turn = np.abs(angle[-1] - angle[0])
        if total_turn < 2 * np.pi:
            periods[name] = None                  # didn't finish one orbit
            continue
        revolutions = total_turn / (2 * np.pi)
        periods[name] = (len(x) * dt) / revolutions
    return periods


def kepler_check(names, periods):
    """Kepler's Third Law: T^2 = a^3 (in years and AU). Print the comparison."""
    print("\nKepler's Third Law check   ( expect  T^2 / a^3  ~ 1.0 )")
    print("-" * 58)
    print(f"{'Planet':<9}{'a [AU]':>9}{'T sim [yr]':>13}{'T^2/a^3':>12}")
    print("-" * 58)
    for name in names:
        if name == "Sun" or periods.get(name) is None:
            continue
        a = BODIES[name][1]
        T = periods[name]
        ratio = T**2 / a**3
        print(f"{name:<9}{a:>9.3f}{T:>13.3f}{ratio:>12.4f}")
    print("-" * 58)
    print("If the last column is ~1.0, your simulator reproduces real physics.")


# ----------------------------------------------------------------------
# 4. Plots
# ----------------------------------------------------------------------

def plot_orbits(names, trajectory, filename="orbits.png"):
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.set_facecolor("#05060f")
    for i, name in enumerate(names):
        colour = BODIES[name][2]
        x = trajectory[:, i, 0]
        y = trajectory[:, i, 1]
        if name == "Sun":
            ax.plot(x[-1], y[-1], "o", color=colour, markersize=BODIES[name][3],
                    label=name, zorder=5)
        else:
            ax.plot(x, y, "-", color=colour, linewidth=0.8, alpha=0.8)
            ax.plot(x[-1], y[-1], "o", color=colour,
                    markersize=BODIES[name][3], label=name, zorder=5)
    ax.set_aspect("equal")
    ax.set_title("Solar System  —  orbits emerging from Newton's law of gravity",
                 color="white")
    ax.set_xlabel("x [AU]", color="white")
    ax.set_ylabel("y [AU]", color="white")
    ax.tick_params(colors="white")
    for s in ax.spines.values():
        s.set_color("#333333")
    ax.legend(loc="upper right", facecolor="#111111", labelcolor="white",
              fontsize=8, framealpha=0.8)
    fig.tight_layout()
    fig.savefig(filename, dpi=130, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"saved {filename}")


def plot_energy(energy, dt, filename="energy_conservation.png"):
    t = np.arange(len(energy)) * dt
    drift = (energy - energy[0]) / abs(energy[0])
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(t, drift, color="#3a6ea5")
    ax.set_title("Energy conservation  (a good integrator keeps this ~flat)")
    ax.set_xlabel("time [years]")
    ax.set_ylabel("relative energy drift")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(filename, dpi=130)
    plt.close(fig)
    print(f"saved {filename}")


# ----------------------------------------------------------------------
# 5. Run it
# ----------------------------------------------------------------------

if __name__ == "__main__":
    # ~165 years so even Neptune completes a full orbit. Small dt keeps
    # fast-moving Mercury accurate. Feel free to change these two numbers!
    DT = 0.002          # time step in years
    N_STEPS = 82500     # 0.002 * 82500 = 165 years

    print("Simulating the Solar System from Newton's law of gravity...")
    names, masses, trajectory, energy = simulate(DT, N_STEPS)

    plot_orbits(names, trajectory)
    plot_energy(energy, DT)

    periods = measure_periods(names, trajectory, DT)
    kepler_check(names, periods)

    print("\nDone. Open orbits.png and energy_conservation.png to see your work.")

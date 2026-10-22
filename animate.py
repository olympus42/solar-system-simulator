"""
animate.py  —  watch the planets move.
======================================

This reuses the exact same physics engine from `solar_system.py` (Newton's law
of gravity + velocity-Verlet) and turns the trajectory into a live animation.

Two ways to run it:

    python3 animate.py            ->  saves an animated GIF (solar_system.gif)
    python3 animate.py --show     ->  opens a live window on your computer
    python3 animate.py --real     ->  starts from where the planets REALLY are
                                      today, using NASA's SPICE data
                                      (combine with --show if you like)

Everything you might want to tweak is in the SETTINGS block below.
"""

import sys
from datetime import datetime, timedelta, timezone

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

from solar_system import simulate, BODIES   # <-- reusing your own engine

# ----------------------------------------------------------------------
# SETTINGS  (play with these!)
# ----------------------------------------------------------------------
YEARS  = 12       # how much time to animate
DT     = 0.005    # simulation time step in years (smaller = more accurate)
FRAMES = 300      # number of animation frames (more = smoother, bigger file)
TRAIL  = 45       # how long a glowing tail each planet leaves behind
VIEW   = 6.0      # half-width of the view in AU. Try 1.8 for the inner planets,
                  # or 32 to fit Neptune's whole orbit.
FPS    = 30       # playback speed of the saved GIF


def main(show=False, real=False):
    n_steps = int(YEARS / DT)
    start_date = datetime.now(timezone.utc)
    if real:
        import nasa_spice   # only needed for --real (pip install spiceypy)
        state, G = nasa_spice.real_state(start_date.strftime("%Y-%m-%d %H:%M"))
        print(f"Starting from NASA's real positions on {start_date:%Y-%m-%d}.")
        print(f"Simulating {YEARS} years ({n_steps} steps)...")
        names, masses, trajectory, energy = simulate(DT, n_steps, state=state, G=G)
    else:
        print(f"Simulating {YEARS} years ({n_steps} steps)...")
        names, masses, trajectory, energy = simulate(DT, n_steps)

    # Choose evenly spaced frames out of the full, high-resolution trajectory.
    frame_idx = np.linspace(0, n_steps - 1, FRAMES).astype(int)
    gap = max(1, n_steps // FRAMES)   # sim-steps between two frames (for trails)

    # --- set up the figure ---
    fig, ax = plt.subplots(figsize=(8, 8))
    fig.patch.set_facecolor("#05060f")
    ax.set_facecolor("#05060f")
    ax.set_xlim(-VIEW, VIEW)
    ax.set_ylim(-VIEW, VIEW)
    ax.set_aspect("equal")
    title = "from today's real positions (NASA SPICE)" if real else "live, from Newton's law of gravity"
    ax.set_title(f"Solar System — {title}", color="white")
    ax.set_xlabel("x [AU]", color="white")
    ax.set_ylabel("y [AU]", color="white")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_color("#333333")

    # One dot + one trailing line per body.
    dots, trails = [], []
    for name in names:
        colour, size = BODIES[name][2], BODIES[name][3]
        (dot,)   = ax.plot([], [], "o", color=colour, markersize=size,
                           zorder=5, label=name)
        (trail,) = ax.plot([], [], "-", color=colour, linewidth=0.8, alpha=0.5)
        dots.append(dot)
        trails.append(trail)

    ax.legend(loc="upper right", facecolor="#111111", labelcolor="white",
              fontsize=7, framealpha=0.7)
    clock = ax.text(0.02, 0.97, "", transform=ax.transAxes, color="white",
                    va="top", fontsize=10, family="monospace")

    def update(i):
        f = frame_idx[i]
        start = max(0, f - TRAIL * gap)
        for k, (dot, trail) in enumerate(zip(dots, trails)):
            dot.set_data([trajectory[f, k, 0]], [trajectory[f, k, 1]])
            trail.set_data(trajectory[start:f + 1, k, 0],
                           trajectory[start:f + 1, k, 1])
        if real:
            clock.set_text((start_date + timedelta(days=(f + 1) * DT * 365.25)).strftime("%Y-%m-%d"))
        else:
            clock.set_text(f"t = {f * DT:6.2f} yr")
        return dots + trails + [clock]

    anim = FuncAnimation(fig, update, frames=FRAMES, interval=1000 // FPS,
                         blit=True)

    if show:
        print("Opening live window — close it to finish.")
        plt.show()
    else:
        out = "solar_system_real.gif" if real else "solar_system.gif"
        anim.save(out, writer=PillowWriter(fps=FPS))
        print(f"saved {out}")


if __name__ == "__main__":
    main(show="--show" in sys.argv, real="--real" in sys.argv)

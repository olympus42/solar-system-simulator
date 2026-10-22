"""
real_solar_system.py  —  your Newton engine vs. NASA.
=====================================================

solar_system.py starts the planets on perfect circles and lets gravity run.
This script does the real thing:

  1. Asks NASA's SPICE Toolkit where every planet ACTUALLY is today
     (JPL's DE440 ephemeris, the same data missions navigate with).
  2. Starts your own Newton engine from those exact positions and speeds.
  3. Runs it forward for years, then checks every planet against where
     NASA says it really is. How far off is Newton, and why?

Run it:

    python3 real_solar_system.py                 # start today, run 20 years
    python3 real_solar_system.py --date 1977-09-05 --years 40
    python3 real_solar_system.py --dt 0.0005     # smaller step, more accurate

First run downloads ~32 MB of NASA data into kernels/.
Outputs: real_sky.png and newton_vs_nasa.png
"""

import argparse
from datetime import datetime, timezone

import matplotlib.pyplot as plt
import numpy as np

import nasa_spice as nasa
from solar_system import BODIES, simulate

AU_KM = nasa.AU_KM


def parse_args():
    p = argparse.ArgumentParser(description="Newton engine vs. NASA JPL ephemeris")
    p.add_argument("--date", default=datetime.now(timezone.utc).strftime("%Y-%m-%d"),
                   help="start date, UTC (default: today). DE440s covers 1849-2150.")
    p.add_argument("--years", type=float, default=20, help="how long to simulate")
    p.add_argument("--dt", type=float, default=0.001, help="time step in years")
    p.add_argument("--checks", type=int, default=400, help="how many times to compare with NASA")
    return p.parse_args()


# ----------------------------------------------------------------------
# 1. Today's sky, straight from NASA
# ----------------------------------------------------------------------

def print_sky(date):
    print(f"\nThe Solar System on {date} (UTC), from NASA JPL's DE440 ephemeris")
    print("-" * 52)
    print(f"{'Planet':<9}{'from Sun':>11}{'from Earth':>14}{'light time':>13}")
    print("-" * 52)
    for name, sun, earth, minutes in nasa.sky_today(date):
        light = f"{minutes:.1f} min" if minutes < 60 else f"{minutes / 60:.2f} h"
        print(f"{name:<9}{sun:>8.3f} AU{earth:>11.3f} AU{light:>13}")
    print("-" * 52)
    print("Light time = how old the light is when you see that planet tonight.")


def plot_real_sky(names, et, date, filename="real_sky.png"):
    """Every planet where NASA says it is on `date`, with its real orbit."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    for ax, view, title in [(axes[0], 1.8, "Inner planets"), (axes[1], 32, "Whole Solar System")]:
        ax.set_facecolor("#05060f")
        for i, name in enumerate(names):
            colour, size = BODIES[name][2], BODIES[name][3]
            if name != "Sun":
                # one full orbit, sampled from NASA's data (Kepler: T = a^1.5 years)
                period = BODIES[name][1] ** 1.5
                times = np.linspace(-period, 0, 240)
                path = nasa.real_positions([name], et, times)[:, 0]
                ax.plot(path[:, 0], path[:, 1], "-", color=colour, lw=0.8, alpha=0.6)
            now = nasa.real_positions([name], et, [0.0])[0, 0]
            ax.plot(now[0], now[1], "o", color=colour, markersize=size, label=name, zorder=5)
        ax.set_xlim(-view, view)
        ax.set_ylim(-view, view)
        ax.set_aspect("equal")
        ax.set_title(f"{title} on {date}", color="white")
        ax.set_xlabel("x [AU]  (toward the March equinox)", color="white")
        ax.set_ylabel("y [AU]", color="white")
        ax.tick_params(colors="white")
        for s in ax.spines.values():
            s.set_color("#333333")
    axes[1].legend(loc="upper right", facecolor="#111111", labelcolor="white", fontsize=8)
    fig.patch.set_facecolor("#05060f")
    fig.suptitle("Real positions from NASA SPICE (JPL DE440)", color="white")
    fig.tight_layout()
    fig.savefig(filename, dpi=120, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    print(f"saved {filename}")


# ----------------------------------------------------------------------
# 2. Newton vs. NASA
# ----------------------------------------------------------------------

def compare(date, years, dt, checks):
    state, G = nasa.real_state(date)
    names = state[0]
    n_steps = int(round(years / dt))
    print(f"\nRunning your Newton engine from NASA's real starting state:")
    print(f"  {years:g} years, dt = {dt:g} yr ({dt * 365.25:.2f} days), {n_steps:,} steps ...")
    _, _, trajectory, energy = simulate(dt, n_steps, state=state, G=G)

    # Compare at `checks` evenly spaced moments (trajectory[k] is time (k+1)*dt)
    idx = np.unique(np.linspace(0, n_steps - 1, checks).astype(int))
    times = (idx + 1) * dt
    truth = nasa.real_positions(names, nasa.to_et(date), times)
    error_km = np.linalg.norm(trajectory[idx] - truth, axis=2) * AU_KM
    return names, times, error_km, energy


def print_report(names, times, error_km):
    marks = [m for m in (1, 5, 10, 20, 50, 100) if m <= times[-1] + 1e-9]
    print("\nHow far is your simulation from where NASA says each planet really is?")
    header = f"{'Planet':<9}" + "".join(f"{f'after {m} yr':>14}" for m in marks)
    print("-" * len(header))
    print(header)
    print("-" * len(header))
    for i, name in enumerate(names):
        if name == "Sun":
            continue
        cells = []
        for m in marks:
            k = int(np.argmin(np.abs(times - m)))
            cells.append(f"{fmt_km(error_km[k, i]):>14}")
        print(f"{name:<9}" + "".join(cells))
    print("-" * len(header))
    print("For scale: the Moon is 384,400 km from Earth. Earth is 12,742 km wide.")


def fmt_km(km):
    return f"{km:,.0f} km" if km < 1e6 else f"{km / 1e6:,.2f} M km"


def plot_comparison(names, times, error_km, date, dt, filename="newton_vs_nasa.png"):
    fig, ax = plt.subplots(figsize=(10, 5.5))
    for i, name in enumerate(names):
        if name == "Sun":
            continue
        ax.plot(times, error_km[:, i], color=BODIES[name][2], lw=1.4, label=name)
    ax.axhline(384_400, color="#888888", ls="--", lw=0.8)
    ax.text(times[-1], 384_400 * 0.7, "Earth–Moon distance", color="#666666",
            fontsize=8, ha="right", va="top")
    ax.set_yscale("log")
    ax.set_ylim(bottom=1)
    ax.set_xlabel(f"years after {date}")
    ax.set_ylabel("distance from NASA's real position [km]")
    ax.set_title(f"Your Newton engine vs. NASA JPL DE440   (dt = {dt * 365.25:.2f} days)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(filename, dpi=130)
    plt.close(fig)
    print(f"saved {filename}")


# ----------------------------------------------------------------------
# 3. Run it
# ----------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()
    nasa.load_kernels()
    et = nasa.to_et(args.date)

    print_sky(args.date)
    plot_real_sky(list(nasa.SPICE_ID), et, args.date)

    names, times, error_km, energy = compare(args.date, args.years, args.dt, args.checks)
    print_report(names, times, error_km)
    plot_comparison(names, times, error_km, args.date, args.dt)

    drift = abs(energy[-1] - energy[0]) / abs(energy[0])
    print(f"\nEnergy drift over the run: {drift:.1e}  (velocity-Verlet keeps this tiny)")
    print("Open real_sky.png and newton_vs_nasa.png to see the results.")

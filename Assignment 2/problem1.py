"""
Problem 1: Solitary wave propagation and grid effects
Runs sbouss 3 times with dx=2, 1, 0.5 and plots eta at t=70.
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os

FOLDER = os.path.dirname(os.path.abspath(__file__))
SBOUSS = os.path.join(FOLDER, "sbouss")

# Parameters
L = 150       # domain length (wave starts at x=20, travels ~73 units in t=70)
alpha = 0.1   # amplitude
x0 = 20       # initial position
t_end = 70    # simulation end time
dt_factor = 0.5  # reduction factor for dt/dx

dx_values = [2, 1, 0.5]


def run_sbouss(dx, output_dir):
    """Run sbouss with given dx, saving output files to output_dir."""
    os.makedirs(output_dir, exist_ok=True)
    N = int(L / dx)  # number of grid points

    # Build the interactive input for sbouss
    inputs = "\n".join([
        "flat",           # bathimetry
        str(float(L)),    # total length
        str(N),           # number of grid points
        "Boussinesq",     # equation type (appropriate for solitary waves)
        "no",             # discrete correction term: off as required
        "soliton",        # generation condition
        str(alpha),       # a/h amplitude
        str(float(x0)),   # initial position
        "no",             # propagate towards decreasing x? no = rightward
        "1",              # number of cycles (1 output at t=t_end)
        str(float(t_end)),# time interval
        str(dt_factor),   # reduction factor for dt/dx
        "all",            # times for printing
        "",               # (extra newline in case needed)
    ])

    result = subprocess.run(
        [SBOUSS],
        input=inputs,
        capture_output=True,
        text=True,
        cwd=output_dir,
    )

    if result.returncode != 0:
        print(f"sbouss error for dx={dx}:")
        print(result.stderr)
    else:
        print(f"dx={dx}: done")


def read_eta(filepath):
    """Read a two-column sbouss output file (x, eta)."""
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, 1]


# --- Run simulations ---
print("Running simulations...")
for dx in dx_values:
    out_dir = os.path.join(FOLDER, "data", "problem1", f"run_dx{dx}".replace(".", "p"))
    run_sbouss(dx, out_dir)

# --- Plot results ---
fig, ax = plt.subplots(figsize=(10, 5))

for dx in dx_values:
    out_dir = os.path.join(FOLDER, "data", "problem1", f"run_dx{dx}".replace(".", "p"))
    eta_file = os.path.join(out_dir, "eta1")  # eta1 = output at t=t_end
    if os.path.exists(eta_file):
        x, eta = read_eta(eta_file)
        ax.plot(x, eta, label=f"$\\Delta x = {dx}$", linewidth=2)
        print(f"dx={dx}: max eta = {eta.max():.6f}")
    else:
        print(f"dx={dx}: output file eta1 not found in {out_dir}")

ax.set_xlabel("x", fontsize=16)
ax.set_ylabel("$\\eta$", fontsize=16)
ax.set_title(f"Solitary wave at $\\mathbf{{t = {t_end}}}$, $\\mathbf{{\\alpha = {alpha}}}$", fontsize=18, fontweight='bold')
ax.legend(fontsize=16)
ax.tick_params(axis='both', labelsize=16)
ax.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(FOLDER, "results", "problem1", "problem1.pdf"), dpi=150)
plt.show()

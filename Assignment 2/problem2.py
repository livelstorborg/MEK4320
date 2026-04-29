"""
Problem 2: Wave dispersion
Initial elevation: eta(x,0) = 2A*cos^2(pi*(x-x0)/lam) for -lam/2 <= x-x0 <= lam/2, else 0
u(x,0) = 0

Parameters: A=0.2, lam=20, x0=0, L=100, t=80, h=1
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os

FOLDER = os.path.dirname(os.path.abspath(__file__))
SBOUSS = os.path.join(FOLDER, "sbouss")

# Common parameters
A = 0.2
lam = 20.0
x0 = 0.0
L = 100.0
t_end = 80.0


def make_eta_in(filepath, dx_fine=0.1):
    """Generate eta.in with the initial cos^2 elevation profile."""
    x = np.arange(0, L + dx_fine, dx_fine)
    eta = np.zeros_like(x)
    mask = np.abs(x - x0) <= lam / 2
    eta[mask] = 2 * A * np.cos(np.pi * (x[mask] - x0) / lam) ** 2
    np.savetxt(filepath, np.column_stack([x, eta]), fmt="%.6f")


def make_u_in(filepath, dx_fine=0.1):
    """Generate u1.in with zero initial velocity."""
    x = np.arange(0, L + dx_fine, dx_fine)
    np.savetxt(filepath, np.column_stack([x, np.zeros_like(x)]), fmt="%.6f")


def run_sbouss(eq_type, run_name, dx=0.5, dt_factor=0.5, t_sim=None):
    """Run sbouss for a given equation type, returns output directory path."""
    if t_sim is None:
        t_sim = t_end

    output_dir = os.path.join(FOLDER, "data", "problem2", f"run_p2_{run_name}")
    os.makedirs(output_dir, exist_ok=True)
    N = int(L / dx)

    for fname in ["eta.in", "u1.in"]:
        src = os.path.join(FOLDER, "data", "problem2", fname)
        dst = os.path.join(output_dir, fname)
        with open(src) as f:
            content = f.read()
        with open(dst, "w") as f:
            f.write(content)

    inputs = "\n".join([
        "flat",
        str(float(L)),
        str(N),
        eq_type,
        "yes",
        "readfromfile",
        "eta.in",
        "u1.in",
        "1",
        str(float(t_sim)),
        str(dt_factor),
        "all",
        "",
    ])

    result = subprocess.run(
        [SBOUSS], input=inputs, capture_output=True, text=True, cwd=output_dir
    )
    if result.returncode != 0:
        print(f"sbouss error ({run_name}): {result.stderr[-300:]}")
    else:
        print(f"{run_name}: done")
    return output_dir


def read_eta(filepath):
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, 1]


def setup_plot(title):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_xlabel("x", fontsize=18)
    ax.set_ylabel("$\\eta$", fontsize=18)
    ax.set_title(title, fontsize=22, fontweight='bold')
    ax.tick_params(axis='both', labelsize=16)
    ax.grid(True)
    return fig, ax


# =============================================================================
def problem2a():
    """
    a) LSW equations with different resolutions.
    For LSW, the result should be two non-dispersive waves travelling in
    opposite directions at speed c=1. Adjust resolution until this is clear.
    """
    print("\n--- Problem 2a: LSW ---")
    dx_values = [1.0, 0.5, 0.25]

    fig, ax = setup_plot(f"LSW at $\\mathbf{{t = {t_end}}}$")
    for dx in dx_values:
        out_dir = run_sbouss("LSW", f"a_dx{str(dx).replace('.','p')}", dx=dx)
        eta_file = os.path.join(out_dir, "eta1")
        if os.path.exists(eta_file):
            x, eta = read_eta(eta_file)
            ax.plot(x, eta, label=f"$\\Delta x = {dx}$", linewidth=2)

    ax.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(FOLDER, "results", "problem2", "problem2a.pdf"))
    plt.show()
    print("Saved problem2a.pdf")


# =============================================================================
def problem2b():
    """
    b) Linearized Boussinesq equations. Plot at t=80.
    Dispersive effects cause wave packet to spread.
    """
    print("\n--- Problem 2b: LBoussinesq ---")
    dx = 0.5
    out_dir = run_sbouss("LBoussinesq", f"b_LBoussinesq_dx{str(dx).replace('.','p')}", dx=dx)

    fig, ax = setup_plot(f"LBoussinesq at $\\mathbf{{t = {t_end}}}$")
    eta_file = os.path.join(out_dir, "eta1")
    if os.path.exists(eta_file):
        x, eta = read_eta(eta_file)
        ax.plot(x, eta, linewidth=2)
    plt.tight_layout()
    plt.savefig(os.path.join(FOLDER, "results", "problem2", "problem2b.pdf"))
    plt.show()
    print("Saved problem2b.pdf")


# =============================================================================
def problem2c():
    """
    c) Full Boussinesq vs linearized Boussinesq. Plot together.
    Nonlinear effects may cause leading wave to steepen into a solitary wave.
    """
    print("\n--- Problem 2c: Boussinesq vs LBoussinesq ---")
    dx = 0.5
    out_dir_b = os.path.join(FOLDER, "data", "problem2", f"run_p2_b_LBoussinesq_dx{str(dx).replace('.','p')}")
    out_dir_c = run_sbouss("Boussinesq", f"c_Boussinesq_dx{str(dx).replace('.','p')}", dx=dx)

    fig, ax = setup_plot(f"Boussinesq vs LBoussinesq at $\\mathbf{{t = {t_end}}}$")
    for out_dir, label in [(out_dir_b, "LBoussinesq"), (out_dir_c, "Boussinesq")]:
        eta_file = os.path.join(out_dir, "eta1")
        if os.path.exists(eta_file):
            x, eta = read_eta(eta_file)
            ax.plot(x, eta, label=label, linewidth=2)

    ax.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(FOLDER, "results", "problem2", "problem2c.pdf"))
    plt.show()
    print("Saved problem2c.pdf")


# =============================================================================
def problem2d():
    """
    d) NLSW with grid refinement. Short features that change with refinement
    are artifacts (invalid solution).
    """
    print("\n--- Problem 2d: NLSW grid refinement ---")
    dx_values = [1.0, 0.5, 0.25]
    alphas = [1.0, 0.75, 0.5]

    fig, ax = setup_plot(f"NLSW at $\\mathbf{{t = {t_end}}}$ — grid refinement")
    for dx, alpha in zip(dx_values, alphas):
        out_dir = run_sbouss("NLSW", f"d_NLSW_dx{str(dx).replace('.','p')}", dx=dx)
        eta_file = os.path.join(out_dir, "eta1")
        if os.path.exists(eta_file):
            x, eta = read_eta(eta_file)
            ax.plot(x, eta, label=f"$\\Delta x = {dx}$", linewidth=2, alpha=alpha)

    ax.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(FOLDER, "results", "problem2", "problem2d.pdf"))
    plt.show()
    print("Saved problem2d.pdf")


# =============================================================================
def plot_max_steepness_2e(times):
    """Plot maximum wave steepness max|dη/dx| vs time for NLSW runs."""
    print("\n--- Problem 2e: max steepness vs time ---")
    max_steepness = []
    for t in times:
        out_dir = run_sbouss("NLSW", f"e_NLSW_t{int(t*10):03d}", dx=0.5, t_sim=t)
        eta_file = os.path.join(out_dir, "eta1")
        if os.path.exists(eta_file):
            x, eta = read_eta(eta_file)
            dx = x[1] - x[0]
            deta_dx = np.abs(np.gradient(eta, dx))
            max_steepness.append(np.max(deta_dx))
        else:
            max_steepness.append(np.nan)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(times, max_steepness, linewidth=2)
    ax.set_xlabel("$t$", fontsize=18)
    ax.set_ylabel(r"$\max\left| \frac{\partial\eta}{\partial x}\right |$", fontsize=18)
    ax.set_title("NLSW — maximum wave steepness", fontsize=22, fontweight='bold')
    ax.tick_params(axis='both', labelsize=16)
    ax.grid(True)
    ax.axvline(x=3, color='gray', linestyle='--', linewidth=1.5)
    # ax.axvline(x=28, color='gray', linestyle='--', linewidth=1.5)
    plt.tight_layout()
    plt.savefig(os.path.join(FOLDER, "results", "problem2", "problem2e_steepness.pdf"))
    plt.show()
    print("Saved problem2e_steepness.pdf")


def problem2e():
    """
    e) NLSW at early times, before artifacts appear. Shows wave steepening
    (shock formation) — physical nonlinear effect.
    """
    print("\n--- Problem 2e: NLSW early times ---")
    early_times = [10.0, 20.0, 30.0]

    fig, ax = setup_plot("NLSW — Early")
    for t in early_times:
        out_dir = run_sbouss("NLSW", f"e_NLSW_t{int(t)}", dx=0.5, t_sim=t)
        eta_file = os.path.join(out_dir, "eta1")
        if os.path.exists(eta_file):
            x, eta = read_eta(eta_file)
            ax.plot(x, eta, label=f"$t = {t}$", linewidth=2)

    ax.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(FOLDER, "results", "problem2", "problem2e.pdf"))
    plt.show()
    print("Saved problem2e.pdf")

    steepness_times = np.arange(0, 35.0, 0.5)
    plot_max_steepness_2e(steepness_times)


# =============================================================================
def problem2_initial():
    """Plot the initial elevation profile."""
    x_plot = np.linspace(0, L, 1000)
    eta_init = np.where(np.abs(x_plot - x0) <= lam / 2,
                        2 * A * np.cos(np.pi * (x_plot - x0) / lam) ** 2, 0)
    fig, ax = setup_plot("Initial elevation $\\mathbf{\\eta(x,0)}$")
    ax.plot(x_plot, eta_init, 'k', linewidth=2)
    ax.set_xlim(0, lam * 2)
    plt.tight_layout()
    plt.savefig(os.path.join(FOLDER, "results", "problem2", "problem2_initial.pdf"))
    plt.show()
    print("Saved problem2_initial.pdf")


if __name__ == "__main__":
    make_eta_in(os.path.join(FOLDER, "data", "problem2", "eta.in"))
    make_u_in(os.path.join(FOLDER, "data", "problem2", "u1.in"))

    # problem2_initial()
    # problem2a()
    # problem2b()
    # problem2c()
    # problem2d()
    problem2e()

"""
Problem 3: Bores
Same initial elevation as Problem 2 but with L=1000, lam=100.
u(x,0) = 0, constant depth h=1 (except 3c).
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os

FOLDER = os.path.dirname(os.path.abspath(__file__))
SBOUSS = os.path.join(FOLDER, "sbouss")
DATA   = os.path.join(FOLDER, "data", "problem3")

# Common parameters
A   = 0.2
lam = 100.0
x0  = 0.0
L   = 1000.0
dt_factor = 0.5


def make_eta_in(dx_fine=0.5):
    """Generate eta.in with cos^2 initial elevation for problem 3."""
    filepath = os.path.join(DATA, "eta.in")
    x = np.arange(0, L + dx_fine, dx_fine)
    eta = np.zeros_like(x)
    mask = np.abs(x - x0) <= lam / 2
    eta[mask] = 2 * A * np.cos(np.pi * (x[mask] - x0) / lam) ** 2
    np.savetxt(filepath, np.column_stack([x, eta]), fmt="%.6f")
    print(f"Written {filepath}")


def make_u_in(dx_fine=0.5):
    """Generate u1.in with zero initial velocity."""
    filepath = os.path.join(DATA, "u1.in")
    x = np.arange(0, L + dx_fine, dx_fine)
    np.savetxt(filepath, np.column_stack([x, np.zeros_like(x)]), fmt="%.6f")
    print(f"Written {filepath}")


def make_depth_file():
    """
    Generate depth file for problem 3c:
    h=1 for x<40, h=0.2 for x>50, linear slope for 40<=x<=50.
    """
    filepath = os.path.join(DATA, "depth.dat")
    x = np.arange(0, L + 0.1, 0.1)
    h = np.ones_like(x)
    slope = (x >= 40) & (x <= 50)
    h[slope] = 1.0 - (1.0 - 0.2) * (x[slope] - 40) / 10.0
    h[x > 50] = 0.2
    np.savetxt(filepath, np.column_stack([x, h]), fmt="%.6f")
    print(f"Written {filepath}")


def run_sbouss(eq_type, run_name, dx=0.5, t_sim=200.0, n_cycles=1,
               flat=True, depth_file=None, use_soliton=False,
               soliton_alpha=0.05, soliton_x0=20.0):
    """Run sbouss, returns output directory path."""
    output_dir = os.path.join(DATA, f"run_p3_{run_name}")
    os.makedirs(output_dir, exist_ok=True)
    N = int(L / dx)

    # Copy input files
    for fname in ["eta.in", "u1.in"]:
        src = os.path.join(DATA, fname)
        dst = os.path.join(output_dir, fname)
        if os.path.exists(src):
            with open(src) as f:
                content = f.read()
            with open(dst, "w") as f:
                f.write(content)

    if depth_file:
        src = os.path.join(DATA, "depth.dat")
        dst = os.path.join(output_dir, "depth.dat")
        with open(src) as f:
            content = f.read()
        with open(dst, "w") as f:
            f.write(content)

    # Build inputs
    if flat:
        bathimetry_lines = ["flat", str(float(L))]
    else:
        bathimetry_lines = ["readfromfile", "depth.dat", str(float(L))]

    if use_soliton:
        generation_lines = ["soliton", str(soliton_alpha), str(float(soliton_x0)), "no"]
    else:
        generation_lines = ["readfromfile", "eta.in", "u1.in"]

    inputs = "\n".join(
        bathimetry_lines +
        [str(N), eq_type, "yes"] +
        generation_lines +
        [str(n_cycles), str(float(t_sim / n_cycles)), str(dt_factor), "all", ""]
    )

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


# =============================================================================
def plot_max_diff_3a(diff_times, eta_data):
    """Plot max absolute difference between NLSW and Boussinesq vs time."""
    for t in diff_times:
        for eq in ["NLSW", "Boussinesq"]:
            if (eq, t) not in eta_data:
                out_dir = run_sbouss(eq, f"a_{eq}_t{int(t)}", dx=0.5, t_sim=t)
                eta_file = os.path.join(out_dir, "eta1")
                if os.path.exists(eta_file):
                    eta_data[(eq, t)] = read_eta(eta_file)

    t_vals, max_diffs = [], []
    for t in diff_times:
        if ("NLSW", t) in eta_data and ("Boussinesq", t) in eta_data:
            _, eta_n = eta_data[("NLSW", t)]
            _, eta_b = eta_data[("Boussinesq", t)]
            t_vals.append(t)
            max_diffs.append(np.max(np.abs(eta_n - eta_b)))

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t_vals, max_diffs, 'o-', linewidth=2, markersize=5, color='k')
    ax.set_xlabel("$t$", fontsize=16)
    ax.set_ylabel("$\\|\\eta_{NLSW} - \\eta_{Bouss}\\|_{\\infty}$", fontsize=16)
    ax.set_title("Max absolute difference: NLSW vs Boussinesq", fontsize=16, fontweight='bold')
    ax.set_yscale('log')
    ax.tick_params(axis='both', labelsize=14)
    ax.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(FOLDER, "results", "problem3", "problem3a_diff.pdf"))
    plt.show()
    print("Saved problem3a_diff.pdf")


def problem3a():
    """
    a) Compare NLSW and Boussinesq for t<200.
    Same color per time, solid=NLSW, dashed=Boussinesq.
    """
    print("\n--- Problem 3a: NLSW vs Boussinesq, t<200 ---")
    times  = [50.0, 100.0, 150.0, 200.0]
    colors = ['#7B2D8B', '#E07B00', '#1A7ABF', '#2A9D4E']

    diff_times = list(np.arange(10.0, 201.0, 10.0))

    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(12, 5))

    styles = {
        "NLSW":       dict(linestyle="-",  linewidth=2, alpha=0.5),
        "Boussinesq": dict(linestyle="--",  linewidth=2, alpha=1.0),
    }

    nlsw_handles, bouss_handles = [], []
    eta_data = {}  # store (x, eta) keyed by (eq, t)
    for t, color in zip(times, colors):
        for eq, handles in [("NLSW", nlsw_handles), ("Boussinesq", bouss_handles)]:
            out_dir = run_sbouss(eq, f"a_{eq}_t{int(t)}", dx=0.05, t_sim=t)
            eta_file = os.path.join(out_dir, "eta1")
            if os.path.exists(eta_file):
                x, eta = read_eta(eta_file)
                ax.plot(x, eta, color=color, **styles[eq])
                eta_data[(eq, t)] = (x, eta)
            handles.append(Line2D([0], [0], color=color, label=f"$t={int(t)}$", **styles[eq]))

    legend1 = ax.legend(handles=nlsw_handles, title="NLSW",        fontsize=14, title_fontsize=15, loc='upper right', bbox_to_anchor=(0.5, 1.0))
    legend2 = ax.legend(handles=bouss_handles, title="Boussinesq", fontsize=14, title_fontsize=15, loc='upper left',  bbox_to_anchor=(0.5, 1.0))
    legend1.get_title().set_fontweight('bold')
    legend2.get_title().set_fontweight('bold')
    ax.add_artist(legend1)

    ax.set_xlim(0, 220)
    ax.set_xlabel("x", fontsize=16)
    ax.set_ylabel("$\\eta$", fontsize=16)
    ax.set_title("NLSW vs Boussinesq, $\\mathbf{t < 200}$", fontsize=18, fontweight='bold')
    ax.tick_params(axis='both', labelsize=14)
    ax.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(FOLDER, "results", "problem3", "problem3a.pdf"))
    plt.show()
    print("Saved problem3a.pdf")

    plot_max_diff_3a(diff_times, eta_data)







# =============================================================================
def problem3b():
    """
    b) Boussinesq simulation until t=800.
    Compare second crest to solitary wave solution.
    """
    print("\n--- Problem 3b: Boussinesq until t=800 ---")
    times  = [200.0, 400.0, 600.0, 800.0]
    colors = ['#7B2D8B', '#E07B00', '#1A7ABF', '#2A9D4E']

    fig, ax = plt.subplots(figsize=(12, 5))

    for t, color in zip(times, colors):
        out_dir = run_sbouss("Boussinesq", f"b_Boussinesq_t{int(t)}", dx=0.5, t_sim=t)
        eta_file = os.path.join(out_dir, "eta1")
        if os.path.exists(eta_file):
            x, eta = read_eta(eta_file)
            ax.plot(x, eta, color=color, linewidth=2, label=f"$t={int(t)}$")

    ax.legend(fontsize=12)
    ax.set_xlabel("x", fontsize=16)
    ax.set_ylabel("$\\eta$", fontsize=16)
    ax.set_title("Boussinesq — Soliton trains", fontsize=18, fontweight='bold')
    ax.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(FOLDER, "results", "problem3", "problem3b.pdf"))
    plt.show()
    print("Saved problem3b.pdf")

    # --- Part 2: second crest vs soliton solution at t=800 ---
    from scipy.signal import find_peaks

    out_dir = run_sbouss("Boussinesq", "b_Boussinesq_t800", dx=0.1, t_sim=800.0)
    eta_file = os.path.join(out_dir, "eta1")
    if os.path.exists(eta_file):
        x, eta = read_eta(eta_file)

        peaks, _ = find_peaks(eta, height=0.01)
        peaks = peaks[np.argsort(x[peaks])[::-1]]  # sort by x, front first

        if len(peaks) >= 2:
            idx    = peaks[1]
            alpha_s = eta[idx]
            x0_s    = x[idx]

            # Solitary wave solution
            k     = np.sqrt(3 * alpha_s / 4)
            c     = 1 + alpha_s / 2
            x_sol = np.linspace(x0_s - 100, x0_s + 50, 1000)
            eta_sol = alpha_s / np.cosh(k * (x_sol - x0_s)) ** 2

            fig2, ax2 = plt.subplots(figsize=(10, 5))
            ax2.plot(x, eta, color='#1A7ABF', linewidth=2, label="Boussinesq $t=800$")
            ax2.plot(x_sol, eta_sol, color='#E07B00', linewidth=2, linestyle='--',
                     label=r"$\zeta = \alpha \, \mathrm{sech}^2\!\left(\left(\dfrac{3\alpha}{4}\right)^{1/2} \psi\right)$")
            ax2.set_xlim(x0_s - 100, x0_s + 50)
            ax2.legend(fontsize=12)
            ax2.set_xlabel("x", fontsize=16)
            ax2.set_ylabel("$\\eta$", fontsize=16)
            ax2.set_title("Second crest vs soliton solution at $\mathbf{t=800}$", fontsize=18, fontweight='bold')
            ax2.grid(True)
            plt.tight_layout()
            plt.savefig(os.path.join(FOLDER, "results", "problem3", "problem3b_soliton.pdf"))
            plt.show()
            print(f"Saved problem3b_soliton.pdf  (alpha={alpha_s:.4f}, x0={x0_s:.1f})")


            eta_sol_x = alpha_s / np.cosh(k * (x - x0_s)) ** 2

            fig, (ax3, ax4) = plt.subplots(1, 2, figsize=(10, 5))

            ax3.plot(x, eta, color='#1A7ABF', linewidth=2, label="Boussinesq")
            ax3.plot(x_sol, eta_sol, color='#E07B00', linewidth=2, linestyle='--',
                     label="Soliton Solution")
            ax3.legend(fontsize=12)
            ax3.set_xlim(920, 940)
            ax3.set_xticks(range(920, 941, 5))
            ax3.set_xlabel("x", fontsize=16)
            ax3.set_ylabel("$\\eta$", fontsize=16)
            ax3.grid(alpha=0.5)

            ax4.plot(x, np.abs(eta - eta_sol_x), color='#7B2D8B', linewidth=2)
            ax4.set_xlim(920, 940)
            ax4.set_ylim(0, 0.01)
            ax4.set_xticks(range(920, 941, 5))
            ax4.set_xlabel("x", fontsize=16)
            ax4.set_ylabel("$|\\eta - \\eta_{\\mathrm{sol}}|$", fontsize=16)
            ax4.grid(alpha=0.5)

            fig.suptitle(r"$\mathbf{Second\ crest\ vs\ soliton\ solution\ at\ t=800}$" + "\n" +
                r"IC: $\eta(x,0) = 2A\cos^2\!\left(\frac{\pi(x-x_0)}{\lambda}\right)$",
                fontsize=16)
            plt.tight_layout()
            plt.savefig(os.path.join(FOLDER, "results", "problem3", "problem3b_zoom.pdf"))
            plt.show()
            print("Saved problem3b_zoom.pdf")


def problem3b_soliton_initial():
    """
    Run Boussinesq with a soliton initial condition and compare the output at
    t=800 to the analytical soliton solution. Shows smaller error than the
    cos^2 initial condition case.
    """
    from scipy.signal import find_peaks

    print("\n--- Problem 3b soliton IC: Boussinesq t=800 ---")
    alpha_ic = A          # use same amplitude as problem
    x0_ic    = 50.0       # start well inside domain

    out_dir = run_sbouss("Boussinesq", "b_soliton_ic_t800", dx=0.1, t_sim=800.0,
                         use_soliton=True, soliton_alpha=alpha_ic, soliton_x0=x0_ic)
    eta_file = os.path.join(out_dir, "eta1")
    if not os.path.exists(eta_file):
        print("Output not found.")
        return

    x, eta = read_eta(eta_file)

    peaks, _ = find_peaks(eta, height=0.01)
    peaks = peaks[np.argsort(x[peaks])[::-1]]

    if len(peaks) == 0:
        print("No peaks found.")
        return

    idx     = peaks[0]
    alpha_s = eta[idx]
    x0_s    = x[idx]

    k       = np.sqrt(3 * alpha_s / 4)
    x_sol   = np.linspace(x0_s - 50, x0_s + 50, 1000)
    eta_sol = alpha_s / np.cosh(k * (x_sol - x0_s)) ** 2
    eta_sol_x = alpha_s / np.cosh(k * (x - x0_s)) ** 2

    x_lo = int(x0_s) - 10
    x_hi = int(x0_s) + 10

    fig, (ax3, ax4) = plt.subplots(1, 2, figsize=(10, 5))

    ax3.plot(x, eta, color='#1A7ABF', linewidth=2, label="Boussinesq")
    ax3.plot(x_sol, eta_sol, color='#E07B00', linewidth=2, linestyle='--', label="Soliton solution")
    ax3.set_xlim(x_lo, x_hi)
    ax3.set_xticks(range(x_lo, x_hi + 1, 5))
    ax3.set_xlabel("x", fontsize=16)
    ax3.set_ylabel("$\\eta$", fontsize=16)
    ax3.legend(fontsize=12)
    ax3.grid(alpha=0.5)

    ax4.plot(x, np.abs(eta - eta_sol_x), color='#7B2D8B', linewidth=2)
    ax4.set_xlim(x_lo, x_hi)
    ax4.set_xticks(range(x_lo, x_hi + 1, 5))
    ax4.set_ylim(bottom=0)
    ax4.set_xlabel("x", fontsize=16)
    ax4.set_ylabel("$|\\eta - \\eta_{\\mathrm{sol}}|$", fontsize=16)
    ax4.grid(alpha=0.5)

    fig.suptitle(r"$\mathbf{Soliton\ crest\ vs\ soliton\ solution\ at\ t=800}$" + "\n" +
             r"IC: $\eta(x,0) = \alpha \, \mathrm{sech}^2\!\left(\left(\frac{3\alpha}{4}\right)^{1/2}(x-x_0)\right)$",
             fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(FOLDER, "results", "problem3", "problem3b_soliton_ic.pdf"))
    plt.show()
    print(f"Saved problem3b_soliton_ic.pdf  (alpha={alpha_s:.4f}, x0={x0_s:.1f})")

# =============================================================================
def problem3c():
    """
    c) Undular bore from shoaling.
    Solitary wave A=0.05 from deep (h=1) to shallow (h=0.2) region.
    """
    print("\n--- Problem 3c: Shoaling, undular bore ---")
    out_dir = run_sbouss(
        "Boussinesq", "c_shoaling",
        dx=0.1,           # fine grid needed for shallow shelf
        t_sim=400.0,
        flat=False,
        depth_file="depth.dat",
        use_soliton=True,
        soliton_alpha=0.05,
        soliton_x0=20.0,
    )

    x_depth, h = read_eta(os.path.join(DATA, "depth.dat"))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                                gridspec_kw={'height_ratios': [3, 1]})

    eta_file = os.path.join(out_dir, "eta1")
    if os.path.exists(eta_file):
        x, eta = read_eta(eta_file)
        ax1.plot(x, eta, linewidth=2, color='#1A7ABF')

    ax1.set_ylabel("$\\eta$", fontsize=16)
    ax1.set_title("Shoaling — undular bore", fontsize=18, fontweight='bold')
    ax1.tick_params(axis='both', labelsize=14)
    ax1.grid(True)

    ax2.fill_between(x_depth, -h, -1.1, color='#8B6914', alpha=0.6)
    ax2.plot(x_depth, -h, color='#8B6914', linewidth=1.5)
    ax2.set_xlabel("x", fontsize=16)
    ax2.tick_params(axis='x', labelsize=14)
    ax2.yaxis.set_visible(False)
    ax2.grid(True)

    # Depth arrows
    x_arrow_deep = 10
    x_arrow_shallow = 200
    h_deep = h[np.argmin(np.abs(x_depth - x_arrow_deep))]
    h_shallow = h[np.argmin(np.abs(x_depth - x_arrow_shallow))]

    ax2.annotate('', xy=(x_arrow_deep, 0),
                xytext=(x_arrow_deep, -h_deep),
                arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
    ax2.text(x_arrow_deep + 2, -h_deep/2, '$h=1$', fontsize=13, va='center')

    ax2.annotate('', xy=(x_arrow_shallow, 0),
                xytext=(x_arrow_shallow, -h_shallow),
                arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
    ax2.text(x_arrow_shallow + 2, -h_shallow/2, '$h=0.2$', fontsize=13, va='center')

    ax1.set_xlim(0, 300)
    plt.tight_layout()
    plt.savefig(os.path.join(FOLDER, "results", "problem3", "problem3c.pdf"))
    plt.show()
    print("Saved problem3c.pdf")


# =============================================================================
def problem3c_animation():
    """Animate the shoaling/undular bore from t=0 to t=500, saved as mp4."""
    from matplotlib.animation import FuncAnimation, PillowWriter

    print("\n--- Problem 3c animation ---")
    frame_times = np.arange(0, 505, 5)

    print("Running sbouss for each frame...")
    frames = []
    for t in frame_times:
        if t == 0:
            frames.append(None)
            continue
        out_dir = run_sbouss(
            "Boussinesq", f"c_anim_t{int(t):04d}",
            dx=0.1, t_sim=float(t),
            flat=False, depth_file="depth.dat",
            use_soliton=True, soliton_alpha=0.05, soliton_x0=20.0,
        )
        eta_file = os.path.join(out_dir, "eta1")
        if os.path.exists(eta_file):
            frames.append(read_eta(eta_file))
        else:
            frames.append(None)

    x_depth, h = read_eta(os.path.join(DATA, "depth.dat"))

    # Use first valid frame to get x grid
    x0_frame = next(f for f in frames if f is not None)[0]
    eta_init = np.zeros_like(x0_frame)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                                    gridspec_kw={'height_ratios': [3, 1]})

    line, = ax1.plot(x0_frame, eta_init, linewidth=2, color='#1A7ABF')
    time_text = ax1.text(0.02, 0.93, '', transform=ax1.transAxes, fontsize=14)

    ax1.set_xlim(0, 300)
    ax1.set_ylim(-0.05, 0.15)
    ax1.set_ylabel("$\\eta$", fontsize=16)
    ax1.set_title("Shoaling — undular bore", fontsize=18, fontweight='bold')
    ax1.tick_params(axis='both', labelsize=14)
    ax1.grid(True)

    ax2.fill_between(x_depth, -h, -1.1, color='#8B6914', alpha=0.6)
    ax2.plot(x_depth, -h, color='#8B6914', linewidth=1.5)
    ax2.set_xlim(0, 300)
    ax2.set_xlabel("x", fontsize=16)
    ax2.tick_params(axis='x', labelsize=14)
    ax2.yaxis.set_visible(False)
    ax2.grid(True)

    x_arrow_deep, x_arrow_shallow = 10, 200
    h_deep    = h[np.argmin(np.abs(x_depth - x_arrow_deep))]
    h_shallow = h[np.argmin(np.abs(x_depth - x_arrow_shallow))]
    ax2.annotate('', xy=(x_arrow_deep, 0), xytext=(x_arrow_deep, -h_deep),
                 arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
    ax2.text(x_arrow_deep + 2, -h_deep / 2, '$h=1$', fontsize=13, va='center')
    ax2.annotate('', xy=(x_arrow_shallow, 0), xytext=(x_arrow_shallow, -h_shallow),
                 arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
    ax2.text(x_arrow_shallow + 2, -h_shallow / 2, '$h=0.2$', fontsize=13, va='center')

    plt.tight_layout()

    def update(i):
        t = frame_times[i]
        if frames[i] is None:
            line.set_ydata(eta_init)
        else:
            _, eta = frames[i]
            line.set_ydata(eta)
        time_text.set_text(f'$t = {int(t)}$')
        return line, time_text

    ani = FuncAnimation(fig, update, frames=len(frame_times), interval=80, blit=True)

    out_path = os.path.join(FOLDER, "results", "problem3", "problem3c_animation.gif")
    writer = PillowWriter(fps=8)
    ani.save(out_path, writer=writer)
    plt.close(fig)
    print(f"Saved {out_path}")


# =============================================================================
if __name__ == "__main__":
    make_eta_in()
    make_u_in()
    make_depth_file()

    # problem3a()
    # problem3b()
    # problem3b_soliton_initial()
    # problem3c()
    # problem3c_animation()

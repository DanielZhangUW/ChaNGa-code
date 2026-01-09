import os, math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import rebound

# ===== constants / ICs =====
OMEGA = 1.0
T     = 2.0*math.pi
DT    = 0.01*T          # 100 steps per orbit
NSTEPS= 100

# Use BASE_X as the reference radius/amplitude; in the paper formula this is R
BASE_X = 1.0

# Default epicycle initial condition for comparisons (x0,y0,vx0,vy0)=(1,0,0,-2)
IC = dict(x=1.0, y=0.0, z=0.0, vx=0.0, vy=-2.0, vz=0.0)

# Visualization options (hodograph disabled for this task)
CENTER_HODOGRAPH = False
HODO_CENTER_Y    = 0.0


def run_sei_points_for_ic(ic_override, dt=DT, nsteps=NSTEPS):
    sim = rebound.Simulation()
    sim.integrator = "sei"
    sim.ri_sei.OMEGA = OMEGA
    sim.gravity = "none"; sim.boundary = "none"; sim.collision = "none"
    sim.add(**ic_override); sim.dt = dt
    X,Y = [], []
    p=sim.particles[0]; X.append(p.x); Y.append(p.y)
    t=0.0
    for _ in range(nsteps):
        t += dt
        sim.integrate(t)
        p=sim.particles[0]
        X.append(p.x); Y.append(p.y)
    return X,Y

def run_sei_points(dt=DT, nsteps=NSTEPS):
    return run_sei_points_for_ic(IC, dt=dt, nsteps=nsteps)


def vy_from_e(x, e, R, Omega, branch="-"):
    """Compute vy from target eccentricity e using
    e^2 = 1/(R^2 Ω^2) (P_x^2 − 4Ω x P_y + Ω^2 x^2 + 4 P_y^2), with P_x=v_x, P_y=v_y+2Ωx.

    We enforce v_x=0 and y=0 at t=0. Solve quadratic for P_y and choose branch.
    Branch "-" ensures e=1 gives the straight-down case vy = -2Ωx when R=x.
    """
    # Solve for P_y (q): 4 q^2 − 4 Ω x q + Ω^2(x^2 − e^2 R^2) = 0
    # Roots: q = (Ω/2)(x ± e R)
    if branch not in ("-", "+"):
        branch = "-"
    q = 0.5*Omega*(x - e*R) if branch == "-" else 0.5*Omega*(x + e*R)
    vy = q - 2.0*Omega*x
    return vy

def ic_from_e_fixed_guiding_center(e, R, Omega, x_gc=0.0):
    """Compute initial conditions (x0, vy) for a given eccentricity e 
    while keeping the guiding center fixed at x_gc.
    
    For X_gc = 0 and phase φ=0 (vx=0, y0=0):
    - x0 = a = e·R (amplitude/semi-major axis)
    - vy = -2Ω·a = -2Ω·x0
    
    Special case: When e=1, this gives x0=R, vy=-2ΩR, which produces
    a straight vertical line (the particle falls straight down).
    
    This ensures:
    1. Guiding center X_gc = x0 + vy/(2Ω) = x0 - x0 = 0 (fixed)
    2. All orbits share the same guiding center
    3. e represents the normalized amplitude
    
    General formula for arbitrary X_gc:
    - x0 = X_gc + a
    - vy = -2Ω·a - (3/2)Ω·X_gc = -2Ω·x0 + (1/2)Ω·X_gc
    """
    if abs(x_gc) < 1e-10:  # X_gc = 0 case
        x0 = e * R
        vy = -2.0 * Omega * x0
        # For e=1, this should give x0=R, vy=-2ΩR, which is the straight-down case
        # Verify: when e=1, x0=R, so vy=-2ΩR, which matches the paper definition
    else:  # General case for non-zero X_gc
        a = e * R
        x0 = x_gc + a
        vy = -2.0 * Omega * a - 1.5 * Omega * x_gc
    return x0, vy


def ic_for_transition_eccentricity(e, x_right, R=BASE_X, Omega=OMEGA):
    """Return IC, guiding center, and epicycle amplitude for fixed rightmost start."""
    a = e * R
    X_gc = x_right - a
    vy = -2.0 * Omega * a - 1.5 * Omega * X_gc
    ic = dict(x=x_right, y=0.0, z=0.0, vx=0.0, vy=vy, vz=0.0)
    return ic, X_gc, a


def ic_vertical_stream(x0, Omega=OMEGA):
    """Initial condition that produces a pure shear (vertical) trajectory at fixed x = x0."""
    vy = -1.5 * Omega * x0  # ensures dx/dt = 0, dy/dt = -1.5 Ω x0
    ic = dict(x=x0, y=0.0, z=0.0, vx=0.0, vy=vy, vz=0.0)
    return ic, x0  # guiding center equals x0 when the epicycle amplitude is zero


def run_sei_states_for_ic(ic_override, dt=DT, nsteps=NSTEPS):
    """Integrate with SEI and return (X, Y, VX, VY) over nsteps+1 outputs."""
    sim = rebound.Simulation()
    sim.integrator = "sei"
    sim.ri_sei.OMEGA = OMEGA
    sim.gravity = "none"; sim.boundary = "none"; sim.collision = "none"
    sim.add(**ic_override); sim.dt = dt
    X,Y,VX,VY = [], [], [], []
    p=sim.particles[0]
    X.append(p.x); Y.append(p.y); VX.append(p.vx); VY.append(p.vy)
    t=0.0
    for _ in range(nsteps):
        t += dt
        sim.integrate(t)
        p=sim.particles[0]
        X.append(p.x); Y.append(p.y); VX.append(p.vx); VY.append(p.vy)
    return X,Y,VX,VY


def resample_equal_arc(xs, ys, segments=10):
    """Resample a polyline (xs, ys) into (segments+1) points with equal arc length.
    Keeps endpoints. Assumes xs, ys are ordered along the trajectory.
    """
    if len(xs) < 2:
        return xs[:segments+1], ys[:segments+1]
    import bisect
    # cumulative arc length
    s = [0.0]
    for i in range(1, len(xs)):
        dx = xs[i] - xs[i-1]
        dy = ys[i] - ys[i-1]
        s.append(s[-1] + (dx*dx + dy*dy) ** 0.5)
    total = s[-1]
    if total == 0:
        return xs[:segments+1], ys[:segments+1]
    targets = [total * k / segments for k in range(segments+1)]
    xr, yr = [], []
    for t in targets:
        j = bisect.bisect_left(s, t)
        if j == 0:
            xr.append(xs[0]); yr.append(ys[0]); continue
        if j >= len(s):
            xr.append(xs[-1]); yr.append(ys[-1]); continue
        # linear interpolate between j-1 and j
        s0, s1 = s[j-1], s[j]
        w = 0.0 if s1 == s0 else (t - s0) / (s1 - s0)
        xr.append(xs[j-1] + w * (xs[j] - xs[j-1]))
        yr.append(ys[j-1] + w * (ys[j] - ys[j-1]))
    return xr, yr

def quinn_step(x,y,z,vx,vy,vz,dt,Omega):
    # Quinn KDK with conserved Py reuse
    vx_q14 = vx - 0.5*dt*(Omega*Omega*x)
    Py     = vy + 2.0*Omega*x
    vx_h   = vx_q14 + dt*Omega*Py
    vy_h   = Py - Omega*x - Omega*(x + dt*vx_h)
    vz_h   = vz + 0.5*dt*(-Omega*Omega*z)
    x1 = x + dt*vx_h;  y1 = y + dt*vy_h;  z1 = z + dt*vz_h
    vx_q34 = vx_h + dt*Omega*Py
    vx1    = vx_q34 - 0.5*dt*(Omega*Omega*x1)
    vy1    = Py - 2.0*Omega*x1
    vz1    = vz_h + 0.5*dt*(-Omega*Omega*z1)
    return x1,y1,z1,vx1,vy1,vz1

def run_quinn_points(dt=DT, nsteps=NSTEPS):
    x,y,z = IC["x"],IC["y"],IC["z"]
    vx,vy,vz = IC["vx"],IC["vy"],IC["vz"]
    X,Y = [x],[y]
    for _ in range(nsteps):
        x,y,z,vx,vy,vz = quinn_step(x,y,z,vx,vy,vz,dt,OMEGA)
        X.append(x); Y.append(y)
    return X,Y

def run_quinn_points_for_ic(ic_override, dt=DT, nsteps=NSTEPS):
    """Integrate with Quinn and return (X, Y) over nsteps+1 outputs."""
    x,y,z = ic_override["x"],ic_override["y"],ic_override["z"]
    vx,vy,vz = ic_override["vx"],ic_override["vy"],ic_override["vz"]
    print(
        f"Quinn IC (code units): x0={x:.6f}, y0={y:.6f}, "
        f"vx0={vx:.6f}, vy0={vy:.6f}"
    )
    X,Y = [x],[y]
    for _ in range(nsteps):
        x,y,z,vx,vy,vz = quinn_step(x,y,z,vx,vy,vz,dt,OMEGA)
        X.append(x); Y.append(y)
    return X,Y

def run_quinn_states_for_ic(ic_override, dt=DT, nsteps=NSTEPS):
    """Integrate with Quinn and return (X, Y, VX, VY) over nsteps+1 outputs."""
    x,y,z = ic_override["x"],ic_override["y"],ic_override["z"]
    vx,vy,vz = ic_override["vx"],ic_override["vy"],ic_override["vz"]
    X,Y,VX,VY = [x],[y],[vx],[vy]
    for _ in range(nsteps):
        x,y,z,vx,vy,vz = quinn_step(x,y,z,vx,vy,vz,dt,OMEGA)
        X.append(x); Y.append(y); VX.append(vx); VY.append(vy)
    return X, Y, VX, VY

# ===== Modified leapfrog integrator =====
def lf_modified_step(x,y,z,vx,vy,vz,dt,Omega):
    """Modified leapfrog: predictor-corrector, second-order, non-symplectic."""
    # Kick (half)
    vx_h = vx + 0.5*dt*( 3*(Omega**2)*x + 2*Omega*vy )
    vy_h = vy + 0.5*dt*( -2*Omega*vx )
    vz_h = vz + 0.5*dt*( -(Omega**2)*z )
    # Predictor
    vx_bar = vx + dt*( 3*(Omega**2)*x + 2*Omega*vy )
    vy_bar = vy + dt*( -2*Omega*vx )
    vz_bar = vz + dt*( -(Omega**2)*z )
    # Drift
    x1 = x + dt*vx_h
    y1 = y + dt*vy_h
    z1 = z + dt*vz_h
    # Kick (half) using predictor
    vx1 = vx_h + 0.5*dt*( 3*(Omega**2)*x1 + 2*Omega*vy_bar )
    vy1 = vy_h + 0.5*dt*( -2*Omega*vx_bar )
    vz1 = vz_h + 0.5*dt*( -(Omega**2)*z1 )
    return x1,y1,z1,vx1,vy1,vz1

def run_leapfrog_states_for_ic(ic_override, dt=DT, nsteps=NSTEPS):
    """Integrate with modified leapfrog and return (X, Y, VX, VY) over nsteps+1 outputs."""
    x,y,z = ic_override["x"],ic_override["y"],ic_override["z"]
    vx,vy,vz = ic_override["vx"],ic_override["vy"],ic_override["vz"]
    X,Y,VX,VY = [x],[y],[vx],[vy]
    for _ in range(nsteps):
        x,y,z,vx,vy,vz = lf_modified_step(x,y,z,vx,vy,vz,dt,OMEGA)
        X.append(x); Y.append(y); VX.append(vx); VY.append(vy)
    return X,Y,VX,VY


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(script_dir, "fig")
    os.makedirs(out_dir, exist_ok=True)

    plot_quinn_transition(out_dir, steps_per_orbit=80, suffix="", header_note="Δt = T/80 (≈0.0785)")
    plot_quinn_transition(out_dir, steps_per_orbit=10, suffix="_step10",
                          header_note="Δt = T/10 (10 steps/orbit)", orbits_closed=1.0, orbits_shear=3.0)
    plot_quinn_transition_blended(out_dir, steps_per_orbit=80, suffix="_blend",
                                  header_note="Δt = T/80 (X_gc blend)")
    plot_quinn_transition_blended(out_dir, steps_per_orbit=10, suffix="_blend_step10",
                                  header_note="Δt = T/10 (blend, 10 steps)", orbits_closed=1.0, orbits_shear=3.0)
    plot_quinn_a_sweep(out_dir, steps_per_orbit=80, suffix="_a_sweep",
                       header_note="Δt = T/80 (a sweep, e = 0)")
    plot_quinn_a_sweep(out_dir, steps_per_orbit=10, suffix="_a_sweep_step10",
                       header_note="Δt = T/10 (a sweep, e = 0, 10 steps)", orbits=3.0)


def plot_quinn_transition(out_dir, steps_per_orbit, suffix, header_note, orbits_closed=6.0, orbits_shear=6.0):
    """Generate Quinn-integrated transition plot for a given timestep."""
    plt.rcParams.update({
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "legend.fontsize": 10,
    })
    fig, ax = plt.subplots(1, 1, figsize=(18.0, 10.0))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x", fontweight="medium")
    ax.set_ylabel("y", fontweight="medium")
    ax.grid(True, alpha=0.18, linestyle="--", linewidth=0.6)

    e_levels = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
    elliptic_es = [e for e in e_levels if e < 1.0]
    R_epicycle = 0.3

    dt_quinn = T / float(steps_per_orbit)
    nsteps_closed = int(math.ceil(orbits_closed * T / dt_quinn))
    nsteps_shear = int(math.ceil(orbits_shear * T / dt_quinn))

    color_cycle = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b"]
    if len(color_cycle) < len(elliptic_es):
        raise ValueError("Not enough distinct colors for the chosen eccentricity levels.")
    color_vertical = "#d62728"

    print("\n" + "=" * 86)
    print(f"Quinn sweep ({header_note}): closed epicycles (X_gc = 0) → vertical shear (e = 1)")
    print(f"Ω = {OMEGA:.1f}, steps/orbit = {steps_per_orbit}, dt = {dt_quinn:.4f}, closed-steps = {nsteps_closed}, shear-steps = {nsteps_shear}")
    print("=" * 86)
    print(f"{'e':>5} {'a':>7} {'x0':>7} {'vy':>10} {'X_gc target':>14} {'X_gc true':>14}")
    print("-" * 86)

    for e in e_levels:
        if e < 1.0:
            x0, vy0 = ic_from_e_fixed_guiding_center(e=e, R=R_epicycle, Omega=OMEGA, x_gc=0.0)
            ic_now = dict(x=x0, y=0.0, z=0.0, vx=0.0, vy=vy0, vz=0.0)
            nsteps = nsteps_closed
            target_Xgc = 0.0
        else:
            ic_now, target_Xgc = ic_vertical_stream(R_epicycle, Omega=OMEGA)
            x0 = ic_now["x"]
            nsteps = nsteps_shear

        X_full, Y_full, _, _ = run_quinn_states_for_ic(ic_now, dt=dt_quinn, nsteps=nsteps)

        color = color_cycle[elliptic_es.index(e)] if e < 1.0 else color_vertical
        lw = 2.3 if abs(e - 1.0) < 1e-9 else 1.4
        line_kwargs = dict(
            linestyle="-",
            color=color,
            lw=lw,
            zorder=2,
            antialiased=False,
        )
        if steps_per_orbit <= 20:
            line_kwargs.update(marker="o", ms=4.0, markevery=1)
        ax.plot(X_full, Y_full, **line_kwargs)

        X_gc_true = (2.0 * ic_now["vy"]) / OMEGA + 4.0 * ic_now["x"]
        print(f"{e:5.1f} {e*R_epicycle:7.3f} {ic_now['x']:7.3f} {ic_now['vy']:10.6f} "
              f"{target_Xgc:14.6f} {X_gc_true:14.6f}")
        print(f"    stored points: {len(X_full)}")

        if e < 1.0:
            ax.scatter([X_full[0]], [Y_full[0]], marker="o", s=70,
                       facecolors="none", edgecolors=color, linewidths=1.4, zorder=5)
        else:
            ax.scatter([X_full[0]], [Y_full[0]], marker="x", s=80,
                       color=color, linewidths=1.8, zorder=5)

    print("-" * 86 + "\n")

    ax.set_title("Quinn integrator: epicycles (X_gc = 0) → vertical shear (e = 1)")
    ax.set_xlim(-25.0, 25.0)
    ax.set_ylim(-25.0, 25.0)

    legend_handles = [
        plt.Line2D([], [], color=color_cycle[i], lw=2.0, label=f"e = {elliptic_es[i]:.1f}")
        for i in range(len(elliptic_es))
    ]
    legend_handles.append(plt.Line2D([], [], color=color_vertical, lw=2.3, label="e = 1.0"))
    ax.legend(handles=legend_handles, loc="upper right", frameon=True, framealpha=0.92)

    plt.tight_layout()
    basename = f"e_sweep_quinn_only{suffix}"
    plt.savefig(os.path.join(out_dir, f"{basename}.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(out_dir, f"{basename}.pdf"), bbox_inches="tight")
    plt.show()


def plot_quinn_transition_blended(out_dir, steps_per_orbit, suffix, header_note,
                                  orbits_closed=6.0, orbits_shear=6.0, blend_exponent=3.0):
    """Same as plot_quinn_transition but smoothly shifts X_gc toward the shear solution as e→1."""
    plt.rcParams.update({
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "legend.fontsize": 10,
    })
    fig, ax = plt.subplots(1, 1, figsize=(18.0, 10.0))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x", fontweight="medium")
    ax.set_ylabel("y", fontweight="medium")
    ax.grid(True, alpha=0.18, linestyle="--", linewidth=0.6)

    e_levels = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
    elliptic_es = [e for e in e_levels if e < 1.0]
    R_epicycle = 0.3

    dt_quinn = T / float(steps_per_orbit)
    nsteps_closed = int(math.ceil(orbits_closed * T / dt_quinn))
    nsteps_shear = int(math.ceil(orbits_shear * T / dt_quinn))

    color_cycle = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b"]
    if len(color_cycle) < len(elliptic_es):
        raise ValueError("Not enough distinct colors for the chosen eccentricity levels.")
    color_vertical = "#d62728"

    print("\n" + "=" * 90)
    print(f"Quinn sweep blended ({header_note}): guided ellipses → shear (X_gc interpolation)")
    print(f"Ω = {OMEGA:.1f}, steps/orbit = {steps_per_orbit}, dt = {dt_quinn:.4f}, closed-steps = {nsteps_closed}, shear-steps = {nsteps_shear}")
    print("=" * 90)
    print(f"{'e':>5} {'blend':>8} {'x0':>7} {'vy':>10} {'X_gc target':>14} {'X_gc true':>14}")
    print("-" * 90)

    for e in e_levels:
        if e < 1.0:
            blend = min(1.0, max(0.0, e ** blend_exponent))
            x0_closed, vy_closed = ic_from_e_fixed_guiding_center(e=e, R=R_epicycle, Omega=OMEGA, x_gc=0.0)
            ic_vert, X_gc_vert = ic_vertical_stream(R_epicycle, Omega=OMEGA)
            x0_vert = ic_vert["x"]
            vy_vert = ic_vert["vy"]

            x0 = (1.0 - blend) * x0_closed + blend * x0_vert
            vy0 = (1.0 - blend) * vy_closed + blend * vy_vert
            target_Xgc = blend * X_gc_vert
            ic_now = dict(x=x0, y=0.0, z=0.0, vx=0.0, vy=vy0, vz=0.0)
            nsteps = nsteps_closed
        else:
            ic_now, target_Xgc = ic_vertical_stream(R_epicycle, Omega=OMEGA)
            x0 = ic_now["x"]
            vy0 = ic_now["vy"]
            blend = 1.0
            nsteps = nsteps_shear

        X_full, Y_full, _, _ = run_quinn_states_for_ic(ic_now, dt=dt_quinn, nsteps=nsteps)

        color = color_cycle[elliptic_es.index(e)] if e < 1.0 else color_vertical
        lw = 2.3 if abs(e - 1.0) < 1e-9 else 1.4
        line_kwargs = dict(linestyle="-", color=color, lw=lw, zorder=2, antialiased=False)
        if steps_per_orbit <= 20:
            line_kwargs.update(marker="o", ms=4.0, markevery=1)
        ax.plot(X_full, Y_full, **line_kwargs)

        X_gc_true = (2.0 * ic_now["vy"]) / OMEGA + 4.0 * ic_now["x"]
        print(f"{e:5.1f} {blend:8.3f} {ic_now['x']:7.3f} {ic_now['vy']:10.6f} "
              f"{target_Xgc:14.6f} {X_gc_true:14.6f}")
        print(f"    stored points: {len(X_full)}")

        if e < 1.0:
            ax.scatter([X_full[0]], [Y_full[0]], marker="o", s=70,
                       facecolors="none", edgecolors=color, linewidths=1.4, zorder=5)
        else:
            ax.scatter([X_full[0]], [Y_full[0]], marker="x", s=80,
                       color=color, linewidths=1.8, zorder=5)

    print("-" * 90 + "\n")

    ax.set_title("Quinn: blended X_gc → shear (e → 1)")
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-1.5, 1.5)

    legend_handles = [
        plt.Line2D([], [], color=color_cycle[i], lw=2.0, label=f"e = {elliptic_es[i]:.1f}")
        for i in range(len(elliptic_es))
    ]
    legend_handles.append(plt.Line2D([], [], color=color_vertical, lw=2.3, label="e = 1.0"))
    ax.legend(handles=legend_handles, loc="upper right", frameon=True, framealpha=0.92)

    plt.tight_layout()
    basename = f"e_sweep_quinn_only{suffix}"
    plt.savefig(os.path.join(out_dir, f"{basename}.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(out_dir, f"{basename}.pdf"), bbox_inches="tight")
    plt.show()


def plot_quinn_a_sweep(out_dir, steps_per_orbit, suffix, header_note,
                       a_levels=None, orbits=3.0, Omega=OMEGA):
    """Sweep semi-major axis (guiding centre) at zero eccentricity."""
    if a_levels is None:
        a_levels = [-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0]

    plt.rcParams.update({
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "legend.fontsize": 10,
    })
    fig, ax = plt.subplots(1, 1, figsize=(10.5, 7.0))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x", fontweight="medium")
    ax.set_ylabel("y", fontweight="medium")
    ax.grid(True, alpha=0.18, linestyle="--", linewidth=0.6)

    dt_quinn = T / float(steps_per_orbit)
    nsteps = int(math.ceil(orbits * T / dt_quinn))

    colors = plt.cm.tab10(np.linspace(0.05, 0.95, len(a_levels)))

    print("\n" + "=" * 82)
    print(f"Quinn a-sweep ({header_note}): e = 0, vary guiding centre (semi-major axis)")
    print(f"Ω = {Omega:.1f}, steps/orbit = {steps_per_orbit}, dt = {dt_quinn:.4f}, steps = {nsteps}")
    print("=" * 82)
    print(f"{'a (X_gc)':>10} {'x0':>10} {'vy':>12} {'X_gc check':>14}")
    print("-" * 82)

    for a_val, color in zip(a_levels, colors):
        ic_now, X_gc_target = ic_vertical_stream(a_val, Omega=Omega)
        X_full, Y_full, _, _ = run_quinn_states_for_ic(ic_now, dt=dt_quinn, nsteps=nsteps)

        line_kwargs = dict(linestyle="-", color=color, lw=2.0, zorder=2, antialiased=False)
        if steps_per_orbit <= 20:
            line_kwargs.update(marker="o", ms=4.0, markevery=1)
        ax.plot(X_full, Y_full, **line_kwargs)

        X_gc_check = ic_now["x"] + 2.0 * ic_now["vy"] / (3.0 * Omega)
        print(f"{a_val:10.3f} {ic_now['x']:10.3f} {ic_now['vy']:12.6f} {X_gc_check:14.6f}")
        print(f"    stored points: {len(X_full)}")

        ax.scatter([X_full[0]], [Y_full[0]], marker="x", s=70,
                   color=color, linewidths=1.6, zorder=5)

    print("-" * 82 + "\n")

    ax.set_title("Quinn integrator: e = 0, semi-major axis sweep")
    a_abs_max = max(abs(a) for a in a_levels) if a_levels else 0.1
    x_limit = 10.0
    ax.set_xlim(-25.0, 25.0)
    ax.set_xticks(np.arange(-25, 26, 5))
    ax.set_ylim(-25.0, 25.0)

    legend_handles = [
        plt.Line2D([], [], color=colors[i], lw=2.0, label=f"a = {a_levels[i]:.2f}")
        for i in range(len(a_levels))
    ]
    ax.legend(handles=legend_handles, loc="upper right", frameon=True, framealpha=0.92)

    plt.tight_layout()
    basename = f"quinn_a_sweep{suffix}"
    plt.savefig(os.path.join(out_dir, f"{basename}.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(out_dir, f"{basename}.pdf"), bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()

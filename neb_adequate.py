from __future__ import annotations

import os
import subprocess
import sys
import time

from lammps import lammps

# ---------------------------------------------------------------------------
# Paths and Constants
# ---------------------------------------------------------------------------

DIR = os.path.dirname(os.path.abspath(__file__))

# Create a dedicated directory for all temporary/output files
TMP_DIR = os.path.join(DIR, "tmp_neb_files")
os.makedirs(TMP_DIR, exist_ok=True)

DATA_M = os.path.join(DIR, "GaN_ortho_rotated.lmp")
DATA_A = os.path.join(DIR, "GaN_ortho_rotated.lmp")
LIB = os.path.join(DIR, "library.meam")
MEAM = os.path.join(DIR, "GaN.meam")

ZVAC = 10.0
FROZ = 2.0
ADATOM_BOND_DIST = 2.2  # Angstroms above the surface to safely place the adatom

# Surface-N anchor (x, y)
M_N = (1.5945, 8.000)
A_N = (1.5945, 54.310)

# Unit-cell steps (x, y)
UC_1120 = (3.189, 0.0)
UC_0001 = (0.0, 5.185)

REF = {
    "mplane_11-20": 0.21,
    "mplane_0001": 0.93,
    "aplane_0001": 0.32,
    "aplane_1-100": 0.63,
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _lmp():
    return lammps()


def _setup(lmp, data):
    lmp.command("clear")
    lmp.command("units metal")
    lmp.command("dimension 3")
    lmp.command("atom_style atomic")
    lmp.command("boundary p p f")

    # Use forward slashes for LAMMPS paths to prevent escape character issues
    lmp.command(f"read_data {data.replace(os.sep, '/')}")
    lmp.command("change_box all z delta 0 %.1f" % ZVAC)
    lmp.command("mass 1 14.007")
    lmp.command("mass 2 69.723")
    lmp.command("pair_style meam")
    lmp.command(
        f"pair_coeff * * {LIB.replace(os.sep, '/')} Ga N {MEAM.replace(os.sep, '/')} Ga N"
    )
    lmp.command("neighbor 2.0 bin")
    lmp.command("neigh_modify delay 0 every 1 check yes")

    # Groups
    lmp.command("region freeze block INF INF INF INF INF %.1f" % FROZ)
    lmp.command("group frozen region freeze")
    lmp.command("group substrate type 1 2")
    lmp.command("group unfrozen subtract substrate frozen")
    lmp.command("fix freeze frozen setforce 0.0 0.0 0.0")
    lmp.command("thermo 0")
    lmp.command("thermo_style custom step pe")


def _get_max_z(lmp):
    """Dynamically finds the highest Z coordinate of the existing substrate."""
    lmp.command("compute maxz substrate property/atom z")
    lmp.command("compute zmax substrate reduce max c_maxz")
    lmp.command("run 0")
    return lmp.extract_compute("zmax", 0, 0)


# ---------------------------------------------------------------------------
# Single-image relaxation (Endpoints)
# ---------------------------------------------------------------------------


def _relax(data, xy, suffix=""):
    """Return (PE, datafile_path, adatom_id, relaxed_xyz) after CG minimisation."""
    lmp = _lmp()
    _setup(lmp, data)

    # 1. Dynamically find the surface height
    z_max = _get_max_z(lmp)

    # 2. Wrap position into periodic box
    lx = lmp.get_thermo("lx")
    ly = lmp.get_thermo("ly")

    wx = ((xy[0] % lx) + lx) % lx
    wy = ((xy[1] % ly) + ly) % ly
    wz = z_max + ADATOM_BOND_DIST  # Safely place exactly on top of the surface

    # 3. Create the single adatom
    lmp.command("create_atoms 2 single %.10f %.10f %.4f units box" % (wx, wy, wz))
    ad_id = lmp.get_natoms()

    # 4. Minimize
    lmp.command("group nebatoms id %d" % ad_id)
    lmp.command("min_style cg")
    lmp.command("minimize 1.0e-8 1.0e-10 10000 10000")
    lmp.command("run 0")
    pe = lmp.get_thermo("pe")

    # 5. Extract relaxed adatom coordinates
    coords = lmp.gather_atoms("x", 1, 3)
    i = (ad_id - 1) * 3
    relaxed_xyz = (coords[i], coords[i + 1], coords[i + 2])

    # Save the output to the temporary folder
    out = os.path.join(TMP_DIR, "relax%s.data" % suffix)
    lmp.command(f"write_data {out.replace(os.sep, '/')}")
    return pe, out, ad_id, relaxed_xyz


# ---------------------------------------------------------------------------
# NEB runner
# ---------------------------------------------------------------------------


def run_neb(label, data, init_xy, final_xy, n_img=12):
    t0 = time.time()
    print("\n%s" % "=" * 60)
    print("NEB  %s" % label)
    print("      init  (x=%.4f, y=%.4f)" % init_xy)
    print("      final (x=%.4f, y=%.4f)" % final_xy)
    print("%s" % "=" * 60)

    # --- endpoints --------------------------------------------------------
    print("  initial ...")
    pe0, init_data, _, _ = _relax(data, init_xy, "_init")
    print("          E = %.6f eV" % pe0)

    print("  final   ...")
    pe1, final_data, ad_id, relaxed_xyz = _relax(data, final_xy, "_final")
    print("          E = %.6f eV   naive dE = %+.6f eV" % (pe1, pe1 - pe0))

    # --- write NEB final coordinate file ----------------------------------
    final_neb = os.path.join(TMP_DIR, "final_%s.neb" % label)
    print(
        "  final adatom id=%d  relaxed coords=(%.6f, %.6f, %.6f)"
        % (ad_id, relaxed_xyz[0], relaxed_xyz[1], relaxed_xyz[2])
    )
    with open(final_neb, "w") as f:
        f.write("1\n")
        f.write(
            "%d %.10f %.10f %.4f\n"
            % (ad_id, relaxed_xyz[0], relaxed_xyz[1], relaxed_xyz[2])
        )

    # --- fix neb via subprocess -------------------------------------------
    print("  fix neb (%d images) ..." % n_img)
    barrier = None
    neb_in = os.path.join(TMP_DIR, "in.%s.neb" % label)
    dump_path = os.path.join(TMP_DIR, f"dump_{label}.*.xyz")

    with open(neb_in, "w") as f:
        f.write("clear\nunits metal\ndimension 3\natom_style atomic\n")
        f.write("boundary p p f\n")
        f.write("atom_modify map array sort 0 0.0\n")
        f.write("read_data %s\n" % init_data.replace(os.sep, "/"))
        f.write("change_box all z delta 0 %.1f\n" % ZVAC)
        f.write("mass 1 14.007\nmass 2 69.723\n")
        f.write("pair_style meam\n")
        f.write(
            "pair_coeff * * %s Ga N %s Ga N\n"
            % (LIB.replace(os.sep, "/"), MEAM.replace(os.sep, "/"))
        )
        f.write("neighbor 2.0 bin\nneigh_modify delay 0 every 1 check yes\n")
        f.write("region freeze block INF INF INF INF INF %.1f\n" % FROZ)
        f.write("group frozen region freeze\n")
        f.write("group substrate type 1 2\n")
        f.write("group unfrozen subtract substrate frozen\n")
        f.write("variable _nad equal count(all)\n")
        f.write("group nebatoms id ${_nad}\n")
        f.write("group nonneb subtract all nebatoms\n")
        f.write("fix freeze frozen setforce 0.0 0.0 0.0\n")
        f.write("fix freeze_nonneb nonneb setforce 0.0 0.0 0.0\n")
        f.write("fix neb1 nebatoms neb 1.0 parallel ideal\n")
        f.write("timestep 0.0005\n")
        f.write("min_style fire\n")
        f.write("thermo 50\n")
        f.write("thermo_style custom step temp pe ke etotal press fmax fnorm\n")
        f.write(
            "dump dmp all custom 50 %s id type xu yu zu\n"
            % dump_path.replace(os.sep, "/")
        )
        f.write("dump_modify dmp element N Ga\n")
        f.write("neb 0.0 0.1 1000 500 1 final %s\n" % final_neb.replace(os.sep, "/"))

    try:
        r = subprocess.run(
            [
                "mpirun",
                "--oversubscribe",
                "--map-by",
                "slot",
                "--bind-to",
                "none",
                "-np",
                str(n_img),
                "/opt/homebrew/bin/lmp_mpi",
                "-partition",
                "%dx1" % n_img,
                "-in",
                neb_in,
            ],
            capture_output=True,
            text=True,
            timeout=600,
            cwd=DIR,
        )
        if r.returncode != 0:
            print("  fix neb LAMMPS error (exit code %d):" % r.returncode)
            for line in r.stderr.strip().splitlines()[-15:]:
                print("    | " + line)
        else:
            pe_max = None
            for line in r.stdout.splitlines():
                p = line.split()
                if len(p) >= 11:
                    try:
                        for i in range(10, len(p), 2):
                            pe = float(p[i])
                            if pe_max is None or pe > pe_max:
                                pe_max = pe
                    except (ValueError, IndexError):
                        pass
            if pe_max is not None:
                barrier = pe_max - pe0
                print("  fix neb barrier = %.6f eV" % barrier)
            else:
                print("  fix neb: could not parse barrier from output")
    except Exception as e:
        print("  fix neb failed: %s" % e)
    finally:
        if os.path.exists(neb_in):
            os.remove(neb_in)

    # --- cleanup ----------------------------------------------------------
    # Removes the initial, final, and neb coordinate files to keep the tmp directory strictly for dumps
    for p in (init_data, final_data, final_neb):
        if os.path.exists(p):
            os.remove(p)

    elapsed = time.time() - t0
    print("  elapsed %.1f s\n" % elapsed)
    return pe0, pe1, barrier, elapsed


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    print("=" * 60)
    print("GaN adatom diffusion NEB")
    print("Lymperakis & Neugebauer PRB 79 241308 (2009)")
    print("Temporary files stored in: %s" % TMP_DIR)
    print("=" * 60)

    # 2D coordinates only (x, y). The Z coordinate is calculated dynamically.
    paths = [
        (
            "mplane_11-20",
            DATA_M,
            (M_N[0], M_N[1]),
            (M_N[0] + UC_1120[0], M_N[1] + UC_1120[1]),
        ),
        (
            "mplane_0001",
            DATA_M,
            (M_N[0], M_N[1]),
            (M_N[0] + UC_1120[0] + UC_0001[0], M_N[1] + UC_0001[1]),
        ),
        (
            "aplane_0001",
            DATA_A,
            (A_N[0], A_N[1]),
            (A_N[0] + UC_0001[0], A_N[1] + UC_0001[1]),
        ),
        ("aplane_1-100", DATA_A, (A_N[0], A_N[1]), (A_N[0] + UC_1120[0], A_N[1])),
    ]

    results = []
    for label, data, init, final in paths:
        pe0, pe1, bar, dt = run_neb(label, data, init, final)
        results.append((label, pe0, pe1, bar, dt))

    # Report
    print("\n" + "=" * 70)
    print(
        "%-20s  %8s  %8s  %8s  %8s"
        % ("Path", "E0(eV)", "Ebar(eV)", "DFT(eV)", "Time(s)")
    )
    print("-" * 70)
    for label, pe0, pe1, bar, dt in results:
        bar_str = f"{bar:8.4f}" if bar is not None else "   ERROR"
        print("%-20s  %8.4f  %8s  %8.3f  %8.1f" % (label, pe0, bar_str, REF[label], dt))
    print("=" * 70)


if __name__ == "__main__":
    main()

"""
GaN adatom diffusion NEB pipeline
Lymperakis & Neugebauer, PRB 79, 241308 (2009)

Flow per path
-------------
  1. lammps: read_data slab → create_atoms 2 single <init> → minimize cg
             → write_data <label>_init.data
  2. lammps: read_data slab → create_atoms 2 single <final> → minimize cg
             → write_data <label>_final.data
  3. Write in.<label>.neb  (fix neb + neb command)
  4. subprocess: mpirun -np 12 lmp -partition 12x1 -in in.<label>.neb
  5. Parse barrier, remove temp files

Surfaces
--------
  m-plane : GaN_ortho_rotated.lmp
    x = [2-1-10]  y = [0001]  z = [-12-10] (normal)
    adatom z = 29.5 A (vacuum)
  a-plane : GaN_mplane.lmp
    x = [2-1-10]  y = [0001]  z = [1-100] (normal)
    adatom y = 62.0 A (vacuum)

Reference barriers (DFT)
------------------------
  m-plane [11-20] : 0.21 eV   m-plane [0001] : 0.93 eV
  a-plane [0001]  : 0.32 eV   a-plane [1-100]: 0.63 eV
"""

from __future__ import annotations

import os
import subprocess
import sys
import time

from lammps import lammps

# ---------------------------------------------------------------------------
# paths and constants
# ---------------------------------------------------------------------------

DIR = os.path.dirname(os.path.abspath(__file__))

DATA_M = os.path.join(DIR, "GaN_ortho_rotated.lmp")
DATA_A = os.path.join(DIR, "GaN_mplane.lmp")
LIB = os.path.join(DIR, "library.meam")
MEAM = os.path.join(DIR, "GaN.meam")

ZVAC = 10.0
FROZ = 2.0
T_AD = 3

# surface-N anchor (x, y) and adatom vacuum coordinate
M_N = (1.5945, 8.000)
A_N = (1.5945, 54.310)
Z_M = 29.5
Z_A = 30.0

# unit-cell steps (x, y)
UC_1120 = (3.189, 0.0)
UC_0001 = (0.0, 5.185)

REF = {
    "mplane_11-20": 0.21,
    "mplane_0001": 0.93,
    "aplane_0001": 0.32,
    "aplane_1-100": 0.63,
}

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _lmp():
    return lammps()


def _setup(lmp, data):
    lmp.command("clear")
    lmp.command("units metal")

    lmp.command("dimension 3")
    lmp.command("atom_style atomic")
    lmp.command("boundary p p f")
    lmp.command("read_data " + data)
    lmp.command("change_box all z delta 0 %.1f" % ZVAC)
    lmp.command("mass 1 14.007")
    lmp.command("mass 2 69.723")
    lmp.command("pair_style meam")
    lmp.command("pair_coeff * * %s Ga N %s Ga N" % (LIB, MEAM))
    lmp.command("neighbor 2.0 bin")
    lmp.command("neigh_modify delay 0 every 1 check yes")
    lmp.command("region freeze block INF INF INF INF INF %.1f" % FROZ)
    lmp.command("group frozen region freeze")
    lmp.command("group substrate type 1 2")
    lmp.command("group unfrozen subtract substrate frozen")
    lmp.command("fix freeze frozen setforce 0.0 0.0 0.0")
    lmp.command("thermo 0")
    lmp.command("thermo_style custom step pe")


# ---------------------------------------------------------------------------
# single-image relaxation  (used for endpoints and fallback)
# ---------------------------------------------------------------------------


def _relax(data, xyz, suffix=""):
    """Return (PE, datafile_path, adatom_id, relaxed_xyz) after CG minimisation."""
    lmp = _lmp()
    _setup(lmp, data)
    # wrap position into periodic box to avoid create_atoms silently failing
    lx = lmp.get_thermo("lx")
    ly = lmp.get_thermo("ly")
    lz = lmp.get_thermo("lz")
    wx = ((xyz[0] % lx) + lx) % lx
    wy = ((xyz[1] % ly) + ly) % ly
    wz = ((xyz[2] % lz) + lz) % lz
    lmp.command("create_atoms 2 single %.10f %.10f %.4f units box" % (wx, wy, wz))
    ad_id = lmp.get_natoms()
    lmp.command("group nebatoms id %d" % ad_id)
    lmp.command("min_style cg")
    lmp.command("minimize 1.0e-8 1.0e-10 10000 10000")
    lmp.command("run 0")
    pe = lmp.get_thermo("pe")
    # extract relaxed adatom coordinates
    coords = lmp.gather_atoms("x", 1, 3)
    i = (ad_id - 1) * 3
    relaxed_xyz = (coords[i], coords[i + 1], coords[i + 2])
    out = os.path.join(DIR, "_tmp_relax%s.data" % suffix)
    lmp.command("write_data " + out)
    return pe, out, ad_id, relaxed_xyz


# ---------------------------------------------------------------------------
# NEB runner
# ---------------------------------------------------------------------------


def run_neb(label, data, init_xyz, final_xyz, n_img=12):
    t0 = time.time()
    print("\n%s" % "=" * 60)
    print("NEB  %s" % label)
    print("      init  (%.4f %.4f %.4f)" % init_xyz)
    print("      final (%.4f %.4f %.4f)" % final_xyz)
    print("%s" % "=" * 60)

    # --- endpoints --------------------------------------------------------
    print("  initial …")
    pe0, init_data, _, _ = _relax(data, init_xyz, "_init")
    print("          E = %.6f eV" % pe0)

    print("  final   …")
    pe1, final_data, ad_id, relaxed_xyz = _relax(data, final_xyz, "_final")
    print("          E = %.6f eV   naive dE = %+.6f eV" % (pe1, pe1 - pe0))

    # --- write NEB final coordinate file (simple format) ------------------
    final_neb = os.path.join(DIR, "_final_%s.neb" % label)
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

    # --- try true fix neb via subprocess ----------------------------------
    print("  fix neb (%d images) ..." % n_img)
    barrier = None
    neb_in = os.path.join(DIR, "in.%s.neb" % label)
    with open(neb_in, "w") as f:
        f.write("clear\nunits metal\ndimension 3\natom_style atomic\n")
        f.write("boundary p p f\n")
        f.write("atom_modify map array sort 0 0.0\n")
        f.write("read_data %s\n" % init_data)
        f.write("change_box all z delta 0 %.1f\n" % ZVAC)
        f.write("mass 1 14.007\nmass 2 69.723\n")
        f.write("pair_style meam\n")
        f.write("pair_coeff * * %s Ga N %s Ga N\n" % (LIB, MEAM))
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
        f.write("dump dmp all custom 50 dump_%s.*.xyz id type xu yu zu\n" % label)
        f.write("dump_modify dmp element N Ga\n")
        f.write("neb 0.0 0.1 1000 500 1 final %s\n" % final_neb)

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
            stderr_lines = r.stderr.strip().splitlines()
            for line in stderr_lines[-15:]:
                print("    | " + line)
            stdout_lines = r.stdout.strip().splitlines()
            if stdout_lines:
                print("  last stdout:")
                for line in stdout_lines[-10:]:
                    print("    | " + line)
        else:
            pe_max = None
            for line in r.stdout.splitlines():
                p = line.split()
                # NEB output: 9 header cols + RD[0] PE[0] RD[1] PE[1] ... (2 per replica)
                # Total = 9 + 2*N (odd). PE values start at index 10, step 2.
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
                print("  last stdout lines:")
                for line in r.stdout.strip().splitlines()[-15:]:
                    print("    | " + line)
    except FileNotFoundError:
        print("  ERROR: 'lmp' executable not found in PATH")
        print("  Install LAMMPS or add it to your PATH")
    except subprocess.TimeoutExpired:
        print("  fix neb timed out after 600 s")
    except Exception as e:
        print("  fix neb failed: %s" % e)
    finally:
        if os.path.exists(neb_in):
            os.remove(neb_in)

    if barrier is None:
        sys.exit(
            "NEB failed (see messages above). "
            "Run manually: mpirun -np %d lmp -partition %dx1 -in in.%s.neb"
            % (n_img, n_img, label)
        )

    # --- cleanup ----------------------------------------------------------
    for p in (init_data, final_data):
        if os.path.exists(p):
            os.remove(p)

    elapsed = time.time() - t0
    print("  elapsed %.1f s\n" % elapsed)
    return pe0, pe1, barrier, elapsed


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main():
    print("=" * 60)
    print("GaN adatom diffusion NEB")
    print("Lymperakis & Neugebauer PRB 79 241308 (2009)")
    print("=" * 60)

    paths = [
        (
            "mplane_11-20",
            DATA_M,
            (M_N[0], M_N[1], Z_M),
            (M_N[0] + UC_1120[0], M_N[1] + UC_1120[1], Z_M),
        ),
        (
            "mplane_0001",
            DATA_M,
            (M_N[0], M_N[1], Z_M),
            (M_N[0] + UC_1120[0] + UC_0001[0], M_N[1] + UC_0001[1], Z_M),
        ),
        (
            "aplane_0001",
            DATA_A,
            (A_N[0], A_N[1], Z_A),
            (A_N[0] + UC_0001[0], A_N[1] + UC_0001[1], Z_A),
        ),
        (
            "aplane_1-100",
            DATA_A,
            (A_N[0], A_N[1], Z_A),
            (A_N[0] + UC_1120[0], A_N[1], Z_A),
        ),
    ]

    results = []
    for label, data, init, final in paths:
        pe0, pe1, bar, dt = run_neb(label, data, init, final)
        results.append((label, pe0, pe1, bar, dt))

    # report
    print("\n" + "=" * 70)
    print(
        "%-20s  %8s  %8s  %8s  %8s"
        % ("Path", "E0(eV)", "Ebar(eV)", "DFT(eV)", "Time(s)")
    )
    print("-" * 70)
    for label, pe0, pe1, bar, dt in results:
        print("%-20s  %8.4f  %8.4f  %8.3f  %8.1f" % (label, pe0, bar, REF[label], dt))
    print("=" * 70)
    print("""
To run true multi-replica fix neb (requires MPI):
  mpirun -np 12 lmp -partition 12x1 -in in.<label>.neb

Reference (DFT, PRB 79 241308 2009):
  M-plane [0001]  : 0.93 eV   [11-20] : 0.21 eV
  A-plane [0001]  : 0.32 eV   [1-100] : 0.63 eV
""")


if __name__ == "__main__":
    main()

import os
import subprocess
import time
from lammps import lammps

# ---------------------------------------------------------------------------
# Paths and Constants
# ---------------------------------------------------------------------------

DIR = os.path.dirname(os.path.abspath(__file__))
TMP_DIR = os.path.join(DIR, "tmp_neb_files")
os.makedirs(TMP_DIR, exist_ok=True)

# 1-100 a-plane specific files and parameters
DATA_A = os.path.join(DIR, "GaN_ortho_rotated.lmp")
LIB = os.path.join(DIR, "library.meam")
MEAM = os.path.join(DIR, "GaN.meam")

ZVAC = 10.0
FROZ = 2.0
# Increased slightly to prevent MEAM overlap explosion before minimization
ADATOM_BOND_DIST = 2.5

# Start and End coordinates for aplane_1-100
START_XY = (1.5945, 54.310)
END_XY = (1.5945 + 3.189, 54.310)

# Number of NEB images (Keep this matching your MPI processors)
N_IMAGES = 8

# ---------------------------------------------------------------------------
# Endpoint Relaxation (Initial & Final)
# ---------------------------------------------------------------------------


def _relax(data, xy, suffix=""):
    """
    Relaxes the ENTIRE unfrozen system (substrate + adatom).
    Returns PE, path to data file, number of atoms, and relaxed coordinates.
    """
    lmp = lammps()

    # 1. Initialization
    lmp.command("clear")
    lmp.command("units metal")
    lmp.command("dimension 3")
    lmp.command("atom_style atomic")
    lmp.command("atom_modify map array sort 0 0.0")
    lmp.command("boundary p p f")

    # 2. System Definition
    lmp.command(f"read_data {data.replace(os.sep, '/')}")
    lmp.command(f"change_box all z delta 0 {ZVAC}")
    lmp.command("mass 1 14.007")
    lmp.command("mass 2 69.723")

    # 3. Settings
    lmp.command("pair_style meam")
    lmp.command(
        f"pair_coeff * * {LIB.replace(os.sep, '/')} Ga N {MEAM.replace(os.sep, '/')} Ga N"
    )
    lmp.command("neighbor 2.0 bin")
    lmp.command("neigh_modify delay 0 every 1 check yes")

    # 4. Atom Creation (Dynamically place the adatom)
    lmp.command("compute maxz all property/atom z")
    lmp.command("compute zmax all reduce max c_maxz")
    lmp.command("run 0")
    z_max = lmp.extract_compute("zmax", 0, 0)

    lx, ly = lmp.get_thermo("lx"), lmp.get_thermo("ly")
    wx = ((xy[0] % lx) + lx) % lx
    wy = ((xy[1] % ly) + ly) % ly
    wz = z_max + ADATOM_BOND_DIST

    lmp.command(f"create_atoms 2 single {wx} {wy} {wz} units box")

    # 5. Define Groups & Fixes
    lmp.command(f"region freeze block INF INF INF INF INF {FROZ}")
    lmp.command("group frozen region freeze")
    lmp.command("group unfrozen subtract all frozen")

    # Freeze the bottom layer so it anchors the crystal
    lmp.command("fix f_freeze frozen setforce 0.0 0.0 0.0")

    # 6. Run (Minimize the whole unfrozen system using FIRE for stability)
    lmp.command("min_style fire")
    lmp.command("minimize 1.0e-8 1.0e-10 10000 10000")
    lmp.command("run 0")

    pe = lmp.get_thermo("pe")
    out_data = os.path.join(TMP_DIR, f"relax_{suffix}.data").replace(os.sep, "/")
    lmp.command(f"write_data {out_data}")

    # Gather relaxed coordinates of ALL atoms to write to final.neb
    natoms = lmp.get_natoms()
    x = lmp.gather_atoms("x", 1, 3)
    coords = [x[i] for i in range(natoms * 3)]

    return pe, out_data, natoms, coords


# ---------------------------------------------------------------------------
# NEB Execution
# ---------------------------------------------------------------------------


def run_simple_neb():
    t0 = time.time()
    print("=" * 60)
    print("Running NEB for aplane_1-100")
    print("=" * 60)

    # --- 1. Relax Endpoints ---
    print("Minimizing initial state...")
    pe0, init_data, natoms_i, _ = _relax(DATA_A, START_XY, "init")
    print(f"  E_init = {pe0:.6f} eV")

    print("Minimizing final state...")
    pe1, _, natoms_f, final_coords = _relax(DATA_A, END_XY, "final")
    print(f"  E_final = {pe1:.6f} eV   (dE = {pe1 - pe0:+.6f} eV)")

    # --- 2. Write final.neb using ALL atoms ---
    # Because the unfrozen surface relaxes differently in the initial vs final state,
    # we must supply the coordinates of ALL atoms so LAMMPS interpolates them smoothly.
    final_neb = os.path.join(TMP_DIR, "final_aplane.neb")
    with open(final_neb, "w") as f:
        f.write(f"{natoms_f}\n")
        for i in range(natoms_f):
            x, y, z = (
                final_coords[i * 3],
                final_coords[i * 3 + 1],
                final_coords[i * 3 + 2],
            )
            # LAMMPS IDs are 1-indexed
            f.write(f"{i + 1} {x:.10f} {y:.10f} {z:.10f}\n")

    # --- 3. Generate NEB Input Script ---
    neb_in = os.path.join(TMP_DIR, "in.neb")
    dump_path = os.path.join(TMP_DIR, "dump_aplane.*.xyz").replace(os.sep, "/")

    with open(neb_in, "w") as f:
        f.write("clear\nunits metal\ndimension 3\n")
        f.write("atom_style atomic\n")
        f.write("atom_modify map array sort 0 0.0\n")
        f.write("boundary p p f\n")

        # Load the fully relaxed initial state
        f.write(f"read_data {init_data}\n")

        f.write("mass 1 14.007\nmass 2 69.723\n")
        f.write("pair_style meam\n")
        f.write(
            f"pair_coeff * * {LIB.replace(os.sep, '/')} Ga N {MEAM.replace(os.sep, '/')} Ga N\n"
        )
        f.write("neighbor 2.0 bin\nneigh_modify delay 0 every 1 check yes\n")

        f.write(f"region freeze block INF INF INF INF INF {FROZ}\n")
        f.write("group frozen region freeze\n")

        # CRITICAL FIX: The entire unfrozen surface must be in the NEB group
        # This prevents it from being locked in place and exploding.
        f.write("group nebatoms subtract all frozen\n")
        f.write("fix freeze frozen setforce 0.0 0.0 0.0\n")

        f.write("fix neb1 nebatoms neb 1.0 parallel ideal\n")
        f.write("timestep 0.001\n")
        f.write("min_style fire\n")
        f.write("thermo 50\n")
        f.write("thermo_style custom step temp pe etotal press fmax fnorm\n")
        f.write(f"dump dmp all custom 50 {dump_path} id type xu yu zu\n")
        f.write("dump_modify dmp element N Ga\n")

        f.write(f"neb 0.0 0.1 1000 500 1 final {final_neb.replace(os.sep, '/')}\n")

    # --- 4. Run NEB via MPI ---
    print(f"Running NEB with {N_IMAGES} images...")
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
                str(N_IMAGES),
                "/opt/homebrew/bin/lmp_mpi",  # <--- UPDATE THIS IF NECESSARY
                "-partition",
                f"{N_IMAGES}x1",
                "-in",
                neb_in,
            ],
            capture_output=True,
            text=True,
            timeout=600,
            cwd=DIR,
        )
        if r.returncode != 0:
            print(f"LAMMPS crashed (Exit code {r.returncode}). Tail of stderr:")
            for line in r.stderr.strip().splitlines()[-15:]:
                print("  | " + line)
        else:
            pe_values = []

            for line in r.stdout.splitlines():
                tokens = line.split()
                # Ensure line isn't a header/log message and contains enough tokens
                if len(tokens) >= 11 and "fire" not in line.lower():
                    try:
                        # Attempt numeric conversion of potential energy at index 10
                        val = float(tokens[10])
                        pe_values.append(val)
                    except ValueError:
                        # Skip thermo headers or non-numeric standard output lines
                        continue

            # Ensure we captured enough image data
            if len(pe_values) >= N_IMAGES:
                # Slice the LAST N_IMAGES representing the final converged state
                final_images_pe = pe_values[-N_IMAGES:]
                pe_max = max(final_images_pe)
                barrier = pe_max - pe0
                print(f"Success! Diffusion Barrier: {barrier:.4f} eV")
            else:
                raise ValueError(
                    f"Expected at least {N_IMAGES} potential energy values, found {len(pe_values)}."
                )
                print(f"Success! Diffusion Barrier: {barrier:.4f} eV")

    except Exception as e:
        print(f"Execution failed: {e}")

    print(f"Elapsed Time: {time.time() - t0:.1f} s\n")


if __name__ == "__main__":
    run_simple_neb()

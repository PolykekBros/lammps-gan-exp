import subprocess

N = 4

IN_FILE = "in.neb"
DUMP_LOAD = "dump.load"
DUMP_MINIMIZE = "dump.minimize"
DUMP_ATOM = "dump.atom"
DUMP_ATOM_MINIMIZE = "dump.atom_minimize"
NEB_FINAL = "final.neb"

a_id = 4001
a_x = 15.945
a_y = 29.776
a_z = 27.277
a_d = 5.1835


def exec_lmp(in_file):
    subprocess.run(
        [
            "mpirun",
            "-n",
            str(N),
            "lmp",
            "-in",
            in_file,
        ],
        text=True,
        check=True,
    )


def relax_stuff(x, y, z, n):
    in_file = f"{IN_FILE}.{n}"
    with open(in_file, mode="w") as f:
        f.write(
            "units metal\n"
            "dimension 3\n"
            "boundary p p f\n"
            "atom_style atomic\n"
            "read_data GaN_ortho_rotated.lmp\n"
            "change_box all z delta 0 10\n"
            f"write_dump all custom {DUMP_LOAD}.{n} id type x y z\n"
            "mass 1 69.723\n"
            "mass 2 14.007\n"
            "pair_style meam\n"
            "pair_coeff * * library.meam Ga N GaN.meam Ga N\n"
            "minimize 1.0e-10 1.0e-10 10000 10000\n"
            "reset_timestep 0\n"
            f"write_dump all custom {DUMP_MINIMIZE}.{n} id type x y z\n"
            f"create_atoms 1 single {x} {y} {z} units box\n"
            f"write_dump all custom {DUMP_ATOM}.{n} id type x y z\n"
            "minimize 1.0e-10 1.0e-10 10000 10000\n"
            "reset_timestep 0\n"
            f"write_dump all custom {DUMP_ATOM_MINIMIZE}.{n} id type x y z\n"
        )
    exec_lmp(in_file)


def extract_coords(dump_file):
    with open(dump_file, mode="r") as f:
        lines = f.readlines()
    for line in lines:
        tokens = line.split()
        if len(tokens) == 0:
            continue
        try:
            c_a_id = int(tokens[0])
        except ValueError:
            continue
        if c_a_id == a_id and len(tokens) == 5:
            x = float(tokens[2])
            y = float(tokens[3])
            z = float(tokens[4])
            return x, y, z
    raise f"no atom with id {a_id}"


def write_final(x, y, z):
    with open(NEB_FINAL, mode="w") as f:
        f.write(f"{1}\n{a_id} {x} {y} {z}\n")


def do_neb(n):
    with open(IN_FILE, mode="w") as f:
        f.write(
            "units metal\n"
            "dimension 3\n"
            "boundary p p f\n"
            "atom_style atomic\n"
            "atom_modify map array sort 0 0.0\n"
            "region 1 block 0 1 0 1 0 1\n"
            "create_box 2 1\n"
            f"read_dump {DUMP_ATOM_MINIMIZE}.{n} 0 id type x y z box yes add yes\n"
            "mass 1 69.723\n"
            "mass 2 14.007\n"
            "write_dump all custom dump.neb_check id type x y z\n"
            "pair_style meam\n"
            "pair_coeff * * library.meam Ga N GaN.meam Ga N\n"
            "min_style fire\n"
            "fix 1 all neb 1.0\n"
            f"neb 0.0 1.0e-10 1000 1000 100 final {NEB_FINAL}\n"
        )
    subprocess.run(
        [
            "mpirun",
            "-n",
            str(N),
            "lmp",
            "-partition",
            f"{N}x1",
            "-in",
            IN_FILE,
        ],
        text=True,
        check=True,
    )


def main():
    relax_stuff(a_x, a_y, a_z, 1)
    relax_stuff(a_x, a_y + a_d, a_z, 2)
    x, y, z = extract_coords(f"{DUMP_ATOM_MINIMIZE}.{2}")
    write_final(x, y, z)
    do_neb(1)


if __name__ == "__main__":
    main()

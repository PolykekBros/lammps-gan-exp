import lammps_mpi4py

NEB_FINAL = "neb.final"

a_type = 1
a_x = 15.945
a_y = 29.776
a_z = 27.277
a_d = 5.1835


def main(lmp):
    lmp.commands_list(
        [
            "units metal",
            "dimension 3",
            "boundary p p f",
            "atom_style atomic",
            "read_data GaN_ortho_rotated.lmp",
            "change_box all z delta 0 10",
            "write_dump all custom dump.neb_load id type x y z",
            "mass 1 69.723",
            "mass 2 14.007",
            "pair_style meam",
            "pair_coeff * * library.meam Ga N GaN.meam Ga N",
            "minimize 1.0e-10 1.0e-10 10000 10000",
            "write_dump all custom dump.neb_minimize id type x y z",
            f"create_atoms {a_type} single {a_x} {a_y} {a_z} units box",
            "write_dump all custom dump.neb_atom id type x y z",
            "minimize 1.0e-10 1.0e-10 10000 10000",
            "write_dump all custom dump.neb_atom_minimize id type x y z",
        ]
    )
    with open(NEB_FINAL, mode="w") as f:
        f.write(f"{1}\n{4001} {a_x} {a_y + a_d} {a_z}")
    lmp.command(f"neb 0.0 1.0e-10 10000 10000 50 final {NEB_FINAL}")


if __name__ == "__main__":
    lammps_mpi4py.run(main)

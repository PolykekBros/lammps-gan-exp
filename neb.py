import lammps_mpi4py

NEB_FINAL = "neb.final"

a_type = 1
a_x = 15.945
a_y = 29.776
a_z = 27.277
a_d = 5.1835


def main(lmp):
    lmp.commands_list([
    "units metal",
    "dimension 3",
    c.boundary("p p f")
    c.atom_style("atomic")
    c.read_data("GaN_ortho_rotated.lmp")
    c.change_box("all z delta", 0, 10)
    c.write_dump("all custom dump.neb_load id type x y z")
    c.mass(1, 69.723)  # Ga
    c.mass(2, 14.007)  # N
    c.pair_style("meam")
    c.pair_coeff("* * library.meam Ga N GaN.meam Ga N")
    c.minimize(1.0e-10, 1.0e-10, 10000, 10000)
    c.write_dump("all custom dump.neb_minimize id type x y z")
    c.create_atoms(a_type, "single", a_x, a_y, a_z, "units box")
    c.write_dump("all custom dump.neb_atom id type x y z")
    c.minimize(1.0e-10, 1.0e-10, 10000, 10000)
    c.write_dump("all custom dump.neb_atom_minimize id type x y z")
    ])
    with open(NEB_FINAL, mode="w") as f:
        f.write(f"{1}\n{4001} {a_x} {a_y + a_d} {a_z}")
    c.neb(0.0, 1.0e-10, 10000, 10000, 50, f"final {NEB_FINAL}")


if __name__ == "__main__":
    lammps_mpi4py.run(main)


# thermo          100
# thermo_style    custom step pe ke etotal temp lx ly lz press

# create_atoms    1 single 15.945 29.776 27.277 units box

# write_dump      all custom 001 id type x y z
# min_style       cg
# minimize        1.0e-8 1.0e-8 10000 10000
# write_dump      all custom 002 id type x y z

# reset_timestep  0
# timestep        0.001     # 1 fs timestep (metal units)

# velocity        all create 1.0e-6 87287 mom yes rot yes

# fix             1 all nve

# thermo          100
# thermo_style    custom step temp pe ke etotal press

# run             5000

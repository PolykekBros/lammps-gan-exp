clear
units           metal
dimension       3
boundary        p p f
atom_style      atomic

read_data       GaN_ortho_rotated.lmp
change_box      all z delta 0 10

mass            1 69.723
mass            2 14.007

pair_style      meam
pair_coeff      * * library.meam Ga N GaN.meam Ga N

neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes

thermo          100
thermo_style    custom step pe ke etotal temp lx ly lz press

min_style       cg
minimize        1.0e-8 1.0e-8 10000 10000

variable        zcut equal 5.0
region          upper block INF INF INF INF ${zcut} INF
group           upper region upper

variable        nx equal 20
variable        ny equal 20
variable        dx equal lx/${nx}
variable        dy equal ly/${ny}
variable        reset_y equal -${ny}*${dy}

print           "x_disp y_disp pe" file pes_surface.dat

variable        i loop ${nx}
label           loopx
variable        j loop ${ny}
label           loopy

run             0
variable        current_pe equal pe
variable        x_disp equal (${i}-1)*${dx}
variable        y_disp equal (${j}-1)*${dy}

print           "${x_disp} ${y_disp} ${current_pe}" append pes_surface.dat
displace_atoms  upper move 0 ${dy} 0

next            j
jump            SELF loopy
variable        j delete

displace_atoms  upper move ${dx} ${reset_y} 0

next            i
jump            SELF loopx

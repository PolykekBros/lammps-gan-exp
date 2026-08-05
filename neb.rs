#!/usr/bin/env -S cargo +nightly -Zscript
---
[package]
edition = "2024"
[dependencies]
anyhow = "1.0.104"
---
use anyhow::{Context, Result};
use std::{
    cell::LazyCell,
    ffi::OsStr,
    fs::File,
    io::{BufWriter, BufReader, BufRead, Write},
    iter,
    path::{Path, PathBuf},
    process::Command,
};

const N: usize = 4;

const IN_FILE: LazyCell<PathBuf> = LazyCell::new(|| PathBuf::from("in.neb"));
const DUMP_LOAD: &str = "dump.load";
const DUMP_MINIMIZE: &str = "dump.minimize";
const DUMP_ATOM: &str = "dump.atom";
const DUMP_ATOM_MINIMIZE: &str = "dump.atom_minimize";
const NEB_FINAL: &str = "final.neb";

const A_ID: usize = 4001;
const A_X: f64 = 15.945;
const A_Y: f64 = 29.776;
const A_Z: f64 = 27.277;
const A_D: f64 = 5.1835;

fn get_lmp_cmd() -> Result<&'static str> {
    Ok("lmp")
}

fn exec_lmp_with_args<P, I, S>(in_file: P, args: I) -> Result<()>
where
    P: AsRef<Path>,
    I: IntoIterator<Item = S>,
    S: AsRef<OsStr>,
{
    let lmp_cmd = get_lmp_cmd()?;
    Command::new("mpirun")
        .args(["-n", &N.to_string(), lmp_cmd, "-in"])
        .arg(in_file.as_ref())
        .args(args)
        .status()?;
    Ok(())
}

fn exec_lmp<P>(in_file: P) -> Result<()>
where
    P: AsRef<Path>,
{
    exec_lmp_with_args(in_file, iter::empty::<String>())
}

fn relax_stuff(coords: [f64; 3], n: usize) -> Result<()> {
    let x = coords[0];
    let y = coords[1];
    let z = coords[2];
    let in_file = IN_FILE.with_added_extension(&n.to_string());
    let file = File::open(&in_file)
        .with_context(|| format!("file: {}", in_file.display()))?;
    let mut w = BufWriter::new(file);
    writeln!(w, "units metal")?;
    writeln!(w, "units metal")?;
    writeln!(w, "dimension 3")?;
    writeln!(w, "boundary p p f")?;
    writeln!(w, "atom_style atomic")?;
    writeln!(w, "read_data GaN_ortho_rotated.lmp")?;
    writeln!(w, "change_box all z delta 0 10")?;
    writeln!(w, "write_dump all custom {DUMP_LOAD}.{n} id type x y z")?;
    writeln!(w, "mass 1 69.723")?;
    writeln!(w, "mass 2 14.007")?;
    writeln!(w, "pair_style meam")?;
    writeln!(w, "pair_coeff * * library.meam Ga N GaN.meam Ga N")?;
    writeln!(w, "minimize 1.0e-10 1.0e-10 10000 10000")?;
    writeln!(w, "reset_timestep 0")?;
    writeln!(w, "write_dump all custom {DUMP_MINIMIZE}.{n} id type x y z")?;
    writeln!(w, "create_atoms 1 single {x} {y} {z} units box")?;
    writeln!(w, "write_dump all custom {DUMP_ATOM}.{n} id type x y z")?;
    writeln!(w, "minimize 1.0e-10 1.0e-10 10000 10000")?;
    writeln!(w, "reset_timestep 0")?;
    writeln!(
        w,
        "write_dump all custom {DUMP_ATOM_MINIMIZE}.{n} id type x y z"
    )?;
    exec_lmp(in_file)
}

fn extract_coords<P>(dump_file: P) -> Result<[f64; 3]>
where P: AsRef<Path> {
    let dump_file = dump_file.as_ref();
    let dump_file = File::open(&dump_file)
        .with_context(|| format!("file {}", dump_file.display()))?;
    let reader = BufReader::new(dump_file); 
    for line in reader.lines() {
       let tokens = line?.whitespace_split();

    }
    /*
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
    */
    Ok([0.0,0.0,0.0])
}

fn main() -> Result<()> {
    relax_stuff([A_X, A_Y, A_Z], 1)?;
    relax_stuff([A_X, A_Y + A_D, A_Z], 2)?;
    Ok(())
}

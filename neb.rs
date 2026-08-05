#!/usr/bin/env -S cargo +nightly -Zscript
---
[package]
edition = "2024"
[dependencies]
anyhow = "1.0.104"
---
#![feature(frontmatter)]
use anyhow::{bail, Context, Result};
use std::{
    cell::LazyCell,
    ffi::OsStr,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    iter,
    path::{Path, PathBuf},
    process::Command,
};

const N: usize = 4;

const IN_FILE: LazyCell<PathBuf> = LazyCell::new(|| PathBuf::from("in.neb"));
const DUMP_LOAD: LazyCell<PathBuf> =
    LazyCell::new(|| PathBuf::from("dump.load"));
const DUMP_MINIMIZE: LazyCell<PathBuf> =
    LazyCell::new(|| PathBuf::from("dump.minimize"));
const DUMP_ATOM: LazyCell<PathBuf> =
    LazyCell::new(|| PathBuf::from("dump.atom"));
const DUMP_ATOM_MINIMIZE: LazyCell<PathBuf> =
    LazyCell::new(|| PathBuf::from("dump.atom_minimize"));
const NEB_FINAL: LazyCell<PathBuf> =
    LazyCell::new(|| PathBuf::from("final.neb"));

const A_ID: usize = 4001;
const A_X: f64 = 15.945;
const A_Y: f64 = 29.776;
const A_Z: f64 = 27.277;
const A_D: f64 = 5.1835;

fn get_lmp_cmd() -> Result<&'static str> {
    Command::new("lmp_mpi")
        .arg("-h")
        .status()
        .map(|_| "lmp_mpi")
        .or_else(|e| match e.kind() {
            std::io::ErrorKind::NotFound => {
                Command::new("lmp").arg("-h").status().map(|_| "lmp")
            }
            _ => Err(e),
        })
        .map_err(anyhow::Error::from)
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
    let dump_load = DUMP_LOAD.with_added_extension(&n.to_string());
    let dump_minimize = DUMP_MINIMIZE.with_added_extension(&n.to_string());
    let dump_atom = DUMP_ATOM.with_added_extension(&n.to_string());
    let dump_atom_minimize =
        DUMP_ATOM_MINIMIZE.with_added_extension(&n.to_string());
    let file = File::create(&in_file)
        .with_context(|| format!("file: {}", in_file.display()))?;
    let mut w = BufWriter::new(file);
    writeln!(w, "units metal")?;
    writeln!(w, "dimension 3")?;
    writeln!(w, "boundary p p f")?;
    writeln!(w, "atom_style atomic")?;
    writeln!(w, "read_data GaN_ortho_rotated.lmp")?;
    writeln!(w, "change_box all z delta 0 10")?;
    writeln!(
        w,
        "write_dump all custom {} id type x y z",
        dump_load.display()
    )?;
    writeln!(w, "mass 1 69.723")?;
    writeln!(w, "mass 2 14.007")?;
    writeln!(w, "pair_style meam")?;
    writeln!(w, "pair_coeff * * library.meam Ga N GaN.meam Ga N")?;
    writeln!(w, "minimize 1.0e-10 1.0e-10 10000 10000")?;
    writeln!(w, "reset_timestep 0")?;
    writeln!(
        w,
        "write_dump all custom {} id type x y z",
        dump_minimize.display()
    )?;
    writeln!(w, "create_atoms 1 single {x} {y} {z} units box")?;
    writeln!(
        w,
        "write_dump all custom {} id type x y z",
        dump_atom.display()
    )?;
    writeln!(w, "minimize 1.0e-10 1.0e-10 10000 10000")?;
    writeln!(w, "reset_timestep 0")?;
    writeln!(
        w,
        "write_dump all custom {} id type x y z",
        dump_atom_minimize.display()
    )?;
    exec_lmp(in_file)
}

fn extract_coords<P>(dump_file: P) -> Result<[f64; 3]>
where
    P: AsRef<Path>,
{
    let dump_file = dump_file.as_ref();
    let dump_file = File::open(&dump_file)
        .with_context(|| format!("file {}", dump_file.display()))?;
    let reader = BufReader::new(dump_file);
    for line in reader.lines() {
        let line = line?;
        let tokens = line.split_whitespace().collect::<Vec<_>>();
        let idx = 0;
        let Some(Ok(a_id)) = tokens
            .get(idx)
            .map(|s| str::parse::<usize>(s))
        else {
            continue;
        };
        if a_id != A_ID {
            continue;
        }
        let begin = idx + 2;
        let end = begin + 3;
        let mut coords = [0f64; 3];
        for (idx, coord) in tokens
            .get(begin..end)
            .with_context(|| format!("columns {} trough {}", begin, end))?
            .iter()
            .enumerate()
        {
            coords[idx] = coord.parse()?;
        }
        return Ok(coords);
    }
    bail!("could not find coords for atom id {A_ID}");
}

fn write_final(coords: [f64; 3]) -> Result<()> {
    let file = File::create(NEB_FINAL.as_path())?;
    let mut w = BufWriter::new(file);
    writeln!(w, "1")?;
    writeln!(w, "{A_ID} {} {} {}", coords[0], coords[1], coords[2])?;
    Ok(())
}

fn do_neb(n: usize) -> Result<()> {
    let file = File::create(IN_FILE.as_path())?;
    let mut w = BufWriter::new(file);
    writeln!(w, "units metal")?;
    writeln!(w, "dimension 3")?;
    writeln!(w, "boundary p p f")?;
    writeln!(w, "atom_style atomic")?;
    writeln!(w, "atom_modify map array sort 0 0.0")?;
    writeln!(w, "region 1 block 0 1 0 1 0 1")?;
    writeln!(w, "create_box 2 1")?;
    writeln!(
        w,
        "read_dump {} 0 id type x y z box yes add yes",
        DUMP_ATOM_MINIMIZE
            .with_added_extension(&n.to_string())
            .display()
    )?;
    writeln!(w, "mass 1 69.723")?;
    writeln!(w, "mass 2 14.007")?;
    writeln!(w, "write_dump all custom dump.neb_check id type x y z")?;
    writeln!(w, "pair_style meam")?;
    writeln!(w, "pair_coeff * * library.meam Ga N GaN.meam Ga N")?;
    writeln!(w, "min_style fire")?;
    writeln!(w, "fix 1 all neb 1.0")?;
    writeln!(
        w,
        "neb 0.0 1.0e-10 1000 1000 100 final {}",
        NEB_FINAL.display()
    )?;
    exec_lmp_with_args(IN_FILE.as_path(), ["-partition", &format!("{N}x1")])
}

fn main() -> Result<()> {
    relax_stuff([A_X, A_Y, A_Z], 1)?;
    relax_stuff([A_X, A_Y + A_D, A_Z], 2)?;
    let coords = extract_coords(
        DUMP_ATOM_MINIMIZE.with_added_extension(&2.to_string()),
    )?;
    write_final(coords)?;
    do_neb(1)?;
    Ok(())
}

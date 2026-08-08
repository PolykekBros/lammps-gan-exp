#!/usr/bin/env -S cargo +nightly -Zscript
---
[package]
edition = "2024"
[dependencies]
anyhow = "1.0.104"
---
use anyhow::{bail, Context, Result};
use std::cell::OnceCell;
use std::{
    cell::LazyCell,
    ffi::OsStr,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    iter,
    path::{Path, PathBuf},
    process::Command,
};

const N: usize = 12;

const IN_FILE: LazyCell<PathBuf> = LazyCell::new(|| PathBuf::from("in.neb"));
const DUMP_LOAD: LazyCell<PathBuf> =
    LazyCell::new(|| PathBuf::from("dump.load"));
const DUMP_MINIMIZE: LazyCell<PathBuf> =
    LazyCell::new(|| PathBuf::from("dump.minimize"));
const DUMP_ATOM: LazyCell<PathBuf> =
    LazyCell::new(|| PathBuf::from("dump.atom"));
const NEB_DATA: LazyCell<PathBuf> = LazyCell::new(|| PathBuf::from("data.neb"));
const NEB_FINAL: LazyCell<PathBuf> =
    LazyCell::new(|| PathBuf::from("final.neb"));

#[derive(Debug, Clone)]
struct InputWriter {
    path: PathBuf,
    file: String,
}

#[derive(Debug, Clone)]
struct LammpsCMDFinder {
    candidates: Vec<&'static str>,
    cmd: OnceCell<Option<&'static str>>,
}

impl Default for LammpsCMDFinder {
    fn default() -> Self {
        Self {
            candidates: vec!["lmp_mpi", "lmp", "lmp_serial"],
            cmd: OnceCell::new(),
        }
    }
}

impl LammpsCMDFinder {
    fn find(&self) -> Result<&'static str> {
        self.cmd
            .get_or_init(|| {
                let args = ["-h", "-screen", "none"];
                for variant in &self.candidates {
                    let err = match Command::new(variant).args(&args).status() {
                        Ok(_) => return Some(variant),
                        Err(err) => err,
                    };
                    match err.kind() {
                        std::io::ErrorKind::NotFound => continue,
                        _ => return Some(variant),
                    }
                }
                None
            })
            .with_context(|| {
                format!(
                    "could not find any LAMMPS executables: {:?}",
                    self.candidates
                )
            })
    }
}

#[derive(Debug, Clone, Default)]
struct LammpsExecutor {
    cmd_finder: LammpsCMDFinder,
}

impl LammpsExecutor {
    fn exec_with_args<P, I, S>(&self, in_file: P, args: I) -> Result<()>
    where
        P: AsRef<Path>,
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        let cmd = self.cmd_finder.find()?;
        let status = Command::new("mpirun")
            .args(["-n", &N.to_string(), cmd, "-in"])
            .arg(in_file.as_ref())
            .args(args)
            .status()?;
        if !status.success() {
            bail!("LAMMPS exited with status {}", status);
        }
        Ok(())
    }

    fn exec<P: AsRef<Path>>(&self, in_file: P) -> Result<()> {
        self.exec_with_args(in_file, iter::empty::<String>())
    }
}

const A_ID: usize = 4001;
const A_X: f64 = 15.945;
const A_Y: f64 = 29.776;
const A_Z: f64 = 27.277;
const A_D: f64 = 5.1835;

fn relax_stuff(coords: [f64; 3], n: usize) -> Result<()> {
    let x = coords[0];
    let y = coords[1];
    let z = coords[2];
    let in_file = IN_FILE.with_added_extension(&n.to_string());
    let dump_load = DUMP_LOAD.with_added_extension(&n.to_string());
    let dump_minimize = DUMP_MINIMIZE.with_added_extension(&n.to_string());
    let dump_atom = DUMP_ATOM.with_added_extension(&n.to_string());
    let neb_data = NEB_DATA.with_added_extension(&n.to_string());
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
    writeln!(w, "write_data {}", neb_data.display())?;
    w.flush()?;
    LammpsExecutor::default().exec(in_file)
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
        if tokens.len() != 8 {
            continue;
        }
        let idx = 0;
        let a_id = tokens
            .get(idx)
            .map(|s| s.parse::<usize>())
            .transpose()?
            .context("no id column")?;
        if a_id != A_ID {
            continue;
        }
        let begin = idx + 2;
        let end = begin + 3;
        let mut coords = [0f64; 3];
        for (idx, coord) in tokens
            .get(begin..end)
            .with_context(|| format!("no columns {} trough {}", begin, end))?
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
    writeln!(
        w,
        "read_data {}",
        NEB_DATA.with_added_extension(&n.to_string()).display()
    )?;
    writeln!(w, "pair_style meam")?;
    writeln!(w, "pair_coeff * * library.meam Ga N GaN.meam Ga N")?;
    writeln!(w, "region fixed block INF INF INF INF INF $(zlo+5.0)")?;
    writeln!(w, "group fixed region fixed")?;
    writeln!(w, "group mobile subtract all fixed")?;
    writeln!(w, "min_style fire")?;
    writeln!(w, "min_modify dmax 0.05")?;
    writeln!(w, "fix fixed fixed setforce 0.0 0.0 0.0")?;
    writeln!(w, "fix neb mobile neb 0.5 parallel ideal")?;
    writeln!(w, "variable P uloop {N}")?;
    writeln!(
        w,
        "dump neb all custom 100 dump.neb.$P id type xu yu zu fx fy fz"
    );
    writeln!(
        w,
        "neb 0.0 0.05 10000 10000 1000 final {}",
        NEB_FINAL.display()
    )?;
    w.flush()?;
    LammpsExecutor::default()
        .exec_with_args(IN_FILE.as_path(), ["-partition", &format!("{N}x1")])
}

fn main() -> Result<()> {
    relax_stuff([A_X, A_Y, A_Z], 1)?;
    relax_stuff([A_X, A_Y + A_D, A_Z], 2)?;
    let coords = extract_coords(NEB_DATA.with_added_extension(&2.to_string()))?;
    write_final(coords)?;
    do_neb(1)?;
    Ok(())
}

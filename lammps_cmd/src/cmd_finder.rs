use std::cell::{LazyCell, OnceCell};
use std::process::{Command, Stdio};

use crate::error::Error;

#[derive(Debug)]
struct LAMMPSCMDChecker {
    cmd: String,
}

impl LAMMPSCMDChecker {
    fn new(cmd: String) -> Self {
        Self { cmd }
    }

    fn check(self) -> Result<String, Error> {
        let status = Command::new(&self.cmd)
            .stdin(Stdio::null())
            .status()
            .map_err(|err| Error::CMDExecError(self.cmd.clone(), err))?;
        status
            .success()
            .ok_or_else(|| Error::CMDExecFailure(self.cmd.clone(), status))
            .map(|_| self.cmd)
    }
}

pub trait CMDFinder {
    fn find(&self) -> Result<String, Error>;
}

#[derive(Debug)]
pub struct LAMMPSCMDFinder {
    candidates: Vec<String>,
    cmd: OnceCell<Option<String>>,
}

impl Default for LAMMPSCMDFinder {
    fn default() -> Self {
        Self::new(vec![
            String::from("lmp_mpi"),
            String::from("lmp"),
            String::from("lmp_serial"),
        ])
    }
}

impl LAMMPSCMDFinder {
    fn new(candidates: Vec<String>) -> Self {
        Self {
            candidates,
            cmd: OnceCell::new(),
        }
    }
}

impl CMDFinder for LAMMPSCMDFinder {
    fn find(&self) -> Result<String, Error> {
        self.cmd
            .get_or_init(|| {
                self.candidates
                    .iter()
                    .cloned()
                    .filter_map(|cmd| LAMMPSCMDChecker::new(cmd).check().ok())
                    .next()
            })
            .clone()
            .ok_or_else(|| Error::NoSuitableCMDs(self.candidates.clone()))
    }
}
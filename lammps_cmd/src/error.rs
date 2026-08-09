use std::{
    io,
    process::{ExitStatus},
};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("Could not find any LAMMPS executable, tried: {0:?}")]
    NoSuitableCMDs(Vec<String>),
    #[error("Could not execute command `{0}`: {1}")]
    CMDExecError(String, io::Error),
    #[error("Command `{0}` failed: {1}")]
    CMDExecFailure(String, ExitStatus),
}
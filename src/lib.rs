pub mod d2q9;
pub mod d3q27;
pub mod global_variables;
pub mod io;
pub mod post;

pub use global_variables::*;

#[derive(Copy, Clone, PartialEq)]
pub enum NodeType {
    Fluid = 0,
    Solid = 1,
}

#[derive(Clone)]
pub struct Residuals {
    pub density: Float,
    pub velocity: Vec<Float>,
}

pub trait Lattice {
    fn compute_residuals();
}

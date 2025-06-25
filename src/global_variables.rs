pub const CASE_NAME: &'static str = "Case Test";

pub const IGNORE_SOLID: bool = true;

pub type Float = f64;

pub const TOLERANCE_DENSITY: Float = 1e-7;

pub const TOLERANCE_VELOCITY_X: Float = 1e-7;

pub const TOLERANCE_VELOCITY_Y: Float = 1e-7;

pub const TOLERANCE_VELOCITY_Z: Float = 1e-7;

pub const MIN_ITER: usize = 1000;

pub const MAX_ITER: usize = 100_000;

pub const DELTA_T: Float = 1.0;

pub const DELTA_X: Float = 1.0;

pub const LATTICE_DENSITY: Float = 1.0;

pub const CS_2: Float = 1.0 / 3.0 * DELTA_X * DELTA_X / DELTA_T / DELTA_T;

pub const CS_2_INV: Float = 3.0;

pub const CS_4_INV: Float = 9.0;

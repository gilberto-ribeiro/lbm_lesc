use super::*;
use crate::global_variables::*;

const Q_BAR: [usize; Q] = [
    0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26,
    25,
];

const Q_EAST: [usize; 9] = [1, 7, 9, 13, 15, 19, 21, 23, 26];

const Q_WEST: [usize; 9] = [2, 8, 10, 14, 16, 20, 22, 24, 25];

const Q_NORTH: [usize; 9] = [3, 7, 11, 14, 17, 19, 21, 24, 25];

const Q_SOUTH: [usize; 9] = [4, 8, 12, 13, 18, 20, 22, 23, 26];

const Q_TOP: [usize; 9] = [5, 9, 11, 16, 18, 19, 22, 23, 25];

const Q_BOTTOM: [usize; 9] = [6, 10, 12, 15, 17, 20, 21, 24, 26];

#[derive(Clone)]
pub enum BoundaryCondition {
    NoSlip,
    Dirichlet {
        density: Float,
        velocity: [Float; D],
    },
    PressureOutlet {
        density: Float,
    },
    Periodic,
}

#[derive(Hash, Eq, PartialEq, Clone, Copy)]
pub enum BoundaryFace {
    East,
    West,
    North,
    South,
    Top,
    Bottom,
}

impl Lattice {
    pub fn apply_boundary_condition(
        &mut self,
        boundary_face: BoundaryFace,
        boundary_condition: BoundaryCondition,
    ) {
        match boundary_condition {
            BoundaryCondition::NoSlip => self.no_slip(boundary_face),
            BoundaryCondition::Dirichlet { density, velocity } => {
                self.dirichlet(boundary_face, density, velocity)
            }
            BoundaryCondition::PressureOutlet { density } => {
                self.pressure_outlet(boundary_face, density)
            }
            BoundaryCondition::Periodic => {}
        }
    }

    fn no_slip(&mut self, boundary_face: BoundaryFace) {
        let q_face = match boundary_face {
            BoundaryFace::East => Q_EAST,
            BoundaryFace::West => Q_WEST,
            BoundaryFace::North => Q_NORTH,
            BoundaryFace::South => Q_SOUTH,
            BoundaryFace::Top => Q_TOP,
            BoundaryFace::Bottom => Q_BOTTOM,
        };
        self.nodes
            .par_iter_mut()
            .filter(|node| {
                node.boundary_faces
                    .as_ref()
                    .unwrap()
                    .contains(&boundary_face)
            })
            .for_each(|node| {
                inner_no_slip(node, q_face);
            });
    }

    fn dirichlet(&mut self, boundary_face: BoundaryFace, density: Float, velocity: [Float; D]) {
        let q_face = match boundary_face {
            BoundaryFace::East => Q_EAST,
            BoundaryFace::West => Q_WEST,
            BoundaryFace::North => Q_NORTH,
            BoundaryFace::South => Q_SOUTH,
            BoundaryFace::Top => Q_TOP,
            BoundaryFace::Bottom => Q_BOTTOM,
        };
        self.nodes
            .par_iter_mut()
            .filter(|node| {
                node.boundary_faces
                    .as_ref()
                    .unwrap()
                    .contains(&boundary_face)
            })
            .for_each(|node| match node.node_type {
                NodeType::Fluid => {
                    for q in q_face {
                        let u_dot_c = velocity
                            .iter()
                            .zip(C[q].iter())
                            .map(|(&u, &c)| u * (c as Float))
                            .sum::<Float>();
                        node.f[Q_BAR[q]] =
                            node.f_star[q] - 2.0 * W[q] * density * CS_2_INV * u_dot_c;
                    }
                }
                NodeType::Solid => {
                    inner_no_slip(node, q_face);
                }
            });
    }

    fn pressure_outlet(&mut self, boundary_face: BoundaryFace, density: Float) {
        match boundary_face {
            BoundaryFace::East => {
                for j in 0..self.ny {
                    for k in 0..self.nz {
                        match self.get_node(&[self.nx - 1, j, k]).node_type {
                            NodeType::Fluid => {
                                let ux = 1.5 * self.get_node(&[self.nx - 1, j, k]).velocity[0]
                                    - 0.5 * self.get_node(&[self.nx - 2, j, k]).velocity[0];
                                let uy = 1.5 * self.get_node(&[self.nx - 1, j, k]).velocity[1]
                                    - 0.5 * self.get_node(&[self.nx - 2, j, k]).velocity[1];
                                let uz = 1.5 * self.get_node(&[self.nx - 1, j, k]).velocity[2]
                                    - 0.5 * self.get_node(&[self.nx - 2, j, k]).velocity[2];
                                let u_2 = ux * ux + uy * uy + uz * uz;
                                for q in Q_EAST {
                                    let cx = C[q][0] as Float;
                                    let cy = C[q][1] as Float;
                                    let cz = C[q][2] as Float;
                                    let u_dot_c = cx * ux + cy * uy + cz * uz;
                                    self.get_node_mut(&[self.nx - 1, j, k]).f[Q_BAR[q]] =
                                        -self.get_node(&[self.nx - 1, j, k]).f_star[q]
                                            + 2.0
                                                * W[q]
                                                * density
                                                * (1.0 + 0.5 * CS_4_INV * u_dot_c * u_dot_c
                                                    - 0.5 * CS_2_INV * u_2);
                                }
                            }
                            NodeType::Solid => {
                                for q in Q_EAST {
                                    self.get_node_mut(&[self.nx - 1, j, k]).f[Q_BAR[q]] =
                                        self.get_node(&[self.nx - 1, j, k]).f_star[q];
                                }
                            }
                        }
                    }
                }
            }
            BoundaryFace::West => {
                for j in 0..self.ny {
                    for k in 0..self.nz {
                        match self.get_node(&[0, j, k]).node_type {
                            NodeType::Fluid => {
                                let ux = 1.5 * self.get_node(&[0, j, k]).velocity[0]
                                    - 0.5 * self.get_node(&[1, j, k]).velocity[0];
                                let uy = 1.5 * self.get_node(&[0, j, k]).velocity[1]
                                    - 0.5 * self.get_node(&[1, j, k]).velocity[1];
                                let uz = 1.5 * self.get_node(&[0, j, k]).velocity[2]
                                    - 0.5 * self.get_node(&[1, j, k]).velocity[2];
                                let u_2 = ux * ux + uy * uy + uz * uz;
                                for q in Q_WEST {
                                    let cx = C[q][0] as Float;
                                    let cy = C[q][1] as Float;
                                    let cz = C[q][2] as Float;
                                    let u_dot_c = cx * ux + cy * uy + cz * uz;
                                    self.get_node_mut(&[0, j, k]).f[Q_BAR[q]] =
                                        -self.get_node(&[0, j, k]).f_star[q]
                                            + 2.0
                                                * W[q]
                                                * density
                                                * (1.0 + 0.5 * CS_4_INV * u_dot_c * u_dot_c
                                                    - 0.5 * CS_2_INV * u_2);
                                }
                            }
                            NodeType::Solid => {
                                for q in Q_WEST {
                                    self.get_node_mut(&[0, j, k]).f[Q_BAR[q]] =
                                        self.get_node(&[0, j, k]).f_star[q];
                                }
                            }
                        }
                    }
                }
            }
            BoundaryFace::North => {
                for i in 0..self.nx {
                    for k in 0..self.nz {
                        match self.get_node(&[i, self.ny - 1, k]).node_type {
                            NodeType::Fluid => {
                                let ux = 1.5 * self.get_node(&[i, self.ny - 1, k]).velocity[0]
                                    - 0.5 * self.get_node(&[i, self.ny - 2, k]).velocity[0];
                                let uy = 1.5 * self.get_node(&[i, self.ny - 1, k]).velocity[1]
                                    - 0.5 * self.get_node(&[i, self.ny - 2, k]).velocity[1];
                                let uz = 1.5 * self.get_node(&[i, self.ny - 1, k]).velocity[2]
                                    - 0.5 * self.get_node(&[i, self.ny - 2, k]).velocity[2];
                                let u_2 = ux * ux + uy * uy + uz * uz;
                                for q in Q_NORTH {
                                    let cx = C[q][0] as Float;
                                    let cy = C[q][1] as Float;
                                    let cz = C[q][2] as Float;
                                    let u_dot_c = cx * ux + cy * uy + cz * uz;
                                    self.get_node_mut(&[i, self.ny - 1, k]).f[Q_BAR[q]] =
                                        -self.get_node(&[i, self.ny - 1, k]).f_star[q]
                                            + 2.0
                                                * W[q]
                                                * density
                                                * (1.0 + 0.5 * CS_4_INV * u_dot_c * u_dot_c
                                                    - 0.5 * CS_2_INV * u_2);
                                }
                            }
                            NodeType::Solid => {
                                for q in Q_NORTH {
                                    self.get_node_mut(&[i, self.ny - 1, k]).f[Q_BAR[q]] =
                                        self.get_node(&[i, self.ny - 1, k]).f_star[q];
                                }
                            }
                        }
                    }
                }
            }
            BoundaryFace::South => {
                for i in 0..self.nx {
                    for k in 0..self.nz {
                        match self.get_node(&[i, 0, k]).node_type {
                            NodeType::Fluid => {
                                let ux = 1.5 * self.get_node(&[i, 0, k]).velocity[0]
                                    - 0.5 * self.get_node(&[i, 1, k]).velocity[0];
                                let uy = 1.5 * self.get_node(&[i, 0, k]).velocity[1]
                                    - 0.5 * self.get_node(&[i, 1, k]).velocity[1];
                                let uz = 1.5 * self.get_node(&[i, 0, k]).velocity[2]
                                    - 0.5 * self.get_node(&[i, 1, k]).velocity[2];
                                let u_2 = ux * ux + uy * uy + uz * uz;
                                for q in Q_SOUTH {
                                    let cx = C[q][0] as Float;
                                    let cy = C[q][1] as Float;
                                    let cz = C[q][2] as Float;
                                    let u_dot_c = cx * ux + cy * uy + cz * uz;
                                    self.get_node_mut(&[i, 0, k]).f[Q_BAR[q]] =
                                        -self.get_node(&[i, 0, k]).f_star[q]
                                            + 2.0
                                                * W[q]
                                                * density
                                                * (1.0 + 0.5 * CS_4_INV * u_dot_c * u_dot_c
                                                    - 0.5 * CS_2_INV * u_2);
                                }
                            }
                            NodeType::Solid => {
                                for q in Q_SOUTH {
                                    self.get_node_mut(&[i, 0, k]).f[Q_BAR[q]] =
                                        self.get_node(&[i, 0, k]).f_star[q];
                                }
                            }
                        }
                    }
                }
            }
            BoundaryFace::Top => {
                for i in 0..self.nx {
                    for j in 0..self.ny {
                        match self.get_node(&[i, j, self.nz - 1]).node_type {
                            NodeType::Fluid => {
                                let ux = 1.5 * self.get_node(&[i, j, self.nz - 1]).velocity[0]
                                    - 0.5 * self.get_node(&[i, j, self.nz - 2]).velocity[0];
                                let uy = 1.5 * self.get_node(&[i, j, self.nz - 1]).velocity[1]
                                    - 0.5 * self.get_node(&[i, j, self.nz - 2]).velocity[1];
                                let uz = 1.5 * self.get_node(&[i, j, self.nz - 1]).velocity[2]
                                    - 0.5 * self.get_node(&[i, j, self.nz - 2]).velocity[2];
                                let u_2 = ux * ux + uy * uy + uz * uz;
                                for q in Q_TOP {
                                    let cx = C[q][0] as Float;
                                    let cy = C[q][1] as Float;
                                    let cz = C[q][2] as Float;
                                    let u_dot_c = cx * ux + cy * uy + cz * uz;
                                    self.get_node_mut(&[i, j, self.nz - 1]).f[Q_BAR[q]] =
                                        -self.get_node(&[i, j, self.nz - 1]).f_star[q]
                                            + 2.0
                                                * W[q]
                                                * density
                                                * (1.0 + 0.5 * CS_4_INV * u_dot_c * u_dot_c
                                                    - 0.5 * CS_2_INV * u_2);
                                }
                            }
                            NodeType::Solid => {
                                for q in Q_TOP {
                                    self.get_node_mut(&[i, j, self.nz - 1]).f[Q_BAR[q]] =
                                        self.get_node(&[i, j, self.nz - 1]).f_star[q];
                                }
                            }
                        }
                    }
                }
            }
            BoundaryFace::Bottom => {
                for i in 0..self.nx {
                    for j in 0..self.ny {
                        match self.get_node(&[i, j, 0]).node_type {
                            NodeType::Fluid => {
                                let ux = 1.5 * self.get_node(&[i, j, 0]).velocity[0]
                                    - 0.5 * self.get_node(&[i, j, 1]).velocity[0];
                                let uy = 1.5 * self.get_node(&[i, j, 0]).velocity[1]
                                    - 0.5 * self.get_node(&[i, j, 1]).velocity[1];
                                let uz = 1.5 * self.get_node(&[i, j, 0]).velocity[2]
                                    - 0.5 * self.get_node(&[i, j, 1]).velocity[2];
                                let u_2 = ux * ux + uy * uy + uz * uz;
                                for q in Q_BOTTOM {
                                    let cx = C[q][0] as Float;
                                    let cy = C[q][1] as Float;
                                    let cz = C[q][2] as Float;
                                    let u_dot_c = cx * ux + cy * uy + cz * uz;
                                    self.get_node_mut(&[i, j, 0]).f[Q_BAR[q]] =
                                        -self.get_node(&[i, j, 0]).f_star[q]
                                            + 2.0
                                                * W[q]
                                                * density
                                                * (1.0 + 0.5 * CS_4_INV * u_dot_c * u_dot_c
                                                    - 0.5 * CS_2_INV * u_2);
                                }
                            }
                            NodeType::Solid => {
                                for q in Q_BOTTOM {
                                    self.get_node_mut(&[i, j, 0]).f[Q_BAR[q]] =
                                        self.get_node(&[i, j, 0]).f_star[q];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

fn inner_no_slip(node: &mut Node, q_face: [usize; 9]) {
    for q in q_face {
        node.f[Q_BAR[q]] = node.f_star[q];
    }
}

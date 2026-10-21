use crate::{Float, P2, V2, configs::NodeConfig, motor::Motor, spring::Spring};
use std::{fmt::Display, sync::RwLock, sync::Weak};

pub type NodeId = usize;
pub const NEAREST_NODE_OVER_SPRING_OFFSET: Float = 4.0;

#[derive(Debug, Clone)]
pub struct Attached {
    pub spring: Weak<RwLock<Spring>>,
    pub other_node_id: NodeId,
}

impl Attached {
    pub fn relative_position(&self, nodes: &[Node], main_node_position: P2) -> V2 {
        let node_ix = nodes
            .binary_search_by_key(&self.other_node_id, |n| n.id)
            .expect("relative_position(...) => missing node");
        nodes[node_ix].position() - main_node_position
    }

    pub fn length(&self, nodes: &[Node], main_node_position: P2) -> Float {
        self.relative_position(nodes, main_node_position).length()
    }

    pub fn is_destruct(&self, nodes: &[Node], main_node_position: P2) -> bool {
        if let Some(spring_ref) = self.spring.upgrade() {
            if let Ok(spring) = spring_ref.read() {
                !((spring.min_len * spring.resting_length_factor)
                    ..=(spring.max_len * spring.resting_length_factor))
                    .contains(&self.length(nodes, main_node_position))
            } else {
                false
            }
        } else {
            false
        }
    }

    pub fn print_nodes(&self, nodes: &[Node], main_node: &Node) {
        let node_ix = nodes
            .binary_search_by_key(&self.other_node_id, |n| n.id)
            .expect("relative_position(...) => missing node");
        print!("Node@({})", nodes[node_ix]);
        print!("Node@({:?})", main_node.position().to_array());
        println!();
    }
}

#[derive(Debug, Clone)]
pub struct Node {
    pub id: NodeId,
    pub radius: Float,
    pub position: P2,
    pub velocity: V2,
    pub acceleration: V2,
    pub weight: Float,
    stationary: bool,
    motor: Option<Motor>,
}

impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Node {
    pub fn new(id: NodeId, position: P2, radius: Float) -> Node {
        Node {
            id,
            radius,
            position,
            velocity: V2::new(0., 0.),
            acceleration: V2::new(0., 0.),
            weight: radius * radius,
            stationary: false,
            motor: None,
        }
    }

    pub fn is_stationary(&self) -> bool {
        self.stationary
    }

    pub fn set_stationary(&mut self, stationary: bool) {
        self.stationary = stationary;
        if !stationary {
            self.motor = None;
        }
    }

    pub fn get_radius(&self) -> Float {
        self.radius
    }

    pub fn set_motor(&mut self, motor: Option<Motor>) {
        self.motor = motor;
        if motor.is_some() {
            self.stationary = true;
        }
    }

    pub fn get_motor(&self) -> Option<Motor> {
        self.motor
    }

    pub fn bounds_collision(&mut self, width: Float, height: Float) {
        match (
            self.position()[0] < self.radius,
            self.position()[0] > (width - self.radius),
        ) {
            (true, _) => {
                if self.velocity[0] < 0. {
                    self.velocity[0] = -self.velocity[0];
                }
                self.position()[0] = self.radius;
            }
            (_, true) => {
                if self.velocity[0] > 0. {
                    self.velocity[0] = -self.velocity[0];
                }
                self.position()[0] = width - self.radius;
            }
            (false, false) => { // actually no collision occurred => do nothing
            }
        }
        match (
            self.position()[1] < self.radius,
            self.position()[1] > (height - self.radius),
        ) {
            (true, _) => {
                if self.velocity[1] < 0. {
                    self.velocity[1] = -self.velocity[1];
                }
                self.position()[1] = self.radius;
            }
            (_, true) => {
                if self.velocity[1] > 0. {
                    self.velocity[1] = -self.velocity[1];
                }
                self.position()[1] = height - self.radius;
            }
            (false, false) => { // actually no collision occurred => do nothing
            }
        }
    }

    pub fn configure(&mut self, config: NodeConfig) {
        self.stationary = config.stationary;
        self.radius = config.node_radius;
        self.motor = config.motor;
        self.weight = self.radius.powi(2);
    }

    pub fn position(&self) -> P2 {
        if self.stationary {
            if let Some(motor) = self.motor {
                motor.position(self.position)
            } else {
                self.position
            }
        } else {
            self.position
        }
    }
}

impl Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "velocity: {:.2}", self.velocity)?;
        writeln!(f, "acc: {:.2}", self.acceleration.length())?;
        writeln!(f, "weight: {:.2}", self.weight)
    }
}

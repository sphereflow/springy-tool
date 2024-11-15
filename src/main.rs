use ::glam::DVec2;
use ::glam::Vec2;
use configs::GridConfig;
use configs::GridType;
use configs::NodeConfig;
use core::f64;
use egui::Slider;
use egui::Ui;
use egui::Window;
use grid::Grid;
use macroquad::math::Rect;
use macroquad::prelude::*;
use motor::Motor;
use std::cmp::Ordering;
pub(crate) use std::{
    cell::RefCell,
    fmt::Display,
    rc::{Rc, Weak},
};
pub mod configs;
pub mod grid;
pub mod motor;

type Float = f64;
type V2 = DVec2;
type P2 = DVec2;
type NodeId = usize;

#[derive(Debug, Clone)]
pub struct Contraption {
    pub nodes: Vec<Node>,
    next_node_id: NodeId,
    pub springs: Vec<Rc<RefCell<Spring>>>,
    // intermediate values for Runge Kutta
    // stored in the struct to avoid reallocations
    kps: [Vec<P2>; 4],
    kvs: [Vec<V2>; 4],
    kaccs: [Vec<V2>; 4],
    g: V2,
}

impl Default for Contraption {
    fn default() -> Self {
        Self::new()
    }
}

impl Contraption {
    pub fn new() -> Contraption {
        Contraption {
            nodes: Vec::new(),
            next_node_id: 0,
            springs: Vec::new(),
            kps: [const { Vec::new() }; 4],
            kvs: [const { Vec::new() }; 4],
            kaccs: [const { Vec::new() }; 4],
            g: V2::new(0., 0.),
        }
    }

    pub fn rk4(&mut self, delta_t: Float) {
        let factors = [0.5, 0.5, 1.0];
        self.rk4_init();

        // Runge Kutta steps

        for (n, factor) in factors.iter().enumerate() {
            for (i, node) in self.nodes.iter().enumerate() {
                if self.nodes[i].stationary {
                    self.kps[n + 1][i] = node.position();
                    continue;
                }
                self.kps[n + 1][i] = node.position + factor * delta_t * self.kvs[n][i];
                self.kvs[n + 1][i] = node.velocity + factor * delta_t * self.kaccs[n][i];
            }

            for (i, _) in self.nodes.iter().enumerate() {
                self.kaccs[n + 1][i] =
                    self.g + self.acceleration(i, &self.kps[n + 1], &self.kvs[n + 1]);
            }
        }

        // updating nodes

        for i in 0..self.nodes.len() {
            if self.nodes[i].stationary {
                continue;
            }
            // calculate position
            let weighted_v =
                (self.kvs[0][i] + self.kvs[1][i] * 2.0 + self.kvs[2][i] * 2.0 + self.kvs[3][i])
                    * (1.0 / 6.0);
            self.nodes[i].position += delta_t * weighted_v;

            // calculate velocity
            let weighted_acc = (self.kaccs[0][i]
                + self.kaccs[1][i] * 2.0
                + self.kaccs[2][i] * 2.0
                + self.kaccs[3][i])
                * (1.0 / 6.0);
            self.nodes[i].velocity += delta_t * weighted_acc;

            // put all the nodes positions into kps[0]
            self.kps[0][i] = self.nodes[i].position();
            // put all the nodes velocities into kvs[0]
            self.kvs[0][i] = self.nodes[i].velocity;
        }

        // calculate acceleration
        for i in 0..self.nodes.len() {
            // use the positions and velocities of self.nodes which where collected in the previous loop
            self.nodes[i].acceleration = self.g + self.acceleration(i, &self.kps[0], &self.kvs[0]);
        }
    }

    pub fn acceleration(&self, node_ix: usize, kps: &[P2], kvs: &[V2]) -> V2 {
        let mut acc = V2::ZERO;
        let node = &self.nodes[node_ix];
        for attached in node.attached.iter() {
            if let Some(spring) = attached.spring.upgrade() {
                if let Some(other_ix) = self.get_node_ix(attached.other_node_id) {
                    let dir = kps[other_ix] - kps[node_ix];
                    let vdiff = kvs[other_ix] - kvs[node_ix];
                    acc += spring.borrow().force(&dir, &vdiff) / self.nodes[node_ix].weight;
                }
            }
        }
        acc
    }

    // initialize kps, kvs, kaccs
    fn rk4_init(&mut self) {
        // make shure k values have enough capacity
        for ps in self.kps.iter_mut() {
            while ps.len() < self.nodes.len() {
                ps.push(P2::new(0.0, 0.0));
            }
        }
        for vs in self.kvs.iter_mut() {
            while vs.len() < self.nodes.len() {
                vs.push(V2::ZERO);
            }
        }
        for accs in self.kaccs.iter_mut() {
            while accs.len() < self.nodes.len() {
                accs.push(V2::ZERO);
            }
        }
        // fill in initial values
        for (i, node) in self.nodes.iter().enumerate() {
            self.kps[0][i] = node.position();
            self.kvs[0][i] = node.velocity;
            self.kaccs[0][i] = node.acceleration;
        }
    }

    pub fn step_motors(&mut self, delta_t: Float) {
        for node in self.nodes.iter_mut() {
            if node.stationary {
                if let Some(motor) = node.motor.as_mut() {
                    motor.step(delta_t);
                }
            }
        }
    }

    pub fn get_node_ix(&self, id: NodeId) -> Option<usize> {
        if self.nodes.is_empty() {
            return None;
        }
        let min_ix = 0;
        let max_ix = id.min(self.nodes.len() - 1);
        // println!("max_ix: {}, nodes.len(): {}", max_ix, self.nodes.len());
        self.get_node_ix_min_max(id, min_ix, max_ix)
    }

    fn get_node_ix_min_max(&self, id: NodeId, min_ix: usize, max_ix: usize) -> Option<usize> {
        if min_ix == max_ix {
            if self.nodes[min_ix].id == id {
                return Some(min_ix);
            } else {
                return None;
            }
        }
        let hix = (max_ix + min_ix) / 2;
        match self.nodes[hix].id.cmp(&id) {
            Ordering::Less => self.get_node_ix_min_max(id, hix + 1, max_ix),
            Ordering::Equal => Some(hix),
            Ordering::Greater => self.get_node_ix_min_max(id, min_ix, hix),
        }
    }

    pub fn insert_spring(
        &mut self,
        node_id_a: NodeId,
        node_id_b: NodeId,
        spring_constant: Float,
        dampening: Float,
    ) {
        let ixa = self.get_node_ix(node_id_a);
        let ixb = self.get_node_ix(node_id_b);
        if let (Some(ixa), Some(ixb)) = (ixa, ixb) {
            let node_distance = self.nodes[ixa].position.distance(self.nodes[ixb].position);
            self.create_spring(node_distance, spring_constant, dampening);
            let aa = Attached {
                spring: Rc::downgrade(self.springs.last().unwrap()),
                other_node_id: self.nodes[ixb].id,
            };
            self.nodes[ixa].attached.push(aa);
            let ab = Attached {
                spring: Rc::downgrade(self.springs.last().unwrap()),
                other_node_id: self.nodes[ixa].id,
            };
            self.nodes[ixb].attached.push(ab);
        }
    }

    pub fn insert_node(&mut self, position: P2, radius: Float) {
        let node = Node {
            id: self.next_node_id,
            radius,
            position,
            velocity: V2::ZERO,
            acceleration: V2::ZERO,
            weight: radius * radius,
            stationary: false,
            motor: None,
            attached: Vec::new(),
        };
        self.next_node_id += 1;
        self.nodes.push(node);
    }

    pub fn spring_exists_between(&self, node_id_a: NodeId, node_id_b: NodeId) -> bool {
        if let Some(ixa) = self.get_node_ix(node_id_a) {
            for attached in self.nodes[ixa].attached.iter() {
                if attached.other_node_id == node_id_b {
                    return true;
                }
            }
        }
        if let Some(ixb) = self.get_node_ix(node_id_b) {
            for attached in self.nodes[ixb].attached.iter() {
                if attached.other_node_id == node_id_a {
                    return true;
                }
            }
        }
        false
    }

    pub fn insert_spring_on_indices(
        &mut self,
        ixa: usize,
        ixb: usize,
        spring_constant: Float,
        spring_dampening: Float,
    ) {
        let node_distance = self.nodes[ixa].position.distance(self.nodes[ixb].position);
        self.create_spring(node_distance, spring_constant, spring_dampening);
        let ab = Attached {
            spring: Rc::downgrade(self.springs.last().unwrap()),
            other_node_id: self.nodes[ixb].id,
        };
        self.nodes[ixa].attached.push(ab);
        let ba = Attached {
            spring: Rc::downgrade(self.springs.last().unwrap()),
            other_node_id: self.nodes[ixa].id,
        };
        self.nodes[ixb].attached.push(ba);
    }

    pub fn fully_connect(&mut self, node_ids: &[NodeId], spring_constant: Float, dampening: Float) {
        for ida in node_ids.iter().copied() {
            for idb in node_ids.iter().copied() {
                if ida != idb && !self.spring_exists_between(ida, idb) {
                    self.insert_spring(ida, idb, spring_constant, dampening);
                }
            }
        }
    }

    pub fn apply_resting_length_factors(&self) {
        for spring in self.springs.iter() {
            spring.borrow_mut().apply_resting_length_factor();
        }
    }

    pub fn add_grid(
        &mut self,
        GridConfig {
            ty,
            node_distance,
            node_radius,
            spring_constant,
            spring_dampening,
            width,
            mut height,
            upper_left,
        }: GridConfig,
    ) {
        let mut nodes = Vec::new();
        // create all the nodes
        let nd2 = node_distance.powi(2);
        let tri_height = (nd2 - nd2 / 4.0).sqrt();
        if ty == GridType::Hex {
            height *= 2;
            height += 1;
        }
        for y in 0..height {
            for x in 0..width {
                let pos = match ty {
                    GridType::Quad | GridType::QuadTriangle | GridType::Cross => {
                        upper_left + node_distance * V2::new(x as Float, y as Float)
                    }
                    GridType::Triangle => {
                        if y % 2 == 0 {
                            upper_left
                                + V2::new(x as Float * node_distance, y as Float * tri_height)
                        } else {
                            upper_left
                                + V2::new(
                                    x as Float * node_distance - (node_distance / 2.0),
                                    y as Float * tri_height,
                                )
                        }
                    }
                    GridType::Hex => {
                        if (y % 2) == 0 {
                            upper_left
                                + V2::new(x as Float * node_distance * 3., y as Float * tri_height)
                        } else {
                            upper_left
                                + V2::new(
                                    x as Float * node_distance * 3. - (node_distance / 2.0),
                                    y as Float * tri_height,
                                )
                        }
                    }
                };
                nodes.push(Node::new(self.next_node_id, pos, node_radius));
                self.next_node_id += 1;
                if ty == GridType::Hex {
                    let pos2 = if (y % 2) == 0 {
                        P2::new(pos.x + node_distance, pos.y)
                    } else {
                        P2::new(pos.x + 2. * node_distance, pos.y)
                    };
                    nodes.push(Node::new(self.next_node_id, pos2, node_radius));
                    self.next_node_id += 1;
                }
            }
        }
        dbg!(nodes.len());
        let base_ix = self.nodes.len();
        self.nodes.extend_from_slice(&nodes);
        dbg!(width, height);
        // link the nodes
        for y in 0..height {
            for x in 0..width {
                // Hex grids are fully linked in the match statement bellow
                if ty != GridType::Hex {
                    // link horizontally
                    if x < (width - 1) {
                        let ixa = base_ix + x + y * width;
                        let ixb = ixa + 1;
                        self.insert_spring_on_indices(ixa, ixb, spring_constant, spring_dampening);
                    }
                    // link vertically
                    if y < (height - 1) {
                        let ixa = base_ix + x + y * width;
                        let ixb = ixa + width;
                        self.insert_spring_on_indices(ixa, ixb, spring_constant, spring_dampening);
                    }
                }
                match ty {
                    GridType::Quad => {}
                    GridType::QuadTriangle => {
                        if (x < (width - 1)) && (y < (height - 1)) {
                            let ixa = base_ix + x + y * width;
                            let ixb = ixa + 1 + width;
                            self.insert_spring_on_indices(
                                ixa,
                                ixb,
                                spring_constant,
                                spring_dampening,
                            );
                        }
                    }
                    GridType::Triangle => {
                        if (x < (width - 1)) && (y < (height - 1)) {
                            let ixa = base_ix + x + y * width;
                            let ixb = if (y % 2) == 0 {
                                ixa + 1 + width
                            } else if x == 0 {
                                continue;
                            } else {
                                ixa - 1 + width
                            };
                            self.insert_spring_on_indices(
                                ixa,
                                ixb,
                                spring_constant,
                                spring_dampening,
                            );
                        }
                    }
                    GridType::Cross => {
                        if (x < (width - 1)) && (y < (height - 1)) {
                            let ixa = base_ix + x + y * width;
                            let ixb = ixa + 1 + width;
                            self.insert_spring_on_indices(
                                ixa,
                                ixb,
                                spring_constant,
                                spring_dampening,
                            );
                        }
                        if (x > 0) && (y < (height - 1)) {
                            let ixa = base_ix + x + y * width;
                            let ixb = ixa - 1 + width;
                            self.insert_spring_on_indices(
                                ixa,
                                ixb,
                                spring_constant,
                                spring_dampening,
                            );
                        }
                    }
                    GridType::Hex => {
                        // left to right
                        let ixa1 = base_ix + 2 * x + y * 2 * width + (y % 2);
                        let ixb1 = ixa1 + 1;
                        // make sure the x coordinate does not overshoot
                        if 2 * x + 1 + (y % 2) < 2 * width {
                            self.insert_spring_on_indices(
                                ixa1,
                                ixb1,
                                spring_constant,
                                spring_dampening,
                            );
                            //diagonally down right
                            if y + 1 < height {
                                let ixa3 = ixb1;
                                let ixb3 = ixa3 + 2 * width;
                                self.insert_spring_on_indices(
                                    ixa3,
                                    ixb3,
                                    spring_constant,
                                    spring_dampening,
                                );
                            }
                        }
                        if y + 1 < height {
                            // diagonally down left
                            let ixa2 = ixa1;
                            let ixb2 = ixa2 + 2 * width;
                            self.insert_spring_on_indices(
                                ixa2,
                                ixb2,
                                spring_constant,
                                spring_dampening,
                            );
                        }
                        if x == 0 && (y % 2) == 1 && (y + 1 < height) {
                            let ixa = y * 2 * width;
                            let ixb = ixa + 2 * width;
                            self.insert_spring_on_indices(
                                ixa,
                                ixb,
                                spring_constant,
                                spring_dampening,
                            );
                        }
                    }
                }
            }
        }
    }

    pub fn bounds_collision(&mut self, width: Float, height: Float, elasticity: Float) {
        for node in self.nodes.iter_mut() {
            node.bounds_collision(width, height, elasticity);
        }
    }

    pub fn nearest_node_ix(&self, pos: &P2) -> Option<(usize, Float)> {
        let mut res = None;
        for (ix, node) in self.nodes.iter().enumerate() {
            let new_dist = pos.distance_squared(node.position());
            if let Some((_, dist)) = res {
                if new_dist < dist {
                    res = Some((ix, new_dist));
                }
            } else {
                res = Some((ix, new_dist));
            }
        }
        res
    }

    /// returns a reference to the spring and its distance
    pub fn nearest_spring(&self, to: &P2) -> Option<(Weak<RefCell<Spring>>, Float)> {
        let mut res = None;
        for node in self.nodes.iter() {
            for attached in node.attached.iter() {
                if let Some(ix) = self.get_node_ix(attached.other_node_id) {
                    let other_node = &self.nodes[ix];
                    let dirab = other_node.position() - node.position();
                    if dirab.length() <= f64::EPSILON {
                        continue;
                    }
                    let dirba = -dirab;
                    // if to is closer to the spring than any of its nodes
                    let a = to - node.position();
                    let b = to - other_node.position();
                    if (dirab.dot(a) > 0.0) && (dirba.dot(b) > 0.0) {
                        // calculate the distance
                        let dirabn = dirab.normalize();
                        let dist = (a.length_squared() - dirabn.dot(a).powi(2)).sqrt();
                        if let Some((_, current_dist)) = res {
                            if dist < current_dist {
                                res = Some((attached.spring.clone(), dist));
                            }
                        } else {
                            res = Some((attached.spring.clone(), dist));
                        }
                    }
                }
            }
        }
        res
    }

    pub fn multi_select(
        &self,
        rect: &Rect,
        previous_nodes: &[NodeId],
    ) -> (Vec<NodeId>, Vec<Weak<RefCell<Spring>>>) {
        let mut node_ids = Vec::from_iter(previous_nodes.iter().copied());
        let mut node_indices = Vec::new();
        for id in node_ids.iter() {
            if let Some(ix) = self.get_node_ix(*id) {
                node_indices.push(ix);
            }
        }
        for (ix, node) in self.nodes.iter().enumerate() {
            let pos =
                macroquad::math::Vec2::new(node.position().x as f32, node.position().y as f32);
            if rect.contains(pos) {
                node_ids.push(node.id);
                node_indices.push(ix);
            }
        }
        node_ids.sort();
        node_ids.dedup();
        node_indices.sort();
        node_indices.dedup();
        let mut springs = Vec::new();
        for node_ix in node_indices.iter() {
            for attached in self.nodes[*node_ix].attached.iter() {
                if node_ids.binary_search(&attached.other_node_id).is_ok() {
                    springs.push(attached.spring.clone());
                }
            }
        }
        (node_ids, springs)
    }

    pub fn create_spring(
        &mut self,
        resting_length: Float,
        spring_constant: Float,
        damping_factor: Float,
    ) {
        let spring = Spring {
            spring_constant,
            resting_length,
            resting_length_factor: 1.0,
            damping_factor,
            min_len: 0.2 * resting_length,
            max_len: 5.0 * resting_length,
        };
        self.springs.push(Rc::new(RefCell::new(spring)));
    }

    pub fn create_spring_between(
        &mut self,
        ixa: usize,
        ixb: usize,
        spring_constant: Float,
        damping_factor: Float,
    ) {
        let resting_length = self.nodes[ixa].position.distance(self.nodes[ixb].position);
        let spring = Spring {
            spring_constant,
            resting_length,
            resting_length_factor: 1.0,
            damping_factor,
            min_len: 0.25 * resting_length,
            max_len: 4.0 * resting_length,
        };
        let spring_ref = Rc::new(RefCell::new(spring));
        self.springs.push(spring_ref.clone());
        let ata = Attached {
            spring: Rc::downgrade(&spring_ref),
            other_node_id: self.nodes[ixb].id,
        };
        self.nodes[ixa].attached.push(ata);
        let atb = Attached {
            spring: Rc::downgrade(&spring_ref),
            other_node_id: self.nodes[ixa].id,
        };
        self.nodes[ixb].attached.push(atb);
    }

    pub fn spring_break(&mut self) {
        let mut delete_indices = Vec::new();
        for (node_ix, node) in self.nodes.iter().enumerate() {
            for (attached_ix, attached) in node.attached.iter().enumerate() {
                if attached.is_destruct(&self.nodes, node.position()) {
                    delete_indices.push((node_ix, attached_ix));
                }
            }
        }
        // iterate in reverse, so the indices remain valid after one or more Attached's have been
        // removed
        for (node_ix, attached_ix) in delete_indices.iter().rev().copied() {
            self.nodes[node_ix].attached.remove(attached_ix);
        }
    }

    pub fn delete_node(&mut self, node_id: NodeId) {
        if let Some(ix) = self.get_node_ix(node_id) {
            self.nodes.remove(ix);
        }
    }

    pub fn delete_spring(&mut self, spring_rc: Weak<RefCell<Spring>>) {
        let mut remove_ix = None;
        for (ix, spring) in self.springs.iter().enumerate() {
            let downgraded = Rc::downgrade(spring);
            if Weak::ptr_eq(&downgraded, &spring_rc) {
                remove_ix = Some(ix);
                break;
            }
        }
        if let Some(ix) = remove_ix {
            self.springs.remove(ix);
        }
    }

    pub fn clear(&mut self) {
        self.springs.clear();
        self.nodes.clear();
    }

    pub fn energy(&self) -> Float {
        let mut energy = 0.0;
        for node in &self.nodes {
            energy += 0.5 * node.weight * node.velocity.length_squared();
            for attached in &node.attached {
                if let Some(spring) = attached.spring.upgrade() {
                    if let Some(other_node_ix) = self.get_node_ix(attached.other_node_id) {
                        let distance = node
                            .position()
                            .distance(self.nodes[other_node_ix].position());
                        energy += 0.25
                            * (spring.borrow().resting_length - distance).powi(2)
                            * spring.borrow().spring_constant;
                    }
                }
            }
        }
        energy
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Spring {
    spring_constant: Float,
    resting_length: Float,
    resting_length_factor: Float,
    damping_factor: Float,
    min_len: Float,
    max_len: Float,
}

impl Spring {
    fn force(&self, dir: &V2, vdiff: &V2) -> V2 {
        if *dir == V2::ZERO {
            V2::ZERO
        } else {
            let diff =
                (1.0 - ((self.resting_length * self.resting_length_factor) / dir.length())) * dir;
            self.spring_constant * diff + self.damping_factor * vdiff.project_onto(*dir)
        }
    }

    pub fn apply_resting_length_factor(&mut self) {
        self.resting_length *= self.resting_length_factor;
        self.resting_length_factor = 1.0;
    }

    pub fn ui(&mut self, ui: &mut Ui) -> bool {
        let mut changed = false;
        ui.horizontal(|ui| {
            ui.label("spring constant");
            changed |= ui
                .add(Slider::new(&mut self.spring_constant, 0.0..=100000.0).logarithmic(true))
                .changed();
        });
        ui.horizontal(|ui| {
            ui.label("resting length factor");
            changed |= ui
                .add(Slider::new(&mut self.resting_length_factor, 0.0..=2.0))
                .changed();
        });
        ui.horizontal(|ui| {
            ui.label("min len");
            changed |= ui
                .add(Slider::new(
                    &mut self.min_len,
                    0.0..=(self.resting_length * self.resting_length_factor),
                ))
                .changed();
        });
        ui.horizontal(|ui| {
            ui.label("max len");
            changed |= ui
                .add(Slider::new(
                    &mut self.max_len,
                    (self.resting_length * self.resting_length_factor)
                        ..=(10.0 * self.resting_length * self.resting_length_factor),
                ))
                .changed();
        });
        ui.horizontal(|ui| {
            ui.label("damping factor");
            changed |= ui
                .add(Slider::new(&mut self.damping_factor, 0.0..=100.0).logarithmic(true))
                .changed();
        });
        changed
    }
}

#[derive(Debug, Clone)]
pub struct Attached {
    spring: Weak<RefCell<Spring>>,
    other_node_id: NodeId,
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
            let spring = spring_ref.borrow();
            !(spring.min_len..=spring.max_len).contains(&self.length(nodes, main_node_position))
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
    id: NodeId,
    radius: Float,
    position: P2,
    velocity: V2,
    acceleration: V2,
    weight: Float,
    stationary: bool,
    motor: Option<Motor>,
    attached: Vec<Attached>,
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
            attached: Vec::new(),
        }
    }

    pub fn bounds_collision(&mut self, width: Float, height: Float, elasticity: Float) {
        match (
            self.position()[0] < self.radius,
            self.position()[0] > (width - self.radius),
        ) {
            (true, _) => {
                if self.velocity[0] < 0. {
                    self.velocity[0] = -self.velocity[0];
                    self.velocity *= elasticity;
                }
                self.position()[0] = self.radius;
            }
            (_, true) => {
                if self.velocity[0] > 0. {
                    self.velocity[0] = -self.velocity[0];
                    self.velocity *= elasticity;
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
                    self.velocity *= elasticity;
                }
                self.position()[1] = self.radius;
            }
            (_, true) => {
                if self.velocity[1] > 0. {
                    self.velocity[1] = -self.velocity[1];
                    self.velocity *= elasticity;
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

pub struct SpringyApp {
    app_mode: AppMode,
    contraption: Contraption,
    gravity_on: bool,
    grid: Grid,
    draw_velocities: bool,
    velocity_arrow_factor: Float,
    draw_accelerations: bool,
    acceleration_arrow_factor: Float,
    multi_selection_rect: Option<(Rect, Vec2)>,
    gravity: V2,
    node_config: NodeConfig,
    grid_config: GridConfig,
    create_spring_constant: Float,
    create_spring_dampening: Float,
}

impl Default for SpringyApp {
    fn default() -> Self {
        Self::new()
    }
}

impl SpringyApp {
    pub fn new() -> Self {
        SpringyApp {
            app_mode: AppMode::NodeSelected { node_id: None },
            contraption: Contraption::new(),
            gravity_on: false,
            grid: Grid::default(),
            draw_velocities: true,
            velocity_arrow_factor: 0.1,
            draw_accelerations: true,
            acceleration_arrow_factor: 0.01,
            multi_selection_rect: None,
            gravity: V2::new(0., 1000.),
            grid_config: GridConfig {
                ty: GridType::Cross,
                node_distance: 30.0,
                node_radius: 4.0,
                spring_constant: 10000.0,
                spring_dampening: 4.0,
                width: 15,
                height: 15,
                upper_left: P2::new(200., 200.),
            },
            node_config: Default::default(),
            create_spring_constant: 10.0,
            create_spring_dampening: 10.0,
        }
    }

    pub fn app_mode_ui(&mut self, ui: &mut Ui) {
        let mut app_mode_changed = false;
        if ui.button("Select").clicked() {
            self.app_mode = AppMode::NodeSelected { node_id: None };
            self.contraption.apply_resting_length_factors();
            app_mode_changed = true;
        }
        if ui.button("MultiSelection").clicked() {
            self.app_mode = AppMode::MultiSelection {
                node_ids: Vec::new(),
                springs: Vec::new(),
            };
            app_mode_changed = true;
        }
        if ui.button("Draw Node").clicked() {
            self.app_mode = AppMode::PlaceNode;
            app_mode_changed = true;
        }
        if ui.button("Draw Joints").clicked() {
            self.app_mode = AppMode::SpringStart;
            app_mode_changed = true;
        }
        if ui.button("Cut Spring").clicked() {
            self.app_mode = AppMode::CutSpring;
            app_mode_changed = true;
        }
        if ui.button("Clear").clicked() {
            self.contraption.clear();
        }
        if app_mode_changed {
            self.contraption.apply_resting_length_factors();
        }
    }

    /// select the nearest node or spring and conveniently return the coresponding AppMode
    pub fn select_nearest(&self) -> AppMode {
        let mouse_pos = self.grid.mouse_position();
        let mouse_pos = P2::new(mouse_pos.0 as Float, mouse_pos.1 as Float);
        let nearest_node = self.contraption.nearest_node_ix(&mouse_pos);
        let nearest_spring = self.contraption.nearest_spring(&mouse_pos);
        match (nearest_node, nearest_spring) {
            (None, None) => AppMode::NodeSelected { node_id: None },
            (None, Some((spring, _))) => AppMode::SpringSelected {
                spring: spring.clone(),
            },
            (Some((node_ix, _)), None) => AppMode::NodeSelected {
                node_id: Some(self.contraption.nodes[node_ix].id),
            },
            (Some((node_ix, node_dist)), Some((spring, spring_dist))) => {
                if node_dist < spring_dist {
                    AppMode::NodeSelected {
                        node_id: Some(self.contraption.nodes[node_ix].id),
                    }
                } else {
                    AppMode::SpringSelected {
                        spring: spring.clone(),
                    }
                }
            }
        }
    }

    pub fn node_ui(&mut self, ui: &mut Ui) {
        let radius_range = 0.1..=10.0;
        match self.app_mode {
            AppMode::NodeSelected { node_id: Some(id) } => {
                if let Some(ix) = self.contraption.get_node_ix(id) {
                    let node = &mut self.contraption.nodes[ix];
                    if ui.checkbox(&mut node.stationary, "stationary").changed() {
                        node.velocity = V2::ZERO;
                        node.acceleration = V2::ZERO;
                    }
                    ui.label(format!("{node}"));
                    ui.horizontal(|ui| {
                        ui.label("node radius");
                        if ui
                            .add(Slider::new(&mut node.radius, radius_range))
                            .changed()
                        {
                            node.weight = node.radius.powi(2);
                        }
                    });
                }
            }
            AppMode::PlaceNode => {
                ui.horizontal(|ui| {
                    ui.label("node radius");
                    ui.add(Slider::new(&mut self.node_config.node_radius, radius_range));
                });
            }
            _ => {}
        }
    }

    pub fn spring_ui(&mut self, ui: &mut Ui) {
        let spring_constant_range = 0.0..=10000.0;
        let dampening_range = 0.0..=100.0;
        match self.app_mode.clone() {
            AppMode::SpringSelected { spring } => {
                if let Some(spring_cell) = spring.upgrade() {
                    let mut spring_ref = spring_cell.borrow_mut();
                    ui.horizontal(|ui| {
                        ui.label("spring constant");
                        ui.add(
                            Slider::new(&mut spring_ref.spring_constant, spring_constant_range)
                                .logarithmic(true),
                        );
                    });
                    ui.horizontal(|ui| {
                        ui.label("spring dampening");
                        ui.add(
                            Slider::new(&mut spring_ref.damping_factor, dampening_range)
                                .logarithmic(true),
                        );
                    });
                }
            }
            AppMode::SpringStart | AppMode::SpringEnd { .. } => {
                ui.horizontal(|ui| {
                    ui.label("spring constant");
                    ui.add(Slider::new(
                        &mut self.create_spring_constant,
                        spring_constant_range,
                    ));
                });
                ui.horizontal(|ui| {
                    ui.label("spring dampening");
                    ui.add(
                        Slider::new(&mut self.create_spring_dampening, dampening_range)
                            .logarithmic(true),
                    );
                });
            }
            _ => {}
        };
    }

    pub fn arrows_ui(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.draw_velocities, "Draw Velocities");
            if self.draw_velocities {
                ui.add(
                    Slider::new(&mut self.velocity_arrow_factor, 0.001..=10.0).logarithmic(true),
                );
            }
        });
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.draw_accelerations, "Draw Accelerations");
            if self.draw_accelerations {
                ui.add(
                    Slider::new(&mut self.acceleration_arrow_factor, 0.001..=10.0)
                        .logarithmic(true),
                );
            }
        });
    }

    pub fn multi_selection_ui(&mut self, ui: &mut Ui) {
        if let AppMode::MultiSelection { node_ids, springs } = &self.app_mode {
            let mut node_ids_iter = node_ids.iter();
            let mut current_id = node_ids_iter.next();
            if self.node_config.ui(ui) {
                for node in self.contraption.nodes.iter_mut() {
                    if Some(&node.id) == current_id {
                        current_id = node_ids_iter.next();
                        node.configure(self.node_config);
                    }
                }
            }
            let master_spring = springs.first();
            if let Some(spring_ref) = master_spring {
                if let Some(spring_rc) = spring_ref.upgrade() {
                    if spring_rc.borrow_mut().ui(ui) {
                        let mut copy_spring = *spring_rc.borrow();
                        for weak_spring in springs.iter() {
                            if let Some(spring_rc) = weak_spring.upgrade() {
                                // resting length should not be modified by the UI
                                // resting_length_factor should be instead
                                copy_spring.resting_length = spring_rc.borrow().resting_length;
                                spring_rc.replace(copy_spring);
                            }
                        }
                    }
                }
            }
            if ui.button("fully connect").clicked() {
                self.contraption.fully_connect(node_ids, 1000.0, 4.0);
            }
        }
    }

    pub fn draw_contraption(&self) {
        for (node_ix, node) in self.contraption.nodes.iter().enumerate() {
            let mut nearest_spring = None;
            let mut nearest_node = None;
            let mpos = self.grid.mouse_position();
            let mpos = DVec2::new(mpos.0 as Float, mpos.1 as Float);
            if let AppMode::NodeSelected { node_id: None } = &self.app_mode {
                match (
                    self.contraption.nearest_spring(&mpos),
                    self.contraption.nearest_node_ix(&mpos),
                ) {
                    (None, None) => {}
                    (None, Some((node_ix, _))) => nearest_node = Some(node_ix),
                    (Some((spring, _)), None) => nearest_spring = Some(spring),
                    (Some((spring, spring_dist)), Some((node_ix, node_dist))) => {
                        if spring_dist < node_dist {
                            nearest_spring = Some(spring);
                        } else {
                            nearest_node = Some(node_ix);
                        }
                    }
                }
            }
            let pos = node.position();
            let mut node_color = GREEN;
            if Some(node_ix) == nearest_node {
                node_color = PURPLE;
            } else if node.stationary {
                node_color = RED;
            }
            draw_circle(pos.x as f32, pos.y as f32, node.radius as f32, node_color);
            if let AppMode::MultiSelection {
                node_ids,
                springs: _,
            } = &self.app_mode
            {
                if node_ids.binary_search(&node.id).is_ok() {
                    let mut color = SKYBLUE;
                    color.a = 0.5;
                    draw_circle(pos.x as f32, pos.y as f32, 2.0 * node.radius as f32, color);
                }
            }
            let pos_a = node.position();
            for attached in &node.attached {
                if let Some(other_ix) = self.contraption.get_node_ix(attached.other_node_id) {
                    let pos_b = self.contraption.nodes[other_ix].position();
                    if let Some(spring) = attached.spring.upgrade() {
                        let dist = pos_a.distance(pos_b);
                        let resting_length = spring.borrow().resting_length;
                        let max_len = spring.borrow().max_len;
                        let stretch_factor = if dist < resting_length {
                            (resting_length - dist) / resting_length
                        } else {
                            (dist - resting_length) / (max_len - resting_length)
                        };
                        let color = colorous::RED_YELLOW_GREEN
                            .eval_continuous((1.0 - stretch_factor).powi(4));
                        let mut spring_color =
                            macroquad::color::Color::from_rgba(color.r, color.g, color.b, 255);
                        if let Some(nearest) = nearest_spring.as_ref() {
                            if nearest.as_ptr() == attached.spring.as_ptr() {
                                spring_color = PURPLE;
                            }
                        }
                        draw_line(
                            pos_a[0] as f32,
                            pos_a[1] as f32,
                            pos_b[0] as f32,
                            pos_b[1] as f32,
                            2.,
                            spring_color,
                        );
                    }
                }
            }
            if self.draw_velocities {
                let mut velocity_color = BEIGE;
                velocity_color.a = 0.5;
                draw_vector(
                    &node.position(),
                    &(node.velocity * self.velocity_arrow_factor),
                    velocity_color,
                );
            }
            if self.draw_accelerations {
                let mut acceleration_color = VIOLET;
                acceleration_color.a = 0.5;
                draw_vector(
                    &node.position(),
                    &(node.acceleration * self.acceleration_arrow_factor),
                    acceleration_color,
                );
            }
        }
    }

    fn pull_to_mouse(&mut self, force: Float) {
        let mpos = self.grid.mouse_position();
        let mpos = P2::new(mpos.0 as Float, mpos.1 as Float);
        for node in self.contraption.nodes.iter_mut() {
            if node.stationary {
                continue;
            }
            let dir = mpos - node.position;
            let dist = dir.length();
            let dir = dir.normalize();
            if dist > 0. {
                node.velocity += dir * (force * 10. / (dist + 10.));
            }
        }
    }
}

#[macroquad::main("Springy")]
async fn main() {
    let mut app = SpringyApp::new();
    let mut mpos;
    let mut mouse_hover = false;
    let mut start = get_time();
    let mut stop;
    let mut delta_t = 0.001;
    let mut step_elapsed;
    let mut paused: bool = false;
    loop {
        let width = screen_width() as Float;
        let height = screen_height() as Float;
        mpos = app.grid.mouse_position();
        clear_background(BLACK);
        egui_macroquad::ui(|ctx| {
            Window::new("Springy").show(ctx, |ui| {
                ui.label(format!("Mouse: ({}, {})", mpos.0, mpos.1));
                ui.label(format!("num nodes: {}", app.contraption.nodes.len()));
                ui.label(format!("energy: {}", app.contraption.energy()));
                if paused {
                    if ui.button("play").clicked() {
                        paused = false;
                    }
                } else if ui.button("pause").clicked() {
                    paused = true;
                }
                disp_app_mode(ui, &app.app_mode);
                match app.app_mode {
                    AppMode::SpringStart { .. } | AppMode::SpringEnd { .. } => {}
                    _ => {}
                }
                ui.horizontal(|ui| {
                    ui.label("delta_t");
                    ui.add(Slider::new(&mut delta_t, 0.0001..=0.1).logarithmic(true));
                });
                ui.checkbox(&mut app.gravity_on, "gravity");
                app.app_mode_ui(ui);
                app.arrows_ui(ui);
                app.grid.ui(ui);
                app.spring_ui(ui);
                app.node_ui(ui);
                app.multi_selection_ui(ui);
                app.grid_config.ui(ui, &mut app.contraption);
                ui.label("spring velocities:");
            });
            mouse_hover = ctx.is_pointer_over_area();
        });

        // MultiSelection: store the mouse position where the selection starts
        if !mouse_hover && is_mouse_button_pressed(MouseButton::Left) {
            if let AppMode::MultiSelection { .. } = app.app_mode {
                let mouse_pos = app.grid.mouse_position();
                app.multi_selection_rect = Some((
                    macroquad::math::Rect::new(mouse_pos.0, mouse_pos.1, 0.0, 0.0),
                    Vec2::new(mouse_pos.0, mouse_pos.1),
                ));
            }
        }

        // draw multiselection rect
        if !mouse_hover && is_mouse_button_down(MouseButton::Left) {
            if let AppMode::MultiSelection { .. } = app.app_mode {
                if let Some((rect, start_pos)) = app.multi_selection_rect.as_mut() {
                    let mut selection_rect_color = macroquad::color::SKYBLUE;
                    selection_rect_color.a = 0.5;
                    let mouse_pos = app.grid.mouse_position();
                    rect.w = (start_pos.x - mouse_pos.0).abs();
                    rect.h = (start_pos.y - mouse_pos.1).abs();
                    rect.x = start_pos.x.min(mouse_pos.0);
                    rect.y = start_pos.y.min(mouse_pos.1);
                    draw_rectangle(rect.x, rect.y, rect.w, rect.h, selection_rect_color);
                }
            }
        }

        if !mouse_hover && is_mouse_button_released(MouseButton::Left) {
            app.contraption.apply_resting_length_factors();
            match &mut app.app_mode {
                AppMode::NodeSelected { node_id: _ } | AppMode::SpringSelected { spring: _ } => {
                    app.app_mode = app.select_nearest();
                }
                AppMode::MultiSelection { node_ids, springs } => {
                    if let Some((rect, _)) = app.multi_selection_rect {
                        if is_key_down(KeyCode::LeftControl) || is_key_down(KeyCode::RightControl) {
                            (*node_ids, *springs) = app.contraption.multi_select(&rect, node_ids);
                            springs.sort_by_key(|weak_spring| weak_spring.as_ptr());
                            springs.dedup_by_key(|weak_spring| weak_spring.as_ptr());
                        } else {
                            (*node_ids, *springs) = app.contraption.multi_select(&rect, &[]);
                        }
                        if let Some(id) = node_ids.first() {
                            if let Some(ix) = app.contraption.get_node_ix(*id) {
                                app.node_config = NodeConfig::from_node(&app.contraption.nodes[ix]);
                            }
                        }
                    }
                }
                AppMode::PlaceNode => {
                    app.contraption.insert_node(
                        P2::new(mpos.0 as Float, mpos.1 as Float),
                        app.node_config.node_radius,
                    );
                }
                AppMode::SpringStart => {
                    let mouse_position = P2::new(mpos.0 as Float, mpos.1 as Float);
                    if let Some((ixa, _)) = app.contraption.nearest_node_ix(&mouse_position) {
                        app.app_mode = AppMode::SpringEnd {
                            node_a_id: app.contraption.nodes[ixa].id,
                        };
                    }
                }
                AppMode::SpringEnd { node_a_id } => {
                    let mouse_position = P2::new(mpos.0 as Float, mpos.1 as Float);
                    if let Some(ixa) = app.contraption.get_node_ix(*node_a_id) {
                        if let Some((ixb, _)) = app.contraption.nearest_node_ix(&mouse_position) {
                            app.contraption.create_spring_between(
                                ixa,
                                ixb,
                                app.create_spring_constant,
                                app.create_spring_dampening,
                            );
                            app.app_mode = AppMode::SpringStart;
                        }
                    }
                }
                AppMode::CutSpring => {}
            }
        }

        if is_key_pressed(KeyCode::P) {
            paused = !paused;
        }

        if paused && is_key_released(KeyCode::S) {
            app.contraption.rk4(delta_t as Float);
            app.contraption.step_motors(delta_t as Float);
        }

        // handle gravity_on
        if app.gravity_on {
            app.contraption.g = app.gravity;
        } else {
            app.contraption.g = V2::new(0., 0.)
        }

        // calculate delta_t
        stop = get_time();
        let elapsed = stop - start;
        start = stop;
        step_elapsed = 0.0;

        if !paused {
            loop {
                // update contraption
                app.contraption.bounds_collision(width, height, 0.9);
                if (step_elapsed + delta_t) < elapsed {
                    app.contraption.rk4(delta_t as Float);
                    step_elapsed += delta_t;
                } else {
                    // final step
                    app.contraption.rk4((elapsed - step_elapsed) as Float);
                    break;
                }
                app.contraption.spring_break();
            }
            app.contraption.step_motors(elapsed);
        }
        if is_mouse_button_down(MouseButton::Right) && !mouse_hover {
            app.pull_to_mouse(500.);
        }
        if is_mouse_button_down(MouseButton::Left) && (app.app_mode == AppMode::CutSpring) {
            let mpos = app.grid.mouse_position();
            if let Some((spring_ref, dist)) = app
                .contraption
                .nearest_spring(&P2::new(mpos.0 as Float, mpos.1 as Float))
            {
                if dist < 2.0 {
                    app.contraption.delete_spring(spring_ref);
                }
            }
        }
        app.draw_contraption();
        // draw the grid if shift is pressed
        if is_key_down(KeyCode::LeftShift) || is_key_down(KeyCode::RightShift) {
            app.grid.draw();
        }
        // draw the frame
        egui_macroquad::draw();
        next_frame().await
    }
}

fn disp_app_mode(ui: &mut Ui, mode: &AppMode) {
    ui.horizontal(|ui| {
        ui.label("App Mode: ");
        match mode {
            AppMode::PlaceNode => ui.label("Draw Nodes"),
            AppMode::SpringSelected { .. } => ui.label("SpringSelected"),
            AppMode::NodeSelected { .. } => ui.label("NodeSelected"),
            AppMode::MultiSelection { .. } => ui.label("MultiSelection"),
            AppMode::SpringStart => ui.label("JointStart"),
            AppMode::SpringEnd { .. } => ui.label("JointEnd"),
            AppMode::CutSpring => ui.label("CutSpring"),
        }
    });
}

fn draw_vector(start_point: &P2, v: &V2, color: Color) {
    let end_point = start_point + v;
    draw_line(
        start_point[0] as f32,
        start_point[1] as f32,
        end_point[0] as f32,
        end_point[1] as f32,
        2.,
        color,
    );
    let s = v[1].signum();
    let mut angle = v.angle_to(V2::new(s, 0.0)).abs();
    if v.y < 0.0 {
        angle += std::f64::consts::PI;
    }
    let angle = 180.0 * angle / std::f64::consts::PI;
    draw_poly(
        end_point[0] as f32,
        end_point[1] as f32,
        3,
        5.,
        angle as f32,
        color,
    );
}

#[derive(Clone, Debug)]
pub enum AppMode {
    MultiSelection {
        node_ids: Vec<NodeId>,
        springs: Vec<Weak<RefCell<Spring>>>,
    },
    NodeSelected {
        node_id: Option<NodeId>,
    },
    SpringSelected {
        spring: Weak<RefCell<Spring>>,
    },
    PlaceNode,
    SpringStart,
    SpringEnd {
        node_a_id: NodeId,
    },
    CutSpring,
}

impl PartialEq for AppMode {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (AppMode::PlaceNode, AppMode::PlaceNode) => true,
            (
                AppMode::NodeSelected { node_id: node_a_id },
                AppMode::NodeSelected { node_id: node_b_id },
            ) => node_a_id == node_b_id,
            (
                AppMode::SpringEnd { node_a_id },
                AppMode::SpringEnd {
                    node_a_id: node_b_id,
                },
            ) => node_a_id == node_b_id,
            (AppMode::SpringStart, AppMode::SpringStart) => true,
            (AppMode::SpringSelected { spring: sa }, AppMode::SpringSelected { spring: sb }) => {
                sa.ptr_eq(sb)
            }
            (AppMode::CutSpring, AppMode::CutSpring) => true,
            _ => false,
        }
    }
}

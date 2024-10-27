use core::f64;
use egui::*;
use macroquad::prelude::*;
use nalgebra::{distance, distance_squared, Point2, UnitVector2, Vector2};
pub(crate) use std::{
    cell::RefCell,
    fmt::Display,
    rc::{Rc, Weak},
};
use std::{cmp::Ordering, ops::Mul};

type Float = f64;
type V2 = Vector2<Float>;
type P2 = Point2<Float>;
type U2 = UnitVector2<Float>;
type NodeId = usize;

pub trait Projectable {
    fn project_on(&self, onto: &V2) -> V2;
}

impl Projectable for V2 {
    fn project_on(&self, onto: &V2) -> V2 {
        let onton = onto.normalize();
        self.dot(&onton) * onton
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum GridType {
    Quad,
    Triangle,
    Cross,
}

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
            self.kps[0][i] = self.nodes[i].position;
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
        let mut acc = V2::zeros();
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
                vs.push(V2::zeros());
            }
        }
        for accs in self.kaccs.iter_mut() {
            while accs.len() < self.nodes.len() {
                accs.push(V2::zeros());
            }
        }
        // fill in initial values
        for (i, node) in self.nodes.iter().enumerate() {
            self.kps[0][i] = node.position;
            self.kvs[0][i] = node.velocity;
            self.kaccs[0][i] = node.acceleration;
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
            let node_distance = distance(&self.nodes[ixa].position, &self.nodes[ixb].position);
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
            velocity: V2::zeros(),
            acceleration: V2::zeros(),
            weight: radius * radius,
            stationary: false,
            attached: Vec::new(),
        };
        self.next_node_id += 1;
        self.nodes.push(node);
    }

    pub fn insert_spring_on_indices(
        &mut self,
        ixa: usize,
        ixb: usize,
        spring_constant: Float,
        spring_dampening: Float,
    ) {
        let node_distance = distance(&self.nodes[ixa].position, &self.nodes[ixb].position);
        self.create_spring(node_distance, spring_constant, spring_dampening);
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

    pub fn add_grid(
        &mut self,
        GridConfig {
            ty,
            node_distance,
            node_radius,
            spring_constant,
            spring_dampening,
            width,
            height,
            upper_left,
        }: GridConfig,
    ) {
        let mut nodes = Vec::new();
        // create all the nodes
        for y in 0..height {
            for x in 0..width {
                let pos = upper_left + node_distance * V2::new(x as Float, y as Float);
                nodes.push(Node::new(self.next_node_id, pos, node_radius));
                self.next_node_id += 1;
            }
        }
        let base_ix = self.nodes.len();
        self.nodes.extend_from_slice(&nodes);
        // link the nodes
        for y in 0..height {
            for x in 0..width {
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
                match ty {
                    GridType::Quad => {}
                    GridType::Triangle => {
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
            let new_dist = distance_squared(pos, &node.position);
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
                    let dirab = other_node.position - node.position;
                    if dirab.norm() <= f64::EPSILON {
                        continue;
                    }
                    let dirba = -dirab;
                    // if to is closer to the spring than any of its nodes
                    let a = to - node.position;
                    let b = to - other_node.position;
                    if (dirab.dot(&a) > 0.0) && (dirba.dot(&b) > 0.0) {
                        // calculate the distance
                        let dirabn = dirab.normalize();
                        let dist = (a.norm_squared() - dirabn.dot(&a).powi(2)).sqrt();
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

    pub fn create_spring(
        &mut self,
        resting_length: Float,
        spring_constant: Float,
        damping_factor: Float,
    ) {
        let spring = Spring {
            spring_constant,
            resting_length,
            damping_factor,
            min_len: 0.25 * resting_length,
            max_len: 4.0 * resting_length,
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
        let resting_length = distance(&self.nodes[ixa].position, &self.nodes[ixb].position);
        let spring = Spring {
            spring_constant,
            resting_length,
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
                if attached.is_destruct(&self.nodes, node.position) {
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
        todo!()
    }

    pub fn clear(&mut self) {
        self.springs.clear();
        self.nodes.clear();
    }

    pub fn energy(&self) -> Float {
        let mut energy = 0.0;
        for node in &self.nodes {
            energy += 0.5 * node.weight * node.velocity.norm().powi(2);
            for attached in &node.attached {
                if let Some(spring) = attached.spring.upgrade() {
                    if let Some(other_node_ix) = self.get_node_ix(attached.other_node_id) {
                        let distance =
                            distance(&node.position, &self.nodes[other_node_ix].position);
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

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GridConfig {
    ty: GridType,
    node_distance: Float,
    node_radius: Float,
    spring_constant: Float,
    spring_dampening: Float,
    width: usize,
    height: usize,
    upper_left: P2,
}

#[derive(Debug, Clone, Copy)]
pub struct Spring {
    spring_constant: Float,
    resting_length: Float,
    damping_factor: Float,
    min_len: Float,
    max_len: Float,
}

impl Spring {
    fn force(&self, dir: &V2, vdiff: &V2) -> V2 {
        if *dir == V2::zeros() {
            V2::zeros()
        } else {
            let diff = (1.0 - (self.resting_length / dir.norm())).mul(5.0).powi(9) * dir;
            self.spring_constant * diff + self.damping_factor * vdiff.project_on(dir)
        }
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
        nodes[node_ix].position - main_node_position
    }

    pub fn length(&self, nodes: &[Node], main_node_position: P2) -> Float {
        self.relative_position(nodes, main_node_position).norm()
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
        print!("Node@({})", main_node.position.xy());
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
            attached: Vec::new(),
        }
    }

    pub fn bounds_collision(&mut self, width: Float, height: Float, elasticity: Float) {
        match (
            self.position[0] < self.radius,
            self.position[0] > (width - self.radius),
        ) {
            (true, _) => {
                if self.velocity[0] < 0. {
                    self.velocity[0] = -self.velocity[0];
                    self.velocity *= elasticity;
                }
                self.position[0] = self.radius;
            }
            (_, true) => {
                if self.velocity[0] > 0. {
                    self.velocity[0] = -self.velocity[0];
                    self.velocity *= elasticity;
                }
                self.position[0] = width - self.radius;
            }
            (false, false) => { // actually no collision occurred => do nothing
            }
        }
        match (
            self.position[1] < self.radius,
            self.position[1] > (height - self.radius),
        ) {
            (true, _) => {
                if self.velocity[1] < 0. {
                    self.velocity[1] = -self.velocity[1];
                    self.velocity *= elasticity;
                }
                self.position[1] = self.radius;
            }
            (_, true) => {
                if self.velocity[1] > 0. {
                    self.velocity[1] = -self.velocity[1];
                    self.velocity *= elasticity;
                }
                self.position[1] = height - self.radius;
            }
            (false, false) => { // actually no collision occurred => do nothing
            }
        }
    }
}

impl Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "velocity: {:.2}", self.velocity)?;
        writeln!(f, "acc: {:.2}", self.acceleration.norm())?;
        writeln!(f, "weight: {:.2}", self.weight)
    }
}

pub struct SpringyApp {
    app_mode: AppMode,
    contraption: Contraption,
    gravity_on: bool,
    gravity: V2,
    node_radius: Float,
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
            node_radius: 10.0,
            create_spring_constant: 10.0,
            create_spring_dampening: 10.0,
        }
    }

    /// select the nearest node or spring and conveniently return the coresponding AppMode
    pub fn select_nearest(&self) -> AppMode {
        let mouse_pos = mouse_position();
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
                        node.velocity = V2::zeros();
                        node.acceleration = V2::zeros();
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
                    ui.add(Slider::new(&mut self.node_radius, radius_range));
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

    fn grid_ui(&mut self, ui: &mut Ui) {
        ui.separator();
        ui.label("Grid:");
        let grid_types = [GridType::Quad, GridType::Triangle, GridType::Cross];
        ui.horizontal(|ui| {
            ComboBox::from_label("grid_type")
                .selected_text(format!("{:?}", self.grid_config.ty))
                .show_ui(ui, |ui| {
                    for gt in grid_types {
                        ui.selectable_value(&mut self.grid_config.ty, gt, format!("{:?}", gt));
                    }
                });
        });
        ui.horizontal(|ui| {
            ui.label("width");
            ui.add(Slider::new(&mut self.grid_config.width, 1..=100));
        });
        ui.horizontal(|ui| {
            ui.label("height");
            ui.add(Slider::new(&mut self.grid_config.height, 1..=100));
        });
        ui.horizontal(|ui| {
            ui.label("node distance");
            ui.add(Slider::new(
                &mut self.grid_config.node_distance,
                1.0..=100.0,
            ));
        });
        ui.horizontal(|ui| {
            ui.label("node radius");
            ui.add(Slider::new(&mut self.grid_config.node_radius, 1.0..=10.0));
        });
        ui.horizontal(|ui| {
            ui.label("spring constant");
            ui.add(Slider::new(
                &mut self.grid_config.spring_constant,
                10.0..=100000.0,
            ));
        });
        ui.horizontal(|ui| {
            ui.label("dampening");
            ui.add(
                Slider::new(&mut self.grid_config.spring_dampening, 0.0..=100.0).logarithmic(true),
            );
        });
        if ui.button("create grid").clicked() {
            self.contraption.add_grid(self.grid_config);
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
        mpos = mouse_position();
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
                if ui.button("Select").clicked() {
                    app.app_mode = AppMode::NodeSelected { node_id: None };
                }
                if ui.button("Draw Node").clicked() {
                    app.app_mode = AppMode::PlaceNode;
                }
                if ui.button("Draw Joints").clicked() {
                    app.app_mode = AppMode::SpringStart;
                }
                if ui.button("Clear").clicked() {
                    app.contraption.clear();
                }
                app.spring_ui(ui);
                app.node_ui(ui);
                app.grid_ui(ui);
                ui.label("spring velocities:");
            });
            mouse_hover = ctx.is_pointer_over_area();
        });

        if !mouse_hover && is_mouse_button_released(MouseButton::Left) {
            match &mut app.app_mode {
                AppMode::NodeSelected { node_id: _ } | AppMode::SpringSelected { spring: _ } => {
                    app.app_mode = app.select_nearest();
                }
                AppMode::PlaceNode => {
                    app.contraption
                        .insert_node(P2::new(mpos.0 as Float, mpos.1 as Float), app.node_radius);
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
            }
        }

        if is_key_pressed(KeyCode::P) {
            paused = !paused;
        }

        if paused && is_key_released(KeyCode::S) {
            app.contraption.rk4(delta_t as Float);
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
        }
        if is_mouse_button_down(MouseButton::Right) && !mouse_hover {
            pull_to_mouse(&mut app.contraption, 50.);
        }

        draw_contraption(&app.contraption);

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
            AppMode::SpringSelected { spring: _ } => ui.label("SpringSelected"),
            AppMode::NodeSelected { node_id: _ } => ui.label("NodeSelected"),
            AppMode::SpringStart => ui.label("JointStart"),
            AppMode::SpringEnd { node_a_id: _ } => ui.label("JointEnd"),
        }
    });
}

fn pull_to_mouse(contraption: &mut Contraption, force: Float) {
    let mpos = mouse_position();
    let mpos = P2::new(mpos.0 as Float, mpos.1 as Float);
    for node in contraption.nodes.iter_mut() {
        if node.stationary {
            continue;
        }
        let dir = mpos - node.position;
        let dist = dir.norm();
        let dir = dir.normalize();
        if dist > 0. {
            node.velocity += dir * (force * 10. / (dist + 10.));
        }
    }
}

fn draw_contraption(contraption: &Contraption) {
    for node in &contraption.nodes {
        let pos = node.position;
        if node.stationary {
            draw_circle(pos.x as f32, pos.y as f32, node.radius as f32, RED);
        } else {
            draw_circle(pos.x as f32, pos.y as f32, node.radius as f32, GREEN);
        }
        let pos_a = node.position;
        for attached in &node.attached {
            if let Some(other_ix) = contraption.get_node_ix(attached.other_node_id) {
                let pos_b = contraption.nodes[other_ix].position;
                draw_line(
                    pos_a[0] as f32,
                    pos_a[1] as f32,
                    pos_b[0] as f32,
                    pos_b[1] as f32,
                    2.,
                    WHITE,
                );
            }
        }
        draw_vector(&node.position, &node.velocity, YELLOW);
        draw_vector(&node.position, &(node.acceleration * 0.1), RED);
    }
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
    let angle = v.angle(&V2::new(s, 0.0));
    let mut angle = 180.0 * angle / std::f64::consts::PI;
    if s < 0.0 {
        angle += 60.0;
    }
    draw_poly(
        end_point[0] as f32,
        end_point[1] as f32,
        3,
        5.,
        angle as f32,
        color,
    );
}

#[derive(Clone)]
pub enum AppMode {
    NodeSelected { node_id: Option<NodeId> },
    SpringSelected { spring: Weak<RefCell<Spring>> },
    PlaceNode,
    SpringStart,
    SpringEnd { node_a_id: NodeId },
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
            _ => false,
        }
    }
}

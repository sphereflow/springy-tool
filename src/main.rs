use egui::*;
use macroquad::prelude::*;
use nalgebra::{distance, Point2, UnitVector2, Vector2};
use std::{
    cell::RefCell,
    rc::{Rc, Weak},
};

type V2 = Vector2<f32>;
type P2 = Point2<f32>;
type U2 = UnitVector2<f32>;

pub struct Contraption {
    pub nodes: Vec<Rc<RefCell<Node>>>,
    pub springs: Vec<Spring>,
}
impl Contraption {
    pub fn new() -> Contraption {
        Contraption {
            nodes: Vec::new(),
            springs: Vec::new(),
        }
    }

    pub fn apply_global_force(&mut self, force: V2) {
        for node in &self.nodes {
            node.borrow_mut().forces.push(force);
        }
    }

    pub fn apply_global_resistance(&mut self, resistance: f32) {
        for rcn in &self.nodes {
            let mut node = rcn.borrow_mut();
            let v_abs = node.velocity.norm();
            if v_abs > 0. {
                let force = v_abs * v_abs * resistance;
                let force_vec = -force * node.velocity.normalize();
                node.forces.push(force_vec);
            }
        }
    }

    pub fn bounds_collision(&mut self, width: f32, height: f32, elasticity: f32) {
        for rcn in &self.nodes {
            rcn.borrow_mut().bounds_collision(width, height, elasticity);
        }
    }

    pub fn nearest_to(&self, pos: &P2) -> Weak<RefCell<Node>> {
        let mut dist = f32::MAX;
        let mut res = Weak::new();
        for rcn in &self.nodes {
            let new_dist = distance(pos, &rcn.borrow().position);
            if new_dist < dist {
                dist = new_dist;
                res = Rc::downgrade(rcn);
            }
        }
        res
    }

    pub fn create_spring(
        &mut self,
        wa: Weak<RefCell<Node>>,
        wb: Weak<RefCell<Node>>,
        spring_constant: f32,
        damping_factor: f32,
    ) {
        let mut dist = 0.;
        if let Some((rca, rcb)) = wa.upgrade().zip(wb.upgrade()) {
            dist = distance(&rca.borrow().position, &rcb.borrow().position);
        }
        self.springs.push(Spring {
            spring_constant,
            resting_length: dist,
            damping_factor,
            min_len: 0.,
            max_len: 0.,
            node_a: wa,
            node_b: wb,
        });
    }

    pub fn update(&mut self, delta_t: f32) {
        // update forces
        for spring in &self.springs {
            if let (Some(node_a), Some(node_b)) = (spring.node_a.upgrade(), spring.node_b.upgrade())
            {
                let abn = (node_b.borrow().position - node_a.borrow().position).normalize();
                let force_a = spring.force();
                node_a.borrow_mut().forces.push(force_a * abn);
                node_b.borrow_mut().forces.push(-force_a * abn);
            }
        }
        // call update on nodes
        for node in &self.nodes {
            node.borrow_mut().update(delta_t);
        }
    }

    pub fn print_nodes(&self) {
        for rcn in &self.nodes {
            let node = rcn.borrow();
            print!("Node@({})", node.position.xy());
        }
        println!();
    }

    pub fn clear(&mut self) {
        self.springs.clear();
        self.nodes.clear();
    }
}
pub struct Spring {
    spring_constant: f32,
    resting_length: f32,
    damping_factor: f32,
    min_len: f32,
    max_len: f32,
    node_a: Weak<RefCell<Node>>,
    node_b: Weak<RefCell<Node>>,
}

impl Spring {
    pub fn length(&self) -> f32 {
        if let Some((rca, rcb)) = self.node_a.upgrade().zip(self.node_b.upgrade()) {
            distance(&rca.borrow().position, &rcb.borrow().position)
        } else {
            self.resting_length
        }
    }

    // which speed is this spring being streched/compressed with
    pub fn spring_velocity(&self) -> f32 {
        if let Some((rca, rcb)) = self.node_a.upgrade().zip(self.node_b.upgrade()) {
            let dir_ab: U2 = U2::new_normalize(rcb.borrow().position - rca.borrow().position);
            // velocity of node a in spring direction
            let va = rca.borrow().velocity_in_direction(&-dir_ab);
            let vb = rcb.borrow().velocity_in_direction(&dir_ab);
            (va + vb).norm()
        } else {
            0.
        }
    }

    pub fn force(&self) -> f32 {
        // damping is a force counteracting expansion / compression velocity
        (self.length() - self.resting_length) * self.spring_constant
            - self.damping_factor * self.spring_velocity()
    }
}

#[derive(PartialEq)]
pub struct Node {
    position: P2,
    radius: f32,
    velocity: V2,
    weight: f32,
    forces: Vec<V2>,
    stationary: bool,
}
impl Node {
    pub fn new(position: P2, radius: f32) -> Node {
        Node {
            position,
            radius,
            velocity: V2::new(0., 0.),
            weight: 1.0,
            forces: Vec::new(),
            stationary: false,
        }
    }

    pub fn bounds_collision(&mut self, width: f32, height: f32, elasticity: f32) {
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
    pub fn net_force(&self) -> V2 {
        self.forces.iter().sum()
    }

    pub fn velocity_in_direction(&self, dir: &U2) -> V2 {
        self.net_force().dot(dir) * dir.into_inner()
    }

    pub fn update(&mut self, delta_t: f32) {
        let net_force = self.net_force();
        // v = dt * a
        if !self.stationary {
            self.velocity += delta_t * net_force / self.weight;
            self.position += delta_t * self.velocity;
        }
        // remove all the forces because in the next update on Contraption they will be
        // recalculated
        self.forces.clear();
    }
}

#[macroquad::main("Springy")]
async fn main() {
    let mut app_mode = AppMode::NodeSelected { node: Weak::new() };
    let mut mpos;
    let mut mouse_hover = false;
    let mut contraption = Contraption::new();
    let mut global_resistance = 0.;
    let gravity = 1000. * V2::new(0., 1.);
    let mut gravity_on = true;
    let mut start = get_time();
    let mut stop;
    let sub_updates = 0.0001;
    let mut step_elapsed;
    let node_radius = 10.;
    loop {
        let width = screen_width();
        let height = screen_height();
        mpos = mouse_position();
        clear_background(BLACK);
        egui_macroquad::ui(|ctx| {
            Window::new("Springy").show(ctx, |ui| {
                ui.label(format!("Mouse: ({}, {})", mpos.0, mpos.1));
                ui.label(format!("num nodes: {}", contraption.nodes.len()));
                disp_app_mode(ui, &app_mode);
                match app_mode {
                    AppMode::SpringStart { .. } | AppMode::SpringEnd { .. } => {}
                    _ => {}
                }
                ui.checkbox(&mut gravity_on, "gravity");
                if ui.button("No Draw").clicked() {
                    app_mode = AppMode::NodeSelected { node: Weak::new() };
                }
                if ui.button("Draw Node").clicked() {
                    app_mode = AppMode::Node;
                }
                if ui.button("Draw Joints").clicked() {
                    app_mode = AppMode::SpringStart;
                }
                if ui.button("Clear").clicked() {
                    contraption.clear();
                }
                // we need to clone here because otherwise rustc decides to move app_mode into the
                // closure
                if let AppMode::NodeSelected { node } = &app_mode {
                    node_ui(ui, node.clone());
                }
                ui.add(Slider::new(&mut global_resistance, 0.0..=0.02));
                ui.label("spring velocities:");
                for spring in &contraption.springs {
                    ui.label(format!("{}", spring.spring_velocity()));
                }
            });
            mouse_hover = ctx.is_pointer_over_area();
        });

        if !mouse_hover && is_mouse_button_released(MouseButton::Left) {
            match app_mode {
                AppMode::NodeSelected { node: _ } => {
                    if let Some(rcell) = contraption.nearest_to(&P2::new(mpos.0, mpos.1)).upgrade()
                    {
                        app_mode = AppMode::NodeSelected {
                            node: Rc::downgrade(&rcell),
                        };
                    }
                }
                AppMode::Node => {
                    contraption.nodes.push(Rc::new(RefCell::new(Node::new(
                        P2::new(mpos.0, mpos.1),
                        node_radius,
                    ))));
                }
                AppMode::SpringStart => {
                    let mouse_position = P2::new(mpos.0, mpos.1);
                    app_mode = AppMode::SpringEnd {
                        node_a: contraption.nearest_to(&mouse_position),
                    };
                }
                AppMode::SpringEnd { node_a } => {
                    let mouse_position = P2::new(mpos.0, mpos.1);
                    contraption.create_spring(
                        node_a,
                        contraption.nearest_to(&mouse_position),
                        1000.0,
                        0.,
                    );
                    app_mode = AppMode::SpringStart;
                }
            }
        }

        // calculate delta_t
        stop = get_time();
        let elapsed = stop - start;
        start = stop;
        step_elapsed = sub_updates;

        loop {
            // update contraption
            if gravity_on {
                contraption.apply_global_force(gravity);
            }
            contraption.bounds_collision(width, height, 0.9);
            contraption.apply_global_resistance(global_resistance);
            if is_mouse_button_down(MouseButton::Right) && !mouse_hover {
                pull_to_mouse(&contraption, 5.);
            }
            if (step_elapsed + sub_updates) < elapsed {
                contraption.update(sub_updates as f32);
            } else {
                // final step
                contraption.update((elapsed - step_elapsed) as f32);
                break;
            }
            step_elapsed += sub_updates;
        }

        draw_contraption(&contraption);

        // draw the frame
        egui_macroquad::draw();
        next_frame().await
    }
}

fn node_ui(ui: &mut Ui, node: Weak<RefCell<Node>>) {
    if let Some(rcell) = node.upgrade() {
        let mut n = rcell.borrow_mut();
        ui.checkbox(&mut n.stationary, "stationary");
    }
}

fn disp_app_mode(ui: &mut Ui, mode: &AppMode) {
    ui.horizontal(|ui| {
        ui.label("App Mode: ");
        match mode {
            AppMode::Node => ui.label("Draw Nodes"),
            AppMode::NodeSelected { node: _ } => ui.label("NodeSelected"),
            AppMode::SpringStart => ui.label("JointStart"),
            AppMode::SpringEnd { node_a: _ } => ui.label("JointEnd"),
        }
    });
}

fn pull_to_mouse(contraption: &Contraption, force: f32) {
    let mpos = mouse_position();
    let mpos = P2::new(mpos.0, mpos.1);
    for rcn in &contraption.nodes {
        let mut node = rcn.borrow_mut();
        let dir = mpos - node.position;
        let dist = dir.norm();
        let dir = dir.normalize();
        if dist > 0. {
            node.velocity += dir * (force * 100. / (dist + 10.));
        }
    }
}

fn draw_contraption(contraption: &Contraption) {
    for rcn in &contraption.nodes {
        let node = rcn.borrow();
        let pos = node.position;
        if node.stationary {
            draw_circle(pos.x, pos.y, node.radius, RED);
        } else {
            draw_circle(pos.x, pos.y, node.radius, GREEN);
        }
    }
    for spring in &contraption.springs {
        if let Some((rca, rcb)) = spring.node_a.upgrade().zip(spring.node_b.upgrade()) {
            let pos_a = rca.borrow().position;
            let pos_b = rcb.borrow().position;
            draw_line(pos_a[0], pos_a[1], pos_b[0], pos_b[1], 2., WHITE);
        }
    }
}

enum AppMode {
    NodeSelected { node: Weak<RefCell<Node>> },
    Node,
    SpringStart,
    SpringEnd { node_a: Weak<RefCell<Node>> },
}

impl PartialEq for AppMode {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (AppMode::Node, AppMode::Node) => true,
            (AppMode::NodeSelected { node: _ }, AppMode::NodeSelected { node: _ }) => true,
            (AppMode::SpringEnd { node_a }, AppMode::SpringEnd { node_a: node_b }) => {
                if let Some((rca, rcb)) = node_a.upgrade().zip(node_b.upgrade()) {
                    *rca.borrow() == *rcb.borrow()
                } else {
                    false
                }
            }
            (AppMode::SpringStart, AppMode::SpringStart) => true,
            _ => false,
        }
    }
}

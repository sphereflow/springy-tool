use ::glam::{DVec2, Vec2};
use configs::{GridConfig, GridType, NodeConfig};
use core::f64;
pub(crate) use core::f64::consts::PI;
use egui::{Slider, Ui, Window};
use grid::Grid;
use macroquad::math::Rect;
use macroquad::prelude::*;
use spring::Spring;
pub(crate) use std::sync::{Arc, RwLock, Weak};
use std::{
    collections::VecDeque,
    time::{Duration, Instant},
};

use crate::{
    contraption::Contraption,
    egui_macroquad::Egui,
    node::{NEAREST_NODE_OVER_SPRING_OFFSET, NodeId},
};

pub mod configs;
pub mod contraption;
pub mod egui_macroquad;
pub mod grid;
pub mod motor;
pub mod node;
pub mod spring;

type Float = f64;
type V2 = DVec2;
type P2 = DVec2;
pub struct SpringyApp {
    app_mode: AppMode,
    contraption: Contraption,
    staging: Contraption,
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
    middle_mouse: MiddleMouse,
    min_fully_connect_distance: Float,
    max_fully_connect_distance: Float,
    pull_to_mouse_force: Float,
    brayon: bool,
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
            staging: Contraption::new(),
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
                node_distance: 15.0,
                node_radius: 4.0,
                spring_constant: 10000.0,
                spring_dampening: 4.0,
                spring_pretension_force: 0.0,
                width: 50,
                height: 50,
                upper_left: P2::new(200., 200.),
            },
            node_config: Default::default(),
            create_spring_constant: 10.0,
            create_spring_dampening: 10.0,
            middle_mouse: Default::default(),
            min_fully_connect_distance: 0.0,
            max_fully_connect_distance: 2000.0,
            pull_to_mouse_force: 200.0,
            brayon: false,
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
        if ui.button("Grab Node").clicked() {
            self.app_mode = AppMode::GrabNode { node_id: None };
            app_mode_changed = true;
        }
        if ui.button("Clear").clicked() {
            self.contraption.clear();
        }
        ui.horizontal(|ui| {
            ui.label("pull force");
            ui.add(Slider::new(&mut self.pull_to_mouse_force, 0.0..=10000.0));
        });
        if app_mode_changed {
            self.contraption.apply_resting_length_factors();
        }
    }

    pub fn rayon_ui(&mut self, ui: &mut Ui) {
        ui.checkbox(&mut self.brayon, "rayon");
    }

    /// select the nearest node or spring and conveniently return the coresponding AppMode
    pub fn select_nearest(&self, mouse_pos: (f32, f32)) -> AppMode {
        let mouse_pos = P2::new(mouse_pos.0 as Float, mouse_pos.1 as Float);
        let nearest_node = self.contraption.nearest_node_ix(&mouse_pos);
        let nearest_spring = self.contraption.nearest_spring(&mouse_pos);
        match (nearest_node, nearest_spring) {
            (None, None) => AppMode::NodeSelected { node_id: None },
            (None, Some((spring, _))) => AppMode::SpringSelected {
                spring: spring.clone(),
            },
            (Some((node_ix, _)), None) => AppMode::NodeSelected {
                node_id: Some(self.contraption.get_node(node_ix).id),
            },
            (Some((node_ix, node_dist)), Some((spring, spring_dist))) => {
                // prefer nodes since they are harder to select
                if node_dist <= (spring_dist + NEAREST_NODE_OVER_SPRING_OFFSET) {
                    AppMode::NodeSelected {
                        node_id: Some(self.contraption.get_node(node_ix).id),
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
                    let node = &mut self.contraption.get_node_mut(ix);
                    let mut stationary = node.is_stationary();
                    if ui.checkbox(&mut stationary, "stationary").changed() {
                        node.velocity = V2::ZERO;
                        node.acceleration = V2::ZERO;
                    }
                    node.set_stationary(stationary);
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
                if let Some(spring_cell) = spring.upgrade()
                    && let Ok(mut spring_ref) = spring_cell.write()
                {
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

    pub fn circle_ui(&mut self, ui: &mut Ui) {
        if ui.button("add circle").clicked() {
            self.app_mode = AppMode::PlaceCircle {
                radius: 50.0,
                num_nodes: 4,
                node_config: NodeConfig::default(),
            };
        }
        if let AppMode::PlaceCircle { node_config, .. } = &mut self.app_mode {
            node_config.ui(ui);
        }
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
                for node in self.contraption.nodes_iter_mut() {
                    if Some(&node.id) == current_id {
                        current_id = node_ids_iter.next();
                        node.configure(self.node_config);
                    }
                }
            }
            let master_spring = springs.first();
            if let Some(master_spring_ref) = master_spring
                && let Some(master_spring_rc) = master_spring_ref.upgrade()
                && let Ok(option_change) = master_spring_rc.write().map(|mut ms| ms.ui(ui))
                && let Some(change) = option_change
                && let Ok(copy_spring) = master_spring_rc.read()
            {
                for weak_spring in springs.iter() {
                    if let Some(copy_to_spring_rc) = weak_spring.upgrade() {
                        // resting length should not be modified by the UI
                        // resting_length_factor should be instead
                        if Arc::ptr_eq(&master_spring_rc, &copy_to_spring_rc) {
                            continue;
                        }
                        if let Ok(mut spring) = copy_to_spring_rc.write() {
                            match change {
                                spring::SpringChange::Constant => {
                                    spring.spring_constant = copy_spring.spring_constant
                                }
                                spring::SpringChange::RestingLengthFactor => {
                                    spring.resting_length_factor = copy_spring.resting_length_factor
                                }
                                spring::SpringChange::PretensionForce => {
                                    spring.pretension_force = copy_spring.pretension_force
                                }
                                spring::SpringChange::DampingFactor => {
                                    spring.damping_factor = copy_spring.damping_factor
                                }
                                spring::SpringChange::MinLength => {
                                    spring.min_len = copy_spring.min_len
                                }
                                spring::SpringChange::MaxLength => {
                                    spring.max_len = copy_spring.max_len
                                }
                            }
                        }
                    }
                }
            }

            if ui.button("fully connect").clicked() {
                self.contraption.fully_connect(
                    node_ids,
                    1000.0,
                    4.0,
                    self.min_fully_connect_distance,
                    self.max_fully_connect_distance,
                );
            }
            if let AppMode::MultiSelection { node_ids, springs } = &self.app_mode
                && ui.button("grab selection").clicked()
            {
                let mut nodes_offsets = Vec::new();
                let mut middle = P2::ZERO;
                let mut middle_count = 0.0;
                for id in node_ids {
                    if let Some(ix) = self.contraption.get_node_ix(*id) {
                        middle += self.contraption.get_node(ix).position();
                        middle_count += 1.0;
                    }
                }
                if middle_count != 0.0 {
                    middle /= middle_count;
                }
                for id in node_ids {
                    if let Some(ix) = self.contraption.get_node_ix(*id) {
                        nodes_offsets
                            .push((*id, self.contraption.get_node(ix).position() - middle));
                    }
                }
                self.app_mode = AppMode::GrabSelection {
                    node_ids_offsets: nodes_offsets,
                    springs: springs.clone(),
                }
            }
            ui.horizontal(|ui| {
                ui.label("min distance");
                ui.add(Slider::new(
                    &mut self.min_fully_connect_distance,
                    0.0..=self.max_fully_connect_distance,
                ));
            });
            ui.horizontal(|ui| {
                ui.label("max distance");
                ui.add(Slider::new(
                    &mut self.max_fully_connect_distance,
                    self.min_fully_connect_distance..=2000.0,
                ));
            });
        }
    }

    pub fn draw_contraption(&self, mpos: (f32, f32)) {
        let mpos = DVec2::new(mpos.0 as Float, mpos.1 as Float);
        let mut nearest_spring = None;
        let mut nearest_node = None;
        if let AppMode::NodeSelected { node_id: None } = &self.app_mode {
            match (
                self.contraption.nearest_spring(&mpos),
                self.contraption.nearest_node_ix(&mpos),
            ) {
                (None, None) => {}
                (None, Some((node_ix, _))) => nearest_node = Some(node_ix),
                (Some((spring, _)), None) => nearest_spring = Some(spring),
                (Some((spring, spring_dist)), Some((node_ix, node_dist))) => {
                    if (spring_dist + NEAREST_NODE_OVER_SPRING_OFFSET) < node_dist {
                        nearest_spring = Some(spring);
                    } else {
                        nearest_node = Some(node_ix);
                    }
                }
            }
        }
        for (node_ix, node) in self.contraption.nodes_iter().enumerate() {
            let pos = node.position();
            let mut node_color = GREEN;
            if (AppMode::NodeSelected {
                node_id: Some(node.id),
            }) == self.app_mode
            {
                node_color = ORANGE;
            } else if Some(node_ix) == nearest_node {
                node_color = PURPLE;
            } else if node.is_stationary() {
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
                    color.a = 0.75;
                    draw_circle(pos.x as f32, pos.y as f32, 2.0 * node.radius as f32, color);
                }
            }
            if let AppMode::GrabSelection {
                node_ids_offsets,
                springs: _,
            } = &self.app_mode
            {
                if node_ids_offsets
                    .binary_search_by_key(&node.id, |id_offset| id_offset.0)
                    .is_ok()
                {
                    let mut color = SKYBLUE;
                    color.a = 0.75;
                    draw_circle(pos.x as f32, pos.y as f32, 2.0 * node.radius as f32, color);
                }
            }
            let pos_a = node.position();
            for attached in &node.attached {
                if let Some(other_ix) = self.contraption.get_node_ix(attached.other_node_id) {
                    let pos_b = self.contraption.get_node(other_ix).position();
                    if let Some(spring) = attached.spring.upgrade() {
                        let dist = pos_a.distance(pos_b);
                        if let Ok(spring) = spring.read() {
                            let max_len = spring.max_len;
                            let stretch_factor = if dist < spring.resting_length {
                                (spring.resting_length - dist) / spring.resting_length
                            } else {
                                (dist - spring.resting_length) / (max_len - spring.resting_length)
                            };
                            let color = colorous::RED_YELLOW_GREEN
                                .eval_continuous((1.0 - stretch_factor).powi(4));
                            let mut spring_color =
                                macroquad::color::Color::from_rgba(color.r, color.g, color.b, 255);
                            if (AppMode::SpringSelected {
                                spring: attached.spring.clone(),
                            }) == self.app_mode
                            {
                                spring_color = ORANGE;
                            } else if let Some(nearest) = nearest_spring.as_ref()
                                && nearest.as_ptr() == attached.spring.as_ptr()
                            {
                                spring_color = PURPLE;
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
            }

            // draw velocity and acceleration arrows
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

        // draw staging area
        for node in self.staging.nodes_iter() {
            draw_circle(
                node.position.x as f32,
                node.position.y as f32,
                node.radius as f32,
                SKYBLUE,
            );
            for attached in node.attached.iter() {
                if let Some(other_ix) = self.staging.get_node_ix(attached.other_node_id) {
                    let other_pos = self.staging.get_node(other_ix).position;
                    draw_line(
                        node.position.x as f32,
                        node.position.y as f32,
                        other_pos.x as f32,
                        other_pos.y as f32,
                        3.0,
                        SKYBLUE,
                    );
                }
            }
        }
    }

    fn pull_to_mouse(&mut self, mpos: (f32, f32)) {
        let mpos = P2::new(mpos.0 as Float, mpos.1 as Float);
        for node in self.contraption.nodes_iter_mut() {
            if node.is_stationary() {
                continue;
            }
            let dir = mpos - node.position;
            let dist = dir.length();
            let dir = dir.normalize();
            if dist > 0. {
                node.velocity += dir * (self.pull_to_mouse_force * 10. / (dist + 10.));
            }
        }
    }

    pub fn update_middle_mouse(&mut self) {
        self.middle_mouse.update();
        if let AppMode::PlaceCircle {
            radius,
            num_nodes,
            node_config: _,
        } = &mut self.app_mode
        {
            match self.middle_mouse.change {
                MiddleMouseChange::Radius => self.middle_mouse.update_float(radius),
                MiddleMouseChange::NumNodes => self.middle_mouse.update_usize(num_nodes),
            }
        }
    }
}

#[macroquad::main("Springy")]
async fn main() {
    let mut app = SpringyApp::new();
    let mut egui = Egui::new();
    let mut mpos_prev;
    let mut mpos = (0.0, 0.0);
    let mut mouse_hover = false;
    let mut start = get_time();
    let mut stop;
    let mut delta_t = 0.01;
    let mut compute_time = 0.;
    let mut step_elapsed;
    let mut paused: bool = false;
    loop {
        let width = screen_width() as Float;
        let height = screen_height() as Float;
        mpos_prev = mpos;
        mpos = app.grid.mouse_position();
        clear_background(BLACK);
        app.update_middle_mouse();
        {
            egui.calc_ui(|ctx| {
                Window::new("Springy").show(ctx, |ui| {
                    ui.label(format!("Mouse: ({}, {})", mpos.0, mpos.1));
                    ui.label(format!("compute time: {:.2} ms", compute_time * 1000.));
                    ui.label(format!("num nodes: {}", app.contraption.num_nodes()));
                    ui.label(format!("num springs: {}", app.contraption.springs.len()));
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
                        AppMode::SpringStart | AppMode::SpringEnd { .. } => {}
                        _ => {}
                    }
                    ui.horizontal(|ui| {
                        ui.label("delta_t");
                        ui.add(Slider::new(&mut delta_t, 0.0001..=0.1).logarithmic(true));
                    });
                    ui.checkbox(&mut app.gravity_on, "gravity");
                    app.app_mode_ui(ui);
                    app.rayon_ui(ui);
                    app.arrows_ui(ui);
                    app.grid.ui(ui);
                    app.spring_ui(ui);
                    app.node_ui(ui);
                    app.circle_ui(ui);
                    app.multi_selection_ui(ui);
                    app.grid_config.ui(ui, &mut app.contraption);
                    ui.label("spring velocities:");
                });
                mouse_hover = ctx.is_pointer_over_area();
            });
        }

        // MultiSelection: store the mouse position where the selection starts
        if !mouse_hover && is_mouse_button_pressed(MouseButton::Left) {
            if let AppMode::MultiSelection { .. } = app.app_mode {
                app.multi_selection_rect = Some((
                    macroquad::math::Rect::new(mpos.0, mpos.1, 0.0, 0.0),
                    Vec2::new(mpos.0, mpos.1),
                ));
            }
        }

        // draw multiselection rect
        if !mouse_hover && is_mouse_button_down(MouseButton::Left) {
            if let AppMode::MultiSelection { .. } = app.app_mode {
                if let Some((rect, start_pos)) = app.multi_selection_rect.as_mut() {
                    let mut selection_rect_color = macroquad::color::SKYBLUE;
                    selection_rect_color.a = 0.5;
                    rect.w = (start_pos.x - mpos.0).abs();
                    rect.h = (start_pos.y - mpos.1).abs();
                    rect.x = start_pos.x.min(mpos.0);
                    rect.y = start_pos.y.min(mpos.1);
                    draw_rectangle(rect.x, rect.y, rect.w, rect.h, selection_rect_color);
                }
            }
        }

        if !mouse_hover && is_mouse_button_released(MouseButton::Left) {
            app.contraption.apply_resting_length_factors();
            match &mut app.app_mode {
                AppMode::NodeSelected { node_id: _ } | AppMode::SpringSelected { spring: _ } => {
                    app.app_mode = app.select_nearest(mpos);
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
                                app.node_config =
                                    NodeConfig::from_node(&app.contraption.get_node(ix));
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
                            node_a_id: app.contraption.get_node(ixa).id,
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
                AppMode::PlaceCircle { .. } => {
                    app.contraption.merge(&mut app.staging);
                }
                AppMode::CutSpring => {}
                AppMode::GrabNode { .. } => {}
                AppMode::GrabSelection { .. } => {}
            }
        }

        if !is_mouse_button_down(MouseButton::Left) {
            app.staging.clear();
        }
        if !mouse_hover && !is_mouse_button_down(MouseButton::Left) {
            if let AppMode::PlaceCircle {
                radius,
                num_nodes,
                node_config,
            } = app.app_mode
            {
                app.staging.add_circle(
                    P2::new(mpos.0 as Float, mpos.1 as Float),
                    radius,
                    num_nodes,
                    node_config,
                );
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
                app.contraption.bounds_collision(width, height);
                if (step_elapsed + delta_t) < elapsed {
                    if app.brayon {
                        app.contraption.rk4_rayon(delta_t as Float);
                    } else {
                        app.contraption.rk4(delta_t as Float);
                    }
                    step_elapsed += delta_t;
                } else {
                    // final step
                    if app.brayon {
                        app.contraption.rk4_rayon((elapsed - step_elapsed) as Float);
                    } else {
                        app.contraption.rk4((elapsed - step_elapsed) as Float);
                    }
                    break;
                }
                app.contraption.spring_break();
            }
            app.contraption.step_motors(elapsed);
        }
        compute_time = get_time() - start;
        if is_mouse_button_down(MouseButton::Right) && !mouse_hover {
            app.pull_to_mouse(mpos);
        }
        if is_mouse_button_down(MouseButton::Left) && (app.app_mode == AppMode::CutSpring) {
            if let Some((spring_ref, dist)) = app
                .contraption
                .nearest_spring(&P2::new(mpos.0 as Float, mpos.1 as Float))
            {
                if dist < 2.0 {
                    app.contraption.delete_spring(spring_ref);
                }
            }

            let ma = P2::new(mpos_prev.0 as Float, mpos_prev.1 as Float);
            let mb = P2::new(mpos.0 as Float, mpos.1 as Float);
            app.contraption.cut_springs(ma, mb);
        }
        if is_mouse_button_down(MouseButton::Left) && !mouse_hover {
            let pos = P2::new(mpos.0 as Float, mpos.1 as Float);
            if let AppMode::GrabNode {
                node_id: Some(node_id),
            } = app.app_mode
            {
                if let Some(ix) = app.contraption.get_node_ix(node_id) {
                    let node = app.contraption.get_node_mut(ix);
                    node.position = pos;
                    node.velocity = V2::ZERO;
                }
            } else if let AppMode::GrabNode { node_id: None } = app.app_mode {
                if let Some((ix, _)) = app.contraption.nearest_node_ix(&pos) {
                    app.app_mode = AppMode::GrabNode {
                        node_id: Some(app.contraption.get_node(ix).id),
                    }
                }
            }
            if let AppMode::GrabSelection {
                node_ids_offsets,
                springs: _,
            } = &app.app_mode
            {
                for (id, offset) in node_ids_offsets {
                    if let Some(ix) = app.contraption.get_node_ix(*id) {
                        let node = app.contraption.get_node_mut(ix);
                        node.position = pos + offset;
                        node.velocity = V2::ZERO;
                    }
                }
            }
        }
        app.draw_contraption(mpos);
        // draw the grid if shift is pressed
        if is_key_down(KeyCode::LeftShift) || is_key_down(KeyCode::RightShift) {
            app.grid.draw();
        }
        // draw the frame
        egui.draw();
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
            AppMode::PlaceCircle { .. } => ui.label("PlaceCircle"),
            AppMode::CutSpring => ui.label("CutSpring"),
            AppMode::GrabNode { .. } => ui.label("GrabNode"),
            AppMode::GrabSelection { .. } => ui.label("GrabSelection"),
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

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum MiddleMouseChange {
    Radius,
    NumNodes,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MiddleMouse {
    last_scroll: Instant,
    scroll_increase_threshold: Duration,
    scroll_step: Float,
    wheel: (f32, f32),
    change: MiddleMouseChange,
    scroll_timings: VecDeque<(Instant, f32)>,
}

impl Default for MiddleMouse {
    fn default() -> Self {
        Self {
            last_scroll: Instant::now(),
            scroll_increase_threshold: Duration::from_millis(100),
            scroll_step: 1.0,
            wheel: mouse_wheel(),
            change: MiddleMouseChange::Radius,
            scroll_timings: VecDeque::new(),
        }
    }
}

impl MiddleMouse {
    pub fn update(&mut self) {
        self.wheel = mouse_wheel();
        if is_mouse_button_released(MouseButton::Middle) {
            match self.change {
                MiddleMouseChange::Radius => self.change = MiddleMouseChange::NumNodes,
                MiddleMouseChange::NumNodes => self.change = MiddleMouseChange::Radius,
            }
        }
        let now = Instant::now();
        let passed = now - self.last_scroll;
        self.scroll_timings.push_back((now, self.wheel.1));
        if self.scroll_timings.len() > 10 {
            self.scroll_timings.pop_front();
        }
        if self.wheel.1 != 0.0 {
            if passed <= self.scroll_increase_threshold {
                self.scroll_step *= 1.3;
            }
            self.last_scroll = now;
        }
        if (self.scroll_step > 1.5) && (passed >= self.scroll_increase_threshold) {
            self.scroll_step /= 1.3;
            self.last_scroll = now;
        }
    }

    fn scroll_density(&self) -> f32 {
        let mut acc = 0.0;
        if let Some((start, _)) = self.scroll_timings.front() {
            for (current_time, val) in self.scroll_timings.iter() {
                let time = (*current_time - *start).as_secs_f32();
                let numerator = *val * 1.5_f32.powf(time);
                acc += numerator;
            }
        }
        acc
    }

    pub fn update_float(&self, val: &mut Float) {
        //*val += self.wheel.1 as Float * self.scroll_step;
        let acc = self.scroll_density();
        let density = (acc.powi(3) - 0.19 * acc.signum()).abs().floor() * acc.signum();
        *val += density as Float / 5.0;
    }

    pub fn update_usize(&self, val: &mut usize) {
        let mut float_val = *val as f64;
        let acc = self.scroll_density();
        let density = (acc - 0.06 * acc.signum()).abs().floor() * acc.signum();
        float_val += density as Float;
        *val = float_val as usize;
    }
}

#[derive(Clone, Debug)]
pub enum AppMode {
    MultiSelection {
        node_ids: Vec<NodeId>,
        springs: Vec<Weak<RwLock<Spring>>>,
    },
    NodeSelected {
        node_id: Option<NodeId>,
    },
    SpringSelected {
        spring: Weak<RwLock<Spring>>,
    },
    PlaceNode,
    SpringStart,
    SpringEnd {
        node_a_id: NodeId,
    },
    PlaceCircle {
        radius: Float,
        num_nodes: usize,
        node_config: NodeConfig,
    },
    CutSpring,
    GrabNode {
        node_id: Option<NodeId>,
    },
    GrabSelection {
        node_ids_offsets: Vec<(NodeId, P2)>,
        springs: Vec<Weak<RwLock<Spring>>>,
    },
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
            (AppMode::GrabNode { .. }, AppMode::GrabNode { .. }) => true,
            (AppMode::GrabSelection { .. }, AppMode::GrabSelection { .. }) => true,
            _ => false,
        }
    }
}

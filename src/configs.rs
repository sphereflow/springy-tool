use crate::{motor::Motor, Contraption, Float, Node, P2};
use egui::{ComboBox, Slider, Ui};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum GridType {
    Quad,
    QuadTriangle,
    Cross,
    Triangle,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GridConfig {
    pub ty: GridType,
    pub node_distance: Float,
    pub node_radius: Float,
    pub spring_constant: Float,
    pub spring_dampening: Float,
    pub width: usize,
    pub height: usize,
    pub upper_left: P2,
}

impl GridConfig {
    pub fn ui(&mut self, ui: &mut Ui, contraption: &mut Contraption) {
        ui.separator();
        ui.label("Grid:");
        let grid_types = [
            GridType::Quad,
            GridType::QuadTriangle,
            GridType::Cross,
            GridType::Triangle,
        ];
        ui.horizontal(|ui| {
            ComboBox::from_label("grid_type")
                .selected_text(format!("{:?}", self.ty))
                .show_ui(ui, |ui| {
                    for gt in grid_types {
                        ui.selectable_value(&mut self.ty, gt, format!("{:?}", gt));
                    }
                });
        });
        ui.horizontal(|ui| {
            ui.label("width");
            ui.add(Slider::new(&mut self.width, 1..=100));
        });
        ui.horizontal(|ui| {
            ui.label("height");
            ui.add(Slider::new(&mut self.height, 1..=100));
        });
        ui.horizontal(|ui| {
            ui.label("node distance");
            ui.add(Slider::new(&mut self.node_distance, 1.0..=100.0));
        });
        ui.horizontal(|ui| {
            ui.label("node radius");
            ui.add(Slider::new(&mut self.node_radius, 1.0..=10.0));
        });
        ui.horizontal(|ui| {
            ui.label("spring constant");
            ui.add(Slider::new(&mut self.spring_constant, 10.0..=100000.0));
        });
        ui.horizontal(|ui| {
            ui.label("dampening");
            ui.add(Slider::new(&mut self.spring_dampening, 0.0..=100.0).logarithmic(true));
        });
        if ui.button("create grid").clicked() {
            contraption.add_grid(*self);
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NodeConfig {
    pub stationary: bool,
    pub node_radius: Float,
    pub motor: Option<Motor>,
}

impl Default for NodeConfig {
    fn default() -> Self {
        Self {
            stationary: false,
            node_radius: 1.0,
            motor: None,
        }
    }
}

impl NodeConfig {
    pub fn from_node(node: &Node) -> Self {
        NodeConfig {
            stationary: node.stationary,
            node_radius: node.radius,
            motor: node.motor,
        }
    }
    pub fn ui(&mut self, ui: &mut Ui) -> bool {
        let radius_range = 0.1..=10.0;
        let mut changed = false;
        changed |= ui.checkbox(&mut self.stationary, "stationary").changed();
        changed |= ui
            .add(Slider::new(&mut self.node_radius, radius_range))
            .changed();
        if self.stationary {
            let mut has_motor = self.motor.is_some();
            changed |= ui.checkbox(&mut has_motor, "Motor?").changed();
            if has_motor {
                ui.separator();
                changed |= self.motor.get_or_insert(Motor::default()).ui(ui);
            } else {
                self.motor = None;
            }
        }
        changed
    }
}

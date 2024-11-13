use std::f64::consts::PI;

use crate::{Float, P2};
use egui::{Slider, Ui};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum Axis {
    X,
    Y,
    #[default]
    XY,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Motor {
    axis: Axis,
    phase: Float,
    frequency: Float,
    amplitude: Float,
}

impl Motor {
    pub fn ui(&mut self, ui: &mut Ui) -> bool {
        let mut changed = false;
        changed |= ui.radio_value(&mut self.axis, Axis::X, "X").changed();
        changed |= ui.radio_value(&mut self.axis, Axis::Y, "Y").changed();
        changed |= ui.radio_value(&mut self.axis, Axis::XY, "XY").changed();
        ui.horizontal(|ui| {
            ui.label("amplitude");
            changed |= ui
                .add(Slider::new(&mut self.amplitude, 0.0..=10.0))
                .changed();
        });
        ui.horizontal(|ui| {
            ui.label("frequency");
            changed |= ui
                .add(Slider::new(&mut self.frequency, 0.0..=10.0))
                .changed();
        });
        ui.horizontal(|ui| {
            ui.label("phase");
            changed |= ui
                .add(Slider::new(&mut self.phase, 0.0..=(2.0 * PI)))
                .changed();
        });
        changed
    }

    pub fn step(&mut self, delta_t: Float) {
        self.phase += 2.0 * PI * self.frequency * delta_t;
        if self.phase > 2.0 * PI {
            self.phase -= 2.0 * PI;
        }
    }

    pub fn position(&self, base_point: P2) -> P2 {
        match self.axis {
            Axis::X => base_point + (P2::new(self.phase.cos(), 0.0) * self.amplitude),
            Axis::Y => base_point + (P2::new(0.0, self.phase.sin()) * self.amplitude),
            Axis::XY => base_point + (P2::new(self.phase.cos(), self.phase.sin()) * self.amplitude),
        }
    }
}

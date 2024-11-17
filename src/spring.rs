use crate::{Float, V2};
use egui::{Slider, Ui};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Spring {
    pub spring_constant: Float,
    pub resting_length: Float,
    pub resting_length_factor: Float,
    pub pretension_force: Float,
    pub damping_factor: Float,
    pub min_len: Float,
    pub max_len: Float,
}

impl Spring {
    pub fn force(&self, dir: &V2, vdiff: &V2) -> V2 {
        if *dir == V2::ZERO {
            V2::ZERO
        } else {
            // magnitude of the diff vector
            let diffmag = 1.0 - ((self.resting_length * self.resting_length_factor) / dir.length());
            let diff = diffmag * dir;
            self.pretension_force * dir.normalize() * diffmag.signum()
                + self.spring_constant * diff
                + self.damping_factor * vdiff.project_onto(*dir)
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
            ui.label("pretension force");
            changed |= ui
                .add(Slider::new(&mut self.pretension_force, 0.0..=100000.0).logarithmic(true))
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

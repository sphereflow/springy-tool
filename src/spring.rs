use std::collections::HashMap;

use crate::{
    Float, P2, V2,
    node::{Node, NodeId},
};
use egui::{Slider, Ui};

pub type SpringId = usize;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpringChange {
    Constant,
    RestingLengthFactor,
    PretensionForce,
    DampingFactor,
    MinLength,
    MaxLength,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Spring {
    pub id: SpringId,
    pub spring_constant: Float,
    pub resting_length: Float,
    pub resting_length_factor: Float,
    pub pretension_force: Float,
    pub damping_factor: Float,
    pub min_len: Float,
    pub max_len: Float,
    pub node_id_a: NodeId,
    pub node_id_b: NodeId,
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
        self.min_len *= self.resting_length_factor;
        self.max_len *= self.resting_length_factor;
        self.resting_length_factor = 1.0;
    }

    pub fn is_destruct(&self, nodes: &[Node], node_map: &HashMap<NodeId, usize>) -> bool {
        // search for the 2 nodes connecting this spring with the ids
        let ixa = node_map[&self.node_id_a];
        let ixb = node_map[&self.node_id_b];
        let length = (nodes[ixa].position - nodes[ixb].position).length();
        !((self.min_len * self.resting_length_factor)..=(self.max_len * self.resting_length_factor))
            .contains(&length)
    }

    pub fn crosses(
        &self,
        nodes: &[Node],
        node_map: &HashMap<NodeId, usize>,
        ma: P2,
        mb: P2,
    ) -> bool {
        let ixa = node_map[&self.node_id_a];
        let ixb = node_map[&self.node_id_b];

        let nap = nodes[ixa].position;
        let nbp = nodes[ixb].position;
        let va = ma - nap;
        let vb = mb - nap;
        let vc = nbp - ma;
        let vs = nbp - nap;
        let vm = mb - ma;
        (vs.perp_dot(va) * vs.perp_dot(vb)) <= 0.0 && (vm.perp_dot(va) * vm.perp_dot(vc)) >= 0.0
    }

    pub fn ui(&mut self, ui: &mut Ui) -> Option<SpringChange> {
        let mut changed = None;
        ui.horizontal(|ui| {
            ui.label("spring constant");
            if ui
                .add(Slider::new(&mut self.spring_constant, 0.0..=100000.0).logarithmic(true))
                .changed()
            {
                changed = Some(SpringChange::Constant)
            }
        });
        ui.horizontal(|ui| {
            ui.label("resting length factor");
            if ui
                .add(Slider::new(&mut self.resting_length_factor, 0.0..=2.0))
                .changed()
            {
                changed = Some(SpringChange::RestingLengthFactor)
            }
        });
        ui.horizontal(|ui| {
            ui.label("pretension force");
            if ui
                .add(Slider::new(&mut self.pretension_force, 0.0..=100000.0).logarithmic(true))
                .changed()
            {
                changed = Some(SpringChange::PretensionForce);
            }
        });
        ui.horizontal(|ui| {
            ui.label("min len");
            if ui
                .add(Slider::new(
                    &mut self.min_len,
                    0.0..=(self.resting_length * self.resting_length_factor),
                ))
                .changed()
            {
                changed = Some(SpringChange::MinLength);
            }
        });
        ui.horizontal(|ui| {
            ui.label("max len");
            if ui
                .add(Slider::new(
                    &mut self.max_len,
                    (self.resting_length * self.resting_length_factor)
                        ..=(10.0 * self.resting_length * self.resting_length_factor),
                ))
                .changed()
            {
                changed = Some(SpringChange::MaxLength);
            }
        });
        ui.horizontal(|ui| {
            ui.label("damping factor");
            if ui
                .add(Slider::new(&mut self.damping_factor, 0.0..=100.0).logarithmic(true))
                .changed()
            {
                changed = Some(SpringChange::DampingFactor);
            }
        });
        changed
    }
}

use macroquad::math::Rect;
use rayon::prelude::*;
use std::{
    cmp::Ordering,
    collections::HashMap,
    sync::{Arc, RwLock, Weak},
};

use crate::{
    Float, P2, V2,
    configs::{GridConfig, GridType, NodeConfig},
    node::{Attached, Node, NodeId},
    spring::Spring,
};

#[derive(Debug, Clone)]
pub struct Contraption {
    nodes: Vec<Node>,
    node_map: HashMap<NodeId, usize>,
    next_node_id: NodeId,
    pub springs: Vec<Arc<RwLock<Spring>>>,
    // intermediate values for Runge Kutta
    // stored in the struct to avoid reallocations
    kps: [Vec<P2>; 4],
    kvs: [Vec<V2>; 4],
    kaccs: [Vec<V2>; 4],
    pub g: V2,
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
            node_map: HashMap::new(),
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
                if self.nodes[i].is_stationary() {
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
            if self.nodes[i].is_stationary() {
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

    pub fn rk4_rayon(&mut self, delta_t: Float) {
        let factors = [0.5, 0.5, 1.0];
        self.rk4_init();

        // Runge Kutta steps

        for (n, factor) in factors.iter().enumerate() {
            let [_, kpsnp1, ..] = &mut self.kps[n..] else {
                panic!("Contraption::kps does not have enough elements")
            };
            let [kvsn, kvsnp1, ..] = &mut self.kvs[n..] else {
                panic!("Contraption::kvs does not have enough elements")
            };
            let [kaccsn, kaccsnp1, ..] = &mut self.kaccs[n..] else {
                panic!("Contraption::kaccs does not have enough elements")
            };

            // calculate positions and velocities
            // this does not involve springs
            (&self.nodes, kpsnp1, kvsn, kvsnp1, kaccsn)
                .into_par_iter()
                .for_each(|(node, kpnp1, kvn, kvnp1, kaccn)| {
                    if node.is_stationary() {
                        *kpnp1 = node.position();
                        return;
                    }
                    *kpnp1 = node.position + factor * delta_t * *kvn;
                    *kvnp1 = node.velocity + factor * delta_t * *kaccn;
                });

            (&self.nodes, kaccsnp1)
                .into_par_iter()
                .enumerate()
                .for_each(|(i, (node, kaccsnp1))| {
                    *kaccsnp1 = self.g
                        + Self::acceleration_par(
                            &self.node_map,
                            i,
                            node,
                            &self.kps[n + 1],
                            &self.kvs[n + 1],
                        );
                });
        }

        // updating nodes
        let [ref mut kp0s, ..] = self.kps;
        let [ref mut kv0s, ref kv1s, ref kv2s, ref kv3s] = self.kvs;
        let [ref kacc0s, ref kacc1s, ref kacc2s, ref kacc3s] = self.kaccs;

        (
            &mut self.nodes,
            kp0s,
            kv0s,
            kv1s,
            kv2s,
            kv3s,
            kacc0s,
            kacc1s,
            kacc2s,
            kacc3s,
        )
            .into_par_iter()
            .for_each(
                |(node, kp0, kv0, kv1, kv2, kv3, kacc0, kacc1, kacc2, kacc3)| {
                    if node.is_stationary() {
                        return;
                    }
                    // calculate position
                    let weighted_v = (*kv0 + kv1 * 2.0 + kv2 * 2.0 + kv3) * (1.0 / 6.0);
                    node.position += delta_t * weighted_v;

                    // calculate velocity
                    let weighted_acc = (kacc0 + kacc1 * 2.0 + kacc2 * 2.0 + kacc3) * (1.0 / 6.0);
                    node.velocity += delta_t * weighted_acc;

                    // put all the nodes positions into kps[0]
                    *kp0 = node.position();
                    // put all the nodes velocities into kvs[0]
                    *kv0 = node.velocity;
                },
            );

        // calculate acceleration
        for i in 0..self.nodes.len() {
            // use the positions and velocities of self.nodes which where collected in the previous loop
            self.nodes[i].acceleration = self.g + self.acceleration(i, &self.kps[0], &self.kvs[0]);
        }
    }

    pub fn acceleration(&self, node_ix: usize, kps: &[P2], kvs: &[V2]) -> V2 {
        let mut forces = V2::ZERO;
        let node = &self.nodes[node_ix];
        for attached in node.attached.iter() {
            if let Some(spring) = attached.spring.upgrade()
                && let Some(other_ix) = self.get_node_ix(attached.other_node_id)
            {
                let dir = kps[other_ix] - kps[node_ix];
                let vdiff = kvs[other_ix] - kvs[node_ix];
                if let Ok(spring) = spring.read() {
                    forces += spring.force(&dir, &vdiff);
                }
            }
        }
        forces / node.weight
    }

    pub fn acceleration_par(
        map: &HashMap<NodeId, usize>,
        node_ix: usize,
        node: &Node,
        kps: &[P2],
        kvs: &[V2],
    ) -> V2 {
        let mut forces = V2::ZERO;
        for attached in node.attached.iter() {
            if let Some(spring) = attached.spring.upgrade()
                && let Some(&other_ix) = map.get(&attached.other_node_id)
            {
                let dir = kps[other_ix] - kps[node_ix];
                let vdiff = kvs[other_ix] - kvs[node_ix];
                if let Ok(spring) = spring.read() {
                    forces += spring.force(&dir, &vdiff);
                }
            }
        }
        forces / node.weight
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
            if node.is_stationary() {
                if let Some(mut motor) = node.get_motor() {
                    motor.step(delta_t);
                    node.set_motor(Some(motor));
                }
            }
        }
    }

    pub fn push_node(&mut self, mut node: Node) {
        node.id = self.next_node_id;
        self.node_map.insert(node.id, self.nodes.len());
        self.nodes.push(node);
        self.next_node_id += 1;
    }

    #[inline(always)]
    pub fn get_node_ix(&self, id: NodeId) -> Option<usize> {
        self.node_map.get(&id).copied()
    }

    pub fn rehash_nodes(&mut self) {
        self.node_map.clear();
        for (ix, node) in self.nodes.iter().enumerate() {
            self.node_map.insert(node.id, ix);
        }
    }

    // pub fn get_node_ix(&self, id: NodeId) -> Option<usize> {
    //     if self.nodes.is_empty() {
    //         return None;
    //     }
    //     let min_ix = 0;
    //     let max_ix = id.min(self.nodes.len() - 1);
    //     // println!("max_ix: {}, nodes.len(): {}", max_ix, self.nodes.len());
    //     self.get_node_ix_min_max(id, min_ix, max_ix)
    // }

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
        pretension_force: Float,
    ) {
        let ixa = self.get_node_ix(node_id_a);
        let ixb = self.get_node_ix(node_id_b);
        if let (Some(ixa), Some(ixb)) = (ixa, ixb) {
            let node_distance = self.nodes[ixa].position.distance(self.nodes[ixb].position);
            self.create_spring(node_distance, spring_constant, dampening, pretension_force);
            let aa = Attached {
                spring: Arc::downgrade(self.springs.last().unwrap()),
                other_node_id: self.nodes[ixb].id,
            };
            self.nodes[ixa].attached.push(aa);
            let ab = Attached {
                spring: Arc::downgrade(self.springs.last().unwrap()),
                other_node_id: self.nodes[ixa].id,
            };
            self.nodes[ixb].attached.push(ab);
        }
    }

    pub fn insert_node(&mut self, position: P2, radius: Float) {
        let node = Node::new(0, position, radius);
        self.push_node(node);
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
        spring_pretension_force: Float,
    ) {
        let node_distance = self.nodes[ixa].position.distance(self.nodes[ixb].position);
        self.create_spring(
            node_distance,
            spring_constant,
            spring_dampening,
            spring_pretension_force,
        );
        let ab = Attached {
            spring: Arc::downgrade(self.springs.last().unwrap()),
            other_node_id: self.nodes[ixb].id,
        };
        self.nodes[ixa].attached.push(ab);
        let ba = Attached {
            spring: Arc::downgrade(self.springs.last().unwrap()),
            other_node_id: self.nodes[ixa].id,
        };
        self.nodes[ixb].attached.push(ba);
    }

    pub fn fully_connect(
        &mut self,
        node_ids: &[NodeId],
        spring_constant: Float,
        dampening: Float,
        min_distance: Float,
        max_distance: Float,
    ) {
        for ida in node_ids.iter().copied() {
            for idb in node_ids.iter().copied() {
                if let Some(ixa) = self.get_node_ix(ida) {
                    if let Some(ixb) = self.get_node_ix(idb) {
                        let node_dist =
                            (self.nodes[ixa].position() - self.nodes[ixb].position()).length();
                        let distance_constraint =
                            (min_distance..=max_distance).contains(&node_dist);
                        if ida != idb
                            && !self.spring_exists_between(ida, idb)
                            && distance_constraint
                        {
                            self.insert_spring(ida, idb, spring_constant, dampening, 0.0);
                        }
                    }
                }
            }
        }
    }

    pub fn apply_resting_length_factors(&self) {
        for spring in self.springs.iter() {
            if let Ok(mut spring_write) = spring.write() {
                spring_write.apply_resting_length_factor();
            }
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
            spring_pretension_force,
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
                if ty == GridType::Hex {
                    let pos2 = if (y % 2) == 0 {
                        P2::new(pos.x + node_distance, pos.y)
                    } else {
                        P2::new(pos.x + 2. * node_distance, pos.y)
                    };
                    nodes.push(Node::new(self.next_node_id, pos2, node_radius));
                }
            }
        }
        let base_ix = self.nodes.len();
        for node in nodes {
            self.push_node(node);
        }
        // link the nodes
        for y in 0..height {
            for x in 0..width {
                // Hex grids are fully linked in the match statement bellow
                if ty != GridType::Hex {
                    // link horizontally
                    if x < (width - 1) {
                        let ixa = base_ix + x + y * width;
                        let ixb = ixa + 1;
                        self.insert_spring_on_indices(
                            ixa,
                            ixb,
                            spring_constant,
                            spring_dampening,
                            spring_pretension_force,
                        );
                    }
                    // link vertically
                    if y < (height - 1) {
                        let ixa = base_ix + x + y * width;
                        let ixb = ixa + width;
                        self.insert_spring_on_indices(
                            ixa,
                            ixb,
                            spring_constant,
                            spring_dampening,
                            spring_pretension_force,
                        );
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
                                spring_pretension_force,
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
                                spring_pretension_force,
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
                                spring_pretension_force,
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
                                spring_pretension_force,
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
                                spring_pretension_force,
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
                                    spring_pretension_force,
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
                                spring_pretension_force,
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
                                spring_pretension_force,
                            );
                        }
                    }
                }
            }
        }
    }

    pub fn add_circle(
        &mut self,
        center: P2,
        radius: Float,
        num_nodes: usize,
        node_config: NodeConfig,
    ) {
        for i in 0..num_nodes {
            let div = 2.0 * crate::PI * i as Float / num_nodes as Float;
            self.push_node(node_config.place(P2::new(div.sin(), div.cos()) * radius + center, 0));
        }
    }

    pub fn merge(&mut self, other: &mut Contraption) {
        for node in other.nodes.drain(..) {
            self.push_node(node);
        }
        for spring in other.springs.drain(..) {
            self.springs.push(spring);
        }
    }

    pub fn bounds_collision(&mut self, width: Float, height: Float) {
        for node in self.nodes.iter_mut() {
            node.bounds_collision(width, height);
        }
    }

    pub fn nearest_node_ix(&self, pos: &P2) -> Option<(usize, Float)> {
        let mut res = None;
        for (ix, node) in self.nodes.iter().enumerate() {
            let new_dist = pos.distance(node.position());
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
    pub fn nearest_spring(&self, to: &P2) -> Option<(Weak<RwLock<Spring>>, Float)> {
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
    ) -> (Vec<NodeId>, Vec<Weak<RwLock<Spring>>>) {
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
        pretension_force: Float,
    ) {
        let spring = Spring {
            spring_constant,
            resting_length,
            resting_length_factor: 1.0,
            pretension_force,
            damping_factor,
            min_len: 0.2 * resting_length,
            max_len: 5.0 * resting_length,
        };
        self.springs.push(Arc::new(RwLock::new(spring)));
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
            pretension_force: 0.0,
            damping_factor,
            min_len: 0.25 * resting_length,
            max_len: 4.0 * resting_length,
        };
        let spring_ref = Arc::new(RwLock::new(spring));
        self.springs.push(spring_ref.clone());
        let ata = Attached {
            spring: Arc::downgrade(&spring_ref),
            other_node_id: self.nodes[ixb].id,
        };
        self.nodes[ixa].attached.push(ata);
        let atb = Attached {
            spring: Arc::downgrade(&spring_ref),
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

    /// ma : previous mouse position, mb : current mouse position
    pub fn cut_springs(&mut self, ma: P2, mb: P2) {
        let mut springs_to_cut = Vec::new();
        for node in self.nodes.iter() {
            let nap = node.position;
            for attached in node.attached.iter() {
                if let Some(ixb) = self.get_node_ix(attached.other_node_id) {
                    let nbp = self.nodes[ixb].position;
                    let va = ma - nap;
                    let vb = mb - nap;
                    let vc = nbp - ma;
                    let vs = nbp - nap;
                    let vm = mb - ma;
                    if (vs.perp_dot(va) * vs.perp_dot(vb)) <= 0.0
                        && (vm.perp_dot(va) * vm.perp_dot(vc)) >= 0.0
                    {
                        springs_to_cut.push(attached.spring.clone());
                    }
                }
            }
        }
        springs_to_cut.sort_unstable_by_key(|s| s.as_ptr());
        springs_to_cut.dedup_by_key(|s| s.as_ptr());
        self.springs.retain(|spring_rc| {
            springs_to_cut
                .binary_search_by_key(&Arc::downgrade(spring_rc).as_ptr(), |weak| weak.as_ptr())
                .is_err()
        });
    }

    pub fn delete_node(&mut self, node_id: NodeId) {
        if let Some(ix) = self.get_node_ix(node_id) {
            self.nodes.remove(ix);
        }
    }

    pub fn delete_spring(&mut self, spring_rc: Weak<RwLock<Spring>>) {
        let mut remove_ix = None;
        for (ix, spring) in self.springs.iter().enumerate() {
            let downgraded = Arc::downgrade(spring);
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
        self.node_map.clear();
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
                        if let Ok(spring) = spring.read() {
                            energy += 0.25
                                * (spring.resting_length - distance).powi(2)
                                * spring.spring_constant;
                        }
                    }
                }
            }
        }
        energy
    }

    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub fn nodes_iter<'a>(&'a self) -> std::slice::Iter<'a, Node> {
        self.nodes.iter()
    }

    pub fn nodes_iter_mut<'a>(&'a mut self) -> std::slice::IterMut<'a, Node> {
        self.nodes.iter_mut()
    }

    pub fn get_node<'a>(&'a self, ix: usize) -> &'a Node {
        &self.nodes[ix]
    }

    pub fn get_node_mut<'a>(&'a mut self, ix: usize) -> &'a mut Node {
        &mut self.nodes[ix]
    }
}

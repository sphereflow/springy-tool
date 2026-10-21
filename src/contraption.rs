use macroquad::math::Rect;
use rayon::prelude::*;
use std::{cmp::Ordering, collections::HashMap};

use crate::{
    Float, P2, V2,
    configs::{GridConfig, GridType, NodeConfig},
    node::{Node, NodeId},
    spring::{Spring, SpringId},
};

// if the structure of the contraption is changed meaning adding/deleting a node/spring
// node_map, springs, spring_map and node_to_springs must be updated accordingly
// also next_node_id and next_spring_id need to be adjusted if a node/spring is added or when merging
// with other contraptions
// essentially the first 7 variables of the struct must be kept in sync
#[derive(Debug, Clone)]
pub struct Contraption {
    nodes: Vec<Node>,
    node_map: HashMap<NodeId, usize>,
    next_node_id: NodeId,
    springs: Vec<Spring>,
    spring_map: HashMap<SpringId, usize>,
    next_spring_id: SpringId,
    node_to_springs: HashMap<NodeId, Vec<SpringId>>,
    // intermediate values for Runge Kutta
    // stored in the struct to avoid reallocations
    kps: [Vec<P2>; 4],
    kvs: [Vec<V2>; 4],
    kaccs: [Vec<V2>; 4],
    pub g: V2,
    thread_local_accs: Vec<Vec<V2>>,
}

impl Default for Contraption {
    fn default() -> Self {
        Self::new()
    }
}

impl Contraption {
    pub fn new() -> Contraption {
        let mut thread_local_accs = Vec::new();
        for _ in 0..rayon::current_num_threads() {
            thread_local_accs.push(Vec::new());
        }
        Contraption {
            nodes: Vec::new(),
            node_map: HashMap::new(),
            next_node_id: 0,
            springs: Vec::new(),
            spring_map: HashMap::new(),
            next_spring_id: 0,
            node_to_springs: HashMap::new(),
            kps: [const { Vec::new() }; 4],
            kvs: [const { Vec::new() }; 4],
            kaccs: [const { Vec::new() }; 4],
            thread_local_accs,
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

            Contraption::acceleration_phase(
                &self.node_map,
                &self.springs,
                &self.nodes,
                &self.kps[n + 1],
                &self.kvs[n + 1],
                &mut self.kaccs[n + 1],
            );
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

        // use the positions and velocities of self.nodes which where collected in the previous loop
        self.kaccs[0].fill(V2::ZERO);
        Contraption::acceleration_phase(
            &self.node_map,
            &self.springs,
            &self.nodes,
            &self.kps[0],
            &self.kvs[0],
            &mut self.kaccs[0],
        );
        // transfer calculated accelerations and add gravity
        for i in 0..self.nodes.len() {
            self.nodes[i].acceleration = self.g + self.kaccs[0][i];
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

            Contraption::acceleration_phase_rayon(
                &self.node_map,
                &self.springs,
                &self.nodes,
                &self.kps[n + 1],
                &self.kvs[n + 1],
                kaccsnp1,
                &mut self.thread_local_accs,
            );
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

        // use the positions and velocities of self.nodes which where collected in the previous loop
        Contraption::acceleration_phase_rayon(
            &self.node_map,
            &self.springs,
            &self.nodes,
            &self.kps[0],
            &self.kvs[0],
            &mut self.kaccs[0],
            &mut self.thread_local_accs,
        );
        // transfer calculated accelerations and add gravity
        for i in 0..self.nodes.len() {
            self.nodes[i].acceleration = self.g + self.kaccs[0][i];
        }
    }

    /// calculate all accelerations in parallel
    /// needs positions and velocities and a target buffer to write to
    pub fn acceleration_phase(
        map: &HashMap<NodeId, usize>,
        springs: &[Spring],
        nodes: &[Node],
        kps: &[P2],
        kvs: &[V2],
        acc_target_buffer: &mut [V2],
    ) {
        for spring in springs {
            let ixa = map[&spring.node_id_a];
            let ixb = map[&spring.node_id_b];
            let dir = kps[ixb] - kps[ixa];
            let vdiff = kvs[ixb] - kvs[ixa];
            let force = spring.force(&dir, &vdiff);
            // rk4_init initializes n + 1 accelerations to 0 so adding onto them is fine
            acc_target_buffer[ixa] += force / nodes[ixa].weight;
            acc_target_buffer[ixb] += -force / nodes[ixb].weight;
        }
    }

    /// calculate all accelerations in parallel
    /// needs positions and velocities and a target buffer to write to
    pub fn acceleration_phase_rayon(
        map: &HashMap<NodeId, usize>,
        springs: &[Spring],
        nodes: &[Node],
        kps: &[P2],
        kvs: &[V2],
        acc_target_buffer: &mut [V2],
        thread_local_accs: &mut Vec<Vec<V2>>,
    ) {
        // zero the thread_local_accs buffer
        thread_local_accs
            .par_iter_mut()
            .for_each(|tla| tla.fill(V2::ZERO));
        acc_target_buffer
            .par_iter_mut()
            .for_each(|acc| *acc = V2::ZERO);
        let num_threads = thread_local_accs.len();
        // do not let chunk_size get zero
        let chunk_size = if springs.len() < num_threads {
            springs.len().max(1)
        } else {
            springs.len().div_ceil(num_threads)
        };
        springs
            .par_chunks(chunk_size)
            .zip(&mut *thread_local_accs)
            .for_each(|(springs, tl_accs)| {
                for spring in springs {
                    let ixa = map[&spring.node_id_a];
                    let ixb = map[&spring.node_id_b];
                    let dir = kps[ixb] - kps[ixa];
                    let vdiff = kvs[ixb] - kvs[ixa];
                    let force = spring.force(&dir, &vdiff);
                    // rk4_init initializes n + 1 accelerations to 0 so adding onto them is fine
                    tl_accs[ixa] += force / nodes[ixa].weight;
                    tl_accs[ixb] += -force / nodes[ixb].weight;
                }
            });
        for (ix, acc) in acc_target_buffer.iter_mut().enumerate() {
            for tl_acc in thread_local_accs.iter() {
                *acc += tl_acc[ix];
            }
        }
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
        let num_threads = rayon::current_num_threads();
        for thread_ix in 0..num_threads {
            while self.thread_local_accs[thread_ix].len() < self.nodes.len() {
                self.thread_local_accs[thread_ix].push(V2::ZERO);
            }
        }
        // fill in initial values
        for (i, node) in self.nodes.iter().enumerate() {
            self.kps[0][i] = node.position();
            self.kvs[0][i] = node.velocity;
            self.kaccs[0][i] = node.acceleration;
        }
        self.kaccs[1].fill(V2::ZERO);
        self.kaccs[2].fill(V2::ZERO);
        self.kaccs[3].fill(V2::ZERO);
    }

    pub fn step_motors(&mut self, delta_t: Float) {
        for node in self.nodes.iter_mut() {
            if node.is_stationary()
                && let Some(mut motor) = node.get_motor()
            {
                motor.step(delta_t);
                node.set_motor(Some(motor));
            }
        }
    }

    pub fn push_node(&mut self, mut node: Node) {
        node.id = self.next_node_id;
        self.node_map.insert(node.id, self.nodes.len());
        self.nodes.push(node);
        self.next_node_id += 1;
    }

    /// NOTE: node_id's must be valid
    pub fn push_spring(&mut self, mut spring: Spring) {
        spring.id = self.next_spring_id;
        self.spring_map.insert(spring.id, self.springs.len());

        // insert the spring into self.node_to_springs
        if let Some(connected_spring_ids) = self.node_to_springs.get_mut(&spring.node_id_a) {
            connected_spring_ids.push(spring.id);
        } else {
            self.node_to_springs
                .insert(spring.node_id_a, vec![spring.id]);
        }
        if let Some(connected_spring_ids) = self.node_to_springs.get_mut(&spring.node_id_b) {
            connected_spring_ids.push(spring.id);
        } else {
            self.node_to_springs
                .insert(spring.node_id_b, vec![spring.id]);
        }

        self.springs.push(spring);
        self.next_spring_id += 1;
    }

    #[inline(always)]
    pub fn get_node_ix(&self, id: NodeId) -> Option<usize> {
        self.node_map.get(&id).copied()
    }

    #[inline(always)]
    pub fn get_spring_ix(&self, id: SpringId) -> Option<usize> {
        self.spring_map.get(&id).copied()
    }

    pub fn rehash_nodes(&mut self) {
        self.node_map.clear();
        for (ix, node) in self.nodes.iter().enumerate() {
            self.node_map.insert(node.id, ix);
        }
    }

    pub fn rehash_springs(&mut self) {
        self.node_map.clear();
        for (ix, spring) in self.springs.iter().enumerate() {
            self.spring_map.insert(spring.id, ix);
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

    /// manual binary search for a NodeId between two indices
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
            let spring = Spring {
                id: self.next_spring_id,
                spring_constant,
                resting_length: node_distance,
                resting_length_factor: 1.0,
                pretension_force,
                damping_factor: dampening,
                min_len: 0.2 * node_distance,
                max_len: 5.0 * node_distance,
                node_id_a,
                node_id_b,
            };
            self.push_spring(spring);
        }
    }

    pub fn insert_node(&mut self, position: P2, radius: Float) {
        let node = Node::new(0, position, radius);
        self.push_node(node);
    }

    pub fn spring_exists_between(&self, node_id_a: NodeId, node_id_b: NodeId) -> bool {
        // first both nodes must exist in self.node_to_springs
        if let Some(springs_a) = self.node_to_springs.get(&node_id_a)
            && let Some(springs_b) = self.node_to_springs.get(&node_id_b)
        {
            if springs_a.is_empty() || springs_b.is_empty() {
                return false;
            }
            // both arrays are sorted so step in 'sync'
            let mut ixsa = 0;
            let mut ixsb = 0;
            loop {
                if ixsa >= springs_a.len() {
                    return false;
                }
                if ixsb >= springs_b.len() {
                    return false;
                }
                match springs_a[ixsa].cmp(&springs_b[ixsb]) {
                    Ordering::Less => ixsa += 1,
                    Ordering::Equal => return true,
                    Ordering::Greater => ixsb += 1,
                }
            }
        } else {
            false
        }
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
        let node_id_a = self.nodes[ixa].id;
        let node_id_b = self.nodes[ixb].id;
        let spring = Spring {
            id: self.next_spring_id,
            spring_constant,
            resting_length: node_distance,
            resting_length_factor: 1.0,
            pretension_force: spring_pretension_force,
            damping_factor: spring_dampening,
            min_len: 0.2 * node_distance,
            max_len: 5.0 * node_distance,
            node_id_a,
            node_id_b,
        };
        self.push_spring(spring);
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
                if let Some(ixa) = self.get_node_ix(ida)
                    && let Some(ixb) = self.get_node_ix(idb)
                {
                    let node_dist =
                        (self.nodes[ixa].position() - self.nodes[ixb].position()).length();
                    let distance_constraint = (min_distance..=max_distance).contains(&node_dist);
                    if ida != idb && !self.spring_exists_between(ida, idb) && distance_constraint {
                        self.insert_spring(ida, idb, spring_constant, dampening, 0.0);
                    }
                }
            }
        }
    }

    pub fn apply_resting_length_factors(&mut self) {
        self.springs.iter_mut().for_each(|spring| {
            spring.apply_resting_length_factor();
        });
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
        // shift NodeId's of other up into a range which is not occupied by selfs NodeId's
        // and do the same for SpringId's
        if let Some(first_node_of_other) = other.nodes.first() {
            // update self.nodes
            let adjust_node_id =
                self.next_node_id.max(first_node_of_other.id) - first_node_of_other.id;
            for mut node in other.nodes.drain(..) {
                // we can not use self.push_node(...) since springs rely on adjust_node_id to be constant
                // so we insert the node and do the bookkeeping manually
                node.id += adjust_node_id;
                self.next_node_id = node.id + 1;
                self.node_map.insert(node.id, self.nodes.len());
                self.nodes.push(node);
            }

            // update self.springs
            for mut spring in other.springs.drain(..) {
                // adjust id's
                spring.node_id_a += adjust_node_id;
                spring.node_id_b += adjust_node_id;

                // do all the needed updates
                // push_spring does the update to spring.id and next_spring_id so we do not need
                // to care about that here
                self.push_spring(spring);
            }
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

    /// returns the spring's id and its distance
    pub fn nearest_spring(&self, to: &P2) -> Option<(SpringId, Float)> {
        let mut res = None;
        for spring in self.springs.iter() {
            let ixa = self.node_map[&spring.node_id_a];
            let ixb = self.node_map[&spring.node_id_b];
            let dirab = self.nodes[ixa].position() - self.nodes[ixb].position();
            if dirab.length() <= f64::EPSILON {
                continue;
            }
            let dirba = -dirab;
            // if to is closer to the spring than any of its nodes
            let a = to - self.nodes[ixa].position();
            let b = to - self.nodes[ixb].position();
            if (dirab.dot(a) > 0.0) && (dirba.dot(b) > 0.0) {
                // calculate the distance
                let dirabn = dirab.normalize();
                let dist = (a.length_squared() - dirabn.dot(a).powi(2)).sqrt();
                if let Some((_, current_dist)) = res {
                    if dist < current_dist {
                        res = Some((spring.id, dist));
                    }
                } else {
                    res = Some((spring.id, dist));
                }
            }
        }
        res
    }

    /// previous_nodes are previously selected nodes which should be part of the
    /// new selection as well
    pub fn multi_select(
        &self,
        rect: &Rect,
        previous_nodes: &[NodeId],
    ) -> (Vec<NodeId>, Vec<SpringId>) {
        let mut node_ids = Vec::from_iter(previous_nodes.iter().copied());
        // collect every node
        for node in &self.nodes {
            let pos =
                macroquad::math::Vec2::new(node.position().x as f32, node.position().y as f32);
            if rect.contains(pos) {
                node_ids.push(node.id);
            }
        }
        node_ids.sort();
        node_ids.dedup();
        let mut spring_ids = Vec::new();
        for spring in &self.springs {
            if node_ids.binary_search(&spring.node_id_a).is_ok()
                && node_ids.binary_search(&spring.node_id_b).is_ok()
            {
                spring_ids.push(spring.id);
            }
        }
        spring_ids.sort();
        spring_ids.dedup();
        (node_ids, spring_ids)
    }

    pub fn create_spring(
        &mut self,
        node_id_a: NodeId,
        node_id_b: NodeId,
        resting_length: Float,
        spring_constant: Float,
        damping_factor: Float,
        pretension_force: Float,
    ) {
        let spring = Spring {
            id: self.next_spring_id,
            spring_constant,
            resting_length,
            resting_length_factor: 1.0,
            pretension_force,
            damping_factor,
            min_len: 0.2 * resting_length,
            max_len: 5.0 * resting_length,
            node_id_a,
            node_id_b,
        };
        self.push_spring(spring);
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
            id: self.next_spring_id,
            spring_constant,
            resting_length,
            resting_length_factor: 1.0,
            pretension_force: 0.0,
            damping_factor,
            min_len: 0.25 * resting_length,
            max_len: 4.0 * resting_length,
            node_id_a: self.nodes[ixa].id,
            node_id_b: self.nodes[ixb].id,
        };
        self.push_spring(spring);
    }

    pub fn spring_break(&mut self) {
        let mut delete_spring_indices = Vec::new();
        for (ix, spring) in self.springs.iter().enumerate() {
            if spring.is_destruct(&self.nodes, &self.node_map) {
                delete_spring_indices.push(ix);
            }
        }
        // delete indices in reverse order since delete_springs swap remove could mess with higher indices
        for ix in delete_spring_indices.iter().rev() {
            self.delete_spring(*ix);
        }
    }

    pub fn get_attached_spring_ids(&self, node_id: NodeId) -> Vec<SpringId> {
        self.node_to_springs
            .get(&node_id)
            .unwrap_or(&vec![])
            .clone()
    }

    /// remove springs between ma and mb
    /// ma : previous mouse position, mb : current mouse position
    pub fn cut_springs(&mut self, ma: P2, mb: P2) {
        let mut delete_spring_indices = Vec::new();
        for (ix, spring) in self.springs.iter().enumerate() {
            if spring.crosses(&self.nodes, &self.node_map, ma, mb) {
                delete_spring_indices.push(ix);
            }
        }
        // delete indices in reverse order since delete_springs swap remove could mess with higher indices
        for ix in delete_spring_indices.iter().rev() {
            self.delete_spring(*ix);
        }
    }

    pub fn delete_node(&mut self, node_id: NodeId) {
        if let Some(ix) = self.get_node_ix(node_id) {
            // first remove all attached springs
            for spring_id in self.node_to_springs[&node_id].clone().iter() {
                self.delete_spring(self.spring_map[spring_id]);
            }
            self.nodes.swap_remove(ix);

            self.reassign_node_id(self.nodes[ix].id, node_id, ix);
        }
    }

    fn reassign_node_id(&mut self, from: NodeId, to: NodeId, node_ix: usize) {
        self.nodes[node_ix].id = to;
        if let Some(spring_ids) = self.node_to_springs.remove(&from) {
            for spring_id in spring_ids.iter() {
                if let Some(spring) = self.springs.get_mut(self.spring_map[spring_id]) {
                    if spring.node_id_a == from {
                        spring.node_id_a = to;
                    }
                    if spring.node_id_b == from {
                        spring.node_id_b = to;
                    }
                }
            }
            self.node_to_springs.insert(to, spring_ids);
        }
        if let Some(last_node) = self.nodes.last() {
            self.next_node_id = last_node.id + 1;
        }
    }

    pub fn delete_spring(&mut self, spring_ix: usize) {
        let is_last_elem = spring_ix == (self.springs.len() - 1);
        let removed = self.springs.swap_remove(spring_ix);
        // remove deleted SpringId from node_to_springs
        if let Ok(node_to_spring_ix) =
            self.node_to_springs[&removed.node_id_a].binary_search(&removed.id)
            && let Some(spring_ids) = self.node_to_springs.get_mut(&removed.node_id_a)
        {
            spring_ids.remove(node_to_spring_ix);
        }
        if let Ok(node_to_spring_ix) =
            self.node_to_springs[&removed.node_id_b].binary_search(&removed.id)
            && let Some(spring_ids) = self.node_to_springs.get_mut(&removed.node_id_b)
        {
            spring_ids.remove(node_to_spring_ix);
        }
        if is_last_elem {
            return;
        }

        self.reassign_spring_id(self.springs[spring_ix].id, removed.id, spring_ix);
    }

    /// preserve order of SpringId's in self.springs
    /// this assumes that the SpringId 'to' is not yet taken
    /// this also assumes that after reassigning the order of SpringId's in self.springs will be
    /// sorted after the operation as it does not sort self.springs
    fn reassign_spring_id(&mut self, from: SpringId, to: SpringId, spring_ix: usize) {
        let node_id_a = self.springs[spring_ix].node_id_a;
        if let Ok(node_to_spring_ix) = self.node_to_springs[&node_id_a].binary_search(&from)
            && let Some(spring_ids) = self.node_to_springs.get_mut(&node_id_a)
        {
            spring_ids[node_to_spring_ix] = to;
            spring_ids.sort_unstable();
        }
        let node_id_b = self.springs[spring_ix].node_id_b;
        if let Ok(node_to_spring_ix) = self.node_to_springs[&node_id_b].binary_search(&from)
            && let Some(spring_ids) = self.node_to_springs.get_mut(&node_id_b)
        {
            spring_ids[node_to_spring_ix] = to;
            spring_ids.sort_unstable();
        }
        self.springs[spring_ix].id = to;
        self.next_spring_id = self.springs.last().map_or(0, |spring| spring.id + 1);
    }

    pub fn clear(&mut self) {
        self.springs.clear();
        self.nodes.clear();
        self.node_map.clear();
        self.spring_map.clear();
        self.node_to_springs.clear();
        self.next_node_id = 0;
        self.next_spring_id = 0;
    }

    // non mutating functions

    pub fn energy(&self) -> Float {
        let mut energy = 0.0;
        // nodes kinetic energy
        for node in &self.nodes {
            energy += 0.5 * node.weight * node.velocity.length_squared();
        }
        // springs potential energy
        for spring in &self.springs {
            let ixa = self.node_map[&spring.node_id_a];
            let ixb = self.node_map[&spring.node_id_b];
            let distance = (self.nodes[ixa].position - self.nodes[ixb].position).length();
            energy += 0.5 * (spring.resting_length - distance).powi(2) * spring.spring_constant;
        }
        energy
    }

    pub fn check_indices_equal_ids(&self) -> bool {
        for (ix, node) in self.nodes.iter().enumerate() {
            if ix != node.id {
                return false;
            }
        }
        for (ix, spring) in self.springs.iter().enumerate() {
            if ix != spring.id {
                return false;
            }
        }
        true
    }

    pub fn check_id_integrity(&self) -> bool {
        for (ix, node) in self.nodes.iter().enumerate() {
            if self.node_map[&node.id] != ix {
                return false;
            }
        }
        for (ix, spring) in self.springs.iter().enumerate() {
            if self.spring_map[&spring.id] != ix {
                return false;
            }
            if !self.node_map.contains_key(&spring.node_id_a) {
                return false;
            }
            if !self.node_map.contains_key(&spring.node_id_b) {
                return false;
            }
        }
        for (&node_id, spring_ids) in &self.node_to_springs {
            for spring_id in spring_ids.iter() {
                let spring = self.springs[self.spring_map[spring_id]];
                if (spring.node_id_a != node_id) && (spring.node_id_b != node_id) {
                    return false;
                }
            }
        }
        true
    }

    // Spring and Node accessors

    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub fn nodes_iter<'a>(&'a self) -> std::slice::Iter<'a, Node> {
        self.nodes.iter()
    }

    pub fn nodes_iter_mut<'a>(&'a mut self) -> std::slice::IterMut<'a, Node> {
        self.nodes.iter_mut()
    }

    pub fn get_node(&self, ix: usize) -> &Node {
        &self.nodes[ix]
    }

    pub fn get_node_mut(&mut self, ix: usize) -> &mut Node {
        &mut self.nodes[ix]
    }

    pub fn num_springs(&self) -> usize {
        self.springs.len()
    }

    pub fn springs_iter<'a>(&'a self) -> std::slice::Iter<'a, Spring> {
        self.springs.iter()
    }

    pub fn springs_iter_mut<'a>(&'a mut self) -> std::slice::IterMut<'a, Spring> {
        self.springs.iter_mut()
    }

    pub fn get_spring(&self, ix: usize) -> &Spring {
        &self.springs[ix]
    }

    pub fn get_spring_mut(&mut self, ix: usize) -> &mut Spring {
        &mut self.springs[ix]
    }
}

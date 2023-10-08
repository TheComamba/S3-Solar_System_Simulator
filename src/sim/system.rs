use crate::{
    sim::body::Body,
    sim::initial_parameters::{Float, InitialParameters, DIMENSIONALITY, G},
};

pub(crate) struct StellarSystem {
    current_time: Float,
    bodies: Vec<Body>,
}

impl StellarSystem {
    pub(crate) fn new(params: InitialParameters) -> StellarSystem {
        let mut bodies = vec![Body::new(
            0,
            0.,
            0.,
            params.total_mass * params.stellar_mass_fraction,
        )];
        for i in 0..params.body_count {
            let mass = params.total_mass * (1. - params.stellar_mass_fraction)
                / params.body_count as Float;
            bodies.push(Body::new(
                i + 1,
                params.position_variance,
                params.velocity_variance,
                mass,
            ));
        }
        StellarSystem {
            current_time: 0.,
            bodies,
        }
    }

    fn get_acceleration(body1: &Body, body2: &Body) -> Vec<Float> {
        let mut r = vec![0.; DIMENSIONALITY];
        for i in 0..DIMENSIONALITY {
            r[i] = body1.position[i] - body2.position[i];
        }
        let r_squared = r.iter().map(|x| x * x).sum::<Float>();
        let r_cubed = r_squared * r_squared.sqrt();
        let r_unit = r
            .iter()
            .map(|x| x / r_squared.sqrt())
            .collect::<Vec<Float>>();
        r_unit
            .iter()
            .map(|x| -G * body1.mass * body2.mass * x / r_cubed)
            .collect::<Vec<Float>>()
    }

    fn do_collisions(&mut self, time_step: Float) {
        let mut bodies_to_merge = vec![];
        for i in 0..self.bodies.len() {
            for j in (i + 1)..self.bodies.len() {
                let body1 = &self.bodies[i];
                let body2 = &self.bodies[j];
                if body1.collides_with(&body2, time_step) {
                    bodies_to_merge.push((body1, body2));
                }
            }
        }

        if !bodies_to_merge.is_empty() {
            let mut new_bodies = self.bodies.clone();
            for (body1, body2) in bodies_to_merge.iter() {
                let mut new_body = (*body1).clone();
                new_body.merge_with(&body2);
                new_bodies.retain(|x| x.index != body1.index && x.index != body2.index);
                new_bodies.push(new_body);
            }
            self.bodies = new_bodies;
        }
    }

    pub(crate) fn evolve(&mut self, time_step: Float) {
        self.do_collisions(time_step);

        let mut accelerations = vec![vec![0.; DIMENSIONALITY]; self.bodies.len()];
        for i in 0..self.bodies.len() {
            for j in 0..self.bodies.len() {
                if i == j {
                    continue;
                }
                let acceleration = Self::get_acceleration(&self.bodies[i], &self.bodies[j]);
                for k in 0..DIMENSIONALITY {
                    accelerations[i][k] += acceleration[k];
                }
            }
        }
        for i in 0..self.bodies.len() {
            for k in 0..DIMENSIONALITY {
                self.bodies[i].position[k] += self.bodies[i].velocity[k] * time_step;
                self.bodies[i].velocity[k] += accelerations[i][k] * time_step;
            }
        }
        self.current_time += time_step;
    }
}

use crate::{
    body::Body,
    initial_parameters::{Float, InitialParameters, DIMENSIONALITY, G},
};

pub(crate) struct StellarSystem {
    current_time: Float,
    bodies: Vec<Body>,
}

impl StellarSystem {
    pub(crate) fn new(params: InitialParameters) -> StellarSystem {
        let mut bodies = vec![Body::new(
            0.,
            0.,
            params.total_mass * params.stellar_mass_fraction,
        )];
        for _ in 0..params.body_count {
            let mass = params.total_mass * (1. - params.stellar_mass_fraction)
                / params.body_count as Float;
            bodies.push(Body::new(
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

    pub(crate) fn evolve(&mut self, time_step: Float) {
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

use crate::{
    sim::body::Body,
    sim::initial_parameters::{Float, InitialParameters, DIMENSIONALITY, G},
};

pub(crate) const MIN_TIMESTEP: Float = 1e-5;

#[derive(Clone, Debug)]
pub(crate) struct StellarSystem {
    pub(crate) current_time: Float,
    pub(crate) bodies: Vec<Body>,
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

    fn get_acceleration(accelerated: &Body, accelerating: &Body) -> Vec<Float> {
        let mut r = vec![0.; DIMENSIONALITY];
        for i in 0..DIMENSIONALITY {
            r[i] = accelerated.position[i] - accelerating.position[i];
        }
        let r_squared = r.iter().map(|x| x * x).sum::<Float>();
        let r_cubed = r_squared * r_squared.sqrt();
        let r_unit = r
            .iter()
            .map(|x| x / r_squared.sqrt())
            .collect::<Vec<Float>>();
        r_unit
            .iter()
            .map(|x| -G * accelerating.mass * x / r_cubed)
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

    /*
       r > v * t
       r / v > t
       r^2 / v^2 > t^2
    */
    fn get_timestep(&self, max: Float) -> Float {
        let mut time_step_sqrd = max * max;
        for i in 0..self.bodies.len() {
            for j in (i + 1)..self.bodies.len() {
                let body1 = &self.bodies[i];
                let body2 = &self.bodies[j];

                assert_eq!(DIMENSIONALITY, 2);
                let r = body1
                    .position
                    .iter()
                    .zip(body2.position.iter())
                    .map(|(x, y)| (x - y).abs())
                    .collect::<Vec<Float>>();
                let v = body1
                    .velocity
                    .iter()
                    .zip(body2.velocity.iter())
                    .map(|(x, y)| (x - y).abs())
                    .collect::<Vec<Float>>();
                let r_sqrd = r.iter().map(|x| x * x).sum::<Float>();
                let mut v_sqrd = v.iter().map(|x| x * x).sum::<Float>();
                if v_sqrd < MIN_TIMESTEP * MIN_TIMESTEP {
                    v_sqrd = MIN_TIMESTEP * MIN_TIMESTEP;
                }
                let candidate = r_sqrd / v_sqrd;
                if candidate < time_step_sqrd {
                    time_step_sqrd = candidate;
                }
            }
        }
        if time_step_sqrd < MIN_TIMESTEP * MIN_TIMESTEP {
            MIN_TIMESTEP
        } else {
            time_step_sqrd.sqrt()
        }
    }

    fn do_evolution_step(&mut self, time_step: Float) {
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

    pub(crate) fn evolve_for(&mut self, time: Float) {
        let target_time = self.current_time + time;
        while self.current_time < target_time {
            let time_step = self.get_timestep(target_time - self.current_time);
            self.do_collisions(time_step);
            self.do_evolution_step(time_step);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn acceleration_is_antisymmetric() {
        let body1 = Body {
            position: vec![0., 0.],
            velocity: vec![0., 0.],
            mass: 1.,
            index: 1,
        };
        let body2 = Body {
            position: vec![1., 1.],
            velocity: vec![0., 0.],
            mass: 1.,
            index: 2,
        };
        let acc_1 = StellarSystem::get_acceleration(&body1, &body2);
        let acc_2 = StellarSystem::get_acceleration(&body2, &body1);
        println!("Acceleration of body 1:\n{:?}", acc_1);
        println!("Acceleration of body 2:\n{:?}", acc_2);
        assert!((acc_1[0] + acc_2[0]).abs() < 1e-5);
        assert!((acc_1[1] + acc_2[1]).abs() < 1e-5);
    }

    #[test]
    fn acceleration_depends_on_inertia() {
        let body1 = Body {
            position: vec![0., 0.],
            velocity: vec![0., 0.],
            mass: 1e5,
            index: 1,
        };
        let body2 = Body {
            position: vec![1., 1.],
            velocity: vec![0., 0.],
            mass: 1e-5,
            index: 2,
        };
        let acc_large_body = StellarSystem::get_acceleration(&body1, &body2);
        let acc_small_body = StellarSystem::get_acceleration(&body2, &body1);
        println!("Acceleration of large body:\n{:?}", acc_large_body);
        println!("Acceleration of small body:\n{:?}", acc_small_body);
        assert!(acc_large_body[0].abs() < 1e-5);
        assert!(acc_large_body[1].abs() < 1e-5);
        assert!(acc_small_body[0].abs() > 1e-5);
        assert!(acc_small_body[1].abs() > 1e-5);
    }

    #[test]
    fn symmetric_bodies_end_up_resting() {
        const TIME_STEP: Float = 1e1;

        let position1 = vec![1., 0.];
        let position2 = vec![-1., 0.];
        let velocity1 = vec![0., 0.];
        let velocity2 = vec![0., 0.];
        let body1 = Body {
            position: position1,
            velocity: velocity1,
            mass: 1.,
            index: 1,
        };
        let body2 = Body {
            position: position2,
            velocity: velocity2,
            mass: 1.,
            index: 2,
        };
        let mut system = StellarSystem {
            current_time: 0.,
            bodies: vec![body1, body2],
        };

        system.evolve_for(TIME_STEP);
        system.evolve_for(TIME_STEP);
        println!("{:?}", system);

        assert!(system.bodies[0].position[0] < 1.);
        assert!(system.bodies[1].position[0] > -1.);
        assert!(system.bodies[0].velocity[0] < 0.);
        assert!(system.bodies[1].velocity[0] > 0.);

        let mut loop_count = 0;
        while system.bodies.len() > 1 && loop_count < 10_000 {
            system.evolve_for(TIME_STEP);
            loop_count += 1;
        }
        println!("Looped {} times.", loop_count);
        println!("{:?}", system);

        assert!(system.bodies.len() == 1);
        assert!(system.bodies[0].position[0].abs() < 1e-5);
        assert!(system.bodies[0].position[1].abs() < 1e-5);
        assert!(system.bodies[0].velocity[0].abs() < 1e-5);
        assert!(system.bodies[0].velocity[1].abs() < 1e-5);
    }

    #[test]
    fn small_body_falls_onto_large_body() {
        const TIME_STEP: Float = 1e1;
        let body1 = Body {
            position: vec![0., 0.],
            velocity: vec![0., 0.],
            mass: 1e5,
            index: 1,
        };
        let body2 = Body {
            position: vec![1., 1.],
            velocity: vec![0., 0.],
            mass: 1e-5,
            index: 2,
        };
        let mut system = StellarSystem {
            current_time: 0.,
            bodies: vec![body1, body2],
        };
        let mut loop_count = 0;
        while system.bodies.len() > 1 && loop_count < 10_000 {
            system.evolve_for(TIME_STEP);
            loop_count += 1;
        }
        println!("Looped {} times.", loop_count);
        println!("{:?}", system);

        assert!(system.bodies.len() == 1);
        assert!(system.bodies[0].position[0].abs() < 1e-5);
        assert!(system.bodies[0].position[1].abs() < 1e-5);
        assert!(system.bodies[0].velocity[0].abs() < 1e-5);
        assert!(system.bodies[0].velocity[1].abs() < 1e-5);
    }

    #[test]
    fn two_body_system_is_stable() {
        const TIME_STEP: Float = 1.;
        let values = vec![-1., 1., 1e2];
        for x in values.iter() {
            for y in values.iter() {
                for v_x in values.iter() {
                    for v_y in values.iter() {
                        println!("x = {}, y = {}, v_x = {}, v_y = {}", x, y, v_x, v_y);

                        let body1 = Body {
                            position: vec![0., 0.],
                            velocity: vec![0., 0.],
                            mass: 1.,
                            index: 1,
                        };
                        let body2 = Body {
                            position: vec![*x, *y],
                            velocity: vec![*v_x, *v_y],
                            mass: 1.,
                            index: 2,
                        };
                        let initial_sam = body1.specific_relative_angular_momentum(&body2);
                        if initial_sam.abs() < 1e-5 {
                            continue;
                        }
                        let mut system = StellarSystem {
                            current_time: 0.,
                            bodies: vec![body1, body2],
                        };

                        println!("{:?}", system);
                        println!("Initial specific angular momentum: {}", initial_sam);
                        for i in 0..10_000 {
                            system.evolve_for(TIME_STEP);

                            println!("Loop {}", i);
                            println!("{:?}", system);
                            assert!(system.bodies.len() == 2);

                            let new_sam = system.bodies[0]
                                .specific_relative_angular_momentum(&system.bodies[1]);

                            println!("Specific angular momentum: {}", new_sam);
                            assert!((new_sam - initial_sam).abs() < 1e-3 * initial_sam.abs());
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn bodies_at_same_position_collide() {
        let v_x_values = vec![-1., 0., 1., 1e5];
        let v_y_values = v_x_values.clone();
        for v_x in v_x_values.iter() {
            for v_y in v_y_values.iter() {
                println!("v_x = {}, v_y = {}", v_x, v_y);
                const TIME_STEP: Float = 1e1;
                let body1 = Body {
                    position: vec![0., 0.],
                    velocity: vec![*v_x, *v_y],
                    mass: 1.,
                    index: 1,
                };
                let body2 = Body {
                    position: vec![0., 0.],
                    velocity: vec![-v_x, -v_y],
                    mass: 1.,
                    index: 2,
                };
                let mut system = StellarSystem {
                    current_time: 0.,
                    bodies: vec![body1, body2],
                };
                system.evolve_for(TIME_STEP);
                println!("{:?}", system);

                assert!(system.bodies.len() == 1);
                assert!(system.bodies[0].position[0].abs() < 1e-5);
                assert!(system.bodies[0].position[1].abs() < 1e-5);
                assert!(system.bodies[0].velocity[0].abs() < 1e-5);
                assert!(system.bodies[0].velocity[1].abs() < 1e-5);
            }
        }
    }

    #[test]
    fn timestep_for_bodies_at_same_position_is_small() {
        let values = vec![-1., 0., 1., 1e5];
        for v_x_1 in values.iter() {
            for v_y_1 in values.iter() {
                for v_x_2 in values.iter() {
                    for v_y_2 in values.iter() {
                        println!(
                            "v_x_1 = {}, v_y_1 = {}, v_x_2 = {}, v_y_2 = {}",
                            v_x_1, v_y_1, v_x_2, v_y_2
                        );
                        const TIME_STEP: Float = 1e1;
                        let body1 = Body {
                            position: vec![0., 0.],
                            velocity: vec![*v_x_1, *v_y_1],
                            mass: 1.,
                            index: 1,
                        };
                        let body2 = Body {
                            position: vec![0., 0.],
                            velocity: vec![*v_x_2, *v_y_2],
                            mass: 1.,
                            index: 2,
                        };
                        let system = StellarSystem {
                            current_time: 0.,
                            bodies: vec![body1, body2],
                        };
                        assert!(system.get_timestep(TIME_STEP) < 1.1 * MIN_TIMESTEP);
                    }
                }
            }
        }
    }
}

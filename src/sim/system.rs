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
        r.iter()
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

    fn leap_frog_kick(&mut self, time_step_half: Float, accelerations: Vec<Vec<Float>>) {
        for i in 0..self.bodies.len() {
            for k in 0..DIMENSIONALITY {
                self.bodies[i].velocity[k] += accelerations[i][k] * time_step_half;
            }
        }
    }

    fn leap_frog_drift(&mut self, time_step: Float) {
        for i in 0..self.bodies.len() {
            for k in 0..DIMENSIONALITY {
                self.bodies[i].position[k] += self.bodies[i].velocity[k] * time_step;
            }
        }
    }

    fn get_all_accelerations(&mut self) -> Vec<Vec<f32>> {
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
        accelerations
    }

    fn do_evolution_step(&mut self, time_step: Float) {
        let accelerations = self.get_all_accelerations();
        self.leap_frog_kick(0.5 * time_step, accelerations);
        self.leap_frog_drift(time_step);
        let accelerations = self.get_all_accelerations();
        self.leap_frog_kick(0.5 * time_step, accelerations);
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

    #[cfg(test)]
    fn total_energy(&self) -> Float {
        let mut total_energy = 0.;
        for i in 0..self.bodies.len() {
            total_energy += self.bodies[i].kinetic_energy();
            for j in (i + 1)..self.bodies.len() {
                total_energy += self.bodies[i].relative_potential_energy(&self.bodies[j]);
            }
        }
        total_energy
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
    fn acceleration_is_isotropic() {
        let mut acceleration_modulus = None;
        for dir in 0..1 {
            for sign in vec![-1., 1.] {
                let body1 = Body {
                    position: vec![0., 0.],
                    velocity: vec![0., 0.],
                    mass: 1.,
                    index: 1,
                };
                let mut pos = vec![0., 0.];
                pos[dir] = sign;
                println!("Position: {:?}", pos);
                let body2 = Body {
                    position: pos,
                    velocity: vec![0., 0.],
                    mass: 1.,
                    index: 2,
                };
                let new_acc = StellarSystem::get_acceleration(&body1, &body2)
                    .iter()
                    .map(|x| x * x)
                    .sum::<Float>()
                    .sqrt();
                match acceleration_modulus {
                    None => acceleration_modulus = Some(new_acc),
                    Some(old_acc) => {
                        println!("Old acceleration:\n{:?}", new_acc);
                        println!("New acceleration:\n{:?}", old_acc);
                        assert!((new_acc - old_acc).abs() < 1e-5);
                    }
                }
            }
        }
    }

    #[test]
    fn acceleration_is_translation_invariant() {
        let values = [-1., 0., 1., 1e5];
        let diff = vec![1., -3.];
        let mut initial_acceleration = None;
        for x in values.iter() {
            for y in values.iter() {
                println!("x = {}, y = {}", x, y);
                let body1 = Body {
                    position: vec![*x, *y],
                    velocity: vec![0., 0.],
                    mass: 1.,
                    index: 1,
                };
                let body2 = Body {
                    position: vec![*x + diff[0], *y + diff[1]],
                    velocity: vec![0., 0.],
                    mass: 1.,
                    index: 2,
                };
                let acceleration = StellarSystem::get_acceleration(&body1, &body2);
                match &initial_acceleration {
                    None => initial_acceleration = Some(acceleration.clone()),
                    Some(old_acc) => {
                        println!("Initial acceleration:\n{:?}", old_acc);
                        println!("New acceleration:\n{:?}", acceleration);
                        assert!((old_acc[0] - acceleration[0]).abs() < 1e-5);
                        assert!((old_acc[1] - acceleration[1]).abs() < 1e-5);
                    }
                }
            }
        }
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
    fn acceleration_vanishes_over_disctance() {
        let body1 = Body {
            position: vec![-1e5, -1e5],
            velocity: vec![0., 1.],
            mass: 1.,
            index: 1,
        };
        let body2 = Body {
            position: vec![-1e5, 1e5],
            velocity: vec![2., -3.],
            mass: 1.,
            index: 2,
        };
        let body3 = Body {
            position: vec![1e5, -1e5],
            velocity: vec![4., 5e5],
            mass: 1.,
            index: 3,
        };
        let body4 = Body {
            position: vec![1e5, 1e5],
            velocity: vec![6e6, 7.],
            mass: 1.,
            index: 4,
        };
        let bodies = vec![body1, body2, body3, body4];
        for body1 in bodies.iter() {
            for body2 in bodies.iter() {
                if body1.index == body2.index {
                    continue;
                }
                let acc_1 = StellarSystem::get_acceleration(&body1, &body2);
                let acc_2 = StellarSystem::get_acceleration(&body2, &body1);
                println!("Acceleration of body 1:\n{:?}", acc_1);
                println!("Acceleration of body 2:\n{:?}", acc_2);
                assert!(acc_1[0].abs() < 1e-5);
                assert!(acc_1[1].abs() < 1e-5);
                assert!(acc_2[0].abs() < 1e-5);
                assert!(acc_2[1].abs() < 1e-5);
            }
        }
    }

    /*
       F = G * m1 * m2 / r^2
    */
    #[test]
    fn unit_masses_at_unit_distance_accelerate_with_g() {
        let body1 = Body {
            position: vec![0., 0.],
            velocity: vec![0., 0.],
            mass: 1.,
            index: 1,
        };
        let body2 = Body {
            position: vec![1., 0.],
            velocity: vec![0., 0.],
            mass: 1.,
            index: 2,
        };
        let acc_1 = StellarSystem::get_acceleration(&body1, &body2)[0].abs();
        let acc_2 = StellarSystem::get_acceleration(&body2, &body1)[0].abs();
        println!("Acceleration of body 1: {}", acc_1);
        println!("Acceleration of body 2: {}", acc_2);
        print!("G: {}", G);
        assert!((acc_1 - G).abs() < 1e-5);
        assert!((acc_2 - G).abs() < 1e-5);
    }

    /*
       a = G * M / r^2
    */
    #[test]
    fn acceleration_is_proportional_to_well_mass() {
        let mut well = Body {
            position: vec![0., 0.],
            velocity: vec![0., 0.],
            mass: 1.,
            index: 1,
        };
        let test = Body {
            position: vec![1., 0.],
            velocity: vec![0., 0.],
            mass: 1.,
            index: 2,
        };
        let acc_1 = StellarSystem::get_acceleration(&test, &well)[0].abs();
        well.mass = 2.;
        let acc_2 = StellarSystem::get_acceleration(&test, &well)[0].abs();
        println!("Acceleration 1: {}", acc_1);
        println!("Acceleration 2: {}", acc_2);
        assert!((2. * acc_1 - acc_2).abs() < 1e-5);
    }

    /*
       a = G * M / r^2
    */
    #[test]
    fn acceleration_is_independent_of_test_mass() {
        let well = Body {
            position: vec![0., 0.],
            velocity: vec![0., 0.],
            mass: 1.,
            index: 1,
        };
        let mut test = Body {
            position: vec![1., 0.],
            velocity: vec![0., 0.],
            mass: 1.,
            index: 2,
        };
        let acc_1 = StellarSystem::get_acceleration(&test, &well)[0].abs();
        test.mass = 2.;
        let acc_2 = StellarSystem::get_acceleration(&test, &well)[0].abs();
        println!("Acceleration 1: {}", acc_1);
        println!("Acceleration 2: {}", acc_2);
        assert!((acc_1 - acc_2).abs() < 1e-5);
    }

    /*
       a = G * M / r^2
    */
    #[test]
    fn acceleration_is_proportional_to_one_over_distance_squared() {
        let well = Body {
            position: vec![0., 0.],
            velocity: vec![0., 0.],
            mass: 1.,
            index: 1,
        };
        let mut test = Body {
            position: vec![1., 0.],
            velocity: vec![0., 0.],
            mass: 1.,
            index: 2,
        };
        let acc_1 = StellarSystem::get_acceleration(&test, &well)[0].abs();
        test.position[0] = 2.;
        let acc_2 = StellarSystem::get_acceleration(&test, &well)[0].abs();
        println!("Acceleration 1: {}", acc_1);
        println!("Acceleration 2: {}", acc_2);
        assert!((acc_1 - 4. * acc_2).abs() < 1e-5);
    }

    #[test]
    fn force_is_negative_gradient_of_potential_energy() {
        let masses = vec![1., 1e2, 1e3];
        for m1 in masses.iter() {
            for m2 in masses.iter() {
                println!("m1 = {}, m2 = {}", m1, m2);
                const D: Float = 1e-5;
                let body1 = Body {
                    position: vec![0., 0.],
                    velocity: vec![0., 0.],
                    mass: *m1,
                    index: 1,
                };
                let mut body2 = Body {
                    position: vec![1., 0.],
                    velocity: vec![0., 0.],
                    mass: *m2,
                    index: 2,
                };
                let force_1_from_acc =
                    body1.mass * StellarSystem::get_acceleration(&body1, &body2)[0].abs();
                let force_2_from_acc =
                    body2.mass * StellarSystem::get_acceleration(&body2, &body1)[0].abs();
                let pot_1 = body1.relative_potential_energy(&body2);
                body2.position[0] += D;
                println!("{:?}", body2);
                let pot_2 = body1.relative_potential_energy(&body2);
                let force_from_pot = ((pot_1 - pot_2) / D).abs();

                println!("Force from acceleration 1: {}", force_1_from_acc);
                println!("Force from acceleration 2: {}", force_2_from_acc);
                println!("Force from potential: {}", force_from_pot);
                assert!((force_1_from_acc - force_from_pot).abs() < 1e-2 * force_1_from_acc);
                assert!((force_2_from_acc - force_from_pot).abs() < 1e-2 * force_from_pot);
            }
        }
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

    //#[test]
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

    #[test]
    fn near_miss_is_irrelevant_for_large_body() {
        const TIME_STEP: Float = 1e1;
        let body1 = Body {
            position: vec![0., 0.],
            velocity: vec![0., 0.],
            mass: 1e5,
            index: 1,
        };
        let body2 = Body {
            position: vec![0.1, -10.],
            velocity: vec![0., 10.],
            mass: 1e-5,
            index: 2,
        };
        let mut system = StellarSystem {
            current_time: 0.,
            bodies: vec![body1, body2],
        };
        let mut loop_count = 0;
        while (system.bodies[1].position[0].powi(2) + system.bodies[1].position[1].powi(2)).sqrt()
            < 1e3
            && loop_count < 10_000
        {
            system.evolve_for(TIME_STEP);
            loop_count += 1;
            assert!(system.bodies.len() == 2);
        }
        println!("Looped {} times.", loop_count);
        println!("{:?}", system);

        assert!(system.bodies[0].position[0].abs() < 1e-5);
        assert!(system.bodies[0].position[1].abs() < 1e-5);
        assert!(system.bodies[0].velocity[0].abs() < 1e-5);
        assert!(system.bodies[0].velocity[1].abs() < 1e-5);
    }

    #[test]
    fn near_miss_does_not_change_total_momenta() {
        const TIME_STEP: Float = 1e3;
        let body1 = Body {
            position: vec![-1., -1e3],
            velocity: vec![0., 1e2],
            mass: 1.,
            index: 1,
        };
        let body2 = Body {
            position: vec![1., 1e3],
            velocity: vec![0., -1e2],
            mass: 1.,
            index: 2,
        };
        let initial_sam = body1.specific_relative_angular_momentum(&body2);
        let initial_v_1 = (body1.velocity[0].powi(2) + body1.velocity[1].powi(2)).sqrt();
        let initial_v_2 = (body2.velocity[0].powi(2) + body2.velocity[1].powi(2)).sqrt();
        let mut system = StellarSystem {
            current_time: 0.,
            bodies: vec![body1, body2],
        };

        let mut loop_count = 0;
        while (system.bodies[1].position[0].powi(2) + system.bodies[1].position[1].powi(2)).sqrt()
            < 1e3
            && loop_count < 10_000
        {
            system.evolve_for(TIME_STEP);
            loop_count += 1;
            assert!(system.bodies.len() == 2);
        }
        println!("Looped {} times.", loop_count);
        println!("{:?}", system);

        let new_sam = system.bodies[0].specific_relative_angular_momentum(&system.bodies[1]);
        let new_v1 =
            (system.bodies[0].velocity[0].powi(2) + system.bodies[0].velocity[1].powi(2)).sqrt();
        let new_v2 =
            (system.bodies[1].velocity[0].powi(2) + system.bodies[1].velocity[1].powi(2)).sqrt();
        assert!((new_sam - initial_sam).abs() < 1e-5);
        assert!((new_v1 - initial_v_1).abs() < 1e-5);
        assert!((new_v2 - initial_v_2).abs() < 1e-5);
    }

    #[test]
    fn bound_three_body_system_remains_bound() {
        const TIME_STEP: Float = 1e0;
        let body1 = Body {
            position: vec![0., 0.],
            velocity: vec![0., 0.],
            mass: 3.,
            index: 1,
        };
        let body2 = Body {
            position: vec![1., 0.],
            velocity: vec![0., 0.01],
            mass: 2.,
            index: 2,
        };
        let body3 = Body {
            position: vec![2., 0.],
            velocity: vec![0., 0.001],
            mass: 1.,
            index: 3,
        };
        let mut system = StellarSystem {
            current_time: 0.,
            bodies: vec![body1, body2, body3],
        };

        let initial_total_energy = system.total_energy();
        println!("Initial total energy: {}", initial_total_energy);

        for _ in 0..100 {
            system.evolve_for(TIME_STEP);
        }
        println!("{:?}", system);
        let new_total_energy = system.total_energy();
        println!("New total energy: {}", new_total_energy);

        assert!(system.bodies.len() == 3);
        assert!(initial_total_energy < 0.);
        assert!(new_total_energy < 0.);
        assert!((new_total_energy - initial_total_energy).abs() < 1e-5);
    }

    #[test]
    fn escaping_a_massive_body_conserves_total_energy() {
        const TIME_STEP: Float = 1e1;
        let body1 = Body {
            position: vec![0., 0.],
            velocity: vec![0., 0.],
            mass: 1e5,
            index: 1,
        };
        let body2 = Body {
            position: vec![0.1, 0.1],
            velocity: vec![0., 20.],
            mass: 1.,
            index: 2,
        };
        let mut system = StellarSystem {
            current_time: 0.,
            bodies: vec![body1, body2],
        };
        let initial_potential_energy =
            system.bodies[1].relative_potential_energy(&system.bodies[0]);
        let initial_kinetic_energy1 = system.bodies[0].kinetic_energy();
        let initial_kinetic_energy2 = system.bodies[1].kinetic_energy();
        let initial_total_energy =
            initial_potential_energy + initial_kinetic_energy1 + initial_kinetic_energy2;
        println!("Initial potential: {}", initial_potential_energy);
        println!("Initial kinetic 1: {}", initial_kinetic_energy1);
        println!("Initial kinetic 2: {}", initial_kinetic_energy2);
        println!("Initial total: {}", initial_total_energy);
        assert!(initial_total_energy > 0.); //body escapes

        let mut loop_count = 0;
        while system.bodies[1].position[1].abs() < 1e3 && loop_count < 10_000 {
            system.evolve_for(TIME_STEP);
            loop_count += 1;
        }
        println!("Looped {} times.", loop_count);
        println!("{:?}", system);

        let new_potential_energy = system.bodies[1].relative_potential_energy(&system.bodies[0]);
        let new_kinetic_energy_1 = system.bodies[0].kinetic_energy();
        let new_kinetic_energy_2 = system.bodies[1].kinetic_energy();
        let new_total_energy = new_potential_energy + new_kinetic_energy_1 + new_kinetic_energy_2;
        println!("New potential: {}", new_potential_energy);
        println!("New kinetic energy 1: {}", new_kinetic_energy_2);
        println!("New kinetic energy 2: {}", new_kinetic_energy_1);
        println!("New total: {}", new_total_energy);

        assert!(new_potential_energy.abs() < 1e-5);
        assert!(new_potential_energy > initial_potential_energy);
        assert!(new_kinetic_energy_2 < initial_kinetic_energy2);
        assert!((initial_total_energy - new_total_energy).abs() < 1e-5);

        assert!(system.bodies[1].velocity[1] < 1.);
    }
}

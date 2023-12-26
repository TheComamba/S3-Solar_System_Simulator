use super::units::{Float, DIMENSIONALITY, G};
use crate::{sim::body::Body, sim::initial_parameters::InitialParameters};

pub(crate) const MIN_TIMESTEP: Float = 1e-8;

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
        Condition 1: Distance does not change too much
        r d > v t
        => r d / v > t
        => t^2 < r^2 d^2 / v^2

        Condition 2: Gravitational field is roughly homogenous
        d a(0) > |a(0) - a(t)|
        => d G M / r^2 > G M / r^2 - G M / (r + v t)^2
        => d > 1 - r^2 / (r + v t)^2
        => r^2 / (r + v t)^2 > 1 - d
        => (r + v t)^2 < r^2 / (1 - d)
        => (1 + v t / r)^2 < 1 / (1 - d)
        t is small, d is small
        => 1 + 2 v t / r < 1 + d
        => t < d r / (2 v)
        => t^2 < d^2 r^2 / (4 v^2) = d'^2 r^2 / v^2
        => Condition 2 implies Condition 1

        Condition 3: Velocity change is small
        d v(0) > |v(0) - v(t)|
        => d v > |v - v - a t| = G M t / r^2
        => t < d v r^2 / (G M)
        => t^2 < v^2 r^4 d^2 / (G M)^2
    */
    fn get_timestep(&self, max: Float) -> Float {
        const D_SQUARED: Float = 1e-8;
        assert_eq!(DIMENSIONALITY, 2);
        let mut time_step_sqrd = max * max;
        for i in 0..self.bodies.len() {
            let d_over_g_m_squared = D_SQUARED / (G * self.bodies[i].mass).powi(2);
            for j in (i + 1)..self.bodies.len() {
                let body1 = &self.bodies[i];
                let body2 = &self.bodies[j];

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
                let r_sqrd = r.iter().map(|x| x.powi(2)).sum::<Float>();
                let mut v_sqrd = v.iter().map(|x| x.powi(2)).sum::<Float>();
                if v_sqrd < MIN_TIMESTEP * MIN_TIMESTEP {
                    v_sqrd = MIN_TIMESTEP * MIN_TIMESTEP;
                }
                let candidate = D_SQUARED * r_sqrd / v_sqrd;
                if candidate < time_step_sqrd {
                    time_step_sqrd = candidate;
                }
                let candidate = v_sqrd * r_sqrd.powi(2) * d_over_g_m_squared;
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
            // println!("Time step: {}", time_step);
            self.do_collisions(time_step);
            self.do_evolution_step(time_step);
        }
    }

    #[cfg(test)]
    fn total_energy(&self) -> Float {
        let mut total_energy = 0.;
        for i in 0..self.bodies.len() {
            total_energy += self.bodies[i].kinetic_energy();
            for j in 0..self.bodies.len() {
                if i == j {
                    continue;
                }
                total_energy += self.bodies[i].relative_potential_energy(&self.bodies[j]);
            }
        }
        total_energy
    }
}

#[cfg(test)]
mod tests {
    use crate::sim::{
        body,
        units::{DISTANCE_TO_M, PI, TIME_TO_SECONDS, VELOCITY_TO_KM_PER_S},
    };

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

    /*
        Initially: v = 0, r = R, m = M
        T = T_1 + T_2 = 0
        V_0 = - 2 G M m / R
        Finally: r = R/2
        => V_f = - 4 G M m / R
        => T_f = V_0 - V_f = 2 G M m / R
               = T_1 + T_2 = 2 T_1 = 2 T_2
        => T_1 = 1/2 m v^2 = T_f / 2 = G M m / R
        => v_1 = sqrt(2 G M / R)
        Passed time:
        t = sqrt(R^3 / 2 G M) (sqrt(1/2 (1 - 1/2)) + arccos sqrt(1/2)) = sqrt(R^3 / 2 G M) (1/2 + Pi/4)
        R = 1., M = 1.
        Time until collision:
        sqrt(R^3 / 2 G M^2) * Pi/2
    */
    //#[test]
    fn relative_velocity_of_two_bodies_falling_half_their_distance() {
        const ACC: Float = 1e-4;
        let values: Vec<Float> = vec![1e0, 1e1, 1e2];
        for r in values.clone() {
            for m1 in values.clone() {
                for m2 in values.clone() {
                    let t = (r.powi(3) / (2. * G * (m1 + m2))).sqrt() * (0.5 + PI / 4.) as Float;
                    let expected_velocity1 = (2. * G * m2 / r).sqrt();
                    let expected_velocity2 = (2. * G * m1 / r).sqrt();

                    println!("r = {}, m1 = {}, m2 = {}", r, m1, m2);
                    println!("Time for half the fall: {}", t);
                    println!("Expected velocity 1: {}", expected_velocity1);
                    println!("Expected velocity 2: {}", expected_velocity2);

                    let body1 = Body {
                        position: vec![-r / 2., 0.],
                        velocity: vec![0., 0.],
                        mass: m1,
                        index: 1,
                    };
                    let body2 = Body {
                        position: vec![r / 2., 0.],
                        velocity: vec![0., 0.],
                        mass: m2,
                        index: 2,
                    };
                    let mut system = StellarSystem {
                        current_time: 0.,
                        bodies: vec![body1, body2],
                    };
                    let initial_potential_energy =
                        2. * system.bodies[0].relative_potential_energy(&system.bodies[1]);
                    let initial_total_energy = system.total_energy();

                    system.evolve_for(t);

                    println!("{:?}", system);

                    let relative_distance =
                        (system.bodies[0].position[0] - system.bodies[1].position[0]).abs();
                    println!("Relative distance: {}", relative_distance);

                    let final_potential_energy =
                        2. * system.bodies[0].relative_potential_energy(&system.bodies[1]);
                    let final_kinetic_energy =
                        system.bodies[0].kinetic_energy() + system.bodies[1].kinetic_energy();
                    let final_total_energy = system.total_energy();
                    println!("Initial total energy: {}", initial_total_energy);
                    println!("Final total energy: {}", final_total_energy);
                    println!("Kinetic energy: {}", final_kinetic_energy);
                    println!(
                        "Potential energy difference: {}",
                        initial_potential_energy - final_potential_energy
                    );

                    let velocity1 = system.bodies[0].velocity[0].abs();
                    let velocity2 = system.bodies[1].velocity[0].abs();
                    println!("Velocity 1: {}", velocity1);
                    println!("Velocity 2: {}", velocity2);

                    assert!(system.bodies[0].position[1].abs() < ACC);
                    assert!(system.bodies[1].position[1].abs() < ACC);
                    assert!(system.bodies[0].velocity[1].abs() < ACC);
                    assert!(system.bodies[1].velocity[1].abs() < ACC);
                    assert!((relative_distance - r / 2.).abs() < ACC);
                    assert!(
                        (initial_total_energy - final_total_energy).abs()
                            < ACC * initial_total_energy
                    );
                    assert!(
                        (final_kinetic_energy + final_potential_energy - initial_potential_energy)
                            .abs()
                            < ACC * initial_total_energy
                    );
                    assert!((velocity1 - expected_velocity1).abs() < ACC * velocity1);
                    assert!((velocity2 - expected_velocity2).abs() < ACC * velocity2);
                }
            }
        }
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
        const ACC: Float = 1e-3;
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

        assert!(new_potential_energy.abs() < ACC * initial_total_energy.abs());
        assert!(new_potential_energy > initial_potential_energy);
        assert!(new_kinetic_energy_2 < initial_kinetic_energy2);
        assert!((initial_total_energy - new_total_energy).abs() < ACC * initial_total_energy.abs());
    }

    /*
        r_0_1 = 0
        r_0_2 = R
        a_0_1 = G m_2 / R^2
        a_0_2 = - G m_1 / R^2
        v_half_1 = a_0_1 * t / 2
        v_half_2 = a_0_2 * t / 2
        r_1_1 = v_half_1 * t
        r_1_2 = R - v_half_2 * t
        a_1_1 = G m_2 / (r_1_2-r_1_1)^2
        a_1_2 = - G m_1 / (r_1_2-r_1_1)^2
        v_1_1 = v_half_1 + a_1_1 * t / 2
        v_1_2 = v_half_2 + a_1_2 * t / 2
    */
    #[test]
    fn leap_frog_for_simple_case() {
        let values: Vec<Float> = vec![1e0, 1e1, 1e2];
        for m1 in values.clone() {
            for m2 in values.clone() {
                for r in values.clone() {
                    for t in values.clone() {
                        println!("m1 = {}, m2 = {}, r = {}, t = {}", m1, m2, r, t);

                        let a_0_1 = G * m2 / r.powi(2);
                        let a_0_2 = -G * m1 / r.powi(2);
                        let v_half_1 = a_0_1 * t / 2.;
                        let v_half_2 = a_0_2 * t / 2.;
                        let r_1_1 = v_half_1 * t;
                        let r_1_2 = r + v_half_2 * t;
                        let acc_sign = if r_1_2 > r_1_1 { 1. } else { -1. };
                        let a_1_1 = acc_sign * G * m2 / (r_1_2 - r_1_1).powi(2);
                        let a_1_2 = -acc_sign * G * m1 / (r_1_2 - r_1_1).powi(2);
                        let v_1_1 = v_half_1 + a_1_1 * t / 2.;
                        let v_1_2 = v_half_2 + a_1_2 * t / 2.;

                        println!("a_0_1 = {}", a_0_1);
                        println!("a_0_2 = {}", a_0_2);
                        println!("v_half_1 = {}", v_half_1);
                        println!("v_half_2 = {}", v_half_2);
                        println!("r_1_1 = {}", r_1_1);
                        println!("r_1_2 = {}", r_1_2);
                        println!("a_1_1 = {}", a_1_1);
                        println!("a_1_2 = {}", a_1_2);
                        println!("v_1_1 = {}", v_1_1);
                        println!("v_1_2 = {}", v_1_2);

                        let body1 = Body {
                            position: vec![0., 0.],
                            velocity: vec![0., 0.],
                            mass: m1,
                            index: 1,
                        };
                        let body2 = Body {
                            position: vec![r, 0.],
                            velocity: vec![0., 0.],
                            mass: m2,
                            index: 2,
                        };
                        let mut system = StellarSystem {
                            current_time: 0.,
                            bodies: vec![body1, body2],
                        };

                        system.do_evolution_step(t);
                        println!("{:?}", system);

                        assert!((system.bodies[0].position[0] - r_1_1).abs() < 1e-5 * r_1_1.abs());
                        assert!((system.bodies[1].position[0] - r_1_2).abs() < 1e-5 * r_1_2.abs());
                        assert!((system.bodies[0].velocity[0] - v_1_1).abs() < 1e-5 * v_1_1.abs());
                        assert!((system.bodies[1].velocity[0] - v_1_2).abs() < 1e-5 * v_1_2.abs());
                    }
                }
            }
        }
    }

    #[test]
    fn sun_earth_system() {
        let sun_mass = 333_000.;
        let earth_mass = 1.;
        let au: Float = 1.;
        let earth_orbital_speed = 2. * PI;

        let accuracy = 1e-3;
        let time_step = 0.25;

        let sun = Body {
            position: vec![0., 0.],
            velocity: vec![0., 0.],
            mass: sun_mass,
            index: 1,
        };
        let earth = Body {
            position: vec![au, 0.],
            velocity: vec![0., earth_orbital_speed],
            mass: earth_mass,
            index: 2,
        };
        //TODO: this is not working
        //earth.set_velocity_for_circular_orbit(&sun);
        let mut system = StellarSystem {
            current_time: 0.,
            bodies: vec![sun, earth],
        };

        println!("\n{:?}\n", system);
        assert!(system.bodies[1].velocity[0].abs() < accuracy);
        assert!((system.bodies[1].velocity[1] - earth_orbital_speed).abs() < accuracy);
        assert!((system.bodies[1].position[0] - au).abs() < accuracy);
        assert!(system.bodies[1].position[1].abs() < accuracy);

        system.evolve_for(time_step);
        println!("\n{:?}\n", system);
        assert!(system.bodies[1].position[0].abs() < accuracy);
        assert!((system.bodies[1].position[1] - au).abs() < accuracy);

        system.evolve_for(time_step);
        println!("\n{:?}\n", system);
        assert!((system.bodies[1].position[0] + au).abs() < accuracy);
        assert!(system.bodies[1].position[1].abs() < accuracy);

        system.evolve_for(time_step);
        println!("\n{:?}\n", system);
        assert!(system.bodies[1].position[0].abs() < accuracy);
        assert!((system.bodies[1].position[1] + au).abs() < accuracy);
    }

    #[test]
    fn earth_moon_system() {
        let earth_mass = 1.;
        let moon_mass = 0.012;
        let distance: Float = 384_748_000. / DISTANCE_TO_M;
        let period = 0.074855;
        let d_moon_barycenter = distance * earth_mass / (moon_mass + earth_mass);
        let moon_orbital_speed = 2. * PI * d_moon_barycenter / period;
        let d_earth_barycenter = distance * moon_mass / (moon_mass + earth_mass);
        let earth_orbital_speed = 2. * PI * d_earth_barycenter / period;
        let time_step = 0.25 * period;
        println!("Period = {} years", period);
        println!(
            "d_moon_barycenter = {} km",
            d_moon_barycenter * DISTANCE_TO_M / 1000.
        );
        println!(
            "d_earth_barycenter = {} km",
            d_earth_barycenter * DISTANCE_TO_M / 1000.
        );
        println!(
            "moon_orbital_speed = {} km/s",
            moon_orbital_speed * VELOCITY_TO_KM_PER_S
        );
        println!(
            "earth_orbital_speed = {} km/s",
            earth_orbital_speed * VELOCITY_TO_KM_PER_S
        );

        let accuracy = 3e-5;

        let earth = Body {
            position: vec![-d_earth_barycenter, 0.],
            velocity: vec![0., -earth_orbital_speed],
            mass: earth_mass,
            index: 1,
        };
        let moon = Body {
            position: vec![d_moon_barycenter, 0.],
            velocity: vec![0., moon_orbital_speed],
            mass: moon_mass,
            index: 2,
        };
        //TODO: this is not working
        //earth.set_velocity_for_circular_orbit(&sun);
        let mut system = StellarSystem {
            current_time: 0.,
            bodies: vec![earth, moon],
        };

        println!("\n{:?}\n", system);
        assert!((system.bodies[0].position[0] + d_earth_barycenter).abs() < accuracy);
        assert!(system.bodies[0].position[1].abs() < accuracy);
        assert!((system.bodies[1].position[0] - d_moon_barycenter).abs() < accuracy);
        assert!(system.bodies[1].position[1].abs() < accuracy);

        system.evolve_for(time_step);
        println!("\n{:?}\n", system);
        assert!(system.bodies[0].position[0].abs() < accuracy);
        assert!((system.bodies[0].position[1] + d_earth_barycenter).abs() < accuracy);
        assert!(system.bodies[1].position[0].abs() < accuracy);
        assert!((system.bodies[1].position[1] - d_moon_barycenter).abs() < accuracy);

        system.evolve_for(time_step);
        println!("\n{:?}\n", system);
        assert!((system.bodies[0].position[0] - d_earth_barycenter).abs() < accuracy);
        assert!(system.bodies[0].position[1].abs() < accuracy);
        assert!((system.bodies[1].position[0] + d_moon_barycenter).abs() < accuracy);
        assert!(system.bodies[1].position[1].abs() < accuracy);

        system.evolve_for(time_step);
        println!("\n{:?}\n", system);
        assert!(system.bodies[0].position[0].abs() < accuracy);
        assert!((system.bodies[0].position[1] - d_earth_barycenter).abs() < accuracy);
        assert!(system.bodies[1].position[0].abs() < accuracy);
        assert!((system.bodies[1].position[1] + d_moon_barycenter).abs() < accuracy);
    }

    //#[test]
    fn specific_stability_test() {
        let mut system = StellarSystem {
            current_time: 110.0,
            bodies: vec![
                Body {
                    index: 0,
                    position: vec![-0.16409463, 0.31138736],
                    velocity: vec![0.20223364, -0.07811384],
                    mass: 0.5,
                },
                Body {
                    index: 1,
                    position: vec![-0.14044371, 0.3014494],
                    velocity: vec![-0.40218723, 0.15970851],
                    mass: 0.25,
                },
                Body {
                    index: 2,
                    position: vec![0.50843185, 1.4344311],
                    velocity: vec![-0.002280068, -0.0034808216],
                    mass: 0.25,
                },
            ],
        };
        println!("\n{:?}\n", system);

        let initial_energy = system.total_energy();
        println!("Initial energy: {}", initial_energy);
        for _ in 0..100 {
            system.evolve_for(0.1);

            println!("\n{:?}\n", system);
            let new_energy = system.total_energy();
            println!("New energy: {}", new_energy);
            assert!((new_energy - initial_energy).abs() < 1e-5 * initial_energy.abs());
        }
    }

    fn other_specific_test_case() {
        let mut system = StellarSystem {
            current_time: 60.0,
            bodies: vec![
                Body {
                    index: 0,
                    position: vec![-0.15938611, -0.040578824],
                    velocity: vec![-0.008027485, -0.0020192228],
                    mass: 0.5,
                },
                Body {
                    index: 1,
                    position: vec![-1.7621905, -1.2192367],
                    velocity: vec![0.0011315977, 0.00090261595],
                    mass: 0.25,
                },
                Body {
                    index: 2,
                    position: vec![-0.37285408, -0.09300329],
                    velocity: vec![0.01492337, 0.0031358302],
                    mass: 0.25,
                },
            ],
        };
    }
}

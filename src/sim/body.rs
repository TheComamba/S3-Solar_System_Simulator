use super::system::MIN_TIMESTEP;
use super::units::{Float, DIMENSIONALITY, G};
use crate::sim::units::{PI, ROCK_DENSITY};
use rand_distr::{Distribution, Normal};

#[derive(Clone, Debug)]
pub(crate) struct Body {
    pub(crate) index: u32,
    pub(crate) position: Vec<Float>,
    pub(crate) velocity: Vec<Float>,
    pub(crate) mass: Float,
}

impl PartialEq for Body {
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index
    }
}

impl Body {
    fn random_vector(variance: Float) -> Vec<Float> {
        let mut vector = vec![0.; DIMENSIONALITY];
        let distribution = Normal::new(0., variance).unwrap();
        for i in 0..DIMENSIONALITY {
            vector[i] = distribution.sample(&mut rand::thread_rng());
        }
        vector
    }

    pub(crate) fn new(
        index: u32,
        position_variance: Float,
        velocity_variance: Float,
        mass: Float,
    ) -> Body {
        Body {
            index,
            position: Self::random_vector(position_variance),
            velocity: Self::random_vector(velocity_variance),
            mass,
        }
    }

    pub(crate) fn radius(&self) -> Float {
        const MASS_TO_VOLUME_FACTOR: Float = 3. / (4. * PI * ROCK_DENSITY);
        (self.mass * MASS_TO_VOLUME_FACTOR).cbrt()
    }

    fn comes_close(&self, other: &Self, time_step: Float) -> bool {
        let mut relative_position = vec![0.; DIMENSIONALITY];
        for i in 0..DIMENSIONALITY {
            relative_position[i] = self.position[i] - other.position[i];
        }
        let mut relative_velocity = vec![0.; DIMENSIONALITY];
        for i in 0..DIMENSIONALITY {
            relative_velocity[i] = self.velocity[i] - other.velocity[i];
        }
        let distance_squared = relative_position.iter().map(|x| x * x).sum::<Float>();
        let relative_speed_squared = relative_velocity.iter().map(|x| x * x).sum::<Float>();
        2. * relative_speed_squared.sqrt() * time_step + MIN_TIMESTEP > distance_squared.sqrt()
    }

    pub(crate) fn specific_relative_angular_momentum(&self, other: &Self) -> Float {
        let mut relative_position = vec![0.; DIMENSIONALITY];
        for i in 0..DIMENSIONALITY {
            relative_position[i] = self.position[i] - other.position[i];
        }
        let mut relative_velocity = vec![0.; DIMENSIONALITY];
        for i in 0..DIMENSIONALITY {
            relative_velocity[i] = self.velocity[i] - other.velocity[i];
        }
        assert!(DIMENSIONALITY == 2);
        relative_position[0] * relative_velocity[1] - relative_position[1] * relative_velocity[0]
    }

    fn gravitational_parameter(&self, other: &Self) -> Float {
        G * (self.mass + other.mass)
    }

    fn minimal_distance(&self, other: &Self) -> Float {
        let specific_relative_angular_momentum = self.specific_relative_angular_momentum(other);
        let gravitational_parameter = self.gravitational_parameter(other);
        specific_relative_angular_momentum * specific_relative_angular_momentum
            / gravitational_parameter
    }

    fn will_collide(&self, other: &Self) -> bool {
        self.minimal_distance(other) < self.radius() + other.radius()
    }

    pub(crate) fn collides_with(&self, other: &Self, time_step: Float) -> bool {
        self.comes_close(other, time_step) && self.will_collide(other)
    }

    pub(crate) fn merge_with(&mut self, other: &Self) {
        let total_mass = self.mass + other.mass;
        for i in 0..DIMENSIONALITY {
            self.position[i] =
                (self.position[i] * self.mass + other.position[i] * other.mass) / total_mass;
            self.velocity[i] =
                (self.velocity[i] * self.mass + other.velocity[i] * other.mass) / total_mass;
        }
        self.mass = total_mass;
    }

    #[cfg(test)]
    pub(crate) fn kinetic_energy(&self) -> Float {
        let relative_speed_squared = self.velocity.iter().map(|x| x * x).sum::<Float>();
        0.5 * self.mass * relative_speed_squared
    }

    #[cfg(test)]
    pub(crate) fn relative_potential_energy(&self, other: &Self) -> Float {
        let mut relative_position = vec![0.; DIMENSIONALITY];
        for i in 0..DIMENSIONALITY {
            relative_position[i] = self.position[i] - other.position[i];
        }
        let distance_squared = relative_position.iter().map(|x| x * x).sum::<Float>();
        -G * self.mass * other.mass / distance_squared.sqrt()
    }

    #[cfg(test)]
    pub(crate) fn set_velocity_for_circular_orbit(&mut self, other: &Self) {
        let mut barycenter = vec![0.; DIMENSIONALITY];
        for i in 0..DIMENSIONALITY {
            barycenter[i] = (self.position[i] * self.mass + other.position[i] * other.mass)
                / (self.mass + other.mass);
        }
        let mut position_relative_to_barycenter = vec![0.; DIMENSIONALITY];
        for i in 0..DIMENSIONALITY {
            position_relative_to_barycenter[i] = self.position[i] - barycenter[i];
        }
        let distance_to_barycenter = position_relative_to_barycenter
            .iter()
            .map(|x| x.powi(2))
            .sum::<Float>()
            .sqrt();
        let reduced_mass = self.mass * other.mass / (self.mass + other.mass);
        let orbital_speed = (G * reduced_mass / distance_to_barycenter).sqrt();
        let mut tangential_velocity = vec![0.; DIMENSIONALITY];
        tangential_velocity[0] =
            -position_relative_to_barycenter[1] * orbital_speed / distance_to_barycenter;
        tangential_velocity[1] =
            position_relative_to_barycenter[0] * orbital_speed / distance_to_barycenter;
        self.velocity = tangential_velocity;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn collide_bodies_with_same_mass_and_opposite_velocities() {
        let position1 = vec![1., 2.];
        let position2 = vec![-1., -2.];
        let velocity1 = vec![3., 4.];
        let velocity2 = vec![-3., -4.];
        let mut body1 = Body {
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

        body1.merge_with(&body2);

        assert!(body1.position[0].abs() < 1e-5);
        assert!(body1.position[1].abs() < 1e-5);
        assert!(body1.velocity[0].abs() < 1e-5);
        assert!(body1.velocity[1].abs() < 1e-5);
    }

    #[test]
    fn collide_large_with_small_body() {
        let position1 = vec![0., 0.];
        let position2 = vec![1., 1.];
        let velocity1 = vec![0., 0.];
        let velocity2 = vec![1., 1.];
        let mut body1 = Body {
            position: position1,
            velocity: velocity1,
            mass: 1e5,
            index: 1,
        };
        let body2 = Body {
            position: position2,
            velocity: velocity2,
            mass: 1e-5,
            index: 2,
        };

        body1.merge_with(&body2);

        assert!(body1.position[0].abs() < 1e-5);
        assert!(body1.position[1].abs() < 1e-5);
        assert!(body1.velocity[0].abs() < 1e-5);
        assert!(body1.velocity[1].abs() < 1e-5);
    }

    #[test]
    fn bodies_at_same_position_come_close_and_will_collide() {
        let values = vec![-1., 0., 1., 1e5];
        for x in values.iter() {
            for y in values.iter() {
                for v_x_1 in values.iter() {
                    for v_y_1 in values.iter() {
                        for v_x_2 in values.iter() {
                            for v_y_2 in values.iter() {
                                println!(
                                    "x: {}, y: {}, v_x_1: {}, v_y_1: {}, v_x_2: {}, v_y_2: {}",
                                    x, y, v_x_1, v_y_1, v_x_2, v_y_2
                                );
                                let body1 = Body {
                                    position: vec![*x, *y],
                                    velocity: vec![*v_x_1, *v_y_1],
                                    mass: 1.,
                                    index: 1,
                                };
                                let body2 = Body {
                                    position: vec![*x, *y],
                                    velocity: vec![*v_x_2, *v_y_2],
                                    mass: 1.,
                                    index: 2,
                                };

                                assert!(body1.comes_close(&body2, 1.));
                                assert!(body1.will_collide(&body2));
                            }
                        }
                    }
                }
            }
        }
    }

    //#[test]
    fn close_but_escaping_bodies_will_not_collide() {
        let values: Vec<Float> = vec![-1., 0., 1., 1e5];
        for v_x_1 in values.iter() {
            for v_y_1 in values.iter() {
                for v_x_2 in values.iter() {
                    for v_y_2 in values.iter() {
                        println!(
                            "v_x_1: {}, v_y_1: {}, v_x_2: {}, v_y_2: {}",
                            v_x_1, v_y_1, v_x_2, v_y_2
                        );
                        if (v_x_1 - v_x_2).abs() < 1e-5 && (v_y_1 - v_y_2).abs() < 1e-5 {
                            continue;
                        }
                        let body1 = Body {
                            position: vec![*v_x_1 / 10., *v_y_1 / 10.],
                            velocity: vec![*v_x_1, *v_y_1],
                            mass: 1.,
                            index: 1,
                        };
                        let body2 = Body {
                            position: vec![*v_x_2 / 10., *v_y_2 / 10.],
                            velocity: vec![*v_x_2, *v_y_2],
                            mass: 1.,
                            index: 2,
                        };
                        println!("body1: {:?}", body1);
                        println!("body2: {:?}", body2);

                        assert!(body1.comes_close(&body2, 1.));
                        assert!(!body1.will_collide(&body2));
                    }
                }
            }
        }
    }

    #[test]
    fn potential_energy_of_unit_masses_is_minus_g() {
        let body1 = Body {
            position: vec![0., 0.],
            velocity: vec![1., 2.],
            mass: 1.,
            index: 1,
        };
        let body2 = Body {
            position: vec![1., 0.],
            velocity: vec![3., 4.],
            mass: 1.,
            index: 2,
        };

        let pot = body1.relative_potential_energy(&body2);
        println!("pot: {}", pot);
        assert!((pot + G).abs() < 1e-5);
    }

    #[test]
    fn potential_energy_is_proportional_to_masses() {
        let mut body1 = Body {
            position: vec![0., 0.],
            velocity: vec![1., 2.],
            mass: 1.,
            index: 1,
        };
        let mut body2 = Body {
            position: vec![1., 0.],
            velocity: vec![3., 4.],
            mass: 1.,
            index: 2,
        };

        let pot1 = body1.relative_potential_energy(&body2);
        body1.mass = 2.;
        body2.mass = 2.;
        let pot2 = body1.relative_potential_energy(&body2);
        println!("pot1: {}", pot1);
        println!("pot2: {}", pot2);
        assert!((4. * pot1 - pot2).abs() < 1e-5);
    }

    #[test]
    fn potential_energy_is_proportional_to_inverse_distance() {
        let body1 = Body {
            position: vec![0., 0.],
            velocity: vec![1., 2.],
            mass: 1.,
            index: 1,
        };
        let mut body2 = Body {
            position: vec![1., 0.],
            velocity: vec![3., 4.],
            mass: 1.,
            index: 2,
        };

        let pot1 = body1.relative_potential_energy(&body2);
        body2.position[0] = 2.;
        let pot2 = body1.relative_potential_energy(&body2);
        println!("pot1: {}", pot1);
        println!("pot2: {}", pot2);
        assert!((pot1 - 2. * pot2).abs() < 1e-5);
    }

    #[test]
    fn kinetic_energy_of_unit_mass_is_one_half() {
        let body = Body {
            position: vec![1e1, 0.],
            velocity: vec![1., 0.],
            mass: 1.,
            index: 2,
        };

        let kin = body.kinetic_energy();
        println!("kin: {}", kin);
        assert!((kin - 0.5).abs() < 1e-5);
    }

    #[test]
    fn kinetic_energy_is_proportional_to_mass() {
        let mut body = Body {
            position: vec![1e1, 0.],
            velocity: vec![1., 0.],
            mass: 1.,
            index: 2,
        };

        let kin1 = body.kinetic_energy();
        body.mass = 2.;
        let kin2 = body.kinetic_energy();
        println!("kin: {}", kin1);
        println!("kin: {}", kin2);
        assert!((2. * kin1 - kin2).abs() < 1e-5);
    }

    #[test]
    fn kinetic_energy_is_proportional_to_velocity_squared() {
        let mut body = Body {
            position: vec![1e1, 0.],
            velocity: vec![1., 0.],
            mass: 1.,
            index: 2,
        };

        let kin1 = body.kinetic_energy();
        body.velocity[0] = 2.;
        let kin2 = body.kinetic_energy();
        println!("kin: {}", kin1);
        println!("kin: {}", kin2);
        assert!((4. * kin1 - kin2).abs() < 1e-5);
    }
}

use crate::initial_parameters::{Float, DIMENSIONALITY, G, ROCK_DENSITY};
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
        const MASS_TO_VOLUME_FACTOR: Float = 3. / (4. * std::f32::consts::PI * ROCK_DENSITY);
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
        relative_speed_squared.sqrt() * time_step > distance_squared.sqrt()
    }

    fn specific_relative_angular_momentum(&self, other: &Self) -> Float {
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
}

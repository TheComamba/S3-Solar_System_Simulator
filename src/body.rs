use crate::initial_parameters::{Float, ROCK_DENSITY};

pub(crate) struct Body {
    position: Vec<Float>,
    velocity: Vec<Float>,
    mass: Float,
}

impl Body {
    pub(crate) fn new(position: Vec<Float>, velocity: Vec<Float>, mass: Float) -> Body {
        Body {
            position,
            velocity,
            mass,
        }
    }

    pub(crate) fn radius(&self) -> Float {
        const MASS_TO_VOLUME_FACTOR: Float = 3. / (4. * std::f32::consts::PI * ROCK_DENSITY);
        (self.mass * MASS_TO_VOLUME_FACTOR).cbrt()
    }
}

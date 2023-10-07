use crate::{
    body::Body,
    initial_parameters::{Float, InitialParameters},
};

pub(crate) struct StellarSystem {
    central_body: Body,
    bodies: Vec<Body>,
}

impl StellarSystem {
    pub(crate) fn new(params: InitialParameters) -> StellarSystem {
        let central_body = Body::new(0., 0., params.total_mass * params.stellar_mass_fraction);
        let bodies = (0..params.body_count)
            .map(|_| {
                let mass = params.total_mass * (1. - params.stellar_mass_fraction)
                    / params.body_count as Float;
                Body::new(params.position_variance, params.velocity_variance, mass)
            })
            .collect();
        StellarSystem {
            central_body,
            bodies,
        }
    }
}

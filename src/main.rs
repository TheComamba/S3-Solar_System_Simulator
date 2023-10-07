use initial_parameters::InitialParameters;
use system::StellarSystem;

mod body;
mod initial_parameters;
mod system;

fn main() {
    let params = InitialParameters::default();
    let system = StellarSystem::new(params);
}

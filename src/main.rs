use gui::Gui;
use iced::Sandbox;

mod gui;
mod sim;

fn main() -> iced::Result {
    let mut window_settings = iced::window::Settings::default();
    window_settings.size = (1820, 980);
    let mut settings = iced::settings::Settings::default();
    settings.window = window_settings;
    settings.antialiasing = true;
    Gui::run(settings)
}

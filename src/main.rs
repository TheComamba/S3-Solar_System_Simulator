use gui::Gui;
use iced::Size;

mod gui;
mod sim;

fn main() -> iced::Result {
    let mut window_settings = iced::window::Settings::default();
    window_settings.size = Size::new(1820., 980.);
    iced::application("S3 - Solar System Simulator", Gui::update, Gui::view)
        .antialiasing(true)
        .window(window_settings)
        .run()
}

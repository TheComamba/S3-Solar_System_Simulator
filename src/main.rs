use gui::Gui;
use iced::{Sandbox, Settings};

mod gui;
mod sim;

fn main() -> iced::Result {
    Gui::run(Settings::default())
}

use std::vec;

use iced::{
    widget::{
        canvas::{self, Path},
        Button, Column, Row, Text,
    },
    Color, Length, Size,
};

use crate::sim::{
    body::Body,
    initial_parameters::InitialParameters,
    system::StellarSystem,
    units::{Float, DISTANCE_TO_M, TIME_TO_SECONDS},
};

pub(crate) struct Gui {
    canvas_state: CanvasState,
    time_step: Float,
}

impl Default for Gui {
    fn default() -> Self {
        Gui {
            canvas_state: CanvasState::new(),
            time_step: 1e-2,
        }
    }
}

impl Gui {
    pub(crate) fn update(&mut self, message: GuiMessage) {
        match message {
            GuiMessage::Evolve => {
                self.canvas_state.system.evolve_for(self.time_step);
                // println!("\n{:?}\n", self.canvas_state.system);
            }
            GuiMessage::ChangeTimeStep(step) => self.time_step = step,
            GuiMessage::ChangeBodySize(fac) => self.canvas_state.body_size_factor = fac,
            GuiMessage::ChangeSystemSize(fac) => self.canvas_state.system_size_factor = fac,
            GuiMessage::ChangeSystemCenterX(x) => self.canvas_state.system_offset.0 = x,
            GuiMessage::ChangeSystemCenterY(y) => self.canvas_state.system_offset.1 = y,
        };
        self.canvas_state.bodies_cache.clear();
    }

    pub(crate) fn view(&self) -> iced::Element<'_, GuiMessage> {
        let header = Row::new()
            .push(self.status_block())
            .push(self.control_block());
        let canvas = iced::widget::canvas(&self.canvas_state)
            .width(Length::Fill)
            .height(Length::Fill);
        Column::new()
            .push(header)
            .push(canvas)
            .align_x(iced::Alignment::Center)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }
}

impl Gui {
    fn time_format(time: Float) -> String {
        let time_in_seconds = time * TIME_TO_SECONDS;
        if time_in_seconds.abs() < 1e2 {
            format!("{:.2} s", time_in_seconds)
        } else if time_in_seconds.abs() < 1e4 {
            format!("{:.2} min", time_in_seconds / 60.)
        } else if time_in_seconds.abs() < 1e5 {
            format!("{:.2} h", time_in_seconds / 3600.)
        } else if time_in_seconds.abs() < 1e8 {
            format!("{:.2} d", time_in_seconds / 86400.)
        } else if time_in_seconds.abs() < 1e11 {
            format!("{:.2} y", time_in_seconds / 31536000.)
        } else {
            format!("{:.2} ky", time_in_seconds / 31536000. / 1e3)
        }
    }

    fn distance_format(distance: Float) -> String {
        let distance_in_m = distance * DISTANCE_TO_M;
        if distance_in_m.abs() < 1e3 {
            format!("{:.2} m", distance_in_m)
        } else if distance_in_m.abs() < 1e6 {
            format!("{:.2} km", distance_in_m / 1e3)
        } else if distance_in_m.abs() < 1e9 {
            format!("{:.2} Mm", distance_in_m / 1e6)
        } else {
            format!("{:.2} AU", distance)
        }
    }

    fn status_block(&self) -> iced::Element<'_, GuiMessage> {
        let time_text = Self::time_format(self.canvas_state.system.current_time);
        let time_text = Text::new(format!("Evolution time: {}", time_text));
        Column::new().push(time_text).into()
    }

    fn value_control_block(
        &self,
        value_name: &'static str,
        value: String,
        decrease: GuiMessage,
        increase: GuiMessage,
    ) -> iced::Element<'_, GuiMessage> {
        let decrease_time_step_button = Button::new(Text::new("<<")).on_press(decrease);
        let increase_time_step_button = Button::new(Text::new(">>")).on_press(increase);
        Row::new()
            .push(Text::new(value_name))
            .push(decrease_time_step_button)
            .push(Text::new(value))
            .push(increase_time_step_button)
            .into()
    }

    fn control_block(&self) -> iced::Element<'_, GuiMessage> {
        let time_step_block = self.value_control_block(
            "Time step: ",
            Self::time_format(self.time_step),
            GuiMessage::ChangeTimeStep(0.5 * self.time_step),
            GuiMessage::ChangeTimeStep(2. * self.time_step),
        );
        let system_center_x = self.value_control_block(
            "System center x-pos: ",
            Self::distance_format(self.canvas_state.system_offset.0),
            GuiMessage::ChangeSystemCenterX(
                self.canvas_state.system_offset.0 - 10. * self.canvas_state.system_size_factor,
            ),
            GuiMessage::ChangeSystemCenterX(
                self.canvas_state.system_offset.0 + 10. * self.canvas_state.system_size_factor,
            ),
        );
        let system_center_y = self.value_control_block(
            "System center y-pos: ",
            Self::distance_format(self.canvas_state.system_offset.1),
            GuiMessage::ChangeSystemCenterY(
                self.canvas_state.system_offset.1 - 10. * self.canvas_state.system_size_factor,
            ),
            GuiMessage::ChangeSystemCenterY(
                self.canvas_state.system_offset.1 + 10. * self.canvas_state.system_size_factor,
            ),
        );
        let system_size_block = self.value_control_block(
            "System size factor: ",
            format!("{:.2}", self.canvas_state.system_size_factor),
            GuiMessage::ChangeSystemSize(0.5 * self.canvas_state.system_size_factor),
            GuiMessage::ChangeSystemSize(2. * self.canvas_state.system_size_factor),
        );
        let body_size_block = self.value_control_block(
            "Body size factor: ",
            format!("{:.2}", self.canvas_state.body_size_factor),
            GuiMessage::ChangeBodySize(0.5 * self.canvas_state.body_size_factor),
            GuiMessage::ChangeBodySize(2. * self.canvas_state.body_size_factor),
        );
        let evolve_button = Button::new(Text::new("Evolve")).on_press(GuiMessage::Evolve);
        Column::new()
            .push(time_step_block)
            .push(system_center_x)
            .push(system_center_y)
            .push(system_size_block)
            .push(body_size_block)
            .push(evolve_button)
            .into()
    }
}

struct CanvasState {
    background_cache: canvas::Cache,
    bodies_cache: canvas::Cache,
    system: StellarSystem,
    body_size_factor: Float,
    system_size_factor: Float,
    system_offset: (Float, Float),
}

impl CanvasState {
    fn new() -> CanvasState {
        let params = InitialParameters::default();
        let system = StellarSystem::new(params);
        CanvasState {
            background_cache: canvas::Cache::default(),
            bodies_cache: canvas::Cache::default(),
            system,
            body_size_factor: 1e5,
            system_size_factor: 1e-2,
            system_offset: (0., 0.),
        }
    }

    fn body_display_radius(&self, body: &Body) -> Float {
        let radius = body.radius() * self.body_size_factor;
        if radius < 1.5 {
            1.5
        } else {
            radius.sqrt()
        }
    }

    fn body_display_position(&self, body: &Body) -> (Float, Float) {
        (
            (body.position[0] - self.system_offset.0) / self.system_size_factor,
            (body.position[1] - self.system_offset.1) / self.system_size_factor,
        )
    }
}

impl<GuiMessage> canvas::Program<GuiMessage> for CanvasState {
    type State = ();

    fn draw(
        &self,
        _state: &Self::State,
        renderer: &iced::Renderer,
        _theme: &iced::theme::Theme,
        bounds: iced::Rectangle,
        _cursor: iced::mouse::Cursor,
    ) -> Vec<canvas::Geometry> {
        let background = self
            .background_cache
            .draw(renderer, bounds.size(), |frame| {
                frame.fill_rectangle(frame.center(), Size::UNIT, Color::BLACK)
            });
        let bodies = self.bodies_cache.draw(renderer, bounds.size(), |frame| {
            let bodies = Path::new(|path_builder| {
                for body in &self.system.bodies {
                    let radius = self.body_display_radius(body);
                    let (x, y) = self.body_display_position(body);
                    let pos = frame.center() + iced::Vector::new(x, y);
                    path_builder.circle(pos, radius);
                }
            });
            frame.fill(&bodies, Color::BLACK);
        });
        vec![bodies]
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum GuiMessage {
    Evolve,
    ChangeTimeStep(Float),
    ChangeBodySize(Float),
    ChangeSystemSize(Float),
    ChangeSystemCenterX(Float),
    ChangeSystemCenterY(Float),
}

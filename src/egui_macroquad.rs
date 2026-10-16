use egui_miniquad::EguiMq;
use macroquad::prelude::*;
use miniquad::{self as mq};

pub use egui;
pub use macroquad;

pub struct Egui {
    gl_ctx: InternalGlContext<'static>,
    egui_mq: EguiMq,
    input_subscriber_id: usize,
}

impl Egui {
    pub fn new() -> Self {
        let gl_ctx = unsafe { get_internal_gl() };
        let egui_mq = EguiMq::new(&mut *gl_ctx.quad_context);
        Self {
            gl_ctx,
            egui_mq,
            input_subscriber_id: macroquad::input::utils::register_input_subscriber(),
        }
    }

    fn ui<F>(&mut self, f: F)
    where
        F: FnOnce(&mut dyn mq::RenderingBackend, &egui::Context),
    {
        macroquad::input::utils::repeat_all_miniquad_input(self, self.input_subscriber_id);

        self.egui_mq.run(self.gl_ctx.quad_context, f);
    }

    pub fn draw(&mut self) {
        // Ensure that macroquad's shapes are not goint to be lost, and draw them now
        self.gl_ctx.flush();
        self.egui_mq.draw(self.gl_ctx.quad_context);
    }

    /// Calculates egui ui. Must be called once per frame.
    pub fn calc_ui<F: FnOnce(&egui::Context)>(&mut self, f: F) {
        self.ui(|_, ctx| f(ctx))
    }

    /// Configure egui without beginning or ending a frame.
    pub fn cfg<F: FnOnce(&egui::Context)>(&mut self, f: F) {
        f(self.egui_mq.egui_ctx());
    }
}

impl Default for Egui {
    fn default() -> Self {
        Self::new()
    }
}

impl mq::EventHandler for Egui {
    fn update(&mut self) {}

    fn draw(&mut self) {}

    fn mouse_motion_event(&mut self, x: f32, y: f32) {
        self.egui_mq.mouse_motion_event(x, y);
    }

    fn mouse_wheel_event(&mut self, dx: f32, dy: f32) {
        self.egui_mq.mouse_wheel_event(dx, dy);
    }

    fn mouse_button_down_event(&mut self, mb: mq::MouseButton, x: f32, y: f32) {
        self.egui_mq.mouse_button_down_event(mb, x, y);
    }

    fn mouse_button_up_event(&mut self, mb: mq::MouseButton, x: f32, y: f32) {
        self.egui_mq.mouse_button_up_event(mb, x, y);
    }

    fn char_event(&mut self, character: char, _keymods: mq::KeyMods, _repeat: bool) {
        self.egui_mq.char_event(character);
    }

    fn key_down_event(&mut self, keycode: mq::KeyCode, keymods: mq::KeyMods, _repeat: bool) {
        self.egui_mq.key_down_event(keycode, keymods);
    }

    fn key_up_event(&mut self, keycode: mq::KeyCode, keymods: mq::KeyMods) {
        self.egui_mq.key_up_event(keycode, keymods);
    }
}

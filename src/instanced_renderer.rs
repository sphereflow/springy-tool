//! InstancedRenderer — draws circles and lines with two instanced draw calls
//! instead of ~N immediate-mode macroquad calls. Built for miniquad 0.4.x.

use miniquad as mq;

// macroquad re-exports miniquad; `unsafe { get_internal_gl() }` hands you
// the same RenderingBackend macroquad itself draws through.
use macroquad::{
    color::Color,
    window::{get_internal_gl, screen_height, screen_width},
};

use crate::V2;

/// One spring: endpoints + color. 32 bytes.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct LineInstance {
    pub start: [f32; 2],
    pub end: [f32; 2],
    pub color: [f32; 4],
}

/// One node: center + radius + color. 28 bytes.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct CircleInstance {
    pub center: [f32; 2],
    pub radius: f32,
    pub color: [f32; 4],
}

/// Field order must match the shader-meta uniform lists: miniquad computes
/// byte offsets sequentially (scale@0, offset@8, px@16, width@20).
#[repr(C)]
#[derive(Clone, Copy)]
struct Uniforms {
    scale: [f32; 2], // world -> NDC
    offset: [f32; 2],
    px: f32,    // world units per screen pixel (for AA)
    width: f32, // line width, world units
}

const LINE_VERT: &str = r#"#version 100
attribute vec2 corner;      // unit quad, (-1..1)
attribute vec2 a_start;     // per-instance
attribute vec2 a_end;       // per-instance
attribute vec4 a_color;     // per-instance
uniform vec2 u_scale;
uniform vec2 u_offset;
uniform float u_px;
uniform float u_width;
varying mediump vec4 v_color;
varying mediump float v_side; // signed distance from the line axis, in pixels
void main() {
    vec2 d = a_end - a_start; // direction vector
    float len = length(d);
    d = len > 0.0 ? d / len : vec2(1.0, 0.0); // normalize direction vector
    vec2 perp = vec2(-d.y, d.x); // negative normalized direction vector
    vec2 world = (a_start + a_end) * 0.5
               + corner.x * d * (len * 0.5)
               + corner.y * perp * (u_width * 0.5);
    gl_Position = vec4(world * u_scale + u_offset, 0.0, 1.0);
    v_color = a_color;
    v_side = corner.y * (u_width * 0.5) / u_px;
}"#;

const LINE_FRAG: &str = r#"#version 100
precision mediump float;
varying mediump vec4 v_color;
varying mediump float v_side;
uniform highp float u_width;   // line width, world units
uniform highp float u_px;      // world units per screen pixel
void main() {
    // feather the outer ~1 px on each side, anchored at the actual edge so
    // wider lines stay wider (feather shrinks to fit for thin lines)
    float half_px = u_width * 0.5 / u_px;
    if (half_px <= 0.0) { discard; }
    float aa = min(1.0, half_px);
    float alpha = (1.0 - smoothstep(half_px - aa, half_px, abs(v_side))) * v_color.a;
    gl_FragColor = vec4(v_color.rgb, alpha);
}"#;

const CIRCLE_VERT: &str = r#"#version 100
attribute vec2 corner;      // per-vertex unit quad
attribute vec2 a_center;    // per-instance
attribute float a_radius;   // per-instance
attribute vec4 a_color;     // per-instance
uniform vec2 u_scale;
uniform vec2 u_offset;
uniform float u_px;
varying mediump vec2 v_local;   // offset from center, world units
varying mediump float v_radius;
varying mediump vec4 v_color;
void main() {
    v_local = corner * a_radius;
    vec2 world = a_center + v_local;
    gl_Position = vec4(world * u_scale + u_offset, 0.0, 1.0);
    v_radius = a_radius;
    v_color = a_color;
}"#;

const CIRCLE_FRAG: &str = r#"#version 100
precision mediump float;
varying mediump vec2 v_local;
varying mediump float v_radius;
varying mediump vec4 v_color;
void main() {
    float d = length(v_local);
    // feather the rim over ~1.5 px
    float alpha = (1.0 - smoothstep(v_radius - 1.5, v_radius, d)) * v_color.a;
    gl_FragColor = vec4(v_color.rgb, alpha);
}"#;

/// quad vertices, plus 2 triangles of indices — shared by both pipelines
const QUAD: [[f32; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
const QUAD_INDICES: [u16; 6] = [0, 1, 2, 0, 2, 3];

pub struct InstancedRenderer {
    line_pipeline: mq::Pipeline,
    line_quad_vb: mq::BufferId,
    line_ib: mq::BufferId,
    line_inst_vb: mq::BufferId,
    line_instances: Vec<LineInstance>,

    circle_pipeline: mq::Pipeline,
    circle_quad_vb: mq::BufferId,
    circle_ib: mq::BufferId,
    circle_inst_vb: mq::BufferId,
    circle_instances: Vec<CircleInstance>,

    view: Uniforms,
    pub line_width: f32, // world units
}

fn line_meta() -> mq::ShaderMeta {
    mq::ShaderMeta {
        images: vec![],
        uniforms: mq::UniformBlockLayout {
            uniforms: vec![
                mq::UniformDesc::new("u_scale", mq::UniformType::Float2),
                mq::UniformDesc::new("u_offset", mq::UniformType::Float2),
                mq::UniformDesc::new("u_px", mq::UniformType::Float1),
                mq::UniformDesc::new("u_width", mq::UniformType::Float1),
            ],
        },
    }
}

fn circle_meta() -> mq::ShaderMeta {
    mq::ShaderMeta {
        images: vec![],
        uniforms: mq::UniformBlockLayout {
            uniforms: vec![
                mq::UniformDesc::new("u_scale", mq::UniformType::Float2),
                mq::UniformDesc::new("u_offset", mq::UniformType::Float2),
                mq::UniformDesc::new("u_px", mq::UniformType::Float1),
            ],
        },
    }
}

fn blending() -> Option<mq::BlendState> {
    Some(mq::BlendState::new(
        mq::Equation::Add,
        mq::BlendFactor::Value(mq::BlendValue::SourceAlpha),
        mq::BlendFactor::OneMinusValue(mq::BlendValue::SourceAlpha),
    ))
}

impl InstancedRenderer {
    /// Call once per app run, after macroquad is initialized (start of the
    /// main loop is fine). Requires instancing support (desktop GL3.2+ ok).
    pub fn new() -> InstancedRenderer {
        unsafe {
            let macroquad::window::InternalGlContext {
                quad_context: ctx, ..
            } = get_internal_gl();

            assert!(ctx.info().features.instancing, "backend lacks instancing");

            let line_quad_vb = ctx.new_buffer(
                mq::BufferType::VertexBuffer,
                mq::BufferUsage::Immutable,
                mq::BufferSource::slice(&QUAD),
            );
            let line_ib = ctx.new_buffer(
                mq::BufferType::IndexBuffer,
                mq::BufferUsage::Immutable,
                mq::BufferSource::slice(&QUAD_INDICES),
            );
            let line_inst_vb = ctx.new_buffer(
                mq::BufferType::VertexBuffer,
                mq::BufferUsage::Dynamic,
                mq::BufferSource::empty::<LineInstance>(1024),
            );
            let line_shader = ctx
                .new_shader(
                    mq::ShaderSource::Glsl {
                        vertex: LINE_VERT,
                        fragment: LINE_FRAG,
                    },
                    line_meta(),
                )
                .expect("line shader");
            let line_pipeline = ctx.new_pipeline(
                &[
                    mq::BufferLayout::default(), // per-vertex quad
                    mq::BufferLayout {
                        stride: 0,
                        step_func: mq::VertexStep::PerInstance,
                        step_rate: 1,
                    },
                ],
                &[
                    mq::VertexAttribute::with_buffer("corner", mq::VertexFormat::Float2, 0),
                    mq::VertexAttribute::with_buffer("a_start", mq::VertexFormat::Float2, 1),
                    mq::VertexAttribute::with_buffer("a_end", mq::VertexFormat::Float2, 1),
                    mq::VertexAttribute::with_buffer("a_color", mq::VertexFormat::Float4, 1),
                ],
                line_shader,
                mq::PipelineParams {
                    color_blend: blending(),
                    depth_write: false,
                    depth_test: mq::Comparison::Always,
                    ..Default::default()
                },
            );

            let circle_quad_vb = ctx.new_buffer(
                mq::BufferType::VertexBuffer,
                mq::BufferUsage::Immutable,
                mq::BufferSource::slice(&QUAD),
            );
            let circle_ib = ctx.new_buffer(
                mq::BufferType::IndexBuffer,
                mq::BufferUsage::Immutable,
                mq::BufferSource::slice(&QUAD_INDICES),
            );
            let circle_inst_vb = ctx.new_buffer(
                mq::BufferType::VertexBuffer,
                mq::BufferUsage::Dynamic,
                mq::BufferSource::empty::<CircleInstance>(1024),
            );
            let circle_shader = ctx
                .new_shader(
                    mq::ShaderSource::Glsl {
                        vertex: CIRCLE_VERT,
                        fragment: CIRCLE_FRAG,
                    },
                    circle_meta(),
                )
                .expect("circle shader");
            let circle_pipeline = ctx.new_pipeline(
                &[
                    mq::BufferLayout::default(),
                    mq::BufferLayout {
                        stride: 0,
                        step_func: mq::VertexStep::PerInstance,
                        step_rate: 1,
                    },
                ],
                &[
                    mq::VertexAttribute::with_buffer("corner", mq::VertexFormat::Float2, 0),
                    mq::VertexAttribute::with_buffer("a_center", mq::VertexFormat::Float2, 1),
                    mq::VertexAttribute::with_buffer("a_radius", mq::VertexFormat::Float1, 1),
                    mq::VertexAttribute::with_buffer("a_color", mq::VertexFormat::Float4, 1),
                ],
                circle_shader,
                mq::PipelineParams {
                    color_blend: blending(),
                    depth_write: false,
                    depth_test: mq::Comparison::Always,
                    ..Default::default()
                },
            );

            InstancedRenderer {
                line_pipeline,
                line_quad_vb,
                line_ib,
                line_inst_vb,
                line_instances: Vec::new(),
                circle_pipeline,
                circle_quad_vb,
                circle_ib,
                circle_inst_vb,
                circle_instances: Vec::new(),
                view: Uniforms {
                    scale: [1.0, 1.0],
                    offset: [0.0, 0.0],
                    px: 1.0,
                    width: 1.0,
                },
                line_width: 1.0,
            }
        }
    }

    /// Configure the world->screen transform: world point `center` maps to the
    /// screen center, `zoom` = screen pixels per world unit, y grows downward
    /// (macroquad convention). Call each frame before pushing instances.
    pub fn set_view(&mut self, center: V2, zoom: f32) {
        let (w, h) = (screen_width(), screen_height());
        self.view = Uniforms {
            scale: [2.0 * zoom / w, -2.0 * zoom / h],
            offset: [
                -2.0 * zoom * center.x as f32 / w,
                2.0 * zoom * center.y as f32 / h,
            ],
            px: 1.0 / zoom,
            width: self.line_width,
        };
    }

    pub fn push_line(&mut self, start: V2, end: V2, color: Color) {
        self.line_instances.push(LineInstance {
            start: [start.x as f32, start.y as f32],
            end: [end.x as f32, end.y as f32],
            color: [color.r, color.g, color.b, color.a],
        });
    }

    pub fn push_circle(&mut self, center: V2, radius: f32, color: Color) {
        self.circle_instances.push(CircleInstance {
            center: [center[0] as f32, center[1] as f32],
            radius,
            color: [color.r, color.g, color.b, color.a],
        });
    }

    /// Render everything pushed since the last call, then clear the queues.
    /// Call after `clear_background` (and after any macroquad draws that should
    /// be underneath), before UI that should sit on top.
    pub fn draw(&mut self) {
        let line_n = self.line_instances.len();
        let circle_n = self.circle_instances.len();
        if line_n == 0 && circle_n == 0 {
            return;
        }
        let mut gl = unsafe { macroquad::window::get_internal_gl() };

        // submit macroquad's queued draws (each in its own pass) so the
        // instanced geometry below lands after them, then open our pass
        gl.flush();
        let ctx = gl.quad_context;

        ctx.begin_default_pass(mq::PassAction::Nothing);

        if line_n > 0 {
            self.ensure_line_capacity(&mut *ctx, line_n);
            ctx.buffer_update(
                self.line_inst_vb,
                mq::BufferSource::slice(&self.line_instances),
            );
            ctx.apply_pipeline(&self.line_pipeline);
            ctx.apply_bindings(&mq::Bindings {
                vertex_buffers: vec![self.line_quad_vb, self.line_inst_vb],
                index_buffer: self.line_ib,
                images: vec![],
            });
            ctx.apply_uniforms(mq::UniformsSource::table(&self.view));
            ctx.draw(0, 6, line_n as i32);
        }

        if circle_n > 0 {
            self.ensure_circle_capacity(&mut *ctx, circle_n);
            ctx.buffer_update(
                self.circle_inst_vb,
                mq::BufferSource::slice(&self.circle_instances),
            );
            ctx.apply_pipeline(&self.circle_pipeline);
            ctx.apply_bindings(&mq::Bindings {
                vertex_buffers: vec![self.circle_quad_vb, self.circle_inst_vb],
                index_buffer: self.circle_ib,
                images: vec![],
            });
            ctx.apply_uniforms(mq::UniformsSource::table(&self.view));
            ctx.draw(0, 6, circle_n as i32);
        }

        ctx.end_render_pass();
        self.line_instances.clear();
        self.circle_instances.clear();
    }

    fn ensure_line_capacity(&mut self, ctx: &mut dyn mq::RenderingBackend, n: usize) {
        let _ = ctx;
        let cap = ctx.buffer_size(self.line_inst_vb) / std::mem::size_of::<LineInstance>();
        if cap < n {
            ctx.delete_buffer(self.line_inst_vb);
            self.line_inst_vb = ctx.new_buffer(
                mq::BufferType::VertexBuffer,
                mq::BufferUsage::Dynamic,
                mq::BufferSource::empty::<LineInstance>(n.next_power_of_two()),
            );
        }
    }

    fn ensure_circle_capacity(&mut self, ctx: &mut dyn mq::RenderingBackend, n: usize) {
        let cap = ctx.buffer_size(self.circle_inst_vb) / std::mem::size_of::<CircleInstance>();
        if cap < n {
            ctx.delete_buffer(self.circle_inst_vb);
            self.circle_inst_vb = ctx.new_buffer(
                mq::BufferType::VertexBuffer,
                mq::BufferUsage::Dynamic,
                mq::BufferSource::empty::<CircleInstance>(n.next_power_of_two()),
            );
        }
    }
}

impl Default for InstancedRenderer {
    fn default() -> Self {
        Self::new()
    }
}

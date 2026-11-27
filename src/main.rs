//! The next.BREP adapter: answers the B-rep protocol's requests on stdin and stdout, so the
//! conformance suite can drive this kernel.

use std::process::ExitCode;

use brep_kernel::{
    AffineTransform, BooleanOperation, BooleanOptions, BrepSolid, KernelRefusal, LawSegment,
    NurbsCurve, OffsetCorners, PointClass, RadiusLaw, RefusalClass, SectionPlacement, SweepPath,
    Vec3, Vec4, boolean_operation, boolean_operation_nary, chamfer_edges_asymmetric,
    classify_point, delete_faces_and_heal, draft_faces, export_step, extrude_profile_brep,
    extrude_profile_with_holes, fillet_edges_conic_setback, fillet_edges_g2, fillet_edges_radii,
    fillet_edges_setback, fillet_edges_setback_g2, fillet_edges_variable_law, fillet_faces,
    full_round, import_step, loft_profile_brep, make_arc, make_box_brep, make_cone_brep,
    make_cylinder_brep, make_line, make_periodic_bspline, make_slab_brep, make_sphere_brep,
    make_torus_brep, offset_shell_rounded, offset_solid, project_point_to_curve,
    resolve_face_by_point, revolve_profile_brep, revolve_profile_with_holes, rotate_faces,
    screw_profile_brep, sew_step, solid_closest_points, solid_display, solid_interference,
    solid_least_box, solid_mass_properties, solid_mass_properties_full, solid_outward_normals,
    solid_principal_curvatures, solid_section, solid_section_area, split_solid_by_plane,
    sweep_profile_along_chain, tessellate_brep_watertight, thicken_faces, transform_brep,
};
use brep_protocol::{Refusal, model, model::PathSeg, model::Primitive, profile};

/// The contract's `10·ε`, within which a point names an edge or a face.
const ADDRESS: f64 = 1e-6;

/// The tolerance a point is classified to, mm: the suite's, as the contract asks of every setting.
const CLASSIFY: f64 = 1e-7;

/// The loosest tolerance governing next.BREP's results, mm: the edge-to-face disagreement its
/// validator accepts in a valid solid (`KernelTolerances::pcurve_consistency`, at least 4e-3).
const DECLARED: f64 = 4e-3;

struct NextBrep;

/// A solid: next.BREP's, or `None` for the empty solid, which next.BREP has no value for.
type Solid = Option<BrepSolid>;

/// next.BREP's refusal as the protocol's, by the class next.BREP gives it: geometry it does not
/// model yet is unsupported, an input it rejects or that has no single answer is invalid, and
/// anything else is a failure.
fn refused(e: KernelRefusal) -> Refusal {
    let reason = e.to_string();
    match e.class {
        RefusalClass::UnsupportedGeometry { .. } => Refusal::unsupported(reason),
        RefusalClass::InvalidInput { .. } | RefusalClass::IllPosed { .. } => {
            Refusal::invalid(reason)
        }
        _ => Refusal::failed(reason),
    }
}

/// next.BREP's untyped error, which says nothing of its kind.
fn failed(e: String) -> Refusal {
    Refusal::failed(e)
}

fn v3(a: model::Vec3) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}

/// A result as a solid: next.BREP says a result is empty by a solid with no faces.
fn solid(s: BrepSolid) -> Solid {
    if s.shells.iter().all(|shell| shell.faces.is_empty()) {
        None
    } else {
        Some(s)
    }
}

/// The map `p ↦ M·p + t` for the rows `[M | t]`, as next.BREP's row-major 4 × 4.
fn affine(m: [[f64; 4]; 3]) -> Result<AffineTransform, Refusal> {
    let mut e = [0.0; 16];
    for i in 0..3 {
        e[4 * i..4 * i + 4].copy_from_slice(&m[i]);
    }
    e[15] = 1.0;
    AffineTransform::new(e).map_err(Refusal::invalid)
}

fn det3(m: [[f64; 4]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// next.BREP's affine map of a whole solid; a map turning space over turns the faces round.
fn mapped(s: &BrepSolid, m: [[f64; 4]; 3]) -> Result<BrepSolid, Refusal> {
    transform_brep(s, affine(m)?, det3(m) < 0.0).map_err(failed)
}

/// The solid built about the world origin, moved onto `frame`.
fn place(s: BrepSolid, f: &model::Frame) -> Result<BrepSolid, Refusal> {
    mapped(&s, [0, 1, 2].map(|i| [f.x[i], f.y[i], f.z[i], f.origin[i]]))
}

/// Where a profile is drawn: its points in space, and the plane's axes.
#[derive(Clone, Copy)]
enum Plane {
    /// The xy plane at height `h`.
    Xy(f64),
    /// A revolution's meridian: x out from the z axis, y along it.
    Meridian,
}

impl Plane {
    fn point(self, p: [f64; 2]) -> Vec3 {
        match self {
            Plane::Xy(h) => Vec3::new(p[0], p[1], h),
            Plane::Meridian => Vec3::new(p[0], 0.0, p[1]),
        }
    }

    fn axes(self) -> (Vec3, Vec3) {
        match self {
            Plane::Xy(_) => (Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)),
            Plane::Meridian => (Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        }
    }
}

/// The suite's profile as next.BREP's curves, the loop's in order. next.BREP's profile
/// constructions take one loop, so a profile with holes is one they have no call for. A clockwise
/// arc is the same arc drawn counter-clockwise about the plane's y axis turned back.
fn sketch(p: &profile::Profile, plane: Plane) -> Result<Vec<NurbsCurve>, Refusal> {
    if p.loops.len() > 1 {
        return Err(Refusal::unsupported(
            "next.BREP's profile constructions take one loop: a profile with holes",
        ));
    }
    p.loops
        .first()
        .map_or(Ok(Vec::new()), |segs| sketch_loop(segs, plane, true))
}

/// A profile with holes as next.BREP's curves: its first loop, the outer one, and its other
/// loops, the holes, each a closed chain.
fn sketch_holes(
    p: &profile::Profile,
    plane: Plane,
) -> Result<(Vec<NurbsCurve>, Vec<Vec<NurbsCurve>>), Refusal> {
    let mut loops = p.loops.iter().map(|segs| sketch_loop(segs, plane, false));
    let outer = loops.next().unwrap_or(Ok(Vec::new()))?;
    Ok((outer, loops.collect::<Result<_, _>>()?))
}

/// One loop of the suite's profile as next.BREP's curves, in order.
fn sketch_loop(
    segs: &[profile::Seg],
    plane: Plane,
    halve_circles: bool,
) -> Result<Vec<NurbsCurve>, Refusal> {
    let (ax, ay) = plane.axes();
    let mut curves = Vec::new();
    for seg in segs {
        let next = match *seg {
            profile::Seg::Line { a, b } => {
                make_line(plane.point(a), plane.point(b)).map_err(failed)?
            }
            profile::Seg::Arc {
                center,
                radius,
                start,
                sweep,
            } => {
                let (y, a0, a1) = if sweep >= 0.0 {
                    (ay, start, start + sweep)
                } else {
                    (ay.scale(-1.0), -start, -start - sweep)
                };
                let c = plane.point(center);
                // A whole circle goes as its two halves: next.BREP's lofts take a section of two
                // curves or more.
                if halve_circles && a1 - a0 >= std::f64::consts::TAU - 1e-12 {
                    let mid = a0 + (a1 - a0) / 2.0;
                    curves.push(make_arc(c, ax, y, radius, a0, mid).map_err(failed)?);
                    make_arc(c, ax, y, radius, mid, a1).map_err(failed)?
                } else {
                    make_arc(c, ax, y, radius, a0, a1).map_err(failed)?
                }
            }
            // next.BREP's periodic B-spline is on evenly spaced knots, from one period's points;
            // the protocol's lists those points and then the first `degree` again.
            profile::Seg::Spline(ref sp) => {
                let steps: Vec<f64> = sp.knots.windows(2).map(|k| k[1] - k[0]).collect();
                let even = steps
                    .iter()
                    .all(|d| (d - steps[0]).abs() <= 1e-12 * steps[0].abs().max(1.0));
                if !even {
                    return Err(Refusal::unsupported(
                        "next.BREP's periodic B-spline takes evenly spaced knots",
                    ));
                }
                let period = sp.ctrl.len().saturating_sub(sp.degree);
                let points: Vec<Vec3> = sp.ctrl[..period].iter().map(|&q| plane.point(q)).collect();
                make_periodic_bspline(sp.degree, &points).map_err(failed)?
            }
            // The rational quadratic through `a`, `apex` and `b`, weighted 1, `w`, 1.
            profile::Seg::Conic { a, apex, b, w } => NurbsCurve::new(
                2,
                vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                vec![
                    Vec4::from_point(plane.point(a), 1.0),
                    Vec4::from_point(plane.point(apex), w),
                    Vec4::from_point(plane.point(b), 1.0),
                ],
            )
            .map_err(failed)?,
            profile::Seg::Gn { .. } => {
                return Err(Refusal::unsupported(
                    "next.BREP profiles have no blend-section segment",
                ));
            }
            _ => {
                return Err(Refusal::unsupported(
                    "the adapter does not know this kind of profile segment",
                ));
            }
        };
        curves.push(next);
    }
    Ok(curves)
}

/// The suite's sweep path as next.BREP's chain of curves, from the origin starting along z: each
/// arc bends toward the current x about an axis along the current y.
fn trace(segs: &[PathSeg]) -> Result<Vec<NurbsCurve>, Refusal> {
    let (mut at, mut x, mut z) = (
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    );
    let mut curves = Vec::new();
    for s in segs {
        match *s {
            PathSeg::Line { length } => {
                let to = at.add(z.scale(length));
                curves.push(make_line(at, to).map_err(failed)?);
                at = to;
            }
            PathSeg::Arc { radius, angle } => {
                let c = at.add(x.scale(radius));
                curves.push(make_arc(c, x.scale(-1.0), z, radius, 0.0, angle).map_err(failed)?);
                let (s, co) = angle.sin_cos();
                at = c.add(x.scale(-radius * co)).add(z.scale(radius * s));
                (x, z) = (x.scale(co).add(z.scale(-s)), z.scale(co).add(x.scale(s)));
            }
        }
    }
    Ok(curves)
}

fn make(p: &Primitive) -> Result<BrepSolid, Refusal> {
    let o = Vec3::new(0.0, 0.0, 0.0);
    let z = Vec3::new(0.0, 0.0, 1.0);
    let (s, frame) = match p {
        Primitive::Box { frame, size } => (
            make_box_brep(o, size[0], size[1], size[2]).map_err(failed)?,
            frame,
        ),
        Primitive::Cylinder {
            frame,
            radius,
            height,
        } => (
            make_cylinder_brep(o, z, *radius, *height).map_err(failed)?,
            frame,
        ),
        Primitive::Sphere { frame, radius } => {
            (make_sphere_brep(o, *radius, z).map_err(failed)?, frame)
        }
        Primitive::Torus {
            frame,
            major,
            minor,
        } => (
            make_torus_brep(o, z, *major, *minor).map_err(failed)?,
            frame,
        ),
        Primitive::Cone {
            frame,
            bottom,
            top,
            height,
        } => (
            make_cone_brep(o, z, *bottom, *top, *height).map_err(failed)?,
            frame,
        ),
        // A profile with holes is next.BREP's extrusion of an outer loop and its holes.
        Primitive::Extrude {
            frame,
            profile,
            height,
        } if profile.loops.len() > 1 => {
            let (outer, holes) = sketch_holes(profile, Plane::Xy(0.0))?;
            (
                extrude_profile_with_holes(&outer, &holes, z, *height).map_err(refused)?,
                frame,
            )
        }
        Primitive::Extrude {
            frame,
            profile,
            height,
        } => (
            extrude_profile_brep(&sketch(profile, Plane::Xy(0.0))?, z, *height).map_err(refused)?,
            frame,
        ),
        Primitive::Revolve {
            frame,
            profile,
            angle,
        } if profile.loops.len() > 1 => {
            let (outer, holes) = sketch_holes(profile, Plane::Meridian)?;
            (
                revolve_profile_with_holes(&outer, &holes, o, z, *angle).map_err(refused)?,
                frame,
            )
        }
        Primitive::Revolve {
            frame,
            profile,
            angle,
        } => (
            revolve_profile_brep(&sketch(profile, Plane::Meridian)?, o, z, *angle)
                .map_err(failed)?,
            frame,
        ),
        // next.BREP's screw: the meridian profile turned right-handed about z through `angle`,
        // rising `rise` per radian.
        Primitive::Screw {
            frame,
            profile,
            rise,
            angle,
        } => (
            screw_profile_brep(&sketch(profile, Plane::Meridian)?, o, z, *rise, *angle)
                .map_err(refused)?,
            frame,
        ),
        // next.BREP's slab: the block under the bicubic height field, on the frame.
        Primitive::Slab {
            frame,
            size,
            knots_u,
            knots_v,
            heights,
        } => (
            make_slab_brep(*size, knots_u, knots_v, heights).map_err(refused)?,
            frame,
        ),
        // The top section is the scaled copy of the profile the protocol sends.
        Primitive::Loft {
            frame,
            profile,
            height,
            top,
            ..
        } => {
            let sections = [
                sketch(profile, Plane::Xy(0.0))?,
                sketch(top, Plane::Xy(*height))?,
            ];
            (loft_profile_brep(&sections).map_err(refused)?, frame)
        }
        // The protocol's sweep carries the profile rigidly: each line extrudes it and each arc
        // revolves it, next.BREP's `Rigid` placement.
        Primitive::Sweep {
            frame,
            profile,
            path,
        } => {
            let chain = SweepPath::new(trace(path)?, Vec::new()).map_err(refused)?;
            let section = sketch(profile, Plane::Xy(0.0))?;
            let s = sweep_profile_along_chain(
                &section,
                &chain,
                0.0,
                None,
                None,
                SectionPlacement::Rigid,
                "",
            )
            .map_err(refused)?;
            (s, frame)
        }
        Primitive::SweepPath {
            frame,
            profile,
            path,
            twist,
        } => {
            let ctrl = path
                .ctrl
                .iter()
                .map(|c| Vec4 {
                    x: c[0],
                    y: c[1],
                    z: c[2],
                    w: 1.0,
                })
                .collect();
            let curve = NurbsCurve::new(path.degree, path.knots.clone(), ctrl).map_err(failed)?;
            // Each loop as the one-loop profile constructions draw it (a whole circle halved).
            let mut loops = profile
                .loops
                .iter()
                .map(|segs| sketch_loop(segs, Plane::Xy(0.0), true));
            let section = loops.next().unwrap_or(Ok(Vec::new()))?;
            let holes = loops.collect::<Result<Vec<_>, _>>()?;
            // The profile as drawn, carried rigidly by the path from its start (the protocol's
            // path starts at the origin along z, square to the profile's plane) by the
            // rotation-minimizing frame with the twist about the path; its stations refined to
            // 1e-7, under the suite's area bar of 1e-6 relative.
            let swept = |section: &[NurbsCurve]| {
                brep_kernel::sweep_profile_twisted_within(
                    section,
                    &curve,
                    *twist,
                    None,
                    brep_kernel::SectionPlacement::Rigid,
                    None,
                    Some(1e-7),
                )
                .map_err(refused)
                .map(|swept| swept.solid)
            };
            // A hole is the same rigid sweep of its loop, taken away: the same frame carries both.
            let mut s = swept(&section)?;
            for hole in &holes {
                s = boolean(BooleanOperation::Subtract, &Some(s), &Some(swept(hole)?))?
                    .ok_or_else(|| Refusal::unsupported("the hole took the whole sweep"))?;
            }
            (s, frame)
        }
    };
    place(s, frame)
}

/// The boolean of `a` and `b`, the empty solid answered as the protocol asks.
fn boolean(op: BooleanOperation, a: &Solid, b: &Solid) -> Result<Solid, Refusal> {
    match (op, a, b) {
        (_, Some(a), Some(b)) => Ok(solid(
            boolean_operation(a, b, op, &BooleanOptions::default()).map_err(refused)?,
        )),
        (BooleanOperation::Union, x, None) | (BooleanOperation::Subtract, x, None) => Ok(x.clone()),
        (BooleanOperation::Union, None, x) => Ok(x.clone()),
        _ => Ok(None),
    }
}

/// The edge a point lies on, by next.BREP's own projection onto each edge, within `10·ε`.
fn edge_at(s: &BrepSolid, p: model::Vec3) -> Result<u64, Refusal> {
    let q = v3(p);
    s.edges
        .iter()
        .filter(|e| !e.degenerate)
        .filter_map(|e| {
            let hit = project_point_to_curve(&e.curve, q).ok()?;
            let within = hit.u >= e.t0.min(e.t1) && hit.u <= e.t0.max(e.t1);
            (within && hit.distance <= ADDRESS).then_some((hit.distance, e.id))
        })
        .min_by(|a, b| a.0.total_cmp(&b.0))
        .map(|(_, id)| id)
        .ok_or_else(|| Refusal::invalid(format!("no edge passes through {p:?}")))
}

/// The face each point lies on. next.BREP's resolver takes the nearest face within a millimetre,
/// so the point is first held to the contract's `10·ε`: next.BREP must classify it on the
/// boundary at that tolerance, and a point on an edge, where several faces meet, is one the
/// adapter does not choose for.
fn faces_at(s: &BrepSolid, points: &[model::Vec3]) -> Result<Vec<u64>, Refusal> {
    points
        .iter()
        .map(|&p| {
            let on = classify_point(v3(p), s, ADDRESS).map_err(failed)?.class;
            if on != PointClass::On {
                return Err(Refusal::invalid(format!("no face passes through {p:?}")));
            }
            if edge_at(s, p).is_ok() {
                return Err(Refusal::unsupported(
                    "a point on an edge or a vertex, where several faces meet",
                ));
            }
            resolve_face_by_point(s, v3(p)).map_err(refused)
        })
        .collect()
}

/// The points that name edges, each held to the contract's `10·ε` first: next.BREP takes an
/// edge by a point and finds the nearest within a millimetre.
fn edge_points(s: &BrepSolid, ps: &[model::Vec3]) -> Result<Vec<Vec3>, Refusal> {
    ps.iter().map(|&p| edge_at(s, p).map(|_| v3(p))).collect()
}

/// The axis and angle of a rotation matrix, and the point of the axis nearest the origin, where
/// the rows `[R | t]` turn about that axis without sliding along it.
fn turn(m: [[f64; 4]; 3]) -> Option<(Vec3, Vec3, f64)> {
    let r = |i: usize, j: usize| m[i][j];
    let angle = ((r(0, 0) + r(1, 1) + r(2, 2) - 1.0) / 2.0)
        .clamp(-1.0, 1.0)
        .acos();
    if angle.abs() < 1e-12 || (std::f64::consts::PI - angle).abs() < 1e-9 {
        return None;
    }
    let n = Vec3::new(r(2, 1) - r(1, 2), r(0, 2) - r(2, 0), r(1, 0) - r(0, 1))
        .normalized()
        .ok()?;
    let t = Vec3::new(m[0][3], m[1][3], m[2][3]);
    if t.dot(n).abs() > 1e-9 * (1.0 + t.length()) {
        return None;
    }
    let p = t
        .add(n.cross(t).scale(1.0 / (angle / 2.0).tan()))
        .scale(0.5);
    Some((p, n, angle))
}

/// The solid an operation acts on, or the refusal for an empty one.
fn some<'a>(s: &'a Solid, what: &str) -> Result<&'a BrepSolid, Refusal> {
    s.as_ref()
        .ok_or_else(|| Refusal::invalid(format!("nothing to {what}")))
}

impl brep_protocol::Kernel for NextBrep {
    type Solid = Solid;

    fn name(&self) -> String {
        "next.BREP".into()
    }

    fn version(&self) -> String {
        "0.32.2".into()
    }

    fn tolerance(&self) -> f64 {
        DECLARED
    }

    fn configuration(&self) -> String {
        format!("classify_point tolerance {CLASSIFY} mm; tolerances 1e-7 mm where a call takes one")
    }

    fn make(&mut self, p: &Primitive) -> Result<Solid, Refusal> {
        make(p).map(Some)
    }

    /// next.BREP's loft through every section, each profile drawn at its height.
    fn loft_through(
        &mut self,
        frame: &model::Frame,
        profiles: &[profile::Profile],
        heights: &[f64],
    ) -> Result<Solid, Refusal> {
        let sections = profiles
            .iter()
            .zip(heights)
            .map(|(p, &h)| sketch(p, Plane::Xy(h)))
            .collect::<Result<Vec<_>, _>>()?;
        let s = loft_profile_brep(&sections).map_err(refused)?;
        place(s, frame).map(Some)
    }

    fn union(&mut self, a: &Solid, b: &Solid) -> Result<Solid, Refusal> {
        boolean(BooleanOperation::Union, a, b)
    }

    fn intersection(&mut self, a: &Solid, b: &Solid) -> Result<Solid, Refusal> {
        boolean(BooleanOperation::Intersect, a, b)
    }

    fn difference(&mut self, a: &Solid, b: &Solid) -> Result<Solid, Refusal> {
        boolean(BooleanOperation::Subtract, a, b)
    }

    fn union_all(&mut self, solids: &[&Solid]) -> Result<Solid, Refusal> {
        let solids: Vec<BrepSolid> = solids.iter().filter_map(|s| (*s).clone()).collect();
        match solids.len() {
            0 => Ok(None),
            1 => Ok(solids.into_iter().next()),
            _ => Ok(solid(
                boolean_operation_nary(&solids, BooleanOperation::Union).map_err(refused)?,
            )),
        }
    }

    fn difference_all(&mut self, target: &Solid, tools: &[&Solid]) -> Result<Solid, Refusal> {
        let Some(target) = target else {
            return Ok(None);
        };
        let mut operands = vec![target.clone()];
        operands.extend(tools.iter().filter_map(|s| (*s).clone()));
        if operands.len() == 1 {
            return Ok(Some(target.clone()));
        }
        Ok(solid(
            boolean_operation_nary(&operands, BooleanOperation::Subtract).map_err(refused)?,
        ))
    }

    fn volume(&mut self, s: &Solid) -> Result<f64, Refusal> {
        let Some(s) = s else { return Ok(0.0) };
        Ok(solid_mass_properties(s).map_err(failed)?.volume)
    }

    fn area(&mut self, s: &Solid) -> Result<f64, Refusal> {
        let Some(s) = s else { return Ok(0.0) };
        Ok(solid_mass_properties(s).map_err(failed)?.surface_area)
    }

    fn centroid(&mut self, s: &Solid) -> Result<Option<model::Vec3>, Refusal> {
        let Some(s) = s else { return Ok(None) };
        let c = solid_mass_properties_full(s).map_err(failed)?.centroid;
        Ok(Some([c.x, c.y, c.z]))
    }

    fn inertia(&mut self, s: &Solid) -> Result<[[f64; 3]; 3], Refusal> {
        let Some(s) = s else {
            return Ok([[0.0; 3]; 3]);
        };
        Ok(solid_mass_properties_full(s).map_err(failed)?.inertia)
    }

    fn check(&mut self, s: &Solid) -> Result<Vec<String>, Refusal> {
        let Some(s) = s else {
            return Ok(Vec::new());
        };
        Ok(s.validate()
            .iter()
            .map(|issue| format!("{issue:?}"))
            .collect())
    }

    fn classify(&mut self, s: &Solid, pts: &[model::Vec3]) -> Result<Vec<model::Side>, Refusal> {
        let Some(s) = s else {
            return Ok(vec![model::Side::Out; pts.len()]);
        };
        // One classifier for all the points: next.BREP's classify_point builds the same one per
        // call.
        // next.BREP's classifier calls a point On within ten times the tolerance it is built with,
        // so it
        // is built with a tenth of CLASSIFY: On means within CLASSIFY of the boundary, as
        // configured.
        let classifier = brep_kernel::SolidClassifier::new(s, CLASSIFY / 10.0).map_err(failed)?;
        pts.iter()
            .map(|&p| {
                Ok(match classifier.classify(v3(p)).map_err(failed)?.class {
                    PointClass::In => model::Side::In,
                    PointClass::On => model::Side::On,
                    PointClass::Out => model::Side::Out,
                })
            })
            .collect()
    }

    fn transform(&mut self, s: &Solid, matrix: [[f64; 4]; 3]) -> Result<Solid, Refusal> {
        let Some(s) = s else { return Ok(None) };
        mapped(s, matrix).map(Some)
    }

    fn affine(&mut self, s: &Solid, matrix: [[f64; 4]; 3]) -> Result<Solid, Refusal> {
        let Some(s) = s else { return Ok(None) };
        mapped(s, matrix).map(Some)
    }

    fn fillet(&mut self, s: &Solid, edges: &[model::Vec3], radius: f64) -> Result<Solid, Refusal> {
        let s = some(s, "fillet")?;
        let found = edge_points(s, edges)?;
        // The suite's fillet follows each picked edge's tangent chain and refuses a pair of
        // balls that could not each touch both faces they roll on.
        let options = brep_kernel::FilletOptions {
            follow_tangent_chains: true,
            refuse_pinch: true,
        };
        Ok(Some(
            brep_kernel::fillet_edges_with_options(s, &found, None, radius, options, None)
                .map_err(refused)?,
        ))
    }

    /// next.BREP's interference of two solids as they stand: their common volume and their gap.
    fn interference(&mut self, a: &Solid, b: &Solid) -> Result<brep_protocol::Meeting, Refusal> {
        let (a, b) = (some(a, "meet")?, some(b, "meet")?);
        let met = solid_interference(a, b, CLASSIFY).map_err(refused)?;
        Ok(brep_protocol::Meeting {
            volume: met.volume,
            gap: met.gap,
        })
    }

    /// next.BREP's least box round the solid, on its exact geometry.
    fn bounding_box(&mut self, s: &Solid) -> Result<[model::Vec3; 2], Refusal> {
        let s = some(s, "box")?;
        let [low, high] = solid_least_box(s, CLASSIFY).map_err(refused)?;
        Ok([[low.x, low.y, low.z], [high.x, high.y, high.z]])
    }

    /// next.BREP's nearest point of the solid, the point itself where the solid holds it.
    fn closest(&mut self, s: &Solid, points: &[model::Vec3]) -> Result<Vec<model::Vec3>, Refusal> {
        let s = some(s, "find a point of")?;
        let points: Vec<Vec3> = points.iter().map(|&p| v3(p)).collect();
        Ok(solid_closest_points(s, &points, CLASSIFY)
            .map_err(refused)?
            .into_iter()
            .map(|p| [p.x, p.y, p.z])
            .collect())
    }

    /// next.BREP's outward unit normal at boundary points: `None` at a crease or off the boundary.
    /// The empty solid has no boundary.
    fn normal(
        &mut self,
        s: &Solid,
        points: &[model::Vec3],
    ) -> Result<Vec<Option<model::Vec3>>, Refusal> {
        let Some(s) = s.as_ref() else {
            return Ok(vec![None; points.len()]);
        };
        let points: Vec<Vec3> = points.iter().map(|&p| v3(p)).collect();
        Ok(solid_outward_normals(s, &points, CLASSIFY)
            .map_err(refused)?
            .into_iter()
            .map(|n| n.map(|n| [n.x, n.y, n.z]))
            .collect())
    }

    /// next.BREP's principal curvatures at the boundary point nearest each point, `None` on a
    /// crease.
    fn curvature(
        &mut self,
        s: &Solid,
        points: &[model::Vec3],
    ) -> Result<Vec<Option<[f64; 2]>>, Refusal> {
        let s = some(s, "measure the curvature of")?;
        let points: Vec<Vec3> = points.iter().map(|&p| v3(p)).collect();
        solid_principal_curvatures(s, &points, CLASSIFY).map_err(refused)
    }

    /// next.BREP's conic round sized by setback, its shoulder `rho` of the way to the edge.
    fn fillet_conic(
        &mut self,
        s: &Solid,
        edges: &[model::Vec3],
        setback: f64,
        rho: f64,
    ) -> Result<Solid, Refusal> {
        let s = some(s, "fillet")?;
        let found = edge_points(s, edges)?;
        Ok(Some(
            fillet_edges_conic_setback(s, &found, None, setback, rho, None).map_err(refused)?,
        ))
    }

    /// next.BREP's display of the solid within `tolerance`: each face's triangles with its normals,
    /// each edge's polyline with the faces it bounds, and the seams, field for field.
    fn display(&mut self, s: &Solid, tolerance: f64) -> Result<brep_protocol::Display, Refusal> {
        let s = some(s, "display")?;
        let shown = solid_display(s, tolerance).map_err(refused)?;
        let p = |q: Vec3| [q.x, q.y, q.z];
        Ok(brep_protocol::Display::new(
            shown
                .faces
                .into_iter()
                .map(|f| {
                    brep_protocol::DisplayFace::new(
                        f.positions.into_iter().map(p).collect(),
                        f.normals.into_iter().map(p).collect(),
                        f.triangles,
                    )
                })
                .collect(),
            shown
                .edges
                .into_iter()
                .map(|e| {
                    brep_protocol::DisplayEdge::new(
                        e.points.into_iter().map(p).collect(),
                        e.faces,
                        None,
                        false,
                    )
                })
                .collect(),
            shown
                .seams
                .into_iter()
                .map(|seam| seam.into_iter().map(p).collect())
                .collect(),
        ))
    }

    /// next.BREP's fillet whose fillets stop `setback` short of the corner, closed by one patch.
    fn fillet_setback(
        &mut self,
        s: &Solid,
        edges: &[model::Vec3],
        radius: f64,
        setback: f64,
    ) -> Result<Solid, Refusal> {
        let s = some(s, "fillet")?;
        let found = edge_points(s, edges)?;
        Ok(Some(
            fillet_edges_setback(s, &found, None, radius, setback, None).map_err(refused)?,
        ))
    }

    /// next.BREP's curvature-continuous set-back corner: each fillet the order-2 section of tangent
    /// offset `radius` at the equal-area thumbweight, stopping `setback` short of the corner, and
    /// the corner closed by one patch with the fillets' and faces' curvature.
    fn fillet_setback_g2(
        &mut self,
        s: &Solid,
        edges: &[model::Vec3],
        radius: f64,
        setback: f64,
    ) -> Result<Solid, Refusal> {
        let s = some(s, "fillet")?;
        let found = edge_points(s, edges)?;
        Ok(Some(
            fillet_edges_setback_g2(s, &found, None, radius, setback, 0.0, None)
                .map_err(refused)?,
        ))
    }

    /// next.BREP's full round: the face `middle` lies on replaced by a round tangent to it and to
    /// the faces `side_a` and `side_b` lie on.
    fn full_round(
        &mut self,
        s: &Solid,
        side_a: model::Vec3,
        middle: model::Vec3,
        side_b: model::Vec3,
    ) -> Result<Solid, Refusal> {
        let s = some(s, "round")?;
        faces_at(s, &[side_a, middle, side_b])?;
        Ok(Some(
            full_round(s, v3(side_a), v3(middle), v3(side_b), None).map_err(refused)?,
        ))
    }

    /// next.BREP's fillet between two faces by a ball of `radius` touching both.
    fn fillet_faces(
        &mut self,
        s: &Solid,
        face_a: model::Vec3,
        face_b: model::Vec3,
        radius: f64,
    ) -> Result<Solid, Refusal> {
        let s = some(s, "fillet")?;
        faces_at(s, &[face_a, face_b])?;
        Ok(Some(
            fillet_faces(s, v3(face_a), v3(face_b), radius, None).map_err(refused)?,
        ))
    }

    /// next.BREP's thickening of the faces the points lie on, along their outward normals.
    fn thicken(
        &mut self,
        s: &Solid,
        faces: &[model::Vec3],
        thickness: f64,
    ) -> Result<Solid, Refusal> {
        let s = some(s, "thicken")?;
        faces_at(s, faces)?;
        let points: Vec<Vec3> = faces.iter().map(|&p| v3(p)).collect();
        Ok(solid(
            thicken_faces(s, &points, thickness).map_err(refused)?,
        ))
    }

    /// next.BREP's hidden-line drawing looking along `direction`, each chord within `tolerance`.
    fn project(
        &mut self,
        s: &Solid,
        direction: model::Vec3,
        tolerance: f64,
    ) -> Result<brep_protocol::Projection, Refusal> {
        let Some(s) = s.as_ref() else {
            return Ok(brep_protocol::Projection::default());
        };
        let drawn = brep_render::sheets::hlr::hidden_line_drawing(s, direction, tolerance)
            .map_err(refused)?;
        Ok(brep_protocol::Projection::new(drawn.visible, drawn.hidden))
    }

    /// next.BREP's section of the boundary by a plane, each piece an edge of the cut sampled within
    /// `tolerance`.
    fn section(
        &mut self,
        s: &Solid,
        origin: model::Vec3,
        normal: model::Vec3,
        tolerance: f64,
    ) -> Result<Vec<Vec<model::Vec3>>, Refusal> {
        let Some(s) = s.as_ref() else {
            return Ok(Vec::new());
        };
        Ok(solid_section(s, v3(origin), v3(normal), tolerance)
            .map_err(refused)?
            .into_iter()
            .map(|piece| piece.into_iter().map(|p| [p.x, p.y, p.z]).collect())
            .collect())
    }

    /// next.BREP's area of the cut by a plane.
    fn section_area(
        &mut self,
        s: &Solid,
        origin: model::Vec3,
        normal: model::Vec3,
    ) -> Result<f64, Refusal> {
        let Some(s) = s.as_ref() else { return Ok(0.0) };
        solid_section_area(s, v3(origin), v3(normal)).map_err(refused)
    }

    /// next.BREP's sewing of a STEP file's faces into one solid.
    fn sew(&mut self, text: &str, tolerance: f64) -> Result<Solid, Refusal> {
        Ok(solid(sew_step(text, tolerance).map_err(refused)?))
    }

    /// next.BREP's fillet under a radius law on the edge's arc length, here one linear segment from
    /// one radius to the other. A single edge's law runs from the edge's own start (`t0`), so the
    /// radius named for the end nearest `from` goes first when that end is the start.
    fn fillet_variable(
        &mut self,
        s: &Solid,
        edge: model::Vec3,
        from: model::Vec3,
        r_from: f64,
        r_to: f64,
    ) -> Result<Solid, Refusal> {
        let s = some(s, "fillet")?;
        let id = edge_at(s, edge)?;
        let record = s
            .edges
            .iter()
            .find(|e| e.id == id)
            .ok_or_else(|| Refusal::invalid("the edge named has no record"))?;
        let start = record.curve.evaluate(record.t0).map_err(failed)?;
        let end = record.curve.evaluate(record.t1).map_err(failed)?;
        let near = v3(from);
        let (first, last) = if near.sub(start).length() <= near.sub(end).length() {
            (r_from, r_to)
        } else {
            (r_to, r_from)
        };
        let law = RadiusLaw::from_segments(&[LawSegment::Linear {
            length: 1.0,
            start_radius: first,
            end_radius: last,
        }])
        .map_err(refused)?;
        Ok(Some(
            fillet_edges_variable_law(s, &[v3(edge)], None, &law, false, None).map_err(refused)?,
        ))
    }

    /// next.BREP's fillet of several edges in one operation, each edge at its own radius.
    fn fillet_radii(
        &mut self,
        s: &Solid,
        edges: &[model::Vec3],
        radii: &[f64],
    ) -> Result<Solid, Refusal> {
        let s = some(s, "fillet")?;
        let found = edge_points(s, edges)?;
        Ok(Some(
            fillet_edges_radii(s, &found, None, radii, false, None).map_err(refused)?,
        ))
    }

    /// next.BREP's curvature-continuous (order 2) fillet: tangent offset and centre radius.
    fn fillet_g2(
        &mut self,
        s: &Solid,
        edges: &[model::Vec3],
        offset: f64,
        center: f64,
    ) -> Result<Solid, Refusal> {
        let s = some(s, "fillet")?;
        let found = edge_points(s, edges)?;
        Ok(Some(
            fillet_edges_g2(s, &found, None, offset, center, None).map_err(refused)?,
        ))
    }

    /// next.BREP's chamfer of several edges in one operation, each edge at its own distance.
    fn chamfer_distances(
        &mut self,
        s: &Solid,
        edges: &[model::Vec3],
        distances: &[f64],
    ) -> Result<Solid, Refusal> {
        let s = some(s, "chamfer")?;
        let found = edge_points(s, edges)?;
        Ok(Some(
            fillet_edges_radii(s, &found, None, distances, true, None).map_err(refused)?,
        ))
    }

    /// next.BREP's two-distance chamfer, each distance set back along its face, both alike.
    fn chamfer(&mut self, s: &Solid, edges: &[model::Vec3], d: f64) -> Result<Solid, Refusal> {
        let s = some(s, "chamfer")?;
        let found = edge_points(s, edges)?;
        let options = brep_kernel::ChamferOptions {
            follow_tangent_chains: true,
        };
        Ok(Some(
            brep_kernel::chamfer_edges_with_options(s, &found, None, d, options, None)
                .map_err(refused)?,
        ))
    }

    /// next.BREP sets `d1` back along the first of the edge's two faces in its own order, so the
    /// named face's distance goes first when the named face is that one.
    fn chamfer_two(
        &mut self,
        s: &Solid,
        edge: model::Vec3,
        face: model::Vec3,
        d1: f64,
        d2: f64,
    ) -> Result<Solid, Refusal> {
        let s = some(s, "chamfer")?;
        let id = edge_at(s, edge)?;
        let named = faces_at(s, &[face])?[0];
        let (_, faces) = brep_blend::fillet::analyze::edge_and_faces(s, id).map_err(refused)?;
        let (a, b) = match faces.as_slice() {
            [first, _] if first.id == named => (d1, d2),
            [_, second] if second.id == named => (d2, d1),
            _ => {
                return Err(Refusal::invalid(
                    "the face named is not one of the two faces at the edge named",
                ));
            }
        };
        Ok(Some(
            chamfer_edges_asymmetric(s, &[v3(edge)], None, a, b, None).map_err(refused)?,
        ))
    }

    /// Every point within `distance` of the solid: next.BREP's offset with rounded corners.
    fn offset(&mut self, s: &Solid, distance: f64) -> Result<Solid, Refusal> {
        let s = some(s, "offset")?;
        Ok(Some(
            offset_solid(s, distance, OffsetCorners::Rounded)
                .map_err(refused)?
                .solid,
        ))
    }

    /// The wall every point of the solid within `t` of its kept faces: next.BREP's shell joined
    /// by arcs, inward for a positive distance.
    fn shell(&mut self, s: &Solid, open: &[model::Vec3], t: f64) -> Result<Solid, Refusal> {
        let s = some(s, "shell")?;
        let ids = faces_at(s, open)?;
        // Open: next.BREP's shell, its walls offset and re-joined sharp (the protocol's walls
        // `t` thick, meeting at the source's own edges). Closed (no opening): the rounded
        // hollow, the points within `t` of the boundary.
        if ids.is_empty() {
            return Ok(Some(offset_shell_rounded(s, &ids, t).map_err(refused)?));
        }
        // Walls meeting across the whole part leave no cavity: refused, as the protocol's shell
        // must leave one.
        let options = brep_kernel::ShellOptions {
            refuse_no_cavity: true,
        };
        Ok(Some(
            brep_kernel::offset_shell_with_options(s, &ids, t, options)
                .map_err(refused)?
                .solid,
        ))
    }

    fn draft(
        &mut self,
        s: &Solid,
        faces: &[model::Vec3],
        origin: model::Vec3,
        pull: model::Vec3,
        angle: f64,
    ) -> Result<Solid, Refusal> {
        let s = some(s, "draft")?;
        let ids = faces_at(s, faces)?;
        Ok(Some(
            draft_faces(s, &ids, v3(pull), v3(origin), angle).map_err(refused)?,
        ))
    }

    fn delete_faces(&mut self, s: &Solid, faces: &[model::Vec3]) -> Result<Solid, Refusal> {
        let s = some(s, "delete faces from")?;
        let ids = faces_at(s, faces)?;
        Ok(Some(delete_faces_and_heal(s, &ids).map_err(refused)?))
    }

    /// next.BREP offsets faces onto their offset carriers in one call (its Push Face push).
    fn offset_faces(
        &mut self,
        s: &Solid,
        faces: &[model::Vec3],
        distance: f64,
    ) -> Result<Solid, Refusal> {
        let s = some(s, "offset faces of")?;
        faces_at(s, faces)?;
        let points: Vec<Vec3> = faces.iter().map(|&p| v3(p)).collect();
        Ok(Some(
            brep_kernel::offset_faces(s, &points, distance).map_err(refused)?,
        ))
    }

    /// next.BREP moves faces by a translation, or turns them about an axis: a motion that does
    /// both, turning about an axis while sliding along it, has no one call.
    fn move_faces(
        &mut self,
        s: &Solid,
        faces: &[model::Vec3],
        m: [[f64; 4]; 3],
    ) -> Result<Solid, Refusal> {
        let s = some(s, "move faces of")?;
        let ids = faces_at(s, faces)?;
        let identity = (0..3).all(|i| (0..3).all(|j| m[i][j] == if i == j { 1.0 } else { 0.0 }));
        if identity {
            let points: Vec<Vec3> = faces.iter().map(|&p| v3(p)).collect();
            return Ok(Some(
                brep_kernel::transform_faces(s, &points, m).map_err(refused)?,
            ));
        }
        let Some((point, axis, angle)) = turn(m) else {
            // A turn and a slide at once: next.BREP solves the whole rigid motion in one call.
            let points: Vec<Vec3> = faces.iter().map(|&p| v3(p)).collect();
            return Ok(Some(
                brep_kernel::transform_faces(s, &points, m).map_err(refused)?,
            ));
        };
        Ok(Some(
            rotate_faces(s, &ids, point, axis, angle).map_err(refused)?,
        ))
    }

    /// next.BREP's plane split gives the piece below the plane, then the one above it: the
    /// protocol's order is the normal's side first.
    fn split(
        &mut self,
        s: &Solid,
        origin: model::Vec3,
        normal: model::Vec3,
    ) -> Result<(Solid, Solid), Refusal> {
        let Some(s) = s else { return Ok((None, None)) };
        let (below, above) = split_solid_by_plane(s, v3(origin), v3(normal)).map_err(failed)?;
        Ok((solid(above), solid(below)))
    }

    fn mesh(&mut self, s: &Solid, tolerance: f64) -> Result<brep_protocol::Mesh, Refusal> {
        let Some(s) = s else {
            return Ok(brep_protocol::Mesh::default());
        };
        let m = tessellate_brep_watertight(s, tolerance).map_err(failed)?;
        Ok(brep_protocol::Mesh::new(
            m.positions
                .as_chunks::<3>()
                .0
                .iter()
                .map(|p| [p[0], p[1], p[2]])
                .collect(),
            m.indices
                .as_chunks::<3>()
                .0
                .iter()
                .map(|t| [t[0], t[1], t[2]])
                .collect(),
        ))
    }

    fn write_step(&mut self, s: &Solid) -> Result<String, Refusal> {
        let s = some(s, "write as STEP")?;
        export_step(
            std::slice::from_ref(s),
            "solid",
            "mm",
            "1970-01-01T00:00:00",
        )
        .map_err(failed)
    }

    fn read_step(&mut self, text: &str) -> Result<Solid, Refusal> {
        let mut solids = import_step(text).map_err(failed)?;
        match solids.len() {
            0 => Err(Refusal::invalid("the file holds no solid")),
            1 => Ok(solid(solids.remove(0))),
            n => Err(Refusal::unsupported(format!(
                "the file holds {n} solids, and next.BREP reads each as a solid of its own"
            ))),
        }
    }

    /// What next.BREP 0.32.2 has one call of the protocol's shape for. It has no such call for ray
    /// casts, distances or Betti numbers, so those are left out with the operations it does not
    /// have.
    fn supports(&self) -> Vec<String> {
        [
            "box",
            "cylinder",
            "sphere",
            "torus",
            "cone",
            "extrude",
            "revolve",
            "loft",
            "loft_through",
            "sweep",
            "sweep_path",
            "union",
            "intersection",
            "difference",
            "union_all",
            "difference_all",
            "fillet",
            "fillet_variable",
            "fillet_radii",
            "fillet_g2",
            "fillet_faces",
            "full_round",
            "fillet_setback",
            "fillet_setback_g2",
            "display",
            "chamfer_distances",
            "closest",
            "interference",
            "slab",
            "bounding_box",
            "normal",
            "curvature",
            "sew",
            "thicken",
            "fillet_conic",
            "screw",
            "section",
            "section_area",
            "project",
            "chamfer",
            "chamfer_two",
            "offset",
            "delete_faces",
            "move_faces",
            "offset_faces",
            "shell",
            "draft",
            "write_step",
            "read_step",
            "volume",
            "area",
            "classify",
            "check",
            "centroid",
            "inertia",
            "split",
            "mesh",
            "transform",
            "affine",
        ]
        .map(String::from)
        .to_vec()
    }
}

fn main() -> ExitCode {
    match brep_protocol::serve(&mut NextBrep) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{e}");
            ExitCode::FAILURE
        }
    }
}

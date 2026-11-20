# brep-conformance-next-brep

The [brep-conformance](https://github.com/springcraft-ai/brep-conformance) adapter for
[next.BREP](https://github.com/mmiscool/next.BREP.io_RUST_BREP_KERNEL), the kernel published on
crates.io as [`BREP_kernel`](https://crates.io/crates/BREP_kernel): it puts next.BREP behind the
suite's protocol, so that the suite can score it.

## Who maintains this

brep-conformance scores a B-rep kernel through an adapter: a small program that answers the suite's
requests with the kernel's own calls. Springcraft, which maintains the suite, wrote this adapter and
maintains it.

- **It is not next.BREP's.** next.BREP's maintainers did not write it or choose the tests it runs.
- **A failure can be the adapter's.** A case next.BREP fails in the results may be this adapter
  driving it wrongly. Corrections are welcome: open an issue or a pull request.
- **It adds nothing to the kernel.** The adapter is pass-through, as the suite's
  [contract](https://github.com/springcraft-ai/brep-conformance/blob/main/doc/contract.md#what-an-adapter-may-do)
  asks: one call of the kernel for each function, and anything the kernel has no one call for is
  listed as unsupported. What is scored is next.BREP, not this code.
- **next.BREP's maintainers are welcome here.** They know best how to drive next.BREP, so their
  review and corrections carry the most weight. Ask in an issue to be added as a reviewer.
- **The suite decides what is scored.** An official run uses only a commit of this adapter that the
  suite's maintainers have reviewed as pass-through, whoever wrote it.

## Overview

[`src/main.rs`](src/main.rs) answers each function with next.BREP's one call for it, as
[the contract](https://github.com/springcraft-ai/brep-conformance/blob/main/doc/contract.md#what-an-adapter-may-do)
asks, and lists the rest as unsupported. It builds the suite's box, cylinder, sphere, torus and cone
with next.BREP's own constructors, and its extrusions, revolutions, lofts and sweeps with its
profile constructions, each about the origin and moved onto its frame by next.BREP's affine map. It
answers the booleans, the booleans of many, the mass properties, the validity check,
classification, rigid and affine maps, fillets, chamfers, offsets, shells, drafts, face deletes and
moves, plane splits, meshes and STEP with next.BREP's own calls.

next.BREP takes an edge or a face by a point and finds the nearest within about a millimetre. The
contract names an edge or a face only within `10·ε`, so the adapter first holds each point to that,
by next.BREP's own projection onto the edges and its own classification of the point on the
boundary, and refuses a point that names nothing. A point on an edge, where several faces meet,
names no one face, and the call is unsupported.

A two-distance chamfer sets its first distance back along the first of the edge's faces in
next.BREP's own order, read from `BREP_blend`, the same release's blend crate; the adapter gives
the named face's distance in that place. A face move is next.BREP's translation of faces, or its
turn of faces about an axis; a move that turns and slides at once is unsupported.

next.BREP says a result is empty by a solid with no faces, and the adapter reads that, and nothing
else, as the empty solid. Its refusals keep their class: geometry next.BREP does not model yet is
unsupported, an input it rejects or that has no single answer is invalid, and anything else is a
failure.

The adapter classifies points at the suite's tolerance, 1e-7 mm, and says so in `hello`. It
declares 4e-3 mm, the edge-to-face disagreement next.BREP's validator accepts in a valid solid.

Listed as unsupported, since next.BREP 0.24.0 has no public call of the protocol's shape for it:

- Face offsets, thickening a solid's faces, and outward normals: next.BREP builds these, but behind
  its feature pipeline or in private modules.
- A variable fillet: next.BREP's runs its radius from the edge's own start, and cannot be told the
  end the protocol names.
- Screws (its helical sweep keeps the section square to the helix), slabs, conic fillets sized by
  setback, fillets and chamfers of several sizes at once, curvature-continuous, setback, face and
  full-round fillets, sewing, the nearest point, curvatures, distance, interference, Betti numbers,
  sections, projections, the least bounding box, displays and ray casts.

next.BREP's profile constructions take one loop, so a profile with holes is unsupported: next.BREP
cuts holes as booleans of their own. A profile's periodic spline is next.BREP's periodic B-spline
on the period's points, which takes evenly spaced knots, so a spline on other knots is
unsupported; a conic is its rational quadratic through the two ends and the apex. A blend-section
segment has no next.BREP curve and is unsupported.

## Setup

The Rust toolchain [`rust-toolchain.toml`](rust-toolchain.toml) pins.

## Usage

The suite runs the adapter as a separate process. Build the adapter here, then run the suite from a
checkout of [brep-conformance](https://github.com/springcraft-ai/brep-conformance) beside this one:

```bash
cargo build --release
cd ../brep-conformance
cargo build --release -p brep-conformance
target/release/brep-conformance \
    --adapter ../brep-conformance-next-brep/target/release/next-brep-adapter --out next-brep.json
```

## Licences

The adapter in this repository is under the MIT licence: see [LICENSE](LICENSE).

That licence covers the adapter and nothing else. It does not cover the kernel, and nothing here
relicenses the kernel.

next.BREP is a separate work under its own licence, shipped as `LICENSE.md` in each of its crates;
its repository is the authority on that, not this page. This adapter uses the crates as published
and changes none of them.

Building the adapter downloads next.BREP's crates from crates.io and links them in, so a built
adapter contains next.BREP and is subject to next.BREP's terms as well as this repository's.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

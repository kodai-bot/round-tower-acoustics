# Agent instruction — full rebuild of the Clondalkin tower model

**This supersedes all previous geometry.** An RIA archive cross-section of Clondalkin has
shown the earlier model used a diameter roughly twice the true internal value — wall
thickness appears to have been omitted. Every frequency-domain result produced to date is
void until this rebuild is complete.

Please treat this as a clean build, not a patch. Archive the old outputs rather than
overwriting them, and tag all new outputs so the two sets can never be confused.

---

## 1. Corrected geometry (bell storey / top chamber)

| Parameter | Imperial | Metric | Source |
|---|---|---|---|
| Internal diameter | 6'-6" | **1.981 m** (r = 0.991 m) | RIA drawing, legible |
| Chamber height | 9'-6" | **2.896 m** | RIA drawing |
| Cone height above chamber | 6'-3" | **1.905 m** | RIA drawing |
| Wall thickness at bell storey | — | **0.745 m** (derived) | see note |
| External diameter at base | — | **4.043 m** | published circumference 12.7 m |
| External diameter at bell storey | — | **3.472 m** (estimated) | taper estimate — see note |
| Window opening height | 3'-7" | **1.092 m** | RIA drawing |
| Window opening width | 2'-0" | **0.610 m** | Dublin CC 3D model (external face) |
| Window reveal depth | = wall thickness | **0.745 m** | derived |
| Window sill above floor | ~4'-0" | **1.219 m** (approx) | estimated from drawing |
| Floor area | — | **3.08 m²** | derived |

**How the external dimensions were derived.** Published sources give Clondalkin a base
circumference of 12.7 m, i.e. an external diameter of **4.043 m** — consistent with the
independently quoted 4.04–4.05 m "widest point." Round towers taper; applying a taper
comparable to Devenish (4.82 m at base to 4.14 m at the top windows, a 14% reduction)
gives an external diameter at the bell storey of about **3.47 m**. Against the RIA
drawing's internal diameter of 1.981 m, that implies a wall thickness of **0.745 m** —
which sits between the 0.64 m secondary estimate and the 0.813 m reading of "2'-8"" as a
reveal depth, and is preferred because it is derived from two independent measurements
rather than from either alone. **Run 0.64 m and 0.813 m as sensitivity cases.**

Note also that published heights vary (25.6 m, 27.5 m) against the drawing's 86'-9"
(26.44 m) — probably different datums. Not critical for the bell-storey model.

Windows: **four, cardinal** (independently confirmed: N/S/E/W). Model each as a
rectangular duct 0.610 × 1.092 m through the full 0.745 m wall depth, not as a thin
aperture.

**Window subtense — quote both faces, they are not the same thing:**

| Face | Radius | Subtense each | Four windows |
|---|---|---|---|
| External | 1.736 m | 20.2° | 80.9° = **22% of circumference** |
| Internal | 0.991 m | 35.8° | 143.3° = **40% of circumference** |

This is the *same* physical aperture: the fraction differs only because the internal
circumference is smaller. It fully explains the observation that the windows "look bigger
from the inside" **without requiring any splay at all.** An earlier version of this
document quoted the 40% figure alone, which overstated how open the tower appears
externally — the photographs show narrow openings, and 22% is the number consistent with
them.

**Splay:** possibly a slight inward taper, but given the above the visual impression is
explained by geometry, so there is no good evidence for splay in either direction. Model
straight-sided and parallel; treat splay as unresolved rather than as a finding.

Note Clondalkin is documented as the most slender surviving tower — do not generalise
this radius to other towers.

---

## 2. Build two mesh variants

The drawing marks corbels at the base of the cone — at the ceiling level, not above the
window lintels — so **Variant A is the primary run.** Build B as the comparison case.

- **Variant A — CEILED (PRIMARY)**: rigid flat ceiling at z = 2.896 m. Cone void excluded
  from the acoustic domain. Run the full pipeline on this.
- **Variant B — OPEN (comparison)**: chamber continuous with the conical void above
  (total internal height 2.896 + 1.905 = 4.801 m at the axis). Frequency sweep only is
  sufficient unless it diverges sharply from A.

If cheap to add, a **Variant A-timber** with a finite-impedance ceiling rather than rigid
would be informative, since a timber floor is partially transmitting and absorptive.

**Note the proportions this produces.** With the ceiling at 2.896 m and window heads at
2.311 m, there is only **0.585 m (about 2 ft) of wall above the window lintels.** The
chamber is a short, wide-mouthed canister: 1.98 m across, 2.90 m tall, with four openings
occupying 40% of its internal wall surface and a hard lid close above them. This is a very
different acoustic object from the tall open shaft previously modelled, and the ceiling is
close enough to the apertures that it should materially affect coupling.

---

## 3. Source

Single monopole on the chamber axis at **z = 1.50 m** (standing adult mouth height).
This falls inside the window band (sill 1.219 m, head 2.311 m), which is a change from
the previous configuration.

Floor area of 3.08 m² permits **one occupant only** — do not model paired or multiple
sources for this tower. The centred-source assumption is now physically justified rather
than assumed.

---

## 4. Validation targets before trusting any output

Analytic rigid-cylinder modes for r = 0.991 m, H = 2.896 m. The BEM interior solution
should show features near these. **If it doesn't, stop and debug before proceeding.**

Transverse (independent of ceiling condition):

| Mode | Frequency |
|---|---|
| (1,0) 1st azimuthal | 101.5 Hz |
| (2,0) 2nd azimuthal | 168.3 Hz |
| (0,1) 1st radial | 211.2 Hz |
| (3,0) | 231.5 Hz |
| (4,0) | 293.0 Hz |

Axial, Variant A (closed both ends): 59.2, 118.5, 177.7, 236.9, 296.1 Hz
Axial, Variant B (quarter-wave approximation): 29.6, 88.8, 148.1 Hz

The difference in the axial series is the cleanest discriminator between the two
variants — worth reporting explicitly.

---

## 5. Flags to resolve or report

Window width and wall thickness are now both resolved to reasonable confidence. Report
subtense at **both** faces in the run log (22% external, 40% internal) — quoting only one
misrepresents the geometry.

Remaining uncertainties, in order of consequence:

1. **Wall thickness** — 0.745 m is derived from a taper estimate, not measured at the
   bell storey. Run 0.64 m and 0.813 m as bounds.
2. **Sill height** — estimated at ~4'-0", not read from the drawing. This sets how much
   of the singer's output falls within the aperture band, so it matters.
3. **Ceiling material** — corbels indicate a floor was there; whether timber or stone is
   unknown. Rigid is the baseline, timber the sensitivity case.
4. **Splay** — unresolved, and the visual impression that prompted the question turns out
   to be explained by circumference geometry alone. Do not treat as evidence either way.

---

## 6. Pipeline to re-run

Everything downstream of geometry, in order:

1. Interior solve — window driving conditions (magnitude **and phase**, as before)
2. Frequency sweep, 34 points as previously, **extended upward to at least 500 Hz** —
   modes have moved up and the old ceiling of 390 Hz may now truncate the region of
   interest
3. Ground projection — **36 azimuths at 10° spacing**, not 8. Include near distances
   **15, 25, 30 m** (real tower-to-church distances at Clonmacnois, Kells, Glendalough)
   alongside the existing 50/100/200/500/1000/2000 m grid
4. Pilgrim ring walk — both rings, 36 points
5. Elevation sweep — report the lower/upper hemisphere split with cos-weighted
   linear-intensity averaging, and state the method in the output header
6. Advantage-zone overlay on the site plan
7. Height comparison (ground / platform / tower)

---

## 7. Outputs

Prefix everything `v2_` and include a header line in each file stating: geometry version,
mesh variant (A / B / A-timber), window width used, and the script that produced it.

Please also produce a short `v2_validation.md` comparing the BEM interior results against
the analytic modes in §4, and a one-page diff summarising how the headline numbers moved
against the old (void) results.

---

## 8. Still outstanding from the previous round

Carry these forward — they were not resolved and remain relevant:

- Reciprocity check on the two boundary-integral formulations
- Direct vs CFIE cross-check at whichever frequency now shows the sharpest feature
- Crowd layer as a finite-impedance boundary (Shabalina equivalent-fluid parameters),
  replacing the borrowed seat-dip curve
- Water as a reflecting plane (method of images first, meshed boundary only if needed)

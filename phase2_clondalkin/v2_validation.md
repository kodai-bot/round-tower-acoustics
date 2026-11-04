# v2 validation gate -- interior-A analytic mode comparison

Probe point: [-0.563815572471545, -0.2052120859954012, 0.8] (off-axis, off-source)

Analytic-target self-check (transcribed vs independently computed): PASS

| Target | Target Hz | Nearest detected peak Hz | Delta Hz | Tolerance Hz | Result |
|---|---|---|---|---|---|
| axial n=1 | 59.2 | 61.0 | 1.8 | 8.0 | PASS |
| (1,0) | 101.5 | 101.0 | 0.5 | 8.0 | PASS |
| axial n=2 | 118.5 | 117.0 | 1.5 | 8.0 | PASS |
| (2,0) | 168.3 | 169.0 | 0.7 | 8.4 | PASS |
| axial n=3 | 177.7 | 177.0 | 0.7 | 8.9 | PASS |
| (0,1) | 211.2 | 205.0 | 6.2 | 10.6 | PASS |
| (3,0) | 231.5 | 237.0 | 5.5 | 11.6 | PASS |
| axial n=4 | 236.9 | 237.0 | 0.1 | 11.8 | PASS |
| (4,0) | 293.0 | 293.0 | 0.0 | 14.7 | PASS |
| axial n=5 | 296.1 | 293.0 | 3.1 | 14.8 | PASS |

**Overall: PASS**

Probe sweep: 45-309Hz, step 4.0Hz, 67 points, mesh v2_mesh_interior_A.msh (r_inner=0.9905m).

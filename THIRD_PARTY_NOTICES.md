# Third-Party Notices

This project is licensed under the MIT License (see `LICENSE`). It also
incorporates or adapts third-party material as listed below. Each component
remains under its own license; those licenses are reproduced here as required.

---

## NumPy (BSD-3-Clause)

The following code is adapted from the NumPy project
(<https://github.com/numpy/numpy>, `numpy/random/src/distributions/`):

- the ziggurat standard-normal sampler (`ran_standard_normal`),
- the ziggurat standard-exponential sampler (`ran_standard_exponential`),
- the gamma sampler (`ran_gamma`, Marsaglia–Tsang / Ahrens–Dieter),
- the ziggurat constant tables in `src/rng.h` (`ki_f64`, `wi_f64`, `fi_f64`,
  `ki_f32`, `wi_f32`, `fi_f32`).

```
Copyright (c) 2005-2024, NumPy Developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
       copyright notice, this list of conditions and the following
       disclaimer in the documentation and/or other materials provided
       with the distribution.

    * Neither the name of the NumPy Developers nor the names of any
       contributors may be used to endorse or promote products derived
       from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

---

## raddebugger (MIT)

The base layer in `src/base.h` (type aliases, the `internal`/`global`/
`read_only` keyword macros, and assert helpers) adapts conventions from
raddebugger by Epic Games / Ryan Fleury
(<https://github.com/EpicGamesExt/raddebugger>).

```
MIT License

Copyright (c) 2024 Epic Games Tools

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## xoshiro256** and splitmix64 (public domain, CC0)

The pseudo-random generator (`ran_u64`) and its seeding (`splitmix64`) in
`src/rng.c` are by David Blackman and Sebastiano Vigna, released to the public
domain (CC0). See <https://prng.di.unimi.it/>.

> To the extent possible under law, the author has dedicated all copyright and
> related and neighboring rights to this software to the public domain
> worldwide. This software is distributed without any warranty.
> See <http://creativecommons.org/publicdomain/zero/1.0/>.

---

## Algorithms implemented from published descriptions

The following are independent implementations written from the cited
publications. Algorithms are not subject to copyright; these references are
provided for attribution and provenance.

- **Uniform double generation** (`ran_f64_downey`): Allen B. Downey,
  "Generating Pseudo-random Floating-Point Values" (2007),
  <https://allendowney.com/research/rand/>.
- **Poisson PTRD** (`ran_poisson`, large-lambda branch): W. Hörmann,
  "The transformed rejection method for generating Poisson random variables",
  *Insurance: Mathematics and Economics* 12 (1993) 39–45.
- **Regularized incomplete gamma** (`gamma_p`, test code only): standard
  series / continued-fraction evaluation per the NIST DLMF, ch. 8
  (<https://dlmf.nist.gov/8>), using the modified Lentz algorithm
  (W. J. Lentz, 1976). Cf. the Cephes library (BSD) `igam`/`igamc`.
- **Gauss–Legendre nodes/weights** (`compute_gauss_legendre_roots_and_weights`):
  Newton iteration on the Legendre three-term recurrence with the Tricomi
  asymptotic initial guess; standard numerical-analysis material.

# Heston Pricing

This is a personal project that I started as a first computationally heavy project in C in order to learn. The goal is to implement different methods of pricing in the heston model, perform rigorous benchmarks and show failure modes.

The heston Model:

```math
\begin{aligned}
dS_t &= (r-q)\,S_t\,dt + \sqrt{v_t}\,S_t\,dW^S_t \\
dv_t &= \kappa(\theta - v_t)\,dt + \sigma\sqrt{v_t}\,dW^v_t \\
d\langle W^S, W^v\rangle_t &= \rho\,dt
\end{aligned}
```

I try to implement various methods for pricing european derivatives with the heston model. Mainly:

- [Lewis (2001)](#lewis-european-call-formula). 🚧
- COS (Fang–Oosterlee 2008). 🚧
- Carr–Madan (1999). 🚧
- Monte Carlo (Andersen QE 2008). 🚧


### Lewis European call formula

$$
C = F - \frac{\sqrt{F K}}{\pi}
\int_0^{\infty}\mathrm{Re}\left[\frac{e^{i u X}\phi_T(u - i/2)}{u^2 + \tfrac14}\right] du, \quad\quad X=\log{F/K}
$$
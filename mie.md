#Mie code Equations
This document lays out the math behind Mätzler's Mie scattering codes.

## mie()
Computation of Mie efficiencies for a single sphere.
### Mätzler's comment
Computation of Mie Efficiencies for give complex refractive-index ration m=m'+im" and size parameter x=k0*a, where k0= wave number in ambient medium, a=sphere radius, using complex Mie Coefficients an and bn for n=1 to nmax,
s. Bohren and Huffman (1983) BEWI:TDD122, p. 103, 119-122,477.
Result: m', m", x, efficiencies for extinction (qext), scattering (qsca), absorption (qabs), backscattering (qb), asymmetry parameter (asy=<costeta>) and (qration=qb/qsca).
Uses the function "Mie_ab" for an and bn, for n=1 to nmax.
C. Mätzler, May 2002, revised July 2002.

### Explanation
Starting with some simple definitions:
$$x=k_0a$$
$$k_0={{2\pi}\over{\lambda}}$$
with a as the radius of the sphere.
Given the complex relative refractive index $m={k \over k_0}={N \over N_0}$ where $N$ and $N_1$ are the refractive indices of the particle and medium (here, vacuum) respectively, the scattering coefficiences are:

$$a_n = {{\mu_0 m^2j_n(mx)\left[xj_n(x)\right]'-\mu j_n(x)\left[mxj_n(mx)\right]'}\over{\mu_0 m^2j_n(mx)\left[xh_{n}^{(1)}(x)\right]'-\mu h_{n}^{(1)}(x)\left[mxj_n(mx)\right]'}}$$

$$b_n={{\mu j_n(mx)\left[xj_n(x)\right]'-\mu_0 j_n(x)\left[mxj_n(mx)\right]'}\over{\mu j_n(mx)\left[xH_{n}^{(1)}(x)\right]'-\mu_0 h_{n}^{(1)}(x)\left[mxj_n(mx)\right]'}}$$

Here, $j_n$ is related to the Bessel function of the first kind ($J_n$) by

$$j_n(\rho)=\sqrt{\pi \over {2\rho}}J_{n+1/2}(\rho)$$

and $h_{n}^{(1)}$ is the spherical Hankel function (also called the spherical Bessel function of the third kind):

$$h_n^{(1)}=j_n(\rho)+iy_n(\rho)$$

which relies on the spherical Bessel function of the second kind ($Y_n$)through

$$y_n(\rho)=\sqrt{\pi \over {2\rho}}Y_{n+1/2}(\rho).$$

$\mu_0$ and $\mu$ are the magnetic permeability of the medium and sphere respectively.
If we take the permeability of the sphere and the medium to be the same (as is the case with a non-magnetic sphere in a vacuum), then

$$a_n={{m\psi_n(mx)\psi'(x)-m\psi_n(x)\psi_n'(mx)}\over{m\psi_n(mx)\xi_n'(x)-m\xi_n(x)\psi_n'(mx)}}$$

$$b_n={{\psi_n(mx)\psi_n'(x)-\psi_n(x)\psi_n'(mx)}\over{\psi_n(mx)\xi_n'(x)-m\xi_n(x)\psi_n'(mx)}}$$

where,

$$\psi_n(\rho)=\rho j_n(\rho)$$

$$\xi_n(\rho)=\rho h_n^{(1)}(\rho).$$

Note that $a_n$ and $b_n$ vanish as m approaches unity, as is expected.

Fortunately, this program expects to receive $a_n$ and $b_n$ from an external source, so we need only show how they are used to calculate the efficiencies.
The first question to answer in doing this is, "how many terms are required for convergence?" Thankfully, Bohren and Huffman provide an answer and a citation to Wiscombe (1979, 1980):

$$n_{max}=x+4x^{1/3}+2.$$

Mätzler differs from B&H's notation using `nmax` instead of `NSTOP`. He also precomputes a number of terms that appear in later calculations.

#### $Q_{ext}$
From B&H eq. 4.62:

$$C_{ext}={{2\pi}\over{k^2}} \sum_{n=1}^\infty{(2n+1)Re\{a_n + b_n\}}$$

and the definition of $Q_{ext}$:

$$Q_{ext}={C_{ext}\over G}$$

where,

$$G=\pi a^2$$

for a sphere of radius $a$.

The term before the summation becomes,
$${{2\pi}\over{k^2 \pi a^2}}.$$
Simplifying, and recalling the definition of $x$ yields,
$$Q_{ext}={2\over x^2}\sum_{n=1}^\infty{(2n+1)Re\{a_n + b_n\}}$$

#### $Q_{sca}$
Similar to $Q_{ext}$ above, B&H eq. 4.61:

$$C_{sca} = {{2\pi}\over{k^2}} \sum_{n=1}^\infty {(2n+1)(|a_n|^2+|b_n|^2)}$$

and the definition of $Q_{sca}$:

$$Q_{sca} = {C_{sca}\over G}.$$

Substituting and simplifying as above,

$$Q_{sca}={2\over x^2}\sum_{n=1}^\infty {(2n+1)(|a_n|^2+|b_n|^2)}$$

#### $Q_{abs}$
Since we already have $Q_{sca}$ and $Q_{ext}$ we can simply calculate this as,

$$Q_{abs} = Q_{ext}-Q_{sca}$$

#### $Q_b$
The top equation on B&H p. 122 provides:

$$Q_b = {1 \over x^2}\left|\sum_n (2n+1)(-1)^n(a_n-b_n) \right|^2$$

#### Asymmetry Parameter
The second equation on B&H p. 120 provides:

$$\left<cos \theta\right > = {4\over {x^2Q_{sca}}} \left[\sum_n{{{n(n+2)}\over{n+1}}Re\{a_n a_{n+1}^* + b_n b_{n+1}^* \}} \\
+ \sum_n {{{2n+1}\over{n(n+1)}}Re\{a_nb_n^* \}}\right]$$

#### Q ratio
This is just $Q_b \over Q_{sca}$.

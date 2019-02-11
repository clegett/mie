function result = Mie_tetado(x,fact)

% Computation of Mie Power Scattering and diffraction functions 
% and gi coefficients of Legendre Polynomial decomposition
% for complex refractive-index ratio m=m'+im", size parameters x=k0*a, 
% according to Bohren and Huffman (1983) BEWI:TDD122
% 1) polar diagram, linear or in dB scale with respect to minimum, with
% SL in upper semicircle, SR in lower semicircle and 3 cartesian diagrams
% 2) same for SL0 and SR0 without diffraction pattern, 
% 3) scattered intensity S (lin or log scale), and degree of polarisation
% 4) scattered intensity without diffraction peak S0 (lin or log scale),
% 5) beam efficiencies of S and S0, diffraction efficiency Qd
% 6) gi-factors (coefficients of Legendre Polynomials of Phase Function).
% C. Mätzler, April 2004.

m1=[1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5];
nj=length(m1);
Q=[];gi=[];g0i=[];
m=m1+0.001i;
for j=1:nj,
    y=mie_teta(m(j),x,fact);
    y.Q=[m1(j),y.Q];
    Q=[Q;y.Q];
    gi=[gi;y.gi];
    g0i=[g0i;y.g0i];
end;

result.Q=Q;
result.gi=gi;
result.g0i=g0i;
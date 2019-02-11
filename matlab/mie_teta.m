function result = Mie_teta(m, x, tetadlim)

% Computation of Mie Efficienicies and gi coefficients of Legendre 
% Polynomial decomposition for complex refractive-index ratio m=m'+im",
% size parameters x=k0*a, scattered intensity with and without 
% diffraction peak. 
% tetadlim (deg): Optional value of integration limit 
% for beam and diffraction efficiencies (default 180°)
% Output.etab: etab(180°), etab0(180°), etab(tetadlim)
% Output.Q: Qext, Qsca, Qabs, Qb, Qd, Qdlim, asy
% Output.gi and Output.g0i Legendre Coefficients
% C. Mätzler, April 2004.

nsteps=round(23*x);     % number of angular steps
if nargin==2,
    tetadlim=180;       % default value 180°
end;
nx=(1:nsteps); dteta=pi/nsteps;
Q=mie(m,x);  Qext=Q(1); Qsca=Q(2); Qabs=Q(3); Qb=Q(4); asy=Q(5);
nmax=round(2+x+4*x^(1/3));
ab=mie_ab(m,x); an=ab(1,:); bn=ab(2,:);
teta=(nx-0.5).*dteta; tetad=teta*180/pi;
u=cos(teta); s=sin(teta);
px=pi*x^2;
st=pi*s*dteta/Qsca;      % Constant factor of angular integrands
for j = 1:nsteps, 
    pt=mie_pt(u(j),nmax);
    pin =pt(1,:);
    tin =pt(2,:);
    n=(1:nmax);
    n2=(2*n+1)./(n.*(n+1));
    pin=n2.*pin; tin=n2.*tin;
    S1=(an*pin'+bn*tin'); 
    S2=(an*tin'+bn*pin');
    xs=x.*s(j);
    if abs(xs)<0.00002   % Diffraction pattern according to BH, p. 110
        S3=x.*x*0.25.*(1+u(j));    % avoiding division by zero
    else
        S3=x.*x*0.5.*(1+u(j)).*besselj(1,xs)./xs;    
    end;
    S4=S1-S3;
    S5=S2-S3;
    SR(j)= real(S1'*S1)/px;
    SL(j)= real(S2'*S2)/px;
    SR0(j)=real(S4'*S4)/px;        
    SL0(j)=real(S5'*S5)/px;        
end;
z=st.*(SL+SR);       % Integrand 
z0=st.*(SL0+SR0);    % Integrand 
etabdlim=cumsum(z);  % Beam Efficiency with limited angle
nj=11;               
for jj=1:nj,         % Phase fct. decomposition in Legendre Polynomials
    xa=legendre(jj-1,u);  % Legendre Function
    x0=xa(1,:);      % Legendre Polynomial
    gi(jj)=x0*z';    % Beam Eff., asymmetry factor and higher gi's 
    g0i(jj)=x0*z0';  % same, but diffraction signal removed 
end;
etab=gi(1);     gi=gi/etab;     gi=gi(2:nj);
etab0=g0i(1);   g0i=g0i/etab0;  g0i=g0i(2:nj);
Qd=Qsca*(etab-etab0);      % Qd = diffraction efficiency
n=max(find(tetad<tetadlim)); 
Qdlim=Qsca*etabdlim(n);    % Qd = limited diffraction efficiency

result.etab=[etab, etab0, etabdlim(n)];
result.Q=[Qext, Qsca, Qabs, Qb, Qd, Qdlim, asy];
result.gi=gi;
result.g0i=g0i;
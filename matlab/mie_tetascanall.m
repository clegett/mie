function result = Mie_tetascanall(m, x, nsteps, nsmooth, type)

% Computation and plot of Mie Power Scattering and diffraction functions 
% for complex refractive-index ratio m=m'+im", size parameters x=k0*a, 
% according to Bohren and Huffman (1983) BEWI:TDD122
% 1) polar diagram, linear or in dB scale with respect to minimum, with
% SL in upper semicircle, SR in lower semicircle and 3 cartesian diagrams
% 2) same for SL0 and SR0 without diffraction pattern, 
% 3) scattered intensity S (lin or log scale), and degree of polarisation
% 4) scattered intensity without diffraction peak S0 (lin or log scale),
% 5) beam efficiencies of S and S0, diffraction efficiency Qd
% 6) gi-factors (coefficients of Legendre Polynomials of Phase Function).
% nsteps: number of scattering angles (for accurate comp. use nsteps=22*x)
% nsmooth: number of values to be averaged in polarisation, S and S0 
% (for 'log' type only)
% type:= 'log' or 'lin' for logarithmic or linear plots
% C. Mätzler, April 2004.

nstart=round(min(0.5*nsteps,nsteps*pi/x+nsmooth));
m1=real(m); m2=imag(m);
nx=(1:nsteps); dteta=pi/nsteps;
Q=mie(m,x);  Qext=Q(1); Qsca=Q(2); Qabs=Q(3); Qb=Q(4); asy=Q(5);
nmax=round(2+x+4*x^(1/3));
ab=mie_ab(m,x); an=ab(1,:); bn=ab(2,:);
teta=(nx-0.5).*dteta; tetad=teta*180/pi;
u=cos(teta); s=sin(teta);
px=pi*x^2;
st=pi*s*dteta/Qsca;
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
    if abs(xs)<0.00002         % Diffraction pattern according to BH, p. 110
        S3=x.*x*0.25.*(1+u(j));    % avoiding division by zero
    else
        S3=x.*x*0.5.*(1+u(j)).*besselj(1,xs)./xs;    
    end;
    S4=S1-S3;
    S5=S2-S3;
    SR(j)= real(S1'*S1)/px;
    SL(j)= real(S2'*S2)/px;
    SD(j)= real(S3'*S3)/px;        
    SR0(j)=real(S4'*S4)/px;        
    SL0(j)=real(S5'*S5)/px;        
end;
z=st.*(SL+SR);
z0=st.*(SL0+SR0);

nj=11;             % Phase fct decomposition in Legendre Polynomials
for jj=1:nj,
    xa=legendre(jj-1,u);
    x0=xa(1,:);
    gi(jj)=x0*z';    % beam efficiency, asymm. factor and higher gi's 
    g0i(jj)=x0*z0';  % same as gi's, but diffraction signal removed 
end;
etab=gi(1);     gi=gi/etab;     gi=gi(2:nj);
etab0=g0i(1);   g0i=g0i/etab0;  g0i=g0i(2:nj);
Qd=Qsca*(1-etab0);           % Qd = diffraction efficiency
z=cumsum(z);                 % Beam Efficiency vs. tetalim
z0=cumsum(z0);

S=(SL+SR); S0=(SL0+SR0);     % Intensity
Ss=smooth(S,nsmooth); S0s=smooth(S0,nsmooth);
dS=(SR-SL)./S;          % Degree of polarisation
dSs=smooth(dS,nsmooth);

figure;
if type=='lin'                   % linear plots
    y=[teta teta+pi;SR SL(nsteps:-1:1)]'; 
    polar(y(:,1),y(:,2)),
    title(sprintf('Mie Scattering Diagram: m=%g+%gi, x=%g',m1,m2,x)),
    xlabel('Scattering Angle'),
figure;
subplot(2,1,1);
    plot(tetad,S,'k-')
    title(sprintf('Mie Angular Scattering: m=%g+%gi, x=%g',m1,m2,x)),
    xlabel('Scattering Angle'),
    ylabel('S');
subplot(2,1,2);
    plot(tetad,dS,'k-')
    xlabel('Scattering Angle'),
    ylabel('Polarisation Degree ');
elseif type=='log',              % logar. plots
    y=[teta teta+pi;10*log10(SR) 10*log10(SL(nsteps:-1:1))]';
    ymin=min(y(:,2));            % Minimum for normalisation of log-polar plot
    y(:,2)=y(:,2)-ymin;                  
    polar(y(:,1),y(:,2)),
    title(sprintf('Mie Scattering Diagram: m=%g+%gi, x=%g, min(dB)= %g',m1,m2,x,ymin)),
    xlabel('Scattering Angle (deg)');
    y=[teta teta+pi;10*log10(SR0) 10*log10(SL0(nsteps:-1:1))]';
    ymin=min(y(:,2));            % Minimum for normalisation of log-polar plot
    y(:,2)=y(:,2)-ymin;                  
figure;
    polar(y(:,1),y(:,2)),
    title(sprintf('No-Peak Scattering Diagram: m=%g+%gi, x=%g, min(dB)= %g',m1,m2,x,ymin)),
    xlabel('Scattering Angle (deg)');
figure;
subplot(2,1,1);
    semilogy(tetad,Ss,'k-'),    
    title(sprintf('Mie Angular Scattering: m=%g+%gi, x=%g',m1,m2,x)),
    xlabel('Scattering Angle'),
    ylabel('S');
subplot(2,1,2);
    plot(tetad,dSs,'k-'),
    xlabel('Scattering Angle'),
    ylabel('Polarisation Degree ');
figure;
    semilogy(tetad(nstart:nsteps),Ss(nstart:nsteps),'r:',tetad,S0s,'k-'),    
    title(sprintf('No-Peak Angular Scattering: m=%g+%gi, x=%g',m1,m2,x)),
    xlabel('Scattering Angle'),
    ylabel('S0');
figure;
    xmin=min(tetad/180);
    semilogx(tetad/180,z,'r-',tetad/180,z0,'k--',tetad/180,z-z0,'b:'),
    xlabel('Maximum Scattering Angle/180°'),
    axis([xmin, 1, 0, 1.1]);
end;
result.s=[tetad',SR',SL',SR0',SL0',SD',dS',z',z0']; 
result.Q=[Qext,Qsca,Qabs,Qb,Qd,asy];
result.gi=gi;
result.g0i=g0i;
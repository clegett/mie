function result = mie_tetascan(m, x, nsteps, type)

% Computation and plot of Mie Power Scattering function for given 
% complex refractive-index ratio m=m'+im", size parameters x=k0*a, 
% according to Bohren and Huffman (1983) BEWI:TDD122
% type ='pol' for polar diagram, else cartesian for intensity (0)
% polarisation (1) of both modes ('b')
% C. Mätzler, May 2002, revised March 2004.

nsteps=nsteps;
m1=real(m); m2=imag(m);
nx=(1:nsteps); dteta=pi/(nsteps-1);
teta=(nx-1).*dteta;
    for j = 1:nsteps, 
        u=cos(teta(j));
        a(:,j)=mie_S12(m,x,u);
        SL(j)= real(a(1,j)'*a(1,j))/(pi*x^2);
        SR(j)= real(a(2,j)'*a(2,j))/(pi*x^2);
    end;
S=0.5*(SL+SR);                       % Intenisty
dS=0.5*(SL-SR)./S;                   % Degree of lienar polarisation
dS10=smooth(dS,10);
y=[teta teta+pi;SL SR(nsteps:-1:1)]'; 
tetad=teta*180/pi;
if type=='pol'
    polar(y(:,1),y(:,2))
    title(sprintf('Mie angular scattering: m=%g+%gi, x=%g',m1,m2,x));
    xlabel('Scattering Angle')
elseif type==0,
    semilogy(tetad,S,'k-')
    title(sprintf('Mie angular scattering: m=%g+%gi, x=%g',m1,m2,x))
    xlabel('Scattering Angle')
    ylabel('Intensity')
elseif type==1,
    plot(tetad,dS,'k-')
    title(sprintf('Mie angular scattering: m=%g+%gi, x=%g',m1,m2,x))
    xlabel('Scattering Angle'),
    ylabel('Degree of Linear Polarisation')
else
    semilogy(tetad,SR,'r-',tetad,SL,'b--')
    title(sprintf('Mie angular scattering: m=%g+%gi, x=%g',m1,m2,x))
    xlabel('Scattering Angle')
    legend('SR','SL')
end;
result=[tetad',cos(tetad)',SR',SL',S',dS']; 

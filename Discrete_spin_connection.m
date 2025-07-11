close all
clear all
clc


%This file provides the code for finding site coordinates, radius of the
%corresponding geodesic curve and the discrete spin connection in the 0th
%and 1st generation of arbitrary tessellation of Poincare disk

%Define the Schlaffli symbols {p,q} and angles corresponding to the regular
%tessellation of the Poincare disk
p=10;
q=3;
alpha=pi/p;
beta=pi/q;

%r0=AB=AD, see Fig. 6 (radius of 0th generation sites to the center of the
%Poincare disk)
r0=sqrt(cos(alpha+beta)/cos(beta-alpha));
%kappa=AC, see Fig. 6
kappa=(cos(beta)-sin(alpha))/sqrt(cos(beta+alpha)*cos(beta-alpha));
%OC=OF=OB=R, see Fig. 6 (radius of the geodesic curve in the 0th
%generation)
R=r0*sin(alpha)/sin(pi/2-alpha-beta);
%Euclidean distance d_E between 0th and 1st generation polygon centers, see Eq. (94)
x=sqrt(1-sin(alpha)^2/cos(beta)^2);

%Define edge of the Poincare disk
c=linspace(-1,1,1000);
%Define auxiliary matrices for site coordinates
phi=zeros(p+1,p); %p+1 for proper defining the 0th generation
r=zeros(p+1,1);
rx=zeros(p+1,p);
ry=zeros(p+1,p);
%Matrices for centers of geodesic curves (xc,yc), type from Eq. (21)
xc=zeros(p+1,p);
yc=zeros(p+1,p);
chi=zeros(p+1,p);
%Auxiliary coordinate when the geodesic curve changes from one site to the
%other
x1=zeros(p);
y1=zeros(p);
r1=zeros(p,1);
%Radius of geodesic curve, type Eq. (21)
Rc=zeros(p,1);
%Auxiliary variables for theta=+/-1/2(mu-nu)+theta0, Eq.(F12)
muj=zeros(p,1);
nuj=zeros(p,1);
mui=zeros(p,1);
nui=zeros(p,1);
mui1=zeros(p,1);
nui1=zeros(p,1);
%Define the integration constant theta0 and a
thetaj=zeros(p,1);
thetai=zeros(p,1);
thetai1=zeros(p,1);
theta0=zeros(p);
theta1=zeros(p);
a1=ones(p+1,p); 
%Parallel transport vector
V1=zeros(2,p+1,p);
V2=zeros(p,1);
%Reference vector
Vr=zeros(2,p,p);
%Tangent vectors to geodesic curves
Vgamma1=zeros(2,p,p); %on the first site on the geodesic
Vgamma=zeros(2,p,p); %on the second site on the geodesic
%beta_j and alpha_i
alpha_g=zeros(p); %pure geometrical
alpha_p=zeros(p); %parallel transport
beta_p=zeros(p); %parallel transport
%theta_{ij} needed to construct the discrete spin connection matrix
thetaij_g=zeros(p); %pure geometrical
thetaij_p=zeros(p); %parallel transport
%Discrete spin connection matrix
Omega_g=zeros(2,2,p); %pure geometrical
Omega_p=zeros(2,2,p); %parallel transport
%Variable to check Eq.(121)
P1=eye(2);
%Gamma matrices
gamma_x=[0 1; 1 0];
gamma_y=[0 -1i; 1i 0];
%Auxiliary to check the tetrad hypothesis
l=zeros(2,2,2*p);

%Set initial values
Rc(p)=R;

%Plot the edge of the Poincare disk
figure (1)
plot(c,sqrt(1-c.^2),'k','Linewidth',1.5);
axis equal;
hold all;
ylim([-1 1]);
plot(c,-sqrt(1-c.^2),'k','Linewidth',1.5);

for m=1:1:p
    for j=1:1:p
        ej=exp(2i*j*alpha);
        em=exp(1i*(2*m-1)*alpha);
        em1=exp(2i*(m-1)*alpha);
        e_m=exp(-1i*(2*m-1)*alpha);
        e_m1=exp(-2i*(m-1)*alpha);
        e_j=exp(-2i*j*alpha);
        %Find coordinates of every site in 0th and 1st generation, see Eq.
        %(83)
        %1st direction
        z=(ej-x^2)/(1-x^2*ej)*(r0*em-x*em1*(ej-1)/(ej-x^2))/(1-x/em1*r0*em*(1-ej)/(1-x^2*ej));
        rx(j,m)=real(z);
        ry(j,m)=imag(z);
        if ry(j,m)>0
            phi(j,m)=angle(z);
        else
            phi(j,m)=2*pi+angle(z);
        end
        r(j)=abs(z);
        %2nd direction
        z1=(e_j-x^2)/(1-x^2*e_j)*(r0*e_m-x*e_m1*(e_j-1)/(e_j-x^2))/(1-x/e_m1*r0*e_m*(1-e_j)/(1-x^2*e_j));
        %Set initial conditions
        yc(p,m)=(kappa+R)*imag(em1);
        xc(p,m)=(kappa+R)*real(em1);
        %Plot sites of the 0th and 1st generation
        figure(1)
        scatter(real(z),imag(z),'filled','b');
        %scatter(xc(p,m),yc(p,m),'filled','r');
    end
    rx(p+1,m)=rx(1,m);
    ry(p+1,m)=ry(1,m);
    phi(p+1,m)=phi(1,m);
end
r(p+1)=r(1);

for m=1:1:p
    for j=1:1:p
        if j==1
            %Reference vector
            Vr(1,p,m)=(rx(j,m)-rx(p,m))/sqrt((rx(j,m)-rx(p,m))^2+(ry(j,m)-ry(p,m))^2);
            Vr(2,p,m)=(ry(j,m)-ry(p,m))/sqrt((rx(j,m)-rx(p,m))^2+(ry(j,m)-ry(p,m))^2);
            %Tangent vectors
            gamma=atan((rx(j,m)-xc(p,m))/(yc(p,m)-ry(j,m)));
            gamma1=atan((rx(p,m)-xc(p,m))/(yc(p,m)-ry(p,m)));
            %Auxiliary variables to fix the tangent vector direction
            tau1=(xc(p,m)-rx(j,m))/sqrt((xc(p,m)-rx(j,m))^2+(yc(p,m)-ry(j,m))^2);
            tau2=(yc(p,m)-ry(j,m))/sqrt((xc(p,m)-rx(j,m))^2+(yc(p,m)-ry(j,m))^2);
            tau3=(xc(p,m)-rx(p,m))/sqrt((xc(p,m)-rx(p,m))^2+(yc(p,m)-ry(p,m))^2);
            tau4=(yc(p,m)-ry(p,m))/sqrt((xc(p,m)-rx(p,m))^2+(yc(p,m)-ry(p,m))^2);
            if tau1*sin(gamma)-tau2*cos(gamma)<0
                Vgamma(1,j,m)=cos(gamma);
                Vgamma(2,j,m)=sin(gamma);
            else
                Vgamma(1,j,m)=-cos(gamma);
                Vgamma(2,j,m)=-sin(gamma);
            end
            if tau3*sin(gamma1)-tau4*cos(gamma1)<0
                Vgamma1(1,p,m)=cos(gamma1);
                Vgamma1(2,p,m)=sin(gamma1);
            else
                Vgamma1(1,p,m)=-cos(gamma1);
                Vgamma1(2,p,m)=-sin(gamma1);
            end
        end
        if j>=2
            %Reference vector
            Vr(1,j-1,m)=(rx(j,m)-rx(j-1,m))/sqrt((rx(j,m)-rx(j-1,m))^2+(ry(j,m)-ry(j-1,m))^2);
            Vr(2,j-1,m)=(ry(j,m)-ry(j-1,m))/sqrt((rx(j,m)-rx(j-1,m))^2+(ry(j,m)-ry(j-1,m))^2);
            %Geodesic curve between two sites type Eq. (20), 1st direction
            if vpa(phi(j-1,m))==vpa(phi(j,m))
                Vgamma(1,j,m)=(rx(j,m)-rx(j-1,m))/sqrt((rx(j,m)-rx(j-1,m))^2+(ry(j,m)-ry(j-1,m))^2);
                Vgamma(2,j,m)=(ry(j,m)-ry(j-1,m))/sqrt((rx(j,m)-rx(j-1,m))^2+(ry(j,m)-ry(j-1,m))^2);
                Vgamma1(1,j-1,m)=(rx(j,m)-rx(j-1,m))/sqrt((rx(j,m)-rx(j-1,m))^2+(ry(j,m)-ry(j-1,m))^2);
                Vgamma1(2,j-1,m)=(ry(j,m)-ry(j-1,m))/sqrt((rx(j,m)-rx(j-1,m))^2+(ry(j,m)-ry(j-1,m))^2);
            %Geodesic curve between two sites type Eq. (21), 1st direction
            else
                syms u v
                %Find the center of the orthogonal circle between points
                %z(j,m) and z(j-1,m)
                [xc(j-1,m),yc(j-1,m)]=vpasolve([(rx(j,m)-u).^2+(ry(j,m)-v).^2==u.^2+v.^2-1,(rx(j-1,m)-u).^2+(ry(j-1,m)-v).^2==u.^2+v.^2-1],[u,v]);
                figure(1)
                %scatter(xc(j-1,m),yc(j-1,m),'filled','r');
                Rc(j-1)=sqrt(xc(j-1,m)^2+yc(j-1,m)^2-1);
                %Determine the angle of the geodesic curve center
                if yc(j-1,m)>=0
                    chi(j-1,m)=angle(xc(j-1,m)+1i*yc(j-1,m));
                else
                    chi(j-1,m)=2*pi+angle(xc(j-1,m)+1i*yc(j-1,m));
                end
                %Find theta for parallel transport, Eq.(F12)
                muj(j-1)=asin((-(2*Rc(j-1)^2+1)+r(j)^2)/(2*Rc(j-1)*sqrt(Rc(j-1)^2+1)));
                nuj(j-1)=asin((-(2*Rc(j-1)^2+1)+1/r(j)^2)/(2*Rc(j-1)*sqrt(Rc(j-1)^2+1)));
                thetaj(j-1)=1/2*(muj(j-1)-nuj(j-1));
                mui(j-1)=asin((-(2*Rc(j-1)^2+1)+r(j-1)^2)/(2*Rc(j-1)*sqrt(Rc(j-1)^2+1)));
                nui(j-1)=asin((-(2*Rc(j-1)^2+1)+1/r(j-1)^2)/(2*Rc(j-1)*sqrt(Rc(j-1)^2+1)));
                thetai(j-1)=1/2*(mui(j-1)-nui(j-1));
                %Tangent vector to the circle in point z(j,m)
                gamma=atan((rx(j,m)-xc(j-1,m))/(yc(j-1,m)-ry(j,m)));
                %Tangent vector to the circle in point z(j-1,m)
                gamma1=atan((rx(j-1,m)-xc(j-1,m))/(yc(j-1,m)-ry(j-1,m)));
                %Auxiliary vectors
                tau1=(xc(j-1,m)-rx(j,m))/sqrt((xc(j-1,m)-rx(j,m))^2+(yc(j-1,m)-ry(j,m))^2);
                tau2=(yc(j-1,m)-ry(j,m))/sqrt((xc(j-1,m)-rx(j,m))^2+(yc(j-1,m)-ry(j,m))^2);
                tau3=(xc(j-1,m)-rx(j-1,m))/sqrt((xc(j-1,m)-rx(j-1,m))^2+(yc(j-1,m)-ry(j-1,m))^2);
                tau4=(yc(j-1,m)-ry(j-1,m))/sqrt((xc(j-1,m)-rx(j-1,m))^2+(yc(j-1,m)-ry(j-1,m))^2);
                if tau1*sin(gamma)-tau2*cos(gamma)>0
                    Vgamma(1,j,m)=cos(gamma);
                    Vgamma(2,j,m)=sin(gamma);
                else
                    Vgamma(1,j,m)=-cos(gamma);
                    Vgamma(2,j,m)=-sin(gamma);
                end
                if tau3*sin(gamma1)-tau4*cos(gamma1)>0
                    Vgamma1(1,j-1,m)=cos(gamma1);
                    Vgamma1(2,j-1,m)=sin(gamma1);
                else
                    Vgamma1(1,j-1,m)=-cos(gamma1);
                    Vgamma1(2,j-1,m)=-sin(gamma1);
                end
            end
        end
    end
end

for m=1:1:p
    %Determine the angle of the geodesic curve center for 0th generation
    if yc(p,m)>=0
        chi(p,m)=angle(xc(p,m)+1i*yc(p,m));
    else
        chi(p,m)=2*pi+angle(xc(p,m)+1i*yc(p,m));
    end
end

%Find theta for parallel transport, Eq.(F12) for 0th generation
muj(p)=asin((-(2*Rc(p)^2+1)+r(p+1)^2)/(2*Rc(p)*sqrt(Rc(p)^2+1)));
nuj(p)=asin((-(2*Rc(p)^2+1)+1/r(p+1)^2)/(2*Rc(p)*sqrt(Rc(p)^2+1)));
thetaj(p)=1/2*(muj(p)-nuj(p));
mui(p)=asin(vpa((-(2*Rc(p)^2+1)+r(p)^2)/(2*Rc(p)*sqrt(Rc(p)^2+1))));
nui(p)=asin(vpa((-(2*Rc(p)^2+1)+1/r(p)^2)/(2*Rc(p)*sqrt(Rc(p)^2+1))));
thetai(p)=1/2*(mui(p)-nui(p));
                
%Set initial conditions for parallely transported vectors
V1(1,1,1)=a1(1,1)*(1-r(1)^2)*Vgamma1(1,1,1);
V1(2,1,1)=a1(1,1)*(1-r(1)^2)*Vgamma1(2,1,1);

%Define the parallely transported vectors
for m=1:1:p
    for j=1:1:p+1
        if j==1 && m>=2
            V1(1,j,m)=V1(1,p,m-1);
            V1(2,j,m)=V1(2,p,m-1);
        end
        if j>=2
            %If the geodesic curve between sites is type Eq. (20)
            if vpa(phi(j-1,m))==vpa(phi(j,m))
                V1(1,j,m)=(a1(j,m)*(1-r(j)^2))/(a1(j-1,m)*(1-r(j-1)^2))*V1(1,j-1,m);
                V1(2,j,m)=(a1(j,m)*(1-r(j)^2))/(a1(j-1,m)*(1-r(j-1)^2))*V1(2,j-1,m);
            %If the geodesic curve between sites is type Eq. (21)
            else
                %Checking if the sign of parallel transport changes in Eq.(F10)
                if cot(vpa(chi(j-1,m)-phi(j-1,m)))*cot(vpa(chi(j-1,m)-phi(j,m)))<0
                    syms u v 'real'
                    assumeAlso(u>-1 & u<1)
                    assumeAlso(v>-1 & v<1)
                    u0=rx(j-1,m);
                    v0=ry(j-1,m);
                    sol=vpasolve([(yc(j-1,m)-v)./(xc(j-1,m)-u).*u-v==0,(u-xc(j-1,m)).^2+(v-yc(j-1,m)).^2==xc(j-1,m).^2+yc(j-1,m).^2-1],[u,v],[u0,v0]);
                    if ~isempty(sol) && isfield(sol,'u') && isfield(sol,'v')
                        if numel(sol.u)>1
                            distances=double((sol.u-rx(j-1,m)).^2+(sol.v-ry(j-1,m)).^2);
                            [~,idx]=min(distances);
                            %Defining the coordinates where the switch
                            %happens
                            x1(j-1,m)=double(sol.u(idx));
                            y1(j-1,m)=double(sol.v(idx));
                        else
                            x1(j-1,m)=double(sol.u);
                            y1(j-1,m)=double(sol.v);
                        end
                    else
                        warning('No valid solution found at j=%d,m=%d',j-1,m);
                        x1(j-1,m)=NaN;
                        y1(j-1,m)=NaN;
                    end
                    r1(j-1,m)=sqrt(x1(j-1,m)^2+y1(j-1,m)^2);
                    if cot(vpa(chi(j-1,m)-phi(j-1,m)))<0
                        %Integration constant theta0
                        theta0(j-1,m)=atan(V1(2,j-1,m)/V1(1,j-1,m))+pi/2+thetai(j-1)-phi(j-1,m);
                        %Auxiliary variables for finding theta, Eq.(F12)
                        mui1(j-1)=asin((-(2*Rc(j-1)^2+1)+r1(j-1)^2)/(2*Rc(j-1)*sqrt(Rc(j-1)^2+1)));
                        nui1(j-1)=asin(vpa((-(2*Rc(j-1)^2+1)+1/r1(j-1)^2)/(2*Rc(j-1)*sqrt(Rc(j-1)^2+1))));
                        thetai1(j-1)=1/2*(mui1(j-1)-nui1(j-1));
                        %Integration constant a, Eq.(F8)
                        a11=(V1(1,j-1,m)*cos(-thetai1(j-1)+pi/2+theta0(j-1,m)+phi(j-1,m))+V1(2,j-1,m)*sin(-thetai1(j-1)+pi/2+theta0(j-1,m)+phi(j-1,m)))/(1-r(j)^2);
                        V(1)=(1-r1(j-1,m)^2)*a11*cos(-thetai1(j-1)+pi/2+theta0(j-1,m)+chi(j-1,m));
                        V(2)=(1-r1(j-1,m)^2)*a11*sin(-thetai1(j-1)+pi/2+theta0(j-1,m)+chi(j-1,m));
                        %Change of a sign in parallel transport, Eq.(F10)
                        theta1(j-1,m)=atan(V(2)/V(1))-pi/2-thetai1(j-1)-chi(j-1,m);
                        a1(j,m)=(V(1)*cos(pi/2+thetaj(j-1)+theta1(j-1,m)+chi(j-1,m))+V(2)*sin(pi/2+thetaj(j-1)+theta1(j-1,m)+chi(j-1,m)))/(1-r1(j-1)^2);
                        V1(1,j,m)=(1-r(j)^2)*a1(j,m)*cos(pi/2+thetaj(j-1)+theta1(j-1,m)+phi(j,m));
                        V1(2,j,m)=(1-r(j)^2)*a1(j,m)*sin(pi/2+thetaj(j-1)+theta1(j-1,m)+phi(j,m));
                    else
                        theta0(j-1,m)=atan(vpa(V1(2,j-1,m)/V1(1,j-1,m)))+pi/2-thetai(j-1)-phi(j-1,m);
                        mui1(j-1)=asin(vpa((-(2*Rc(j-1)^2+1)+r1(j-1)^2)/(2*Rc(j-1)*sqrt(Rc(j-1)^2+1))));
                        nui1(j-1)=asin(vpa((-(2*Rc(j-1)^2+1)+1/r1(j-1)^2)/(2*Rc(j-1)*sqrt(Rc(j-1)^2+1))));
                        thetai1(j-1)=1/2*(mui1(j-1)-nui1(j-1));
                        a11=(V1(1,j-1,m)*cos(pi/2+thetai1(j-1)+theta0(j-1,m)+phi(j-1,m))+V1(2,j-1,m)*sin(pi/2+thetai1(j-1)+theta0(j-1,m)+phi(j-1,m)))/(1-r(j)^2);
                        V(1)=(1-r1(j-1,m)^2)*a11*cos(pi/2+thetai1(j-1)+theta0(j-1,m)+chi(j-1,m));
                        V(2)=(1-r1(j-1,m)^2)*a11*sin(pi/2+thetai1(j-1)+theta0(j-1,m)+chi(j-1,m));
                        %Change of a sign in parallel transport, Eq.(F10)
                        theta1(j-1,m)=atan(vpa(V(2)/V(1)))+pi/2+thetai1(j-1)-chi(j-1,m);
                        a1(j,m)=(V(1)*cos(pi/2-thetaj(j-1)+theta1(j-1,m)+chi(j-1,m))+V(2)*sin(pi/2-thetaj(j-1)+theta1(j-1,m)+chi(j-1,m)))/(1-r1(j-1)^2);
                        V1(1,j,m)=(1-r(j)^2)*a1(j,m)*cos(pi/2-thetaj(j-1)+theta1(j-1,m)+phi(j,m));
                        V1(2,j,m)=(1-r(j)^2)*a1(j,m)*sin(pi/2-thetaj(j-1)+theta1(j-1,m)+phi(j,m));
                    end
                else
                    %If the sign doesn't switch in parallel transport,
                    %Eq.(F10), check which sign is appropriate
                    if cot(vpa(chi(j-1,m)-phi(j,m)))<0
                        theta0(j-1,m)=atan(vpa(V1(2,j-1,m)/V1(1,j-1,m)))+pi/2+thetai(j-1)-phi(j-1,m);
                        a1(j,m)=(V1(1,j-1,m)*cos(pi/2-thetaj(j-1)+theta0(j-1,m)+phi(j-1,m))+V1(2,j-1,m)*sin(pi/2-thetaj(j-1)+theta0(j-1,m)+phi(j-1,m)))/(1-r(j)^2);
                        V1(1,j,m)=(1-r(j)^2)*a1(j,m)*cos(pi/2-thetaj(j-1)+theta0(j-1,m)+phi(j,m));
                        V1(2,j,m)=(1-r(j)^2)*a1(j,m)*sin(pi/2-thetaj(j-1)+theta0(j-1,m)+phi(j,m));
                    else
                        theta0(j-1,m)=atan(vpa(V1(2,j-1,m)/V1(1,j-1,m)))+pi/2-thetai(j-1)-phi(j-1,m);
                        a1(j,m)=(V1(1,j-1,m)*cos(pi/2+thetaj(j-1)+theta0(j-1,m)+phi(j-1,m))+V1(2,j-1,m)*sin(pi/2+thetaj(j-1)+theta0(j-1,m)+phi(j-1,m)))/(1-r(j)^2);
                        V1(1,j,m)=(1-r(j)^2)*a1(j,m)*cos(pi/2+thetaj(j-1)+theta0(j-1,m)+phi(j,m));
                        V1(2,j,m)=(1-r(j)^2)*a1(j,m)*sin(pi/2+thetaj(j-1)+theta0(j-1,m)+phi(j,m));
                    end
                end
            end
        end
    end
end

%Check on figure if everything is calculated okay
%for j=1:1:p+1
%    tr=linspace(-1,1,1000);
%    figure(1)
%    plot(tr,V1(2,j,1)/V1(1,j,1)*(tr-rx(j,1))+ry(j,1),'r');
%end

for m=1:1:p
    for j=1:1:p
        %Check on figure if everything is calculated okay
        %tr=linspace(-1,1,1000);
        %figure(1)
        %plot(tr,V1(2,j,1)/V1(1,j,1)*(tr-rx(j,1))+ry(j,1),'r');
        %plot(tr,Vgamma1(2,j,1)/Vgamma1(1,j,1)*(tr-rx(j,1))+ry(j,1),'g');
        %plot(tr,Vgamma(2,j,1)/Vgamma(1,j,1)*(tr-rx(j,1))+ry(j,1),'m');
        %Define angle alpha_i^g, geometrical approach
        if vpa(Vgamma1(1,j,m)*Vr(2,j,m))>=vpa(Vgamma1(2,j,m)*Vr(1,j,m))
            alpha_g(j,m)=acos(vpa((Vgamma1(1,j,m)*Vr(1,j,m)+Vgamma1(2,j,m)*Vr(2,j,m))));
        else
            alpha_g(j,m)=-acos(vpa((Vgamma1(1,j,m)*Vr(1,j,m)+Vgamma1(2,j,m)*Vr(2,j,m))));
        end
        thetaij_g(j,m)=2*vpa(alpha_g(j,m));
        if j==1
            if vpa(V1(1,p+j,m)*Vr(2,p,m))>=vpa(V1(2,p+j,m)*Vr(1,p,m)) && vpa(V1(1,p,m)*Vr(2,p,m))>=vpa(V1(2,p,m)*Vr(1,p,m))
                beta_p(p,m)=vpa(acos(vpa((V1(1,p+j,m)*Vr(1,p,m)+V1(2,p+j,m)*Vr(2,p,m))/(abs(a1(p+j,m))*(1-r(p+j)^2)))));
                alpha_p(p,m)=vpa(acos(vpa((V1(1,p,m)*Vr(1,p,m)+V1(2,p,m)*Vr(2,p,m))/(abs(a1(p,m))*(1-r(p)^2)))));
            else
                beta_p(p,m)=vpa(2*pi-acos(vpa((V1(1,p+j,m)*Vr(1,p,m)+V1(2,p+j,m)*Vr(2,p,m))/(abs(a1(p+j,m))*(1-r(p+j)^2)))));
                alpha_p(p,m)=vpa(2*pi-acos(vpa((V1(1,p,m)*Vr(1,p,m)+V1(2,p,m)*Vr(2,p,m))/(abs(a1(p,m))*(1-r(p)^2)))));
                if vpa(V1(1,p+j,m)*Vr(2,p,m))>=vpa(V1(2,p+j,m)*Vr(1,p,m))
                    beta_p(p,m)=vpa(acos(vpa((V1(1,p+j,m)*Vr(1,p,m)+V1(2,p+j,m)*Vr(2,p,m))/(abs(a1(p+j,m))*(1-r(p+j)^2)))));
                    alpha_p(p,m)=vpa(2*pi-acos(vpa((V1(1,p,m)*Vr(1,p,m)+V1(2,p,m)*Vr(2,p,m))/(abs(a1(p,m))*(1-r(p)^2)))));
                end
                if vpa(V1(1,p,m)*Vr(2,p,m))>=vpa(V1(2,p,m)*Vr(1,p,m))
                    beta_p(p,m)=vpa(2*pi-acos(vpa((V1(1,p+j,m)*Vr(1,p,m)+V1(2,p+j,m)*Vr(2,p,m))/(abs(a1(p+j,m))*(1-r(p+j)^2)))));
                    alpha_p(p,m)=vpa(acos(vpa((V1(1,p,m)*Vr(1,p,m)+V1(2,p,m)*Vr(2,p,m))/(abs(a1(p,m))*(1-r(p)^2)))));
                end
            end
        end
        %Define angle alpha_i and beta_j from parallel transport approach
        if j>=2
            if vpa(V1(1,j,m)*Vr(2,j-1,m))>=vpa(V1(2,j,m)*Vr(1,j-1,m)) && vpa(V1(1,j-1,m)*Vr(2,j-1,m))>=vpa(V1(2,j-1,m)*Vr(1,j-1,m))
                beta_p(j-1,m)=vpa(acos(vpa((V1(1,j,m)*Vr(1,j-1,m)+V1(2,j,m)*Vr(2,j-1,m))/(abs(a1(j,m))*(1-r(j)^2)))));
                alpha_p(j-1,m)=vpa(acos(vpa((V1(1,j-1,m)*Vr(1,j-1,m)+V1(2,j-1,m)*Vr(2,j-1,m))/(abs(a1(j-1,m))*(1-r(j-1)^2)))));
            else
                beta_p(j-1,m)=vpa(2*pi-acos(vpa((V1(1,j,m)*Vr(1,j-1,m)+V1(2,j,m)*Vr(2,j-1,m))/(abs(a1(j,m))*(1-r(j)^2)))));
                alpha_p(j-1,m)=vpa(2*pi-acos(vpa((V1(1,j-1,m)*Vr(1,j-1,m)+V1(2,j-1,m)*Vr(2,j-1,m))/(abs(a1(j-1,m))*(1-r(j-1)^2)))));
                if vpa(V1(1,j,m)*Vr(2,j-1,m))>=vpa(V1(2,j,m)*Vr(1,j-1,m))
                    beta_p(j-1,m)=vpa(acos(vpa((V1(1,j,m)*Vr(1,j-1,m)+V1(2,j,m)*Vr(2,j-1,m))/(abs(a1(j,m))*(1-r(j)^2)))));
                    alpha_p(j-1,m)=vpa(2*pi-acos(vpa((V1(1,j-1,m)*Vr(1,j-1,m)+V1(2,j-1,m)*Vr(2,j-1,m))/(abs(a1(j-1,m))*(1-r(j-1)^2)))));
                end
                if vpa(V1(1,j-1,m)*Vr(2,j-1,m))>=vpa(V1(2,j-1,m)*Vr(1,j-1,m))
                    beta_p(j-1,m)=vpa(2*pi-acos(vpa((V1(1,j,m)*Vr(1,j-1,m)+V1(2,j,m)*Vr(2,j-1,m))/(abs(a1(j,m))*(1-r(j)^2)))));
                    alpha_p(j-1,m)=vpa(acos(vpa((V1(1,j-1,m)*Vr(1,j-1,m)+V1(2,j-1,m)*Vr(2,j-1,m))/(abs(a1(j-1,m))*(1-r(j-1)^2)))));
                end
            end
            thetaij_p(j-1,m)=vpa(alpha_p(j-1,m)-beta_p(j-1,m));
            if thetaij_p(j-1,m)>=pi
                thetaij_p(j-1,m)=-(2*pi-thetaij_p(j-1,m));
            end
            if thetaij_p(j-1,m)<-pi
                thetaij_p(j-1,m)=-2*pi-thetaij_p(j-1,m);
            end
        end
    end
    thetaij_p(p,m)=vpa(alpha_p(p,m)-beta_p(p,m));
    if thetaij_p(p,m)>=pi
        thetaij_p(p,m)=-(2*pi-thetaij_p(p,m));
    end
    if thetaij_p(p,m)<-pi
        thetaij_p(p,m)=-2*pi-thetaij_p(p,m);
    end
end

%Check Eq.(116)
for j=1:1:p
    Omega_g(:,:,j)=[exp(1i*thetaij_g(j,1)/2) 0; 0 exp(-1i*thetaij_g(j,1)/2)];
    Omega_p(:,:,j)=[exp(1i*thetaij_p(j,1)/2) 0; 0 exp(-1i*thetaij_p(j,1)/2)];
    P1=P1*Omega_g(:,:,j);
end

%Auxiliary checking
%for m=1:1:p
    for j=1:1:p
        if j>=2
         V2(j)=Vgamma1(1,j-1,1)*Vgamma1(1,j,1)+Vgamma1(2,j-1,1)*Vgamma1(2,j,1);
        end
    end
%end

%Check Eq.(104)
for j=1:1:p
    if j==1
        l(:,:,p)=vpa(Vgamma1(1,p)*gamma_x+Vgamma1(2,p)*gamma_y);
        l(:,:,2*p)=-Omega_p(:,:,p)*(vpa(-Vgamma(1,j)*gamma_x-Vgamma(2,j)*gamma_y))*ctranspose(Omega_p(:,:,p)); 
    end
    if j>=2
        l(:,:,j-1)=vpa(Vgamma1(1,j-1)*gamma_x+Vgamma1(2,j-1)*gamma_y);
        l(:,:,j+p-1)=-Omega_p(:,:,j-1)*(vpa(-Vgamma(1,j)*gamma_x-Vgamma(2,j)*gamma_y))*ctranspose(Omega_p(:,:,j-1));
    end
end
%Plot the geodesics curve between every two points of 0th and 1st
%generation
for j=1:1:p
    for m=1:1:p
        if j>=2 && vpa(phi(j-1,m))==vpa(phi(j,m))
            if phi(j-1,m)==0 || phi(j-1,m)==pi
                t=linspace(rx(j-1,m),rx(j,m),1000);
                x11=(ry(j,m)-ry(j-1,m))/(rx(j,m)-rx(j-1,m))*t;
                figure(1)
                plot(t,x11,'k','LineWidth',1.5);
            else
                t=linspace(ry(j-1,m),ry(j,m),1000);
                y11=(rx(j,m)-rx(j-1,m))/(ry(j,m)-ry(j-1,m))*t;
                figure(1)
                plot(y11,t,'k','LineWidth',1.5);
            end
        end
        if j>=2
            if rx(j-1,m)>=xc(j-1,m)
                t1=atan((ry(j-1,m)-yc(j-1,m))/(rx(j-1,m)-xc(j-1,m)));
            else
                t1=pi+atan((ry(j-1,m)-yc(j-1,m))/(rx(j-1,m)-xc(j-1,m)));    
            end
            if rx(j,m)>=xc(j-1,m)
                t2=atan((ry(j,m)-yc(j-1,m))/(rx(j,m)-xc(j-1,m)));
            else
                t2=pi+atan((ry(j,m)-yc(j-1,m))/(rx(j,m)-xc(j-1,m)));
            end
            if abs(t1-t2)<=pi
                t=linspace(t1,t2,1000);
            else
                t=linspace(t2,t1+2*pi,1000);
            end
            figure (1)
            plot(xc(j-1,m)+Rc(j-1)*cos(t),yc(j-1,m)+Rc(j-1)*sin(t),'k','LineWidth',1.5);
        else
            if rx(p,m)>=xc(p,m)
                t1=atan((ry(p,m)-yc(p,m))/(rx(p,m)-xc(p,m)));
            else
                t1=pi+atan((ry(p,m)-yc(p,m))/(rx(p,m)-xc(p,m)));    
            end
            if rx(j,m)>=xc(p,m)
                t2=atan((ry(j,m)-yc(p,m))/(rx(j,m)-xc(p,m)));
            else
                t2=pi+atan((ry(j,m)-yc(p,m))/(rx(j,m)-xc(p,m)));
            end
            if abs(t1-t2)<=pi
                t=linspace(t1,t2,1000);
            else
                t=linspace(t1,t2+2*pi,1000);
            end
            figure (1)
            plot(xc(p,m)+Rc(p)*cos(t),yc(p,m)+Rc(p)*sin(t),'k','LineWidth',1.5);
        end    
    end
end        
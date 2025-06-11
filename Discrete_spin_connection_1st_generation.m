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
%Auxiliary variable
A=1/sqrt(r0^2+(kappa+R)^2-2*r0*(kappa+R)*cos(alpha));

%Define edge of the Poincare disk
c=linspace(-1,1,1000);
%Define auxiliary matrices for site coordinates
phi=zeros(2*p,p);
r=zeros(2*p,1);
rx=zeros(2*p,p);
ry=zeros(2*p,p);
%Matrices for centers of geodesic curves (xc,yc), type from Eq. (21)
xc=zeros(2*p,p);
yc=zeros(2*p,p);
%Radius of geodesic curve, type from Eq. (21)
Rc=zeros(p,1);
%Define the integration constant theta0 and a
theta0=zeros(p,p);
a1=ones(p,p);
%Parallel transport vector
V1=zeros(2,p,p);
V2=zeros(p,p);
%Tangent vectors to geodesic curves
Vgamma=zeros(2,p,p);
Vgamma1=zeros(2,p,p);
%beta_j and alpha_i
b=zeros(p,p);
a=zeros(p,p);
%theta_{ij}=alpha_i-beta_j
theta=zeros(p,p);
%Discrete spin connection matrix
Omega=zeros(2,2,p);
P1=eye(2);
%Gamma matrices
gamma_x=[0 1; 1 0];
gamma_y=[0 -1i; 1i 0];
%Auxiliary to check the tetrad hypothesis
l=zeros(2,2,2*p);

%Set initial values
Rc(p)=R;

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
        phi(j,m)=angle(z);
        r(j)=abs(z);
        rx(j,m)=real(z);
        ry(j,m)=imag(z);
        %2nd direction
        z1=(e_j-x^2)/(1-x^2*e_j)*(r0*e_m-x*e_m1*(e_j-1)/(e_j-x^2))/(1-x/e_m1*r0*e_m*(1-e_j)/(1-x^2*e_j));
        phi(j+p,m)=angle(z1);
        r(j+p)=abs(z1);
        rx(j+p,m)=real(z1);
        ry(j+p,m)=imag(z1);
        %Set initial conditions
        yc(p,m)=(kappa+R)*imag(em1);
        xc(p,m)=(kappa+R)*real(em1);
        yc(2*p,m)=(kappa+R)*imag(e_m1);
        xc(2*p,m)=(kappa+R)*real(e_m1);
        %Plot sites of the 0th and 1st generation
        figure(1)
        scatter(real(z),imag(z),'filled','b');
        %scatter(xc(p,m),yc(p,m),'filled','r');
    end
end

for m=1:1:p
    for j=1:1:p
        %if m==1 && j==1
        if j==1
            gamma=atan((rx(j,m)-xc(p,m))/(yc(p,m)-ry(j,m)));
            gamma1=atan((rx(p,m)-xc(p,m))/(yc(p,m)-ry(p,m)));
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
        %if m==1 && j>=2
        if j>=2
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
                Rc(j-1)=sqrt(xc(j-1,m)^2+yc(j-1,m)^2-1);
                %Tangent vector to the circle in point z(j,m)
                gamma=atan((rx(j,m)-xc(j-1,m))/(yc(j-1,m)-ry(j,m)));
                %Tangent vector to the circle in point z(j-1,m)
                gamma1=atan((rx(j-1,m)-xc(j-1,m))/(yc(j-1,m)-ry(j-1,m)));
                %Auxiliary vectors
                tau1=(xc(j-1,m)-rx(j,m))/sqrt((xc(j-1,m)-rx(j,m))^2+(yc(j-1,m)-ry(j,m))^2);
                tau2=(yc(j-1,m)-ry(j,m))/sqrt((xc(j-1,m)-rx(j,m))^2+(yc(j-1,m)-ry(j,m))^2);
                tau3=(xc(j-1,m)-rx(j-1,m))/sqrt((xc(j-1,m)-rx(j-1,m))^2+(yc(j-1,m)-ry(j-1,m))^2);
                tau4=(yc(j-1,m)-ry(j-1,m))/sqrt((xc(j-1,m)-rx(j-1,m))^2+(yc(j-1,m)-ry(j-1,m))^2);
                if tau1*sin(gamma)-tau2*cos(gamma)<0
                    Vgamma(1,j,m)=cos(gamma);
                    Vgamma(2,j,m)=sin(gamma);
                else
                    Vgamma(1,j,m)=-cos(gamma);
                    Vgamma(2,j,m)=-sin(gamma);
                end
                if tau3*sin(gamma1)-tau4*cos(gamma1)<0
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

%Set initial conditions for parallely transported vectors
V1(1,1,1)=1-r(1)^2;
V1(2,1,1)=0;

%Define the parallely transported vectors
for m=1:1:p
    for j=1:1:p
        if j==1 && m>=2
            V1(1,j,m)=V1(1,p,m-1);
            V1(2,j,m)=V1(2,p,m-1);
        end
        if j>=2
            %1st direction
            if vpa(phi(j-1,m))==vpa(phi(j,m))
                V1(1,j,m)=(a1(j,m)*(1-r(j)^2))/(a1(j-1,m)*(1-r(j-1)^2))*V1(1,j-1,m);
                V1(2,j,m)=(a1(j,m)*(1-r(j)^2))/(a1(j-1,m)*(1-r(j-1)^2))*V1(2,j-1,m);
            else
                theta0(j-1,m)=phi(j-1,m)+atan((V1(1,j-1,m)*(r(j-1)^2+1)-V1(2,j-1,m)*sqrt(4*r(j-1)^2*(Rc(j-1)^2+1)-(r(j-1)^2+1)^2))/(V1(2,j-1,m)*(r(j-1)^2+1)+V1(1,j-1,m)*sqrt(4*r(j-1)^2*(Rc(j-1)^2+1)-(r(j-1)^2+1)^2)));
                a1(j,m)=2*r(j-1)*sqrt(Rc(j-1)^2+1)/((1-r(j-1)^2)*(sqrt(4*r(j-1)^2*(Rc(j-1)^2+1)-(r(j-1)^2+1)^2)*cos(theta0(j-1,m)-phi(j-1,m))+(r(j-1)^2+1)*sin(theta0(j-1,m)-phi(j-1,m))));
                V1(1,j,m)=a1(j,m)*(1-r(j)^2)/(2*r(j)*sqrt(Rc(j-1)^2+1))*(sqrt(4*r(j)^2*(Rc(j-1)^2+1)-(r(j)^2+1)^2)*cos(theta0(j-1,m)-phi(j,m))+(r(j)^2+1)*sin(theta0(j-1,m)-phi(j,m)));
                V1(2,j,m)=a1(j,m)*(1-r(j)^2)/(2*r(j)*sqrt(Rc(j-1)^2+1))*(-sqrt(4*r(j)^2*(Rc(j-1)^2+1)-(r(j)^2+1)^2)*sin(theta0(j-1,m)-phi(j,m))+(r(j)^2+1)*cos(theta0(j-1,m)-phi(j,m)));
            end
        end
        %Define angle beta_j, 1st direction
        if Vgamma(1,j,m)*V1(2,j,m)>=Vgamma(2,j,m)*V1(1,j,m)
            b(j,m)=acos((Vgamma(1,j,m)*V1(1,j,m)+Vgamma(2,j,m)*V1(2,j,m))/(a1(j,m)*(1-r(j)^2)));
        else
            b(j,m)=-acos((Vgamma(1,j,m)*V1(1,j,m)+Vgamma(2,j,m)*V1(2,j,m))/(a1(j,m)*(1-r(j)^2)));
        end
        if V1(1,j,m)*Vgamma1(2,j,m)>=V1(2,j,m)*Vgamma1(1,j,m)
            a(j,m)=acos((Vgamma1(1,j,m)*V1(1,j,m)+Vgamma1(2,j,m)*V1(2,j,m))/(a1(j,m)*(1-r(j)^2)));
        else
            a(j,m)=-acos((Vgamma1(1,j,m)*V1(1,j,m)+Vgamma1(2,j,m)*V1(2,j,m))/(a1(j,m)*(1-r(j)^2)));
        end
        if j>=2
            theta(j-1,m)=a(j-1,m)-b(j,m);
            if j==p
                theta(p,m)=a(p,m)-b(1,m);
            end
        end
    end
end

%Check Eq.(116), *Eq. (58)
for j=1:1:p
    Omega(:,:,j)=[exp(1i*theta(j,1)/2) 0; 0 exp(-1i*theta(j,1)/2)];
    P1=P1*Omega(:,:,j);
end

%Normalizability
for m=1:1:p
    for j=1:1:p
        V2(j,m)=(V1(1,j,m)^2+V1(2,j,m)^2)/(a1(j,m)*(1-r(j)^2))^2;
    end
end

%Check *Eq. (57)
for j=1:1:p
    if j==1
        l(:,:,p)=Vgamma1(1,p)*gamma_x+Vgamma1(2,p)*gamma_y;
        l(:,:,2*p)=-Omega(:,:,p)*ctranspose((-Vgamma(1,j)*gamma_x-Vgamma(2,j)*gamma_y)*Omega(:,:,p));  
    end
    if j>=2
        l(:,:,j-1)=Vgamma1(1,j-1)*gamma_x+Vgamma1(2,j-1)*gamma_y;
        l(:,:,j+p-1)=-Omega(:,:,j-1)*ctranspose((-Vgamma(1,j)*gamma_x-Vgamma(2,j)*gamma_y)*Omega(:,:,j-1));
    end
end
%Plot the geodesics curve between every two points of 0th and 1st
%generation
for j=1:1:p
    for m=1:1:p
        if j>=2 && vpa(phi(j-1,m))==vpa(phi(j,m))
            if phi(j-1,m)==0 || phi(j-1,m)==pi
                t=linspace(rx(j-1,m),rx(j,m),1000);
                x1=(ry(j,m)-ry(j-1,m))/(rx(j,m)-rx(j-1,m))*t;
                figure(1)
                plot(t,x1,'k','LineWidth',1.5);
            else
                t=linspace(ry(j-1,m),ry(j,m),1000);
                y1=(rx(j,m)-rx(j-1,m))/(ry(j,m)-ry(j-1,m))*t;
                figure(1)
                plot(y1,t,'k','LineWidth',1.5);
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
        
        
   

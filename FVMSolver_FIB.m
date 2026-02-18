%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [qL, qR] = reconstruct(Q, q_exact, order, limiter, dx)

end

function source = computeSource(q, gam, dx, x)
    np     = size(q,2);
    source = zeros(size(q));
    for i = 1:np
        rho  = q(1,i); u = q(2,i)/rho; E = q(3,i);
        p    = (gam-1)*(E - 0.5*rho*u^2);
        A    = area(x(i));
        dAdx = (area(x(i)+dx/2) - area(x(i)-dx/2))/dx;
        source(1,i) = -rho*u*dAdx/A;
        source(2,i) = -rho*u^2*dAdx/A;
        source(3,i) = -u*(E+p)*dAdx/A;
    end
end

function fluxes = computeFluxes(qL, qR, gam, dx)
   
end

function aa = area(x)
    if x > 0
        aa = 0.536572 - 0.198086*exp(-log(2)*(x/0.6)^2);
    else
        aa = 1.0 - 0.661514*exp(-log(2)*(x/0.6)^2);
    end
end

function q = ExactSolu(x, flow)
    global x_shock p0 rho0 astar p02 astar2
    gam=1.4; gamm1=gam-1.; q=zeros(3,1);
    warning('off','all');
    if flow == 2
        T0 = p0/rho0;
        if x <= x_shock
            A=area(x); aas2=(A/astar)^2;
            if aas2 < 1.01
                if x>=0; machs=fzero(@(m)MachArea(m,aas2),1.1);
                else;     machs=fzero(@(m)MachArea(m,aas2),0.99); end
            else
                if x>0;  machs=fzero(@(m)MachArea(m,aas2),2);
                else;     machs=fzero(@(m)MachArea(m,aas2),0.1); end
            end
            temp=1+0.5*gamm1*machs^2; p1=p0/temp^(gam/gamm1);
            rho1=rho0/temp^(1/gamm1); c=sqrt(gam*p1/rho1); u1=abs(machs)*c;
        else
            aas2=(area(x)/astar2)^2; mache=fzero(@(m)MachArea(m,aas2),0.05);
            temp=1+0.5*gamm1*mache^2; p1=p02/temp^(gam/gamm1);
            T1=T0/temp; rho1=p1/T1; c=sqrt(gam*T1); u1=abs(mache)*c;
        end
        q(1)=rho1; q(2)=rho1*u1; q(3)=p1/gamm1+0.5*rho1*u1^2;
    elseif flow == 1
        temp=(area(x)/astar)^2; mach=fzero(@(m)MachArea(m,temp),0.05);
        temp=1+0.5*gamm1*mach^2; rho=rho0/temp^(1/gamm1);
        p=p0/temp^(gam/gamm1); c=sqrt(gam*p/rho); u=abs(mach)*c;
        q(1)=rho; q(2)=rho*u; q(3)=p/gamm1+0.5*rho*u^2;
    end
    warning('on','all');
end

function [] = Initialize(flow)
    global x_shock p0 rho0 astar p02 astar2
    gam=1.4; gamm1=gam-1.; gamp1=gam+1.;
    if flow == 2
        mach_i=0.2006533; rho_i=1.; p_i=1./gam; p_e=0.6071752;
        temp=1+0.5*gamm1*mach_i^2;
        p0=p_i*temp^(gam/gamm1); rho0=rho_i*temp^(1/gamm1);
        x_shock=0.5; x_shock_small=0.; x_shock_large=1.; pe=0.; astar=area(0);
        while abs(p_e-pe) > 0.00001
            A=area(x_shock); aas2=(A/astar)^2;
            machs=fzero(@(m)MachArea(m,aas2),2.);
            temp=1+0.5*gamm1*machs^2; p1=p0/temp^(gam/gamm1);
            p2=p1*(1+2*gam/gamp1*(machs^2-1));
            mach2=sqrt((1+0.5*gamm1*machs^2)/(gam*machs^2-0.5*gamm1));
            aas2=1/mach2^2*(2/gamp1*(1+0.5*gamm1*mach2^2))^(gamp1/gamm1);
            astar2=sqrt(A^2/aas2);
            temp=1+0.5*gamm1*mach2^2; p02=p2*temp^(gam/gamm1);
            aas2=(area(4)/astar2)^2; mache=fzero(@(m)MachArea(m,aas2),0.05);
            temp=1+0.5*gamm1*mache^2; pe=p02/temp^(gam/gamm1);
            if pe>p_e; x_shock_small=x_shock; x_shock=0.5*(x_shock+x_shock_large);
            else;       x_shock_large=x_shock; x_shock=0.5*(x_shock+x_shock_small); end
        end
        fprintf('Shock location: x = %.4f\n', x_shock);
    elseif flow == 1
        area_e=area(5.0); mach_e=0.4; rho_e=1.; p_e=1./gam;
        temp=1+0.5*gamm1*mach_e^2;
        p0=p_e*temp^(gam/gamm1); rho0=rho_e*temp^(1/gamm1);
        temp=1/mach_e^2*(2/gamp1*(1+0.5*gamm1*mach_e^2))^(gamp1/gamm1);
        astar=area_e/sqrt(temp);
    end
end

function aas = MachArea(mach, aas2)
    gam=1.4; gamp1=gam+1; gamm1=gam-1;
    aas = 1/mach^2*(2/gamp1*(1+0.5*gamm1*mach^2))^(gamp1/gamm1) - aas2;
end
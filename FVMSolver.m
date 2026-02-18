clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Transonic Nozzle Flow - 2nd Order FV, Rusanov Flux, SSP-RK2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PARAMETERS
flow          = 1;
order         = 2;       % 1 = first order, 2 = second order
limiter       = 1;       % 1 = no limiter, 2 = limiter on
np_values     = [100, 200, 400, 800];
tolerance     = 1e-6;
max_iterations = 10000;
CFL           = 1;
xmin          = -4;
xmax          = 4;
gam           = 1.4;
gamm1         = gam - 1;
gamp1         = gam + 1;

%% INITIALIZE
Initialize(flow);

area_e = area(5.0); mach_e = 0.4; rho_e = 1.0; p_e = 1.0/gam;
temp   = 1 + 0.5*gamm1*mach_e^2;
p0     = p_e*temp^(gam/gamm1);
rho0   = rho_e*temp^(1/gamm1);
temp   = 1.0/mach_e^2*(2/gamp1*(1 + 0.5*gamm1*mach_e^2))^(gamp1/gamm1);
astar  = area_e/sqrt(temp);

green_colors = [linspace(0.5,0,length(np_values))', ...
                linspace(1,0.5,length(np_values))', ...
                linspace(0.5,0,length(np_values))'];

%% MESH LOOP
figure;
for idx = 1:length(np_values)
    np = np_values(idx);
    dx = (xmax - xmin)/(np - 1);
    x  = linspace(xmin, xmax, np);

    % Exact solution
    q_exact = zeros(3, np);
    for i = 1:np
        q_exact(:,i) = ExactSolu(x(i), flow);
    end

    % Initial guess (inlet state everywhere)
    Q = repmat(q_exact(:,1), 1, np);
    Q(:,1)   = q_exact(:,1);
    Q(:,end) = q_exact(:,end);

    %% ITERATION LOOP (SSP-RK2)
    residual  = inf;
    iteration = 0;

    while residual > tolerance && iteration < max_iterations
        Q(:,1)       = q_exact(:,1);
        Q(:,end)     = q_exact(:,end);
        Q(:,end-2:end) = q_exact(:,end-2:end);

        rho = Q(1,:); u = Q(2,:)./rho;
        p   = gamm1*(Q(3,:) - 0.5*Q(2,:).*u);
        c   = sqrt(gam*p./rho);
        dt  = CFL*min(dx./(abs(u) + c));

        % RK2 Stage 1
        [qL1, qR1] = reconstruct(Q, q_exact, order, limiter, dx);
        F1         = computeFluxes(qL1, qR1, gam, dx);
        S1         = computeSource(Q, gam, dx, x);
        Q_star     = Q + dt*(-F1 + S1);
        Q_star(:,1) = q_exact(:,1); Q_star(:,end-2:end) = q_exact(:,end-2:end);

        % RK2 Stage 2
        [qL2, qR2] = reconstruct(Q_star, q_exact, order, limiter, dx);
        F2         = computeFluxes(qL2, qR2, gam, dx);
        S2         = computeSource(Q_star, gam, dx, x);
        Q_new      = 0.5*(Q + Q_star + dt*(-F2 + S2));
        Q_new(:,1) = q_exact(:,1); Q_new(:,end-2:end) = q_exact(:,end-2:end);

        residual  = max(abs(Q_new - Q), [], 'all');
        Q         = Q_new;
        iteration = iteration + 1;
    end

    fprintf('np=%d  iter=%d  residual=%.3e\n', np, iteration, residual);

    % Derived quantities for plotting
    vel_n = Q(2,:)./Q(1,:);
    vel_e = q_exact(2,:)./q_exact(1,:);

    color = green_colors(idx,:);
    subplot(4,1,1);
    plot(x, Q(1,:), 'Color', color, 'DisplayName', ['np=', num2str(np)]); hold on;
    ylabel('Density'); title('Density'); grid on;

    subplot(4,1,2);
    plot(x, Q(2,:), 'Color', color, 'DisplayName', ['np=', num2str(np)]); hold on;
    ylabel('Momentum'); title('Momentum'); grid on;

    subplot(4,1,3);
    plot(x, Q(3,:), 'Color', color, 'DisplayName', ['np=', num2str(np)]); hold on;
    ylabel('Energy'); title('Energy'); grid on;

    subplot(4,1,4);
    plot(x, vel_n, 'Color', color, 'DisplayName', ['np=', num2str(np)]); hold on;
    ylabel('Velocity'); title('Velocity'); xlabel('x'); grid on;
end

% Overlay exact solution (using finest mesh x)
subplot(4,1,1); plot(x, q_exact(1,:), 'k--', 'DisplayName', 'Exact'); legend();
subplot(4,1,2); plot(x, q_exact(2,:), 'k--', 'DisplayName', 'Exact'); legend();
subplot(4,1,3); plot(x, q_exact(3,:), 'k--', 'DisplayName', 'Exact'); legend();
subplot(4,1,4); plot(x, q_exact(2,:)./q_exact(1,:), 'k--', 'DisplayName', 'Exact'); legend();
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [qL, qR] = reconstruct(Q, q_exact, order, limiter, dx)
    np = size(Q, 2);
    qL = zeros(3, np+1);
    qR = zeros(3, np+1);

    % Exact boundary states
    qL(:,1)   = q_exact(:,1);   qR(:,1)   = q_exact(:,1);
    qL(:,end) = q_exact(:,end); qR(:,end) = q_exact(:,end);

    if order == 1
        for i = 2:np
            qL(:,i) = Q(:,i-1);
            qR(:,i) = Q(:,i);
        end

    elseif order == 2
        slope = zeros(3, np);
        slope(:,1)   = (Q(:,2)   - Q(:,1))  /dx;
        slope(:,end) = (Q(:,end) - Q(:,end-1))/dx;

        if limiter == 2
            for j = 1:3
                for i = 2:np-1
                    if Q(j,i) < Q(j,i+1) && Q(j,i) < Q(j,i-1)
                        slope(j,i) = 0;
                    elseif Q(j,i) > Q(j,i+1) && Q(j,i) > Q(j,i-1)
                        slope(j,i) = 0;
                    elseif Q(j,i) == Q(j,i+1) || Q(j,i) == Q(j,i-1)
                        slope(j,i) = 0;
                    else
                        slope(j,i) = min(abs(Q(j,i)-Q(j,i-1)), abs(Q(j,i+1)-Q(j,i)));
                    end
                end
            end
        else
            for i = 2:np-1
                slope(:,i) = (Q(:,i+1) - Q(:,i-1))/(2*dx);
            end
        end

        for i = 2:np
            if limiter == 1 && i >= 3 && i <= np-1
                qL(:,i) = Q(:,i-1) + 0.25*(Q(:,i)   - Q(:,i-2));
                qR(:,i) = Q(:,i)   - 0.25*(Q(:,i+1) - Q(:,i-1));
            else
                qL(:,i) = Q(:,i-1) + 0.5*dx*slope(:,i-1);
                qR(:,i) = Q(:,i)   - 0.5*dx*slope(:,i);
            end
        end
    end
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
    nfaces = size(qL, 2);   % np+1 faces
    ncells = nfaces - 1;

    % Compute Rusanov flux at each face
    F_face = zeros(3, nfaces);
    for i = 1:nfaces
        rhoL = qL(1,i); uL = qL(2,i)/rhoL; EL = qL(3,i);
        pL   = (gam-1)*(EL - 0.5*rhoL*uL^2); cL = sqrt(gam*pL/rhoL);
        FL   = [rhoL*uL; rhoL*uL^2+pL; uL*(EL+pL)];

        rhoR = qR(1,i); uR = qR(2,i)/rhoR; ER = qR(3,i);
        pR   = (gam-1)*(ER - 0.5*rhoR*uR^2); cR = sqrt(gam*pR/rhoR);
        FR   = [rhoR*uR; rhoR*uR^2+pR; uR*(ER+pR)];

        smax       = max(abs(uL)+cL, abs(uR)+cR);
        F_face(:,i) = 0.5*(FL+FR) + 0.5*smax*(qL(:,i)-qR(:,i));
    end

    % Flux divergence per cell
    fluxes = zeros(3, ncells);
    for i = 1:ncells
        fluxes(:,i) = (F_face(:,i+1) - F_face(:,i)) / dx;
    end
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
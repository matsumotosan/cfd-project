function [x,y] = create_mesh(xi,eta,plotmesh,plot_bc)
%CREATE_MESH 

% Create meshgrid from eta and xi ranges
[xixi,ee] = meshgrid(xi,eta);

% Calculate Cartesian coordiantes of ellipse
x = cosh(ee) .* cos(xixi);
y = sinh(ee) .* sin(xixi);

% Plot mesh in x-y plane
if plotmesh
    figure
    plot(x,y,'-k'); hold on;    % plot constant eta
    plot(x',y','-k')            % plot constant xi
    title('Elliptical mesh in x-y plane','FontSize',18)
    xlabel('x','FontSize',14)
    ylabel('y','FontSize',14)
    axis equal tight
end

% Mesh with highlighted boundaries
if plot_bc
    figure
    subplot(1,2,1)
    plot(x,y,'-k'); hold on;    % plot constant eta
    plot(x',y','-k')            % plot constant xi
    plot(x(:,1),y(:,1),'-r','LineWidth',2); hold on
    plot(x(:,end),y(:,end),'-g','LineWidth',2)
    plot(x(1,:),y(1,:),'-c','LineWidth',2)
    plot(x(end,:),y(end,:),'-b','LineWidth',2)
    xlabel('x','FontSize',14)
    ylabel('y','FontSize',14)
    title('Boundaries in x-y plane','FontSize',18)
    axis equal tight

    subplot(1,2,2)
    plot(xixi,ee,'-k'); hold on;    % plot constant eta
    plot(xixi',ee','-k')            % plot constant xi
    plot(xixi(:,1),ee(:,1),'-r','LineWidth',2); hold on
    plot(xixi(:,end),ee(:,end),'-g','LineWidth',2)
    plot(xixi(1,:),ee(1,:),'-c','LineWidth',2)
    plot(xixi(end,:),ee(end,:),'-b','LineWidth',2)
    xlabel('\xi','FontSize',14)
    ylabel('\eta','FontSize',14)
    title('Boundaries in \xi-\eta plane','FontSize',18)
    axis equal tight
end

end


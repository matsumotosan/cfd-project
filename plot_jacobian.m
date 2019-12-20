function plot_jacobian(x,y,xi,eta,J,opt)
%PLOT_JACOBIAN Plot Jacobian of grid transformation on x-y and xi-eta
%planes.

if opt == 1
    str = '(analytical)';
else
    str = '(numerical)';
end

figure
subplot(2,1,1)
contourf(x,y,J,15); hold on
title(['Jacobian of grid transformation on x-y plane ' str],'FontSize',14)
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
colorbar

% Plot Jacobian (eta-xi plane)
subplot(2,1,2)
contourf(xi,eta,J,15); hold on
title(['Jacobian of grid transformation on \xi-\eta plane ' str],'FontSize',14)
xlabel('\xi','FontSize',14)
ylabel('\eta','FontSize',14)
colorbar

end


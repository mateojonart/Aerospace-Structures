function PlotSection(x_path,Tn_path,sigma,units)

    % Precomputations
    X = x_path(:,1);
    Y = x_path(:,2);
    smin = min(sigma(1,:));
    smax = max(sigma(1,:));
    
    % Open plot window
    figure; box on; hold on; axis equal;
    % Plot deformed structure with colorbar for stresses 
    patch(X(Tn_path'),Y(Tn_path'),[sigma(1,:);sigma(1,:)],'EdgeColor','interp','LineWidth',4);
    % Colorbar settings
    clims = get(gca,'clim');
    clim(max(abs(clims))*[-1,1]);
    n = 128; % Number of rows
    c1 = 2/3; % Blue
    c2 = 0; % Red
    s = 0.85; % Saturation
    c = hsv2rgb([c1*ones(1,n),c2*ones(1,n);1:-1/(n-1):0,1/n:1/n:1;s*ones(1,2*n)]');
    colormap(c); 
    cb = colorbar;
    
    % Add labels
    title(sprintf('sigma_{min} = %.3g %s | sigma_{max} = %.3g %s', smin, units, smax, units)); 
    xlabel('x (m)'); 
    ylabel('y (m)');
    cb.Label.String = sprintf('Stress (%s)',units); 

end

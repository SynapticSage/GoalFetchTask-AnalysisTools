
%------------------------------------------------------------
%  ____  _       _       
% |  _ \| | ___ | |_ ___ 
% | |_) | |/ _ \| __/ __|
% |  __/| | (_) | |_\__ \
% |_|   |_|\___/ \__|___/
%------------------------------------------------------------
% Consider moving this to a new function!
%                        
% If allowed positive and negative modulation, doesn't appear to
% track rewarded areas!

fig('Change of angle, no goal to goal'); 
tiledlayout(5,2); 
C = cat(2, conditionals{:});
for i = 2:numel(RC)
    nexttile; 
    q=RC{i}-RC{1}; 
    histogram(q(:,2)); 
    location = C(i,:);
    location(2) = location(2) + 1;
    location = num2cell(location);
    skew(location{:}) = skewness(q(:,2)); 
    title(...
        sprintf("%d " , C(i, :)) + newline + ...
        sprintf("skewness = %2.1f", skew(location{:})))
    xlim([-pi, pi]); 
end

fig('Change of modulation, no goal to goal'); 
tiledlayout(5,2); 
C = cat(2, conditionals{:});
for i = 2:numel(RC)
    nexttile; 
    q=RC{i}-RC{1}; 
    prop = q(:,1);
    histogram(prop); 
    location = C(i,:);
    location(2) = location(2) + 1;
    location = num2cell(location);
    skew(location{:}) = skewness(prop); 
    title(...
        sprintf("%d " , C(i, :)) + newline + ...
        sprintf("skewness = %2.1f", skew(location{:})))
    xlim([-0.2, 0.2]); 
end

fig('Change of x,y, no goal to goal'); 
tiledlayout(5,2); 
C = cat(2, conditionals{:});
clear skew
for i = 2:numel(RC)
    nexttile; 
    q    = RC{i}
    prop = q(:,3:4);
    P = num2cell(prop, 1);
    plot(beh.pos(:,1), beh.pos(:,2), 'k:', 'linewidth', 0.01)
    hold on
    scatter(P{:}, '*'); 
    location = C(i,:);
    location(2) = location(2) + 1;
    location = num2cell(location);
    location = [location, ':'];
    skew(location{:}) = skewness(prop);
    axis square
    title(sprintf("%d " , C(i, :)))
    xlim([-20, 200]);
    ylim([-20, 100]);
end


%% TO POINTS OF INTEREST? BARRIERS?
%% --------------------------------


%% Possible comparisons
%% --------------------
% Just listing out comparisons that would be nice. 
%   - Mean intensitiy of goal versus non
%   - Different symmetries
%       - Egocentric
%           - Collapse ego coding to certain sites?
%               - mirror symmetry
%               - rotational symmmetry
%       - Allocentricc
%           - Mirror symmetry
%           - Rotational symmetry
%   - Comparison to shuffle
%   - Comparison to head direction coding

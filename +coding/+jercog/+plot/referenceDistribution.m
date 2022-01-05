function g = referenceDistribution(params, beh, varargin) % 
% Standardizes plots distribution of reference

ip = inputParser;
ip.addParameter('pos', 'back');
ip.addParameter('bounds', 'behavior');
ip.addParameter('subset', []);
ip.addParameter('gargs', {});
ip.addParameter('fwargs', {});
ip.addParameter('fgargs', {});
ip.addParameter('alpha_behavior',     0.5);
ip.addParameter('alpha_distribution', 0.75);
ip.addParameter('padding', 40);
ip.parse(varargin{:})
Opt = ip.Results;

% Cast into proper table if need be
if isstruct(params)
    tab = coding.vectorCells.jercog.table.params(params);
else
    tab = params;
end

% Process arg consequences
if ~isempty(Opt.bounds)
    bounds = gbehavior.mazebounds(beh);
    subset = util.constrain.minmax(tab.reference_point_x, bounds.x + [-Opt.padding, Opt.padding]) & ...
             util.constrain.minmax(tab.reference_point_y, bounds.y + [-Opt.padding,Opt.padding]) & ...
             tab.amplitude > 0;
end
if ~isempty(Opt.subset)
    subset = subset & Opt.subset;
end

% Create table
g = gramm('x', tab.reference_point_x,...
          'y', tab.reference_point_y,...
          'subset', subset, ...
        Opt.gargs{:});

if ~isempty(Opt.fwargs)
    Opt.fwargs = cellfun(@(x) tab.(x), Opt.fwargs, 'UniformOutput', false);
        g = g.facet_wrap(Opt.fwargs{:});
end
if ~isempty(Opt.fgargs)
    Opt.fgargs = cellfun(@(x) tab.(x), Opt.fgargs, 'UniformOutput', false);
    g = g.facet_grid(Opt.fgargs{:});
end


g = g.stat_bin2d('geom', 'image');
g.axe_property('PlotBoxAspectRatio', [1,1,1], ...,
'xlim', bounds.x+[-Opt.padding, Opt.padding], ...
'ylim', bounds.y+[-Opt.padding, Opt.padding]);
g = g.draw();
g = callback(g, beh, Opt);


%SizeChangedFcn = get(gcf, 'SizeChangedFcn');
%if ~isempty(SizeChangedFcn)
%    if ~iscell(SizeChangedFcn)
%        SizeChangedFcn = {SizeChangedFcn};
%    end
%else
%    SizeChangedFcn = {};
%end
%SizeChangedFcn{end+1} = @() callback(g, Opt);


function g = callback(g, beh, Opt)
    if ~isempty(Opt.pos) && ~isempty(g.facet_axes_handles)

        for ax = g.facet_axes_handles(:)'

            hold(ax, 'on'); 
            s = scatter(ax, beh.pos(:,1), beh.pos(:,2), '.k');

            if strcmpi(Opt.pos, 'front')
                C = 100;
            elseif strcmpi(Opt.pos, 'back')
                C = -100;
            end
            set(s,'ZData', C*ones(size(s.XData)))
        end

        if ~isempty(Opt.alpha_behavior)
            set(s, 'MarkerEdgeAlpha', 1, 'MarkerFaceAlpha', Opt.alpha_behavior);
        end
        if ~isempty(Opt.alpha_distribution)
            alpha(findobj(g.facet_axes_handles, 'Type', 'Patch'), Opt.alpha_distribution);
        end

    end



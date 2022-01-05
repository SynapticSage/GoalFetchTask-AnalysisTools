config.clean_string = '$area == "PFC" | $area == "CA1"';

% ==================
% JERCOG: ALL TIMES
% ==================

params = coding.jercog.table.params(spikes.jercog, 'addCellTable', spikes.cellTable);
params = util.table.query(params, config.clean_string);
params = util.table.string2categorical(params);

% Amplitudes
g = gramm('x', params.amplitude, 'color', params.area);
gfig('Jercog.Amplitude', 'grammObj', g, 'clf', true)
g = g.stat_bin('geom','line', 'normalization', 'cdf', 'fill', 'transparent')
g.axe_property('xlim',[-1, 1.5]);
g.draw()
% text(g.facet_axes_handles, 0.1, 0.8, "Slight biase to higher" + newline + "for pfc", 'Units', 'Normalized')

% Angles
g = gramm('x', params.theta_preference, 'color', params.area);
gfig('Jercog angle preference', 'grammObj', g, 'clf', true)
g = g.stat_bin('geom','bar', 'normalization', 'probability', 'fill', 'transparent')
g = g.set_text_options('interpreter', 'latex')
g=g.axe_property('xlim',[-1.5, 1.5]);
g =g.set_names('x', '$\theta$')
%g.set_polar();
g.draw()
text(g.facet_axes_handles, 0.1, 0.8, "Both biased straight", 'Units', 'Normalized')

% Place field centroid to reference point
g(1,1) = gramm('x', params.ref_to_cent, 'color', params.area, 'subset', abs(params.ref_to_cent) < 200 & abs(params.amplitude) > 0);
gfig('Jercog Reference to Center', 'grammObj', g, 'clf', true)
g(1,1) = g(1,1).stat_bin('geom','line', 'normalization', 'pdf', 'fill', 'transparent', 'nbins', 10)
g(1,2) = gramm('x', params.ref_to_cent, 'color', params.area, 'subset', abs(params.ref_to_cent) < 200);
g(1,2) = g(1,2).stat_bin('geom','line', 'normalization', 'cdf', 'fill', 'all', 'nbins', 10)
g.draw()

% Fraction of cells with ref modulation?
res_all = splitapply(@numel, params.neurons, findgroups(params.area));
filt = abs(params.amplitude) > 0.05;
res_amp = splitapply(@numel, params(filt,:).neurons, findgroups(params(filt,:).area));
clear g
g = gramm('x', categorical(["CA1", "PFC"]), 'y', res_amp./res_all)
gfig("Percent of cells with modulation > 5%", 'grammObj', g, 'clf', true)
g.geom_bar()
g.axe_property('ylim', [0, 1]);
g.set_color_options('chroma', 0, 'lightness', 0);
g.set_names('x', 'Areas', 'y', 'Percentage Above Thresh');
g.draw()

% Relation of cent_to_ref and amplitude?
clear g
g(1) = gramm('x', abs(params.amplitude), 'y', params.ref_to_cent, 'color', params.area, 'subset', abs(params.ref_to_cent) < 200 & abs(params.amplitude) > 0.1);
gfig("Jercog Reference Dist from Place Field" + newline + "vs Amplitude", 'grammObj', g, 'clf', true)
g(1).facet_wrap(params.area)
g(1) = g(1).geom_point()
g(1).axe_property('PlotBoxAspectRatio', [1,1,1]);
g(1) = g(1).stat_glm()
g.set_names('x', 'Absolute Amplitude', 'y', 'Reference point to Place Field', 'column', '')
g.draw()

% Let's plot the reference point distributions
gfig("Ref point dist");clf
coding.jercog.plot.referenceDistribution(params, beh);

gfig("Ref point dist-area");clf
coding.jercog.plot.referenceDistribution(params, beh, 'fwargs', {'area'}, 'padding', 100);


% Let's plot the reference point scatters
gfig("Ref point dist");clf
coding.jercog.plot.referenceDistribution(params, beh);

gfig("Ref point dist-area");clf
coding.jercog.plot.referenceDistribution(params, beh, 'fwargs', {'area'}, 'padding', 100);

% Place field centroid to reference point
%g(1,1) = gramm('x', params.ref_to_cent_x, 'color', params.area, 'subset', abs(params.ref_to_cent_x) < 200);
%gfig('Jercog Reference to Center', 'grammObj', g, 'clf', true)
%g(1,1) = g(1,1).stat_bin('geom','line', 'normalization', 'pdf', 'fill', 'transparent', 'nbins', 10)
%g(1,1) = g(1,1).set_names('x', 'X_{reference to center}');
%g(1,2) = gramm('x', params.ref_to_cent_x, 'color', params.area, 'subset', abs(params.ref_to_cent_x) < 200);
%g(1,2) = g(1,2).stat_bin('geom','line', 'normalization', 'cdf', 'fill', 'all', 'nbins', 10)
%g(1,2) = g(1,2).set_names('x', 'X_{reference to center}');
%g(2,1) = gramm('x', params.ref_to_cent_y, 'color', params.area, 'subset', abs(params.ref_to_cent_y) < 200);
%g(2,1) = g(2,1).stat_bin('geom','line', 'normalization', 'pdf', 'fill', 'transparent', 'nbins', 10)
%g(2,1) = g(2,1).set_names('x', 'Y_{reference to center}');
%g(2,2) = gramm('x', params.ref_to_cent_y, 'color', params.area, 'subset', abs(params.ref_to_cent_y) < 200);
%g(2,2) = g(2,2).stat_bin('geom','line', 'normalization', 'cdf', 'fill', 'all', 'nbins', 10)
%g(2,2) = g(2,2).set_names('x', 'Y_{reference to center}');
%g = g.set_text_options('interpreter', 'latex')
%g.draw()

% =====================
% JERCOG: Split by goal
% =====================

params = coding.jercog.table.params(spikes.jercog.byGoal, 'addCellTable', spikes.cellTable);
params = util.table.query(params, config.clean_string);
params = util.table.string2categorical(params);
params = params(params.correct == true,:);

% Amplitudes
g = gramm('x', params.amplitude, 'color', params.area);
g.facet_wrap(params.stopWell)
gfig('Jercog.GoalWise.Amplitude', 'grammObj', g, 'clf', true)
g = g.stat_bin('geom','line', 'normalization', 'pdf', 'fill', 'transparent')
%g.axe_property('xlim',[-1, 1.5]);
g.draw()

% Angles
g = gramm('x', params.theta_preference, 'color', params.area);
gfig('Jercog angle preference', 'grammObj', g, 'clf', true)
g.facet_wrap(params.stopWell, 'ncols', 2)
g = g.stat_bin('geom','bar', 'normalization', 'probability', 'fill', 'transparent')
g = g.set_text_options('interpreter', 'latex')
g=g.axe_property('xlim',[-1, 1.5]);
g =g.set_names('x', '$\theta$')
%g.set_polar();
g.draw()
text(g.facet_axes_handles, 0.1, 0.8, "Both biased straight", 'Units', 'Normalized')

for g = 1:5
    gfig("Ref point dist-area-splitByGoal" + g);
    clf
    coding.jercog.plot.referenceDistribution(params(params.stopWell==g & params.amplitude>0.1,:), beh(beh.stopWell==g,:), 'fwargs', {'area'}, 'padding', 5);
end


function fig = new_figure(cfg, name, w, h)
% NEW_FIGURE  Create a figure (visible or hidden depending on cfg.show_figures).
if nargin < 3, w = 1000; end
if nargin < 4, h = 600;  end
vis = 'off';
if cfg.show_figures, vis = 'on'; end
fig = figure('Name', name, 'Visible', vis, 'Color', 'w', ...
             'Position', [100 100 w h]);
end

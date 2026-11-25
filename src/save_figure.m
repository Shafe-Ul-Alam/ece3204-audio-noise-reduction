function save_figure(fig, cfg, fname)
% SAVE_FIGURE  Save a figure as PNG into results/figures and close it if hidden.
file = fullfile(cfg.fig_dir, [fname '.png']);
print(fig, file, '-dpng', '-r150');
fprintf('   saved figure  %s\n', file);
if ~cfg.show_figures
    close(fig);
end
end

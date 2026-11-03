function sg_save_figure(fig, name, outdir, dpi)
%SG_SAVE_FIGURE Save vector PDF and raster PNG
if nargin < 4, dpi = 180; end
if ~exist(outdir,'dir'), mkdir(outdir); end
set(fig, 'PaperPositionMode', 'auto');
print(fig, fullfile(outdir, [name '.pdf']), '-dpdf', '-painters');
print(fig, fullfile(outdir, [name '.png']), '-dpng', sprintf('-r%d',dpi));
end

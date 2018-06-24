function screen2tiff(filename)

if nargin < 1

     error('Not enough input arguments!')

end

oldscreenunits = get(gcf,'Units');

oldpaperunits = get(gcf,'PaperUnits');

oldpaperpos = get(gcf,'PaperPosition');

set(gcf,'Units','pixels');

scrpos = get(gcf,'Position');

newpos = scrpos/100;

set(gcf,'PaperUnits','inches','PaperPosition',newpos)

print('-dtiff', filename, '-r100');

drawnow

set(gcf,'Units',oldscreenunits,'PaperUnits',oldpaperunits, 'PaperPosition',oldpaperpos)
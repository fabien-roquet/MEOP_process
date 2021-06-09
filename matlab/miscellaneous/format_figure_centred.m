% script format figure
function format_figure_centred(PP)

set(gcf,'PaperUnits','normalized');
set(gcf,'PaperType','usletter');
orient tall

if exist('PP') & PP~=0
set(gcf,'PaperUnits','centimeters');
PP2=[10.5-PP(1)/2 15-PP(2)/2 PP(1) PP(2)];
set(gcf,'PaperPosition',PP2);
else
orient tall
end
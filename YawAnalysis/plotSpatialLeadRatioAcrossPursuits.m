function S=plotSpatialLeadRatioAcrossPursuits(opts)
% plotSpatialLeadRatioAcrossPursuits  Supplementary spatial-lead ratio plot.
%
% The horizontal head-forward ray intersects the time-oriented target path at
% signed spatial lead d. This figure shows the dimensionless geometric ratio
% d/R, where R is current head-origin-to-target distance. Framewise points and
% run medians use the nominal reconstruction. Capped intervals propagate the
% complete joint pose/lure Monte Carlo reconstructions.

arguments
    opts struct=struct()
end

thisDir=fileparts(mfilename('fullpath'));
resultDir=fullfile(fileparts(thisDir),'outputs','yaw','spatial');
opts=applyDefaults(opts,struct( ...
    'nominalCsv',fullfile(resultDir,'yaw_spatial_lookahead_by_frame.csv'), ...
    'monteCarloMat',fullfile(resultDir,'MonteCarloResults', ...
        'yaw_spatial_lookahead_monte_carlo_analysis.mat'), ...
    'outputDir',resultDir, ...
    'showFigure',false, ...
    'figureResolution',600));

assert(isfile(opts.nominalCsv),'Missing nominal spatial-lead table: %s', ...
    opts.nominalCsv);
assert(isfile(opts.monteCarloMat),'Missing Monte Carlo spatial-lead result: %s', ...
    opts.monteCarloMat);
if ~exist(opts.outputDir,'dir'),mkdir(opts.outputDir);end

runNames=["julesFlick1_20190309";"julesTopRun11_20170829"; ...
    "julesTopRun1_20171209";"lilyRun_20190309"; ...
    "phantomRun_20190303"];
runLabels=["Jules Flick 1";"Jules Top Run 11";"Jules Top Run 1"; ...
    "Lily Run";"Phantom Run"];
colors=[0 114 178;213 94 0;0 158 115;204 121 167;230 159 0]/255;

T=readtable(opts.nominalCsv);
required={'dataset','absFrame','analysisMask','leadDistanceM', ...
    'headOriginToLureDistanceM'};
assert(all(ismember(required,T.Properties.VariableNames)), ...
    'Nominal table does not contain the required spatial-lead fields.');

M=load(opts.monteCarloMat,'R');
assert(isfield(M,'R')&&isfield(M.R,'runs'), ...
    'Invalid Monte Carlo spatial-lead result: %s',opts.monteCarloMat);

nRun=numel(runNames);
ratioByRun=cell(nRun,1);
nFrame=zeros(nRun,1);
q25Ratio=nan(nRun,1);
medianRatio=nan(nRun,1);
q75Ratio=nan(nRun,1);
mcMedianRatio=nan(nRun,1);
mcLower95Ratio=nan(nRun,1);
mcUpper95Ratio=nan(nRun,1);

for j=1:nRun
    keep=string(T.dataset)==runNames(j)&logical(T.analysisMask)& ...
        isfinite(T.leadDistanceM)&isfinite(T.headOriginToLureDistanceM)& ...
        T.headOriginToLureDistanceM>0;
    if runNames(j)=="phantomRun_20190303"
        keep=keep&T.absFrame>=100&T.absFrame<=216;
    end
    ratio=double(T.leadDistanceM(keep)./T.headOriginToLureDistanceM(keep));
    assert(~isempty(ratio),'No valid spatial-lead ratios for %s.',runNames(j));
    ratioByRun{j}=ratio;
    nFrame(j)=numel(ratio);
    q25Ratio(j)=percentile(ratio,25);
    medianRatio(j)=percentile(ratio,50);
    q75Ratio(j)=percentile(ratio,75);

    qIndex=find(cellfun(@(q)string(q.dataset)==runNames(j),M.R.runs),1);
    assert(~isempty(qIndex),'Monte Carlo result is missing %s.',runNames(j));
    Q=M.R.runs{qIndex};
    frameKeep=logical(Q.frameTable.analysisMask);
    if runNames(j)=="phantomRun_20190303"
        frameKeep=frameKeep&Q.frameTable.absFrame>=100& ...
            Q.frameTable.absFrame<=216;
    end
    leadDraw=double(Q.leadDrawM(:,frameKeep));
    rangeDraw=double(Q.headLureDistanceDrawM(:,frameKeep));
    ratioDraw=leadDraw./rangeDraw;
    runMedianDraw=median(ratioDraw,2,'omitnan');
    mcMedianRatio(j)=percentile(runMedianDraw,50);
    mcLower95Ratio(j)=percentile(runMedianDraw,2.5);
    mcUpper95Ratio(j)=percentile(runMedianDraw,97.5);
end

summaryTable=table(runNames,runLabels,nFrame,q25Ratio,medianRatio,q75Ratio, ...
    mcMedianRatio,mcLower95Ratio,mcUpper95Ratio, ...
    'VariableNames',{'dataset','displayLabel','nFrames','nominalQ25_dOverR', ...
    'nominalMedian_dOverR','nominalQ75_dOverR','mcMedian_dOverR', ...
    'mcLower95_dOverR','mcUpper95_dOverR'});
writetable(summaryTable,fullfile(opts.outputDir, ...
    'spatial_lead_ratio_summary_by_run.csv'));

f=drawFigure(ratioByRun,summaryTable,colors,opts);
base=fullfile(opts.outputDir,'supplementary_spatial_lead_ratio');
exportgraphics(f,[base '.pdf'],'ContentType','vector');
print(f,[base '.png'],'-dpng',['-r' num2str(opts.figureResolution)]);
if ~opts.showFigure,close(f);end

writeCaption([base '_CAPTION.txt'],summaryTable);
S=struct('summaryTable',summaryTable,'ratioByRun',{ratioByRun}, ...
    'opts',opts,'nominalMedianRange',[min(medianRatio),max(medianRatio)]);
save(fullfile(opts.outputDir,'supplementary_spatial_lead_ratio.mat'),'S');

fprintf(['Spatial lead ratio: nominal run medians %.3f--%.3f; ' ...
    'saved supplementary figure to %s\n'], ...
    min(medianRatio),max(medianRatio),opts.outputDir);
end

function f=drawFigure(ratioByRun,S,colors,opts)
f=figure('Color','w','Visible',onoff(opts.showFigure), ...
    'Position',[100 100 1120 650]);
ax=axes(f);hold(ax,'on');
yline(ax,0,':','Color',[.35 .35 .35],'LineWidth',.9, ...
    'HandleVisibility','off');

nRun=height(S);
for j=1:nRun
    ratio=ratioByRun{j};
    n=numel(ratio);
    % Deterministic low-discrepancy jitter keeps every frame visible without
    % making the exported figure depend on the global random-number state.
    jitter=.18*(2*mod((1:n)'*0.618033988749895,1)-1);
    scatter(ax,j+jitter,ratio,28,mixWhite(colors(j,:),.34),'filled', ...
        'MarkerEdgeColor','w','LineWidth',.35,'HandleVisibility','off');

    % Nominal interquartile range and median.
    plot(ax,[j j],[S.nominalQ25_dOverR(j),S.nominalQ75_dOverR(j)],'-', ...
        'Color',colors(j,:),'LineWidth',6,'HandleVisibility','off');

    % Complete-draw 95% reconstruction interval for the run median.
    errorbar(ax,j,S.mcMedian_dOverR(j), ...
        S.mcMedian_dOverR(j)-S.mcLower95_dOverR(j), ...
        S.mcUpper95_dOverR(j)-S.mcMedian_dOverR(j), ...
        'Color',[.18 .18 .18],'LineStyle','none','LineWidth',1.15, ...
        'CapSize',10,'HandleVisibility','off');

    plot(ax,j,S.nominalMedian_dOverR(j),'d','Color',colors(j,:), ...
        'MarkerFaceColor','w','MarkerSize',8,'LineWidth',1.8, ...
        'HandleVisibility','off');
end

pointHandle=scatter(ax,nan,nan,28,[.55 .55 .55],'filled', ...
    'MarkerEdgeColor','w','DisplayName','Frames (nominal reconstruction)');
iqrHandle=plot(ax,nan,nan,'-','Color',[.35 .35 .35],'LineWidth',6, ...
    'DisplayName','Nominal interquartile range');
medianHandle=plot(ax,nan,nan,'d','Color',[.25 .25 .25], ...
    'MarkerFaceColor','w','MarkerSize',8,'LineWidth',1.8, ...
    'DisplayName','Nominal median');
intervalHandle=errorbar(ax,nan,nan,nan,nan,'Color',[.18 .18 .18], ...
    'LineStyle','none','LineWidth',1.15,'CapSize',10, ...
    'DisplayName','95% reconstruction interval of median');

xlim(ax,[.55,nRun+.55]);
allRatio=vertcat(ratioByRun{:});
yMin=min([allRatio;S.mcLower95_dOverR]);
yMax=max([allRatio;S.mcUpper95_dOverR]);
yPad=.08*(yMax-yMin);
ylim(ax,[yMin-yPad,yMax+yPad]);
numberY=yMin-.045*(yMax-yMin);
for j=1:nRun
    text(ax,j,numberY,sprintf('%.2f',S.nominalMedian_dOverR(j)), ...
        'Color',colors(j,:),'FontSize',10,'FontWeight','bold', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'BackgroundColor','w','Margin',1);
end
xticks(ax,1:nRun);xticklabels(ax,compose('Run %d',1:nRun));
ylabel(ax,'Spatial lead relative to head-to-target distance, d/R', ...
    'Interpreter','tex');
title(ax,'Spatial expression of predictive head orientation', ...
    'FontWeight','bold');
subtitle(ax,['Numbers at the bottom show nominal medians; ' ...
    'Phantom Run is restricted to frames 100--216'], ...
    'FontWeight','normal');
legend(ax,[pointHandle,iqrHandle,medianHandle,intervalHandle], ...
    'Location','northwest','Box','off','FontSize',9);
grid(ax,'on');box(ax,'off');
set(ax,'FontName','Arial','FontSize',11,'LineWidth',.85, ...
    'TickDir','out','GridColor',[.84 .84 .84],'GridAlpha',.45, ...
    'Layer','top');
end

function writeCaption(path,S)
fid=fopen(path,'w');
assert(fid>=0,'Could not write caption: %s',path);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,['Supplementary Figure | Spatial expression of predictive head ' ...
    'orientation. Each point is the ratio between signed spatial lead, d, and ' ...
    'current head-origin-to-target distance, R, for one analysed frame in the ' ...
    'nominal reconstruction. Thick coloured lines show nominal interquartile ' ...
    'ranges, diamonds show nominal per-pursuit medians, and capped black lines ' ...
    'show 95%% reconstruction intervals for the per-pursuit median propagated ' ...
    'through 500 complete joint pose/target reconstruction draws. Phantom Run ' ...
    'was restricted to absolute frames 100--216, matching the main model ' ...
    'comparison. Runs 1--5 are Jules Flick 1, Jules Top Run 11, Jules Top ' ...
    'Run 1, Lily Run, and Phantom Run, respectively. ' ...
    'Nominal median d/R ranged from %.2f to %.2f across pursuits.\n'], ...
    min(S.nominalMedian_dOverR),max(S.nominalMedian_dOverR));
end

function p=percentile(x,q)
x=sort(x(isfinite(x)));
if isempty(x),p=NaN;return;end
u=1+(numel(x)-1)*q/100;
lo=floor(u);hi=ceil(u);
if lo==hi,p=x(lo);else,p=x(lo)+(u-lo)*(x(hi)-x(lo));end
end

function c=mixWhite(c,a)
c=c+a*(1-c);
end

function s=onoff(tf)
if tf,s='on';else,s='off';end
end

function opts=applyDefaults(opts,defaults)
names=fieldnames(defaults);
for i=1:numel(names)
    if ~isfield(opts,names{i})||isempty(opts.(names{i}))
        opts.(names{i})=defaults.(names{i});
    end
end
end

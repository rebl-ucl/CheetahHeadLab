function R = yawSpatialLookAheadMonteCarloAnalysis(opts)
% yawSpatialLookAheadMonteCarloAnalysis  Reconstruction uncertainty in d(t).
%
% Each successful joint pose/lure Monte Carlo draw is propagated through
% the same horizontal head-ray/full-lure-path intersection used by the
% nominal spatial look-ahead analysis. Outputs are written beneath
% outputs/yaw/spatial/MonteCarloResults. Only the two Figure S15 examples
% are plotted. All five lateral pursuits remain in the numerical analysis.

arguments
    opts struct = struct()
end

thisDir=fileparts(mfilename('fullpath'));
pursuitDir=fileparts(thisDir);
repoRoot=fileparts(thisDir);
opts=applyDefaults(opts,struct( ...
    'dataRoot',fullfile(repoRoot,'outputs','analysis_inputs','yaw_tracking'), ...
    'nominalResultDir',fullfile(repoRoot,'outputs','yaw','spatial'), ...
    'outputDir',fullfile(repoRoot,'outputs','yaw','spatial','MonteCarloResults'), ...
    'maxMC',500,'intervalPercent',[2.5 97.5], ...
    'extendEndpoint',true,'endpointWindowSec',0.100,'plotOnly',false, ...
    'lowIntersectionSupport',0.95,'showFigures',false, ...
    'figureResolution',600,'individualTimeEndExclusive',struct('julesTopRun11_20170829',.64), ...
    'plotCombined',false, ...
    'individualDatasets',["julesTopRun1_20171209","julesTopRun11_20170829"]));

% Enforce manuscript membership for direct calls as well as the wrapper.
allowedDatasets=["julesTopRun1_20171209","julesTopRun11_20170829"];
assert(~opts.plotCombined,'release:nonManuscriptPlot', ...
    'The combined five-pursuit spatial time-series plot is not in the manuscript.');
if isempty(opts.individualDatasets), opts.individualDatasets=allowedDatasets; end
opts.individualDatasets=string(opts.individualDatasets);
assert(all(ismember(opts.individualDatasets,allowedDatasets)), ...
    'release:nonManuscriptPlot','Spatial time-series exports are restricted to Figure S15b,c.');

assert(numel(opts.intervalPercent)==2 && opts.intervalPercent(1)>=0 && ...
    opts.intervalPercent(2)<=100 && diff(opts.intervalPercent)>0, ...
    'intervalPercent must be an increasing two-element percentile range.');
assert(opts.maxMC>=1 && opts.lowIntersectionSupport>=0 && ...
    opts.lowIntersectionSupport<=1,'Invalid Monte Carlo options.');
if ~exist(opts.outputDir,'dir'), mkdir(opts.outputDir); end

if opts.plotOnly
    saved=load(fullfile(opts.outputDir,'yaw_spatial_lookahead_monte_carlo_analysis.mat'),'R');
    R=saved.R;
    assert(R.opts.extendEndpoint==opts.extendEndpoint && ...
        R.opts.endpointWindowSec==opts.endpointWindowSec,'Plot options do not match saved analysis.');
    plotFigures(R.runs,R.configuration,opts);
    return
end

C=steadyLateralConfiguration();
nominalPath=fullfile(opts.nominalResultDir,'yaw_spatial_lookahead_by_frame.csv');
assert(isfile(nominalPath),'Run the nominal spatial look-ahead analysis first: %s',nominalPath);
nominalAll=readtable(nominalPath,'TextType','string');

fprintf('\n==== yawSpatialLookAheadMonteCarloAnalysis ====\n');
fprintf('Output directory: %s\n',opts.outputDir);
runs=cell(height(C),1); frameRows=cell(height(C),1); summaryRows=cell(height(C),1);
for i=1:height(C)
    dataset=char(C.dataset(i)); animal=C.animal(i);
    fprintf('  %-30s ',dataset);
    runs{i}=propagateRun(dataset,animal,nominalAll,opts);
    frameRows{i}=runs{i}.frameTable;
    summaryRows{i}=runs{i}.summaryTable;
    fprintf('%d draws; median CI width %.3f m; median support %.3f\n', ...
        runs{i}.nDraw,runs{i}.summaryTable.medianLeadCIWidthM, ...
        runs{i}.summaryTable.medianIntersectionSupport);
end

frameTable=vertcat(frameRows{:});
summaryTable=vertcat(summaryRows{:});
writetable(frameTable,fullfile(opts.outputDir, ...
    'yaw_spatial_lookahead_monte_carlo_by_frame.csv'));
writetable(summaryTable,fullfile(opts.outputDir, ...
    'yaw_spatial_lookahead_monte_carlo_summary_by_run.csv'));
plotFigures(runs,C,opts);
writeReadme(C,opts);

R=struct('configuration',C,'runs',{runs},'frameTable',frameTable, ...
    'summaryTable',summaryTable,'opts',opts);
save(fullfile(opts.outputDir,'yaw_spatial_lookahead_monte_carlo_analysis.mat'), ...
    'R','-v7.3');
fprintf('Saved Monte Carlo spatial look-ahead analysis.\n');
end

function C=steadyLateralConfiguration()
dataset=["julesFlick1_20190309";"julesTopRun11_20170829"; ...
    "julesTopRun1_20171209";"lilyRun_20190309";"phantomRun_20190303"];
animal=["Jules";"Jules";"Jules";"Lily";"Phantom"];
C=table(dataset,animal);
end

function Q=propagateRun(dataset,animal,nominalAll,opts)
folder=fullfile(opts.dataRoot,dataset);
mcPath=fullfile(folder,'MonteCarloResults',['mcMain_' dataset '.mat']);
assert(isfile(mcPath),'Missing Monte Carlo file: %s',mcPath);
S=load(mcPath,'poseAll','lureXYZAll','nominalLureXYZ','sampleSucceeded', ...
    'absFrames','opts','provenance');
required={'poseAll','lureXYZAll','nominalLureXYZ','sampleSucceeded','absFrames'};
for i=1:numel(required), assert(isfield(S,required{i}),'%s lacks %s.',mcPath,required{i}); end
if isfield(S,'opts') && isfield(S.opts,'reconstructionMode')
    assert(strcmpi(char(S.opts.reconstructionMode),'trajectory'), ...
        '%s is not a trajectory-reconstruction posterior.',dataset);
end

ok=find(logical(S.sampleSucceeded(:)));
nDraw=min(numel(ok),floor(opts.maxMC)); ok=ok(1:nDraw);
assert(nDraw>0,'No successful Monte Carlo draws for %s.',dataset);
absFrame=double(S.absFrames(:)); N=numel(absFrame);
assert(size(S.poseAll,2)==N && size(S.poseAll,3)>=6, ...
    'Unexpected poseAll dimensions for %s.',dataset);
assert(size(S.lureXYZAll,2)==N && size(S.lureXYZAll,3)==3, ...
    'Unexpected lureXYZAll dimensions for %s.',dataset);

resultPath=fullfile(folder,'lureAnalysis_results_trajectory.mat');
A=load(resultPath,'runResults'); A=A.runResults;
assert(isequal(double(A.commonFrames(:)),absFrame), ...
    'Monte Carlo and nominal frames differ for %s.',dataset);
if isfield(A,'fs') && isfinite(A.fs(1)) && A.fs(1)>0, fs=double(A.fs(1)); else, fs=90; end
if isfield(A,'yawAnalysisMask') && numel(A.yawAnalysisMask)==N
    analysisMask=logical(A.yawAnalysisMask(:));
else
    analysisMask=true(N,1);
end
timeSec=(absFrame-absFrame(1))/fs;

[fullFrame,fullNominalXYZ]=loadFullLureTrajectory(folder);
[hasCurrent,currentLoc]=ismember(absFrame,fullFrame);
assert(all(hasCurrent),'Full lure trajectory does not cover all MC frames for %s.',dataset);

nominal=nominalAll(nominalAll.dataset==string(dataset),:);
[tf,loc]=ismember(absFrame,double(nominal.absFrame));
assert(all(tf),'Nominal spatial result does not cover all MC frames for %s.',dataset);
nominal=nominal(loc,:);

leadDraw=nan(nDraw,N,'single'); rangeDraw=nan(nDraw,N,'single');
extensionDraw=false(nDraw,N);
for m=1:nDraw
    pose=double(squeeze(S.poseAll(ok(m),:,:)));
    lure=double(squeeze(S.lureXYZAll(ok(m),:,:)));
    drawPath=deformFullPath(fullFrame,fullNominalXYZ,absFrame, ...
        double(S.nominalLureXYZ),lure,currentLoc);
    [lead,rangeM,extension]=spatialLookAheadForDraw(pose,lure,drawPath,currentLoc,analysisMask,fullFrame,fs,opts);
    extensionDraw(m,:)=extension;
    leadDraw(m,:)=single(lead);
    rangeDraw(m,:)=single(rangeM);
end

[leadMedian,leadLo,leadHi,intersectionSupport,probAhead]= ...
    summarisePosterior(leadDraw,opts.intervalPercent);
[rangeMedian,rangeLo,rangeHi]=summariseFiniteColumns(rangeDraw,opts.intervalPercent);
leadWidth=leadHi-leadLo;
nominalLead=double(nominal.leadDistanceM);
nominalRange=double(nominal.headOriginToLureDistanceM);
medianShift=leadMedian-nominalLead;

frameTable=table(repmat(string(dataset),N,1),repmat(string(animal),N,1), ...
    absFrame,timeSec,repmat(fs,N,1),analysisMask,nominalLead,leadMedian, ...
    leadLo,leadHi,leadWidth,medianShift,intersectionSupport,probAhead, ...
    nominalRange,rangeMedian,rangeLo,rangeHi, ...
    'VariableNames',{'dataset','animal','absFrame','timeSec','fps', ...
    'analysisMask','nominalLeadDistanceM','mcMedianLeadDistanceM', ...
    'mcLeadLowerM','mcLeadUpperM','mcLeadIntervalWidthM', ...
    'mcMedianMinusNominalM','mcIntersectionSupport','mcProbabilityAhead', ...
    'nominalHeadLureDistanceM','mcMedianHeadLureDistanceM', ...
    'mcHeadLureLowerM','mcHeadLureUpperM'});

frameTable.mcEndpointExtensionFraction=mean(extensionDraw,1)';
frameTable.nominalUsesEndpointExtension=logical(nominal.intersectionUsesEndpointExtension);
valid=analysisMask & isfinite(leadMedian);
nominalIQR=percentile(nominalLead(valid),75)-percentile(nominalLead(valid),25);
inside=valid & isfinite(nominalLead) & nominalLead>=leadLo & nominalLead<=leadHi;
signCertain=valid & (probAhead>=0.95 | probAhead<=0.05);
summaryTable=table(string(dataset),string(animal),nDraw,N,nnz(valid), ...
    median(intersectionSupport(analysisMask),'omitnan'), ...
    mean(intersectionSupport(analysisMask)>=opts.lowIntersectionSupport), ...
    median(leadWidth(valid),'omitnan'),percentile(leadWidth(valid),75), ...
    median(abs(medianShift(valid)),'omitnan'), ...
    median(leadWidth(valid),'omitnan')/max(nominalIQR,eps), ...
    nnz(inside)/max(1,nnz(valid & isfinite(nominalLead))), ...
    nnz(signCertain)/max(1,nnz(valid)), ...
    'VariableNames',{'dataset','animal','nDraw','nFrames','nFramesWithMCIntersection', ...
    'medianIntersectionSupport','fractionFramesIntersectionSupport95', ...
    'medianLeadCIWidthM','q75LeadCIWidthM','medianAbsoluteMedianShiftM', ...
    'medianCIWidthToNominalIQR','fractionNominalInsideMCInterval', ...
    'fractionFramesWith95SignCertainty'});

Q=struct('dataset',string(dataset),'animal',string(animal),'fs',fs, ...
    'nDraw',nDraw,'frameTable',frameTable,'summaryTable',summaryTable, ...
    'leadDrawM',leadDraw,'headLureDistanceDrawM',rangeDraw,'endpointExtensionDraw',extensionDraw);
end

function pathXYZ=deformFullPath(fullFrame,fullNominalXYZ,mcFrame,nominalMC,lureDraw,currentLoc)
% MC lure samples exist only on analysed frames. Apply their displacement
% from the MC nominal to the full nominal path; outside the sampled range,
% hold the nearest endpoint displacement constant rather than extrapolating.
delta=lureDraw-nominalMC;
pathXYZ=fullNominalXYZ;
for j=1:3
    good=isfinite(mcFrame)&isfinite(delta(:,j));
    if nnz(good)==0, continue; end
    x=mcFrame(good); y=delta(good,j);
    if numel(x)==1
        displacement=repmat(y,size(fullFrame));
    else
        displacement=interp1(x,y,fullFrame,'linear',NaN);
        displacement(fullFrame<x(1))=y(1);
        displacement(fullFrame>x(end))=y(end);
    end
    pathXYZ(:,j)=pathXYZ(:,j)+displacement;
end
% Guarantee exact agreement with the saved MC lure at analysed frames.
pathXYZ(currentLoc,:)=lureDraw;
end

function [lead,rangeM,extension]=spatialLookAheadForDraw(pose,lure,pathXYZ,currentLoc,analysisMask,pathFrame,fs,opts)
N=size(pose,1); lead=nan(N,1);extension=false(1,N);
u=spatialTerminalDirection(pathFrame,pathXYZ,fs,opts.endpointWindowSec);
eye=pose(:,4:6);
x=pose(:,1); y=pose(:,2); z=pose(:,3);
hx=cos(z).*cos(y)-sin(z).*sin(x).*sin(y);
hy=-sin(z).*cos(y)-cos(z).*sin(x).*sin(y);
headYaw=atan2(hy,hx);
rangeM=sqrt(sum((lure-eye).^2,2));
segmentLength=sqrt(sum(diff(pathXYZ,1,1).^2,2));
pathS=[0;cumsum(segmentLength)];
currentS=pathS(currentLoc);
valid=analysisMask & all(isfinite(eye),2) & all(isfinite(lure),2) & ...
    isfinite(headYaw) & isfinite(currentS) & isfinite(rangeM);
for k=find(valid).'
    lead(k)=firstForwardLead(eye(k,1:2),[cos(headYaw(k)),sin(headYaw(k))], ...
        pathXYZ,pathS,currentS(k));
    if ~isfinite(lead(k)) && opts.extendEndpoint
        lead(k)=spatialForwardExtension(eye(k,1:2),[cos(headYaw(k)),sin(headYaw(k))], ...
            pathXYZ(end,:),u,pathS(end),currentS(k));
        extension(k)=isfinite(lead(k));
    end
end
end

function lead=firstForwardLead(origin,direction,xyz,pathS,currentS)
lead=NaN; d=direction/norm(direction);
A=xyz(1:end-1,1:2); E=diff(xyz(:,1:2),1,1); Q=A-origin;
den=d(1)*E(:,2)-d(2)*E(:,1);
crossQD=Q(:,1)*d(2)-Q(:,2)*d(1);
scaleTol=1e-10*max([ones(size(E,1),1),sqrt(sum(E.^2,2)),sqrt(sum(Q.^2,2))],[],2);
rho=nan(size(den)); alpha=nan(size(den));
nonparallel=abs(den)>scaleTol;
rho(nonparallel)=(Q(nonparallel,1).*E(nonparallel,2)-Q(nonparallel,2).*E(nonparallel,1))./den(nonparallel);
alpha(nonparallel)=crossQD(nonparallel)./den(nonparallel);

parallel=~nonparallel & abs(crossQD)<=scaleTol;
for j=find(parallel).'
    rhoA=dot(A(j,:)-origin,d); rhoB=dot(A(j,:)+E(j,:)-origin,d);
    if max(rhoA,rhoB)<-scaleTol(j), continue; end
    if min(rhoA,rhoB)>=0
        if rhoA<=rhoB, rho(j)=rhoA; alpha(j)=0; else, rho(j)=rhoB; alpha(j)=1; end
    elseif abs(rhoB-rhoA)>eps
        rho(j)=0; alpha(j)=-rhoA/(rhoB-rhoA);
    end
end

valid=isfinite(rho)&isfinite(alpha)&rho>=-scaleTol&alpha>=-scaleTol&alpha<=1+scaleTol;
if ~any(valid), return; end
rho=max(rho,0); alpha=min(max(alpha,0),1);
candidate=find(valid);
candidateS=pathS(candidate)+alpha(candidate).*(pathS(candidate+1)-pathS(candidate));
minRho=min(rho(candidate)); tied=find(abs(rho(candidate)-minRho)<=1e-8);
if numel(tied)>1
    [~,jj]=min(abs(candidateS(tied)-currentS)); chosen=tied(jj);
else
    chosen=tied;
end
lead=candidateS(chosen)-currentS;
end

function [med,lo,hi,support,probAhead]=summarisePosterior(X,intervalPercent)
[med,lo,hi]=summariseFiniteColumns(X,intervalPercent);
n=size(X,1); support=sum(isfinite(X),1)'/n; probAhead=nan(size(support));
for k=1:size(X,2)
    v=double(X(:,k)); v=v(isfinite(v));
    if ~isempty(v), probAhead(k)=mean(v>0); end
end
end

function [med,lo,hi]=summariseFiniteColumns(X,intervalPercent)
N=size(X,2); med=nan(N,1); lo=nan(N,1); hi=nan(N,1);
for k=1:N
    v=double(X(:,k)); v=v(isfinite(v));
    if isempty(v), continue; end
    med(k)=percentile(v,50); lo(k)=percentile(v,intervalPercent(1));
    hi(k)=percentile(v,intervalPercent(2));
end
end

function plotFigures(Q,C,opts)
n=height(C);
for i=1:n
    if ~ismember(C.dataset(i),opts.individualDatasets), continue; end
    f=figure('Color','w','Visible',onoff(opts.showFigures),'Position',[100 100 1150 470]);
    ax=axes(f); plotRun(ax,Q{i},opts,true,false);
    base=fullfile(opts.outputDir,[char(C.dataset(i)) '_spatial_lookahead_monte_carlo_timeseries']);
    exportFigure(f,base,opts.figureResolution); close(f);
end
end

function plotRun(ax,Q,opts,showLegend,isCombined)
leadColor=[0 114 178]/255; rangeColor=[213 94 0]/255;
F=Q.frameTable;
% Display-only crop: never changes saved frames, draws or summary statistics.
timeEnd=Inf;
if ~isCombined && isfield(opts.individualTimeEndExclusive,char(Q.dataset))
    timeEnd=opts.individualTimeEndExclusive.(char(Q.dataset));
    assert(isscalar(timeEnd) && isfinite(timeEnd),'Invalid display time endpoint.');
    F=F(F.timeSec<timeEnd,:);
    assert(height(F)>=2,'Display window must contain at least two frames.');
end
x=F.timeSec; hold(ax,'on');

yyaxis(ax,'right'); hold(ax,'on');
hRangeBand=shadeBand(ax,x,F.mcHeadLureLowerM,F.mcHeadLureUpperM,rangeColor,.09, ...
    '95% MC head-lure interval');
hRange=plot(ax,x,F.mcMedianHeadLureDistanceM,'-','Color',rangeColor,'LineWidth',1.45, ...
    'DisplayName','MC median head-lure distance');
finiteRange=[F.mcHeadLureLowerM;F.mcHeadLureUpperM]; finiteRange=finiteRange(isfinite(finiteRange));
if ~isempty(finiteRange)
    span=range(finiteRange); pad=max(.05*span,.1);
    ylim(ax,[max(0,min(finiteRange)-pad),max(finiteRange)+pad]);
end
ylabel(ax,'Head-origin-to-lure distance (m)');

yyaxis(ax,'left'); hold(ax,'on');
yline(ax,0,':','Color',[.35 .35 .35],'LineWidth',.8,'HandleVisibility','off');
hLeadBand=shadeBand(ax,x,F.mcLeadLowerM,F.mcLeadUpperM,leadColor,.16,'95% MC lead interval');
extended=F.nominalUsesEndpointExtension | F.mcEndpointExtensionFraction>=.5;
observedMedian=F.mcMedianLeadDistanceM;observedMedian(extended)=NaN;
extendedMedian=F.mcMedianLeadDistanceM;extendedMedian(~extended)=NaN;
hMedian=plot(ax,x,observedMedian,'-','Color',leadColor,'LineWidth',1.15, ...
    'DisplayName','MC median: recorded path');
hExtended=plot(ax,x,extendedMedian,'--','Color',leadColor,'LineWidth',1.3, ...
    'DisplayName','MC median: endpoint extension');
hRaw=scatter(ax,x,F.nominalLeadDistanceM,isCombined*8+(~isCombined)*28, ...
    mixWithWhite(leadColor,.03),'filled','MarkerEdgeColor','white','LineWidth',.35, ...
    'DisplayName','Nominal reconstruction');
low=F.mcIntersectionSupport<opts.lowIntersectionSupport & isfinite(F.mcMedianLeadDistanceM);
hLow=plot(ax,x(low),F.mcMedianLeadDistanceM(low),'x','Color',[.35 .35 .35], ...
    'MarkerSize',4,'LineWidth',.7,'DisplayName','<95% intersection support');
local=[F.nominalLeadDistanceM;F.mcLeadLowerM;F.mcLeadUpperM]; local=local(isfinite(local));
if isempty(local), limits=[-1 1]; else
    lo=min(0,min(local)); hi=max(.5,max(local)); pad=.04*(hi-lo);
    limits=[lo-pad hi+pad];
end
if ~isCombined && Q.dataset=="julesTopRun11_20170829"
    limits=[0 40];
    yticks(ax,0:5:40);
    set(hLeadBand,'DisplayName','95% MC lead interval (clipped at 0–40 m)');
end
ylim(ax,limits);
if isfinite(timeEnd), xlim(ax,[x(1),timeEnd]); else, xlim(ax,[x(1),x(end)]); end
ylabel(ax,'Signed look-ahead, d (m)');
if ~isCombined, xlabel(ax,'Time from first reconstruction frame (s)'); end
S=Q.summaryTable;
if isCombined
    title(ax,sprintf('%s  |  median 95%% width = %.2f m  |  median intersection support = %.2f', ...
        char(Q.dataset),S.medianLeadCIWidthM,S.medianIntersectionSupport), ...
        'Interpreter','none','FontWeight','normal');
else
    title(ax,sprintf('%s: Monte Carlo uncertainty in spatial look-ahead', ...
        char(Q.dataset)),'Interpreter','none','FontWeight','normal');
end
grid(ax,'on'); box(ax,'off');
set(ax,'FontName','Arial','FontSize',9+~isCombined,'LineWidth',.75,'TickDir','out', ...
    'GridColor',[.82 .82 .82],'GridAlpha',.45,'Layer','top');
% Explicit title size after axes styling avoids inherited title scaling.
set(ax.Title,'FontUnits','points','FontSize',10+2*~isCombined, ...
    'FontName','Arial','FontWeight','normal','HorizontalAlignment','center');
ax.YAxis(1).Color=leadColor; ax.YAxis(2).Color=rangeColor;
if showLegend
    handles=[hRaw,hMedian,hExtended,hLeadBand,hLow,hRange,hRangeBand];
    handles=handles(isgraphics(handles));
    legendLocation='southwest';
    if ~isCombined && Q.dataset=="julesTopRun11_20170829"
        legendLocation='northwest';
    end
    lgd=legend(ax,handles,'Location',legendLocation,'NumColumns',3,'Box','on');
    lgd.Color='white'; lgd.EdgeColor=[.55 .55 .55]; lgd.FontSize=7+~isCombined;
end
end

function h=shadeBand(ax,x,lo,hi,color,alphaValue,label)
x=x(:); lo=lo(:); hi=hi(:); valid=isfinite(x)&isfinite(lo)&isfinite(hi);
h=gobjects(0); edges=diff([false;valid;false]); starts=find(edges==1); stops=find(edges==-1)-1;
for b=1:numel(starts)
    ix=starts(b):stops(b);
    if numel(ix)==1
        dx=max(eps,1e-4*max(1,range(x))); xx=[x(ix)-dx;x(ix)+dx];
        ll=[lo(ix);lo(ix)]; hh=[hi(ix);hi(ix)];
    else
        xx=x(ix); ll=lo(ix); hh=hi(ix);
    end
    p=fill(ax,[xx;flipud(xx)],[ll;flipud(hh)],color,'FaceAlpha',alphaValue, ...
        'EdgeColor','none');
    if isempty(h), p.DisplayName=label; h=p; else, p.HandleVisibility='off'; end
end
if isempty(h)
    h=fill(ax,NaN,NaN,color,'FaceAlpha',alphaValue,'EdgeColor','none','DisplayName',label);
end
end

function exportFigure(f,base,resolution)
pdfPath=[base '.pdf']; pngPath=[base '.png'];
centreExportTitles(f);
exportgraphics(f,pdfPath,'ContentType','vector');
cmd=sprintf('env -u LD_LIBRARY_PATH pdftoppm -f 1 -singlefile -png -r %g %s %s', ...
    resolution,shellQuote(pdfPath),shellQuote(base));
[status,message]=system(cmd);
if status~=0 || ~isfile(pngPath)
    warning('yawSpatialLookAheadMC:pdfRasterisationFailed', ...
        'pdftoppm failed (%s); using MATLAB print.',strtrim(message));
    print(f,pngPath,'-dpng',['-r' num2str(resolution)]);
end
end

function centreExportTitles(f)
% This headless renderer underestimates text Extent. Use glyph metrics to
% centre titles and horizontal/vertical axis labels before vector export.
[status,fontPath]=system('fc-match -f "%{file}" Arial');
assert(status==0 && isfile(strtrim(fontPath)),'Cannot locate title font.');
font=java.awt.Font.createFont(0,java.io.File(strtrim(fontPath)));
frc=java.awt.font.FontRenderContext([],true,true);
drawnow;
axesList=findall(f,'Type','axes');
for k=1:numel(axesList)
    ax=axesList(k);
    % Both yyaxis rulers have independent labels.
    labels={ax.Title,ax.XLabel};
    for j=1:numel(ax.YAxis), labels{end+1}=ax.YAxis(j).Label; end
    axPx=getpixelposition(ax);
    for j=1:numel(labels)
        t=labels{j};
        if isempty(t.String), continue; end
        set(t,'Interpreter','none','FontUnits','points');
        sizedFont=font.deriveFont(single(t.FontSize));
        bounds=sizedFont.getStringBounds(char(t.String),frc);
        widthPx=bounds.getWidth()*get(groot,'ScreenPixelsPerInch')/72;
        set(t,'Units','normalized'); pos=t.Position;
        if j<=2
            pos(1)=.5-widthPx/(2*axPx(3));
            if j==2
                % Fixed physical offset keeps the label close to the ticks.
                pos(2)=-36*get(groot,'ScreenPixelsPerInch')/72/axPx(4);
                set(t,'VerticalAlignment','top');
            end
        else
            % A left-anchored 90-degree label extends upward from pos(2).
            set(t,'Rotation',90);
            pos(2)=.5-widthPx/(2*axPx(4));
        end
        set(t,'Position',pos,'HorizontalAlignment','left');
    end
end
drawnow;
end

function [frame,xyz]=loadFullLureTrajectory(folder)
path=fullfile(folder,'3DLurePos_full.mat');
assert(isfile(path),'Missing full lure trajectory: %s',path);
S=load(path,'Lure3D'); L=S.Lure3D; frame=double(L.Index(:)); xyz=nan(numel(frame),3);
for i=1:numel(frame)
    mark=L.Mark{i}; if ~isempty(mark)&&numel(mark)>=3, xyz(i,:)=double(mark(1:3)); end
end
valid=isfinite(frame)&all(isfinite(xyz),2); frame=frame(valid); xyz=xyz(valid,:);
[frame,order]=sort(frame); xyz=xyz(order,:);
assert(numel(unique(frame))==numel(frame),'Duplicate full lure frames in %s.',path);
end

function writeReadme(C,opts)
path=fullfile(opts.outputDir,'README.txt'); fid=fopen(path,'w');
assert(fid>=0,'Could not write %s.',path); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'MONTE CARLO SPATIAL LOOK-AHEAD\n\n');
fprintf(fid,'Endpoint continuation enabled: %d; trailing direction window %.3f s.\n',opts.extendEndpoint,opts.endpointWindowSec);
fprintf(fid,'The terminal direction is re-estimated per draw; beyond-path lead assumes continued horizontal motion at constant height. No distance cap is imposed.\n');
fprintf(fid,['Every successful joint pose/lure draw from each mcMain file is propagated through ' ...
    'the horizontal head-ray/full-lure-polyline intersection.\n']);
fprintf(fid,'Intervals are %.1f--%.1f percentiles across draws with a valid intersection.\n', ...
    opts.intervalPercent(1),opts.intervalPercent(2));
fprintf(fid,['mcIntersectionSupport reports the fraction of all draws with an intersection; ' ...
    'intervals are conditional on intersection and should be read with this support.\n\n']);
fprintf(fid,['Saved lure draws cover analysed frames only. To retain the full route, each draw''s ' ...
    'lure displacement from its MC nominal is linearly interpolated onto 3DLurePos_full.mat. ' ...
    'The nearest endpoint displacement is held constant outside the sampled frame range, and ' ...
    'saved lure samples are imposed exactly at analysed frames.\n\n']);
fprintf(fid,'Included datasets (%d):\n',height(C));
for i=1:height(C), fprintf(fid,'  %s\n',C.dataset(i)); end
end

function p=percentile(x,q)
x=sort(x(isfinite(x))); if isempty(x), p=NaN; return; end
u=1+(numel(x)-1)*q/100; lo=floor(u); hi=ceil(u);
if lo==hi, p=x(lo); else, p=x(lo)+(u-lo)*(x(hi)-x(lo)); end
end

function c=mixWithWhite(c,amount), c=c+amount*(1-c); end
function s=onoff(tf), if tf, s='on'; else, s='off'; end, end
function q=shellQuote(s), q=['''' strrep(char(s),'''','''"''"''') '''']; end
function opts=applyDefaults(opts,defaults)
names=fieldnames(defaults);
for i=1:numel(names)
    if ~isfield(opts,names{i})||isempty(opts.(names{i})), opts.(names{i})=defaults.(names{i}); end
end
end

function F=plotHeadTorsoHandoff(R,opts)
%PLOTHEADTORSOHANDOFF  Render the two retained head-torso handoff figures.
%
%   plotHeadTorsoHandoff();      % reads outputs/yaw/head_torso/
%   plotHeadTorsoHandoff(R);     % or takes the returned struct directly
%
% Writes into outputs/yaw/head_torso/:
%   figure_6_closure_vs_decoupling          frame-level closure versus yaw separation
%   figure_7_top_view_decoupling_example    top-view reconstructed exemplar

arguments
    R=[]
    opts struct=struct()
end
thisDir=fileparts(mfilename('fullpath'));
repoRoot=fileparts(thisDir);
defaultDataRoot=fullfile(repoRoot,'outputs','analysis_inputs','yaw_tracking');
if ~isempty(R)&&isstruct(R)&&isfield(R,'opts')&&isfield(R.opts,'dataRoot')
    defaultDataRoot=R.opts.dataRoot;
end
opts=withDefaults(opts,struct('resultDir',fullfile(repoRoot,'outputs','yaw','head_torso'), ...
    'dataRoot',defaultDataRoot,'topViewExampleRun',"phantomRun_20190303", ...
    'showFigures',false,'figureResolution',600,'closureHullCoverage',.90));
assert(isfolder(opts.resultDir),'Run headTorsoHandoff first: %s',opts.resultDir);

runOrder=["julesFlick1_20190309";"julesTopRun11_20170829"; ...
    "julesTopRun1_20171209";"lilyRun_20190309";"phantomRun_20190303"];
labels=["Jules Flick 1";"Jules Top Run 11";"Jules Top Run 1"; ...
    "Lily Run";"Phantom Run"];
% Okabe-Ito colours, held fixed across the two figures.
runColors=[0 114 178;230 159 0;0 158 115;204 121 167;100 100 100]/255;
C=struct('head',[117 79 160]/255,'torso',[0 114 178]/255, ...
    'travel',[213 94 0]/255);

if isempty(R)
    BF=readtable(fullfile(opts.resultDir,'nominal_by_frame.csv'),'TextType','string');
else
    BF=R.byFrame;
end

F=struct();
F.closureVsDecoupling=plotClosureVsDecoupling(BF,runOrder,labels,runColors,opts);
F.topViewExample=plotTopViewExample(runOrder,C,opts);
end

% =========================================================================
function files=plotClosureVsDecoupling(BF,order,labels,runColors,opts)
% Does the cheetah keep closing even though head and torso yaw are separated?
% This preserves the frame-level (eta, c) pairing.  A curve or binned summary
% across eta would be misleading because the runs barely overlap in separation.
% Jules Top Run 11 is excluded a priori because its lure pulls away.
excludedRun="julesTopRun11_20170829";
keep=order~=excludedRun;
order=order(keep);labels=labels(keep);runColors=runColors(keep,:);
nR=numel(order);

cover=opts.closureHullCoverage;
assert(isfinite(cover)&&isscalar(cover)&&cover>0&&cover<=1, ...
    'closureHullCoverage must lie in (0, 1].');

valid=BF.valid==1&isfinite(BF.gazeOffsetDeg)&isfinite(BF.closingSpeedMps);
valid=valid&ismember(BF.run,order);
assert(any(valid),'No valid frames available for closure-vs-decoupling plot.');
eta=cell(nR,1);clos=cell(nR,1);
for j=1:nR
    sub=BF(BF.run==order(j)&valid,:);
    eta{j}=abs(sub.gazeOffsetDeg);   % |wrap(psiHead - psiTorso)|
    clos{j}=sub.closingSpeedMps;
end
etaAll=vertcat(eta{:});closAll=vertcat(clos{:});
nAll=numel(etaAll);
pctClosing=100*nnz(closAll>0)/nAll;

f=figure('Color','w','Visible',onoff(opts.showFigures), ...
    'Position',[40 60 900 560]);
ax=axes(f);hold(ax,'on');
xMax=10*ceil(max(etaAll)/10);
yLo=min(closAll)-.9;
yHi=ceil(linearQuantile(closAll,.99))+2;
patch(ax,[0 xMax xMax 0],[yLo yLo 0 0],[.93 .93 .93],'EdgeColor','none', ...
    'HandleVisibility','off');
yline(ax,0,'-','Color',[.15 .15 .15],'LineWidth',1.2,'HandleVisibility','off');

nClipped=0;
for j=1:nR
    inAxis=clos{j}<=yHi;
    scatter(ax,eta{j}(inAxis),clos{j}(inAxis),13,runColors(j,:),'filled', ...
        'MarkerFaceAlpha',.32,'MarkerEdgeColor','none','HandleVisibility','off');
    if any(~inAxis)
        scatter(ax,eta{j}(~inAxis),repmat(yHi-.15,nnz(~inAxis),1),26, ...
            runColors(j,:),'^','filled','MarkerEdgeColor','w','LineWidth',.6, ...
            'HandleVisibility','off');
        nClipped=nClipped+nnz(~inAxis);
    end
    P=coverageHull(eta{j},clos{j},cover);
    patch(ax,P(:,1),P(:,2),runColors(j,:),'FaceAlpha',.05, ...
        'EdgeColor',runColors(j,:),'LineWidth',1.6,'HandleVisibility','off');
    plot(ax,median(eta{j}),median(clos{j}),'d','MarkerSize',9, ...
        'MarkerFaceColor',runColors(j,:),'MarkerEdgeColor','w','LineWidth',1.2, ...
        'HandleVisibility','off');
end

% Direct labels retain run identity without a legend box over the scatter.
anchorRun=["julesFlick1_20190309";"julesTopRun1_20171209"; ...
    "lilyRun_20190309";"phantomRun_20190303"];
anchorXY=[19.5 .6 -1;52 6.7 1;4 6.4 1;76 -3.2 -1];
for j=1:nR
    a=anchorXY(anchorRun==order(j),:);
    if isempty(a),a=[median(eta{j}) median(clos{j})+1.2 0];end
    directLabel(ax,a(1),a(2),labels(j),runColors(j,:),a(3));
end
text(ax,.5,yLo+.35,'receding  (c < 0)','FontSize',9,'Color',[.42 .42 .42]);
if nClipped>0
    text(ax,xMax-.5,yHi-.45,sprintf('%d frames above axis (triangles)',nClipped), ...
        'FontSize',8,'Color',[.42 .42 .42],'HorizontalAlignment','right');
end
xlim(ax,[0 xMax]);ylim(ax,[yLo yHi]);
xlabel(ax,'head-torso yaw separation,  \eta = |wrap(\psi_H - \psi_T)|  (deg)');
ylabel(ax,'closing speed,  c = -dR/dt  (m/s)');
tidy(ax);
title(ax,{'Closing speed stays positive across the full range of head-torso yaw', ...
    sprintf('separation seen in these runs  (%.0f%% of %d valid frames closing)',pctClosing,nAll)}, ...
    'FontWeight','normal');

files=exportPair(f,fullfile(opts.resultDir,'figure_6_closure_vs_decoupling'),opts);
end

function directLabel(ax,x,y,txt,col,align)
switch align
    case -1,h='right';dx= 1.7;
    case  1,h='left'; dx=-1.7;
    otherwise,h='center';dx=0;
end
text(ax,x,y,txt,'FontSize',9.5,'Color',[.15 .15 .15],'HorizontalAlignment',h, ...
    'VerticalAlignment','middle','BackgroundColor','w','Margin',.5);
if dx~=0
    plot(ax,x+dx,y,'o','MarkerSize',6,'MarkerFaceColor',col, ...
        'MarkerEdgeColor','w','LineWidth',.8,'HandleVisibility','off');
end
end

function P=coverageHull(x,y,frac)
% Convex hull of the central frames: descriptive support, not a confidence set.
x=x(:);y=y(:);
sx=std(x);sy=std(y);
if sx<=0,sx=1;end
if sy<=0,sy=1;end
d=hypot((x-median(x))/sx,(y-median(y))/sy);
k=d<=linearQuantile(d,frac);
if nnz(k)<3,k=true(size(x));end
xk=x(k);yk=y(k);
if size(unique([xk yk],'rows'),1)<3,P=[xk yk];return;end
h=convhull(xk,yk,'Simplify',true);
P=[xk(h) yk(h)];
end

% =========================================================================
function files=plotTopViewExample(order,C,opts)
% One reconstructed closing pursuit in world X-Y.  Five interior snapshots show
% the torso axis tangent to locomotor travel while the head faces the lure.
exampleRun=string(opts.topViewExampleRun);
assert(ismember(exampleRun,order),'Unknown topViewExampleRun: %s',exampleRun);
folder=fullfile(opts.dataRoot,char(exampleRun));
mcPath=fullfile(folder,'MonteCarloResults',['mcMain_' char(exampleRun) '.mat']);
assert(isfile(mcPath),'Missing top-view reconstruction input: %s',mcPath);
M=load(mcPath,'nominalPose','nominalLureXYZ','absFrames');
A=load(fullfile(folder,'lureAnalysis_results_trajectory.mat'),'runResults');
pose=double(M.nominalPose);lure=double(M.nominalLureXYZ);
assert(size(pose,1)==size(lure,1)&&size(pose,2)>=12&&size(lure,2)>=2, ...
    'Unexpected pose/lure dimensions for %s.',exampleRun);
if isfield(A.runResults,'yawAnalysisMask')&& ...
        numel(A.runResults.yawAnalysisMask)==size(pose,1)
    valid=logical(A.runResults.yawAnalysisMask(:));
else
    valid=true(size(pose,1),1);
end

head=pose(:,4:5);torso=pose(:,10:11);lure=lure(:,1:2);
headForward=segmentForwardXY(pose(:,1),pose(:,2),pose(:,3));
torsoForward=segmentForwardXY(pose(:,7),pose(:,8),pose(:,9));
valid=valid&all(isfinite(head),2)&all(isfinite(torso),2)& ...
    all(isfinite(lure),2)&all(isfinite(headForward),2)& ...
    all(isfinite(torsoForward),2);
good=find(valid);
assert(numel(good)>=10,'%s has too few valid frames for a five-snapshot plot.',exampleRun);

snapshotPercent=linspace(15,85,5);
q=round(1+(snapshotPercent/100)*(numel(good)-1));
snapshotIx=good(q);
assert(numel(unique(snapshotIx))==5,'Snapshot selection produced duplicate frames.');

colHead=C.head;colTorso=C.torso;colLure=C.travel;
colHeadTrace=mixWhite(colHead,.67);
colLos=[.72 .72 .72];
xy=[head(good,:);torso(good,:);lure(good,:)];
span=max(max(xy)-min(xy));
axisLength=.085*span;
pad=.08*span;

f=figure('Color','w','Visible',onoff(opts.showFigures), ...
    'Position',[40 60 900 560]);
ax=axes(f);hold(ax,'on');
hTorsoTrace=plot(ax,torso(good,1),torso(good,2),'-','Color',colTorso, ...
    'LineWidth',2.6,'DisplayName','torso-origin trace');
hHeadTrace=plot(ax,head(good,1),head(good,2),'-','Color',colHeadTrace, ...
    'LineWidth',1.7,'DisplayName','head-origin trace');
scatter(ax,head(good,1),head(good,2),9,colHead,'filled', ...
    'MarkerFaceAlpha',.12,'MarkerEdgeAlpha',0,'HandleVisibility','off');
hLureTrace=plot(ax,lure(good,1),lure(good,2),'-','Color',colLure, ...
    'LineWidth',2.3,'DisplayName','lure trace');
plot(ax,torso(good(end),1),torso(good(end),2),'>','Color',colTorso, ...
    'MarkerFaceColor',colTorso,'MarkerSize',7,'HandleVisibility','off');
plot(ax,lure(good(end),1),lure(good(end),2),'>','Color',colLure, ...
    'MarkerFaceColor',colLure,'MarkerSize',7,'HandleVisibility','off');

hLos=gobjects(1);hTorsoAxis=gobjects(1);hHeadAxis=gobjects(1);
rangeM=hypot(lure(snapshotIx,1)-torso(snapshotIx,1), ...
    lure(snapshotIx,2)-torso(snapshotIx,2));
for k=1:numel(snapshotIx)
    i=snapshotIx(k);
    h=plot(ax,[head(i,1) lure(i,1)],[head(i,2) lure(i,2)],':', ...
        'Color',colLos,'LineWidth',1.1, ...
        'HandleVisibility',ternary(k==1,'on','off'), ...
        'DisplayName','simultaneous head-to-lure bearing');
    if k==1,hLos=h;end
    scatter(ax,torso(i,1),torso(i,2),52,colTorso,'s','filled', ...
        'MarkerEdgeColor','w','LineWidth',.8,'HandleVisibility','off');
    scatter(ax,head(i,1),head(i,2),48,colHead,'o','filled', ...
        'MarkerEdgeColor','w','LineWidth',.8,'HandleVisibility','off');
    scatter(ax,lure(i,1),lure(i,2),64,colLure,'o','filled', ...
        'MarkerEdgeColor','w','LineWidth',.9,'HandleVisibility','off');
    h=quiver(ax,torso(i,1),torso(i,2), ...
        axisLength*torsoForward(i,1),axisLength*torsoForward(i,2),0, ...
        'Color',colTorso,'LineWidth',2.2,'MaxHeadSize',.45, ...
        'HandleVisibility',ternary(k==1,'on','off'),'DisplayName','torso forward axis');
    if k==1,hTorsoAxis=h;end
    h=quiver(ax,head(i,1),head(i,2), ...
        axisLength*headForward(i,1),axisLength*headForward(i,2),0, ...
        'Color',colHead,'LineWidth',2.2,'MaxHeadSize',.45, ...
        'HandleVisibility',ternary(k==1,'on','off'),'DisplayName','head forward axis');
    if k==1,hHeadAxis=h;end
    text(ax,torso(i,1),torso(i,2)+.018*span,sprintf('T%d',k), ...
        'Color',mixWhite(colTorso,.08),'FontSize',8.5, ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'BackgroundColor','w','Margin',1);
    text(ax,head(i,1)-.014*span,head(i,2)-.014*span,sprintf('H%d',k), ...
        'Color',mixWhite(colHead,.08),'FontSize',8.5, ...
        'HorizontalAlignment','right','VerticalAlignment','top', ...
        'BackgroundColor','w','Margin',1);
    text(ax,lure(i,1),lure(i,2)-.025*span,sprintf('L%d',k), ...
        'Color',mixWhite(colLure,.08),'FontSize',8.5, ...
        'HorizontalAlignment','center','VerticalAlignment','top', ...
        'BackgroundColor','w','Margin',1);
end

axis(ax,'equal');
xlim(ax,[min(xy(:,1))-pad max(xy(:,1))+pad]);
ylim(ax,[min(xy(:,2))-pad max(xy(:,2))+pad]);
xlabel(ax,'world X (m)');ylabel(ax,'world Y (m)');
legend(ax,[hTorsoTrace hHeadTrace hLureTrace hTorsoAxis hHeadAxis hLos], ...
    'Location','northeast','NumColumns',2,'Box','on','Color','w', ...
    'EdgeColor',[.65 .65 .65],'FontSize',8);
tidy(ax);
grid(ax,'off');
title(ax,'Head-torso decoupling during closure','FontWeight','normal','FontSize',13);

absFrames=double(M.absFrames(:));
if isfield(A.runResults,'fs'),fs=double(A.runResults.fs(1));else,fs=NaN;end
headYaw=atan2(headForward(:,2),headForward(:,1));
torsoYaw=atan2(torsoForward(:,2),torsoForward(:,1));
headBearing=atan2(lure(:,2)-head(:,2),lure(:,1)-head(:,1));
snapshotTable=table((1:5)',snapshotPercent(:),snapshotIx(:), ...
    absFrames(snapshotIx),(absFrames(snapshotIx)-absFrames(good(1)))/fs, ...
    torso(snapshotIx,1),torso(snapshotIx,2),head(snapshotIx,1),head(snapshotIx,2), ...
    lure(snapshotIx,1),lure(snapshotIx,2),rangeM(:), ...
    rad2deg(abs(wrapPiLocal(headYaw(snapshotIx)-torsoYaw(snapshotIx)))), ...
    rad2deg(abs(wrapPiLocal(headYaw(snapshotIx)-headBearing(snapshotIx)))), ...
    'VariableNames',{'snapshot','validSequencePercent','relativeIndex', ...
    'absFrame','timeSec','torsoX','torsoY','headX','headY','lureX','lureY', ...
    'torsoLureRangeM','headTorsoSeparationDeg','headLureOffsetDeg'});
writetable(snapshotTable,fullfile(opts.resultDir, ...
    'figure_7_top_view_decoupling_example_snapshots.csv'));

files=exportPair(f,fullfile(opts.resultDir,'figure_7_top_view_decoupling_example'),opts);
end

% =========================================================================
function tidy(ax)
set(ax,'Box','off','TickDir','out','LineWidth',.9,'FontSize',10, ...
    'XColor',[.25 .25 .25],'YColor',[.25 .25 .25]);
grid(ax,'on');ax.GridColor=[.85 .85 .85];ax.GridAlpha=1;ax.Layer='bottom';
end

function q=linearQuantile(x,p)
x=sort(x(isfinite(x)));
assert(~isempty(x)&&isscalar(p)&&p>=0&&p<=1,'Invalid quantile input.');
z=1+(numel(x)-1)*p;
i=floor(z);j=ceil(z);
if i==j,q=x(i);else,q=x(i)+(z-i)*(x(j)-x(i));end
end

function xy=segmentForwardXY(tx,ty,tz)
x=cos(tz).*cos(ty)-sin(tz).*sin(tx).*sin(ty);
y=-sin(tz).*cos(ty)-cos(tz).*sin(tx).*sin(ty);
n=hypot(x,y);xy=[x./n,y./n];xy(n<eps,:)=NaN;
end

function y=wrapPiLocal(x),y=mod(x+pi,2*pi)-pi;end
function c=mixWhite(c,amount),c=(1-amount)*c+amount*[1 1 1];end
function v=ternary(c,a,b),if c,v=a;else,v=b;end,end

function files=exportPair(f,base,opts)
exportgraphics(f,[base '.png'],'Resolution',opts.figureResolution);
exportgraphics(f,[base '.pdf'],'ContentType','vector');
files={[base '.png'],[base '.pdf']};
if ~opts.showFigures,close(f);end
end

function s=onoff(tf),if tf,s='on';else,s='off';end,end

function opts=withDefaults(opts,defaults)
n=fieldnames(defaults);
for i=1:numel(n)
    if ~isfield(opts,n{i})||isempty(opts.(n{i})),opts.(n{i})=defaults.(n{i});end
end
end

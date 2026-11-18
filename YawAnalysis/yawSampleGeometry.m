function G = yawSampleGeometry(pose,lure,t,valid,fs,cfg)
    N = size(pose,1);
    assert(size(lure,1)==N && numel(t)==N,'Pose/lure/time length mismatch.');
    valid = logical(valid(:)) & all(isfinite(pose),2) & all(isfinite(lure),2);

    % Head-forward vector under the exact mcMainAngles rotation convention:
    % Rz(z)*Rx(x)*Ry(y)*[1;0;0].
    x = pose(:,1); y = pose(:,2); z = pose(:,3);
    hx = cos(z).*cos(y) - sin(z).*sin(x).*sin(y);
    hy = -sin(z).*cos(y) - cos(z).*sin(x).*sin(y);
    headYawRaw = atan2(hy,hx);
    eye = pose(:,4:6);
    bearingRaw = atan2(lure(:,2)-eye(:,2),lure(:,1)-eye(:,1));
    gammaRaw = wrapPi(headYawRaw-bearingRaw);

    gammaFilt = wrapPi(despikeSeries(unwrap(gammaRaw),cfg.hampelK,cfg.hampelSigma));
    bearingUnw = unwrap(bearingRaw);
    headYaw = bearingUnw + unwrap(gammaFilt);

    Q = localPathMetric(lure(:,1:2),eye(:,1:2),headYaw,t,valid,fs,cfg);
    G = Q;
    G.valid = valid & isfinite(gammaFilt) & isfinite(Q.pathLeadM);
    G.gammaDeg = rad2deg(gammaFilt);
    G.headYawRad = headYaw;
    G.bearingCurrentRad = bearingUnw;
    G.eyeXYZ = eye;
    G.lureXYZ = lure;
end

function Q = localPathMetric(lureXY,eyeXY,headYaw,t,valid,fs,cfg)
    N = size(lureXY,1);
    baseHalf = max(2,round(cfg.pathWindowSec*fs/2));
    maxHalf = max(baseHalf,round(cfg.pathMaxWindowSec*fs/2));
    tangent = nan(N,2); reliable = false(N,1); borrowed = false(N,1);

    for c = 1:N
        if ~valid(c), continue; end
        chosen = [NaN NaN];
        for h = baseHalf:maxHalf
            ix = max(1,c-h):min(N,c+h);
            keep = valid(ix) & isfinite(t(ix)) & all(isfinite(lureXY(ix,:)),2);
            ix = ix(keep);
            if numel(ix)<3, continue; end
            tt = t(ix)-t(c);
            X = [ones(numel(ix),1),tt(:)];
            bx = X\lureXY(ix,1); by = X\lureXY(ix,2);
            v = [bx(2),by(2)]; speed = norm(v);
            chosen = v/max(speed,eps);
            xx = lureXY(ix,1); yy = lureXY(ix,2);
            extent = hypot(max(xx)-min(xx),max(yy)-min(yy));
            C = [xx-mean(xx),yy-mean(yy)]; [~,Ss,~] = svd(C,0); ss=diag(Ss);
            if isempty(ss), linearity=NaN; else, linearity=ss(1)^2/max(sum(ss.^2),eps); end
            if extent>=cfg.pathMinDisplacementM && ...
                    linearity>=cfg.pathMinLinearity && speed>1e-9
                reliable(c)=true;
                break;
            end
        end
        tangent(c,:) = chosen;
    end

    good = find(reliable & all(isfinite(tangent),2));
    assert(~isempty(good),'No reliable lure direction in posterior draw.');
    for c = find(valid & ~reliable).'
        [~,j] = min(abs(t(good)-t(c)));
        tangent(c,:) = tangent(good(j),:); borrowed(c)=true;
    end

    range = hypot(lureXY(:,1)-eyeXY(:,1),lureXY(:,2)-eyeXY(:,2));
    gaze = eyeXY + range.*[cos(headYaw),sin(headYaw)];
    z = gaze-lureXY;
    lead = sum(z.*tangent,2);
    cross = -z(:,1).*tangent(:,2)+z(:,2).*tangent(:,1);
    lead(~valid)=NaN; cross(~valid)=NaN;
    Q = struct('pathLeadM',lead,'pathCrossM',cross,'pathTangent',tangent, ...
        'directionBorrowed',borrowed);
end

function y=despikeSeries(x,k,nSigma)
    x=x(:);y=x;fin=isfinite(x);if nnz(fin)<2*k+1,return;end
    xf=x(fin);done=false;
    if exist('hampel','file')==2 || exist('hampel','builtin')==2
        try,xf=hampel(xf,k,nSigma);done=true;catch,end
    end
    if ~done
        yf=xf;
        for i=1:numel(xf)
            w=xf(max(1,i-k):min(numel(xf),i+k));med=median(w);
            sig=1.4826*median(abs(w-med));
            if sig>0 && abs(xf(i)-med)>nSigma*sig,yf(i)=med;end
        end
        xf=yf;
    end
    y(fin)=xf;
end

function y=wrapPi(x),y=mod(x+pi,2*pi)-pi;end

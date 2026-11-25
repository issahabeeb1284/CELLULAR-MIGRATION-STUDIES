clear; close all; clc;

%% SETTINGS
dataFolder = pwd;
fileName = 'WT VS COSMC HYALURONIC ACID1.csv';

pixelSize_um = 0.365;
dt_hr = 5/60;
minFrames = 20;

maxAllowedJump_um = 50;
maxEuclidean_um = 300;
maxAccumulated_um = 1500;

groupOrder = ["U87 WT","COSMC-KO"];

outFolder = fullfile(dataFolder,'WT_COSMC_COMPLETE_MIGRATION_STATS');
if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

%% LOAD DATA
T = readtable(fullfile(dataFolder,fileName));
T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

T.Group = string(T.Group);
T.Group = strtrim(T.Group);

T.Group(T.Group=="U87-WT") = "U87 WT";
T.Group(T.Group=="WT") = "U87 WT";
T.Group(T.Group=="COSMC KO") = "COSMC-KO";
T.Group(T.Group=="COSMC_KO") = "COSMC-KO";
T.Group(T.Group=="KO") = "COSMC-KO";

AllMetrics = table();
AllCoords = table();
RemovedTracks = table();

%% ANALYZE TRACKS
for g = 1:length(groupOrder)

    groupName = groupOrder(g);
    Tg = T(T.Group == groupName,:);

    trackIDs = unique(Tg.Track_No);

    for i = 1:length(trackIDs)

        id = trackIDs(i);
        tr = Tg(Tg.Track_No == id,:);
        tr = sortrows(tr,'Slice_No');

        if height(tr) < minFrames
            RemovedTracks = [RemovedTracks; table(groupName,id,"Too few frames",height(tr),NaN,NaN,NaN, ...
                'VariableNames',{'Group','TrackID','Reason','N_Frames','MaxJump_um','EuclideanDistance_um','AccumulatedDistance_um'})];
            continue
        end

        frame = tr.Slice_No;

        x_um = tr.X * pixelSize_um;
        y_um = tr.Y * pixelSize_um;

        x0 = x_um - x_um(1);
        y0 = y_um - y_um(1);

        dx = diff(x_um);
        dy = diff(y_um);
        stepDist_um = sqrt(dx.^2 + dy.^2);

        accumulatedDistance_um = sum(stepDist_um);

        euclideanDistance_um = sqrt( ...
            (x_um(end)-x_um(1))^2 + ...
            (y_um(end)-y_um(1))^2 );

        totalTime_hr = (length(x_um)-1) * dt_hr;
        speed_um_hr = accumulatedDistance_um / totalTime_hr;

        if accumulatedDistance_um == 0
            migrationIndex = NaN;
        else
            migrationIndex = euclideanDistance_um / accumulatedDistance_um;
        end

        maxJump_um = max(stepDist_um);

        if maxJump_um > maxAllowedJump_um
            RemovedTracks = [RemovedTracks; table(groupName,id,"Jump too large",length(x_um),maxJump_um,euclideanDistance_um,accumulatedDistance_um, ...
                'VariableNames',{'Group','TrackID','Reason','N_Frames','MaxJump_um','EuclideanDistance_um','AccumulatedDistance_um'})];
            continue
        end

        if euclideanDistance_um > maxEuclidean_um
            RemovedTracks = [RemovedTracks; table(groupName,id,"Euclidean displacement too large",length(x_um),maxJump_um,euclideanDistance_um,accumulatedDistance_um, ...
                'VariableNames',{'Group','TrackID','Reason','N_Frames','MaxJump_um','EuclideanDistance_um','AccumulatedDistance_um'})];
            continue
        end

        if accumulatedDistance_um > maxAccumulated_um
            RemovedTracks = [RemovedTracks; table(groupName,id,"Accumulated distance too large",length(x_um),maxJump_um,euclideanDistance_um,accumulatedDistance_um, ...
                'VariableNames',{'Group','TrackID','Reason','N_Frames','MaxJump_um','EuclideanDistance_um','AccumulatedDistance_um'})];
            continue
        end

        motilityCoefficient = calculateMotilityCoefficient(x0,y0,dt_hr);
        persistenceTime_hr = calculatePersistenceTime(x_um,y_um,dt_hr);

        uniqueID = groupName + "_Track_" + string(id);

        coordTable = table();
        coordTable.Group = repmat(groupName,length(x_um),1);
        coordTable.TrackID = repmat(id,length(x_um),1);
        coordTable.UniqueTrackID = repmat(uniqueID,length(x_um),1);
        coordTable.Frame = frame;
        coordTable.Time_hr = (frame-min(frame))*dt_hr;
        coordTable.X_um = x_um;
        coordTable.Y_um = y_um;
        coordTable.X_centered_um = x0;
        coordTable.Y_centered_um = y0;
        coordTable.NetDisplacement_um = sqrt(x0.^2 + y0.^2);

        AllCoords = [AllCoords; coordTable];

        metricTable = table();
        metricTable.Group = groupName;
        metricTable.TrackID = id;
        metricTable.UniqueTrackID = uniqueID;
        metricTable.N_Frames = length(x_um);
        metricTable.TotalTime_hr = totalTime_hr;
        metricTable.Speed_um_hr = speed_um_hr;
        metricTable.AccumulatedDistance_um = accumulatedDistance_um;
        metricTable.EuclideanDistance_um = euclideanDistance_um;
        metricTable.MigrationIndex = migrationIndex;
        metricTable.MotilityCoefficient_um2_hr = motilityCoefficient;
        metricTable.PersistenceTime_hr = persistenceTime_hr;
        metricTable.MaxJump_um = maxJump_um;

        AllMetrics = [AllMetrics; metricTable];

    end
end

%% PANEL A: TRAJECTORIES WITH IDENTICAL AXES
WT_color = [0.00 0.3 0.70];
KO_color = [0.5 0.00 0.15];

% Give every track a unique, reproducible color across both panels.
% Golden-ratio hue spacing separates consecutive track colors.
allTrackIDs = unique(AllCoords.UniqueTrackID);
nTracks = numel(allTrackIDs);
hues = mod((0:nTracks-1)' * 0.618033988749895, 1);
trackColors = hsv2rgb([hues, 0.85*ones(nTracks,1), 0.75*ones(nTracks,1)]);

figure('Color','w','Units','inches','Position',[1 1 8 4]);

globalX = AllCoords.X_centered_um;
globalY = AllCoords.Y_centered_um;

maxRange = max(abs([globalX; globalY]));
maxRange = ceil(maxRange/25)*25;

for g = 1:length(groupOrder)

    subplot(1,2,g); hold on;

    groupName = groupOrder(g);

    tracks = unique(AllCoords.UniqueTrackID(AllCoords.Group==groupName));

    for i = 1:length(tracks)
        idx = AllCoords.UniqueTrackID==tracks(i);
        colorIndex = find(allTrackIDs == tracks(i),1);
        colorNow = trackColors(colorIndex,:);

        plot(AllCoords.X_centered_um(idx), ...
             AllCoords.Y_centered_um(idx), ...
             'Color',colorNow, ...
             'LineWidth',1.2);
    end

    xline(0,'--','Color',[0.8 0.8 0.8],'LineWidth',1.0);
    yline(0,'--','Color',[0.8 0.8 0.8],'LineWidth',1.0);

    axis equal
    xlim([-maxRange maxRange]);
    ylim([-maxRange maxRange]);

    box off
    xlabel('X-direction (\mum)');
    ylabel('Y-direction (\mum)');
    title(groupName,'FontName','Arial','FontSize',18,'FontWeight','bold');

    set(gca,'FontName','Arial','FontSize',18,'LineWidth',1.2,'TickDir','out');
end

saveas(gcf,fullfile(outFolder,'Panel_A_Trajectories_UnifiedAxes.png'));

%% STATISTICAL ANALYSIS
metricsToTest = {
    'Speed_um_hr',                 'Speed_um_per_hr';
    'AccumulatedDistance_um',      'AccumulatedDistance_um';
    'EuclideanDistance_um',        'EuclideanDisplacement_um';
    'MigrationIndex',              'MigrationIndex';
    'MotilityCoefficient_um2_hr',  'MotilityCoefficient_um2_per_hr';
    'PersistenceTime_hr',          'PersistenceTime_hr'
};

StatsTable = table();

for m = 1:size(metricsToTest,1)

    metricName = metricsToTest{m,1};
    cleanName  = metricsToTest{m,2};

    WT = AllMetrics.(metricName)(AllMetrics.Group=="U87 WT");
    KO = AllMetrics.(metricName)(AllMetrics.Group=="COSMC-KO");

    WT = WT(~isnan(WT));
    KO = KO(~isnan(KO));

    [p,tstat,df,ci] = welchTTest_manual(WT,KO);

    meanWT = mean(WT);
    meanKO = mean(KO);
    semWT = std(WT)/sqrt(length(WT));
    semKO = std(KO)/sqrt(length(KO));

    pooledSD = sqrt(((length(WT)-1)*var(WT) + (length(KO)-1)*var(KO)) / ...
        (length(WT)+length(KO)-2));

    CohenD = (meanWT - meanKO) / pooledSD;

    sig = significanceStars(p);

    row = table();
    row.Metric = string(cleanName);
    row.N_WT = length(WT);
    row.N_COSMC_KO = length(KO);
    row.WT_mean = meanWT;
    row.WT_sem = semWT;
    row.COSMC_KO_mean = meanKO;
    row.COSMC_KO_sem = semKO;
    row.Mean_Difference_WT_minus_KO = meanWT - meanKO;
    row.t_statistic = tstat;
    row.df = df;
    row.p_value = p;
    row.Significance = string(sig);
    row.CI95_lower = ci(1);
    row.CI95_upper = ci(2);
    row.CohenD = CohenD;

    StatsTable = [StatsTable; row];
end

%% PANELS B-G: BOXPLOTS WITH STATS
plotMetric(AllMetrics,'Speed_um_hr','Speed (\mum/h)',outFolder,'Panel_B_Speed');
plotMetric(AllMetrics,'AccumulatedDistance_um','Accumulated distance (\mum)',outFolder,'Panel_C_AccumulatedDistance');
plotMetric(AllMetrics,'EuclideanDistance_um','Euclidean displacement (\mum)',outFolder,'Panel_D_EuclideanDisplacement');
plotMetric(AllMetrics,'MigrationIndex','Migration index',outFolder,'Panel_E_MigrationIndex');
plotMetric(AllMetrics,'MotilityCoefficient_um2_hr','Motility coefficient (\mum^2/h)',outFolder,'Panel_F_MotilityCoefficient');
plotMetric(AllMetrics,'PersistenceTime_hr','Persistence time (h)',outFolder,'Panel_G_PersistenceTime');

%% DISPLACEMENT OVER TIME
figure('Color','w'); hold on;

for g = 1:length(groupOrder)

    groupName = groupOrder(g);

    if groupName == "U87 WT"
        colorNow = WT_color;
    else
        colorNow = KO_color;
    end

    Tg = AllCoords(AllCoords.Group == groupName,:);
    timePoints = unique(Tg.Time_hr);

    meanDisp = zeros(length(timePoints),1);
    semDisp = zeros(length(timePoints),1);

    for t = 1:length(timePoints)
        vals = Tg.NetDisplacement_um(Tg.Time_hr == timePoints(t));
        vals = vals(~isnan(vals));
        meanDisp(t) = mean(vals);
        semDisp(t) = std(vals)/sqrt(length(vals));
    end

    upper = meanDisp + semDisp;
    lower = meanDisp - semDisp;

    fill([timePoints; flipud(timePoints)], ...
         [upper; flipud(lower)], ...
         colorNow, ...
         'FaceAlpha',0.15, ...
         'EdgeColor','none');

    plot(timePoints,meanDisp,'Color',colorNow,'LineWidth',2);
end

box off;
xlabel('Time (h)');
ylabel('Net displacement (\mum)');
legend(cellstr(groupOrder),'Box','off');
set(gca,'FontName','Arial','FontSize',18,'LineWidth',1,'TickDir','out');

saveas(gcf,fullfile(outFolder,'Panel_H_DisplacementOverTime.png'));

%% SUMMARY TABLE
summaryTable = table();

for g = 1:length(groupOrder)

    groupName = groupOrder(g);
    Tg = AllMetrics(AllMetrics.Group == groupName,:);

    row = table();
    row.Group = groupName;
    row.N = height(Tg);
    row.Speed_mean = mean(Tg.Speed_um_hr,'omitnan');
    row.Speed_sem = std(Tg.Speed_um_hr,'omitnan')/sqrt(height(Tg));
    row.Accumulated_mean = mean(Tg.AccumulatedDistance_um,'omitnan');
    row.Accumulated_sem = std(Tg.AccumulatedDistance_um,'omitnan')/sqrt(height(Tg));
    row.Euclidean_mean = mean(Tg.EuclideanDistance_um,'omitnan');
    row.Euclidean_sem = std(Tg.EuclideanDistance_um,'omitnan')/sqrt(height(Tg));
    row.MigrationIndex_mean = mean(Tg.MigrationIndex,'omitnan');
    row.MigrationIndex_sem = std(Tg.MigrationIndex,'omitnan')/sqrt(height(Tg));
    row.MotilityCoefficient_mean = mean(Tg.MotilityCoefficient_um2_hr,'omitnan');
    row.MotilityCoefficient_sem = std(Tg.MotilityCoefficient_um2_hr,'omitnan')/sqrt(height(Tg));
    row.PersistenceTime_mean = mean(Tg.PersistenceTime_hr,'omitnan');
    row.PersistenceTime_sem = std(Tg.PersistenceTime_hr,'omitnan')/sqrt(height(Tg));

    summaryTable = [summaryTable; row];
end

%% SAVE OUTPUTS
writetable(summaryTable,fullfile(outFolder,'Summary_Table.xlsx'));
writetable(StatsTable,fullfile(outFolder,'WT_vs_COSMC_KO_Statistical_Analysis.xlsx'));
writetable(AllMetrics,fullfile(outFolder,'All_Track_Metrics.xlsx'));
writetable(AllCoords,fullfile(outFolder,'All_Track_Coordinates.xlsx'));
writetable(RemovedTracks,fullfile(outFolder,'Removed_Tracks_QC.xlsx'));

disp('Analysis complete.');
disp(['Final tracks analyzed: ', num2str(height(AllMetrics))]);
disp(['Bad tracks removed: ', num2str(height(RemovedTracks))]);
disp(['Results saved in: ', outFolder]);

disp(' ');
disp('Summary table:');
disp(summaryTable);

disp(' ');
disp('WT vs COSMC-KO statistical analysis:');
disp(StatsTable);

%% FUNCTION: MOTILITY COEFFICIENT
function motilityCoefficient = calculateMotilityCoefficient(x0,y0,dt_hr)

    n = length(x0);
    maxLag = min(n-1,round(5/dt_hr));

    if maxLag < 3
        motilityCoefficient = NaN;
        return
    end

    lags = (1:maxLag)';
    msd = zeros(length(lags),1);

    for k = 1:length(lags)
        lag = lags(k);
        dx = x0(1+lag:end) - x0(1:end-lag);
        dy = y0(1+lag:end) - y0(1:end-lag);
        msd(k) = mean(dx.^2 + dy.^2);
    end

    lagTime = lags * dt_hr;
    useIdx = lagTime <= min(2,max(lagTime));

    if sum(useIdx) < 3
        useIdx = false(size(lagTime));
        useIdx(1:min(5,length(lagTime))) = true;
    end

    p = polyfit(lagTime(useIdx),msd(useIdx),1);
    motilityCoefficient = p(1)/4;
end

%% FUNCTION: PERSISTENCE TIME
function persistenceTime_hr = calculatePersistenceTime(x_um,y_um,dt_hr)

    dx = diff(x_um);
    dy = diff(y_um);

    step = sqrt(dx.^2 + dy.^2);
    valid = step > 0;

    if sum(valid) < 5
        persistenceTime_hr = NaN;
        return
    end

    ux = dx(valid)./step(valid);
    uy = dy(valid)./step(valid);

    n = length(ux);
    maxLag = min(n-1,round(5/dt_hr));

    if maxLag < 2
        persistenceTime_hr = NaN;
        return
    end

    lags = (1:maxLag)';
    autocorrValue = zeros(length(lags),1);

    for k = 1:length(lags)
        lag = lags(k);
        autocorrValue(k) = mean(ux(1:end-lag).*ux(1+lag:end) + ...
                                uy(1:end-lag).*uy(1+lag:end));
    end

    lagTime = lags * dt_hr;
    threshold = 1/exp(1);

    idx = find(autocorrValue <= threshold,1,'first');

    if isempty(idx)
        persistenceTime_hr = lagTime(end);
    else
        persistenceTime_hr = lagTime(idx);
    end
end

%% FUNCTION: MANUAL BOXPLOT + POINTS + STATS
function plotMetric(T,metricName,yLabelText,outFolder,figName)

    figure('Color','w','Units','inches','Position',[1 1 3.2 4]);
    hold on;

    WT_vals = T.(metricName)(T.Group=="U87 WT");
    KO_vals = T.(metricName)(T.Group=="COSMC-KO");

    WT_vals = WT_vals(~isnan(WT_vals));
    KO_vals = KO_vals(~isnan(KO_vals));

   WT_color = [0.00 0.3 0.70];
   KO_color = [0.5 0.00 0.15];

    drawManualBox(1,WT_vals,WT_color);
    drawManualBox(2,KO_vals,KO_color);

    rng(1);

    jitterWT = 0.16*(rand(length(WT_vals),1)-0.5);
    jitterKO = 0.16*(rand(length(KO_vals),1)-0.5);

    scatter(1+jitterWT,WT_vals,38,WT_color,'filled', ...
        'MarkerEdgeColor','k','LineWidth',0.4,'MarkerFaceAlpha',0.80);

    scatter(2+jitterKO,KO_vals,38,KO_color,'filled', ...
        'MarkerEdgeColor','k','LineWidth',0.4,'MarkerFaceAlpha',0.80);

    WT_mean = mean(WT_vals);
    KO_mean = mean(KO_vals);

    WT_sem = std(WT_vals)/sqrt(length(WT_vals));
    KO_sem = std(KO_vals)/sqrt(length(KO_vals));

    errorbar(1,WT_mean,WT_sem,'k','LineWidth',1.4,'CapSize',10);
    errorbar(2,KO_mean,KO_sem,'k','LineWidth',1.4,'CapSize',10);

    plot(1,WT_mean,'kd','MarkerFaceColor','k','MarkerSize',6);
    plot(2,KO_mean,'kd','MarkerFaceColor','k','MarkerSize',6);

    allVals = [WT_vals; KO_vals];

    [p,~,~,~] = welchTTest_manual(WT_vals,KO_vals);

    if p < 0.0001
        sigText = '****';
        pText = 'P < 0.0001';
    elseif p < 0.001
        sigText = '***';
        pText = sprintf('P = %.1e',p);
    elseif p < 0.01
        sigText = '**';
        pText = sprintf('P = %.3f',p);
    elseif p < 0.05
        sigText = '*';
        pText = sprintf('P = %.3f',p);
    else
        sigText = 'ns';
        pText = sprintf('P = %.3f',p);
    end

    yMax = max(allVals);
    yMin = min(allVals);
    yRange = yMax - yMin;

    if yRange == 0
        yRange = 1;
    end

    yLine = yMax + 0.12*yRange;
    yStar = yMax + 0.18*yRange;
    yPval = yMax + 0.35*yRange;

    plot([1 1 2 2], ...
         [yLine-0.03*yRange yLine yLine yLine-0.03*yRange], ...
         'k','LineWidth',1.2);

    text(1.5,yStar,sigText, ...
        'HorizontalAlignment','center', ...
        'FontName','Arial', ...
        'FontSize',18, ...
        'FontWeight','bold');

    text(1.5,yPval,pText, ...
        'HorizontalAlignment','center', ...
        'FontName','Arial', ...
        'FontSize',18);

    ylim([yMin-0.08*yRange yMax+0.38*yRange]);

    xlim([0.5 2.5]);
    xticks([1 2]);
    xticklabels({'WT','COSMC-KO'});

    ylabel(yLabelText);

    box off;

    set(gca,'FontName','Arial','FontSize',16,'LineWidth',1.2,'TickDir','out');

    saveas(gcf,fullfile(outFolder,[figName '_boxplot_with_statistics.png']));
end

%% FUNCTION: DRAW MANUAL BOXPLOT
function drawManualBox(xpos,vals,colorNow)

    vals = vals(~isnan(vals));

    q1 = prctile_manual(vals,25);
    q2 = median(vals);
    q3 = prctile_manual(vals,75);

    iqrVal = q3 - q1;

    lowerWhisker = min(vals(vals >= q1 - 1.5*iqrVal));
    upperWhisker = max(vals(vals <= q3 + 1.5*iqrVal));

    boxWidth = 0.38;

    rectangle('Position',[xpos-boxWidth/2,q1,boxWidth,q3-q1], ...
        'FaceColor',[colorNow 0.18], ...
        'EdgeColor','k', ...
        'LineWidth',1.2);

    plot([xpos-boxWidth/2 xpos+boxWidth/2],[q2 q2], ...
        'k','LineWidth',1.5);

    plot([xpos xpos],[lowerWhisker q1], ...
        'k','LineWidth',1.2);

    plot([xpos xpos],[q3 upperWhisker], ...
        'k','LineWidth',1.2);

    plot([xpos-0.10 xpos+0.10],[lowerWhisker lowerWhisker], ...
        'k','LineWidth',1.2);

    plot([xpos-0.10 xpos+0.10],[upperWhisker upperWhisker], ...
        'k','LineWidth',1.2);
end

%% FUNCTION: MANUAL WELCH T-TEST
function [p,tstat,df,ci] = welchTTest_manual(x,y)

    x = x(~isnan(x));
    y = y(~isnan(y));

    nx = length(x);
    ny = length(y);

    meanX = mean(x);
    meanY = mean(y);

    varX = var(x);
    varY = var(y);

    se = sqrt(varX/nx + varY/ny);

    tstat = (meanX - meanY) / se;

    df = (varX/nx + varY/ny)^2 / ...
        ((varX/nx)^2/(nx-1) + (varY/ny)^2/(ny-1));

    p = tTwoTailP_manual(tstat,df);

    tcrit = tCritical_manual(0.975,df);

    diffMean = meanX - meanY;

    ci = [diffMean - tcrit*se;
          diffMean + tcrit*se];
end

%% FUNCTION: TWO-TAILED T P-VALUE
function p = tTwoTailP_manual(t,df)

    t = abs(t);
    x = df / (df + t^2);

    p = betainc(x,df/2,0.5);

    if p < 0
        p = 0;
    end

    if p > 1
        p = 1;
    end
end

%% FUNCTION: MANUAL T CRITICAL VALUE
function tcrit = tCritical_manual(prob,df)

    targetP = 2*(1-prob);

    f = @(t) tTwoTailP_manual(t,df) - targetP;

    try
        tcrit = fzero(f,[0 100]);
    catch
        tcrit = 1.96;
    end
end

%% FUNCTION: SIGNIFICANCE STARS
function sig = significanceStars(p)

    if p < 0.0001
        sig = '****';
    elseif p < 0.001
        sig = '***';
    elseif p < 0.01
        sig = '**';
    elseif p < 0.05
        sig = '*';
    else
        sig = 'ns';
    end
end

%% FUNCTION: MANUAL PERCENTILE
function q = prctile_manual(x,p)

    x = sort(x(:));
    n = length(x);

    if n == 1
        q = x;
        return
    end

    pos = 1 + (p/100)*(n-1);
    lower = floor(pos);
    upper = ceil(pos);

    if lower == upper
        q = x(lower);
    else
        q = x(lower) + (pos-lower)*(x(upper)-x(lower));
    end
end
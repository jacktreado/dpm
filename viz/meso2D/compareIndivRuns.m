%% Script to compare multiple individual files of mesoEnthalpyMin sims

clear;
close all;
clc;

% list of files to compare
simloc = 'local/mesoHMin2D_data';
simlist = {'mesoHMin2D_N32_n32_ca1.14_kb0.2_be200_h0.5_da0.2_dl1.5_cL0.5_cB1_t0m0.5_P1e-6_seed12.posctc',...
    'mesoHMin2D_N32_n32_ca1.14_kb0.2_be200_h0.4_da0.2_dl1.5_cL0.5_cB1_t0m0.5_P1e-6_seed12.posctc',...
    'mesoHMin2D_N32_n32_ca1.14_kb0.2_be200_h0.3_da0.2_dl1.5_cL0.5_cB1_t0m0.5_P1e-6_seed12.posctc',...
    'mesoHMin2D_N32_n32_ca1.14_kb0.2_be200_h0.2_da0.2_dl1.5_cL0.5_cB1_t0m0.5_P1e-6_seed12.posctc'};
NSIMS = length(simlist);

% get cell list of data
mesoDataCell = cell(NSIMS,1);
for ss = 1:NSIMS
    fprintf('\n*** Sim %d/%d\n',ss,NSIMS);
    mesoDataCell{ss} = readMesoNetworkCTCS2D([simloc '/' simlist{ss}]);
end

% limit phi
phiMin = 0.48;

% plotting info
lineCLR = jet(NSIMS);
mrkrSYM = 'so<>^dvph';
NSYM = length(mrkrSYM);

% error bar face alpha
faceAlpha = 0.2;

% helpful supercells
phipatch = cell(NSIMS,1);
for ss = 1:NSIMS
    % struct
    mesoDataTmp = mesoDataCell{ss};
    
    % get data
    NFRAMES = mesoDataTmp.NFRAMES;
    NCELLS = mesoDataTmp.NCELLS;
    x = mesoDataTmp.x;
    y = mesoDataTmp.y;
    r = mesoDataTmp.r;
    L = mesoDataTmp.L;

    % compute patch phi
    phipatchtmp = zeros(NFRAMES,1);
    for ff = 1:NFRAMES
        xf = x(ff,:);
        yf = y(ff,:);
        rf = r(ff,:);
        apatch = zeros(NCELLS,1);
        for nn = 1:NCELLS
            xtmp = xf{nn};
            ytmp = yf{nn};
            rtmp = rf{nn};
            cx = mean(xtmp);
            cy = mean(ytmp);
            rx = xtmp - cx;
            ry = ytmp - cy;
            rads = sqrt(rx.^2 + ry.^2);
            xtmp = xtmp + sin(pi/3)*rtmp.*(rx./rads);
            ytmp = ytmp + sin(pi/3)*rtmp.*(ry./rads);
            apatch(nn) = polyarea(xtmp,ytmp);
        end
        phipatchtmp(ff) = sum(apatch)/(L(ff,1)*L(ff,2));
    end
    phipatch{ss} = phipatchtmp;
end

%% Make plots, compare

% calA vs poro
figure(1), clf, hold on, box on;
for ss = 1:NSIMS
    % get sim data
    phipatchtmp = phipatch{ss};
    mesoDataTmp = mesoDataCell{ss};
    
    % set idx
    idx = phipatchtmp > phiMin;
    NFRAMES = sum(idx);
    
    % extract data
    p = mesoDataTmp.p(idx,:);
    a = mesoDataTmp.a(idx,:);
    calA = p.^2./(4.0*pi*a);
    calAMean = mean(calA,2);
    calAStd = std(calA,0,2);
    
    % define scaled porosity
    simporo = max(phipatchtmp) - phipatchtmp;
    
    % plot
    clrtmp = lineCLR(ss,:);
    symtmp = mrkrSYM(mod((ss-1),NSYM)+1);
    linestr = ['-k' symtmp];
    patchErrorBar(simporo,calAMean,calAStd,linestr,10,clrtmp,clrtmp,faceAlpha);
end

% Get Ambrose Data
ambroseData = load('/Users/jacktreado/Jamming/Flowers/structure/plant/ambroseMesoCells/ambroseData.mat');
totalData = load('/Users/jacktreado/Jamming/Flowers/structure/plant/ambroseMesoCells/totalData.mat');
porosity = totalData.porosity;
calAMean = totalData.calAMean;
calAMin = totalData.calAMin;
calAMax = totalData.calAMax;
calAMeas = totalData.calAMeas;
dCalAMin = totalData.dCalAMin;
dCalAMax = totalData.dCalAMax;

ct01CalAMean = totalData.ct01CalAMean;
ct01CalAStd = totalData.ct01CalAStd;
ct01Porosity = totalData.ct01Porosity;

ct17CalAMean = totalData.ct17CalAMean;
ct17CalAStd = totalData.ct17CalAStd;
ct17Porosity = totalData.ct17Porosity;

ct65CalAMean = totalData.ct65CalAMean;
ct65CalAStd = totalData.ct65CalAStd;
ct65Porosity = totalData.ct65Porosity;

% Plot Ambrose Data
errorbar(porosity,calAMeas,dCalAMin,dCalAMax,'k>','markersize',14,'markerfacecolor','r');
errorbar(ct01Porosity,ct01CalAMean,ct01CalAStd,'ko','markersize',12,'markerfacecolor','b');
errorbar(ct17Porosity,ct17CalAMean,ct17CalAStd,'ks','markersize',12,'markerfacecolor','g');
errorbar(ct65Porosity,ct65CalAMean,ct65CalAStd,'kd','markersize',12,'markerfacecolor','r');

% label plot
ylabel('$\mathcal{A}$','Interpreter','latex');
xlabel('$\phi_{\rm max}-\phi$','Interpreter','latex');
ax = gca;
ax.FontSize = 24;
ax.XLim = [0 0.5];






% Contact number over time
figure(2), clf, hold on, box on;
for ss = 1:NSIMS
    % get sim data
    phipatchtmp = phipatch{ss};
    mesoDataTmp = mesoDataCell{ss};
    
    % set idx
    idx = phipatchtmp > phiMin;
    NFRAMES = sum(idx);
    
    % extract data
    zc = mesoDataTmp.zc(idx,:);
    zcMean = mean(zc,2);
    zcStd = std(zc,0,2);
    
    % define scaled porosity
    simporo = max(phipatchtmp) - phipatchtmp;
    
    % plot
    clrtmp = lineCLR(ss,:);
    symtmp = mrkrSYM(mod((ss-1),NSYM)+1);
    linestr = ['-k' symtmp];
    patchErrorBar(simporo,zcMean,zcStd,linestr,10,clrtmp,clrtmp,faceAlpha);
end
ylabel('$z$','Interpreter','latex');
xlabel('$\phi_{\rm max}-\phi$','Interpreter','latex');
ax = gca;
ax.FontSize = 24;














%% Setup
% set current directory to folder where data and surfstat functions are stored
cd '/Volumes/Lacie2/Alisa_Zoltowski/CorticalMorphology/Code/cascio' 
addpath(genpath('/Volumes/Lacie2/Alisa_Zoltowski/CorticalMorphology/Code/surfstat')); % path to surfstat lib
addpath(genpath('/Volumes/Lacie2/Alisa_Zoltowski/CorticalMorphology/Code/SurfStatView')); % path to new figure code
cmeasure = 'lgi';   % EDIT PER MODEL, measure of interest: lgi, sd or ct

%% Read data
load(sprintf('data/%s.mat',cmeasure));  % surface measure
load('env/environment.mat');  % template surfaces, parcellation boundary, etc.

% load full demographic dataset with all intended variables
load('data/demographics_updated_ados.mat');  % final dataset for arz analyses
demographic = demographic_ados2;



%% Data filtering
% filter demographic/structural data based on included subjects
% pick one of the below, based on contrast

% for diagnostic group models, find subjects with valid Dx code (should be all)
subj = find(demographic.Dx == 0 | demographic.Dx == 1); 

% for ADOS models, find subjects with non-missing ADOS score
%subj = find(~isnan(demographic.ADOS)); 

% for age by diagnostic group models, remove influential points
%subj = find(demographic.Age < 50); % most age by dx analyses (LGI, CT)
%subj= [1:70,72:161,163:255,257:573]; % for sulcal depth TD>ASD by age, outliers: 71, 162, 256 

% filter corresponding structural data
Y=Y0(subj,:);

%% Model term initialization

% convert 0,1 categorical variables to corresponding labels (dx and sex)
dx = cell(length(demographic.Dx),1);
dx(demographic.Dx == 0) = {'TD'};
dx(demographic.Dx == 1) = {'ASD'};
sex = cell(length(demographic.Sex),1);
sex(demographic.Sex == 0) = {'F'};
sex(demographic.Sex == 1) = {'M'};

% create variable for each model covariates
tbv = demographic.TBV(subj);
age = demographic.Age(subj);
%ados = demographic.ADOS(subj);  % only run this line for ADOS models
sex = sex(subj);
dx = dx(subj);
scanner = demographic.Scan(subj);
subject = demographic.Subject(subj);

% convert covariates to formatted terms for surfstat
TBV = term( tbv );
Age = term( age );
Sex = term( sex );
%ADOS = term( ados );    % only run this line for ADOS models
Dx = term( dx );
Scanner = term( scanner );
Subject = term( subject );

%% Model setup/fitting

% define me: contrast and model (run only one contrast, M)

% Diagnostic group contrast and model
contrast = Dx.TD - Dx.ASD;  % EDIT: Specify direction (TD-ASD or ASD-TD)
M=1+Age+Sex+Dx+random(Scanner)+random(Subject)+I;   % Dx model


% Age by diagnostic group contrast
%contrast = age .* Dx.ASD - age .* Dx.TD;  % EDIT: Specify direction
%(age .* Dx.ASD - age .* Dx.TD or age .* Dx.TD - age .* Dx.ASD)
%M=1+Age+Sex+Dx+Age*Dx+random(Scanner)+random(Subject)+I;   % Age by dx model

% ados (continuous, no contrast needed)
%M=1+Age+Sex+ADOS+random(Scanner)+random(Subject)+I;   % ADOS model

Y(:,sum(abs(Y))==0) = rand(size(Y(:,sum(abs(Y))==0)))*eps;   % prevent numerical instability
slm = SurfStatLinMod(Y,M,surfwhite);  % run linear regression models, defined above
slm = SurfStatT( slm, contrast ); % t stat for models with specified contrast (dx, age by dx)
%slm = SurfStatT( slm, -ados);      % t stat for ados models, -ados for reverse
% identify clusters and p-values:
[ pval, peak, clus ] = SurfStatP( slm, mask, 0.01); % cluster raw-p

%% Save and reload models, once already ran (internal purposes)
% Examples for labeling:
% dx example: lgi_td-asd 
% age example: lgi_td-asd_byage_no54 or _nooutlier
% ados example: lgi_ados_pos or _neg

% EDIT: which measure and contrast, following examples above
contrastLabel = 'ct_td-asd'; 

% save result files per contrast
%save(sprintf('FinalModels/slm/%s.mat',contrastLabel),'slm'); %model output
%save(sprintf('FinalModels/clus/%s_clus.mat',contrastLabel),'clus'); %clusters
%save(sprintf('FinalModels/clus/%s_peak.mat',contrastLabel),'peak'); %peaks
%save(sprintf('FinalModels/clus/%s_pval.mat',contrastLabel),'pval'); %p

% load whichever model you need, following conventions above:
%load(sprintf('FinalModels/slm/%s.mat',contrastLabel),'slm'); %model
%load(sprintf('FinalModels/clus/%s_clus.mat',contrastLabel),'clus'); %clusters
%load(sprintf('FinalModels/clus/%s_peak.mat',contrastLabel),'peak'); %peaks
%load(sprintf('FinalModels/clus/%s_pval.mat',contrastLabel),'pval'); %p

%% Primary visualization (Figures 2-4): signficant regions after RFT
pval.mask(:)=true; % apply p-value mask to plot
% plot using edited SurfStatViewer (edited by IL)
figure; SurfStatView2( pval, surfinfl, 'Cluster p-value, corrected' );

%% Calculate residual (adjusted Y) if needed for inset and/or supplemental plots

% only run one of the below lines
X=Dx; % contrast: group diff
%X=Age*Dx; % constrast: group by covariate
%X=ADOS; % contrast: ADOS
% note, can analyze this for non-contrast covariates as well for plotting
% (e.g. Sex or Age)

G=term(1);
XG=X+G+X*G;
Yhat=slm.X*slm.coef;
XGmat=double(XG);
Yhatadj=XGmat*pinv(XGmat)*Yhat;
c=mean(Y)-mean(Yhatadj);
Yhatadj=Yhatadj+c;
Yadj=Yhatadj+Y-Yhat;
Yadj(:,~mask) = 0;

%% Figure 3 inset plot: diagnostic group by age interaction

% Interaction plot for specified cluster
clusid=1; % EDIT: define cluster number
cpval = clus.P(clusid); % select vertices within specified cluster
figure; SurfStatPlot( age, double(mean(Y(:, pval.C == cpval),2)), 1, dx ); ylabel(sprintf('%s',cmeasure)); % plot 

% Plots not published- influence plots to determine outlying data points by age
% Initial analyses determined significant clusters with whole group included, per model and contrast 
% Influential data points were then assessed within each cluster
clus_avg = double(mean(Y(:, pval.C == cpval),2)); % average structural index within cluster, number specified above
tbl = table(age,dx,clus_avg,'VariableNames',{'age','dx','clus_avg'}); % select variables for age by diagnostic group regression model
mdl = fitlm(tbl,'clus_avg~age*dx'); % model age by diagnostic group interaction
figure; plotDiagnostics(mdl,'cookd') % calculate Cook's d to determine outliers

%% Figure 4, inset plot: ADOS visualization
clusid=1;   % EDIT: define cluster number
cpval = clus.P(clusid); % select data points per cluster
color = parula(5); % define color scheme
figure; SurfStatPlot( ados,mean(Yadj( :, pval.C == cpval ),2), Age+Sex+Dx, [], 'Color',color(1,:) ); xlabel('ADOS'); ylabel(sprintf('%s',cmeasure)); % plot

%% Tables 2-4 information: coordinate location (of peaks), cluster size, maximum coefficient (effect)
% add path to where files with vertex-to-region-id files are stored
addpath(genpath('/Volumes/Lacie2/Alisa_Zoltowski/CorticalMorphology/Code/')); 
% read in vertex numbers to id's tables
lh = readtable('lh.parc.txt','ReadVariableNames',false);
rh = readtable('rh.parc.txt','ReadVariableNames',false);
bilat = vertcat(lh,rh);

% for each cluster, get related vertex numbers and translate to id's
n_clus = sum(clus.P < 0.05); % find out how many clusters below <0.05 threshold
clus_info = table(0,"",0,'VariableNames',{'NumVertices','IDs','MaxEffect'});
for ii = 1:n_clus
    clusid=ii; %cluster id
    % from peak vertex list, identify corresponding vertex numbers, id
    % codes
    vert_num = peak.vertid(peak.clusid==clusid);
    vert_id = bilat.Var1(vert_num);
    verts=unique(vert_id); % get unique region codes from these vertices
    clus_info.IDs(clusid)=string(sprintf('%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d',verts)); %fill in up to 21 region ID codes, edit if any have more
    % number of vertices per cluster, from clus data structure
    clus_info.NumVertices(clusid)=clus.nverts(clusid);
    % max effect (coefficient for single term or difference for contrast)
    ef_verts = slm.ef(vert_num); 
    clus_info.MaxEffect(clusid)=max(ef_verts); 
    
end

% final note- this loop created framework for data tables, but regions were labeled
% with additional csv file offline

%% Supplemental Figure 1 (main): Correlations between CT, SD with LGI in specific clusters
% load other two indices (LGI should be Y0)
Y0 = load('data/lgi.mat'); 
ct = load('data/ct.mat'); 
sd = load('data/sd.mat'); 

% identify cluster of interest and get vertex column #s
clusid=3; % EDIT per cluster
cpval = clus.P(clusid);

% average each index per cluster
ct_clusI = ct.Y0(:,pval.C == cpval);
ct_clusAvg = mean(ct_clusI,2);

sd_clusI = sd.Y0(:,pval.C == cpval);
sd_clusAvg = mean(sd_clusI,2);

lgi_clusI = Y0(:,pval.C == cpval);
lgi_clusAvg = mean(lgi_clusI,2);

% plot LGI with SD and CT
color = parula(5); % define color scheme
figure; SurfStatPlot2( sd_clusAvg, lgi_clusAvg, 1, 1, dx); xlabel('Sulcal Depth (unadjusted)'); ylabel('LGI (unadjusted)');
figure; SurfStatPlot2( ct_clusAvg, lgi_clusAvg, 1, 1, dx); xlabel('Cortical Thickness (unadjusted)'); ylabel('LGI (unadjusted)');

%% Supplemental Figure 1 (left panel): box plot for significant regions (in group comparison models)
clusid=1;   % EDIT per cluster
cpval = clus.P(clusid);

% define color scheme
color = parula(5);
cc = color([3,1],:);

% initialize box plot and labels
figure('Position', [10 10 300 700]);
boxplot( mean(Yadj( :, pval.C == cpval ),2), dx, 'Notch','off', 'Colors', cc, 'Labels',{'NT','AUT'});
ylabel(sprintf('%s',cmeasure));
xlabel(sprintf('cluster p: %f',cpval));

% add box plot lines
h = findobj(gca,'tag','Median');
for j=1:length(h)
    set(h(j),'linestyle','-.','Color',[0 0 0],{'linew'},{2});
end
set(findobj(gca,'tag','Outliers'),'MarkerEdgeColor',[0 0 0]);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    c = cc(3-j,:);
    patch(get(h(j),'XData'),get(h(j),'YData'), c,'FaceAlpha',.81,'EdgeColor',cc(3-j,:)); % .51 alpha
end


%% Supplemental Figure 2, visualizing measure trends by age
figure; SurfStatPlot2( age, mean(Y,2),Sex+Dx+Scanner,1,dx); xlabel('Age'); ylabel(sprintf('%s',cmeasure));


%% Supplemental Figure 3, random effects by scanner protocol
% main supplemental plot for random effects
figure; SurfStatView( slm.r(1,:).*mask_b, surfinfl, 'Correlation within protocol' ); SurfStatColLim( [0 1] );

% not published: option to plot within-subject random effects
figure; SurfStatView( slm.r(2,:).*mask_b, surfinfl, sprintf('Correlation within subject (r=%.2f)',mean(slm.r(2,:))) ); SurfStatColLim( [0 1] );

%% Supplemental Figure 3, plot of biological sex differences

% Main supplemental plot by biological sex
% calculate group difference across all brain regions and then plot via surfstat
diff = mean(Yadj(contains(sex,'M'),:))-mean(Yadj(contains(sex,'F'),:));
figure; SurfStatView( diff, surfinfl, 'adjusted difference (M - F)' );
cmap = [0 0 0; jet];
SurfStatColormap( jet );
SurfStatColLim( [-2 2] );

% Not published- you can create a similar plot for Dx group differences
% as continuous estimate, not just significant regions
diff = mean(Yadj(contains(dx,'TD'),:))-mean(Yadj(contains(dx,'ASD'),:));
figure; SurfStatView( diff, surfinfl, 'adjusted difference (TD - ASD)' );
cmap = [0 0 0; jet];
SurfStatColormap( jet );
SurfStatColLim( [-2 2] );



%% Overall measure mean per region- optional visualization, not published:
meansubj = mean( double( Yadj ) );
figure; SurfStatView( meansubj, surfinfl, 'adjusted mean' ); SurfStatColLim( [1 5] );   % 1-5: ct, 0-30: sd, 1-15: lgi

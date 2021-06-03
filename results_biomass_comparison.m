files = dir('Results\Results_irreversibility_considered_biomass\*.mat');
bio=[];
mrc_any=[];mrc_ma=[];
species=[];

for f=1:length(files)
    R=load(strcat(files(f).folder,'/',files(f).name));

    species=[species;R.species_list(R.s)];
    
    if isempty(R.Solma.f)
        R.Solma.f=nan;
    end

    bio=[bio;R.Sol2.f R.Solany.f R.Solma.f];
    mrc_any=[mrc_any;(size(R.model.S,1)-size(R.model_any.S,1))/size(R.model.S,1) ...
            (size(R.model.S,2)-size(R.model_any.S,2))/size(R.model.S,2) ...
            (size(R.model.A,1)-size(R.model_any.A,1))/size(R.model.A,1)];

    mrc_ma=[mrc_ma;(size(R.model.S,1)-size(R.model_ma.S,1))/size(R.model.S,1) ...
            (size(R.model.S,2)-size(R.model_ma.S,2))/size(R.model.S,2) ...
            (size(R.model.A,1)-size(R.model_ma.A,1))/size(R.model.A,1)];
end
figure
subplot(1,2,1)
labels={'\itA. niger iMA871'
'\itA. thaliana AraCore'
'\itC. reinhardtii iCre1355'
'\itE. coli iJO1366'
'\itM. acetivorans iMB745'
'\itM. barkeri iAF692'
'\itM. musculus'
'\itM. tuberculosis iNJ661m'
'\itN. pharaonis'
'\itP. putida iJN746'
'\itT. maritima iLJ478'
'\itS. cerevisiae Yeast8'};

plot(abs(bio(:,1)),abs(bio(:,2)),'k.','MarkerSize',10)
hold on
plot(abs(bio([3 8],1)),abs(bio([3 8],2)),'r.','MarkerSize',10)
set(gca,'XScale','Log','YScale','Log',...
    'XTick',[1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3],'YTick',[1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3])
ylabel('optimal biomass reduced network [h^{-1}]')
xlabel('optimal biomass original network [h^{-1}]')
ylim([1e-3 1e3])
xlim([1e-3 1e3])

axes('Position',[.2 .2 .2 .2])
box on
plot(abs(bio([3 8],1)),abs(bio([3 8],2)),'r.','MarkerSize',10)
ylim([0.05 0.055])
xlim([0.05 0.055])

% subplot(2,2,2)
% plot(abs(bio(:,1)),abs(bio(:,3)),'k.','MarkerSize',10)
% set(gca,'XScale','Log','YScale','Log')
% ylabel('optimal biomass reduced network [h^{-1}]')
% xlabel('optimal biomass original network [h^{-1}]')
subplot(1,2,2)
bar(1:12,mrc_any(:,1)*100)
set(gca,'XTick',1:12,'XTickLabels',labels,'XTickLabelRotation',45)
ylabel('species reduction [%]')
% subplot(2,2,4)
% bar(1:12,mrc_ma(:,1)*100)
% set(gca,'XTick',1:12,'XTickLabels',labels,'XTickLabelRotation',45)

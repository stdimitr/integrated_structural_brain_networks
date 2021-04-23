function [isbn CIJtree coef]=integrated_structural_brain_network(diffusion_sbn)

%% An Graph-based Integration of different Diffusion Metrics 
%% Into a Single Structural Brain Network
%%% INPUT : a 3D matrix with dimensions equal to = number of diffusion
%%%         metrics x rois x rois
%%%OUTPUT : isbn = The integrated weighted structural brain network with dimension
%%%                 equals to 90 x 90
%%%        CIJtree = The topological filtered with OMST method integrated structural brain network with dimension
%%%                 equals to 90 x 90
%%%            coef = coefficients of contribution of every diffusion-metric based structural brain network to the integrated weighted structural brain network

%%% For further details see the paper:
%% Paper: Dimitriadis et al., 2017.Improving the Reliability of Network Metrics in Structural Brain Networks by Integrating Different Network Weighting Strategies into a Single Graph. Front. Neurosci., 
%%19 December 2017 | https://doi.org/10.3389/fnins.2017.00694

[no_dmetrics rois rois]=size(diffusion_sbn);

thresholded=zeros(no_dmetrics ,rois,rois);

%%% topologically filtered of brain networks with OMST

for k=1:no_dmetrics 
    [nCIJtree CIJtree mdeg  globalcosteffmax costmax E]=threshold_omst_gce_wu(squeeze(diffusion_sbn(k,:,:)),0);
    thresholded(k,:,:)=CIJtree;
end

%%distance matrix with graph diffusion distance metric

dist=zeros(no_dmetrics ,no_dmetrics );

for k=1:no_dmetrics 
    for l=(k+1):no_dmetrics 
        A1=squeeze(thresholded(k,:,:));
        A2=squeeze(thresholded(l,:,:));
        [gdd,t,t_upperbound]=compute_gdd(A1,A2);
        dist(k,l)=gdd;
        dist(l,k)=dist(k,l);
    end
end

%% take the sum of rows
sum1=sum(dist);

%normalize and keep the coefficients of the linear combination of
%individualized brain networks
coef=sum1./sum(sum1);

%%integrated brain network
isbn=zeros(90,90);

for k=1:no_dmetrics 
    isbn=isbn + coef(k).*squeeze(thresholded(k,:,:));
end

%% normalize the integrated brain network
isbn_norm=isbn/max(isbn);

%% topological filtering of integrated brain network with OMST
 [nCIJtree CIJtree mdeg  globalcosteffmax costmax E]=threshold_omst_gce_wu(isbn_norm,1);
 
 %% PLOT THE COEFFICIENTS THAT DEMONSTRATE THE CONTRIBUTION OF EACH DIFFUSION-BASED
 %% STRUCTURAL BRAIN NETWORK TO THE INTEGRATED SBN
 figure(2),imagesc(coef) ; colorbar
           title('Coefficients of the linear combination of brain networks')
 %%%PLOT THE INTEGRATED SBN          
 figure(3), imagesc(CIJtree)           
           
           

clc
 close all
 clear all
 
 %%
 % data generation
 data = randgen([4 5],[9 0;0 9],[10 18],[6 0;0 6],[15 8],[4 0;0 4]);
 N = size(data,1);
 
 %%
 % real-valued similarities between data points
 % calculate s(i,k)
 s = zeros(N,N);
 for i=1:N
     for k=1:N
         s(i,k) = -norm((data(i,:)-data(k,:)),2);
     end
 end
 
 %%
 % Preferences
 % init_pref = median(s(s~=0));    % set the preference value as a common value of median of s
 init_pref = min(s(s~=0));    % set the preference value as a common value of median of s
 init_pref = init_pref * log(N);
 
 %%
 [idx,netsim,dpsim,expref] = apcluster(s,init_pref);
 % [idx,netsim,dpsim,expref] = apclusterSparse(s,init_pref);
 
 %%
 I = unique(idx);
 c = idx;
 K = length(I);
 for i=1:K
     c(idx==I(i)) = i;
 end
 
 %%
 color_idx = 'yrgbmc';
 marker_idx = 'ox+*dv^<>ph';
 
 figure(100)
 clf
 for i=1:K
     plot(data(c==i,1),data(c==i,2),'Color',color_idx(mod(i,6)+1),'Marker',marker_idx(mod(i,11)+1));
     hold on
     plot(data(I(i),1),data(I(i),2),'Color','k','Marker','s', 'MarkerSize',10);
 end

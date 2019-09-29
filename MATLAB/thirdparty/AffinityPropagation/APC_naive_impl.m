clc
 close all
 clear all
 
 %%
 % data generation
 % N = 20;
 % data = GenerateData(1,N);
 
 data = randgen([4 5],[9 0;0 9],[10 18],[6 0;0 6],[15 8],[4 0;0 4]);
 N = size(data,1);
 
 %%
 % Parameter setting
 N_iter = 100;
 lambda = 0.9;
 
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
 for i=1:N
     s(i,i) = init_pref;
 end
 
 %%
 a = zeros(N,N);
 r = zeros(N,N);
 
 for iter = 1:N_iter   
     % responsibility r(i,k)
     %   : sent from data point i to candidate exemplar point k
     %   : reflects the accumulated evidence for how well-suited point k is to
     %   serve as the exemplar for point i
     r_prev = r;
     a_plus_s = a+s;
     for i = 1:N
         for k = 1:N
             r(i,k) = s(i,k) - max(a_plus_s(i,[1:k-1,k+1:N]));
         end
     end
     %%
     r = (1-lambda)*r + lambda*r_prev;        
     
     %%
     % availability a(i,k)
     %   : sent from candidate exemplar point k to point i.
     %   : reflects the accumulated evidence for how appropriate it would be for
     %   point i to choose point k as its exemplar
     a_prev = a;
     r_over_zero = zeros(N,N);
     r_over_zero_idx = find(r>0);
     r_over_zero(r_over_zero_idx) = r(r_over_zero_idx);
     
     for i=1:N
         for k=1:N
             r_over_zero_buf = r_over_zero;
             r_over_zero_buf(i,k) = 0;
             r_over_zero_buf(k,k) = 0;
             sum_buf = sum(r_over_zero_buf(:,k));
             
             if i == k
                 a(i,k) = sum_buf;
             else
                 a(i,k) = min(0,r(k,k)+sum_buf);
             end
         end
     end    
     %%
     a = (1-lambda)*a + lambda*a_prev;
 
     %%
     E = r+a;
     I = find(diag(E)>0);
     K = length(I);  % Indices of exemplars
     [tmp c] = max(s(:,I),[],2);
     c(I) = 1:K;
     idx = I(c);
     
     %%
     color_idx = 'yrgbmc';
     marker_idx = 'ox+*dv^<>ph';
     
     figure(103)
     clf
     for i=1:K
         plot(data(c==i,1),data(c==i,2),'Color',color_idx(mod(i,6)+1),'Marker',marker_idx(mod(i,11)+1));
         hold on        
         plot(data(I(i),1),data(I(i),2),'Color','k','Marker','s', 'MarkerSize',10);
     end
     pause(0.1);
     
     %%
     [iter, K]
 end
 
 % %%
 % E = r+a;
 % I = find(diag(E)>0);
 % K = length(I)  % Indices of exemplars
 % [tmp c] = max(s(:,I),[],2);
 % c(I) = 1:K;
 % idx = I(c);
 % 
 % %%
 % color_idx = 'rgbkmcyrgbkcmy';
 % figure(100)
 % for i=1:K
 %     plot(data(c==i,1),data(c==i,2),'Color',color_idx(i),'Marker','o');
 %     hold on
 % end
 

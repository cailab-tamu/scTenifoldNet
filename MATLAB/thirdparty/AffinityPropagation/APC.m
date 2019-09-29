S = s;
 
 N=size(S,1); 
 A=zeros(N,N); 
 R=zeros(N,N); % Initialize messages
 
 S=S+1e-12*randn(N,N)*(max(S(:))-min(S(:))); % Remove degeneracies
 lam=0.9; % Set damping factor
 
 for iter=1:300
     % Compute responsibilities
     Rold=R;
     AS=A+S; 
     [Y,I]=max(AS,[],2);
     
     for i=1:N 
         AS(i,I(i))=-realmax; 
     end;
     
     [Y2,I2]=max(AS,[],2);
     R=S-repmat(Y,[1,N]);
     
     for i=1:N 
         R(i,I(i))=S(i,I(i))-Y2(i); 
     end;
     
     R=(1-lam)*R+lam*Rold; % Dampen responsibilities
     
     % Compute availabilities
     Aold=A;
     Rp=max(R,0); 
     for k=1:N 
         Rp(k,k)=R(k,k); 
     end;
     A=repmat(sum(Rp,1),[N,1])-Rp;
     dA=diag(A); A=min(A,0); 
     
     for k=1:N 
         A(k,k)=dA(k); 
     end;
     A=(1-lam)*A+lam*Aold; % Dampen availabilities
     
     %%
     E = R+A;
     I = find(diag(E)>0);
     K = length(I);  % Indices of exemplars
     [tmp c] = max(s(:,I),[],2);
     c(I) = 1:K;
     idx = I(c);    
     
     color_idx = 'yrgbmc';
     marker_idx = 'ox+*dv^<>ph';
     
     figure(105)
     clf
     for i=1:K
         plot(data(c==i,1),data(c==i,2),'Color',color_idx(mod(i,6)+1),'Marker',marker_idx(mod(i,11)+1));
         hold on        
         plot(data(I(i),1),data(I(i),2),'Color','k','Marker','s', 'MarkerSize',10);
     end
     pause(0.05);    
     
     [iter, K]    
 end;
 
 E=R+A; % Pseudomarginals
 I=find(diag(E)>0); 
 K=length(I); % Indices of exemplars
 [tmp c]=max(S(:,I),[],2); 
 c(I)=1:K; idx=I(c); % Assignments

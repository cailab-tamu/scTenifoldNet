function [labels,NC,labelid,i,Sp,Slam,NCfix,netsim,dpsim,expref,idx] = ...
    adapt_apcluster(data,dtype,pvalues,folds,adapt,varargin)

adapt = adapt+1;            % 0: turning to original AP
if adapt < 2
if strcmp(dtype, 'euclidean') || (isnumeric(dtype) && dtype == 1)
   [Dist, dmax] = similarity_euclid(data,2);
elseif strcmp(dtype, 'correlation') || (isnumeric(dtype) && dtype == 2)
   Dist = 1-(1+similarity_pearson(data'))/2;
   dmax = 1;
end

nrow = size(Dist,1);
nap = nrow*nrow-nrow;
s = zeros(nap,3);
j=1;

for i=1:nrow
   for k = [1:i-1,i+1:nrow]
     s(j,1) = i;
     s(j,2) = k; 
     s(j,3) = -Dist(i,k);
     j = j+1;
   end;
 end;
Dist = [];
else
  s = data;
  data =[];
  %pmin = -min(s(:,3));
end

pfixed = 0;
dn = find(s(:,3)>-realmax);
pmedian = median(s(dn,3));                     % Set preference to median similarity
pstep = folds*pmedian;                        % decreasing step of p
if length(pvalues) < 1
    pvalues = pmedian*0.5;
else
    pfixed = 1;
end

% Handle arguments to function
if nargin<2
    error('Too few input arguments');
else
    maxits=500; convits=50; lam=0.5; plt=0; details=0; nonoise=0;
    i=1;
    while i<=length(varargin)
        if strcmp(varargin{i},'plot')
            plt=1; i=i+1;
        elseif strcmp(varargin{i},'details')
            details=1; i=i+1;
		elseif strcmp(varargin{i},'sparse')
			[idx,netsim,dpsim,expref]=apcluster_sparse(s,pvalues,varargin{:});
			return;
        elseif strcmp(varargin{i},'nonoise')
            nonoise=1; i=i+1;
        elseif strcmp(varargin{i},'maxits')
            maxits=varargin{i+1};
            i=i+2;
            if maxits<=0 
                error('maxits must be a positive integer'); 
            end;
        elseif strcmp(varargin{i},'convits')
            convits=varargin{i+1};
            i=i+2;
            if convits<=0
                error('convits must be a positive integer'); 
            end;
        elseif strcmp(varargin{i},'dampfact')
            lam=varargin{i+1};
            i=i+2;
            if (lam<0.5)||(lam>=1)
                error('dampfact must be >= 0.5 and < 1');
            end;
        else i=i+1;
        end;
    end;
end;

if lam>0.9
    fprintf('\n*** Warning: Large damping factor in use. Turn on plotting\n');
    fprintf('    to monitor the net similarity. The algorithm will\n');
    fprintf('    change decisions slowly, so consider using a larger value\n');
    fprintf('    of convits.\n\n');
end;

% Check that standard arguments are consistent in size
if length(size(s))~=2 
    error('s should be a 2D matrix');
elseif length(size(pvalues))>2 
    error('pvalues should be a vector or a scalar');
elseif size(s,2)==3
    tmp=max(max(s(:,1)),max(s(:,2)));
    if length(pvalues)==1 
        N=tmp; 
    else
        N=length(pvalues);
    end;
    if tmp>N
        error('data point index exceeds number of data points');
    elseif min(min(s(:,1)),min(s(:,2)))<=0
        error('data point indices must be >= 1');
    end;
elseif size(s,1)==size(s,2)
    N=size(s,1);
    if (length(pvalues)~=N) && (length(pvalues)~=1)
        error('pvalues should be scalar or a vector of size N');
    end;
else error('s must have 3 columns or be square'); 
end;

% Construct similarity matrix
if N>3000
    fprintf('\n*** Warning: Large memory request. Consider activating\n');
    fprintf('    the sparse version of APCLUSTER.\n\n');
end;
if size(s,2)==3
    S=-Inf*ones(N,N); 
    for j=1:size(s,1)
        S(s(j,1),s(j,2))=s(j,3); 
    end;
else
    S=s;
end;
s = [];

% In case user did not remove degeneracies from the input similarities,
% avoid degenerate solutions by adding a small amount of noise to the
% input similarities
if ~nonoise
    rns=randn('state'); 
    randn('state',0);
    S=S+(eps*S+realmin*100).*rand(N,N);
    randn('state',rns);
end;

% Place preferences on the diagonal of S
if length(pvalues)==1 
    for i=1:N 
        S(i,i)=pvalues; 
    end;
else
    for i=1:N 
        S(i,i)=pvalues(i); 
    end;
end;

% Allocate space for messages, etc
dS=diag(S); 
A=zeros(N,N); 
R=zeros(N,N); 
t=1;
if plt 
    netsim=zeros(1,maxits+1); 
end;
if details
    idx=zeros(N,maxits+1);
    netsim=zeros(1,maxits+1); 
    dpsim=zeros(1,maxits+1); 
    expref=zeros(1,maxits+1); 
end;



% Initialization
dn=0; i=0;
stoptimes = max([maxits/10 2000]); % Stop if K's unchanging times >= stoptimes
if pfixed
    stoptimes = convits;
end
Hstop = zeros(N,stoptimes);  % recording unchanging times for stop condition
Hconvits = zeros(N,convits); % recording convergence at each K
Tdelay = 10;                            % delay time for detecting variation of K
Hdelay = Tdelay;                    % counting delay
Hconverg = 0;                      % whether K examplars convergence
nhalf = round(0.3*convits);
Hdelay2 = Tdelay;
Hconvhalf = zeros(N,nhalf);  % recording K at half convergence
Hsavehalf = 0;
Hn1 = 0; Hn2 = 0;
Kmean=0; Kdown=0; 
Wstart = max([100 round(convits/2)]); % start oscillation monitoring
wsize = 40;                  % monitoring window size for detecting oscillations
Kunchange=0;                 % recording variation of K
Kocil = ones(1,wsize);              % recording K's non-oscillations
Noscil = wsize+10;               % statistic of non-oscillations, initial value
Svib = 0;                               % counting potential oscillations
Hvib = 10; Tvib = 2;                    % counting oscillations, initial values
Hguid = 1;                            % label to start solution recording
Sprefer = pvalues;                     % asigned to diag of S for realizing p decreasing
astep = pstep;                      % adjustable step instead of fixed pstep
Kset = []; Kold = 0;              % recording K
Kfix = 0; nKfix = 0;             % used for speeding up
Kmax = 0; nfix = 0;

% Execute parallel affinity propagation updates
while ~dn
    i=i+1;

    % Compute responsibilities
    AS=A+S; 
    [Y,I]=max(AS,[],2); 
    for k=1:N 
        AS(k,I(k))= -realmax;
    end;
    [Y2,I2]=max(AS,[],2);
    AS = [];
    Rold=R;
    R=S-repmat(Y,[1,N]);
    for k=1:N 
        R(k,I(k))=S(k,I(k))-Y2(k); 
    end;
%     R=(1-lam)*R+lam*Rold;
    Rold = lam*Rold;
    R = (1-lam)*R;
    R = R+Rold;
    Rold = [];

    % Compute availabilities
    Rp=max(R,0); 
    for k=1:N 
        Rp(k,k)=R(k,k); 
    end; 
    Aold=A;
%     A=repmat(sum(Rp,1),[N,1])-Rp;
    A = sum(Rp,1);
    A = repmat(A,[N,1]);
    A = A-Rp;
    Rp = [];
    dA=diag(A);
    A=min(A,0);
    for k=1:N
        A(k,k)=dA(k);
    end;
%     A=(1-lam)*A+lam*Aold;
    Aold = lam*Aold;
    A = (1-lam)*A;
    A = A+Aold;
    Aold = [];

    % Check for convergence
    E = ((diag(A)+diag(R))>0);
    Hconvits(:,mod(i-1,convits)+1) = E;
    K=sum(E);
    Kset(i) = K;
    newp = Sprefer(1);
    newlam = lam;
    Hstop(:,mod(i-1,stoptimes)+1) = E;
    Hconvhalf(:,mod(i-1,nhalf)+1) = E;
    if mod(i,100) == 1 || i == maxits
        fprintf('** running at iteration %d, K = %d\n', i,K);
    end
    
    Hsave = 0; Hsave1 = 0; Hsave2 = 0; Hsave3 = 0;
    if i>=Wstart || i>=maxits
        se = sum(Hconvits,2);
        se1 = sum(se==convits);
        se2 = sum(se==0);
        unconverged = (se1+se2) ~= N;
        Hconverg = ~unconverged; 
        se = sum(Hstop,2);
        se1 = sum(se==stoptimes);
        se2 = sum(se==0);
        if (se1+se2) == N || i == maxits
            dn=1;
            if (se1+se2) == N
               Hsave1 = 1;
            end
        end
        se = sum(Hconvhalf,2);
        se1 = sum(se==nhalf);
        se2 = sum(se==0);
        Hsavehalf = (se1+se2) == N;
        Hsavehalf = Hsavehalf && Hguid == 2;
    end
    
 if adapt   
    if i > 5
      Kmean(i) = mean(Kset(i-5:i));        % covering 6 points
      Kdown(i) = Kmean(i)-Kmean(i-1) < 0;  % covering 7 points
      if Hguid == 2
          Kdown(i) = Kdown(i) && K <= Kold;
      end
      Kunchange(i) = sum(abs(Kset(i)-Kset(i-5:i-1)));
      Kocil(:,mod(i-1,wsize)+1) = Kdown(i) || Kunchange(i) == 0;
      Noscil = sum(Kocil);                     % non- oscillations
    end
    
    % reducing parameter pvalues to yield smaller NC when unchanging
    if Hconverg
      Hdelay = Hdelay+1;
      if Hdelay >= Tdelay
        Hsave1 = 1;        % starting pvalues reduction & result saving
        Hdelay = 0;
        Hn1 = Hn1+1;
        if K == Kfix
            nKfix = nKfix+1;
        else
            nKfix = 0;
        end
        Kfix = K;        
        stepfold = sqrt(K+50)/10;
        pstep = folds*pmedian/stepfold;
        if nKfix > 1
            astep = nKfix*pstep;                 % speeding up
        else
            astep = pstep;
        end
      end
    elseif Hsavehalf
        Hdelay2 = Hdelay2+1;
        if Hdelay2 >= Tdelay
            Hsave2 = 1;        % starting pvalues reduction & result saving
            Hdelay2 = 0;
            Hn2 = Hn2+1;
        end
    end
    if ~Hconverg
        Hn1 = 0;
        Hdelay = 0;
    end
    if ~Hsavehalf
        Hn2 = 0;
        Hdelay2 = 0;
    end
 
    if ((K == 1 || K == 2) && Hsave1)   % avoiding influence of oscillation
       dn = 1;
       unconverged = 0;
    end

    if Hguid == 1 && Hsave1
        Hguid = 2;                               % starting guidance
        labels = zeros(N,K);
        labelid = zeros(N,K);
        NC = zeros(1,K);
        NCfix = zeros(1,K);
        Sp = zeros(1,K);
        Slam = zeros(1,K);
        Kmax = K;
        stepfold = sqrt(Kmax+50)/10;
        pstep = folds*pmedian/stepfold;
    end
      
    if Hsave1                                     % K >= Kold unchanging, K rise
      Svib = 0;
      if ~pfixed
          Sprefer = Sprefer+astep;           % reducing pvalues
          if length(pvalues)==1
              for k=1:N
                  S(k,k) = Sprefer;
              end
          else
              for k=1:N
                  S(k,k) = Sprefer(k);
              end
          end
      end
      
    else
       Svib = Svib+1;                             % counting potential oscillations
       HSvib = (Svib > wsize && Noscil < 0.66*wsize) || Svib > 150;
       HSvib = HSvib && i > Wstart;
       if HSvib                          % = & ! should be 2/3, otherwise oscillations
          Hvib = Hvib+1;
          if Hvib > 10
            lam = max([0.7 lam]);
          elseif Hvib >= 1
            if Tvib >= 3
                if lam >= 0.9
                    lam = min([0.98 0.025+lam]);
                    if lam >= 0.95 && mod(i,9) == 2
                        rns=randn('state');
                        randn('state',0);
                        S=S+(eps*S+realmin*1000).*rand(N,N);
                        randn('state',rns);
                        fprintf(' # A small amount of noise is added\n');
                    end
                end
            else
            lam = min([0.9 0.05+lam]);
            end
            if lam >= 0.85
                Tvib = Tvib+1;
                if ~pfixed
                    if Hguid == 2 && Kold
                       if Kmax <1
                          stepfold = 2;
                       else
                          stepfold = max([3/(sqrt(Kmax)/10+0.4) 1]);
                       end
                       Kvar = 2*sqrt(std(Kset(i-49:i)));
                       astep = min([0.8*Kvar+0.2*Tvib stepfold])*pstep;
                       %0.5*((Kvar+min([Tvib 4]))*pstep); max([sqrt(abs(K-Kold))])
                    else
                       astep = min([Tvib 2])*pstep;
                    end
                    Sprefer = Sprefer+astep;  % escaping oscillations
                    if length(pvalues)==1
                        for k=1:N
                            S(k,k) = Sprefer;
                        end
                    else
                        for k=1:N
                            S(k,k) = Sprefer(k);
                        end
                    end
                    fprintf(' # Escaping oscillation turns on\n');
                end
            end
          end
          Hvib = 0;
          Svib = 0;
          fprintf(' # Damping factor is increased to %g\n', lam);
       else
           Tvib = max([Tvib-0.002 0.98]);
           if lam > 0.9 && Tvib < 1
              lam = max([lam-0.0001 0.5]);
           end
       end
       
    end
 end
 
   if Hguid >= 2 && (K < Kold && K > 1 || Kmean(i) == Kmean(i-1))
      Hsave3 = 1;              % catching decreasing K
   end
   if Hsave1 || Hsave2 || Hsave3
       Hsave = 1;
   end
   Kold = K;

    % Handle plotting and storage of details, if requested
    if plt || details || Hsave
        if K==0
            tmpnetsim=nan; tmpdpsim=nan; tmpexpref=nan; tmpidx=nan;
        else
            I=find(E); 
            [tmp c]=max(S(:,I),[],2); 
            c(I)=1:K;
            if Hsave || dn == 1
               if Hsave1
                  nfix = Hn1*Tdelay+convits;
               elseif Hsave2
                  nfix = Hn2*Tdelay+nhalf;
               elseif Kmean(i) == Kmean(i-1)
                  nfix = 6;
                  if Kmean(i) == Kmean(i-5)
                    nfix = 10;
                  end
               else
                  nfix = 1;
               end
               
               if (K <= Kmax && nfix > NCfix(K)) || (K > Kmax && nfix >= 10)
                   NCfix(K) = nfix;
                   labels(:,K) = c;
                   labelid(:,K) = I(c);
                   NC(K) = K;
                   Sp(:,K) = newp; 
                   Slam(:,K) = newlam;
               end
               if K > Kmax
                   Kmax = K;
                   if length(NCfix) < K
                       NCfix(K) = 0;
                   end
               end
            else
              tmpidx=I(c);
              tmpnetsim=sum(S((tmpidx-1)*N+[1:N]'));
              tmpexpref=sum(dS(I)); 
              tmpdpsim=tmpnetsim-tmpexpref;
            end
        end
    end;
    if details
        netsim(i)=tmpnetsim; dpsim(i)=tmpdpsim; expref(i)=tmpexpref;
        idx(:,i)=tmpidx;
    end;
    if plt
        netsim(i)=tmpnetsim;
        figure(1); 
        tmp=1:i; 
        tmpi=find(~isnan(netsim(1:i)));
        subplot(2,1,1)
        plot(tmp(tmpi),netsim(tmpi),'r-');
        xlabel('# Iterations');
        ylabel('Fitness (net similarity)');
        subplot(2,1,2)
        tmp=max([1 i-3*wsize]):i; 
        plot(tmp,Kset(tmp),'b-');
        xlabel('# Iterations');
        ylabel('Number of clusters');
        drawnow; 
    end;
end;


fprintf(' # Programs run over at K= %d\n', K);
I=find(diag(A+R)>0); 
K=length(I);                                      % Identify exemplars
if K>0
    [tmp c]=max(S(:,I),[],2); 
    c(I)=1:K;                                       % Identify clusters
    % Refine the final set of exemplars and clusters and return results
    for k=1:K 
        ii=find(c==k);
        [y j]=max(sum(S(ii,ii),1)); 
        I(k)=ii(j(1)); 
    end
    [tmp c]=max(S(:,I),[],2); 
    c(I)=1:K; 
    tmpidx=I(c);
    tmpnetsim=sum(S((tmpidx-1)*N+[1:N]')); 
    tmpexpref=sum(dS(I));
    labels(:,K) = c;
    labelid(:,K) = tmpidx;
    NC(K) = K;
    NCfix(K) = nfix;
    Sp(:,K) = newp;
    Slam(:,K) = newlam;
else
    tmpidx=nan*ones(N,1); 
    tmpnetsim=nan; 
    tmpexpref=nan;
end;
if details
    netsim(i+1)=tmpnetsim; netsim=netsim(1:i+1);
    dpsim(i+1)=tmpnetsim-tmpexpref; dpsim=dpsim(1:i+1);
    expref(i+1)=tmpexpref; expref=expref(1:i+1);
    idx(:,i+1)=tmpidx; idx=idx(:,1:i+1);
else
    netsim=tmpnetsim; 
    dpsim=tmpnetsim-tmpexpref;
    expref=tmpexpref; 
    idx=tmpidx;
end;

if length(NC) > 1
   NC(1) = 0;
end
if length(NC) < 1
  Sp = [];
  Slam = [];
  NC = 0;
  NCfix = 0;
  labels = ones(nrow,1);
  labelid = ones(nrow,1);
else
  S = find(NC);
  Sp = Sp(:,S);
  Slam = Slam(:,S);
  NC = NC(:,S);
  labels = labels(:,S);
  labelid = labelid(:,S);
  NCfix = NCfix(:,S);
end

if plt || details
    fprintf('\nNumber of identified clusters: %d\n',K);
    fprintf('Fitness (net similarity): %f\n',tmpnetsim);
    fprintf('  Similarities of data points to exemplars: %f\n',dpsim(end));
    fprintf('  Preferences of selected exemplars: %f\n',tmpexpref);
    fprintf('Number of iterations: %d\n\n',i);
end;
if unconverged
    fprintf('\n*** Warning: Algorithm did not converge at K = %g ! \n', NC(1));
    fprintf('    The similarities may contain degeneracies - add noise to\n');
    fprintf('    the similarities to remove degeneracies. To monitor the net\n');
    fprintf('    similarity, activate plotting. Also, consider increasing\n');
    fprintf('    maxits and if necessary dampfact.\n\n');
end;

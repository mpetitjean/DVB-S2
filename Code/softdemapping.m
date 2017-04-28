function cHat = softdemapping(y,H,maxit)

y = reshape(y,1,size(H,2),[]);
H = repmat(full(H),[1 1 size(y,3)]);
[M,N,R] = size(H);
subrow = repmat(1:M,1,R);
sub3thd = repmat(1:R,M,1);
sub3thd = sub3thd(:).';
Lqij = 2*y.*H;
Lci = -y;
it = 0;
while it < maxit
    alphaij = sign(Lqij);
    alphaij(alphaij==0)=1;
    betaij = full(abs(Lqij));
    betaij(betaij==0) = inf;
    [val,mincols] = min(betaij,[],2);
    mincols = mincols(:).';
    % generate matrix made of minimum value of each row
    minbetaij = repmat(val,[1 N 1]);
    % find indexes of minimum values
    minidxs = sub2ind([M N R],subrow,mincols,sub3thd);
    % replace minimum values with infs
    betaij(minidxs) = inf;
    % find next minimum values
    val = min(betaij,[],2);
    % set original minimum elements to next minimum values
    minbetaij(minidxs) = val;
    Lrji = prod(alphaij,2).*alphaij.*minbetaij.*H;
    colsum = sum(Lrji);
    Lqij = (Lci + colsum- Lrji).* H;
    LQi = Lci + colsum;
    cHat = LQi < 0;
    it = it + 1;
end
cHat = cHat(:,M+1:end,:);
cHat = cHat(:).';
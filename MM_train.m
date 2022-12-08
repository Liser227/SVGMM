function model = MM_train(X,  prior, model, setting)
%%
% X is PolSAR data, its dimention is: N x D^2
% N : number of all samples
% D : element number of your scattering vector
% L : number of looks
% K : component number
tol             = setting.tol;
maxitr          = setting.maxitr; 
%% load date
t               = 0;
Converged       = 0;          
llh_old         = inf;       
llh_itr         = [];

while ~Converged && t< maxitr   
    t           = t+1;
    model       = VBM(X,  model, prior, setting);
    model       = VBE(X,  model, setting);
    
    llh         = LB(X, model, prior, setting);
    llh_new     = llh;
    gap         = abs((llh_new-llh_old)./llh_old);
    R           = model.R;
    N_cla       = sum(R,1)./sum(sum(R));
    Converged   = gap<tol;
         
    llh_itr     = [llh_itr,llh_new]  ;
    llh_old     = llh_new;
    
end
model.llh_itr   = llh_itr;
if Converged
    fprintf('Converged in %d M-steps.\n',t-1);
else
    fprintf('Not converged in %d M-steps.\n',maxitr);
end
end


function model = VBM(X, model,prior,setting)
%%
% W = [1 2 3 4 5 6 7 8 9]';
%%
lnDetX          = model.lnDetX;
X               = model.X;

% Prior.
Beta0          = prior. Beta0;          
Alpha_u0       = prior.Alpha_u0 * ones(1,K);        
Alpha_v0       = prior.Alpha_v0 * ones(1,K);
Gama_s0        = prior.Gama_s0 * ones(1,K);
Gama_k0        = prior.Gama_k0 * ones(1,K);
eta0            = prior.eta0* ones(1,K);           
W0             = prior.W0;

% Intermedia Equations
Ez              = model.R;
EAlpha          = bsxfun(@rdivide,model.Alpha_v, model.Alpha_u);  

EGama           = model.Gama_k./model.Gama_s;
ElnAlpha        = psi(0,model.Alpha_v)-log(model.Alpha_u);
ElnGama         = psi(0,model.Gama_k)-log(model.Gama_s);

L_Tr_Gama       = bsxfun(@plus,bsxfun(@times,real(X*model.ELambda),L),EGama); % N x K

%% Update
% for q(PI)
Beta            = Beta0 + Ez;

eta             = eta0 + L.*sum(Ez,1);
tmp1            = bsxfun(@times, L, X.'*Ez);
tmp2            = bsxfun(@rdivide,bsxfun(@times,bsxfun(@plus,EAlpha,L.*D),tmp1),sum(L_Tr_Gama,1));   
W               = bsxfun(@rdivide,bsxfun(@plus, bsxfun(@times,eta0, W0)...
                                                   ,tmp2) ,eta);
DetW                = zeros(1,K);
ELambda             = zeros(D^2, K);
for i_k = 1:K 
    W_tmp           = reshape(W(:,i_k), D,[]).';
    
    DetW(1,i_k)     = real(det(W_tmp));
    ELambda(:,i_k)  = reshape(inv(W_tmp), [], 1);
end                                 

% for q(Alpha)
Alpha_v             =  bsxfun(@plus,Alpha_v0 ,bsxfun(@times,sum(Ez,1),D.*L));
tmp1                =  bsxfun(@minus, bsxfun(@minus,psi(0,bsxfun(@plus,(D.*L),EAlpha)),psi(0,EAlpha)),(D.*L).*(1./EAlpha));   
Alpha_u             =  bsxfun(@minus,Alpha_u0,      sum(bsxfun(@times,Ez,bsxfun(@minus,tmp1+ElnGama,log(L_Tr_Gama))),1));  

% for q(Gama)
Gama_k              =  bsxfun(@plus,Gama_k0 , bsxfun(@times, sum(Ez,1),EAlpha));
Gama_s              =  bsxfun(@plus,Gama_s0, sum(bsxfun(@times,bsxfun(@plus,EAlpha, D.*L),bsxfun(@rdivide,Ez,L_Tr_Gama)),1) );  
                               
%% Output
model.Beta          = Beta;
model.Alpha_v       = Alpha_v;
model.Alpha_u       = Alpha_u;
model.Gama_k        = Gama_k;
model.Gama_s        = Gama_s;
model.eta           = eta;
model.ELambda       = ELambda;
model.DetW          = DetW;
model.lnDetX        = lnDetX;
model.L_Tr_Gama=L_Tr_Gama;
end


function model = VBE(X,  model,setting)
%%
lnDetX          = model.lnDetX;
ELambda         = model.ELambda;
ElnDetLambda    = -log(model.DetW)-D.*log(model.eta)+MulPolyGamma(model.eta,D,0);
EAlpha          = bsxfun(@rdivide,model.Alpha_v,model.Alpha_u);
ElnAlpha        = psi(0,model.Alpha_v)-log(model.Alpha_u);
ElnGama         = psi(0,model.Gama_k)-log(model.Gama_s);
EGama           = bsxfun(@rdivide,model.Gama_k,model.Gama_s);
ElnPI           = bsxfun(@minus,psi(0,model.Beta),psi(0,sum(model.Beta,2)));
%% Update Ez
tmp1            = bsxfun(@plus,(L.*D.*log(L) - lnMulGamma(L, D)),bsxfun(@times, L-D, model.lnDetX)) ;
tmp2            = log(gamma(EAlpha+D.*L))-log(gamma(EAlpha))-bsxfun(@times,D.*L,log(EAlpha))+bsxfun(@times,D.*L,ElnAlpha);
tmp3            = bsxfun(@times,EAlpha,ElnGama);
L_Tr_Gama       = model.L_Tr_Gama;
tmp4            = bsxfun(@plus,-bsxfun(@times,bsxfun(@plus,EAlpha,D.*L),log(L_Tr_Gama)),L.*ElnDetLambda);  
lnRho_tmp       = bsxfun(@plus,tmp1,tmp2+tmp3)+tmp4+ElnPI;
maxlnrho        = max(lnRho_tmp,[],2);
lnRho           = bsxfun(@minus, lnRho_tmp, maxlnrho);
rho             = exp(lnRho);
R               = bsxfun(@rdivide, rho, sum(rho,2));

%% Output
model.R        = R;
model.EAlpha   = EAlpha;
model.ElnAlpha = ElnAlpha;
model.ElnGama  = ElnGama;
model.EGama    = EGama;
model.ElnPI    = ElnPI;
model.lnDetX   = lnDetX;
model.ELambda  = ELambda;

model.L_Tr_Gama= L_Tr_Gama;
end
function llh = LB(X, model,prior,setting)
%%
L               = model.L;          % for scaled Complex Wishart Component
N               = size(X,1);
D               = model.D;
K               = size(L,2);
Beta          = model.Beta;
Alpha_v       = model.Alpha_v;
Alpha_u       = model.Alpha_u;
Gama_k        = model.Gama_k;
Gama_s        = model.Gama_s;

eta           = model.eta;
ELambda       = model.ELambda;
EAlpha        = model.EAlpha;
ElnAlpha      = model.ElnAlpha;
ElnGama       = model.ElnGama;
EGama         = model.EGama;
ElnPI         = model.ElnPI;
Ez            = model.R;


Beta0              = prior.Beta0;      
Alpha_u0            = prior.Alpha_u0 ;    
Alpha_v0            = prior.Alpha_v0; 
Gama_s0             = prior.Gama_s0;      
Gama_k0             = prior.Gama_k0; 
eta0                = prior.eta0.*ones(1,size(L,2));
W0                  = prior.W0;
lnDetX              = model.lnDetX;
DetW                = model.DetW;
ElnDetLambda        = -log(DetW)-D.*log(eta)+MulPolyGamma(eta,D,0);

%% Update
Tr              = bsxfun(@plus,bsxfun(@times,L,real(X*ELambda)),EGama);
tmp1           = bsxfun(@plus,L.*D.*log(L),log(gamma(L.*D+EAlpha)));
tmp2           = bsxfun(@minus, bsxfun(@minus,bsxfun(@times,EAlpha,ElnGama),log(gamma(EAlpha))) ,lnMulGamma(L, D));
tmp3           = bsxfun(@times,-EAlpha-L.*D, log(Tr));
tmp4           = bsxfun(@plus,bsxfun(@times, L-D, lnDetX),L.*ElnDetLambda);   %2021.4.16
tmp             = bsxfun(@plus,tmp3,tmp1+tmp2)+tmp4;
lnPX            = sum(sum(Ez.*tmp,1),2);

lnPZ            = sum(sum(bsxfun(@times, Ez, ElnPI),1),2);
tmp_1           = gammaln(K.* Beta0)- K.*gammaln(Beta0) + sum((Beta0-1).*ElnPI,2);
lnPPI           = sum(tmp_1, 1);

tmp_1            = bsxfun(@times, Alpha_v0, log(Alpha_u0))-gammaln(Alpha_v0);
tmp_2            = bsxfun(@times, Alpha_v0-1,ElnAlpha) - bsxfun(@times,Alpha_u0,EAlpha);
lnPAlpha         = sum(tmp_1+tmp_2,2);

tmp_1            = bsxfun(@times, Gama_k0, log(Gama_s0))-gammaln(Gama_k0);
tmp_2            = bsxfun(@times, Gama_k0-1,ElnGama) - bsxfun(@times,Gama_s0,EGama);
lnPGama          = sum(tmp_1+tmp_2,2);

lnDetW0         = log(real(det(reshape(W0,D,D) )));
Tr              = W0.'*ELambda;
tmp             = eta0.*D.*log(eta0)-lnMulGamma(eta0, D) + (eta0-D).*ElnDetLambda...
                  +eta0.*lnDetW0 - eta0.*Tr;
lnPLambda       = sum(tmp,2);

lnqZ            = sum(sum(log(Ez.^Ez),1),2);

tmp_1           = gammaln(sum(Beta,2))- sum(gammaln(Beta),2) + sum((Beta-1).*ElnPI,2);
lnqPI           = sum(tmp_1,1);

tmp_1            = bsxfun(@times, Alpha_v, log(Alpha_u))-gammaln(Alpha_v);
tmp_2            = bsxfun(@times, Alpha_v-1,ElnAlpha) - bsxfun(@times,Alpha_u,EAlpha);
lnqAlpha         = sum(tmp_1+tmp_2,2);

tmp_1            = bsxfun(@times, Gama_k, log(Gama_s))-gammaln(Gama_k);
tmp_2            = bsxfun(@times, Gama_k-1,ElnGama) - bsxfun(@times,Gama_s,EGama);
lnqGama          = sum(tmp_1+tmp_2,2);

Tr              = D;  
tmp             = eta.*D.*log(eta)-lnMulGamma(eta, D) + (eta-D).*ElnDetLambda...
                  +eta.*log(DetW) - eta.*Tr;
lnqLambda       = sum(tmp,2);

%% Output
llh         = 1./(N*D).*(lnPX+lnPZ +lnPPI+lnPLambda+lnPAlpha+lnPGama - lnqZ-lnqPI-lnqLambda-lnqAlpha-lnqGama );

end

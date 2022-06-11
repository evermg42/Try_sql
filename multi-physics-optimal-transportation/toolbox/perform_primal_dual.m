function [x,R,y] = perform_primal_dual(x, K,  KS, ProxFS, ProxG, options)

% perform_admm - preconditionned ADMM method
%
%    [x,R,y] = perform_admm(x, K,  KS, ProxFS, ProxG, options);
%
%   Solves 
%       min_x F(K*x) + G(x)
%   where F and G are convex proper functions with an easy to compute proximal operator, 
%   and where K is a linear operator
%
%   Uses the Preconditioned Alternating direction method of multiplier (ADMM) method described in 
%       Antonin Chambolle, Thomas Pock, 
%       A first-order primal-dual algorithm for convex problems with applications to imaging, 
%       Preprint CMAP-685	
%
%   INPUTS:
%   ProxFS(y,sigma) computes Prox_{sigma*F^*}(y)
%   ProxG(x,tau) computes Prox_{tau*G}(x)
%   K(y) is a linear operator.
%   KS(y) compute K^*(y) the dual linear operator.
%   options.sigma and options.tau are the parameters of the
%       method, they shoudl satisfy sigma*tau*norm(K)^2<1
%   options.theta=1 for the ADMM, but can be set in [0,1].
%   options.verb=0 suppress display of progression.
%   options.niter is the number of iterations.
%   options.report(x) is a function to fill in R.
%
%   OUTPUTS:
%   x is the final solution.
%   R(i) = options.report(x) at iteration i.
%
%   Copyright (c) 2010 Gabriel Peyre





options.null = 0;
report       = getoptions(options, 'report', @(x)0);
niter        = getoptions(options, 'niter', 100);
theta        = getoptions(options, 'theta', 1.);
verb         = getoptions(options, 'verb', 1);
plot_hauteur = getoptions(options, 'plot_hauteur', 1);
if(isnumeric(K))
  K          = @(x) K*x;
end  
y            = getoptions(options, 'dualvar', K(x));
if(isnumeric(KS))
  KS         = @(x) KS*x;
end      

%%%% ADMM parameters %%%%
sigma        = getoptions(options, 'sigma', -1);
tau          = getoptions(options, 'tau', -1);
if(sigma<0 || tau<0)
  [L,e]      = compute_operator_norm(@(x)KS(K(x)),randn(size(x)));
  sigma      = 10;
  tau        = .9/(sigma*L);
end
%y            = K(x);
x1           = x;
clear R;
 


      

for i = 1:niter    
  if(verb)
    progressbar(i,niter);
  end
  % record energies
  tt = sum(size(x));
  if tt == 2
  R(i)       = report(x,y);
  else
  R = 0;
  end
  
  % update
  xold       = x;
  
  tt = sigma*K(x1);
  [sum(y(:,:,1,3),'all') sum(tt(:,:,1,3),'all')]
  y          = ProxFS(y + sigma*K(x1),sigma);
  
  out1 = [sum(y(:,:,1,3),'all') sum(y(:,:,2,3),'all')]
  x          = ProxG(x - tau*KS(y),tau);
  
  ttt = x - KS(y);
  out2 = [sum(ttt.M{3}(:,:,1),'all') sum(ttt.M{3}(:,:,2),'all')]
  x1         = x + theta*(x-xold);
% X=interp(x1);
% figure(3);
% subplot(2,2,1)
% quiver(X(:,:,1,2),X(:,:,1,1))
% subplot(2,2,2)
% quiver(X(:,:,floor(end/3),2),X(:,:,floor(end/3),1))
% subplot(2,2,3)
% quiver(X(:,:,floor(end*2/3),2),X(:,:,floor(end*2/3),1))
% subplot(2,2,4)
% quiver(X(:,:,end,2),X(:,:,end,1))
% figure(2);
% subplot(3,4,1)
% contour(X(:,:,1,3))
% subplot(3,4,2)
% contour(X(:,:,4,3))
% subplot(3,4,3)
% contour(X(:,:,7,3))
% subplot(3,4,4)
% contour(X(:,:,10,3))
% subplot(3,4,5)
% contour(X(:,:,13,3))
% subplot(3,4,6)
% contour(X(:,:,16,3))
% subplot(3,4,7)
% contour(X(:,:,19,3))
% subplot(3,4,8)
% contour(X(:,:,22,3))
% subplot(3,4,9)
% contour(X(:,:,25,3))
% subplot(3,4,10)
% contour(X(:,:,28,3))
% subplot(3,4,11)
% contour(X(:,:,30,3))
% subplot(3,4,12)
% contour(X(:,:,32,3))

%% To see the evolution of errors  

% Jlist  = s2v(R,'J');
% Constr = s2v(R,'Constr');
% MinVal = s2v(R,'Min');
% figure(1);
%   subplot(3,1,1);
% plot(Jlist(10:end), '.-'); axis tight;
% title('J');
% subplot(3,1,2);
% plot((Constr(10:end)), '.-'); axis tight;
% title('div=0 violation');
% subplot(3,1,3);
% plot((MinVal(10:end)), '.-'); axis tight;
% title('Minimum value of f');
% drawnow;
end

  

end
 
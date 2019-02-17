%
% Cahn-Hilliard integrator.  Solves the CH equation
%   with natural and no-flux boundary conditions.
%   Nonlinear term: f(u) = u - u^3
%
% Numerical scheme;
%
%    Eyre's linearly stabilized CH integration scheme
%
% Notes;
%
%    1. Uses a fixed time step.
%
%    2. Uses discrete cosine transform to invert
%       the linear update matrix by performing
%       a spectral decomposition of the matrix.
%
%    3. Calculates spinodal decomposition from
%       random initial conditions with mean value
%       close to zero.
%
%    4. With the given parameters, there are roughly
%       10 grid points in a transition layer.
%
%    5. The timestep probably should be shorted
%       initially and lengthed later for a more
%       accurate solution.
%       
%    6. The intermediate output is time and max(|U|)
%

clear
format long 

% spatial dimensions -- adjust N and M to increase decrease
% the size of the computed solution.

N = 128;  M = 128;
delx = 1/(M-1);
delx2 = delx^2;
x = (0:delx:1)';

% graphics parameters

visual_update = 10;
type_update = 10;

% time parameters -- adjust ntmax to take more time
% steps, and delt to take longer time steps.

t = 0;
delt = 0.00005;
ntmax = 250;

% CH parameters -- to see more layering, decrease
% epsilon, but be careful, you might need to increase
% the size of the grid.

epsilon = 0.01;
eps2 = epsilon^2;

% time-step parameter used in Eyre's scheme

a = 2;

% time marching update parameters

lam1 = delt/delx2;
lam2 = eps2*lam1/delx2;

% unscaled eigenvalues of the laplacian (nuemann bc)

Leig  = (((2*cos(pi*(0:N-1)'/(N-1)))-2)*ones(1,M)) + ...
         (ones(N,1)*((2*cos(pi*(0:M-1)/(M-1)))-2));

% scaled eigenvalues of stabilized CH update matrix

CHeig = ones(N,M) - (a*lam1*Leig) + (lam2*Leig.*Leig);

% scaled eigenvalues of the laplacian

Seig = lam1*Leig;

% random initial conditions

U = (rand(N,M)-0.5)*0.01;
hat_U = dct2(U);

% main loop

it = 0; j=0;
t = 0.0;
while it < ntmax

% intermediate output

  if rem(it,type_update) == 0
    [it t max(max(abs(U)))]
  end

% plotting and movie

  if rem(it,visual_update) == 0
    subplot(2,2,1)
    pcolor(U), shading interp, ...
    axis('off'), axis('equal'), title('U');
    subplot(2,2,2)
    mesh(U), axis([0 N 0 M -1 1]), title('U');
    vhat = abs(hat_U); vhat(1,1)=0; 
    vhat = vhat/max(max(abs(vhat))); 
    subplot(2,2,3)
    plot(x',U(M/2,:)), axis([0 1 -1 1]), ...
    title('cross section of U');
    subplot(2,2,4)
    pcolor(vhat), ...
    shading interp, axis('off'), axis('equal'),
    axis([0 M/4 0 N/4]), title('| Uhat |'), pause(1.0);
    if it == 0
      %mov = moviein(ntmax/visual_update,gcf);
    end
    j=j+1; %mov(:,j) = getframe(gcf);
  end

% Update the solution

  it = it+1;
  t = t+delt;

% compute the shifted nonlinear term

  fU = (U.*U.*U) - ((1+a)*U);

% compute the right hand side in tranform space

  hat_rhs = hat_U + (Seig.*dct2(fU));

% compute the updated solution in tranform space

  hat_U = hat_rhs./CHeig;

% invert the cosine transform

  U = idct2(hat_U);

end

% play the movie

%movie(gcf,mov,10,3)


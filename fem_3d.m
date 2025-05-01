tetra10coords={{'cylinder/coords12.csv', 'cylinder/coords13.csv', 'cylinder/coords14.csv'} ...
               {'sphere/coords12.csv', 'sphere/coords12.5.csv', 'sphere/coords12.75.csv'} ...
               {'cone/nosmooth/coords12.8.csv', 'cone/smooth/coords12.8.csv'}};
tetra10nodes={{'cylinder/nodes12.csv', 'cylinder/nodes13.csv', 'cylinder/nodes14.csv'} ...
              {'sphere/nodes12.csv', 'sphere/nodes12.5.csv', 'sphere/nodes12.75.csv'} ...
              {'cone/nosmooth/nodes12.8.csv', 'cone/smooth/nodes12.8.csv'}};
tetra10bdryn={{'cylinder/bdryn12.csv', 'cylinder/bdryn13.csv', 'cylinder/bdryn14.csv'} ...
              {'sphere/bdryn12.csv', 'sphere/bdryn12.5.csv', 'sphere/bdryn12.75.csv'} ...
              {'cone/nosmooth/bdryn12.8.csv', 'cone/smooth/bdryn12.8.csv'}};

ifplot=false;
% ifplot=true;

function u_vals=u_node_vals(nodes, u)
N = size(u, 1);
u_vals = zeros(N, 1);
node_count = zeros(N, 1);

for el = 1:size(nodes, 1)
    idx = nodes(el, :);
    u_local = u(idx);
    
    % Lagrange basis functions at triangle vertices = eye(10)
    % Contribution to each vertex is just the corresponding u value
    for i = 1:10
        u_vals(idx(i)) = u_vals(idx(i)) + u_local(i);
        node_count(idx(i)) = node_count(idx(i)) + 1;
    end
end

u_vals = u_vals ./ node_count;
end

function val = g(x, y, z)
    x0 = -7.1;
    y0 = -7.2;
    z0 = -7.3;
    x = x - x0;
    y = y - y0;
    z = z - z0;
    val = -1 * log(sqrt(((x * x) + (y * y) + (z * z))));
    % val = x + y + z;
end

gm = 1;
MM = 3;
u_h_eps = zeros(MM, 2);

for mm=1:MM
    [coords,nodes,bdryn,hmax]=loadmesh(tetra10coords{gm}{mm},tetra10nodes{gm}{mm}, ...
        tetra10bdryn{gm}{mm},ifplot);

    [bdryn_len,~] = size(bdryn);
    bdryv = zeros(bdryn_len, 1);
    for i=1:bdryn_len
        bd_crd = coords(bdryn(i),:);
        bdryv(i) = g(bd_crd(1), bd_crd(2), bd_crd(3));
    end
    
    u = femsolv(coords,nodes,bdryn,bdryv,@funcoefs_3D,@funrhs_3D);
    u_vals = u_node_vals(nodes, u);

    max_val = 0.0;
    for i=1:bdryn_len
        bd_crd = coords(bdryn(i),:);
        val = g(bd_crd(1), bd_crd(2), bd_crd(3)) - u_vals(bdryn(i));
        if val > max_val
            max_val = val;
        end
    end
    
    max_val
    
    N = size(u, 1);
    max_val = 0.0;
    for i=1:N
        val = abs(g(coords(i,1), coords(i,2), coords(i,3)) - u_vals(i));
        if val > max_val
            max_val = val;
        end
    end
    
    max_val
    
    u_h_eps(mm, 1) = hmax;
    u_h_eps(mm, 2) = max_val;
end

u_h_eps

loglog(u_h_eps(:,1), u_h_eps(:,2), "o"); hold on;
loglog([u_h_eps(MM,1), u_h_eps(1,1)], [u_h_eps(MM,2), u_h_eps(1,2)], "--*");

for i = 1:MM
    text(u_h_eps(i,1), u_h_eps(i,2) * 0.8, sprintf('h_%d', i), ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 8, ...
        'Rotation', 0);
end

% Labels and grid
xlabel('h'); ylabel('eps^h');
title('Log-Log Plot of h vs eps^h');

return;


function [coords,nodes,bdryn,hmax]=loadmesh(coordsf,nodesf,bdrynf,ifplot)
% Load a mesh from some files.
%
%                       Input arguments:
%
%  coordsf - a csv file containing the coordinates of the nodes.
%  nodesf - a csv file containing the nodes for each element.
%  bdrynf - a csv file containing the boundary nodes.
%  ifplot - a boolean determining whether or not plots will be created.
%
%                       Output arguments:
%
%  coords - a (nnodes,ndim) array containing the coordinates of the nodes, 
%       where ndim is the dimension and nnodes is the number of nodes.
%  nodes - the connectivity array; an (nels,nnels) array, containing the 
%       indices of the nodes in each element, where nels is the number
%       of elements and nnels is the number of nodes per element.

coords = readmatrix(coordsf);
[nnodes,~]=size(coords);
nodes = readmatrix(nodesf);
[nels,~]=size(nodes);
bdryn = readmatrix(bdrynf);
hmax=0;

MAX_FACES = 150;
faces = [1 2 3;
         1 2 4;
         1 3 4;
         2 3 4];

if (ifplot)
    figure(1);
    hold on;
    axis equal;
    view(3);
    xlabel('X'); ylabel('Y'); zlabel('Z');
end

for i=1:nels
    % order: z1 ... z10
    en=nodes(i,:);
    crns=en(1,1:4);
    vrts=coords(crns,:);

    if (ifplot && i < MAX_FACES)
        for f = 1:4
            patch('Vertices', vrts, ...
                  'Faces', faces(f,:), ...
                  'FaceColor', 'cyan', ...
                  'FaceAlpha', 0.1, ...
                  'EdgeColor', 'k');
        end
    end
    

    p1=vrts(1,:);
    p2=vrts(2,:);
    p3=vrts(3,:);
    p4=vrts(4,:);

    h=max([norm(p1-p2),norm(p1-p3),norm(p1-p4),norm(p2-p3),norm(p2-p4),norm(p3-p4)]);
    if (h > hmax)
        hmax=h;
    end
end

if (ifplot)
    coordsin=coords(setdiff(1:nnodes,bdryn),:);
    % scatter3(coordsin(:,1), ...
    %      coordsin(:,2), ...
    %      coordsin(:,3), ...
    %      20, 'b', 'filled');
    coordsbd=coords(bdryn,:);
    % scatter3(coordsbd(:,1), ...
    %      coordsbd(:,2), ...
    %      coordsbd(:,3), ...
    %      20, 'r', 'filled');
    hold off;
end    

return;
end


function u=femsolv(coords,nodes,bcdof,bcval,funcoefs,funrhs)
% Solve an equation using the finite element method.
%
%                       Input arguments:
%
%  coords - a (nnodes,ndim) array containing the coordinates of the nodes, 
%       where ndim is the dimension and nnodes is the number of nodes.
%  nodes - the connectivity array; an (nels,nnels) array, containing the 
%       indices of the nodes in each element, where nels is the number
%       of elements and nnels is the number of nodes per element.
%  bcdof - a vector containing the constrained degrees of freedom.
%  bcval - a vector containing the corresponding constrained values.
%  funcoefs - a function evaluating the coefficients.
%  funrhs - a function evaluating the right hand side.
%
%                       Output arguments:
%
%  u - the solution at the system degrees of freedom
%

[nels,~]=size(nodes);
[nnodes,~]=size(coords);
ndof=ones(nnodes,1);
sdof=sum(ndof);

kk=sparse(sdof,sdof);
ff=zeros(sdof,1);

nids=nodeindex(ndof);

mrkrs = [0, 0.25, 0.5, 0.75, 1];
mrkr_ind = 1;

for i=1:nels
    nd=nodes(i,:);
    xs=coords(nd,:);

    ids=eltindex(nd,nids,ndof);
    k=eltmat(xs,funcoefs);
    f=eltrhs(xs,funrhs);

    [kk,ff]=assemble(kk,ff,k,f,ids);

    if(i/nels >= mrkrs(mrkr_ind))
        fprintf('%.2f%% elements complete.\n', mrkrs(mrkr_ind) * 100);
        mrkr_ind = mrkr_ind + 1;
    end
end

[kk,ff]=applybcs(kk,ff,bcdof,bcval);

u=kk\ff;

return;
end


function k=eltmat(xs,funcoefs)
% Construct the element matrix for the element determined by the
% user-specified coordinates.
%
%                       Input arguments:
%
%  xs - an (nnels,ndim) vector of coordinates of the nodes associated with
%       the element, where ndim is the dimension and nnels is the number of
%       nodes per element.
%  funcoefs - a function handle accepting points as its argument; the
%       coefficients of the problem to be solved.
%
%                       Output arguments:
%
%  kk - a matrix of size (edof,edof), where edof is the number of degrees of
%       freedom per element; the element matrix.

[nnels,~]=size(xs);
k = zeros(nnels, nnels);

a = 0.58541020;
b = 0.13819660;
intg_pts = [
    b, b, b;
    a, b, b;
    b, a, b;
    b, b, a
];

for i=1:nnels
    for j=i:nnels
        k(i,j)=funcoefs(i, j, xs, intg_pts);
        k(j,i)=k(i,j);
    end
end

return;
end


function f=eltrhs(xs,funrhs)
% Construct the right hand side vector for the element determined by the 
% user-specified coordinates.
%
%                       Input arguments:
%
%  xs - an (nnels,ndim) vector of coordinates of the nodes associated with
%       the element, where ndim is the dimension and nnels is the number of
%       nodes per element.
%  funrhs - a function handle accepting points as its argument; the right
%       hand side of the problem to be solved.
%
%                       Output arguments:
%
%  ff - a vector of length edof, where edof is the number of degrees of freedom
%       per element; the element right hand side.

[nnels,~]=size(xs);
f = zeros(nnels, 1);
for i=1:nnels
    f(i)=funrhs(i, xs);
end

return;
end


function [kk,ff]=assemble(kk,ff,k,f,ids)
% Assemble the element matrices and right hand side vectors into the system
% matrix and right hand side vector.
%
%                       Input arguments:
%
%  kk - the system matrix, of size (sdof,sdof).
%  ff - the right hand side vector, of length sdof.
%  k - the element matrix, of size (edof,edof).
%  f - the element right hand side vector, of length edof.
%  ids - the degrees of freedom associated with the element.
%
%                       Output arguments:
%
%  kk,ff - the system matrix and right hand side vector, after assembly.
%

edof=length(ids);
for i=1:edof
    ii=ids(i);
    ff(ii)=ff(ii)+f(i);
    for j=1:edof
        jj=ids(j);
        kk(ii,jj)=kk(ii,jj)+k(i,j);
    end
end

return;
end


function [kk, ff] = applybcs(kk, ff, bcdof, bcval)
    % Apply Dirichlet boundary conditions correctly (non-destructive, vectorized).

    for k = 1:length(bcdof)
        id = bcdof(k);
        val = bcval(k);

        ff = ff - val * kk(:, id);   % Subtract DOF influence
        kk(id, :) = 0;               % Zero row
        kk(:, id) = 0;               % Zero column
        kk(id, id) = 1;              % Set diagonal
        ff(id) = val;                % Set RHS value
    end
end


function nids=nodeindex(ndof)
% Compute the number of previous system degrees of freedom for each node.
%
%                       Input arguments:
%
%  ndof - a vector of the number of degrees of freedom for each node.
%
%                       Output arguments:
%
%  nids - the number of previous system degrees of freedom for each node.

[nnodes,~]=size(ndof);
nids=zeros(nnodes,1);

ijk=0;
for i=1:nnodes
    nids(i)=ijk;
    ijk=ijk+ndof(i);
end

end


function ids=eltindex(nd,nids,ndof)
% Compute the system degrees of freedom associated with each element.
%
%                       Input arguments:
%
%  nd - node numbers whose degrees of freedom are to be determined.
%  nids - the number of previous system degrees of freedom for each node.
%  ndof - number of degrees of freedom per node.
%
%                       Output arguments:
%
%  ids - the index vector; an array of the indices of the system degrees
%       of freedom associated with the nodes in nd.
%
nnels=length(nd);
edof=sum(ndof(nd));
ids=zeros(edof,1);

ijk=0;

for i=1:nnels
    start=nids(nd(i));
    for j=1:ndof(nd(i))
        ijk=ijk+1;
        ids(ijk)=start+j;
    end
end

return;
end


function cf=funrhs_3D(~, ~)
    cf = 0;
end

function cf=funcoefs_3D(i, j, xs, pts)

function val = compute_integrand(xi1, xi2, xi3)
    % xs: 10x3 array of physical coordinates

    dphi = zeros(10,3);

    % Derivatives of phi1 = xi1*(2*xi1 - 1)
    dphi(1,:) = [4*xi1 - 1, 0, 0];
    
    % Derivatives of phi2 = xi2*(2*xi2 - 1)
    dphi(2,:) = [0, 4*xi2 - 1, 0];
    
    % phi3 = (1 - xi1 - xi2 - xi3)*(1 - 2*xi1 - 2*xi2 - 2*xi3)
    l1 = 4*(xi1 + xi2 + xi3) - 3;
    dphi(3,:) = [l1, l1, l1];
    
    % phi4 = xi3*(2*xi3 - 1)
    dphi(4,:) = [0, 0, 4*xi3 - 1];
    
    % phi5 = 4*xi1*xi2
    dphi(5,:) = [4*xi2, 4*xi1, 0];
    
    % phi6 = 4*xi2*(1 - xi1 - xi2 - xi3)
    dphi(6,1) = -4*xi2;
    dphi(6,2) = 4*(1 - xi1 - 2*xi2 - xi3);
    dphi(6,3) = -4*xi2;
    
    % phi7 = 4*xi1*(1 - xi1 - xi2 - xi3)
    dphi(7,1) = 4*(1 - 2*xi1 - xi2 - xi3);
    dphi(7,2) = -4*xi1;
    dphi(7,3) = -4*xi1;
    
    % phi8 = 4*xi1*xi3
    dphi(8,:) = [4*xi3, 0, 4*xi1];
    
    % phi9 = 4*xi2*xi3
    dphi(9,:) = [0, 4*xi3, 4*xi2];
    
    % phi10 = 4*xi3*(1 - xi1 - xi2 - xi3)
    dphi(10,1) = -4*xi3;
    dphi(10,2) = -4*xi3;
    dphi(10,3) = 4*(1 - xi1 - xi2 - 2*xi3);

    J = dphi' * xs;

    function grad = grad_phi(ind)
        % Gradient in physical coordinates
        grad = J \ dphi(ind,:)';
    end

    val = abs(det(J)) * dot(grad_phi(i), grad_phi(j));
end

    cf = 0;
    for k=1:4
        cf = cf + compute_integrand(pts(k,1), pts(k,2), pts(k,3));
    end
    cf = cf / 24;

    return;
end

function s1 = update(s);

s1 = s;

imax = s.I;
jmax = s.J;

dx = s.dx;
dy = s.dy;

% Left boundary condition is tau = 1
M = eye(jmax);  % Linear matrix will be imax x jmax
D = ones(jmax, 1);       % Vector of constants


for i=2:imax-1

    % A, B, C matrices are jmax x jmax at each value of i
    A = sparse(zeros(jmax,jmax)); 
    B = sparse(zeros(jmax,jmax));
    C = sparse(zeros(jmax,jmax));
    
    omx = omega(s.x(i)) / s.a;      % Streamwise grid clustering function
    ompx = omegap(s.x(i)) / s.a;    % Its derivative
   
    om1 =  omega(s.y(1)) / s.b;
    omjmax =  omega(s.y(jmax)) / s.b;
    
    for j=2:jmax-1 % Replace value at 1 and jmax
      
       omy = omega(s.y(j)) / s.b;   % Wall-normal grid clustering function
       ompy = omegap(s.y(j)) / s.b; % Its derivative
       
       Ac = omy^2.0 / dy^2.0;   % Constant
      
      % Backflow switch for the upwind scheme.
      if (s.u(i,j) >= 0)
          nu =1;
      else
          nu =0;
      end
      
     % c parameters, see Cassel et al. (1995)
     c(j) = -2.0 * Ac;
     cmin(j) = Ac - omy / (2.0 * dy) * (ompy - s.v(i,j)) ;
     cplus(j) = 2.0*Ac - cmin(j);
     
     d(j) =  0.0; 


     % ci parameters, see Cassel et al. (1995)
     ci(j) =  c(j)  - omx/dx * s.u(i,j) * ( nu - (1-nu) );
     ciplus(j) =    - omx * (1-nu) * s.u(i,j) / dx;
     cimin(j) =     + omx * nu * s.u(i,j) / dx;
     
     % Define matrices
     A(j,j) = cimin(j); 
     B(j,j) = ci(j);    B(j,j-1) = cmin(j); B(j,j+1) = cplus(j);
     C(j,j) = ciplus(j); 
  

    end
     
% end
    
    % Interaction law
    % Central for dtau / dx
    A(1,1) =        omx * ompx * dy / (4.0 * dx) / om1 - omx^2.0 * dy / (2.0 * dx^2.0) /om1;
    A(1,jmax) =     omx * ompx * dy / (4.0 * dx) /omjmax - omx^2.0 * dy / (2.0 * dx^2.0) /omjmax;
    
    B(1,1) =    omx^2.0 * (dy/dx^2.0)/om1 + om1 / dy;
    B(1,jmax) = omx^2.0 * (dy/dx^2.0) / omjmax;
    
    C(1,1) =    - omx * ompx * dy / (4.0 * dx) /om1 - omx^2.0 * dy / (2.0 * dx^2.0) / om1;
    C(1,jmax) = - omx * ompx * dy / (4.0 * dx) /omjmax - omx^2.0 * dy / (2.0 * dx^2.0) /omjmax;
    
    for j=2:jmax-1 
        omy = omega(s.y(j)) / s.b;
        
        A(1,j) =    omx * ompx * dy / (4.0 * dx) *2.0/ omy - omx^2.0 * dy / (2.0 * dx^2.0) * 2.0 / omy; 
        B(1,j) =    omx^2.0 * (dy/dx^2.0) * 2.0 / omy;
        if (j ==2) 
            B(1,j) = B(1,j) - om1 / dy; % d tau / dy first order discretisation
        end
        C(1,j) =    - omx * ompx * dy / (4.0 * dx) *2.0 / omy - omx^2.0 * dy / (2.0 * dx^2.0) * 2.0 / omy;
    end
    d(1) = - surfd2(s,s.xr(i));

     A(jmax,jmax)   = 0.0;
     B(jmax,jmax-1) = 0.0;
     B(jmax,jmax)   = 1.0;
     C(jmax,jmax)   = 0.0;
    
     d(jmax) = 1.0 ;
     
    % Join matrices
    di1 = (i-1)*jmax+1; di2 = i*jmax; % Row corresponding to i
    dj1a = di1 - 1*jmax; dj2a = di2 - 1*jmax; 
    dj1b = di1 - 0*jmax; dj2b = di2 - 0*jmax;
    dj1c = di1 + 1*jmax; dj2c = di2 + 1*jmax;
    
    D( di1:di2, 1) = transpose(d);
    
    % [A B C]
    M( di1:di2, dj1a:dj2a) = A;
    M( di1:di2, dj1b:dj2b) = B;
    M( di1:di2, dj1c:dj2c) = C;
   
    M = sparse(M);
    D = sparse(D);
    
end

% Right BC is tau = 1
tauI = ones(1,jmax);
D = [D; transpose([tauI])]; 
M((imax-1)*jmax+1 : imax*jmax, 1:end ) =  zeros(jmax, imax * jmax);
M((imax-1)*jmax+1 : imax*jmax, (imax-1)*jmax + 1 : imax*jmax) = eye(jmax);

% Solve the linear system
M = sparse(M);
D = sparse(D);

tau = M \ D; % Direct linear solver

s1.tau(1,:) = ones(1,jmax);
for i=2:imax
    s1.tau(i,:) = tau((i-1)*jmax + 1: i*jmax);
end

% Compute u, psi, v and p from tau.
s1.u    = find_u(s1); 
s1.psi  = find_psi(s1);
s1.v    = find_v(s1);
s1.p    = find_pressure(s1);
end

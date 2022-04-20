function G_transient = Sweby(G_00, w, k, dx, domain, G_0 , horizon)
% G_00: initial profile of OC
% w: sed rate [cm][yr]-1
% k: reaction rate constant [yr]-1
% dx: spatial resolution [cm]
% domain: spatial domain [cm]
% G_0: OC at sediment-water interface 
% horizon: depth of oxidant depletion

% Boudreau 1997 p.347

dt = dx/w * 0.8;     % Courant number requirement w/ (arbitrary) safety factor = 0.8
n_t = round(1/dt)+1; % number of grids in time
dt = 1/(n_t-1);      % adjust dt, dt*(n_t-1) = 1 for 1-year simulation

n_x = round(domain/dx); % number of grids in space

G_m = NaN (n_x,n_t);    % define matrix of OC

% initialize with proper initial/boundary conditions
G_m (1,:) = G_0;
G_m (:,1) = G_00;

% define matrices for numerical iteration
F   = NaN (n_x,n_t);      
v   = NaN (n_x,n_t);
r   = NaN (n_x,n_t);
phi = NaN (n_x,n_t);

t   = linspace(0,1,n_t); % time at each gridpoint 
           
for j = 1: n_t
    G_m(1,j) = G_0./exp(-k*t(j));
    F(:,j)   = w .* G_m(:,j);
    for i = 1: n_x-1
        if G_m(i+1,j)-G_m(i,j) ~= 0
            v(i,j) = 1 - dt/dx * ( (F(i+1,j)-F(i,j))/(G_m(i+1,j)-G_m(i,j)) );  
        else
            v(i,j) = 1;
        end
        
    end
    v(n_x,j) = 1;
    for i = 2: n_x-1
        r(i,j) = v(i-1,j)*(F(i,j)-F(i-1,j)) / ( v(i,j) * (F(i+1,j)-F(i,j)) );
    end
    r(1,j)  = 0;
    r(n_x,j)= 0;
    for i = 1: n_x
        if r(i,j) < 0
            phi(i,j) = 0;
        else
            if r(i,j) <= 1
                phi(i,j) = r(i,j);
            else
            phi(i,j) = 1;
            end
        end
    end
    
    for i = 2: n_x-1    
        G_m(i,j+1) = G_m(i,j)- dt/dx * (F(i,j)-F(i-1,j))+...
                            dt/(2*dx)* (phi(i,j)  * v(i,j)  * (F(i+1,j)-F(i,j))-...
                                        phi(i-1,j)* v(i-1,j)* (F(i,j)  -F(i-1,j)));
    end
    G_m(n_x,j+1) =  G_m(n_x,j)- dt/dx * (F(n_x,j)-F(n_x-1,j))+...
                            dt/(2*dx)* (phi(n_x,j)  * v(n_x,j)  * (0-F(n_x,j))-...
                                        phi(n_x-1,j)* v(n_x-1,j)* (F(n_x,j) -F(n_x-1,j)));   
end

% for x < horizon, k = k
% for x > horizon, k = 0
G_transient(1:horizon*(1/dx)+1) = G_m(1:horizon*(1/dx)+1,n_t) .* exp (-k*t(j));
G_transient(horizon*(1/dx)+2:n_x) = G_m(horizon*(1/dx)+2:n_x,n_t);   

end

                                    
                                    
                                    
                                    
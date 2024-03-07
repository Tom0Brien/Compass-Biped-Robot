classdef PortHamiltonianSystem < handle
    properties
        n_x; % Dimension of the full state vector
        n_r; % Dimension of the reduced momentum vector 
        M_0; % Mass matrix function handle
        D_0; % Damping matrix function handle
        G_0; % Input matrix function handle
        V; % Potential energy function handle
        G_c; % Constraint matrix function handle
        Q; % Transposed left-annihilator of the constraint matrix
        M; % Reduced mass matrix
        D; % Reduced damping matrix
        G; % Reduced input matrix
        H; %  Hamiltonian
        dHdq; % Gradient of H w.r.t q
        dHdp; % Gradient of H w.r.t p
    end
    
    methods
        function obj = PortHamiltonianSystem(mass_matrix, damping_matrix, input_matrix, potential_energy, n_x)
            obj.n_x = n_x;
            obj.M_0 = mass_matrix;
            obj.D_0 = damping_matrix;
            obj.G_0 = input_matrix;
            obj.V = potential_energy;
            obj.G_c = [];
            
            syms p [n_x/2 1];
            syms q [n_x/2 1];
            
            H_0 = 0.5 * (p.' / obj.M_0(q) * p) + potential_energy(q);
            obj.H = H_0;
            obj.dHdq = matlabFunction(gradient(H_0, q), 'Vars', {q, p});
            obj.dHdp = matlabFunction(gradient(H_0, p), 'Vars', {q, p});
        end
        
        function set_constraints(obj, G_c)
            syms q [obj.n_x/2 1];

            obj.G_c = G_c;
            obj.Q = matlabFunction(null(G_c(q).'), 'Vars', {q});
            obj.M = @(q) obj.Q(q).' * obj.M_0(q) * obj.Q(q);
            obj.D = @(q) obj.Q(q).' * obj.D_0(q) * obj.Q(q);
            obj.G = @(q) obj.Q(q).' * obj.G_0(q);
            obj.n_r = size(obj.M(q), 1);
            syms p [obj.n_r 1];
            
            obj.H = 0.5 * p.' / obj.M(q) * p + obj.V(q);
            obj.dHdq = matlabFunction(gradient(obj.H, q), 'Vars', {q, p});
            obj.dHdp = matlabFunction(gradient(obj.H, p), 'Vars', {q, p});
        end
        
        function dx = forward_dynamics(obj, state, input)
            if isempty(obj.G_c)
                q = state(1:length(state)/2);
                p = state(length(state)/2+1:end);
                dx = double([obj.dHdp(q, p); -obj.dHdq(q, p) - obj.D_0(q) * obj.dHdp(q, p) + obj.G_0(q) * input]);
            else
                q = state(1:length(state)/2);
                p = state(length(state)/2+1:length(state)/2+obj.n_r);
                dx = double([obj.Q(q) * obj.dHdp(q, p); ...
                             -obj.Q(q).' * obj.dHdq(q, p) - obj.D(q) * obj.dHdp(q, p) + obj.G(q) * input; ...
                             zeros(obj.n_x/2 - obj.n_r, 1)]);
            end
        end
        
        function post_impact_state = impact_mapping(obj, pre_impact_state, new_constraints)
            q_minus = pre_impact_state(1:length(pre_impact_state)/2);
            p_minus = pre_impact_state(length(pre_impact_state)/2+1:length(pre_impact_state)/2+obj.n_r);
            Q_minus = obj.Q(q_minus);
            M_minus = obj.M(q_minus);
            
            obj.set_constraints(new_constraints);

            Q_plus = obj.Q(q_minus);
            p_plus = (Q_plus.' * obj.M_0(q_minus) * Q_minus) / M_minus * p_minus;
            post_impact_state = [q_minus; p_plus; zeros(obj.n_x - length(q_minus) - length(p_plus), 1)];
        end
    end
end

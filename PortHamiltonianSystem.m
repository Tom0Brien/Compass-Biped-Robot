classdef PortHamiltonianSystem < handle
    properties
        n_x;    % Dimension of the full state vector
        n_r;    % Dimension of the reduced momentum vector 
        M_0;    % Mass matrix function handle
        D_0;    % Damping matrix function handle
        G_0;    % Input matrix function handle
        V;      % Potential energy function handle
        J_s;    % Constraint matrix function handle
        J_null; % Transposed left-annihilator of the constraint matrix
        M;      % Reduced mass matrix
        D;      % Reduced damping matrix
        G;      % Reduced input matrix
        H;      % Hamiltonian
        dHdq;   % Gradient of H w.r.t q
        dHdp;   % Gradient of H w.r.t p
    end
    
    methods
        function obj = PortHamiltonianSystem(mass_matrix, damping_matrix, input_matrix, potential_energy, n_x)
            obj.n_x = n_x;
            obj.M_0 = mass_matrix;
            obj.D_0 = damping_matrix;
            obj.G_0 = input_matrix;
            obj.V = potential_energy;
            obj.J_s = [];
            
            syms p [n_x/2 1];
            syms q [n_x/2 1];
            
            H_0 = 0.5 * (p.' / obj.M_0(q) * p) + potential_energy(q);
            obj.H = H_0;
            obj.dHdq = matlabFunction(gradient(H_0, q), 'Vars', {q, p});
            obj.dHdp = matlabFunction(gradient(H_0, p), 'Vars', {q, p});
        end
        
        function set_constraints(obj, J_s)
            syms q [obj.n_x/2 1];

            obj.J_s = J_s;
            obj.J_null = matlabFunction(null(J_s(q).'), 'Vars', {q});
            obj.M = @(q) obj.J_null(q).' * obj.M_0(q) * obj.J_null(q);
            obj.D = @(q) obj.J_null(q).' * obj.D_0(q) * obj.J_null(q);
            obj.G = @(q) obj.J_null(q).' * obj.G_0(q);
            obj.n_r = size(obj.M(q), 1);
            syms p [obj.n_r 1];
            
            obj.H = 0.5 * p.' / obj.M(q) * p + obj.V(q);
            obj.dHdq = matlabFunction(gradient(obj.H, q), 'Vars', {q, p});
            obj.dHdp = matlabFunction(gradient(obj.H, p), 'Vars', {q, p});
        end
        
        function dx = forward_dynamics(obj, state, input)
            if isempty(obj.J_s)
                q = state(1:length(state)/2);
                p = state(length(state)/2+1:end);
                dx = double([obj.dHdp(q, p); -obj.dHdq(q, p) - obj.D_0(q) * obj.dHdp(q, p) + obj.G_0(q) * input]);
            else
                q = state(1:length(state)/2);
                p = state(length(state)/2+1:length(state)/2+obj.n_r);
                dx = double([obj.J_null(q) * obj.dHdp(q, p); ...
                             -obj.J_null(q).' * obj.dHdq(q, p) - obj.D(q) * obj.dHdp(q, p) + obj.G(q) * input; ...
                             zeros(obj.n_x/2 - obj.n_r, 1)]);
            end
        end
        
        function post_impact_state = impact_mapping(obj, pre_impact_state, new_constraints)
            q_minus = pre_impact_state(1:length(pre_impact_state)/2);
            p_minus = pre_impact_state(length(pre_impact_state)/2+1:length(pre_impact_state)/2+obj.n_r);
            J_null_minus = obj.J_null(q_minus);
            M_minus = obj.M(q_minus);
            
            obj.set_constraints(new_constraints);

            J_null_plus = obj.J_null(q_minus);
            p_plus = (J_null_plus.' * obj.M_0(q_minus) * J_null_minus) / M_minus * p_minus;
            post_impact_state = [q_minus; p_plus; zeros(obj.n_x - length(q_minus) - length(p_plus), 1)];
        end
    end
end

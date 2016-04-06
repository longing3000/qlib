classdef GFNMDecohSuperOperator < model.phy.QuantumOperator.SpinOperator.DecoherenceSuperOperator.AbstractDecoherenceOperator
    %GFNMDecohSuperOperator generate a superoperator to describe the decoherence of the multi-spin system.
    % The master equation is generated from a gradient field noise, i.e.,
    % all the spin experience the same noise field with different
    % amplitude: h_j(t)= B^z(x_j)*f(t). Here, B^z(x_j) is the amplitude of
    % the noise field along the z direction, f(t) is a stochastic function, whoes spectrum is normalized:
    % S(omega)=\int_{-\infty}^{\infty} ds <f(t)f(t-s)> exp(-i*omega*s);
    % If we assume <f(t)f(t-s)> = exp(-|s|/ tau_z)/(2*pi), the corresponding spectrum is
    % S(omega) = (1/pi)*tau_z/(omega^2*tau_z^2+1); % tau_z is the correlation time of the noise field.
    % \int_{-\infty}^{\infty} domega S(omega)=1; The integral of this spectrum is 1.
 
    % In this case, the nose fields exerted on different spins are
    % correlated, i.e., <h_i(t)h_j(t')> \neq = 0;
    % the total decoherence operator in the Liouville space is the summation of the decoherence operator 
    % for each single spin and all pairs. But if the spin is not spin half, the corresponding decoherence operator is 
    % replaced by the zero matrix. 

    %   The decoherence operator of the single spin is of the Lindblad form:
    %  L rho(t) = \sum_{j} Gamma_{j}*[sigma_{j}^{z}*rho*sigma_{j}^{z} - rho(t)] % the diagonal part
    %              \sum_{i<j} Gamma_{ij}*[sigma_{i}^{z}*rho(t)*sigma_{j}^{z} + sigma_{j}^{z}*rho(t)*sigma_{i}^{z}
    %              -sigma_{i}^{z}*sigma_{j}^{z}*rho(t) - rho(t)*sigma_{i}^{z}*sigma_{j}^{z}] 
    % where sigma is the Pauli operator of the spin, 
    % Gamma_{j}= [gamm_j*B^z(x_j)]^2 * S(0)/4 = |B^z(x_j)|^2 * tau_z/pi/4;
    % Gamma_{ij}= [gamm_i*B^z(x_i)*gamm_j*B^z(x_j)] * S(0)/4 = [gamm_i*B^z(x_i)*gamm_j*B^z(x_j)] * tau_z/pi/4;
    % gamma_j  is the gyromagnetic ratio of the jth spin
    
    properties
        correlation_time
        
    end
    
    methods
        function obj=GFNMDecohSuperOperator(spin_collection,correlation_time)
            obj@model.phy.QuantumOperator.SpinOperator.DecoherenceSuperOperator.AbstractDecoherenceOperator(spin_collection);
            obj.correlation_time=correlation_time;%The correlation time of the noise field;
            obj.name='decoherence_super_operator';            
        end
        
        function generate_matrix(obj)
            obj.matrix=obj.calculate_matrix();
            obj.hasMatrix=1;
        end
        
        function L_matrix = calculate_matrix(obj)
            dim_tot=obj.spin_collection.getDim;
            nspin=obj.spin_collection.getLength;
            if norm(dim_tot-2^nspin)>1e-10
               disp('This class can only generate the decoherence operator for spin-half system.');
               L_matrix=0;
               return;
            end
            
            L_single=obj.calculate_single_spin_liouvillian;
            L_pair=obj.calculate_spin_pair_liouvillian;            
            L_matrix=L_single+L_pair;

        end
        
        function L_matrix=calculate_single_spin_liouvillian(obj)
            nspin=obj.spin_collection.getLength;            
            L_matrix=0;
            
            for kk=1:nspin
                spin=obj.spin_collection.spin_list{kk};
                if spin.dim>2
                    L_kk=0;
                    L_matrix=L_matrix+L_kk;
                elseif spin.dim==2
                    L_kk=obj.gen_superoperator_single(spin,kk);

                    L_matrix=sparse(L_matrix+L_kk);

                end
            end            
        end
        
        function L_matrix=calculate_spin_pair_liouvillian(obj)
            nspin=obj.spin_collection.getLength;
            L_matrix=0;
            if nspin<2
               % do nothing in this case.
            else
                iter_class=@model.phy.SpinCollection.Iterator.PairSpinIterator;
                iter=iter_class(obj.spin_collection);
                nPair=iter.getLength;
                for kk=1:nPair
                    spins=iter.currentItem();
                    spin_idx=iter.currentIndex();
                    L_pair=obj.gen_superoperator_pair(spins,spin_idx);
                    L_matrix=L_matrix+L_pair;
                    iter.nextItem;
                end
            end
        end
        
        function L_single=gen_superoperator_single(obj,spin,kk)
            % gamma_{j}= [B^z(x_j)]^2 * S(0) = |B^z(x_j)|^2 * tau_z/pi;
            gradient_field=spin.local_field;
            gammaj=spin.gamma;%the gyromagnetic ratio
            Gammaj=gradient_field*gradient_field'*(gammaj^2)*obj.correlation_time/pi;% the parallel decay rate
            disp(['The pure dephasing rate of the ' num2str(kk) 'th spin is ' num2str(Gammaj)]);

            dim=spin.dim;
            if dim ~=2
               error('This class can only generate the decoherence operator for spin-half system.'); 
            end
                                  
            dim_tot=obj.spin_collection.getDim;
            spin_collection=obj.spin_collection;
            % generate ladder operators
            import model.phy.QuantumOperator.SpinOperator.Observable
            import model.phy.QuantumOperator.MatrixStrategy.FromKronProd

            strategy=FromKronProd();
%             Sm=Observable(spin_collection,strategy, 'sigma-', ['1.0 * mat([0,0;1,0])_' num2str(kk)]);
%             sigma_m=full(Sm.getMatrix);
%             Sp=Observable(spin_collection,strategy, 'sigam+', ['1.0 * mat([0,1;0,0])_' num2str(kk)]);
%             sigma_p=full(Sp.getMatrix);
            Sz=Observable(spin_collection,strategy, 'sigmaz', ['1.0 * mat([1,0;0,-1])_'  num2str(kk)]);
            sigma_z=full(Sz.getMatrix);


                % generate supperoperator for the current spin
                %Here, we just consider the pure dephasing case
%                 if Gamma_v>0            
%                     L_v=0.5*Gamma_v*(2*kron(sigma_p',sigma_m)-kron(speye(dim_tot),sigma_p*sigma_m)-kron((sigma_p*sigma_m)',speye(dim_tot)))...
%                         +0.5*Gamma_v*(2*kron(sigma_m',sigma_p)-kron(speye(dim_tot),sigma_m*sigma_p)-kron((sigma_m*sigma_p)',speye(dim_tot)));
%                     %In the kron product, we have used the Hermitian conjugation to replace the transpose operation. because all the matrix is real
% 
%                 else
%                    L_v=0; 
%                 end
                
              L_p=Gammaj*(kron(sigma_z',sigma_z)-speye(dim_tot^2));
                
              L_single=sparse(L_p);

        end
        function L_single=gen_superoperator_pair(obj,spins,index)  

            % gamma_{ij}=[B^z(x_i)*B^z(x_j)] * tau_z/pi;
            g_field1=spins{1}.local_field;
            g_field2=spins{2}.local_field;
            gamma1=spins{1}.gamma;%the gyromagnetic ratio
            gamma2=spins{2}.gamma;%the gyromagnetic ratio
            Gammaij=g_field1*g_field2'*gamma1*gamma2*obj.correlation_time/pi/4;% the parallel decay rate
            disp(['The pure dephasing rate of the pair [' num2str(index(1) ) ' ' ... 
                num2str(index(2) ) '] is ' num2str( Gammaij )]);
                                  
            dim_tot=obj.spin_collection.getDim;
            spin_collection=obj.spin_collection;
            % generate ladder operators
            import model.phy.QuantumOperator.SpinOperator.Observable
            import model.phy.QuantumOperator.MatrixStrategy.FromKronProd

            strategy=FromKronProd();
            S1z=Observable(spin_collection,strategy, 'sigmaz', ['1.0 * mat([1,0;0,-1])_'  num2str(index(1))]);
            sigma_1z=full(S1z.getMatrix);
            S2z=Observable(spin_collection,strategy, 'sigmaz', ['1.0 * mat([1,0;0,-1])_'  num2str(index(2))]);
            sigma_2z=full(S2z.getMatrix);
            S12z=Observable(spin_collection,strategy, 'sigmaz', ['1.0 * mat([1,0;0,-1])_'  num2str( index(1) )...
                ' * mat([1,0;0,-1])_'  num2str(index(2))]);
            sigma_12z=full(S12z.getMatrix);

                
              L_p=Gammaij*(kron(sigma_1z',sigma_2z)+kron(sigma_2z',sigma_1z)...
                  -kron(sigma_12z',speye(dim_tot))-kron(speye(dim_tot),sigma_12z) );
                
              L_single=sparse(L_p);

        end
    end
    
end


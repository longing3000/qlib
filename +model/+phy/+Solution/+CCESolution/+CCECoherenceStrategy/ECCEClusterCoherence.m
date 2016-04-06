classdef ECCEClusterCoherence < model.phy.Solution.CCESolution.CCECoherenceStrategy.AbstractClusterCoherence
    %ECCECLUSTERCOHERENCE 
                %calculate the coherence of single cluster for ensemble CCE 
    properties
      coherence_tilde
    end
    
    methods
        function obj=ECCEClusterCoherence(cluster_spin_index,cluster_parameters)
            obj@model.phy.Solution.CCESolution.CCECoherenceStrategy.AbstractClusterCoherence();
            if nargin >0
               obj.generate(cluster_spin_index,cluster_parameters); 
            end                    
        end
        
        function coh=calculate_cluster_coherence(obj,evolution_para,varargin)
            obj.npulse=evolution_para.npulse;
            center_spin_states=evolution_para.center_spin_states;
            is_secular=evolution_para.is_secular;
            obj.timelist=evolution_para.timelist;
            
            %generate the spin_collection for this cluster including the central spin
            obj.spin_collection= model.phy.SpinCollection.SpinCollection();
            obj.spin_collection.spin_source=model.phy.SpinCollection.Strategy.FromSpinList([{obj.center_spin},obj.cluster_bath_spin]);
            obj.spin_collection.generate();
             
            hamiCell=obj.gen_reduced_hamiltonian(center_spin_states,is_secular);
            [h_list,hami_prefactor]=obj.gen_hami_list(hamiCell);
            
            [bath_cluster_sc,denseMat,initial_state_type]=obj.set_initial_state;
            
            %Observable
            obs=model.phy.QuantumOperator.SpinOperator.Observable(bath_cluster_sc,'IdentityMatrix');
            obs.setMatrix(1);
            
            coh=obj.calculate_coherence_hilbert(h_list,hami_prefactor,obs,denseMat,initial_state_type);
        end
        function [bath_cluster_sc,denseMat,initial_state_type]=set_initial_state(obj)
            %generate a SpinCollection of bath spins in this cluster
            bath_cluster_sc= model.phy.SpinCollection.SpinCollection();
            bath_cluster_sc.spin_source=model.phy.SpinCollection.Strategy.FromSpinList(obj.cluster_bath_spin);
            bath_cluster_sc.generate();
                        
            % DensityMatrix
            denseMat=model.phy.QuantumOperator.SpinOperator.DensityMatrix(bath_cluster_sc,'IdentityMatrix');
            dim=denseMat.dim;            
            denseMat.setMatrix(eye(dim)/dim);
            
            initial_state_type='MixedState';
        end              
    end
    
end

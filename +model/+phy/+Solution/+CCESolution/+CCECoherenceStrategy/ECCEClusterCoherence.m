classdef ECCEClusterCoherence < model.phy.Solution.CCESolution.CCECoherenceStrategy.AbstractClusterCoherence
    %ECCECLUSTERCOHERENCE 
                %calculate the coherence of single cluster for ensemble CCE 
    properties
      coherence_tilde
    end
    
    methods
        function obj=ECCEClusterCoherence(cluster)
            obj@model.phy.Solution.CCESolution.CCECoherenceStrategy.AbstractClusterCoherence();
            if nargin >0
               obj.generate(cluster); 
            end                    
        end
        
        function coh=calculate_cluster_coherence(obj,center_spin_states,timelist,varargin)
            p = inputParser;
            addRequired(p,'center_spin_states');
            addRequired(p,'timelist');
            addOptional(p,'npulse',0,@isnumeric);
            addOptional(p,'is_secular',0,@isnumeric); 

            parse(p,center_spin_states,timelist,varargin{:});

            obj.npulse=p.Results.npulse;
            center_spin_states=p.Results.center_spin_states;
            is_secular=p.Results.is_secular;
            timelist=p.Results.timelist;

            hamiCell=obj.gen_reduced_hamiltonian(center_spin_states,is_secular);
            [~]=obj.gen_hami_list(hamiCell);
            
            bath_spins=obj.spin_collection.spin_list(2:end);
            bath_cluster= model.phy.SpinCollection.SpinCollection();
            bath_cluster.spin_source=model.phy.SpinCollection.Strategy.FromSpinList(bath_spins);
            bath_cluster.generate();
            
            % DensityMatrix
            denseMat=model.phy.QuantumOperator.SpinOperator.DensityMatrix(bath_cluster,'IdentityMatrix');
            dim=denseMat.dim;
            denseMat.setMatrix(eye(dim)/dim);
            

            if dim<2000
                coh=obj.calculate_coherence_hilbert(bath_cluster,denseMat,'MixedState',timelist);
            else
                coh=obj.calculate_coherence_liouville(bath_cluster,denseMat,'MixedState',timelist);
            end
            obj.coherence=coh;
        end
                       
    end
    
end

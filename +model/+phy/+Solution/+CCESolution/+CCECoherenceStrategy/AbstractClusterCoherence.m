classdef AbstractClusterCoherence < handle
    %ABSTRACTCLUSTERCOHERENCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties        
      center_spin% a Spin 
      spin_collection% a SpinCollection: the first spin is the central spin, the rest are bath spins of this cluster
      cluster_bath_spin% a cell including bath spins of this cluster
      bath_spin_collection% a SpinCollection including all bath spins for CCE
      cluster_spin_index% the index of the bath spins of this cluster in bath_spin_collection, e.g. [1,2,3] 
      
      coherence
      hamiltonian_cell
      hami_list      
      preoperator_factor
      npulse
      timelist
    end
    
    methods
        function obj=AbstractClusterCoherence(cluster_spin_index,cluster_parameters)   
            % The input "bath_spins" is a SpinCollection, including the  all bath spins. 
            % The input "cspin" is a Spin 
            % The input "cluster_spin_index" is the index of  bath spin of this cluster in bath_spins            
            if nargin>0
                obj.generate(cluster_spin_index,cluster_parameters);
            end
        end 
        function generate(obj,clst_indx,cluster_parameters)
            obj.center_spin=cluster_parameters.center_spin;
            obj.bath_spin_collection=cluster_parameters.bath_spin_collection;
            obj.cluster_spin_index=clst_indx;
            spin_cell=cell(1, length(clst_indx) );
            for kk=1:length(clst_indx)
                idx=clst_indx(kk);
                spin_cell{kk}=obj.bath_spin_collection.spin_list{idx};
            end
            obj.cluster_bath_spin=spin_cell;
        end
        %  generate reduced hamiltonian for the given central spin states
        function reduced_hami = gen_reduced_hamiltonian(obj,center_spin_state,is_secular)
            cluster=obj.spin_collection;
            hami_cluster=model.phy.QuantumOperator.SpinOperator.Hamiltonian(cluster);
            zee_interaction=model.phy.SpinInteraction.ZeemanInteraction(cluster);
            dip_interaction=model.phy.SpinInteraction.DipolarInteraction(cluster);
            hami_cluster.addInteraction(zee_interaction);
            hami_cluster.addInteraction(dip_interaction);
            ts=model.phy.QuantumOperator.SpinOperator.TransformOperator(cluster,{'1.0 * eigenVectors()_1'});
            hami_cluster.transform(ts);

            hami1=hami_cluster.project_operator(1, center_spin_state(1));
            hami2=hami_cluster.project_operator(1, center_spin_state(2));
            hami1.remove_identity();
            hami2.remove_identity();

            if is_secular && hami1.spin_collection.getLength>1
                approx1=model.phy.SpinApproximation.SpinSecularApproximation(hami1.spin_collection);
                approx2=model.phy.SpinApproximation.SpinSecularApproximation(hami2.spin_collection);
                hami1.apply_approximation(approx1);
                hami2.apply_approximation(approx2);
            end

            reduced_hami=cell(1,2);
            reduced_hami{1,1}=hami1;
            reduced_hami{1,2}=hami2;
        end
        
        % Generate hamiltonian list for the evolution of ensemble CCE with a given pulse number
        function  [hami_list,hami_prefactor]= gen_hami_list(obj,hamiCell)
            time_seq=obj.time_ratio_seq();
            len_time_seq=length(time_seq);
            hami_list=cell(1,len_time_seq);
            p_list=(1:len_time_seq)+obj.npulse;
            parity_list=rem(p_list,2);
            for m=1:len_time_seq %initial the core matrice
                 parity=parity_list(m); 
                 switch parity
                        case 0
                            hami=hamiCell{1};                       
                        case 1
                            hami=hamiCell{2};
                        otherwise
                            error('wrong parity of the index of the hamiltonian sequence.');
                 end
                 hami_list{m}=hami;
            end
            hami_prefactor=1i*time_seq;
        end

        function time_seq = time_ratio_seq(obj)
            if obj.npulse==0
                time_seq=[-1,1];
            elseif obj.npulse>0
                nsegment=obj.npulse+1;
                step=1/obj.npulse/2;
                seq=zeros(1,nsegment);
                for n=1:nsegment
                    if n==1
                        seq(1,n)=step;
                    elseif n==nsegment
                        seq(1,n)=step;
                    else
                        seq(1,n)=2*step;
                    end
                end
                time_seq=[-1*seq,seq];
            end     
        end
            
        function coh=calculate_coherence_hilbert(obj,h_list,hami_prefactor,obs,state,statetype)
            time_list=obj.timelist;
            % Evolution
            switch statetype
                case 'MixedState' 
                    d_mat_evolution=model.phy.Dynamics.EvolutionKernel.DensityMatrixEvolution(h_list,statetype,hami_prefactor);
                    dynamics=model.phy.Dynamics.QuantumDynamics(d_mat_evolution);
                    dynamics.set_initial_state(state,'Hilbert');
                    dynamics.set_time_sequence(time_list);
                    dynamics.addObervable({obs});
                    dynamics.calculate_mean_values();
                case 'PureState'
                    wavefun_evolution=model.phy.Dynamics.EvolutionKernel.WavefunctionEvolution(h_list,statetype,hami_prefactor);
                    dynamics=model.phy.Dynamics.QuantumDynamics(wavefun_evolution);
                    dynamics.set_initial_state(state);
                    dynamics.set_time_sequence(time_list);
                    dynamics.addObervable({obs});
                    dynamics.calculate_mean_values(state');
                otherwise
                    error('No corresponding EvolutionKernel for  such kind of quantum state in Dynamics.');                    
            end
                    
            coh=dynamics.observable_values;
            obj.coherence=coh;
        end

        function coh=calculate_coherence_liouville(obj,liouList,prefactor,obs,state,statetype)
            time_list=obj.timelist;
            d_mat_evolution=model.phy.Dynamics.EvolutionKernel.MatrixVectorEvolution(liouList,statetype,prefactor);
            dynamics=model.phy.Dynamics.QuantumDynamics(d_mat_evolution);
            dynamics.set_initial_state(state,'Liouville');
            dynamics.set_time_sequence(time_list);
            dynamics.addObervable({obs});
            dynamics.calculate_mean_values();
            coh=dynamics.observable_values;
            obj.coherence=coh;
        end
    
        function [coh,coh_tilde]=calculater_cluster_coherence_tilde(obj,center_spin_states,timelist,varargin)
            %   Calculate the tilde coherence of a given cluster
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
            
            coh=obj.calculate_cluster_coherence(center_spin_states,timelist,'npulse',obj.npulse,'is_secular',is_secular);
            subclusters=obj.generate_subclusters();
            len_subclus=length(subclusters);
            if  len_subclus==0
                coh_tilde=coh;
            elseif len_subclus>0 
                 coh_tilde=coh;
                 for m=1:len_subclus
                    cluster=subclusters{m};
                    coh_tilde_sub=cluster.calculater_cluster_coherence_tilde(center_spin_states,timelist,'npulse',obj.npulse,'is_secular',is_secular);
                    coh_tilde=coh_tilde./coh_tilde_sub;
                end
            end

        end
        
        function subclusters=generate_subclusters(obj)
            cluster_sc=obj.spin_collection;
            nspin=cluster_sc.getLength;
            
            %generate index for all possible subclusters
            if nspin==2
                index_list=[];
            elseif nspin==3
                index_list=mat2cell([1,2;1,3],[1,1],2);
            elseif nspin>3               
                idx_base=2:nspin;
                index_list=mat2cell([ones(nspin-1,1),idx_base'],ones(1,nspin-1),2);
                for kk=2:nspin-2
                    mat=nchoosek(idx_base,kk);
                    len=size(mat,1);
                    mat=[ones(len,1),mat];
                    idx2add=mat2cell(mat,ones(1,len),kk+1);
                    index_list=[index_list;idx2add];
                end
            else
                error('This cluster is not a well defined cluster.');
            end
            
            %generate all possible subclusters
            if isempty(index_list)
                subclusters=[];
            else
                iter=model.phy.SpinCollection.Iterator.SpinIterator(cluster_sc,index_list);
                nsubcluster=iter.getLength();
                subclusters=cell(1,nsubcluster);
                
                classname=class(obj);
                for n=1:nsubcluster
                   spins=iter.getItem(n);
                   subcluster_sc=model.phy.SpinCollection.SpinCollection();% spin collection of a cluster
                   subcluster_sc.spin_source=model.phy.SpinCollection.Strategy.FromSpinList(spins);
                   subcluster_sc.generate();

                   subcluster=eval([classname, '(subcluster_sc)']);                   
                   subclusters{1,n}=subcluster;
                end
            end
            
        end
    end
    
        
    
    methods (Abstract)
        coh=calculate_cluster_coherence(obj,para);
    end
    
end


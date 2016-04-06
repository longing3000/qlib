classdef SingleSampleCCESolution < model.phy.Solution.CCESolution.AbstractCCESolution
    %ENSEMBLECCESOLUTION Summary of this class goes here
    %   EnsembleCCESolution needs the following input paramters:
    %   1. parameters.SpinCollectionStrategy
    %   2. parameters.InputFile
    %   3. parameters.SetBathSpins
    %   4. parameters.SetCentralSpin 
    %   5. parameters.MagneticField
    %   6.parameters.CutOff
    %   7.parameters.MaxOrder    
    %   9. parameters.IsSecularApproximation
    %   6. parameters.NPulse
    %   10. parameters.NTime
    %   11. parameters.TMax
    %   12. parameters.TimeList

   
    properties
    end
    
    methods
        function obj=SingleSampleCCESolution(xml_file)
            obj@model.phy.Solution.CCESolution.AbstractCCESolution(xml_file);
        end
        
         function perform(obj)
           Condition=model.phy.LabCondition.getCondition;              
           Condition.setValue('magnetic_field',obj.parameters.MagneticField);  
             
           %  Generate Spin Collection FromSpinList and generate clusters 
           cluster_iterator=obj.generate_cluster_iterator();
           [evolution_para,cluster_para]=obj.pre_calculation(cluster_iterator);
           obj.calculate_total_coherence(evolution_para,cluster_para,cluster_iterator); 
           disp('Calculation of this solution finishes.');
        end
 %%       
        function [evolution_parameter,cluster_parameter]=pre_calculation(obj,cluster_iterator)
           import model.phy.PhysicalObject.NV

           evolution_parameter.center_spin_states=obj.parameters.SetCentralSpin.CentralSpinStates;
           evolution_parameter.timelist=obj.parameters.TimeList;
           evolution_parameter.npulse=obj.parameters.NPulse;
           evolution_parameter.is_secular=obj.parameters.IsSecularApproximation;
           evolution_parameter.MagneticField=obj.parameters.MagneticField; 
           evolution_parameter.strategy_name=obj.parameters.CCEStrategy;           
                      
           center_spin_name=obj.parameters.SetCentralSpin.name;
           para_central_spin=obj.parameters.SetCentralSpin; 
           center_spin=eval(strcat(center_spin_name,'(','para_central_spin',')'));
           obj.keyVariables('center_spin')=center_spin;
           
           cluster_parameter.center_spin=center_spin.espin;
           cluster_parameter.bath_spin_collection=cluster_iterator.spin_collection;
           cluster_parameter.bath_spin_state=obj.generate_bath_spin_state(cluster_iterator);
        end
        function bs_state=generate_bath_spin_state(obj,cluster_iterator)
            seed=obj.parameters.seed;
            nspin=cluster_iterator.spin_collection.getLength;
            dim_list=cluster_iterator.spin_collection.getDimList;
            rng(seed);rand_numbers=randi([1,100],1,nspin);
            
            bs_state=zeros(1,nspin);
            for kk=1:nspin
               dim=dim_list(kk); 
               bs_state(1,kk)=mod(rand_numbers(kk),dim)+1; 
            end
        end
        
%%        
        function calculate_total_coherence(obj, evolu_para,clst_para,cluster_iter) 
%          evolution parameters  
%            evolu_para.npulse:: the pulse number
%            evolu_para.is_secular:: determine whether the secular approximation is taken or not
%            evolu_para.MagneticField:: the magnetic field vector             
%            evolu_para.center_spin_states:: the states of the center spin involved in this calculatioin
%            evolu_para.strategy_name:: CCE strategy name
%            evolu_para.timelist:: the evolution time list

%         cluster parameters
%           clst_para.bath_spin_collection:: a SpinCollection, including the  all bath spins
%           clst_para.center_spin:: a Spin
%           clst_para.bath_spin_state:: a list gives out the specific configuration of the bath
%                                       spin state for SingleSampleCCE,
%                                       e.g. [1,1,1,1,1,1,1,...,1] for the ground state.                                      


           ncluster=cluster_iter.getLength;
           ntime=length(evolu_para.timelist);           
           cluster_index_list=cluster_iter.index_list;
           MagneticField=evolu_para.MagneticField;
           
           CoherenceMatrix=zeros(ncluster,ntime);
           disp('calculate the cluster-coherence matrix ...');
           tic
           
          %In order to eliminate the parfor warning, I have to arrange the
          %swith...case... in this form. Beside, I try to add a field "clst_index" in the structure data "clst_para" in every loop first. 
          % But this action is forbidden in parfor circulation. So I have to change the way to construct 
          % the AbstractClusterCoherence class. This is pretty ulgy, but I have to do this. 

          parpool();
          parfor n=1:ncluster 
              disp(['calculating the ' num2str(n) 'th cluster coherence....']);
              Condition=model.phy.LabCondition.getCondition;              
              Condition.setValue('magnetic_field',MagneticField);

              %calculate cluster coherence              
              clst_index=cluster_index_list{n,1};
              clst_coh=model.phy.Solution.CCESolution.CCECoherenceStrategy.SSCCEClusterCoherence(clst_index,clst_para);
              CoherenceMatrix(n,:)=clst_coh.calculate_cluster_coherence(evolu_para);
              delete(clst_coh);
          end
          delete(gcp('nocreate'));

           toc
           disp('calculation of the cluster-coherence matrix finished.');          

           obj.calculate_tilde_coh_matrix(CoherenceMatrix,cluster_iter);
            
           if ncluster<20000
                obj.keyVariables('coherence_matrix')=CoherenceMatrix;
           else
                timeTag=datestr(clock,'yyyymmdd_HHMMSS');
                save([OUTPUT_FILE_PATH, 'coherence_matrix', timeTag, '.mat'],'CoherenceMatrix');
                clear CoherenceMatrix;
           end
         end
        
    end
    
end


classdef WavefunctionEvolution < model.phy.Dynamics.AbstractEvolutionKernel
    %WAVEFUNCTIONEVOLUTION Summary of this class goes here
    % This class is used to calculate the time evolutuion of the
    % wavefuncton in the Hilbert Spase
    
    properties
        timelist
        initial_state_type
        matrix_prefactor
        matrixList
        result
    end
    
    methods
        function obj=WavefunctionEvolution(qoperatorList,initial_state_type, prefactor)
              noperator=length(qoperatorList);               
              obj.matrixList=cell(1,noperator);
              for n=1:noperator
                   obj.matrixList{1,n}=qoperatorList{n}.getMatrix();
              end
              
              obj.initial_state_type=initial_state_type;
              
              if nargin>2
                   obj.matrix_prefactor=prefactor;
              else
                   obj.matrix_prefactor=-1i;
              end             
              obj.result=0;
        end
          
        function state_out = calculate_evolution(obj, state_in, timelist)
               obj.timelist=timelist;
               ntime=length(timelist);
               dt=timelist(2)-timelist(1);

               noperator=length(obj.matrixList);
               core_mat_list=cell(1,noperator);
               evolution_mat_list=cell(1,noperator);
               for m=1:noperator %initial the core matrice                  
                   core_mat_list{m}=expm(obj.matrix_prefactor(m)*dt*obj.matrixList{1,m});
                   evolution_mat_list{m}=1;
               end
               state_out=cell(1,ntime);
               state_out{1,1}=state_in;
               for n=2:ntime                   
                   for m=1:noperator                      
                     evolution_mat_list{1,m}=evolution_mat_list{1,m}*core_mat_list{1,m};
                   end
                   wave_vec=state_in;
                  for m=1:noperator
                       wave_vec=evolution_mat_list{1,m}*wave_vec;
                  end
                  state_out{1,n}=wave_vec;
               end
                obj.result=state_out;
        end
                   

        function  mean_val=mean_value(obj, obs_list,varargin)
              ntime=length(obj.timelist);
              len_obs=length(obs_list);
              mean_val=zeros(len_obs,ntime);
              
              for n=1:len_obs
                  mat=obs_list{n}.getMatrix;
                  if nargin>2
                      left_vec=varargin{1}{1};
                      mean_val(n,:)=cellfun(@(s) left_vec*mat*s, obj.result);
                  else
                      mean_val(n,:)=cellfun(@(s) s'*mat*s, obj.result);
                  end
              end
          end
          
    end
    
end


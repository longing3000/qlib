classdef MatrixVectorEvolution < model.phy.Dynamics.AbstractEvolutionKernel
    %MATRIXVECTOREVOLUTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        matrix
        matrix_cell
        result
        matrix_prefactor
        initial_state_type
    end
    
    methods
        function obj=MatrixVectorEvolution(qOperators, initial_state_type, prefactors)
            try
                if iscell(qOperators)
                    noperators=length(qOperators);               
                    mat_cell=cell(1,noperators);
                    for kk=1:noperators
                        mat_cell{kk}=qOperators{kk}.getMatrix();              
                    end
                    obj.matrix_cell=mat_cell;
                    obj.matrix=[];
                else
                    obj.matrix=qOperators.getMatrix;
                    obj.matrix_cell=[];
                end
            catch
                error([class(qOperators), 'does not have a property of matrix']);
            end
            obj.result=[];
            
            obj.initial_state_type=initial_state_type;
            
            if nargin > 2
                obj.matrix_prefactor= prefactors;
            else
                obj.matrix_prefactor= -1.j;
            end
        end
        
        function state_out=evolve_step(obj, state_in, dt)
            if isempty(obj.matrix)
                noperator=length(obj.matrix_cell);
                nsection=length(obj.matrix_prefactor);
                state=state_in;
                for kk=1:nsection
                   op_idx=mod(kk,noperator);
                   if op_idx==0
                      op_idx=noperator; 
                   end
                   state=expv(obj.matrix_prefactor(kk)*dt, obj.matrix_cell{op_idx}, state); 
                end
                state_out=state;
            else
                state_out=expv(obj.matrix_prefactor*dt, obj.matrix, state_in);
            end
        end
        
        function state_out=calculate_evolution(obj, state_in, time_list)
            nInterval=length(time_list)-1;
            
            state=state_in;
            for k=1:nInterval
                fprintf('calculating evolution from t=%e to t=%e...\n', time_list(k), time_list(k+1));
                dt=time_list(k+1)-time_list(k);
                state_out=obj.evolve_step(state, dt);                
                state=state_out;
                
                obj.result=[obj.result, state_out];
            end
        end
        
        function mean_val=mean_value(obj, obs_list)
            [len_res_dim, len_res]=size(obj.result);
            len_obs=length(obs_list);
            
            if strcmp(obj.initial_state_type, 'MixedState')
                obs_mat=zeros(len_obs, len_res_dim);
                for k=1:len_obs
                    obs_mat(k,:)=obs_list{k}.getVector()';
                end
                mean_val=obs_mat*obj.result;
                
            elseif strcmp(obj.initial_state_type, 'PureState')
                mean_val=zeros(len_obs, len_res);
                for ii=1:len_obs
                    mat=obs_list{ii}.getMatrix();
                    mat_on_state=mat*obj.result;
                    mean_val=sum(conj(obj.result).*mat_on_state);
                end
            else
                error('not supported.');
            end
        end
        
    end
    
end


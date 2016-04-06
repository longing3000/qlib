classdef AbstractDynamics < handle
    %ABSTRACTDYNAMICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        is_evolved
        observable_list
        observable_keys
        observable_values
        time_list
        state_in
        state_out
        
        kernel
        render
    end
    
    methods
        function obj=AbstractDynamics(kernel)
            obj.kernel=kernel;
            obj.observable_list={};
            obj.observable_keys={};
            obj.is_evolved=0;
            
            obj.render=view.SolutionRender.DynamicsRender(obj);
        end
        
        function evolve(obj)
            obj.state_out=obj.kernel.calculate_evolution(obj.state_in, obj.time_list);
            obj.is_evolved=1;
        end
        
        function calculate_mean_values(obj,varargin)
            if ~obj.is_evolved
                obj.evolve();
            end
            obj.observable_values=obj.kernel.mean_value(obj.observable_list,varargin);
        end
        
        
        function set_initial_state(obj, state,varargin)
            if nargin>2
                Space=varargin{1};
                switch Space
                    case 'Liouville'
                        obj.state_in=state.getVector();
                        obj.kernel.result=state.getVector();
                    case 'Hilbert'
                        obj.state_in=state.getMatrix();
                    otherwise
                        error([class(state), 'does not have a property of vector.']);
                end
            else
                obj.state_in=state;
            end
        end
        
        function set_time_sequence(obj, time_list)
            obj.time_list=time_list;
        end
        
        function addObervable(obj, obs_list)
            obj.observable_list=[obj.observable_list,obs_list];
            obs_num=length(obj.observable_list);
            
            obj.observable_keys=cell(1, obs_num);
            for k=1:obs_num
                obj.observable_keys{k}=obs_list{k}.name;
            end
        end
        
        
    end
    
end


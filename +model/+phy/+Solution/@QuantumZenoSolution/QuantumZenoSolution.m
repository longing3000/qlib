classdef QuantumZenoSolution < model.phy.Solution.AbstractSolution
    %XYMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj=QuantumZenoSolution(xml_file)
            obj@model.phy.Solution.AbstractSolution(xml_file);
        end
        
        function get_parameters(obj, p)
            get_parameters@model.phy.Solution.AbstractSolution(obj, p);
            
            obj.parameters.SpinCollectionStrategy = p.get_parameter('SpinCollection', 'Source');
            obj.parameters.spinType = p.get_parameter('SpinCollection', 'SpinType');
            obj.parameters.nspin  = p.get_parameter('SpinCollection', 'SpinNum');
            
            obj.parameters.AddOnSite1 = p.get_parameter('Interaction1', 'AddOnSite');
            obj.parameters.onSite1 = p.get_parameter('Interaction1', 'OnSite');
            obj.parameters.AddDqtInt1 = p.get_parameter('Interaction1', 'AddDqtInt');
            obj.parameters.dqtInt1 = p.get_parameter('Interaction1', 'DqtInt');
            obj.parameters.AddXYInt1 = p.get_parameter('Interaction1', 'AddXYInt');
            obj.parameters.xyInt1 = p.get_parameter('Interaction1', 'XYInt');
            obj.parameters.AddDipInt1 = p.get_parameter('Interaction1', 'AddDipInt');
            obj.parameters.dipInt1 = p.get_parameter('Interaction1', 'DipInt');
            
            obj.parameters.AddOnSite2 = p.get_parameter('Interaction2', 'AddOnSite');
            obj.parameters.onSite2 = p.get_parameter('Interaction2', 'OnSite');
            obj.parameters.AddDqtInt2 = p.get_parameter('Interaction2', 'AddDqtInt');
            obj.parameters.dqtInt2 = p.get_parameter('Interaction2', 'DqtInt');
            obj.parameters.AddXYInt2 = p.get_parameter('Interaction2', 'AddXYInt');
            obj.parameters.xyInt2 = p.get_parameter('Interaction2', 'XYInt');
            obj.parameters.AddDipInt2 = p.get_parameter('Interaction2', 'AddDipInt');
            obj.parameters.dipInt2 = p.get_parameter('Interaction2', 'DipInt');
            
            obj.parameters.TimeList = p.get_parameter('Dynamics',  'TimeList');
            obj.parameters.nsection = p.get_parameter('Dynamics',  'SectionNum');
            obj.parameters.proportion = p.get_parameter('Dynamics',  'Proportion');
            
            %%iniital state
            tp=p.get_parameter('InitialState', 'Type');
            switch tp
                case 'MixedState'
                    st=p.get_parameter('InitialState', 'DensityMatrix');
                case 'PureState'
                    st=p.get_parameter('InitialState', 'StateVector');
                otherwise
                    error('InitialStateType "%s" is not supported.', tp);
            end
            obj.parameters.InitialStateType = tp;
            obj.parameters.InitialState = st;
            
            %Clustering states
                                    
            obj.parameters.load_cluster_iter=p.get_parameter('Clustering','LoadCluterIterator');
            if obj.parameters.load_cluster_iter
                CluterIteratorName=p.get_parameter('Clustering','CluterIteratorName');
                data=load([OUTPUT_FILE_PATH,CluterIteratorName]);
                obj.keyVariables('cluster_iterator')=data.cluster_iterator;
                disp('cluster_iterator loaded.');
            else
                obj.parameters.CutOff=p.get_parameter('Clustering', 'CutOff');
                obj.parameters.MaxOrder=p.get_parameter('Clustering', 'MaxOrder');
            end
            
            
            
            %%Observables
            nObs=p.get_parameter('Observable', 'ObservableNumber');
            obs_name=cell(1,nObs);  obs_str=cell(1,nObs);
            for k=1:nObs
                obs_name{k} = p.get_parameter('Observable', ['ObservableName',num2str(k)]);
                obs_str{k} = p.get_parameter('Observable',  ['ObservableString',num2str(k)]);
            end
            
            obj.parameters.ObservableNumber=nObs;
            obj.parameters.ObservableName=obs_name;
            obj.parameters.ObservableString=obs_str;
            
        end
        
        
        function perform(obj)
            import model.phy.QuantumOperator.MatrixStrategy.FromKronProd
            
            spin_collection=obj.GetSpinList();
            obj.keyVariables('spin_collection')=spin_collection;
            
            matrix_strategy=FromKronProd();
            [hami1,hami2, liou1,liou2] = obj.GetHamiltonianLiouvillian(spin_collection,matrix_strategy);
            initial_state              = obj.GetInitialState(spin_collection,matrix_strategy);
            observables                = obj.GetObservables(spin_collection,matrix_strategy);
            obj.StoreKeyVariables(spin_collection, hami1,hami2,liou2,liou2, initial_state)
                                    
            dynamics                   = obj.StateEvolve(hami1,hami1, liou1,liou2, initial_state);
            mean_values                = obj.GetMeanValues(dynamics, observables);
            obj.StoreKeyVariables(dynamics,mean_values);
            
            final_states=dynamics.kernel.result;         
            [~] = obj.GetStateInfo(spin_collection,final_states);
        end
    end
    
end


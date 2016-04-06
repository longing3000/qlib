classdef AbstractCCESolution < model.phy.Solution.AbstractSolution
    %ABSTRACTCCESOLUTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

    end
    
    methods
        function obj=AbstractCCESolution(xml_file)
            obj@model.phy.Solution.AbstractSolution(xml_file);
        end
        function get_parameters(obj, p)
            get_parameters@model.phy.Solution.AbstractSolution(obj, p);
            obj.parameters.LoadSpinCollection=p.get_parameter('SpinCollection','LoadSpinCollection');
            if obj.parameters.LoadSpinCollection
                obj.parameters.SpinCollectionStrategy=p.get_parameter('SpinCollection', 'Source');
                obj.parameters.InputFile=p.get_parameter('SpinCollection', 'FileName');
                obj.BathSpinParameters(p);
                disp('spin collection loaded.');
            else
                obj.BathSpinParameters(p);
            end

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

            obj.CentralSpinParameters(p);                      
            obj.parameters.MagneticField=p.get_parameter('Condition', 'MagneticField');
            obj.parameters.IsSecularApproximation=p.get_parameter('Interaction', 'IsSecular');

            obj.parameters.CCEStrategy=p.get_parameter('Dynamics','CCEStrategy');
            NTime=p.get_parameter('Dynamics', 'NTime');
            TMax=p.get_parameter('Dynamics', 'TMax');
            dt=TMax/(NTime-1);
            obj.parameters.TimeList=0:dt:TMax;
            obj.parameters.NTime=NTime;
            obj.parameters.TMax=TMax;
            obj.parameters.NPulse=p.get_parameter('Dynamics', 'NPulse');
            
            % get the special parameters for different CCE strategies
            obj.get_CCE_special_parameter(p);             
        end
        function BathSpinParameters(obj, p)
             %set bath spin
             % get the parameters for setting the bath spins e.g. ZFS, eta, principle_axis, coordinate
                obj.parameters.SetBathSpins.SetSpin=p.get_parameter('SetBathSpins','SetSpin');
                if obj.parameters.SetBathSpins.SetSpin
                    SpeciesNumber=p.get_parameter('SetBathSpins','SpeciesNumber');
                    BathSpinsSettingCell=cell(1,4);
                    NameList=p.get_parameter('SetBathSpins','Name');
                    ZFSList=p.get_parameter('SetBathSpins','ZFS');
                    etaList=p.get_parameter('SetBathSpins','eta');
                    paxisList=p.get_parameter('SetBathSpins','principle_axis');
                    for k=1:SpeciesNumber
                        BathSpinsSettingCell{k}.name=NameList(k,:);
                        BathSpinsSettingCell{k}.ZFS=ZFSList(k,:);
                        if ~strcmp(etaList,'default')                   
                            BathSpinsSettingCell{k}.eta=etaList(k,:);
                        end
                        if ~strcmp(paxisList,'default')                   
                            BathSpinsSettingCell{k}.principle_axis=paxisList(k,:);
                        end                  
                    end
                    obj.parameters.SetBathSpins.BathSpinsSettingCell=BathSpinsSettingCell;               
                end
                
                %set the gradient field options
                obj.parameters.add_gradient_field=p.get_parameter('SetBathSpins','AddGradientField');
                if obj.parameters.add_gradient_field
                    obj.parameters.field_gradient=p.get_parameter('SetBathSpins','FieldGradient');
                end
                
        end
        
         function CentralSpinParameters(obj,p)
             %set central spin
             % get the parameters for setting the central spin e.g. ZFS, eta, principle_axis, coordinate
                obj.parameters.SetCentralSpin.name=p.get_parameter('SetCentralSpin', 'Name');
                obj.parameters.SetCentralSpin.CentralSpinStates=p.get_parameter('SetCentralSpin', 'CentralSpinStates');
                obj.parameters.SetCentralSpin.SetSpin=p.get_parameter('SetCentralSpin', 'SetSpin');
                Orientation=p.get_parameter('SetCentralSpin', 'Orientation');
                Isotope=p.get_parameter('SetCentralSpin', 'Isotope');
                Coordinate=p.get_parameter('SetCentralSpin', 'Coordinate');
                if obj.parameters.SetCentralSpin.SetSpin
                    if ~strcmp(Orientation,'default') 
                        obj.parameters.SetCentralSpin.orientation=Orientation;
                    end
                    if ~strcmp(Isotope,'default') 
                        obj.parameters.SetCentralSpin.isotope=Isotope;
                    end
                    if ~strcmp(Coordinate,'default') 
                        obj.parameters.SetCentralSpin.coordinate=Coordinate;
                    end
                end
         end
         
 
         
    end
    
end


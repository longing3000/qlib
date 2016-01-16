classdef ChainNeighbouringIterator < model.phy.SpinCollection.SpinCollectionIterator
    %CHAINNEIGHBOURITERATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj=ChainNeighbouringIterator(spin_collection,index_generator)
            obj@model.phy.SpinCollection.SpinCollectionIterator(spin_collection);
            if nargin>1
                obj.index_generator=index_generator;
                obj.index_list=obj.ring_index_gen();
            end
        end
        function res = index_gen(obj)
            nspin=obj.spin_collection.getLength();
            res=cell(nspin-1, 1);
            for ii=1:nspin-1
                res{ii}=[ii, ii+1];
            end
        end        
        function res = ring_index_gen(obj)
            nspin=obj.spin_collection.getLength();
            res=cell(nspin, 1);
            for ii=1:nspin-1
                res{ii}=[ii, ii+1];
            end
            res{nspin}=[nspin,1];
        end
    end
    
end


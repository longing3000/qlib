classdef Spin < model.phy.PhysicalObject.PhysicalObject
    %SPIN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        gamma
        ZFS=0;% zero field splitting result from the inhomogeneous electric field ,i.e., electric quadrupole interaction, ZFS for NV etc
        eta=0;   %asymmetry parameter of ZFS
        principle_axis=[0,0,1];%the orientation of the principal axis
        
        coordinate=[];
        local_field=[0 0 0];
        qAxis=[1 0 0; 0 1 0; 0 0 1];
        
        S
        S2
        
        sx_mat=[]
        sy_mat=[]
        sz_mat=[]
        self_hamiltonian
        eigen_vect
        eigen_val
    end
    
    methods
        function obj = Spin(name, coord)
            obj.name=name;
            [obj.dim, obj.gamma] = model.phy.data.NMRData.get_spin(name);
            
            if nargin > 1
                obj.coordinate=coord;
            end
          
            obj.S= 0.5*(obj.dim-1);
            obj.S2= obj.S*(obj.S+1);
            
            obj.sx_mat=Sx(obj.dim);
            obj.sy_mat=Sy(obj.dim);
            obj.sz_mat=Sz(obj.dim);
        end
        
        function set_spin(obj,para)
            if nargin>1
                fname=fieldnames(para);
                n_fname=length(fname);
                for k=1:n_fname
                   obj.(fname{k})=para.(fname{k}); 
                end
            end
            obj.self_hamiltonian=obj.selfHamiltonian();
            [obj.eigen_vect, obj.eigen_val]=eig(full(obj.self_hamiltonian.getMatrix));
        end
        
        function xmat=sx(obj)
%             if isempty(obj.sx_mat)
%                 obj.sx_mat=Sx(obj.dim);
%                 xmat=obj.sx_mat;
%             else
                xmat=obj.sx_mat;
%             end
        end
        
        function ymat=sy(obj)
%             if isempty(obj.sy_mat)
%                 obj.sy_mat=Sy(obj.dim);
%                 ymat=obj.sy_mat;
%             else
                ymat=obj.sy_mat;
%             end
        end
        
        function zmat=sz(obj)
%             if isempty(obj.sz_mat)
%                 obj.sz_mat=Sz(obj.dim);
%                 zmat=obj.sz_mat;
%             else
                zmat=obj.sz_mat;
%             end
        end
        
        function projMat=p(obj, k)
            projMat=zeros(obj.dim);
            projMat(k, k)=1;
        end
        
        function generalMat=mat(obj, m)
            if nargin < 2
               generalMat=speye(obj.dim)/obj.dim; 
               return;
            end
            
            if length(m)==obj.dim
                generalMat=m;
            else
                error('dimension mismatch. matrix_dim=%d is assigned, but spin_dimension=%d is needed.', length(m), obj.dim);
            end
            
        end
        
        function expM=expMat(obj, m)
            if length(m)==obj.dim
                expM=expm(m);
            else
                error('dimension mismatch. matrix_dim=%d is assigned, but spin_dimension=%d is needed.', length(m), obj.dim);
            end
        end
        
        function mat=eigenVectors(obj)
            mat=obj.eigen_vect;
        end
        
        function mat=eigenValues(obj)
            mat=obj.eigen_val;
        end

        
        function hami=selfHamiltonian(obj)
            import model.phy.SpinInteraction.ZeemanInteraction
            
            sc=model.phy.SpinCollection.SpinCollection();
            sc.spin_source = model.phy.SpinCollection.Strategy.FromSpinList({obj});
            sc.generate();

            hami=model.phy.QuantumOperator.SpinOperator.Hamiltonian(sc);
            hami.addInteraction( ZeemanInteraction(sc) );
            hami.generate_matrix();
        end
        
        function IST_state=IST(obj,state) %irreducible spherical tensors
            basis=SphericalTensor(obj.dim);
            IST_state=basis{state+1};
            
        end
        
        function ISTmat=ISTmatrix(obj)
            dim=obj.dim;
            ISTmat=zeros(dim^2,dim^2);
            basis=SphericalTensor(obj.dim);
            for kk=1:dim^2
                IST_state=basis{kk};
                ISTmat(kk,:)=IST_state(:)';                
            end
            ISTmat=sparse(ISTmat);
        end
        
    end
    
end


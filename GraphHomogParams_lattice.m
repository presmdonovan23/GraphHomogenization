classdef GraphHomogParams_lattice

    properties
        
        dim
        geometry
        D0
        rho
        m
        rate
        rateCoeffs
        diagJumps % dont need to supply, default is 0
        ctr
        
    end
        
    properties (Dependent)
        R
        h
    end
    
    methods
        
        function obj = GraphHomogParams_lattice(varargin)
            
            if length(varargin) == 1
                p = varargin{1};
                
                dim = p.dim;
                geometry = p.geometry;
                D0 = p.D0;
                rho = p.rho;
                m = p.m;
                rate = p.rate;
                rateCoeffs = p.rateCoeffs;
                diagJumps = p.diagJumps;
                ctr = p.ctr;
                
            else
                
                dim = varargin{1};
                geometry = varargin{2};
                D0 = varargin{3};
                rho = varargin{4};
                m = varargin{5};
                rate = varargin{6};
                rateCoeffs = varargin{7};
                diagJumps = varargin{8};
                if nargin < 9 || isempty(varargin{9})
                    ctr = .5;
                else
                    ctr = varargin{9};
                end
                
            end
            
            obj.dim = dim;
            obj.geometry = geometry;
            obj.D0 = D0;
            obj.rho = rho;
            obj.m = m;
            obj.rate = rate;
            obj.rateCoeffs = rateCoeffs;
            obj.diagJumps = diagJumps;
            obj.ctr = ctr;
            
        end
        
        function R = get.R(obj)
            R = obj.rho/2;
        end
        
        function h = get.h(obj)
            h = 1./obj.m;
        end
    end
    
end
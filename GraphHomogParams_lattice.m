classdef GraphHomogParams_lattice
    % An auxiliary class for creating a GraphHomogInput object when the
    % following conditions are satisfied:
    %   -Dimension is 2 or 3
    %   -Obstruction is a square or a circle
    %   -Graph invariant under integer vector scaling (periodicity 1)
    % Inputs:
    %   1) dim (scalar) = dimension (2 or 3)
    %   2) geometry (string) = options:
    %       "circle", 
    %       "circleSlowdown", 
    %       "square", 
    %       "squareBonding", 
    %       "squareBdyRepel", 
    %       "squareBdyAttract", 
    %       "squareBdySlow"
    %   3) D0 (scalar) = free diffusivity (typically 1)
    %   4) rho (scalar) = 2*(a+R)/L where a = particle radius, R = obstruction
    %   radius, and L = obstruction spacing
    %   5) m (scalar) = 1/(mesh spacing) = 1/(path length)
    %   6) rate = currently obsolete field
    %   7) diagJumps (scalar) = 0 if no diagonal jumps, 1 if diagonal jumps
    %   8) ctr = center of obstruction (center of unit cell is (.5,...,.5)
    %
    % Other fields:
    %   1) rateCoeffs (object) = scalar fields containing certain
    %   rate function parameters:
    %       .dist = specifies distance for the squareBonding,
    %       squareBdyRepel, squareBdyAttract, and squareBdySlow geometries
    %       .alpha = multiplier on drift term for rate (0 = no drift, 1 =
    %       repulsion, -1 = attraction
    %       .K1 = K1 coefficient for drift
    %       .K2 = K2 coefficient for drift
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
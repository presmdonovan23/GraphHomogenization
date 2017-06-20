classdef GraphHomogInput

    properties
        
        L
        nodes
        edges
        edgeRates
        edgeJumps
        nodeInds
        
    end
    
    methods
        
        function obj = GraphHomogInput(varargin)
            
            if nargin == 0
                
                obj.L = [];
                obj.nodes = [];
                obj.edges = [];
                obj.edgeRates = [];
                obj.edgeJumps = [];
                obj.nodeInds = [];
                
            elseif isa(varargin{1},'GraphHomogParams_lattice')
                
                [L,nodes,edges,edgeRates,edgeJumps,nodeInds] = ...
                    latticeSetup(varargin{1});
                
                obj.L = L;
                obj.nodes = nodes;
                obj.edges = edges;
                obj.edgeRates = edgeRates;
                obj.edgeJumps = edgeJumps;
                obj.nodeInds = nodeInds;
                
            elseif length(varargin) == 1
                p = varargin{1};
                obj.L = p.L;
                obj.nodes = p.nodes;
                obj.edges = p.edges;
                obj.edgeRates = p.edgeRates;
                obj.edgeJumps = p.edgeJumps;
                obj.nodeInds = p.nodeInds;
                
            elseif length(varargin) == 6
                obj.L = varargin{1};
                obj.nodes = varargin{2};
                obj.edges = varargin{3};
                obj.edgeRates = varargin{4};
                obj.edgeJumps = varargin{5};
                obj.nodeInds = varargin{6};
            else
                error('Incorrect inputs.');
            end
            
        end
        
    end
    
end
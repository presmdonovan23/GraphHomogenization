classdef GraphHomogInput
    % An object containing the minimum necessary inputs to the
    % getDeff_homog function.
    % Inputs:
    %   1) L = an s x s square "rate" matrix where s = number of nodes in 
    %   quotient node set. L(i,j) contains the jump rate from node i to 
    %   node j.
    %   2) nodes = an s x d array where d is the dimension. nodes(i,:)
    %   contains the position of the ith node.
    %   3) edges = an e x 2 array where e = number of edges in the quotient
    %   edge set. edges(i,1) = starting node index of the ith edge.
    %   edges(i,2) = ending node index of the ith edge.
    %   4) edgeRates = an e x 1 array where edgeRates(i) = jump rate of the
    %   ith edge.
    %   5) edgeJumps = an e x d array where edgeJumps(i) = the jump of the
    %   ith edge. Equivalently, edgeJumps(i) = nodes(edges(i,2),:) -
    %   nodes(edges(i,1),:).
    
    properties
        
        L
        nodes
        edges
        edgeRates
        edgeJumps
        
    end
    
    methods
        
        function obj = GraphHomogInput(varargin)
            
            if nargin == 0
                % Return empty object
                obj.L = [];
                obj.nodes = [];
                obj.edges = [];
                obj.edgeRates = [];
                obj.edgeJumps = [];
                
            elseif isa(varargin{1},'GraphHomogParams_lattice')
                % Argument is a GraphHomogParams_lattice object
                [L,nodes,edges,edgeRates,edgeJumps] = ...
                    latticeSetup(varargin{1});
                
                obj.L = L;
                obj.nodes = nodes;
                obj.edges = edges;
                obj.edgeRates = edgeRates;
                obj.edgeJumps = edgeJumps;
                
            elseif length(varargin) == 1
                % Argument is an object containing necessary fields.
                p = varargin{1};
                obj.L = p.L;
                obj.nodes = p.nodes;
                obj.edges = p.edges;
                obj.edgeRates = p.edgeRates;
                obj.edgeJumps = p.edgeJumps;
                
            elseif length(varargin) == 5
                % Five arguments passed (not as an object)
                obj.L = varargin{1};
                obj.nodes = varargin{2};
                obj.edges = varargin{3};
                obj.edgeRates = varargin{4};
                obj.edgeJumps = varargin{5};
            else
                error('Incorrect inputs.');
            end
            
        end
        
    end
    
end
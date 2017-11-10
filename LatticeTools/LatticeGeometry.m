classdef LatticeGeometry < handle
    properties
        dim
        m
        name
        obRad
        
        obCtr
        diagJumps

        specialSetting
        
        driftMult
        driftDecay
        
        obSlowdownFctr
        bdyDist
    end
    
    properties (Dependent)
        h
        sideLen
        isValid
    end
    
    properties (Access = protected)
        validDims = [2 3];
        validNames = {'circle','square'};
        validDiagJumps = [0 1 2];
        validSpecialSettings = {'none',...
                                'slowdown',...
                                'bdyBonding',...
                                'bdyAttractRepel',...
                                'm2_slowOneSite',...
                                'bdySlow'};
        % settings where obSlowdownFctr can be defined
        validSlowdownSettings = {'slowdown',...
                                 'bdyBonding',...
                                 'bdyAttractRepel',...
                                 'm2_slowOneSite'};
        
        % settings where bdyDist can be defined
        validBdyDistSettings = {'bdyBonding'};
        
        errorCount
    end
    
    methods
        function value = get.isValid(this)
            value = this.errorCount == 0;
        end
        
        function value = get.h(this)
            value = 1/this.m;
        end
        
        function value = get.sideLen(this)
            if strcmpi(this.name,'square')
                value = 2*this.obRad;
            else
                value = [];
            end
        end 
    end
    
    methods
        function this = LatticeGeometry(dim, m, name, obRad, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist)
            
            if nargin < 4
                error('Must supply dim, m, name, and obstruction radius.');
            end
            
            if nargin < 5 || isempty(obCtr)
                if m == 2
                    obCtr = .75*ones(1,dim);
                else
                    obCtr = .5*ones(1,dim);
                end
            end
            
            if nargin < 6 || isempty(diagJumps)
                diagJumps = 0;
            end
            
            if nargin < 7 || isempty(specialSetting)
                specialSetting = 'none';
            end
            
            if nargin < 8 || isempty(driftMult)
                driftMult = 0;
            end
            
            if nargin < 9 || isempty(driftDecay)
                if driftMult ~= 0
                    driftDecay = 10;
                else
                    driftDecay = 0;
                end
            end
            
            if nargin < 10 || isempty(obSlowdownFctr)
                if ismember(specialSetting,this.validSlowdownSettings)
                    obSlowdownFctr = .5;
                else
                    obSlowdownFctr = [];
                end
            end
            
            if nargin < 11 || isempty(bdyDist)
                if ismember(specialSetting,this.validBdyDistSettings)
                    bdyDist = 1/m;
                else
                    bdyDist = [];
                end
            end
            
            if length(obCtr) == 1
                obCtr = obCtr*ones(1,dim);
            end
            if length(obCtr) ~= dim
                error('obCtr must be a scalar or have correct dimension.')
            end
            
            this.dim = dim;
            this.m = m;
            this.name = name;
            this.obRad = obRad;
            this.obCtr = obCtr;
            this.diagJumps = diagJumps;
            this.specialSetting = specialSetting;

            this.driftMult = driftMult;
            this.driftDecay = driftDecay;

            this.obSlowdownFctr = obSlowdownFctr;
            this.bdyDist = bdyDist;
            
            this.validate;
        end
        
        function this = validate(this)
            errors = [];
            
            % dim
            if ~ismember(this.dim,this.validDims)
                errors{end+1} = 'Invalid dim.';
            end
            % m
            if this.m < 2 || rem(this.m,1) ~= 0
                errors{end+1} = 'm must be an integer >= 2.';
            end
            % name
            if ~ismember(this.name,this.validNames)
                errors{end+1} = 'Invalid name.';
            end
            % obRad
            if this.obRad < 0
                errors{end+1} = 'obRad must be non-negative.';
            end
            if this.m == 2 && this.obRad >= 1/2
                errors{end+1} = 'obRad must be < 1/2 when m = 2.';
            end
            if strcmpi(this.specialSetting,'m2_slowOneSIte') && (this.obRad <= 0 || this.obRad >= 1/2)
                errors{end+1} = 'obRad must be in (0,1/2) when specialSetting = ''m2_slowOneSite''.';
            end
            % obCtr
            if any(this.obCtr < 0) || any(this.obCtr > 1)
                errors{end+1} = 'obCtr must be in [0,1]^d.';
            end
            if this.m == 2 && norm(this.obCtr - .75) > eps
                errors{end+1} = 'obCtr must be (3/4,3/4) when m = 2.';
            end
            % diagJumps
            if ~ismember(this.diagJumps,this.validDiagJumps)
                errors{end+1} = 'diagJumps must be 0, 1, or 2.';
            end
            if this.diagJumps == 2 && (~strcmpi(this.name,'square') || this.dim ~= 2)
                errors{end+1} = 'Corrected diagonal jumps only implemented for square geometry in 2D.';
            end
            % specialSetting
            if ~ismember(this.specialSetting,this.validSpecialSettings)
                errors{end+1} = 'Invalid specialSetting.';
            end
            if contains(this.specialSetting,'m2') && this.m ~= 2
                errors{end+1} = 'Can only use m2 specialSetting if m = 2.';
            end
            if ismember(this.specialSetting,this.validBdyDistSettings) && ...
                    ~strcmpi(this.name,'square')
                errors{end+1} = 'Can only use square geometry for this specialSetting.';
            end
            % driftMult
            if this.driftMult ~= 0 && (this.diagJumps >= 1 || strcmpi(this.name,'square'))
                errors{end+1} = 'Drift only implemented for circle geometry with no diagonal jumps.';
            end
            if this.driftMult ~= 0 && ~strcmpi(this.specialSetting,'none')
                errors{end+1} = 'Drift should be 0 for special settings. Might work fine, haven''t tested.';
            end
            % driftDecay
            if this.driftDecay < 0
                errors{end+1} = 'driftDecay must be non-negative.';
            end
            % obSlowdownFctr
            if this.obSlowdownFctr <= 0
                errors{end+1} = 'obSlowdownFctr must be positive.';
            end
            if ~isempty(this.obSlowdownFctr) && ~ismember(this.specialSetting,this.validSlowdownSettings)
                errors{end+1} = 'obSlowdownFctr must be empty for this specialSetting.';
            end
            if isempty(this.obSlowdownFctr) && ismember(this.specialSetting,this.validSlowdownSettings)
                errors{end+1} = 'obSlowdownFctr must be defined for this specialSetting.';
            end
            % bdyDist
            if this.bdyDist <= 0
                errors{end+1} = 'bdyDist must be positive.';
            end
            if ~isempty(this.bdyDist) && ~ismember(this.specialSetting,this.validBdyDistSettings)
                errors{end+1} = 'bdyDist must be empty for this specialSetting.';
            end
            if isempty(this.bdyDist) && ismember(this.specialSetting,this.validBdyDistSettings)
                errors{end+1} = 'bdyDist must be defined for this specialSetting.';
            end
            
            this.errorCount = length(errors);
            if this.errorCount > 0
                fprintf('%d errors:\n',this.errorCount);
                for i = 1:this.errorCount
                    fprintf('\t%s\n',errors{i});
                end
            end
        end
        
    end
end
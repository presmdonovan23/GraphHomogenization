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
                                'm2_blockOneSite',...
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
                obCtr = .5*ones(1,dim);
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

            count = 0;
            % dim
            if ~ismember(this.dim,this.validDims)
                count = count+1;
                errors{count} = 'Invalid dim.';
            end
            % m
            if this.m < 2 || rem(this.m,1) ~= 0
                count = count+1;
                errors{count} = 'm must be an integer >= 2.';
            end
            % name
            if ~ismember(this.name,this.validNames)
                count = count+1;
                errors{count} = 'Invalid name.';
            end
            % obRad
            if this.obRad < 0
                count = count+1;
                errors{count} = 'obRad must be non-negative.';
            end
            % obCtr
            if any(this.obCtr < 0) || any(this.obCtr > 1)
                count = count+1;
                errors{count} = 'obCtr must be in [0,1]^d.';
            end
            % diagJumps
            if ~ismember(this.diagJumps,this.validDiagJumps)
                count = count+1;
                errors{count} = 'diagJumps must be 0, 1, or 2.';
            end
            if this.diagJumps == 2 && (~strcmpi(this.name,'square') || this.dim ~= 2)
                count = count+1;
                errors{count} = 'Corrected diagonal jumps only implemented for square geometry in 2D.';
            end
            % specialSetting
            if ~ismember(this.specialSetting,this.validSpecialSettings)
                count = count+1;
                errors{count} = 'Invalid specialSetting.';
            end
            if contains(this.specialSetting,'m2') && this.m ~= 2
                count = count+1;
                errors{count} = 'Can only use m2 specialSetting if m = 2.';
            end
            if ismember(this.specialSetting,this.validBdyDistSettings) && ...
                    ~strcmpi(this.name,'square')
                count = count+1;
                errors{count} = 'Can only use square geometry for this specialSetting.';
            end
            % driftMult
            if this.driftMult ~= 0 && (this.diagJumps >= 1 || strcmpi(this.name,'square'))
                count = count+1;
                errors{count} = 'Drift only implemented for circle geometry with no diagonal jumps.';
            end
            if this.driftMult ~= 0 && ~strcmpi(this.specialSetting,'none')
                count = count+1;
                errors{count} = 'Drift should be 0 for special settings. Might work fine, haven''t tested.';
            end
            % driftDecay
            if this.driftDecay < 0
                count = count+1;
                errors{count} = 'driftDecay must be non-negative.';
            end
            % obSlowdownFctr
            if this.obSlowdownFctr <= 0
                count = count+1;
                errors{count} = 'obSlowdownFctr must be positive.';
            end
            if ~isempty(this.obSlowdownFctr) && ~ismember(this.specialSetting,this.validSlowdownSettings)
                count = count+1;
                errors{count} = 'obSlowdownFctr must be empty for this specialSetting.';
            end
            if isempty(this.obSlowdownFctr) && ismember(this.specialSetting,this.validSlowdownSettings)
                count = count+1;
                errors{count} = 'obSlowdownFctr must be defined for this specialSetting.';
            end
            % bdyDist
            if this.bdyDist <= 0
                count = count+1;
                errors{count} = 'bdyDist must be positive.';
            end
            if ~isempty(this.bdyDist) && ~ismember(this.specialSetting,this.validBdyDistSettings)
                count = count+1;
                errors{count} = 'bdyDist must be empty for this specialSetting.';
            end
            if isempty(this.bdyDist) && ismember(this.specialSetting,this.validBdyDistSettings)
                count = count+1;
                errors{count} = 'bdyDist must be defined for this specialSetting.';
            end
            
            this.errorCount = count;
            if this.errorCount > 0
                fprintf('%d errors:\n',this.errorCount);
                for i = 1:this.errorCount
                    fprintf('\t%s\n',errors{i});
                end
            end
        end
        
    end
end
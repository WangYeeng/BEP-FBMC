%==========================================================================
% Ying Wang, wangyingstu@163.com
% (c) 2024 by Ying Wang 
%==========================PAM Constellation===============================
classdef SignalConstellation < handle 
  properties (SetAccess = private)
       Method
       ModulationOrder
       BitMapping
       SymbolMapping
       Implementation
  end
  
  methods
  	function obj = SignalConstellation(varargin)      
        obj.ModulationOrder = varargin{1};
        obj.Method          = varargin{2};      
        %% Gray coded bitmapping
        if strcmp(obj.Method,'PAM')
            obj.BitMapping = [ones(obj.ModulationOrder/2,1);zeros(obj.ModulationOrder/2,1)];
            for i_temp     = 2:log2(obj.ModulationOrder)
                BinaryTemp               = obj.BitMapping(1:2:end,i_temp-1);
                obj.BitMapping(:,i_temp) = [BinaryTemp;BinaryTemp(end:-1:1)];
            end
            obj.SymbolMapping = (2*(1:obj.ModulationOrder)-obj.ModulationOrder-1).';
            obj.SymbolMapping = obj.SymbolMapping/sqrt(mean(abs(obj.SymbolMapping).^2));
        else
           error('Signal constellation method must be PAM!');
        end
        [~,SortOrder]     = sort(bi2de(obj.BitMapping),'ascend');
        obj.SymbolMapping = obj.SymbolMapping(SortOrder);
        obj.BitMapping    = obj.BitMapping(SortOrder,:);
    end   
  end
end
      
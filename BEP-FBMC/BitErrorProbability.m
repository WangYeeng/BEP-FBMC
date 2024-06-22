%==========================================================================
% Ying Wang, wangyingstu@163.com
% (c) 2024 by Ying Wang 
%================Bit Error Probability============
% This function is used to calculate the bit error probability for 
% both QAM and PAM constellations.
%==========================================================================
function BitErrorProbability = BitErrorProbability(...
    SNR_dB, ...        % SNR .
    SymbolAlphabet, ...% The Symbol Alphabet with mean(SymbolAlphabet.*conj(SymbolAlphabet))==1.
    BitMapping)        % The bitmapping corresponding to the symbol Alphabet.
    
   
    % For the decision regions we assume a rectangular regular grid!
    % Calculate the half interval of the decision region
    HalfDecisionInterval = min(abs(real(SymbolAlphabet)));
    % Extracting the real and imaginary parts of the symbol mapping
    realPart = real(SymbolAlphabet);
    imagPart = imag(SymbolAlphabet);
    % Initialize the decision region matrix
    DecisionRegions = zeros(length(SymbolAlphabet), 4);
    % Calculate the decision region
    for i = 1:length(SymbolAlphabet)
        DecisionRegions(i, :) = [realPart(i) - HalfDecisionInterval, ...
                                 realPart(i) + HalfDecisionInterval, ...
                                 imagPart(i) - HalfDecisionInterval, ...
                                 imagPart(i) + HalfDecisionInterval];
    end

    % Deal with the boundary situation
    DecisionRegions(realPart == min(realPart), 1) = -inf;
    DecisionRegions(realPart == max(realPart), 2) = +inf;
    DecisionRegions(imagPart == min(imagPart), 3) = -inf;
    DecisionRegions(imagPart == max(imagPart), 4) = +inf;

    BitErrorProbability = nan(length(SNR_dB),1);
    for i_SNR = 1:length(SNR_dB)
        Pn = 10^(-SNR_dB(i_SNR)/10);
        ProbabilityMatrix = nan(size(SymbolAlphabet,1),size(SymbolAlphabet,1));
        for i_symbol = 1:size(SymbolAlphabet,1)
            x = SymbolAlphabet(i_symbol);
            Ey2 = abs(x).^2+Pn;
            Eh2 = 1;
            Eyh = x;
            ProbabilityMatrix(:,i_symbol)=ProbabilityRectangularboundary(Ey2,Eh2,Eyh,DecisionRegions(:,1),DecisionRegions(:,2),DecisionRegions(:,3),DecisionRegions(:,4));
        end 
        ErrorProbability = nan(2,size(BitMapping,2));
        for i_bit= 1:size(BitMapping,2)
            for i_zero_one = [0 1]
                index_x = (BitMapping(:,i_bit)==i_zero_one);
                ErrorProbability(i_zero_one+1,i_bit) = mean(sum(ProbabilityMatrix(not(index_x),index_x)));
            end   
        end
        BitErrorProbability(i_SNR) = mean(mean(ErrorProbability));
    end
end

%==========================================================================
% Calculates the Probability that the one-tap ZF y/h within "(zRlower zIlower] x (zRupper zIupper]". 
%==========================================================================
function Probability = ProbabilityRectangularboundary(...
    Ey2,...             % Expectation{|y|^2}
    Eh2,...             % Expectation{|h|^2}
    Eyh,...             % Expectation{y*conj(h)}
    zRlower,...         % Determines the rectangular region
    zRupper,...
    zIlower,...
    zIupper)

    CDF_Upperboundary = CDF(Ey2,Eh2,Eyh,zRupper,zIupper);
    CDF_Lowerboundary = CDF(Ey2,Eh2,Eyh,zRlower,zIlower);
    CDF_Rightboundary = CDF(Ey2,Eh2,Eyh,zRlower,zIupper);
    CDF_Leftboundary  = CDF(Ey2,Eh2,Eyh,zRupper,zIlower);
    
    Probability = CDF_Upperboundary+CDF_Lowerboundary-CDF_Rightboundary-CDF_Leftboundary;
end

% This function calculates the CDF of the ZF y/h.
function CDF = CDF(...
    Ey2,...             % Expectation{|y|^2}
    Eh2,...             % Expectation{|h|^2}
    Eyh,...             % Expectation{y*conj(h)}
    zR,...              % Real part of the CDF
    zI)                 % Imaginary part of the CDF

Index0      = (zR == -inf)| (zI == -inf);
Index1      = (zR == inf) & (zI ==  inf);
IndexReal   = (zI == inf) & isfinite(zR);
IndexImag   = (zR == inf) & isfinite(zI);
IndexNormal = isfinite(zR)& isfinite(zI);
% Calculates the CDF
CDF_Real   = 1/2-(real(Eyh/Eh2)-zR(IndexReal))./  (2*sqrt((real(Eyh/Eh2)-zR(IndexReal)).^2+(Ey2/Eh2)-abs(Eyh/Eh2).^2));
CDF_Imag   = 1/2-(imag(Eyh/Eh2)-zI(IndexImag))./  (2*sqrt((imag(Eyh/Eh2)-zI(IndexImag)).^2+(Ey2/Eh2)-abs(Eyh/Eh2).^2));
CDF_Normal = 1/4+(zR(IndexNormal)-real(Eyh/Eh2)).*(2*atan((zI(IndexNormal)-imag(Eyh/Eh2))./sqrt((zR(IndexNormal)-real(Eyh/Eh2)).^2+(Ey2/Eh2)-abs(Eyh/Eh2).^2))+pi)./...
                                                  (4*pi*sqrt((zR(IndexNormal)-real(Eyh/Eh2)).^2+(Ey2/Eh2)-abs(Eyh/Eh2).^2))...
                +(zI(IndexNormal)-imag(Eyh/Eh2)).*(2*atan((zR(IndexNormal)-real(Eyh/Eh2))./sqrt((zI(IndexNormal)-imag(Eyh/Eh2)).^2+(Ey2/Eh2)-abs(Eyh/Eh2).^2))+pi)./...
                                                  (4*pi*sqrt((zI(IndexNormal)-imag(Eyh/Eh2)).^2+(Ey2/Eh2)-abs(Eyh/Eh2).^2));
CDF              = nan(size(zR));
CDF(Index0)      = 0;
CDF(Index1)      = 1;
CDF(IndexReal)   = CDF_Real;
CDF(IndexImag)   = CDF_Imag;
CDF(IndexNormal) = CDF_Normal;

end
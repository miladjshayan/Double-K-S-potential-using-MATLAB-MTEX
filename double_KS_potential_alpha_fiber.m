%% Double K-S Potential Along the Alpha-Fiber

% Before running this script, an EBSD variable named `ebsd` should
% already be available in the MATLAB workspace, with the BCC ferrite
% phase accessible as `ebsd('BCC')`.

%% ========================================================================
% Calculate the BCC orientation distribution function
% ==========================================================================

% The ODF is calculated from the measured BCC ferrite orientations.
odf_alpha = calcDensity( ...
    ebsd('BCC').orientations, ...
    'halfwidth', 5 * degree);


%% ========================================================================
% Define the Kurdjumov-Sachs orientation relationship
% ==========================================================================

KS = orientation.map( ...
    Miller(1,1,1,FCC), Miller(0,1,1,BCC), ...
    Miller(-1,0,1,FCC), Miller(-1,-1,1,BCC) );


%% ========================================================================
% User parameters
% ==========================================================================

% Standard deviation of the angular deviation from the ideal K-S
% orientation
KS_spread = 3.33 * degree;

% Number of random orientations sampled around each candidate alpha2
Ncloud = 100;

% Minimum misorientation between alpha1 and alpha2
minMisori = 10;

% Step size along the alpha-fiber
fiberStep = 1;

% Output folder
outFolder = 'results';


%% ========================================================================
% Define the alpha-fiber
% <110> || RD
% ==========================================================================

Phi_list = (0:fiberStep:90)';

ori_alpha1_list = orientation.byEuler( ...
    zeros(size(Phi_list)) * degree, ...
    Phi_list * degree, ...
    45 * degree, ...
    BCC);


%% ========================================================================
% Initialize the double K-S potential
% ==========================================================================

DoubleKS_Potential = ...
    zeros(length(ori_alpha1_list), 1);


%% ========================================================================
% Create output folder
% ==========================================================================

if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end


%% ========================================================================
% Calculate the double K-S potential
% ==========================================================================

fprintf('Calculating double K-S potential along the alpha-fiber...\n');

for k = 1:length(ori_alpha1_list)

    % Reference alpha1 orientation
    ori_alpha1 = ori_alpha1_list(k);


    % Possible gamma orientations related to alpha1
    ori_gamma = ...
        ori_alpha1 * inv(KS.parents);


    % Candidate alpha2 orientations
    ori_alpha2_all = ...
        variants(KS, ori_gamma);


    % Remove duplicate candidate alpha2 orientations
    [ori_alpha2_unique, ~, ~] = unique( ...
        ori_alpha2_all, ...
        'tolerance', 2 * degree, ...
        'stable');


    % Remove candidate alpha2 orientations that are too close
    % to the reference alpha1 orientation
    misori = ...
        angle(ori_alpha2_unique, ori_alpha1) / degree;

    valid = ...
        misori >= minMisori;

    ori_alpha2_unique = ...
        ori_alpha2_unique(valid);


    % Skip if no valid candidate alpha2 orientations remain
    if isempty(ori_alpha2_unique)

        DoubleKS_Potential(k) = 0;
        continue;

    end


    %% Evaluate the ODF around each candidate alpha2 orientation

    meanDensity = ...
        zeros(length(ori_alpha2_unique), 1);

    stdDensity = ...
        zeros(length(ori_alpha2_unique), 1);

    maxDensity = ...
        zeros(length(ori_alpha2_unique), 1);


    for i = 1:length(ori_alpha2_unique)

        densityCloud = ...
            zeros(Ncloud, 1);


        for n = 1:Ncloud

            % Random rotation axis
            ax = vector3d.rand(1);

            % Gaussian angular deviation from the ideal K-S orientation
            ang = randn * KS_spread;

            R = ...
                rotation.byAxisAngle(ax, ang);

            % Perturbed candidate alpha2 orientation
            alpha2_dev = ...
                R * ori_alpha2_unique(i);

            % Evaluate the experimental ODF
            densityCloud(n) = ...
                eval(odf_alpha, alpha2_dev);

        end


        meanDensity(i) = ...
            mean(densityCloud);

        stdDensity(i) = ...
            std(densityCloud);

        maxDensity(i) = ...
            max(densityCloud);

    end


    %% Calculate the double K-S potential

    DoubleKS_Potential(k) = ...
        sum(meanDensity);


    %% Save candidate alpha2 information

    T = table;

    T.Rank = ...
        (1:length(ori_alpha2_unique))';

    T.phi1_deg = ...
        round(ori_alpha2_unique.phi1 / degree, 2);

    T.Phi_deg = ...
        round(ori_alpha2_unique.Phi / degree, 2);

    T.phi2_deg = ...
        round(ori_alpha2_unique.phi2 / degree, 2);

    T.MeanDensity_mrd = ...
        round(meanDensity, 3);

    T.StdDensity_mrd = ...
        round(stdDensity, 3);

    T.MaxDensity_mrd = ...
        round(maxDensity, 3);

    T.Misori_to_alpha1_deg = ...
        round( ...
            angle(ori_alpha2_unique, ori_alpha1) / degree, ...
            2);


    %% Reference alpha1 orientation

    phi1_now = ...
        round(ori_alpha1.phi1 / degree, 1);

    Phi_now = ...
        round(ori_alpha1.Phi / degree, 1);

    phi2_now = ...
        round(ori_alpha1.phi2 / degree, 1);


    T.Alpha1_phi1_deg = ...
        repmat(phi1_now, height(T), 1);

    T.Alpha1_Phi_deg = ...
        repmat(Phi_now, height(T), 1);

    T.Alpha1_phi2_deg = ...
        repmat(phi2_now, height(T), 1);

    T.DoubleKS_Potential = ...
        repmat(DoubleKS_Potential(k), height(T), 1);


    %% Sort candidate alpha2 orientations

    T = sortrows( ...
        T, ...
        'MeanDensity_mrd', ...
        'descend');


    %% Save candidate alpha2 results

    fileName = sprintf( ...
        'alpha1_phi1_%04.1f_Phi_%04.1f_phi2_%04.1f.csv', ...
        phi1_now, ...
        Phi_now, ...
        phi2_now);

    writetable( ...
        T, ...
        fullfile(outFolder, fileName));

end


%% ========================================================================
% Save the double K-S potential
% ==========================================================================

Result = table;

Result.Phi_deg = ...
    Phi_list;

Result.DoubleKS_Potential = ...
    round(DoubleKS_Potential, 3);


summaryFile = fullfile( ...
    outFolder, ...
    'AlphaFiber_DoubleKS_Potential.csv');

writetable( ...
    Result, ...
    summaryFile);


%% ========================================================================
% Display maximum double K-S potential
% ==========================================================================

[MaxPotential, idxMax] = ...
    max(DoubleKS_Potential);

BestPhi = ...
    Phi_list(idxMax);


fprintf('\n');
fprintf('---------------------------------------------\n');
fprintf('Maximum double K-S potential = %.3f\n', ...
    MaxPotential);

fprintf('Reference alpha-fiber Phi = %.1f deg\n', ...
    BestPhi);

fprintf('---------------------------------------------\n');


%% ========================================================================
% Plot the double K-S potential
% ==========================================================================

figure;

plot( ...
    Phi_list, ...
    DoubleKS_Potential, ...
    '-o', ...
    'LineWidth', 2);

xlabel('\Phi (deg)');
ylabel('Double K-S Potential');

title('Double K-S Potential Along \alpha-Fiber');

grid on;
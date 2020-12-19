% Martha Gizaw
% Engineering Mathematics
% Final Project
% Determining MSEs Between Optimized NN-Based Chromosomes
% December 17, 2020

%% Producing Handwritten Number Image Samples
% Load two csv files to the MATLAB workspace.
tr = csvread('train.csv', 1, 0);
sub = csvread('test.csv', 1, 0);

% Produce grayscale subplots of individual numbers.
figure                                          % plot images
colormap(gray)                                  % set to grayscale
for i = 1:25                                    % preview first 25 samples
    subplot(5,5,i)                              % plot them in 6 x 6 grid
    digit = reshape(tr(i, 2:end), [28,28])';    % row = 28 x 28 image
    imagesc(digit)                              % show the image
    title(num2str(tr(i, 1)))                    % show the label
end

%% Storing Chromosome Strands and MSEs

% Initialize the following empty arrays.
chromsMSE = [];             % optimal chromosome MSEs

storeOptimumChrom = [];     % chromosome strands after 100 generations
storeOptimumChrom2 = [];    % chromosome strands after 1 minute

% Specify any number of trials you want.
numTrials = 5;

%% Start the loop!
for i = 1:numTrials
%% Section 1: Neural Network GUI
    n = size(tr, 1);                    % number of samples in the dataset
    targets  = tr(:,1);                 % 1st column is |label|
    targets(targets == 0) = 10;         % use '10' to present '0'
    targetsd = dummyvar(targets);       % convert label into a dummy variable
    inputs = tr(:,2:end);               % the rest of columns are predictors

    inputs = inputs';                   % transpose input
    targets = targets';                 % transpose target
    targetsd = targetsd';               % transpose dummy variable

    rng(1);                             % for reproducibility
    c = cvpartition(n,'Holdout',n/3);   % hold out 1/3 of the dataset

    Xtrain = inputs(:, training(c));    % 2/3 of the input for training
    Ytrain = targetsd(:, training(c));  % 2/3 of the target for training
    Xtest = inputs(:, test(c));         % 1/3 of the input for testing
    Ytest = targets(test(c));           % 1/3 of the target for testing
    Ytestd = targetsd(:, test(c));      % 1/3 of the dummy variable for testing

%     Ypred = myNNfun(Xtest);             % predicts probability for each label
%     Ypred(:, 1:5)                       % display the first 5 columns
%     [~, Ypred] = max(Ypred);            % find the indices of max probabilities
%     sum(Ytest == Ypred) / length(Ytest) % compare the predicted vs. actual

    sweep = [10,50:50:200];                 % parameter values to test
    scores = zeros(length(sweep), 1);       % pre-allocation
    models = cell(length(sweep), 1);        % pre-allocation
    x = Xtrain;                             % inputs
    t = Ytrain;                             % targets
    trainFcn = 'trainscg';                  % scaled conjugate gradient
    for i = 1:length(sweep)
        hiddenLayerSize = sweep(i);         % number of hidden layer neurons
        net = patternnet(hiddenLayerSize);  % pattern recognition network
        net.divideParam.trainRatio = 70/100;% 70% of data for training
        net.divideParam.valRatio = 15/100;  % 15% of data for validation
        net.divideParam.testRatio = 15/100; % 15% of data for testing
        net = train(net, x, t);             % train the network
        models{i} = net;                    % store the trained network
        p = net(Xtest);                     % predictions
        [~, p] = max(p);                    % predicted labels
        scores(i) = sum(Ytest == p) /...    % categorization accuracy
            length(Ytest);
    end

    figure
    plot(sweep, scores, '.-')
    xlabel('number of hidden neurons')
    ylabel('categorization accuracy')
    title('Number of hidden neurons vs. accuracy')

    %% Section 2: Generation-Based Chromosomes
    rng('shuffle')  %initializes the random number generator

    %specify whether to minimize or maximize the objective function
    %maximize: minmax = 1;
    %minimize: minmax = -1;
    minmax = -1;

    %constants of the genetic algorithm
    numChroms = 60;
    numGenerations = 100;
    reproduced = 0.25;  %percentage of the elite population that is reproduced
    operator1 = 0.25;   %percentage of the new population produced by the first operator
    BLXa = 0.35;        %alpha of the BLXa cross-over operator
    operator2 = 0.25;   %percentage of the new population produced by the second operator
    numB = 1;           %parameter (b) of the non-uniform mutation operator
    randGen = 0.25;     %percentage of the new population produced randomly
    %(reproduced + operator1 + operator2 + randGen = 1)

    %identify the model
    modelScript = 'RungeKuttaExample3ga.m';

    %identify any additional model data or parameters needed ----------
    Yxs = 0.65;

    %--------------------

    %check reproduction/generation percentages
    if reproduced + operator1 + operator2 + randGen ~= 1
        fprintf(['Error. The percentages from reproduction, operator1, operator2, and random generation must add up to 1.' '\n']);
        return
    end

    %input data and run an initial model calculation
    gaInput = transpose(scores);    %for a 2-gene chromosome [umax Ks] based on llsa
    limits(1,:) = 0.75.*gaInput;
    limits(2,:) = 1.25.*gaInput;
    limits(3,:) = limits(2,:) - limits(1,:);

    %initiate chromosomes
    chromLength = size(gaInput,2);
    chrom = zeros(numChroms, chromLength);
    chrom(1,:) = gaInput;

    for i = 2:numChroms

        for j = 1:chromLength

            chrom(i,j) = limits(1,j)+rand.*(limits(3,j));

        end

    end

    %evaluate all of the inital chromosomes (first generation)
    generation = 1;
    chromResult = zeros(numChroms,1);

    for m = 1:size(chrom,1)

        %run the model script, evaluate the objective function
        run(modelScript);

        chromResult(m,1) = result;

    end

    %the evaluated chromosomes are now "old" chromosomes, combine with results
    chromOld = [chrom chromResult];

    if minmax == 1

        chromOld = sortrows(chromOld,-(chromLength+1));

    else

        chromOld = sortrows(chromOld,(chromLength+1));

    end

    %start the Genetic Algorithm loop
    while generation <= numGenerations

        %initiate the new chromosome variable
        chrom = zeros(numChroms, chromLength);

        %reproduction
        %elite chromosomes
        for n = 1:ceil(numChroms.*reproduced)
            chrom(n,:) = chromOld(n,1:chromLength);
        end
        n = n + 1;

        %operator1: BLXa cross-over (a=0.35)
        %strategy: (i) choose two of the reproduced chromosomes, (ii) perform
        %operation, (iii) add them both to the chrom list

        while n < numChroms.*reproduced + numChroms.*operator1

            %choose the reproduced chromosomes and save them as "temporary"
            %chromosomes
            BLXtemp1 = chrom(ceil(rand.*numChroms.*reproduced),:);
            BLXtemp2 = chrom(ceil(rand.*numChroms.*reproduced),:);

            %choose the gene (column) to manipulate
            BLXcol = ceil(rand.*chromLength);

            %BLXa operators
            cmax = max(BLXtemp1(1,BLXcol),BLXtemp2(1,BLXcol));
            cmin = min(BLXtemp1(1,BLXcol),BLXtemp2(1,BLXcol));

            Ic = cmax - cmin;

            if Ic == 0
                cmax = cmax + (rand.*0.5.*cmax);
                cmin = cmin - (rand.*0.5.*cmin);
                Ic = cmax - cmin;
            end

            BLXmin = cmin - Ic.*BLXa;
            BLXmax = cmax + Ic.*BLXa;

            %perform the operation
            BLXtemp1(1,BLXcol) = BLXmin + rand.*(BLXmax - BLXmin);
            BLXtemp2(1,BLXcol) = BLXmin + rand.*(BLXmax - BLXmin);

            %check that these values are within defined limits, if not reset to
            %limit value

            if BLXtemp1(1,BLXcol) < limits(1,BLXcol)
                BLXtemp1(1,BLXcol) = limits(1,BLXcol);
            end
            if BLXtemp1(1,BLXcol) > limits(2,BLXcol)
                BLXtemp1(1,BLXcol) = limits(2,BLXcol);
            end  
            if BLXtemp2(1,BLXcol) < limits(1,BLXcol)
                BLXtemp2(1,BLXcol) = limits(1,BLXcol);
            end
            if BLXtemp2(1,BLXcol) > limits(2,BLXcol)
                BLXtemp2(1,BLXcol) = limits(2,BLXcol);
            end 

            %insert into the chrom list
            chrom(n,:) = BLXtemp1;
            n = n + 1;
            chrom(n,:) = BLXtemp2;
            n = n + 1;
        end

        %operator2: Non-uniform mutation
        while n < numChroms.*reproduced + numChroms.*operator1 + numChroms.*operator2

            %choose a reproduced chromosome and save it as a temporary
            %chromosome
            numTemp1 = chrom(ceil(rand.*numChroms.*reproduced),:);

            %choose the gene (column) to manipulate
            numCol = ceil(rand.*chromLength);

            %non-uniform mutation operators and stochastic variables
            tau = round(rand);

            numYup = limits(2,numCol) - numTemp1(1,numCol);
            numYdown = numTemp1(1,numCol) - limits(1,numCol);

            numTup = numYup.*(1 - rand.^((1-generation./numGenerations).^numB));
            numTdown = numYdown.*(1-rand.^((1-generation./numGenerations).^numB));

            if tau == 0
                numTemp1(1,numCol) = numTemp1(1,numCol) - numTdown;
            else
                numTemp1(1,numCol) = numTemp1(1,numCol) + numTup;
            end

            %check limits
            if numTemp1(1,numCol) < limits(1,numCol)
                numTemp1(1,numCol) = limits(1,numCol);
            end
            if numTemp1(1,numCol) > limits(2,numCol)
                numTemp1(1,numCol) = limits(2,numCol);
            end

            %insert into chrom list
            chrom(n,:) = numTemp1;
            n = n + 1;
        end

        %fill-in the rest of the available chrom spaces with randomly-generated
        %chromosomes

        while n <= numChroms

            %initiate the chromosome
            chrom(n,:) = zeros(1,chromLength);

            %generate random genes around limits
            for j = 1:chromLength

                chrom(n,j) = limits(1,j) + rand.*limits(3,j);

            end
            n = n + 1;

        end

        %evaluate all chromosomes (this generation)
        chromResult = zeros(numChroms,1);

        for m = 1:size(chrom,1)

            %run the model script, evaluate the objective function
            run(modelScript);

            chromResult(m,1) = result;

        end

        %the evaluated chromosomes are now "old" chromosomes with results
        chromOld = [];
        chromOld = [chrom chromResult];

        %sort results
        if minmax == 1
            chromOld = sortrows(chromOld,-(chromLength + 1));
        else
            chromOld = sortrows(chromOld,(chromLength + 1));
        end

        %start the new generation
        fprintf(['Computing generation ' num2str(generation) ' of ' num2str(numGenerations) ...
            '. Optimum result: ' num2str(chromOld(1,chromLength + 1)) '\n']);

        generation = generation + 1;

    end

    optimumChrom = chromOld(1,1:chromLength);
    optimumResult = chromOld(1,chromLength+1);

    fprintf(['Optimum chromosome: ' num2str(optimumChrom) '\n']);
    fprintf(['Optimum result: ' num2str(optimumResult) '\n']);
    fprintf(['Done!' '\n']);

    %% Section 3: Time-Based Chromosomes
    rng('shuffle')  %initializes the random number generator
    tic

    %specify whether to minimize or maximize the objective function
    %maximize: minmax = 1;
    %minimize: minmax = -1;
    minmax = -1;

    %constants of the genetic algorithm
    totalMinutes = 1;   %the total run time (minutes) for the genetic algorithm
    numChroms = 60;
    numGenerations = 100;
    reproduced = 0.25;  %percentage of the elite population that is reproduced
    operator1 = 0.25;   %percentage of the new population produced by the first operator
    BLXa = 0.35;        %alpha of the BLXa cross-over operator
    operator2 = 0.25;   %percentage of the new population produced by the second operator
    numB = 1;           %parameter (b) of the non-uniform mutation operator
    randGen = 0.25;     %percentage of the new population produced randomly
    %(reproduced + operator1 + operator2 + randGen = 1)

    %identify the model
    modelScript = 'RungeKuttaExample3ga.m';

    %identify any additional model data or parameters needed ----------
    Yxs = 0.65;

    %--------------------

    %check reproduction/generation percentages
    if reproduced + operator1 + operator2 + randGen ~= 1
        fprintf(['Error. The percentages from reproduction, operator1, operator2, and random generation must add up to 1.' '\n']);
        return
    end

    %input data and run an initial model calculation
    gaInput = transpose(scores);    %for a 2-gene chromosome [umax Ks] based on llsa
    limits(1,:) = 0.75.*gaInput;
    limits(2,:) = 1.25.*gaInput;
    limits(3,:) = limits(2,:) - limits(1,:);

    %initiate chromosomes
    chromLength = size(gaInput,2);
    chrom = zeros(numChroms, chromLength);
    chrom(1,:) = gaInput;

    for i = 2:numChroms

        for j = 1:chromLength

            chrom(i,j) = limits(1,j)+rand.*(limits(3,j));

        end

    end

    %evaluate all of the inital chromosomes (first generation)
    generation = 1;
    chromResult = zeros(numChroms,1);

    for m = 1:size(chrom,1)

        %run the model script, evaluate the objective function
        run(modelScript);
        chromResult(m,1) = result;
        runTimeMinutes = toc./60;   %run time in minutes

    end

    %the evaluated chromosomes are now "old" chromosomes, combine with results
    chromOld = [chrom chromResult];

    if minmax == 1

        chromOld = sortrows(chromOld,-(chromLength+1));

    else

        chromOld = sortrows(chromOld,(chromLength+1));

    end

    while runTimeMinutes < totalMinutes

        %initiate the new chromosome variable
        chrom = zeros(numChroms, chromLength);

        %reproduction
        %elite chromosomes
        for n = 1:ceil(numChroms.*reproduced)
            chrom(n,:) = chromOld(n,1:chromLength);
        end
        n = n + 1;

        %operator1: BLXa cross-over (a=0.35)
        %strategy: (i) choose two of the reproduced chromosomes, (ii) perform
        %operation, (iii) add them both to the chrom list

        while n < numChroms.*reproduced + numChroms.*operator1

            %choose the reproduced chromosomes and save them as "temporary"
            %chromosomes
            BLXtemp1 = chrom(ceil(rand.*numChroms.*reproduced),:);
            BLXtemp2 = chrom(ceil(rand.*numChroms.*reproduced),:);

            %choose the gene (column) to manipulate
            BLXcol = ceil(rand.*chromLength);

            %BLXa operators
            cmax = max(BLXtemp1(1,BLXcol),BLXtemp2(1,BLXcol));
            cmin = min(BLXtemp1(1,BLXcol),BLXtemp2(1,BLXcol));

            Ic = cmax - cmin;

            if Ic == 0
                cmax = cmax + (rand.*0.5.*cmax);
                cmin = cmin - (rand.*0.5.*cmin);
                Ic = cmax - cmin;
            end

            BLXmin = cmin - Ic.*BLXa;
            BLXmax = cmax + Ic.*BLXa;

            %perform the operation
            BLXtemp1(1,BLXcol) = BLXmin + rand.*(BLXmax - BLXmin);
            BLXtemp2(1,BLXcol) = BLXmin + rand.*(BLXmax - BLXmin);

            %check that these values are within defined limits, if not reset to
            %limit value

            if BLXtemp1(1,BLXcol) < limits(1,BLXcol)
                BLXtemp1(1,BLXcol) = limits(1,BLXcol);
            end
            if BLXtemp1(1,BLXcol) > limits(2,BLXcol)
                BLXtemp1(1,BLXcol) = limits(2,BLXcol);
            end  
            if BLXtemp2(1,BLXcol) < limits(1,BLXcol)
                BLXtemp2(1,BLXcol) = limits(1,BLXcol);
            end
            if BLXtemp2(1,BLXcol) > limits(2,BLXcol)
                BLXtemp2(1,BLXcol) = limits(2,BLXcol);
            end 

            %insert into the chrom list
            chrom(n,:) = BLXtemp1;
            n = n + 1;
            chrom(n,:) = BLXtemp2;
            n = n + 1;
        end

        %operator2: Non-uniform mutation
        while n < numChroms.*reproduced + numChroms.*operator1 + numChroms.*operator2

            %choose a reproduced chromosome and save it as a temporary
            %chromosome
            numTemp1 = chrom(ceil(rand.*numChroms.*reproduced),:);

            %choose the gene (column) to manipulate
            numCol = ceil(rand.*chromLength);

            %non-uniform mutation operators and stochastic variables
            tau = round(rand);

            numYup = limits(2,numCol) - numTemp1(1,numCol);
            numYdown = numTemp1(1,numCol) - limits(1,numCol);

            numTup = numYup.*(1 - rand.^((1-generation./numGenerations).^numB));
            numTdown = numYdown.*(1-rand.^((1-generation./numGenerations).^numB));

            if tau == 0
                numTemp1(1,numCol) = numTemp1(1,numCol) - numTdown;
            else
                numTemp1(1,numCol) = numTemp1(1,numCol) + numTup;
            end

            %check limits
            if numTemp1(1,numCol) < limits(1,numCol)
                numTemp1(1,numCol) = limits(1,numCol);
            end
            if numTemp1(1,numCol) > limits(2,numCol)
                numTemp1(1,numCol) = limits(2,numCol);
            end

            %insert into chrom list
            chrom(n,:) = numTemp1;
            n = n + 1;
        end

        %fill-in the rest of the available chrom spaces with randomly-generated
        %chromosomes

        while n <= numChroms

            %initiate the chromosome
            chrom(n,:) = zeros(1,chromLength);

            %generate random genes around limits
            for j = 1:chromLength

                chrom(n,j) = limits(1,j) + rand.*limits(3,j);

            end
            n = n + 1;

        end

        %evaluate all chromosomes (this generation)
        chromResult = zeros(numChroms,1);

        for m = 1:size(chrom,1)

            %run the model script, evaluate the objective function
            run(modelScript);
            chromResult(m,1) = result;
            runTimeMinutes = toc./60;   %run time in minutes

        end

        %the evaluated chromosomes are now "old" chromosomes with results
        chromOld = [];
        chromOld = [chrom chromResult];

        %sort results
        if minmax == 1
            chromOld = sortrows(chromOld,-(chromLength + 1));
        else
            chromOld = sortrows(chromOld,(chromLength + 1));
        end

        %start the new generation
        fprintf(['Computing for ' num2str(runTimeMinutes) ' of ' num2str(totalMinutes) ...
            ' minutes. Optimum result: ' num2str(chromOld(1,chromLength + 1)) '\n']);

        generation = generation + 1;

    end

    optimumChrom2 = chromOld(1,1:chromLength);
    optimumResult2 = chromOld(1,chromLength+1);

    fprintf(['Optimum chromosome: ' num2str(optimumChrom2) '\n']);
    fprintf(['Optimum result: ' num2str(optimumResult2) '\n']);
    fprintf(['Done!' '\n']);

%% Find the MSE values between the generation and time optimized chromosomes.
    between_chroms = immse(optimumChrom, optimumChrom2);
 
    storeOptimumChrom = [storeOptimumChrom, optimumChrom];
    storeOptimumChrom2 = [storeOptimumChrom2, optimumChrom2];
    chromsMSE = [chromsMSE, between_chroms];
end
%% Plot the MSEs values per trial.
plot([1:numTrials], chromsMSE, 'o-')
xlabel('Trial')
ylabel('Chromosome MSE')
title('MSEs Between Generation and Time Chromosomes')

%% Sources

% [1] “Digit Recognizer.” (2012). Kaggle.
% https://www.kaggle.com/c/digit-recognizer/data (accessed December 17,
% 2020).

% [2] L. Shure. “Artificial Neural Networks for Beginners.” (August
% 4, 2015). MathWorks.
% https://blogs.mathworks.com/loren/2015/08/04/artificial-neural-networks-for-beginners/
% (accessed December 17, 2020).

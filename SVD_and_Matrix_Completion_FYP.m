%% Initiation

clc; clf; clear; rng(2);


%% Data extraction and structuring
%might change data & structure
Data_extraction_and_structuring; %calls script which extracts and structures the data in a meaningfull way

A_original =[0, 0, 5, 4; ...
    0, 1, 4, 0; ...
    4, 5, 2, 0; ...
    0, 4, 2, 1;...
    4, 0, 1, 2; ...
    1, 2, 0, 5];

%% Intial guess
[rows, columns] = size(A_original);
averages = size(A_original);
for i = 1:rows
    for j = 1:columns
      n = numel(nonzeros([A_original(:,j); A_original(i,:)'])) - 1 * (A_original(i,j)~=0); %number of non-zero entries in a row and column        
      averages(i, j) = (sum([A_original(:,j); A_original(i,:)']) - A_original(i, j)) / n; %average of the non-zero rows and columns   
    end
end


A = A_original; %copies original data

for i = 1:rows
    for j = 1:columns
        if A_original(i,j) == 0
            A(i,j) = averages(i,j); %replaces zero values with mean of each that row and column
        end
    end
end

%% SVD

Ac = A - averages; %centered matrix of A

[U, S, V] = svd(Ac); %computes the U matrix, eigenvalues S and eigen vectors V


%% Rank
eigen_val = diag(S); %eigen values in vector form
rank = numel(eigen_val); %number of eigen values

variance = cumsum(eigen_val)/sum(eigen_val); % gives an indication of the varies explained by the lower dimensional matrix

k = find( ( variance > 0.8 ), 1); % chooses the lower rank k, for which 80% of the variance is explained

Ac_recon = U(:, 1:k)*S(1:k, 1:k)*V(:, 1:k)'; %lower k-rank recontructed centered matrix of A

%% Reverse centering operation
A_new = Ac_recon + averages; %new matrix from SVD

%% itterations
sigma = 1; %stopping criterion
itterations = 0; %number of itterations
elements = 0; %counter used to track convergence of values

while sigma > 1e-10
    elements = elements + 1;
    counter = 0; %counter used to track the convergence of values
    for i = 1:rows
        for j = 1:columns
            if A_original(i,j) == 0
                A(i,j) = A_new(i,j); %replaces the unknown values by the itterated values
                counter = counter + 1;
                track_itterations(counter, elements) = A_new(i,j);
            end
        end
    end
    
    %Convergence
    itterations = itterations + 1;
    
    SSl(itterations) = sum(( A(A_original == 0) ).^2); % measure of convergence
    
    
    %repeat centering
    
    for i = 1:rows
        for j = 1:columns
          n = numel([A(:,j); A(i,:)']) - 1; %number of non-zero entries in a row and column        
          averages(i, j) = (sum([[A(:,j); A(i,:)']]) - A(i, j)) / n; %average of the non-zero rows and columns   
        end
    end
    
    
    Ac = A - averages;
    
    %repeat SVD
    
    [U, S, V] = svd(Ac);
    
    %repeat lower rank matrix completion
    
    Ac_recon = U(:, 1:k)*S(1:k, 1:k)*V(:, 1:k)'; %lower k-rank recontructed centered matrix of A
    
    %repeat reverse centering
    A_new = Ac_recon + averages; %new matrix from SVD
    
    %stopping criterion
    if itterations == 1
        sigma = 1;
    else
        sigma = abs( SSl(itterations) - SSl(itterations-1)) / SSl(itterations);
    end
    
    User_input = []; %used for maximum itterations
    
    if rem(itterations,5000) == 0 %some fancy stuff for insuring the loop does not continue indefinitely
        % asks user after every 5000 itterations whether to stop or
        % continue itterations
        
        User_input = 'o';
        
        while User_input == 'o'
            
            fprintf('Maximum number itterations exceeded. Program stopped pre-maturely. Results might be incorrect\n\n');
            User_input = input('Continue itterations? Y-Yes, N-No: \n\n', "s");   
        
            if User_input == 'N' | User_input == 'n'
                break;
                
            elseif User_input == 'Y' | User_input == 'y'
            
            else
                fprintf('Invalid selection. Please enter valid selection\n\n');
                User_input = 'o'; 
            end
        end  
       
    end %end of fancy stuff
    
    if User_input == 'N' | User_input == 'n' %exist the itteration loop
        break
        
    elseif isempty(User_input)
        
    end
end %end of itterations

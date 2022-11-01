%%Initation
clc; clf; rng(1); %clear; 


%%Data extraction
%Data_extraction_and_structuring;
functional_groups = table2struct(readtable("GR Hermanus (22796002).xlsx", "Sheet", "Properties"), 'ToScalar', true);
functional_groups.Name = string(functional_groups.Name);
functional_groups.FunctionalGroup = string(functional_groups.FunctionalGroup);

index_smaller_subset = (functional_groups.FunctionalGroup == "Alkane" | ...
     functional_groups.FunctionalGroup == "Primary Alcohol" | ...
     functional_groups.FunctionalGroup == "Secondary Alcohol" |...
     functional_groups.FunctionalGroup == "Diol" | ...
     functional_groups.FunctionalGroup == "Branched Alcohol" | ...
     functional_groups.FunctionalGroup == "Ether" | ... 
     functional_groups.FunctionalGroup == "Ester"); %index for the smaller subset of the data
index_smaller_subset = ones(size(index_smaller_subset)) == 1;


%%LOOCV and optimal rank
A_original2 = tensor_4D(index_smaller_subset, index_smaller_subset, 1, 1);
[number_rows, number_columns] = size(A_original2);
size_A = min([number_rows, number_columns]);
options = ["c", "r", "all", "rc", "JR", "ideal"];

%below code runs the LOOCV on the data set ans gives the MSE and MAE of the
%data with size of l j k x; where l = composition, j = rank, k = option for
%centering, x = option for the initial guess 
ii = index;
index = 0;

for x = 5:6
    option_initial_guess = options(x);
    

    for k = 1:2
         option_centering = options(k);
    
        
        for l = 1: number_composition    
            
            A_original2 = tensor_4D(index_smaller_subset, index_smaller_subset,l,:);

            for j = 3 : 8 
                rank = j;
                temp = A_original2.*( eye(size_A, size_A) == 0);
                known_indices = find(A_original2.*(eye(size_A,size_A)==0) ~= 0);
                unknown_indices = (A_original2 == 0);

                counter = 0;

                for i = 1: numel(known_indices)

                    if index < ii
                        index = index + 1;
                    else


                        temp = unknown_indices; %temporary variable
                        temp(known_indices(i)) = 1; %sets element i equals to one, i.e. removing that spesific element of the original matrix

                        if l == 5

                            temp = (temp ~= unknown_indices)' + temp;

                        end

                        A_original = A_original2.*(temp == 0); %extracts new non-zero elements of the original matrix

                        if ( (sum(A_original, 1) ~= 0) & (sum(A_original, 2) ~= 0) ) %checks wether the rows and columns has all non-zero entries
                            counter = counter + 1;
                            %if both the rows and columns does not have all zero entries repeat
                            %itterations
                            %% Intial guess

                            A = fun_intial_guess(A_original, option_initial_guess, mole_fraction_Component1(l)); %intial guess

                            %% itterations
                            sigma = 10; %stopping criterion
                            itterations = 0; %number of itterations

                            while sigma > 1

                                A_new = fun_svd_itterations(A, rank, option_centering, mole_fraction_Component1(l)); %centering, SVD, and reverse centering
                                A = A_original + A_new.*temp; %replaces the unknown values by the itterated values

                                %Convergence
                                itterations = itterations + 1;

                                %stopping criterion
                                if itterations == 1
                                    sigma = 10;
                                else
                                    sigma = sum (abs( ( A_new(temp == 1) - A_prev(temp == 1) ) ));
                                end %end of stopping criterion

                                User_input = []; %used for maximum itterations

                                if rem(itterations,1000000) == 0 %some fancy stuff for insuring the loop does not continue indefinitely
                                    % asks user after every 5000 itterations whether to stop or
                                    % continue itterations

                                    User_input = 'o'; %sets user_input to 'o' used to control while loop

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

                                %additional code to exit while loop if user chooses
                                if User_input == 'N' | User_input == 'n' %exist the itteration loop
                                    break

                                elseif isempty(User_input)

                                end %end of additional code
                                %idx = setdiff(1:numel(known_indices), counter);
                                %subplot(3, 1, 1);
                                %plot(itterations, sum( (A_new(known_indices(idx)) - A_original(known_indices(idx))).^2), '.k'); hold on
                                %subplot(3, 1, 2);
                                %plot(itterations, sigma, '.r'); hold on;
                                %subplot(3, 1, 3);
                                %plot(itterations, abs(A_original2(known_indices(counter))-A(known_indices(counter))), '.b'); hold on
                                %pause(0.1);
                                if sigma > 0.5  
                                    A_prev = A_new;
                                end

                            end %end of itterations
                            if l == 5 %for composition = 0.5 matrix should be symmetrical
                                A = (A+A')/2; %makes the matrix symmetrical by taking average;
                            end
                            index = index + 1;
                            LOOCV.Original (index, :) = A_original2(known_indices(i));
                            LOOCV.Itterations (index, :) = A(known_indices(i));
                            LOOCV.Initial (index, :) = options(x);
                            LOOCV.Centering (index, :) = options(k);
                            LOOCV.Composition (index, :) = l;
                            LOOCV.Rank (index, :) = rank;
                            LOOCV.NumberOfItterations (index, :) = itterations;
                            clc;
                            
                            if rem(index,100) == 0
                                save("Results_LOOCV_subset_All");
                            end
                        end

                    end

                end


            end
        
        end
    end
end

function A_guess = fun_intial_guess (A, option, composition)
%this function calculates the intial guesses based on the option type, row
%averages, columns averages and a combination of row and column averages
option = convertCharsToStrings(option); %makes sure that a string is used
A_guess = A;
%[unknown_rows, unknwon_columns] = find(A_guess == 0); %returns the indices of the unknown rows and columns 
check = 0; %intilises to zero

[number_rows, number_columns] = size(A);

while check == 0
    
    if ( (option == "r") | (option == "c") | (option == "rc") | (option == "all") | (option == "JR") | (option == "ideal"))
        
       check = 1;
       
    elseif check == 0
        
        option = input('Please enter a valid selection for "option".\nThe options as "r", "c" or "rc"\n', "s");
       
    end
    
    if option == "r"
        
        for i = 1:number_rows
            
            for j = 1:number_columns
                
                n = numel( nonzeros(A(i, :)) ); %number of non-zerp elements in the row
            
                if n == 0
                    fprintf('Warning! No entries in the row. No calculations were done'); %protective coding
                    break; % goes to next row
                
                else
                    average = sum(A(i, :), 2)/n;
                    
                    if A(i, j) == 0
                        
                        A_guess(i, j) = average; %assignes zero values to row average
                        
                    end
                    
                end
                
            end % end of assigning one average to row i and column j
            
        end % end of assigning average to row
        
        
    %end of rows averages
    
    elseif option == "c"
        
        for i = 1:number_rows
            
            for j = 1:number_columns
                
                n = numel( nonzeros(A(:, j)) ); %number of non-zero elements in the column
            
                if n == 0
                    fprintf('Warning! No entries in the column. No calculations were done'); %protective coding
                    break; % goes to next column
                
                else
                    average = sum(A(:, j), 1)/n;
                    
                    if A(i, j) == 0
                        
                        A_guess(i, j) = average; %assignes zero values to column average
                        
                    end
                    
                end
                
            end % end of assigning one average to row i and column j
            
        end % end of assigning average to a column
        
        
     %end of column averages
    
    elseif option == "rc"
        
        for i = 1:number_rows
            
            for j = 1:number_columns
                
              n = numel(nonzeros([A(:,j); A(i,:)'])) - 1 * (A(i,j)~=0); %number of non-zero entries in a row and column
              
              if n == 0
                    fprintf('Warning! No entries in the row and column. No calculations were done'); %protective coding
                    break; % goes to next entry
                
              else
              
                  average = (sum([A(:,j); A(i,:)']) - A(i, j)) / n; %average of the non-zero rows and columns  
                  
                  if A(i, j) == 0
                      
                      A_guess(i, j) = average;
                      
                  end
                  
                  
              end
              
              
              
            end
        end
        
     elseif option == "all"
    
            average = mean(nonzeros(A), "all");
            unknown = (A == 0);
            A_guess = A_guess + average*unknown;
            
        
    elseif option == "ideal"
        pure_comp = diag(A);
        
        for i = 1 : number_rows
            for j = 1: number_columns
                
                if A(i, j) == 0
                    if (pure_comp(i) == 0) | (pure_comp(j) == 0)
                        
                        A_guess(i, j) = mean(nonzeros(A(:, j))); %if pure component not known replace with column mean
                       
                        
                    else
                        
                        A_guess(i, j) = composition*pure_comp(i) + (1-composition)*pure_comp(j);
                        
                    end
                end
                
            end
        end
        
    elseif option == "JR"
        T = readtable("GR Hermanus (22796002).xlsx", "Sheet", "Densities and molecular weight");
        A_guess = A;
        pure_comp = diag(A);

        for i = 1:number_rows
            for j = 1:number_columns
                
                if A(i,j) == 0
                
                    if pure_comp(i) == 0 | pure_comp(j) == 0
                        
                        A_guess(i,j) = mean( nonzeros(A(i,:)), "all" );
                        
                    else 
                        
                        A_guess(i,j) = ( composition*T.MolecularWeight_kg_mol_(i)...
                            /T.Density_kg_m3_(i)+(1-composition)*T.MolecularWeight_kg_mol_(j)/T.Density_kg_m3_(j) )/...
                            ( (composition*T.MolecularWeight_kg_mol_(i))+((1-composition)*T.MolecularWeight_kg_mol_(j)) )^(1/2) ...
                            * ( ( composition*T.MolecularWeight_kg_mol_(i)/T.Density_kg_m3_(i)+(1-composition)*T.MolecularWeight_kg_mol_(j)/...
                            T.Density_kg_m3_(j) )/( composition*T.Density_kg_m3_(i)*pure_comp(i)^2+...
                            (1-composition)*T.Density_kg_m3_(j)*pure_comp(j)^2) )^(-1/2);
                        
                    end
                
                end
            end
        end
    end
end %end of while loop
    
    
    
end % end of intial guess function


function [Ac, average] = fun_centering(A, option)
%this function calculates the centered matrix based on 3 different methods.
%Centering based on columns means - r, centering based on row means - c,
%centering based on the rows and averages - rc, and centering based on all
%the entries - all
%the entries

option = convertCharsToStrings(option); %converts chars to strings
check = 0;

[number_rows, number_columns] = size(A);

while check == 0
    
     if ( (option == "r") | (option == "c") | (option == "rc") | (option =="all") )
        
       check = 1;
       
    elseif check == 0
        
        option = input('Please enter a valid selection for "option".\nThe options as "r", "c",  "rc" or "all"\n', "s");
       
     end
     
     if option == "r"
         
         for i = 1:number_rows
             
             if sum(A(i,:)) == 0
                 average(i, :) = zeros(1, number_columns);
                 
             else
                 average(i, :) = mean(A(i, :))*ones(1, number_columns);
             end
             
         end
         Ac = A - average;
         
     elseif option == "c"
         
         for j = 1:number_columns
             
             if sum(A(:,j)) == 0
                 average(:, j) = zeros(number_rows, 1);
                 
             else
                 average(:, j) = mean(A(:, j))*ones(number_rows, 1);
             end
             
         end
         Ac = A - average;
         
     elseif option == "rc"
         
         for i = 1:number_rows
            
            for j = 1:number_columns
                
              n = numel(nonzeros([A(:,j); A(i,:)'])) - 1 * (A(i,j)~=0); %number of non-zero entries in a row and column
              
              if n == 0
                   
                  average(i,j) = 0;
                
              else
              
                  average(i,j) = (sum([A(:,j); A(i,:)']) - A(i, j)) / n; %average of the non-zero rows and columns
                   
              end
              
            end
            
         end
         
         Ac = A - average;
         
     elseif option == "all"
         
         average = mean(A, "all")*ones(size(A));
         Ac = A - average;
         
         
     end
end



end %end of the centering function 


function A_new = fun_svd_itterations(A, rank, option, composition)
%this function calculates the svd itterations 
    [Ac, average] = fun_centering(A, option); %centering
    [U, S, V] = svd(Ac); %svd on centerd date
    Ac_recon = U(:,1:rank)*S(1:rank, 1:rank)*V(:,1:rank)'; %reconstructed centered data
    A_new = Ac_recon + average; %new martix
    
    if composition == 0.5
        A_new = (A_new + A_new')/2;
    end
    
end

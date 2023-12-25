function MC_Main()
%% load data
addpath('Data');
addpath('APG');
HI.Matrix=readmatrix('Data\HI_1999_2022_2.xlsx');
HI.Num_Sera=173;
HI.Virus_Year=readmatrix('Data\virus_year_2.xlsx');
HI.Sera_Year=readmatrix("Data\sera_year_2.xlsx");
save HI HI;
%% Process HI matrix 
Matrix = HI.Matrix;
Matrix = abs(Matrix);
known = (Matrix ~= 0);
unknown = (Matrix == 0);
Matrix(known) = log2(Matrix(known));
Max_ij = max(max(Matrix));
Max_j = max(Matrix);
for j = 1:HI.Num_Sera
    Matrix(:, j) = Max_ij - Max_j(j) + Matrix(:, j);
end
Matrix(unknown) = 0;
%% MC for the entire HI matrix
lambda = 0.2;
Full_MC_Mat = APG_outer(Matrix, known, lambda);
Full_MC_Mat(known) = Matrix(known);
%% window size
Win = 12;%[12 16 20 24];
Start_Year = min(HI.Virus_Year);
End_Year = max(HI.Virus_Year);
Count = zeros(size(Matrix));
Final_Mat = zeros(size(Matrix));
for i = 1:length(Win)
    W = Win(i);
    j = 0;
    while true
        A = Start_Year + j;
        B = A + W - 1;
        if B > End_Year
            break;
        end
        fprintf('From %d to %d.\n', A, B);
        [Sub_Mat, row_idx, col_idx] = Get_submat(Matrix, A, B, HI);
        Sub_known = known(row_idx, col_idx);
        Sub_Mat_MC = APG_outer(Sub_Mat, Sub_known, lambda);
        Final_Mat(row_idx, col_idx) = Final_Mat(row_idx, col_idx) + Sub_Mat_MC;
        Count(row_idx, col_idx) = Count(row_idx, col_idx) + 1;
        j = j+1;
    end
    Final_Mat( Final_Mat ~= 0 ) = Final_Mat( Final_Mat ~= 0 ) ./ Count( Final_Mat ~= 0 );
    Missed_idx = ( Final_Mat == 0 );
    Final_Mat(Missed_idx) = Full_MC_Mat(Missed_idx);
    Final_Mat(known) = Matrix(known);
    %% Matlab Kernal MDS test
    %MDS_test(Matrix, HI, 2);
    %MDS_test(Full_MC_Mat, HI, 2);
    %MDS_test(Final_Mat, HI, 2);
end
HI_MC = HI;
HI_MC.Matrix = Final_Mat;
save('HI_MC', 'HI_MC');
xlswrite('HI_MC_1999_2022_2.xlsx',HI_MC.Matrix);
[m, n] = size(Final_Mat);

%Y = zeros( m*(m-1)/2, 1 );
%r=1;
%for i= 1:m
    %for j=i+1:m
         %Y(r) = norm( Final_Mat(i,:) - Final_Mat(j,:), 1 ) / n;
         %r=r+1;
    %end
%end

%xlswrite('HI_distance_new.xlsx',Y);

function coor = MDS_test(Matrix, HI, dim)
%% cityblock = Manhattan Distance
D = pdist(Matrix, 'cityblock');
D = D / HI.Num_Sera;
%%
coor = mdscale(D, dim, 'criterion','metricsstress');
save HI_coor coor;
%xlswrite('HI_coor.xlsx',coor);

%coor = cmdscale(D, dim);
%size(coor)
Interval = 4;
Start_Year = min(HI.Virus_Year);
End_Year = max(HI.Virus_Year);
Set_Num = fix( (End_Year-Start_Year+1)/Interval ) + 1;
Idx_Set = cell(1, Set_Num);
figure;
Shape = 'o*d';
for i = 1:Set_Num
    Year_A = Start_Year + (i - 1)*Interval;
    Year_B = Year_A + Interval - 1;
    Idx_Set{i} = intersect( find(HI.Virus_Year >= Year_A), find(HI.Virus_Year <= Year_B) );
    if dim == 2
        plot(coor(Idx_Set{i},1), coor(Idx_Set{i},2), ...
            Shape( mod(i - 1, length(Shape))+1 ) );
        hold on;
    elseif dim == 3
        scatter3(coor(Idx_Set{i},1), coor(Idx_Set{i},2), coor(Idx_Set{i},3), ...
            Shape( mod(i - 1, length(Shape))+1 ) );
        hold on;
    else
        error('Unknow dim!\n');
    end
end

function [X, row_idx, col_idx] = Get_submat(Matrix, Year_A, Year_B, HI)
row_idx = intersect( find(HI.Virus_Year >= Year_A), find(HI.Virus_Year <= Year_B) );
col_idx = intersect( find(HI.Sera_Year >= Year_A), find(HI.Sera_Year <= Year_B) );
X = Matrix(row_idx, col_idx);








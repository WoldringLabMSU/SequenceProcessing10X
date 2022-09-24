% Importing data from the 10x spreadsheet with concatenated Frw and CDR
% regions (i.e., fwr1-cdr1-fwr2-cdr2-fw3-cdr3-fwr4) Amino acid sequences
% used
sheet = readtable('10XGEN.xlsx');
clonotypes = sheet(:,1);
clonotypes1 = table2array(clonotypes); % The number of the clonotype (i.e., clonotype 1, 2, 3, etc.)
sequences = sheet(:,2);
sequences1 = table2array(sequences);   % The concatenated CDR and Framework regions from the consensus annotations tab of the 10x spreadsheet. This is to make the chimera.
reads = sheet(:,9);
reads1 = table2array(reads);  % The reads given by 10x for each of the chains
chains = sheet(:,10);
chains1 = table2array(chains);  % The identity of the chain (IGH, IGK, or IGL)
vec_1 = zeros(size(sequences1));
Cell_1 = num2cell(vec_1);
vec_2 = zeros(size(sequences1));
Cell_2 = num2cell(vec_2);
reads2 = num2cell(reads1);

for i = 1:9937 
    if contains(chains1(i),'IGH') && contains(chains1(i+1),'IGK') && contains(chains1(i+2),'IGH') % Most common pattern to iterate through 1 clonotype. If H, K are followed by another H, this indicates that 1 clonotype has passed and the final H is the start of a new clonotype 
        Cell_1(i) = sequences1(i); % Sequence for the heavy chain of a clonotype
        Cell_2(i) = sequences1(i+1); % Sequence for the kappa light chain of the same clonotype
    elseif contains(chains1(i),'IGH') && contains(chains1(i+1),'IGL') && contains(chains1(i+2),'IGH') %See previous comment, except in this case the K light chain is an L light chain (lambda lineage)
        Cell_1(i) = sequences1(i);% Sequence for the heavy chain of a clonotype
        Cell_2(i) = sequences1(i+1); % Sequence for the lambda light chain of the same clonotype
    elseif contains(chains1(i),'IGH') && contains(chains1(i+1),'IGK') && contains(chains1(i+2),'IGK') && contains(chains1(i+3),'IGH')  % For the clonotypes that have 2 kappa light chains (by far the most common double light chain case)
        if reads1(i+1) > reads1(i+2)   % If two kappa light chains are present for a clonotype, the one with more reads is selected with this code.
            Cell_1(i) = sequences1(i);
            Cell_2(i) = sequences1(i+1);
        elseif reads1(i+2) > reads1(i+1)
            Cell_1(i) = sequences1(i);
            Cell_2(i) = sequences1(i+2);
        end
    end
end
writecell(Cell_2,'Lightchain.xlsx');    %Printing  the results to two separate spreadsheets. These now contain the VL and VH sequences which can be substituted into the whole Fab sequence for humanization
writecell(Cell_1,'Heavychain.xlsx');
        
 
     
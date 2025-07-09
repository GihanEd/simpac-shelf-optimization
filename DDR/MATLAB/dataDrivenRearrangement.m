% Function: dataDrivenRearrangement
% -------------------------------------------------------------------------
%
%   Performs layout simulations using randomized starts and 3-step shuffling
%   across multiple periods. Calculates expected profit, habitual behavior,
%   and cross-selling metrics for each simulated layout using AMPL models.
%
% Inputs:
%   NumberofTrials    - Number of trials to simulate (e.g., 1000)
%   NumberofPeriods   - Number of periods per trial (e.g., 20)
%   solver1           - Solver name for Step 1 (e.g., 'BARON')
%   solver2           - Solver name for Step 2 (e.g., 'CPLEX')
%   solver3           - Solver name for Step 3 (e.g., 'CPLEX')
%   step1ModPath      - Full path to Step 1 .mod file
%   step2ModPath      - Full path to Step 2 .mod file (single-shelf version)
%   step2AltModPath   - Full path to Step 2 alternate .mod file (two-shelf version)
%   step3ModPath      - Full path to Step 3 .mod file
%
% Outputs (Globals):
%   Analysis4_MetricsOverTime
%   Analysis4_AislesAssigned
%   Analysis4_FinalShelfAllocation
%
% Usage Example:
%   dataDrivenRearrangement(1000, 20, 'BARON', 'CPLEX', 'CPLEX', ...
%       'C:\Models\Step1_Grouping.mod', ...
%       'C:\Models\Step2_SingleShelf.mod', ...
%       'C:\Models\Step2_TwoShelf.mod', ...
%       'C:\Models\Step3_AisleAssign.mod');
%
% Required Files:
%   • AssociationData.mat
%   • InitialLayouts.mat
%   • Subcategory_FP_Growth.mat
%   • Profitability_Matrix.mat
%   • RawTransactionData.mat
%   • Department_FP_Growth.mat
%
% Notes:
%   - The .dat file for AMPL implementation is not required; all input data is generated in-code.
%   - This version allows flexible model assignment via external file paths.
%
% -------------------------------------------------------------------------

function[]=dataDrivenRearrangement(NumberofTrials,NumberofPeriods,solver1,solver2,solver3,step1ModPath, step2ModPath, step2AltModPath, step3ModPath)
tic

% ---- INPUT GLOBALS ----
global ProductGroupingDetails Profitability_SUBCATEGORYwise DepartmentList  % the inputs from 'AssociationData'
global SUBCATEGORY_SUPPORT_MATRIX                                           % the inputs from the relevant FP Growth data file (subcategory) 
global Aisles_and_Shelves Aisles_Inter_Segment_DISTANCES SUBCATEGORY_SpaceReq AisleList    % the inputs from 'AssociationData'
global RawTransactionData                                   % the inputs from the relevant RawTransactionData file
global ORIGINAL_AislesAssigned ORIGINAL_FinalShelfAllocation ORIGINAL_AisleANDShelfAssignments   % the inputs from 'AssociationData'
global DEPARTMENT_SUPPORT_MATRIX                            % the inputs from the relevant FP Growth data file
global TWENTYPERIOD_Step1_VirtualCategoryASSIGNMENT TWENTYPERIOD_Step1_VirtualCategorySUPPORTS TWENTYPERIOD_Step1_Profits TWENTYPERIOD_Step2_FeasibleIndicators TWENTYPERIOD_Step2_PSI TWENTYPERIOD_Step2_VCallocation4eachGROUPandAISLE     % Saved Grouping & PSI values for 3Step Shuffling
global FULL1997SUMMARY_Profitability_SUBCATEGORYwise FULL1997SUMMARY_SUBCATEGORY_SUPPORT_MATRIX FULL1997SUMMARY_RawTransactionData FULL1997SUMMARY_DEPARTMENT_SUPPORT_MATRIX TWENTYPERIOD_ProfitabilityData TWENTYPERIOD_SUBCATEGORY_SUPPORT_MATRIX TWENTYPERIOD_RawTransactionData TWENTYPERIOD_DEPARTMENT_SUPPORT_MATRIX      % New variable as a result of the multi-period approach

% ---- INTERMEDIATE GLOBALS ----
global Step1_VirtualCategoryASSIGNMENT Step1_VirtualCategorySUPPORTS Step1_Profits AllX Alli_cat Obj % outputs of this program
global Step2_FeasibleIndicators Step2_PSI Step2_VCallocation4eachGROUPandAISLE AllYpe AMPLiVC % outputs of this program
global Step3_Objective Step3_AislesAssigned Step3_FinalShelfAllocation AllT AMPLi AMPLFEASIBLECASES AMPLINFEASIBLECASES    % intermediate outputs of this program
global Analysis2_AisleANDShelfAssignments Analysis2_Metrics % SUPPRESS ONCE TAKEN CARE OF
global SEED_Random_AisleANDShelfAssignments    % Intermediate variable storing the 1000 or so original (seed) random allocations    

% ---- OUTPUT GLOBALS ----
global Analysis4_MetricsOverTime Analysis4_AislesAssigned Analysis4_FinalShelfAllocation     % Outputs of this program
%   - Analysis4_MetricsOverTime: Cell array where each cell stores a 3x1 array of [Metric1; Metric2; Metric3] for a period.
%   - Analysis4_AislesAssigned: Cell array where each cell stores a [DeptIndex, AisleIndex] matrix.
%   - Analysis4_FinalShelfAllocation: Cell array where each cell stores [ShelfSeg, Subcat, AllocatedSpace, k, c].



% INSTRUCTIONS FOR USE
% -----------------------------------------------------------
% PRELIMINARY SETUP:
%   Before loading data or running the analysis, initialize the AMPL API environment.
%   In MATLAB, run the following script (found in the Amplspi → Examples folder):
%       >> setupOnce
%
% STEP 1: Load data in the following order (ensure compatibility)
%   a) AssociationData.mat
%   b) InitialLayouts.mat
%   d) Subcategory_FP_Growth.mat (based on base transaction set)
%   e) Profitability_Matrix.mat (based on base transaction set)
%   f) RawTransactionData.mat (based on base transaction set)
%   g) DEPARTMENT_FP_Growth.mat (based on base transaction set)
%
% STEP 2: Run the main analysis script with the following call:
%   dataDrivenRearrangement(1000, 20, 'BARON', 'CPLEX', 'CPLEX', 'OFFICE')
%
% STEP 3: Save the following outputs manually for each completed run:
%   Analysis4_MetricsOverTime.mat
%   Analysis4_AislesAssigned.mat
%   Analysis4_FinalShelfAllocation.mat
 

%------------------ Terminology Differences ------------------------------%
% TERMINOLOGY NOTE: 'Department' in MATLAB is referred to as 'Group' in AMPL;
% 'Subcategory' in MATLAB corresponds to 'Category' in AMPL.
%-------------------------------------------------------------------------%


% A cell specifying all allocation details of the 1000 or so original (seed) random allocations
SEED_Random_AisleANDShelfAssignments=ORIGINAL_AisleANDShelfAssignments;

% A cell specifying the Analysis 4's metrics, for each peiod of each trial considered
% Column1: index of the trial
% Column2: [Metric1; Metric2; Metric3(neighbor)] for period 1
% Column3: [Metric1; Metric2; Metric3(neighbor)] for period 2
% .
% .
% .
% Column(NumberofPeriods+1): [Metric1; Metric2; Metric3(neighbor)] for period 'NumberofPeriods'
Analysis4_MetricsOverTime={};

% Stores the FINAL DEPARTMENT-AISLE ASSIGNMENT for each period
% Column1: index of the trial
% Column2 onwards: 'AislesAssinged' for each period
% [FIND 'Step3_AislesAssigned' BELOW FOR FORMAT WITHIN EACH CELL]
Analysis4_AislesAssigned={};

% Stores the FINAL SHELF SEGMENT-SUBCATEGORY ASSIGNMENT
% Column1: index of the trial
% Column2 onwards: 'FinalShelfAllocation' for each period
% [FIND 'Step3_FinalShelfAllocation' BELOW FOR FORMAT WITHIN EACH CELL]
Analysis4_FinalShelfAllocation={};

% A cell Specifying the Analysis 2's DEPARTMENT-AISLE ASSIGNMENTs in the store, for each initial layout
% Column1: index of the initial layout
% Column2: the initial layout's DEPARTMENT-AISLE ASSIGNMENT - [Column1: group (department) indices ; Column2: the aisle to which the group (department) was ORIGINALLY assigned]
% Column3: the initial layout's SHELF SEGMENT-SUBCATEGORY ASSIGNMENT - [Column 1: Shelf segments, Column 2: Assigned PURE subcategory, Column 3: Space allocated for the subcategory on this shelf segment, Column 4: 'k value' of the shelf segment, Column 5: 'c value (capacity)' of the shelf segment]
%Analysis2_AisleANDShelfAssignments={};

% A matrix specifying the Analysis 2's metrics, for each original layout considered
% Column1: index of the initial layout
% Column2: Metric#1 - The expected profit value
% Column3: Metric#2 - The historical expected profit (not applicable at the initial layout stage)
% Column4: Metric#3a - "Effort" (subcategorywise distance-association)
% Column5: Metric#3b - "Effort" (subcategorywise distance-association for each product's MOST CLOSELY ASSOCIATED PRODUCT)
% Column6: Metric#4 - The number of departments that changed from their original aisle (not applicable at the initial layout stage)
%Analysis2_Metrics=zeros((size(ORIGINAL_AisleANDShelfAssignments,1)),6);






tic
for iTrial=[1:NumberofTrials]    
% for each trial, we'll have many periods (for loop below)
% for each trial (the first two steps do not change based on the initial layout)
% They change only if the dataset (support values), the shelf layout (not allocation) and subcategory l,u values change
Analysis4_MetricsOverTime(iTrial,1)={iTrial};
Analysis4_AislesAssigned(iTrial,1)={iTrial};
Analysis4_FinalShelfAllocation(iTrial,1)={iTrial};

        for iPeriod=[1:NumberofPeriods]
        % for iLayout=[1:(size(ORIGINAL_AisleANDShelfAssignments,1))] 
        % Period 1: Original assignment
        % Periods 2 to NumberofPeriods: result of a random shuffle
        % i.e. {1 Flamand run} + {(NumberofPeriods-1) random shuffles}
        

        %Analysis2_Metrics(iLayout,1)=iLayout;
        CurrentPeriodMetrics=zeros(3,1); % Temporary holder for Metric 1–3 before saving to overall trial-period array 'Analysis4_MetricsOverTime'

        
          if iPeriod==1 % In this case, the first period arrangement is done based on Flamand
            
                
                % The values for the period would practically be all historical data till that point. 
                % However in this case we use the summary data of the full first year (1997) as the data
                % for the first year.
    
                Profitability_SUBCATEGORYwise=FULL1997SUMMARY_Profitability_SUBCATEGORYwise;
                SUBCATEGORY_SUPPORT_MATRIX=FULL1997SUMMARY_SUBCATEGORY_SUPPORT_MATRIX;
                RawTransactionData=FULL1997SUMMARY_RawTransactionData;
                DEPARTMENT_SUPPORT_MATRIX=FULL1997SUMMARY_DEPARTMENT_SUPPORT_MATRIX;
              
                        
                % A cell Specifying the ORIGINAL DEPARTMENT-AISLE ASSIGNMENT in the store
                % Column1: index of the initial layout 
                % Column2: the initial layout's DEPARTMENT-AISLE ASSIGNMENT - [Column1: group (department) indices ; Column2: the aisle to which the group (department) was ORIGINALLY assigned]
                % Column3: the initial layout's SHELF SEGMENT-SUBCATEGORY ASSIGNMENT - [Column 1: Shelf segments, Column 2: Assigned PURE subcategory, Column 3: Space allocated for the subcategory on this shelf segment, Column 4: 'k value' of the shelf segment, Column 5: 'c value (capacity)' of the shelf segment]
                ORIGINAL_AisleANDShelfAssignments={};

                % A matrix specifying the metrics of each original layout
                % Column1: index of the initial layout
                % Column2: Metric#1 - The expected profit value
                % Column3: Metric#2 - The historical expected profit value (not applicable at the initial layout stage)
                % Column4: Metric#3a - "Effort" (subcategorywise distance-association for all combinations within aisle)
                % Column5: Metric#3b - "Effort" (subcategorywise distance-association for each product's MOST CLOSELY ASSOCIATED PRODUCT)
                % Column6: Metric#4 - The number of departments that changed from their original aisle (not applicable at the initial layout stage)
                %ORIGINAL_Metrics_InitialLayouts=zeros(NumberofLayouts,6);

                

                % Assign all results to 'ORIGINAL_AisleANDShelfAssignments' - columns 2 & 3
                %ORIGINAL_AisleANDShelfAssignments(iLayout,3)={FinalSegmentAllocationforLayout}; % Column 3 (segment allocations) - see top for description
                %ORIGINAL_AisleANDShelfAssignments(iLayout,2)={AisleDeptAssignment};             % Column 2 (dept-aisle assignments) - see top for description
                %Step3_AislesAssigned=AisleDeptAssignment;
                Step3_AislesAssigned=SEED_Random_AisleANDShelfAssignments{iTrial,2};


                %%% --------- Forming the Metrics relevant to the layouts --------- %%%

                %FinalShelfAllocationforMetrics=FinalSegmentAllocationforLayout;   % Input for the 'Metric1_ExpectedProfit' program
                FinalShelfAllocationforMetrics=SEED_Random_AisleANDShelfAssignments{iTrial,3};   % Input for the 'Metric1_ExpectedProfit' program
                
                
                
                % The following values for each period are extracted from the SAME PERIOD now (as opposed to the previous period for setting this period's layout)
                Profitability_SUBCATEGORYwise=TWENTYPERIOD_ProfitabilityData{(iPeriod)}{5};
                SUBCATEGORY_SUPPORT_MATRIX=TWENTYPERIOD_SUBCATEGORY_SUPPORT_MATRIX{(iPeriod)};
                RawTransactionData=TWENTYPERIOD_RawTransactionData{(iPeriod)};
                DEPARTMENT_SUPPORT_MATRIX=TWENTYPERIOD_DEPARTMENT_SUPPORT_MATRIX{(iPeriod)};
                
                
                %------------------------ Metric#1 Calculation -----------------------%

                Metric1=0;

                for iSubCategory=(unique(FinalShelfAllocationforMetrics(:,2)))'  % for each subcategory

                    pP=Profitability_SUBCATEGORYwise(find(Profitability_SUBCATEGORYwise(:,1)==iSubCategory),2);     % The subcategory profitability
                    iP=SUBCATEGORY_SUPPORT_MATRIX(iSubCategory,iSubCategory);   % The purchase probability

                    vP=0;

                    for p=(find(FinalShelfAllocationforMetrics(:,2)==iSubCategory))'

                        vP=vP + (  (FinalShelfAllocationforMetrics(p,4))*(FinalShelfAllocationforMetrics(p,3))/(FinalShelfAllocationforMetrics(p,5))  ) ; % The visibility value

                    end

                    Metric1=Metric1+(pP*iP*vP);

                end

                CurrentPeriodMetrics(1,1)=Metric1; % First metric of this period stored in the correct format
                %---------------------------------------------------------------------%


                %------------------------ Metric#2 Calculation -----------------------%

                % Metric#2 deals with 'historical' expected cost

                Metric2=0;

                CurrentPeriodMetrics(2,1)=Metric2; % Second metric of this period stored in the correct format
                %---------------------------------------------------------------------%


                %------------------------ Metric#3 Calculation -----------------------%

                % Not applicable here, as Metric#4 deals with 'historical' expected
                % cost. But since here we simply look at initial layouts (and no
                % shuffling whatsoever).

                Metric3=0;

                for iSubCategory=(unique(FinalShelfAllocationforMetrics(:,2)))'  % for each subcategory


                    % Figure out the closest subcategories to iSubCategory

                    RelevantDept = unique(ProductGroupingDetails((find(ismember(ProductGroupingDetails(:,2), iSubCategory))),4));		% The department corresponding to iSubCategory
                    RelevantSubCategories = unique(ProductGroupingDetails((find(ismember(ProductGroupingDetails(:,4), RelevantDept))),2));			% The intra-department subcategories corresponding to iSubCategory

                    for SubCat1=iSubCategory

                     SUBCATandMinDIST=[];

                        for SubCat2=RelevantSubCategories'

                            if SubCat1==SubCat2

                                SUBCATandMinDIST=SUBCATandMinDIST;

                            else

                                Support=SUBCATEGORY_SUPPORT_MATRIX(SubCat1,SubCat2);
                                AllDistances=[];

                                for SubCat1seg=(FinalShelfAllocationforMetrics(find(FinalShelfAllocationforMetrics(:,2)==SubCat1),1))'
                                    for SubCat2seg=(FinalShelfAllocationforMetrics(find(FinalShelfAllocationforMetrics(:,2)==SubCat2),1))'

                                        AllDistances=[AllDistances Aisles_Inter_Segment_DISTANCES(SubCat1seg,SubCat2seg)];

                                    end
                                end

                                MinDistance=min(AllDistances);

                                % Forming a matrix with each (relevant) subcategory index and their MinDistance to iSubCategory
                                SUBCATandMinDIST=[SUBCATandMinDIST;SubCat2 MinDistance];

                            end

                        end
                    end

                    iSubCategory;
                    ClosestSUBCATs=SUBCATandMinDIST(find(SUBCATandMinDIST(:,2)==min(SUBCATandMinDIST(:,2))),1);		
                    % This may comprise SUBCATs that are directly opposite of iSubCategory, AND NOT IN THE SAME SHELF
                    % The following part takes care of this by isolating only those (if any) in the same shelf as iSubCategory

                    iSubCatSEGMENT=min(FinalShelfAllocationforMetrics(find(FinalShelfAllocationforMetrics(:,2)==iSubCategory),1)) ;
                    iSubCatSHELF=Aisles_and_Shelves(find(Aisles_and_Shelves(:,3)==iSubCatSEGMENT),2) ;
                    SEGMENTSiniSubCatSHELF=Aisles_and_Shelves(find(Aisles_and_Shelves(:,2)==iSubCatSHELF),3) ;

                    SUBCATSiniSubCatSHELF=[];
                    for i=SEGMENTSiniSubCatSHELF'

                        SUBCATSiniSubCatSHELF=[SUBCATSiniSubCatSHELF ; (FinalShelfAllocationforMetrics(find(FinalShelfAllocationforMetrics(:,1)==i),2))];

                    end
                    SUBCATSiniSubCatSHELF=unique(SUBCATSiniSubCatSHELF);

                    ClosestSUBCATs=intersect(SUBCATSiniSubCatSHELF,ClosestSUBCATs);		% ASSUMPTION: selecting only intra-shelf items as a possible 'ClosestSUBCATs'



                        if numel(ClosestSUBCATs)==0

                            Metric3=Metric3+0; % In case there's no neighboring subcategory within the same aisle, we ASSUME there no Cross Selling to account for

                        else

                            CrossSellingTEMP=0;

                                for iNeighbor=ClosestSUBCATs'		
                                % for the closest subcategories, do the following and get the AVERAGE (adding would give unnecessary emphasize in random shuffle) them all up	

                                    pP=Profitability_SUBCATEGORYwise(find(Profitability_SUBCATEGORYwise(:,1)==iNeighbor),2);     % The subcategory profitability
                                    iP=SUBCATEGORY_SUPPORT_MATRIX(iSubCategory,iNeighbor);   % The association between the two

                                    vP=0;

                                    % Here we add up the 'visibility' values of the neighbor AND iSubcategory together(because it is practical)
                                    NeighborANDiSubCategory=[( find(FinalShelfAllocationforMetrics(:,2)==iSubCategory) );( find(FinalShelfAllocationforMetrics(:,2)==iNeighbor) )];

                                    for p=NeighborANDiSubCategory'

                                        vP=vP + (  (FinalShelfAllocationforMetrics(p,4))*(FinalShelfAllocationforMetrics(p,3))/(FinalShelfAllocationforMetrics(p,5))  ) ; % The visibility value

                                    end

                                    CrossSellingTEMP=CrossSellingTEMP+(pP*iP*vP); % Accumulating each neighbor's 'cross selling metric'

                                end

                            Metric3=Metric3+(CrossSellingTEMP/(numel(ClosestSUBCATs)));  % Average of the 'cross selling with closest neighbor' values corresponding to current iSubCategory

                        end

                end

                %Analysis2_Metrics(iLayout,2)=Metric1;
                CurrentPeriodMetrics(3,1)=Metric3; % Third metric of this period stored in the correct format                    
                %---------------------------------------------------------------------%                    

                %%% --------------------------------------------------------------- %%%
              


                        %%--- Assigning this period's three metrics to 'Analysis4_MetricsOverTime'---%%

                        Analysis4_MetricsOverTime{iTrial,(1+iPeriod)}=CurrentPeriodMetrics;

                        %%---------------------------------------------------------------------------%%

                        %% Storing the Aisle & Shelf Allocations for Visualization Purposes
                        
                        Analysis4_AislesAssigned{iTrial,(1+iPeriod)}=Step3_AislesAssigned;
                        Analysis4_FinalShelfAllocation{iTrial,(1+iPeriod)}=Step3_FinalShelfAllocation;
                        
                        

                        %%% --------------------------------------------------------------- %%%

%                         %% Completion
% 
%                         % Notify the completion of Step 3
%                         tTemp=(toc/3600);
%                         fprintf('Step_3 Successfully Completed: %d hours.\n',  tTemp)  
                
                ORIGINAL_AislesAssigned=Step3_AislesAssigned; % Updating the 'original' aisle assignment with the current result, for the use of the next period         
                        
                        
          else  % In this case, for Periods 2 to NumberofPeriods - we do our 3 Step shuffling

                    
                   % The values for each period (2 and beyond) is extracted from our saved database as follows:
                   Profitability_SUBCATEGORYwise=TWENTYPERIOD_ProfitabilityData{(iPeriod-1)}{5};
                   SUBCATEGORY_SUPPORT_MATRIX=TWENTYPERIOD_SUBCATEGORY_SUPPORT_MATRIX{(iPeriod-1)};
                   RawTransactionData=TWENTYPERIOD_RawTransactionData{(iPeriod-1)};
                   DEPARTMENT_SUPPORT_MATRIX=TWENTYPERIOD_DEPARTMENT_SUPPORT_MATRIX{(iPeriod-1)};
              
              
                if (iTrial==1 && iPeriod==2) % We'll carry out Steps1 & 2 of our process only ONCE (i.e the first time we do it: Period #2 of Trial #1)
                              % There onwards we simply repeat Step3, based on historical layout AND Step 2 results 
                    
                        % A cell Specifying the Analysis 2's DEPARTMENT-AISLE ASSIGNMENTs in the store, for each initial layout
                        % Column1: index of the initial layout
                        % Column2: the initial layout's DEPARTMENT-AISLE ASSIGNMENT - [Column1: group (department) indices ; Column2: the aisle to which the group (department) was ORIGINALLY assigned]
                        % Column3: the initial layout's SHELF SEGMENT-SUBCATEGORY ASSIGNMENT - [Column 1: Shelf segments, Column 2: Assigned PURE subcategory, Column 3: Space allocated for the subcategory on this shelf segment, Column 4: 'k value' of the shelf segment, Column 5: 'c value (capacity)' of the shelf segment]
                        %Analysis2_AisleANDShelfAssignments={};

                        % A matrix specifying the Analysis 2's metrics, for each original layout considered
                        % Column1: index of the initial layout
                        % Column2: Metric#1  - The expected profit value
                        % Column3: Metric#2  - The historical expected profit value (not applicable at the initial layout stage)
                        % Column4: Metric#3a - "Effort" (subcategorywise distance-association)
                        % Column5: Metric#3b - "Effort" (subcategorywise distance-association for each product's MOST CLOSELY ASSOCIATED PRODUCT)
                        % Column6: Metric#4  - The number of departments that changed from their original aisle (not applicable at the initial layout stage)
                        %Analysis2_Metrics=zeros((size(ORIGINAL_AisleANDShelfAssignments,1)),6);




                        %% Step 1 - Creation of pairwise Virtual Categories

                        % Stores the OBJECTIVE VALUE of each groupwise run
                        % Column1: group indices ; Column2: corresponding objective value 
                        Step1_Profits=zeros((size(DepartmentList,1)),2);
                        Step1_Profits(:,1)=DepartmentList;

                        % Stores the SUB-CATEGORY ASSIGNMENT to VIRTUAL CATEGORIES after each groupwise run
                        % Column1: group indices ; Column2: Virtual category indices; Column3: sub-categories corresponding to each virtual category 
                        Step1_VirtualCategoryASSIGNMENT={}; % Cell

                        % Stores the inter-subcategory supports corresponding to each VIRTUAL CATEGORY of each groupwise run
                        % This is the intersection; only a result - no further use. In Step 2, we calculate the union for allocation purposes.
                        % Column1: group indices ; Column2: Virtual category indices; Column3: inter-subcategory supports corresponding to each virtual category
                        Step1_VirtualCategorySUPPORTS=[];


                        Step1_TOTALVirtualCategoryCOUNT=0; % Initialization of the TOTALVirtualCategoryCOUNT

                        for iDept=1:size(DepartmentList,1)
                        %for iDept=1:5

                        %--------------------- Start of AMPL API ---------------------------------%

                        format long g

                                % Create an AMPL instance
                                ampl = AMPL;

                                % Change to solver of choice (BARON)
                                ampl.setOption('solver', solver1)
                                %ampl.setOption('baron_options','maxtime=-1') % Setting unlimited max computation time (for a more accurate result)
                                %ampl.setOption('solver', 'BARON')



                                %------- N data, with Zg set -------% 

                                % Create appropriate entities in AMPL for the upcoming assignmentS
                                % (similar to the initial writing we do)
                                ampl.eval('set AMPLZg; param AMPLN {AMPLZg} >= 0;');


                                % Assign data to the Set 'AMPLZg'
                                % Note: The assignment must be stored in a cell array, not a standard matrix

                                % Product subcategories included under the current Department
                                RelevantSubCategories = unique(ProductGroupingDetails((find(ismember(ProductGroupingDetails(:,4), iDept))),2));
                                RelevantSubCategoriesCOUNT=size(RelevantSubCategories,1);

                                if RelevantSubCategoriesCOUNT == 1      % If only one subcategory in Dept. --> form only one virtual category 
                                    Q=0;
                                    R=1;
                                    VirtualCategoryCOUNT=Q+R;

                                elseif RelevantSubCategoriesCOUNT == 2  % If only two subcategories in Dept. --> form two, one-subcategory virtual categories, so that each could go in each shelf
                                    Q=0;
                                    R=2;
                                    VirtualCategoryCOUNT=Q+R;

                                else
                                    Q=fix(RelevantSubCategoriesCOUNT/2);   % Quotient
                                    R=rem(RelevantSubCategoriesCOUNT,2);   % Remainder
                                    VirtualCategoryCOUNT=Q+R;              % The number of duple pairing (if odd, an additional one is introduced);
                                end

                        %       AMPLZg = [100 200 300 400 500];
                                AMPLZg=[(Step1_TOTALVirtualCategoryCOUNT+1):(Step1_TOTALVirtualCategoryCOUNT+VirtualCategoryCOUNT)];
                                CurrentVirtualCategories=AMPLZg;
                                Step1_TOTALVirtualCategoryCOUNT=Step1_TOTALVirtualCategoryCOUNT+VirtualCategoryCOUNT;

                                % Assign data to the Param 'AMPLN'
                        %       AMPLN = [2 2 2 2 1];
                                AMPLN = [repmat(2,1,Q) repmat(1,1,R)];

                                % Create dataframe for Set 'AMPLZg' and Param 'AMPLN'
                                df = DataFrame(1, 'AMPLZg', 'AMPLN');
                                df.setColumn('AMPLZg', num2cell(AMPLZg));
                                df.setColumn('AMPLN', AMPLN);


                                % Set the values to the parameter
                                ampl.setData(df);

                                % Further assign Set data for 'AMPLCART'
                                ampl.setData(df, 'AMPLZg');

                        % %     % Check that data for the above Params are correctly assigned
                        % %     ampl.display('AMPLc')

                                %------------------------------------------% 




                                %------- P_prod data, with PRODg set -------% 

                                % Create appropriate entities in AMPL for the upcoming assignmentS
                                % (similar to the initial writing we do)
                                ampl.eval('set AMPLPRODg; param AMPLP_prod {AMPLPRODg} >= 0;');


                                % Assign data to the Set 'AMPLPRODg'
                                % Note: The assignment must be stored in a cell array, not a standard matrix
                        %       AMPLPRODg = [4 5 6 7 8 9 10 11 12];
                                AMPLPRODg = RelevantSubCategories';

                                % Assign data to the Param 'AMPLP_prod'
                        %       AMPLP_prod = [1.25 6.35 4.80 5.34 4.25 6.35 4.80 5.34 15.34];
                        %         AMPLP_prod = (SUBCATEGORY_Profitability(RelevantSubCategories))';
                                AMPLP_prod = (Profitability_SUBCATEGORYwise(RelevantSubCategories,2))'*1000000; % 1000 is used here as a scaling factor due to AMPL's sensitivity to small values; adjusted later in the process

                                % Create dataframe for Set 'AMPLPRODg' and Param 'AMPLP_prod'
                                df = DataFrame(1, 'AMPLPRODg', 'AMPLP_prod');
                                df.setColumn('AMPLPRODg', num2cell(AMPLPRODg));
                                df.setColumn('AMPLP_prod', AMPLP_prod);


                                % Set the values to the parameter
                                ampl.setData(df);

                                % Further assign Set data for 'AMPLCART'
                                ampl.setData(df, 'AMPLPRODg');

                        % %     % Check that data for the above Params are correctly assigned
                        % %     ampl.display('AMPLc')

                                %------------------------------------------%




                                %------- ASSOCIATION_prod data, with CUSTOMERS and CART sets -------% 

                                % Create appropriate entities in AMPL for the upcoming assignments
                                % (repeating sets are not written here)
                                ampl.eval('param AMPLASSOCIATION_prod {AMPLPRODg, AMPLPRODg} >= 0;');

                                indAMPLPRODg=repmat(AMPLPRODg,1,size(AMPLPRODg,2)); 
                                % creating a row vector with the AMPLPRODg matrix repeated (number of product categories) times

                                tempPRODG1=repmat(AMPLPRODg,[size(AMPLPRODg,2) 1]);
                                indAMPLPRODG1=tempPRODG1(:)';
                                % creating a row vector with each element of AMPLCUSTOMERS repeated (number of product categories) times

                                % Assign data to the Param 'AMPLASSOCIATION_prod'
                        %       AMPLASSOCIATION_prod = [0.81472	0.64543	0.04981	0.56609	0.30971	0.01161	0.19469	0.44532	0.67690; 0.64543	0.90579	0.08324	0.28723	0.28177	0.04861	0.24811	0.13318	0.53013; 0.04981	0.08324	0.12699	0.12067	0.08207	0.09361	0.12182	0.11800	0.06981; 0.56609	0.28723	0.12067	0.91338	0.44857	0.03320	0.15240	0.19140	0.83774; 0.30971	0.28177	0.08207	0.44857	0.63236	0.05709	0.03861	0.10751	0.18075; 0.01161	0.04861	0.09361	0.03320	0.05709	0.09754	0.01456	0.02449	0.07386; 0.19469	0.24811	0.12182	0.15240	0.03861	0.01456	0.27850	0.17157	0.20991; 0.44532	0.13318	0.11800	0.19140	0.10751	0.02449	0.17157	0.54688	0.20806; 0.67690	0.53013	0.06981	0.83774	0.18075	0.07386	0.20991	0.20806	0.95751];
                                AMPLASSOCIATION_prod = SUBCATEGORY_SUPPORT_MATRIX(RelevantSubCategories,RelevantSubCategories);

                                % Converting 'AMPLASSOCIATION_prod' to a long row vector in correct order
                                indAMPLASSOCIATION_prod=reshape((AMPLASSOCIATION_prod)',1,[]);

                                % Create dataframe for Set 'AMPLCUSTOMERS', 'AMPLPRODg' and Param 'AMPLASSOCIATION_prod'
                                df = DataFrame(2, 'AMPLPRODg1', 'AMPLPRODg2', 'AMPLASSOCIATION_prod');
                                df.setColumn('AMPLPRODg1', num2cell(indAMPLPRODG1));
                                df.setColumn('AMPLPRODg2', num2cell(indAMPLPRODg));
                                df.setColumn('AMPLASSOCIATION_prod', indAMPLASSOCIATION_prod);

                                % Set the values to the parameter
                                ampl.setData(df); 

                        %         % Check that data for the above Params are correctly assigned
                        %         ampl.display('AMPLASSOCIATION_prod')

                                %------------------------------------------%




                                %------- Executing AMPL Program -------% 

                                 % Load model from file the ampl model (.mod file)
                                 % The .dat file is irrelevant (since we create all the assignments above)

                                 ampl.read(step1ModPath);        % Step 1: Grouping – Virtual Category creation

                                 % Solve
                                 ampl.solve;

                                %------------------------------------------%



                                %%------------------------------- Results -------------------------------%%


                                %------- Retrieving variable X results -------% 

                                % Get the values of the variable Buy in a dataframe object
                                X = ampl.getVariable('X');
                                df = X.getValues;
                                % Print them
                                df;


                                AllXtable = table;
                                head = df.getHeaders();
                                for i=1:length(head)
                                     AllXtable.(char(head(i)))= df.getColumn(head(i));
                                end

                                AllXtable;                         % Converts the dataframe containing results for 'S' to a table
                                AllX = cell2mat(table2array(AllXtable));     % Converts the table ALLStable to a matrix (for easier indexing)
                                AllX=round(AllX); % Rounding to the nearest integer - 0 or 1

                        %       AllPurchases([(((irun-1)*M)+1):(irun*M)],[3:5])=vec2mat(AllY(:,3),3,M); % Adding the results to the planned format

                                %------------------------------------------%  

                                % Storing virtual category assignment results
                                for v=CurrentVirtualCategories

                                    Step1_VirtualCategoryASSIGNMENT((size(Step1_VirtualCategoryASSIGNMENT,1)+1),1)={iDept};
                                    Step1_VirtualCategoryASSIGNMENT((size(Step1_VirtualCategoryASSIGNMENT,1)),2)={v};
                                    Step1_VirtualCategoryASSIGNMENT((size(Step1_VirtualCategoryASSIGNMENT,1)),3)={(nonzeros(AllX(find(AllX(:,1)==v),2).*AllX(find(AllX(:,1)==v),3)))'};

                                end


                                %------- Retrieving variable i_cat results -------% 

                                % Get the values of the variable Buy in a dataframe object
                                i_cat = ampl.getVariable('i_cat');
                                df = i_cat.getValues;
                                % Print them
                                df;


                                Alli_cattable = table;
                                head = df.getHeaders();
                                for i=1:length(head)
                                     Alli_cattable.(char(head(i)))= df.getColumn(head(i));
                                end

                                Alli_cattable;                         % Converts the dataframe containing results for 'S' to a table
                                Alli_cat = cell2mat(table2array(Alli_cattable));     % Converts the table ALLStable to a matrix (for easier indexing)
                                Alli_cat;

                        %       AllPurchases([(((irun-1)*M)+1):(irun*M)],[3:5])=vec2mat(AllY(:,3),3,M); % Adding the results to the planned format

                                %------------------------------------------%   

                                % Storing virtual category inter-subcategory support values
                                for v=CurrentVirtualCategories

                                    Step1_VirtualCategorySUPPORTS((size(Step1_VirtualCategorySUPPORTS,1)+1),1)=iDept;
                                    Step1_VirtualCategorySUPPORTS((size(Step1_VirtualCategorySUPPORTS,1)),2)=v;
                                    Step1_VirtualCategorySUPPORTS((size(Step1_VirtualCategorySUPPORTS,1)),3)=Alli_cat(find(Alli_cat(:,1)==v),2);

                                end


                                %------- Retrieving objective value -------% 

                                tempobj = ampl.getObjective('Obj');
                                Obj = tempobj.value;
                                %Obj=round(Obj*100)/100 % Rounding to the nearest cent

                        %       AllProfits(irun,2)=Obj;
                                Step1_Profits(iDept,2)=Obj; % Storing the objective value of the corresponding Department

                                %------------------------------------------% 



                                % Close the AMPL object
                                ampl.close();

                        %         % Print the level of progress as a percentage, and the time elapsed so far
                        %         fprintf('Step1 - Current level of completion: %d percent.\n', round((iDept/(size(DepartmentList,1))*100)))
                        %         tTemp=(toc/3600);
                        %         fprintf('Elapsed time: %d hours.\n',  tTemp)        
                        %         remTemp=((toc/3600)/(iDept)*((size(DepartmentList,1))-iDept));
                        %         fprintf('Remaining time: %d hours.\n',  remTemp)

                        fprintf('Step_1 in progress...')

                        end     


                        % Notify the completion of Step 1
                        tTemp=(toc/3600);
                        fprintf('Step_1 Successfully Completed: %d hours.\n',  tTemp)  




                        %% Step 2 - Space Allocation for Virtual Categories


                        % Stores the OBJECTIVE VALUE of each (group,aisle) set, as the PSI value
                        % Rows: Group (i.e. Department) ; Columns: Aisles
                        Step2_PSI=zeros((size(DepartmentList,1)),(size(AisleList,1)));

                        % Stores the FEASIBILITY (1 iff feasible or 0 iff not feasible) of each (group,aisle) set, as the PSI value
                        % Rows: Group (i.e. Department) ; Columns: Aisles
                        Step2_FeasibleIndicators=zeros((size(DepartmentList,1)),(size(AisleList,1)));

                        % Stores the VIRTUAL CATEGORY SHELF ALLOCATION corresponding to each feasible (g,a) combination
                        % Column1: GROUP indices ; Column2: AISLE indices; Column3: [(VC indices)  (shelf indices to which each VC is assigned)]
                        Step2_VCallocation4eachGROUPandAISLE={}; % Cell


                        for iDept=1:size(DepartmentList,1)
                        %for iDept=1:1
                        %for iDept=1:5

                          for iAisle=1:size(AisleList,1)
                          %for iAisle=1:1
                          %for iAisle=1:5

                          iDept;
                          iAisle;

                            % Aisles_and_Shelves includes details about Aisles, shelves and shelf segments
                            % Column1: Aisle indices ; Column2: Corresponding Shelf indices ; Column3: Corresponding Shelf Segment indices ; Column4: Corresponding Shelf Segment CAPACITY values ; Column5: Corresponding Shelf Segment 'K' values
                            % Initially we look at a simple case of 5 aisles - 9 shelves (4 2shelf & 1 1shelf aisles) - 80 shelf segments - 80 capacity values - 80 k values

                            % SUBCATEGORY_SpaceReq includes details about Subcategories - their 'l', 'u' and 'PHI' values
                            % Column1: 'l' values ; Column2: , 'u' values ; Column3: 'PHI' values
                            % Initially we look at a simple case - but a lot thought about setting 'l' & 'u' values


                            % Partial calculations needed for feasibility checks
                            % If the sum of all 'l' values of the department's subcategories are smaller than the sum of all 'c' values of the aisle's shelf segments
                            Min_Dept_Space_Req=sum(SUBCATEGORY_SpaceReq(unique(ProductGroupingDetails((find(ismember(ProductGroupingDetails(:,4), iDept))),2)),1));
                            Available_Aisle_Space=sum(Aisles_and_Shelves((find(Aisles_and_Shelves(:,1)==iAisle)),4));

                            % The number of shelves is the currently considered aisle
                            % This is useful later to define values and direct to the correct AMPL file based on if the aisle is a 2-shelf one or 1-shelf one
                            NumberofShelvesinAisle=numel(unique(Aisles_and_Shelves((find(Aisles_and_Shelves(:,1)==iAisle)),2)));

                            
                            
                            
                            
                            %%%---- Establishing Feasiblity of Dept for Aisle ----%%%
 
                            Min_Dept_Space_Req=sum(SUBCATEGORY_SpaceReq(unique(ProductGroupingDetails((find(ismember(ProductGroupingDetails(:,4), iDept))),2)),1));
                            Available_Aisle_Space=sum(Aisles_and_Shelves((find(Aisles_and_Shelves(:,1)==iAisle)),4));

                            % The number of shelves is the currently considered aisle
                            % This is useful later to define values and direct to the correct AMPL file based on if the aisle is a 2-shelf one or 1-shelf one
                            NumberofShelvesinAisle=numel(unique(Aisles_and_Shelves((find(Aisles_and_Shelves(:,1)==iAisle)),2)));


                            if NumberofShelvesinAisle==1

                                
                                if (Min_Dept_Space_Req > Available_Aisle_Space) 

                                    CaseFeasibility=0;

                                else

                                    CaseFeasibility=1;		

                                end

                                
                            else    % The number of shelves in the aisle is 2

                                
                                AMPLZg = (unique(cell2mat(Step1_VirtualCategoryASSIGNMENT((find(ismember(cell2mat(Step1_VirtualCategoryASSIGNMENT(:,1)), iDept))),2))))';

                                AMPLl = [];     % Obtained by adding 'l' values corresponding to the subcategories in each VIRTUAL CATEGORY
                                for iTemp=AMPLZg

                                    Subcategories_of_VC=(cell2mat(Step1_VirtualCategoryASSIGNMENT((find(ismember(cell2mat(Step1_VirtualCategoryASSIGNMENT(:,2)), iTemp))),3)));
                                    l_of_VC = sum(SUBCATEGORY_SpaceReq((Subcategories_of_VC),1));
                                        AMPLl=[AMPLl l_of_VC]; 

                                end
                                AMPLl;


                                % Assigning data to the set Eb
                %               AMPLEb = [1 2];
                                AMPLEb = (unique(Aisles_and_Shelves((find(ismember(Aisles_and_Shelves(:,1), iAisle))),2)))';

                                % Assigning data to the set Eb1
                %               AMPLEb1 = [17 18 19 20 21 22 23 24];
                                AMPLEb1 = (unique(Aisles_and_Shelves((find(ismember(Aisles_and_Shelves(:,2), AMPLEb(1)))),3)))';

                                % Assigning data to the set Eb2
                %               AMPLEb2 = [25 26 27 28 29 30 31 32];
                                AMPLEb2 = (unique(Aisles_and_Shelves((find(ismember(Aisles_and_Shelves(:,2), AMPLEb(2)))),3)))';


                                AMPLc1 = sum(((Aisles_and_Shelves((find(ismember(Aisles_and_Shelves(:,3), AMPLEb1))),4)))');
                                AMPLc2 = sum(((Aisles_and_Shelves((find(ismember(Aisles_and_Shelves(:,3), AMPLEb2))),4)))');
                                
                                
                                
                                ShelfCount1=ceil(numel(AMPLZg)/NumberofShelvesinAisle);
                                ShelfCount2=floor(numel(AMPLZg)/NumberofShelvesinAisle);

                                AllPossibleLOWERCases=perms(AMPLl);

                                TwoCaseFeasibility=zeros((size(AllPossibleLOWERCases,1)),2);
                                for iPossible=[1:(size(AllPossibleLOWERCases,1))]


                                    % Feasibility of adding more VCs to the first shelf and less to the second
                                    if (  ( sum(AllPossibleLOWERCases(iPossible,[1:ShelfCount1]))<=AMPLc1 ) && ( sum(AllPossibleLOWERCases(iPossible,[(ShelfCount1+1):end]))<=AMPLc2 )  )

                                        TwoCaseFeasibility(iPossible,1)=1;

                                    else

                                        TwoCaseFeasibility(iPossible,1)=0;

                                    end


                                    % Feasibility of adding more VCs to the second shelf and less to the first
                                    if (  ( sum(AllPossibleLOWERCases(iPossible,[1:ShelfCount1]))<=AMPLc2 ) && ( sum(AllPossibleLOWERCases(iPossible,[(ShelfCount1+1):end]))<=AMPLc1 )  )

                                        TwoCaseFeasibility(iPossible,2)=1;

                                    else

                                        TwoCaseFeasibility(iPossible,2)=0;

                                    end

                                end


                                [rowIdcs, colIdcs] = find(TwoCaseFeasibility~=0);


                                if ( size(rowIdcs,1)==0 )

                                    CaseFeasibility=0;

                                else

                                    CaseFeasibility=1;		

                                end



                            end

                            %%%----------------------------------------%%%

                            if CaseFeasibility==0       % (g,a) is not feasible; set relevant PSI value to zero

                                  % set Feasibleindicator to zero
                                  Step2_FeasibleIndicators(iDept,iAisle)=0;
                                  Step2_PSI(iDept,iAisle)=0;


                            else % if feasible, run the Step2 code and store the relevant PSI value.

                                  % set Feasibleindicator to one
                                  Step2_FeasibleIndicators(iDept,iAisle)=1;


                                  %--------------------- Start of AMPL API ---------------------------------%

                                  format long g

                                  % Create an AMPL instance
                                  ampl = AMPL;

                                  % Change to solver of choice (CPLEX)
                                  ampl.setOption('solver', solver2)
                                  %ampl.setOption('solver', 'CPLEX')



                                  %------- Defining the Eb, Eb1, Eb2 Sets -------%
                                  %-------- For a 2-shelf per aisle case --------%

                                  % Applies only when 2-shelf aisles are considered
                                  if NumberofShelvesinAisle==2

                                          %%%%%%-------- See AMPL for detailed descriptions---------%%%%%%

                                          % Create appropriate entities in AMPL for the upcoming assignments
                                          % (repeating sets are not written here)
                                          ampl.eval('set AMPLEb;');
                                          ampl.eval('set AMPLEb1;');
                                          ampl.eval('set AMPLEb2;');


                                          % Assigning data to the set Eb
                        %                 AMPLEb = [1 2];
                                          AMPLEb = (unique(Aisles_and_Shelves((find(ismember(Aisles_and_Shelves(:,1), iAisle))),2)))';
                                          ampl.eval(strcat('data; set AMPLEb := ',num2str(AMPLEb),'; model;'));

                                          % Assigning data to the set Eb1
                        %                 AMPLEb1 = [17 18 19 20 21 22 23 24];
                                          AMPLEb1 = (unique(Aisles_and_Shelves((find(ismember(Aisles_and_Shelves(:,2), AMPLEb(1)))),3)))';
                                          ampl.eval(strcat('data; set AMPLEb1 := ',num2str(AMPLEb1),'; model;'));

                                          % Assigning data to the set Eb2
                        %                 AMPLEb2 = [25 26 27 28 29 30 31 32];
                                          AMPLEb2 = (unique(Aisles_and_Shelves((find(ismember(Aisles_and_Shelves(:,2), AMPLEb(2)))),3)))';
                                          ampl.eval(strcat('data; set AMPLEb2 := ',num2str(AMPLEb2),'; model;'));

                                  end        

                                  %----------------------------------------------------------%





                                  %------- Defining the ALPHAa, BETAa, FIRSTSHELFa, SECONDSHELFa Params -------% 
                                  %--------- For a 2-shelf per aisle case ----------% 

                                  %%%%%%-------- See AMPL for detailed descriptions---------%%%%%%

                                  % Create appropriate entities in AMPL for the upcoming assignments
                                  % (repeating sets are not written here)
                                  ampl.eval('param AMPLALPHAa >= 0; param AMPLBETAa >= 0;');

                                  % All shelf segment indices in the current aisle (this assigned on AMPL later)
                                  AMPLEa=(unique(Aisles_and_Shelves((find(ismember(Aisles_and_Shelves(:,1), iAisle))),3)))';

                                  % Assigning data to the param AMPLALPHAa
                        %         AMPLALPHAa = [17];
                                  AMPLALPHAa = min(AMPLEa);
                                  ampl.eval(strcat('data; param AMPLALPHAa := ',num2str(AMPLALPHAa),'; model;'));

                                  % Assigning data to the param AMPLBETAa
                        %         AMPLBETAa = [32];
                                  AMPLBETAa = max(AMPLEa);
                                  ampl.eval(strcat('data; param AMPLBETAa := ',num2str(AMPLBETAa),'; model;'));


                                  % Applies only when 2-shelf aisles are considered
                                  if NumberofShelvesinAisle==2

                                          % Create appropriate entities in AMPL for the upcoming assignments
                                          % (repeating sets are not written here)
                                          ampl.eval('param AMPLFIRSTSHELFa >= 0; param AMPLSECONDSHELFa >= 0;');


                                          % Assigning data to the param AMPLFIRSTSHELFa
                        %                 AMPLFIRSTSHELFa = [1];
                                          AMPLFIRSTSHELFa = min(AMPLEb);
                                          ampl.eval(strcat('data; param AMPLFIRSTSHELFa := ',num2str(AMPLFIRSTSHELFa),'; model;'));

                                          % Assigning data to the param AMPLSECONDSHELFa
                        %                 AMPLSECONDSHELFa = [2];
                                          AMPLSECONDSHELFa = max(AMPLEb);
                                          ampl.eval(strcat('data; param AMPLSECONDSHELFa := ',num2str(AMPLSECONDSHELFa),'; model;'));  

                                  end

                                  %----------------------------------------------------------%





                                  %------- c data, with Ea set -------% 

                                  %%%%%%-------- See AMPL for detailed descriptions---------%%%%%%

                                  % Create appropriate entities in AMPL for the upcoming assignmentS
                                  % (similar to the initial writing we do)
                                  ampl.eval('set AMPLEa; param AMPLc {AMPLEa} >= 0;');

                                  % Would have to extract from Step1 results
                        %         AMPLEa = [17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
                                  AMPLEa = AMPLEa;

                                  % Assign data to the Param 'AMPLc'
                        %         AMPLc = [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4];
                                  AMPLc = ((Aisles_and_Shelves((find(ismember(Aisles_and_Shelves(:,3), AMPLEa))),4)))';

                                  % Create dataframe for Set 'AMPLZg' and Param 'AMPLN'
                                  df = DataFrame(1, 'AMPLEa', 'AMPLc');
                                  df.setColumn('AMPLEa', num2cell(AMPLEa));
                                  df.setColumn('AMPLc', AMPLc);


                                  % Set the values to the parameter
                                  ampl.setData(df);

                                  % Further assign Set data for 'AMPLCART'
                                  ampl.setData(df, 'AMPLEa');

                        % %       % Check that data for the above Params are correctly assigned
                        % %       ampl.display('AMPLc')

                                  %------------------------------------------%





                                  %------- k data, with Ea set -------% 

                                  %%%%%%-------- See AMPL for detailed descriptions---------%%%%%%

                                  % Create appropriate entities in AMPL for the upcoming assignmentS
                                  % (similar to the initial writing we do)
                                  ampl.eval('param AMPLk {AMPLEa} >= 0;');

                                  % AMPLEa is defined above

                                  % Assign data to the Param 'AMPLk'
                        %         AMPLk = [0.80  0.35  0.45  0.10  0.15  0.55  0.40  0.90  0.80  0.35  0.45  0.10  0.15  0.55  0.40  0.90];
                                  AMPLk = ((Aisles_and_Shelves((find(ismember(Aisles_and_Shelves(:,3), AMPLEa))),5)))';

                                  % Create dataframe for Set 'AMPLEa' and Param 'AMPLk'
                                  df = DataFrame(1, 'AMPLEa', 'AMPLk');
                                  df.setColumn('AMPLEa', num2cell(AMPLEa));
                                  df.setColumn('AMPLk', AMPLk);


                                  % Set the values to the parameter
                                  ampl.setData(df);

                        % %       % Check that data for the above Params are correctly assigned
                        % %       ampl.display('AMPLc')

                                  %------------------------------------------%





                                  %------- P data, with Zg set -------% 

                                  %%%%%%- Applied based on VIRTUAL CATERGORIES, not subcategories-%%%%%%

                                  % Create appropriate entities in AMPL for the upcoming assignmentS
                                  % (similar to the initial writing we do)
                                  ampl.eval('set AMPLZg; param AMPLP {AMPLZg} >= 0;');

                                  % Exacting the VIRTUAL CATEGORIES corresponding iDept, from Step1 results
                        %         AMPLZg = [4 5 6 7 8 9 10 11];
                                  AMPLZg = (unique(cell2mat(Step1_VirtualCategoryASSIGNMENT((find(ismember(cell2mat(Step1_VirtualCategoryASSIGNMENT(:,1)), iDept))),2))))';

                                  % Assign data to the Param 'AMPLP' (would've to use a version of the code for AMPLl, AMPLu)
                        %         AMPLP = [4.25 6.35 4.80 5.34 4.25 6.35 4.80 5.34]; %%%%% UNDECIDED YET %%%%%

                                  AMPLP = [];

                                  for iVC=AMPLZg % for each VIRTUAL CATEGORY

                                        Subcategories_of_VC=(cell2mat(Step1_VirtualCategoryASSIGNMENT((find(ismember(cell2mat(Step1_VirtualCategoryASSIGNMENT(:,2)), iVC))),3)));

                                        TotalProfit_of_VC=0;
                                        %TotalSales_of_VC=0;

                                        iColumn=4;  % Directs to the correct column in 'RawTransactionData' - FOLLOW THE RIGHT FORMAT ALWAYS                
                                        for iSubCategory=Subcategories_of_VC    % Considers each subcategory in the relevant virtual category

                                            for iRow=(find(RawTransactionData(:,iColumn)==iSubCategory))'	% Considers each instance the relevant subcategory appears

                                                CurrentTransactionProfit=( RawTransactionData(iRow,8)-RawTransactionData(iRow,10) );
                                                %CurrentTransactionSales=RawTransactionData(iRow,8);

                                                TotalProfit_of_VC=(TotalProfit_of_VC+CurrentTransactionProfit);
                                                %TotalSales_of_VC=(TotalSales_of_VC+CurrentTransactionSales);   

                                            end

                                        end                

                                  TotalSales=sum(RawTransactionData(:,8)); % This is the total sales volume -  accounts for the VOLUME OF PROFITS brought in by the item in consideration    
                                        
                                  P_of_VC = (( TotalProfit_of_VC/TotalSales )*1000);    % This metric reflects profit margin only; purchase frequency is incorporated separately through support values;
                                                                                        % 1000 is used here as a scaling factor due to AMPL's sensitivity to small values; adjusted later in the process
                                  AMPLP=[AMPLP P_of_VC]; 

                                  end
                                  AMPLP;



                                  % Create dataframe for Set 'AMPLZg' and Param 'AMPLN'
                                  df = DataFrame(1, 'AMPLZg', 'AMPLP');
                                  df.setColumn('AMPLZg', num2cell(AMPLZg));
                                  df.setColumn('AMPLP', AMPLP);


                                  % Set the values to the parameter
                                  ampl.setData(df);

                                  % Further assign Set data for 'AMPLZg'
                                  ampl.setData(df, 'AMPLZg');

                        % %       % Check that data for the above Params are correctly assigned
                        % %       ampl.display('AMPLc')

                                  %------------------------------------------%





                                  %------- i data, with Zg set -------% 

                                  %%%%%%- Applied based on VIRTUAL CATERGORIES, not subcategories-%%%%%%

                                  % Create appropriate entities in AMPL for the upcoming assignmentS
                                  % (similar to the initial writing we do)
                                  ampl.eval('param AMPLi {AMPLZg} >= 0;');

                                  % AMPLZg is defined above

                                  % Assign data to the Param 'AMPLi'
                        %         AMPLi = [0.13  0.005  0.21  0.03  0.23  0.015  0.61  0.78];

                                  % The 'i' of each VC is assumed to be the probability of the union of the subcategories in it
                                  % P (A U B) = P(Buying only A)+P(Buying only B)+P(Buying both A & B)

                                  AMPLiVC = [];
                                  for iVCat=AMPLZg % for each VIRTUAL CATEGORY

                                        Subcategories_of_VC=(cell2mat(Step1_VirtualCategoryASSIGNMENT((find(ismember(cell2mat(Step1_VirtualCategoryASSIGNMENT(:,2)), iVCat))),3)));

                                        %   if the VC has only one subcategory, directly the support of that subcategory is used
                                        if size(Subcategories_of_VC,2)==1

                                            i_of_VC=SUBCATEGORY_SUPPORT_MATRIX(Subcategories_of_VC,Subcategories_of_VC);

                                        %   if the VC has two subcategories, the 'i' of each VC is assumed to be the probability of the union of the subcategories in it
                                        elseif  size(Subcategories_of_VC,2)==2

                        %                   i_of_VC=SUBCATEGORY_SUPPORT_MATRIX(Subcategories_of_VC(1),Subcategories_of_VC(2));
                                            i_of_VC=(SUBCATEGORY_SUPPORT_MATRIX(Subcategories_of_VC(1),Subcategories_of_VC(1)))+(SUBCATEGORY_SUPPORT_MATRIX(Subcategories_of_VC(2),Subcategories_of_VC(2)))-(SUBCATEGORY_SUPPORT_MATRIX(Subcategories_of_VC(1),Subcategories_of_VC(2)));

                                        end

                                  AMPLiVC=[AMPLiVC i_of_VC];      

                                  end

                                  AMPLiVC;

                                  % Create dataframe for Set 'AMPLZg' and Param 'AMPLN'
                                  df = DataFrame(1, 'AMPLZg', 'AMPLi');
                                  df.setColumn('AMPLZg', num2cell(AMPLZg));
                                  df.setColumn('AMPLi', AMPLiVC);


                                  % Set the values to the parameter
                                  ampl.setData(df);

                        % %       % Check that data for the above Params are correctly assigned
                        % %       ampl.display('AMPLc')

                                  %------------------------------------------%





                                  %------- l data, with Zg set -------% 

                                  % Create appropriate entities in AMPL for the upcoming assignmentS
                                  % (similar to the initial writing we do)
                                  ampl.eval('param AMPLl {AMPLZg} >= 0;');

                                  % AMPLZg is defined above

                                  % Assign data to the Param 'AMPLi'
                        %         AMPLl = [6  2  6  3  6  2  6  3];
                                  AMPLl = [];     % Obtained by adding 'l' values corresponding to the subcategories in each VIRTUAL CATEGORY
                                  for iTemp=AMPLZg

                                        Subcategories_of_VC=(cell2mat(Step1_VirtualCategoryASSIGNMENT((find(ismember(cell2mat(Step1_VirtualCategoryASSIGNMENT(:,2)), iTemp))),3)));
                                        l_of_VC = sum(SUBCATEGORY_SpaceReq((Subcategories_of_VC),1));
                                        AMPLl=[AMPLl l_of_VC]; 

                                  end
                                  AMPLl;

                                  % Create dataframe for Set 'AMPLZg' and Param 'AMPLN'
                                  df = DataFrame(1, 'AMPLZg', 'AMPLl');
                                  df.setColumn('AMPLZg', num2cell(AMPLZg));
                                  df.setColumn('AMPLl', AMPLl);


                                  % Set the values to the parameter
                                  ampl.setData(df);

                        % %       % Check that data for the above Params are correctly assigned
                        % %       ampl.display('AMPLc')

                                  %------------------------------------------%





                                  %------- u data, with Zg set -------% 

                                  % Create appropriate entities in AMPL for the upcoming assignmentS
                                  % (similar to the initial writing we do)
                                  ampl.eval('param AMPLu {AMPLZg} >= 0;');

                                  % AMPLZg is defined above

                                  % Assign data to the Param 'AMPLu'
                        %         AMPLu = [10  2  14  6  10  2  14  6];
                                  AMPLu = [];     % Obtained by adding 'u' values corresponding to the subcategories in each VIRTUAL CATEGORY
                                  for iTemp=AMPLZg

                                        Subcategories_of_VC=(cell2mat(Step1_VirtualCategoryASSIGNMENT((find(ismember(cell2mat(Step1_VirtualCategoryASSIGNMENT(:,2)), iTemp))),3)));
                                        u_of_VC = sum(SUBCATEGORY_SpaceReq((Subcategories_of_VC),2));
                                        AMPLu=[AMPLu u_of_VC]; 

                                  end
                                  AMPLu;

                                  % Create dataframe for Set 'AMPLZg' and Param 'AMPLN'
                                  df = DataFrame(1, 'AMPLZg', 'AMPLu');
                                  df.setColumn('AMPLZg', num2cell(AMPLZg));
                                  df.setColumn('AMPLu', AMPLu);


                                  % Set the values to the parameter
                                  ampl.setData(df);

                        % %       % Check that data for the above Params are correctly assigned
                        % %       ampl.display('AMPLc')

                                  %------------------------------------------%





                                  %------- PHI data, with Zg set -------% 

                                  % Create appropriate entities in AMPL for the upcoming assignmentS
                                  % (similar to the initial writing we do)
                                  ampl.eval('param AMPLPHI {AMPLZg} >= 0;');

                                  % AMPLZg is defined above

                                  % Assign data to the Param 'AMPLPHI'
                        %         AMPLPHI = [2  1  2  1  2  1  2  1];
                                  AMPLPHI = [];     % Since we should have space for at least one of each; this is obtained by adding 'PHI' values corresponding to the subcategories in each VIRTUAL CATEGORY
                                  for iTemp=AMPLZg

                                        Subcategories_of_VC=(cell2mat(Step1_VirtualCategoryASSIGNMENT((find(ismember(cell2mat(Step1_VirtualCategoryASSIGNMENT(:,2)), iTemp))),3)));
                                        PHI_of_VC = sum(SUBCATEGORY_SpaceReq((Subcategories_of_VC),3));
                                        AMPLPHI=[AMPLPHI PHI_of_VC]; 

                                  end
                                  AMPLPHI;

                                  % Create dataframe for Set 'AMPLZg' and Param 'AMPLN'
                                  df = DataFrame(1, 'AMPLZg', 'AMPLPHI');
                                  df.setColumn('AMPLZg', num2cell(AMPLZg));
                                  df.setColumn('AMPLPHI', AMPLPHI);


                                  % Set the values to the parameter
                                  ampl.setData(df);

                        % %       % Check that data for the above Params are correctly assigned
                        % %       ampl.display('AMPLc')

                                  %------------------------------------------%





                                  %------- Executing AMPL Program -------% 

                                  % When considering 1-shelf aisles
                                  if NumberofShelvesinAisle==1          

                                          % Load model from file the ampl model (.mod file)
                                          % The .dat file is irrelevant (since we create all the assignments above)

                                          ampl.read(step2ModPath);        % Step 2: Space allocation for single-shelf aisles

                                  end


                                  % When considering 2-shelf aisles
                                  if NumberofShelvesinAisle==2          

                                          % Load model from file the ampl model (.mod file)
                                          % The .dat file is irrelevant (since we create all the assignments above)

                                          ampl.read(step2AltModPath);     % Step 2 (Alt): Space allocation for 2-shelf aisles

                                  end       


                                  % Solve
                                  ampl.solve;

                                  %------------------------------------------%



                                  %%------------------------------- Results -------------------------------%%


                        %%%%%%%%%%%%% Changes incorporated for pure-category outputs  %%%%%%%%%%%%%%%%%%%%%

                                  %------- Retrieving variable Ype results -------% 

                                  % Get the values of the variable Buy in a dataframe object
                                  Ype = ampl.getVariable('Ype');
                                  df = Ype.getValues;
                                  % Print them
                                  df;


                                  AllYpetable = table;
                                  head = df.getHeaders();
                                  for i=1:length(head)
                                     AllYpetable.(char(head(i)))= df.getColumn(head(i));
                                  end

                                  AllYpetable;                         % Converts the dataframe containing results for 'S' to a table
                                  AllYpe = cell2mat(table2array(AllYpetable));     % Converts the table ALLStable to a matrix (for easier indexing)
                                  AllYpe=round(AllYpe); % Rounding to the nearest integer - 0 or 1

                        %         AllPurchases([(((irun-1)*M)+1):(irun*M)],[3:5])=vec2mat(AllY(:,3),3,M); % Adding the results to the planned format          

                                  %------------------------------------------% 

                                  % Organizing the Ype matrix for storage - see line 36 for description of column 3 (which tempYs will fill)
                                  tempYs=[];
                                  tempYs(:,1)=AllYpe(:,1);
                                  tempYs(:,2)=AllYpe(:,2).*AllYpe(:,3);
                                  tempYs=tempYs((find(tempYs(:,2))),:);




                                  %------- Retrieving variable Spe results -------% 

                                  % Get the values of the variable Buy in a dataframe object
                                  Spe = ampl.getVariable('Spe');
                                  df = Spe.getValues;
                                  % Print them
                                  df;


                                  AllSpetable = table;
                                  head = df.getHeaders();
                                  for i=1:length(head)
                                     AllSpetable.(char(head(i)))= df.getColumn(head(i));
                                  end

                                  AllSpetable;                         % Converts the dataframe containing results for 'S' to a table
                                  AllSpe = cell2mat(table2array(AllSpetable));     % Converts the table ALLStable to a matrix (for easier indexing)
                                  AllSpe=round(AllSpe); % Rounding to the nearest integer - 0 or 1

                        %         AllPurchases([(((irun-1)*M)+1):(irun*M)],[3:5])=vec2mat(AllY(:,3),3,M); % Adding the results to the planned format          

                                  %------------------------------------------% 

                                  % Organizing the Ype matrix for storage - see line 36 for description of column 3 (which tempYs will fill)
                                  tempSs=[];
                                  tempSs(:,1)=AllSpe(:,1);
                                  tempSs(:,2)=AllSpe(:,3).*AllYpe(:,3);
                                  tempSs=tempSs((find(tempSs(:,2))),:);
                                  tempYs(:,3)=tempSs(:,2);

                        %           % Storing virtual category shelf allocation results
                        %           Step2_VCallocation4eachGROUPandAISLE((size(Step2_VCallocation4eachGROUPandAISLE,1)),4)={tempSs};


                                  % Storing virtual category shelf allocation results
                                  Step2_VCallocation4eachGROUPandAISLE((size(Step2_VCallocation4eachGROUPandAISLE,1)+1),1)={iDept};
                                  Step2_VCallocation4eachGROUPandAISLE((size(Step2_VCallocation4eachGROUPandAISLE,1)),2)={iAisle};
                                  Step2_VCallocation4eachGROUPandAISLE((size(Step2_VCallocation4eachGROUPandAISLE,1)),3)={tempYs};



                                  %--- Converting VC allocations to original subcategory allocations ---%           

                                  SubCatAllocationsofGroup=[];

                                  for iVCat=AMPLZg % for each VIRTUAL CATEGORY in the current group

                                        Subcategories_of_VC=(cell2mat(Step1_VirtualCategoryASSIGNMENT((find(ismember(cell2mat(Step1_VirtualCategoryASSIGNMENT(:,2)), iVCat))),3)));
                                        SubCatAllocationsofVC=[];

                                        %   if the VC has only one subcategory, ................
                                        if size(Subcategories_of_VC,2)==1

                                           SubCatAllocationsofVC(:,2)=tempYs(( find(tempYs(:,1)==iVCat) ),2);
                                           SubCatAllocationsofVC(:,3)=tempYs(( find(tempYs(:,1)==iVCat) ),3);
                                           SubCatAllocationsofVC(:,1)=Subcategories_of_VC(1);

                                        %   if the VC has two subcategories, ................
                                        elseif  size(Subcategories_of_VC,2)==2

                                            % Temp matrix storing the subcategories of the current VC and their profitabilities
                                            tempSCprofitability=[];
                                            tempSCprofitability(:,1)=Subcategories_of_VC';
                                            tempSCprofitability(:,2)=Profitability_SUBCATEGORYwise(Subcategories_of_VC',2);

                                            HigherProfitSC=tempSCprofitability(find(tempSCprofitability(:,2)==max(tempSCprofitability(:,2))),1);
                                            LowerProfitSC=setdiff(Subcategories_of_VC',HigherProfitSC);


                                            % VCspacereservation holds the segments, their k values and the allocated space on each segment...
                                            % ... FOR THE CURRENT VC
                                            VCspacereservation=[];
                                            VCspacereservation(:,1)=tempYs(find(tempYs(:,1)==iVCat),2);
                                            VCspacereservation(:,3)=tempYs(find(tempYs(:,1)==iVCat),3);

                                            for seg=(VCspacereservation(:,1))'
                                                VCspacereservation( find(VCspacereservation(:,1)==seg) ,2)=Aisles_and_Shelves( find(Aisles_and_Shelves(:,3)==seg) ,5);
                                            end



                                            % VCsubcatspaces holds the VC's subcategories, their l,u and phi values...
                                            % ... FOR THE CURRENT VC
                                            VCsubcatspaces=[];
                                            VCsubcatspaces(:,1)=Subcategories_of_VC';

                                            for isubcat=(VCsubcatspaces(:,1))'
                                                VCsubcatspaces( find(VCsubcatspaces(:,1)==isubcat) ,2 )=SUBCATEGORY_SpaceReq(isubcat,1) ;
                                                VCsubcatspaces( find(VCsubcatspaces(:,1)==isubcat) ,3 )=SUBCATEGORY_SpaceReq(isubcat,2) ;
                                                VCsubcatspaces( find(VCsubcatspaces(:,1)==isubcat) ,4 )=SUBCATEGORY_SpaceReq(isubcat,3) ;
                                            end



                                            % We will allocate the less profitable item first. It will be aloocated the lowest possible amount (l)
                                            % also since this is the lesser profitable item, we allocate it to the groups of adjacent cells having the lowest 'k' values
                                            % For this we will check two options: ascending segment indices and descending segment indices of the reserved (by Step 2) area.
                                            % Once done, the balance segments will be assigned for the more profitable item.

                                            lLessProfitable=VCsubcatspaces(find(VCsubcatspaces(:,1)==LowerProfitSC),2);


                                            % OPTION 1
                                            % Evaluating the 'k' values, if the less profitable item is allocated starting from the lowest indexed segment

                                            ReservationOption1 = sortrows(VCspacereservation,1);
                                            cumsum1=cumsum(ReservationOption1(:,3));


                                            for ikval=1:size(cumsum1,1)

                                                if cumsum1(ikval)>lLessProfitable
                                                break
                                                end

                                            end
                                            ikval;


                                            LowerAllocationOption1=[];

                                            if ikval==1

                                                LowerAllocationOption1(:,2)=ReservationOption1([1:(ikval)],1);
                                                LowerAllocationOption1(:,1)=LowerProfitSC;
                                                LowerAllocationOption1(ikval,3)=lLessProfitable;

                                                LocationWeightedKvalues1=sum((ReservationOption1([1:(ikval)],2)).*(LowerAllocationOption1(:,3)));

                                            elseif cumsum1(ikval-1)==lLessProfitable

                                                LowerAllocationOption1(:,2)=ReservationOption1([1:(ikval-1)],1);
                                                LowerAllocationOption1(:,1)=LowerProfitSC;
                                                LowerAllocationOption1(:,3)=ReservationOption1([1:(ikval-1)],3);

                                                LocationWeightedKvalues1=sum((ReservationOption1([1:(ikval-1)],2)).*(LowerAllocationOption1(:,3)));

                                            elseif cumsum1(ikval-1)<lLessProfitable

                                                LowerAllocationOption1(:,2)=ReservationOption1([1:(ikval)],1);
                                                LowerAllocationOption1(:,1)=LowerProfitSC;
                                                LowerAllocationOption1([1:(ikval-1)],3)=ReservationOption1([1:(ikval-1)],3);
                                                LowerAllocationOption1(ikval,3)=lLessProfitable-sum(ReservationOption1([1:(ikval-1)],3));

                                                LocationWeightedKvalues1=sum((ReservationOption1([1:(ikval)],2)).*(LowerAllocationOption1(:,3)));

                                            end

                                            LowerAllocationOption1;
                                            LocationWeightedKvalues1;




                                            % OPTION 2
                                            % Evaluating the 'k' values, if the less profitable item is allocated starting from the HIGHEST indexed segment

                                            ReservationOption2 = sortrows(VCspacereservation,-1);
                                            cumsum2=cumsum(ReservationOption2(:,3));


                                            for ikval=1:size(cumsum2,1)

                                                if cumsum2(ikval)>lLessProfitable
                                                break
                                                end

                                            end
                                            ikval;


                                            LowerAllocationOption2=[];

                                            if ikval==1

                                                LowerAllocationOption2(:,2)=ReservationOption2([1:(ikval)],1);
                                                LowerAllocationOption2(:,1)=LowerProfitSC;
                                                LowerAllocationOption2(ikval,3)=lLessProfitable;

                                                LocationWeightedKvalues2=sum((ReservationOption2([1:(ikval)],2)).*(LowerAllocationOption2(:,3)));

                                            elseif cumsum2(ikval-1)==lLessProfitable

                                                LowerAllocationOption2(:,2)=ReservationOption2([1:(ikval-1)],1);
                                                LowerAllocationOption2(:,1)=LowerProfitSC;
                                                LowerAllocationOption2(:,3)=ReservationOption2([1:(ikval-1)],3);

                                                LocationWeightedKvalues2=sum((ReservationOption2([1:(ikval-1)],2)).*(LowerAllocationOption2(:,3)));

                                            elseif cumsum2(ikval-1)<lLessProfitable

                                                LowerAllocationOption2(:,2)=ReservationOption2([1:(ikval)],1);
                                                LowerAllocationOption2(:,1)=LowerProfitSC;
                                                LowerAllocationOption2([1:(ikval-1)],3)=ReservationOption2([1:(ikval-1)],3);
                                                LowerAllocationOption2(ikval,3)=lLessProfitable-sum(ReservationOption2([1:(ikval-1)],3));

                                                LocationWeightedKvalues2=sum((ReservationOption2([1:(ikval)],2)).*(LowerAllocationOption2(:,3)));

                                            end

                                            LowerAllocationOption2;
                                            LocationWeightedKvalues2;




                                            % Evaluating which of the above options is the best fit for the less profitable item (the one with the lower 'LocationWeightedKvalue' is always the best option, as we allocate 'l' for the less profitable one on the less attractive consecutive shelves)
                                            % Then, the more profitable item is allocated accordingly, to the remaining space 'reserved' as per our Step 2 solution

                                            HigherAllocation=[];
                                            tempLowerAllocation=[];

                                            if LocationWeightedKvalues1 <= LocationWeightedKvalues2     % (i.e. Option 1 is the 'best' allocation for the less profitable item)

                                                HigherAllocation(:,2)=ReservationOption1(:,1);
                                                HigherAllocation(:,1)=HigherProfitSC;

                                                tempLowerAllocation=[(LowerAllocationOption1(:,3)); zeros(size(ReservationOption1,1)-size(LowerAllocationOption1,1),1) ];
                                                HigherAllocation(:,3)=(ReservationOption1(:,3)-tempLowerAllocation) ;

                                                HigherAllocation=HigherAllocation((find(HigherAllocation(:,3))),:) ;

                                                LowerAllocationOption1;
                                                HigherAllocation ;                       
                                                SubCatAllocationsofVC=[LowerAllocationOption1; HigherAllocation];

                                            elseif LocationWeightedKvalues1 > LocationWeightedKvalues2  % (i.e. Option 2 is the 'best' allocation for the less profitable item)

                                                HigherAllocation(:,2)=ReservationOption2(:,1);
                                                HigherAllocation(:,1)=HigherProfitSC;

                                                tempLowerAllocation=[(LowerAllocationOption2(:,3)); zeros(size(ReservationOption2,1)-size(LowerAllocationOption2,1),1) ];
                                                HigherAllocation(:,3)=(ReservationOption2(:,3)-tempLowerAllocation) ;

                                                HigherAllocation=HigherAllocation((find(HigherAllocation(:,3))),:) ;

                                                LowerAllocationOption2;
                                                HigherAllocation;
                                                SubCatAllocationsofVC=[LowerAllocationOption2; HigherAllocation];

                                            end

                                        end

                                  SubCatAllocationsofGroup = [SubCatAllocationsofGroup; SubCatAllocationsofVC];

                                  end          

                                  SubCatAllocationsofGroup=sortrows(SubCatAllocationsofGroup,2);

                                  Step2_VCallocation4eachGROUPandAISLE((size(Step2_VCallocation4eachGROUPandAISLE,1)),4)={SubCatAllocationsofGroup};

                                  %--------------------------------------------------------------%           




                        %%%%%%%%%%%%%%% Change 'MERGED' file TO this point  %%%%%%%%%%%%%%%%%%%%%%% 



                                  %------- Retrieving objective value -------% 

                                  tempobj = ampl.getObjective('Obj');
                                  Obj = tempobj.value;
                                  %Obj=round(Obj*100)/100 % Rounding to the nearest cent

                        %         AllProfits(irun,2)=Obj;
                                  Step2_PSI(iDept,iAisle)=Obj; % Storing the objective value AS Step2_PSI
                                  % This is Step2_PSI(g,a)

                                  %------------------------------------------%          


                                  % Close the AMPL object
                                  ampl.close();


                            end

                            fprintf('Step_2 in progress...')

                          end
                        end

                        %%%%%-- SAVING Steps 1 & 2 3Step results for future calling by Step 3 --%%%%%

                        Trial1_3Step_Step1_VirtualCategoryASSIGNMENT=Step1_VirtualCategoryASSIGNMENT;
                        Trial1_3Step_Step1_VirtualCategorySUPPORTS=Step1_VirtualCategorySUPPORTS; 
                        Trial1_3Step_Step1_Profits=Step1_Profits;

                        Trial1_3Step_Step2_FeasibleIndicators=Step2_FeasibleIndicators;
                        Trial1_3Step_Step2_PSI=Step2_PSI;
                        Trial1_3Step_Step2_VCallocation4eachGROUPandAISLE=Step2_VCallocation4eachGROUPandAISLE;

                        %%%%%---------------------------------------------------------------------%%%%%

                       

                        % Notify the completion of Step 2
                        tTemp=(toc/3600);
                        fprintf('Step_2 Successfully Completed: %d hours.\n',  tTemp)  


                end

                %%%%%-- Step 1 and Step 2 results are reused from Trial 1 since these steps are not repeated across trials --%%%%%

                Step1_VirtualCategoryASSIGNMENT=Trial1_3Step_Step1_VirtualCategoryASSIGNMENT;
                Step1_VirtualCategorySUPPORTS=Trial1_3Step_Step1_VirtualCategorySUPPORTS; 
                Step1_Profits=Trial1_3Step_Step1_Profits;

                Step2_FeasibleIndicators=Trial1_3Step_Step2_FeasibleIndicators;
                Step2_PSI=Trial1_3Step_Step2_PSI;
                Step2_VCallocation4eachGROUPandAISLE=Trial1_3Step_Step2_VCallocation4eachGROUPandAISLE;

                %%%%%-------------------------------------------------------------------------------%%%%% 
                
                %for iLayout=[1:(size(ORIGINAL_AisleANDShelfAssignments,1))]    
                % for each initial layout (the first two steps do not change based on yhtinitial layout)
                % They change only if the dataset (support values), the shelf layout (not allocation) and subcategory l,u values change

                        %Analysis2_Metrics(iLayout,1)=iLayout;

                        
                        %% Step 3 - Assignment of Aisles           
                        
                        
                        % Stores the FINAL DEPARTMENT-AISLE ASSIGNMENT
                        % Column1: group (department) indices ; Column2: the aisle to which the group (department) is newly assigned
                        Step3_AislesAssigned;

                        % Stores the FINAL SHELF SEGMENT-SUBCATEGORY ASSIGNMENT
                        % Format: [Column 1: Shelf segments, Column 2: Assigned PURE subcategory, Column 3: Space allocated for the subcategory on this shelf segment, Column 4: 'k value' of the shelf segment, Column 5: 'c value (capacity)' of the shelf segment]
                        Step3_FinalShelfAllocation;

                        % Specifies the ORIGINAL DEPARTMENT-AISLE ASSIGNMENT in the store
                        % Column1: group (department) indices ; Column2: the aisle to which the group (department) was ORIGINALLY assigned
                        %ORIGINAL_AislesAssigned=cell2mat(ORIGINAL_AisleANDShelfAssignments(iLayout,2));

                        % Specifies the ORIGINAL DEPARTMENT-AISLE ASSIGNMENT in the store
                        % Format: [Column 1: Shelf segments, Column 2: Assigned PURE subcategory, Column 3: Space allocated for the subcategory on this shelf segment, Column 4: 'k value' of the shelf segment, Column 5: 'c value (capacity)' of the shelf segment]
                        %ORIGINAL_FinalShelfAllocation=cell2mat(ORIGINAL_AisleANDShelfAssignments(iLayout,3));

                        AMPLFEASIBLECASES=[];
                        AMPLINFEASIBLECASES=[];

                        for i=DepartmentList'
                        %for i=1:5
                            for j=AisleList'
                            %for j=1:5

                        %	% Step 2 feasibility checks here ==> have a resultant binary (1 iff feaible, 0 iff infeasible)
                        %	% FeasibleIndicator=1 or 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %    Step2_FeasibleIndicators=datasample([0 1],1,'Weights',[0.1 0.9]); % Randomly generaetd for now 90% feasible

                             if Step2_FeasibleIndicators(i,j)==1
                        %     if Step2_FeasibleIndicators==1

                                TempPair=[i j];
                                allOneString = sprintf('%.0f,' , TempPair);
                                allOneString = allOneString(1:end-1);
                                s1='   (';
                                s2 = ')   ';
                                TempString= strcat(s1,allOneString,s2);

                                AMPLFEASIBLECASES = [AMPLFEASIBLECASES TempString];

                            elseif Step2_FeasibleIndicators(i,j)==0
                        %	elseif Step2_FeasibleIndicators==0

                                TempPair=[i j];
                                allOneString = sprintf('%.0f,' , TempPair);
                                allOneString = allOneString(1:end-1);
                                s1='   (';
                                s2 = ')   ';
                                TempString= strcat(s1,allOneString,s2);

                                AMPLINFEASIBLECASES = [AMPLINFEASIBLECASES TempString];

                            end

                           end
                        end

                        fprintf('Step_3 in progress...')

                        %--------------------- Start of AMPL API ---------------------------------%

                        format long g

                                  % Create an AMPL instance
                                  ampl = AMPL;

                                  % Change to solver of choice (CPLEX)
                                  ampl.setOption('solver', solver3)
                                  %ampl.setOption('solver', 'CPLEX')


                                  %------- PSI data, with GROUPS and AISLES sets -------% 

                                  % Create appropriate entities in AMPL for the upcoming assignments
                                  % (repeating sets are not written here)
                                  ampl.eval('set AMPLAISLES; set AMPLGROUPS; param AMPLPSI {AMPLGROUPS, AMPLAISLES} >= 0;');

                                  % Assign data to the Set 'AMPLAISLES'
                                  % Note: The assignment must be done through explicit writing, which is why concatenation is used below
                                  % (SEE NOTES ON FILE)
                        %           AMPLAISLES = [1 2 3 4 5];
                                  AMPLAISLES = AisleList';
                                  ampl.eval(strcat('data; set AMPLAISLES := ',num2str(AMPLAISLES),'; model;'));

                                  % Assign data to the Set 'AMPLGROUPS'
                                  % Note: The assignment must be done through explicit writing, which is why concatenation is used below
                                  % (SEE NOTES ON FILE)
                        %           AMPLGROUPS = [10 20 30 40 50];
                                  AMPLGROUPS = DepartmentList';
                                  ampl.eval(strcat('data; set AMPLGROUPS := ',num2str(AMPLGROUPS),'; model;'));

                                  indAMPLAISLES=repmat(AMPLAISLES,1,size(AMPLGROUPS,2)); 
                                  % creating a row vector with the AMPLGROUPS matrix repeated (number of aisles) times

                                  tempGROUPS=repmat(AMPLGROUPS,[size(AMPLAISLES,2) 1]);
                                  indAMPLGROUPS=tempGROUPS(:)';
                                  % creating a row vector with each element of AMPLAISLES repeated (number of groups) times

                                  % Assign data to the Param 'AMPLPSI'
                        %         AMPLPSI = [200	150		100		50		60;	220		120		250		150		600; 120		130		140		150		160; 140		170		160		150		160; 20		180		120		350		360];
                                  AMPLPSI = Step2_PSI;

                                  % Converting 'AMPLPSI' to a long row vector in correct order
                                  indAMPLPSI=reshape((AMPLPSI)',1,[]);

                                  % Create dataframe for Set 'AMPLAISLES', 'AMPLGROUPS' and Param 'AMPLPSI'
                                  df = DataFrame(2, 'AMPLGROUPS', 'AMPLAISLES', 'AMPLPSI');
                                  df.setColumn('AMPLGROUPS', num2cell(indAMPLGROUPS));
                                  df.setColumn('AMPLAISLES', num2cell(indAMPLAISLES));
                                  df.setColumn('AMPLPSI', indAMPLPSI);
                                  df;

                                  % Set the values to the parameter
                                  ampl.setData(df);

                                %   % Further assign Set data for 'AMPLCUSTOMERS'
                                %   ampl.setData(df, {'AMPLCUSTOMERS', 'AMPLCART'});  

                                %   % Check that data for the above Params are correctly assigned
                                %   ampl.display('AMPLPSI')

                                  %------------------------------------------%




                                  %------- i data, with GROUPS and AISLES sets -------% 

                                  % Create appropriate entities in AMPL for the upcoming assignments
                                  % (repeating sets are not written here)
                                  ampl.eval('param AMPLi {AMPLGROUPS, AMPLAISLES} >= 0;');

                                  indAMPLAISLES=repmat(AMPLAISLES,1,size(AMPLGROUPS,2)); 
                                  % creating a row vector with the AMPLGROUPS matrix repeated (number of aisles) times

                                  tempGROUPS=repmat(AMPLGROUPS,[size(AMPLAISLES,2) 1]);
                                  indAMPLGROUPS=tempGROUPS(:)';
                                  % creating a row vector with each element of AMPLAISLES repeated (number of groups) times

                                  % Assign data to the Param 'AMPLi': historically relevant association matrix
                                  % [see notes in AMPPL model for exact definition - The purchase impulse of assigning group 'g' to aisle 'a' (this is based on the group located in the same aisle, the previous time peiod)]

                                  AMPLi=zeros((size(DepartmentList,1)),(size(AisleList,1))); % Initialization

                                  for tempD=DepartmentList'
                                      for tempA=AisleList'

                                      HistoricallyRelevantDepartment=ORIGINAL_AislesAssigned( (find(ORIGINAL_AislesAssigned(:,2)==tempA)),1 );  

                                          if tempD==HistoricallyRelevantDepartment

                                              AMPLi(tempD,tempA)=0; % If certain items are unmovable we can set this value to a very high value to force stop movement

                                          else

                                              AMPLi(tempD,tempA)=DEPARTMENT_SUPPORT_MATRIX(tempD,HistoricallyRelevantDepartment); 

                                          end 

                                      end
                                  end

                        %           AMPLi = [0.3783	0.5548	0.0179	0.9198	0.3453; 0.5841	0.4001	0.8886	0.6841	0.1993; 0.7543	0.3056	0.7374	0.3116	0.8586; 0.8183	0.7516	0.1572	0.0101	0.8002; 0.5059	0.7174	0.2269	0.5907	0.9106];

                                  % Converting 'AMPLPSI' to a long row vector in correct order
                                  indAMPLi=reshape((AMPLi)',1,[]);

                                  % Create dataframe for Set 'AMPLAISLES', 'AMPLGROUPS' and Param 'AMPLi'
                                  df = DataFrame(2, 'AMPLGROUPS', 'AMPLAISLES', 'AMPLi');
                                  df.setColumn('AMPLGROUPS', num2cell(indAMPLGROUPS));
                                  df.setColumn('AMPLAISLES', num2cell(indAMPLAISLES));
                                  df.setColumn('AMPLi', indAMPLi);
                                  df;

                                  % Set the values to the parameter
                                  ampl.setData(df);

                                %   % Further assign Set data for 'AMPLCUSTOMERS'
                                %   ampl.setData(df, {'AMPLCUSTOMERS', 'AMPLCART'});  

                                  % Check that data for the above Params are correctly assigned
                                  %ampl.display('AMPLi')

                                  %------------------------------------------%




                                  %------------- Defining the AMPLFEASIBLECASES Set -------------%
                                  %----- [The AMPLFEASIBLECASES string is created in Step 2] ----% 

                                  % Create appropriate entities in AMPL for the upcoming assignments
                                  % (repeating sets are not written here)
                                  ampl.eval('set AMPLFEASIBLECASES within {AMPLGROUPS,AMPLAISLES};');

                                  % The 'Feasible Cases Coding' section constructs the 'AMPLFEASIBLECASES' string matrix using a for-loop in Step 2

                        %         AMPLFEASIBLECASES = '(10,1)   (10,2)   (10,3)   (10,4)   (10,5)   (20,1)   (20,3)   (20,5)   (30,2)   (30,4)   (40,1)   (50,4)   (50,5)' ;

                                  ampl.eval(strcat('data; set AMPLFEASIBLECASES := ',num2str(AMPLFEASIBLECASES),'; model;'));
                                  % ampl.eval('data; set AMPLCUSTOMERS := 1 2 3 4 5 6 7 8 9 10; model;')

                                  %----------------------------------------------------------%




                                  %------------- Defining the AMPLINFEASIBLECASES Set -------------%
                                  %----- [The AMPLINFEASIBLECASES string is created in Step 2] ----% 

                                  % Create appropriate entities in AMPL for the upcoming assignments
                                  % (repeating sets are not written here)
                                  ampl.eval('set AMPLINFEASIBLECASES within {AMPLGROUPS,AMPLAISLES};');

                                  % The 'Feasible Cases Coding' section constructs the 'AMPLFEASIBLECASES' string matrix using a for-loop in Step 2

                        %         AMPLINFEASIBLECASES = '(20,2)   (20,4)   (30,1)   (30,3)   (30,5)   (40,2)   (40,3)   (40,4)   (40,5)   (50,1)   (50,2)   (50,3)' ;

                                  ampl.eval(strcat('data; set AMPLINFEASIBLECASES := ',num2str(AMPLINFEASIBLECASES),'; model;'));
                                  % ampl.eval('data; set AMPLCUSTOMERS := 1 2 3 4 5 6 7 8 9 10; model;')

                                  %----------------------------------------------------------%




                                  %------- Executing AMPL Program -------% 

                                  % Load model from file the ampl model (.mod file)
                                  % The .dat file is irrelevant (since we create all the assignments above)

                                  ampl.read(step3ModPath);        % Step 3: Aisle assignment with feasibility


                                  % Solve
                                  ampl.solve;

                                  %------------------------------------------%



                        %%------------------------------- Results -------------------------------%%


                                  %------- Retrieving variable T results -------% 

                                  % Get the values of the variable Buy in a dataframe object
                                  T = ampl.getVariable('T');
                                  df = T.getValues;
                                  % Print them
                                  df;


                                  AllTtable = table;
                                  head = df.getHeaders();
                                  for i=1:length(head)
                                     AllTtable.(char(head(i)))= df.getColumn(head(i));
                                  end

                                  AllTtable;                         % Converts the dataframe containing results for 'S' to a table
                                  AllT = cell2mat(table2array(AllTtable));     % Converts the table ALLStable to a matrix (for easier indexing)
                                  AllT=round(AllT); % Rounding to the nearest integer - 0 or 1

                        %         AllPurchases([(((irun-1)*M)+1):(irun*M)],[3:5])=vec2mat(AllY(:,3),3,M); % Adding the results to the planned format          

                                  %------------------------------------------% 


                                  % Stores the FINAL DEPARTMENT-AISLE ASSIGNMENT
                                  % Column1: group (department) indices ; Column2: the aisle to which the group (department is assigned)
                                  tempCOLONE=AllT(:,1);
                                  tempCOLTWO=AllT(:,2);
                                  tempCOLTHREE=AllT(:,3);

                                  tempMulti(:,1)=tempCOLONE.*tempCOLTHREE;
                                  tempMulti(:,2)=tempCOLTWO.*tempCOLTHREE;

                                  % Given the presence of many zeros, we multiply the appropriate column of 'AllT' to retain only non-zero (assigned) cases
                                  Step3_AislesAssigned=tempMulti(find(tempMulti(:,1)),:); 



                                  %------- Retrieving objective value -------% 

                                  tempobj = ampl.getObjective('Obj');
                                  Obj = tempobj.value;
                                  %Obj=round(Obj*100)/100 % Rounding to the nearest cent

                        %       AllProfits(irun,2)=Obj;
                                  Step3_Objective=Obj; % Storing the objective value

                                  %------------------------------------------% 


                                  % Close the AMPL object
                                  ampl.close();


                        %------- Forming 'Step3_FinalShelfAllocation' -------% 

                        Step3_FinalShelfAllocation=[];

                        for icase=1:size(Step3_AislesAssigned,1)

                            icasegroup=Step3_AislesAssigned(icase,1);
                            icaseaisle=Step3_AislesAssigned(icase,2);

                            tempA=find(cell2mat(Step2_VCallocation4eachGROUPandAISLE(:,1))==icasegroup);
                            tempB=find(cell2mat(Step2_VCallocation4eachGROUPandAISLE(:,2))==icaseaisle);
                            icaserow=intersect(tempA,tempB); % The row of 'Step2_VCallocation4eachGROUPandAISLE' containing data fitting the Departmetn-Aisle assignment case

                            DeptFinalAssign=cell2mat(Step2_VCallocation4eachGROUPandAISLE(icaserow,4));
                            %%% DeptFinalAssign=cell2mat(Step2_VCallocation4eachGROUPandAISLE(icaserow,3));


                            Step3_FinalShelfAllocation=[Step3_FinalShelfAllocation; DeptFinalAssign];

                        end


                        % Processing the resultant matrix
                        % Format: [Column 1: Shelf segments, Column 2: Assigned PURE subcategory, Column 3: Space allocated for the subcategory on this shelf segment, Column 4: 'k value' of the shelf segment, Column 5: 'c value (capacity)' of the shelf segment]

                        Step3_FinalShelfAllocation(:,[1 2]) = Step3_FinalShelfAllocation(:,[2 1]); 	% Flipping the first two columns to fit the above format
                        Step3_FinalShelfAllocation=sortrows(Step3_FinalShelfAllocation,1);		% Sorting the matrix by the shelf segments in ascending order          

                        for isegmentrow=1:size(Step3_FinalShelfAllocation,1)

                            isegment=Step3_FinalShelfAllocation(isegmentrow,1);     % The shelf segment in this row

                            Step3_FinalShelfAllocation(isegmentrow,4)=Aisles_and_Shelves(find(Aisles_and_Shelves(:,3)==isegment),5);  % Storing the 'k value' (traffic) of the shelf segment
                            Step3_FinalShelfAllocation(isegmentrow,5)=Aisles_and_Shelves(find(Aisles_and_Shelves(:,3)==isegment),4);  % Storing the 'c value' (capacity) of the shelf segment

                        end

                        %------------------------------------------%          

                        % Assign all results to 'Analysis2_AisleANDShelfAssignments' - columns 2 & 3
                        %Analysis2_AisleANDShelfAssignments(iLayout,3)={Step3_FinalShelfAllocation}; % Column 3 (segment allocations) - see top for description
                        %Analysis2_AisleANDShelfAssignments(iLayout,2)={Step3_AislesAssigned};                
              
              
              
              
                    %%% --------- Forming the Metrics relevant to the layouts --------- %%%

                    % Assign all results to 'ORIGINAL_AisleANDShelfAssignments' - columns 2 & 3
                    %ORIGINAL_AisleANDShelfAssignments(iLayout,3)={FinalSegmentAllocationforLayout}; % Column 3 (segment allocations) - see top for description
                    %ORIGINAL_AisleANDShelfAssignments(iLayout,2)={AisleDeptAssignment};             % Column 2 (dept-aisle assignments) - see top for description
                    Step3_AislesAssigned=Step3_AislesAssigned;

                    
                    FinalShelfAllocationforMetrics=Step3_FinalShelfAllocation;   % Input for the 'Metric1_ExpectedProfit' program

                    
                    % The following values for each period are extracted from the SAME PERIOD now (as opposed to the previous period for setting this period's layout)
                    Profitability_SUBCATEGORYwise=TWENTYPERIOD_ProfitabilityData{(iPeriod)}{5};
                    SUBCATEGORY_SUPPORT_MATRIX=TWENTYPERIOD_SUBCATEGORY_SUPPORT_MATRIX{(iPeriod)};
                    RawTransactionData=TWENTYPERIOD_RawTransactionData{(iPeriod)};
                    DEPARTMENT_SUPPORT_MATRIX=TWENTYPERIOD_DEPARTMENT_SUPPORT_MATRIX{(iPeriod)};
                    
                    
                    %------------------------ Metric#1 Calculation -----------------------%

                    Metric1=0;

                    for iSubCategory=(unique(FinalShelfAllocationforMetrics(:,2)))'  % for each subcategory

                        pP=Profitability_SUBCATEGORYwise(find(Profitability_SUBCATEGORYwise(:,1)==iSubCategory),2);     % The subcategory profitability
                        iP=SUBCATEGORY_SUPPORT_MATRIX(iSubCategory,iSubCategory);   % The purchase probability

                        vP=0;

                        for p=(find(FinalShelfAllocationforMetrics(:,2)==iSubCategory))'

                            vP=vP + (  (FinalShelfAllocationforMetrics(p,4))*(FinalShelfAllocationforMetrics(p,3))/(FinalShelfAllocationforMetrics(p,5))  ) ; % The visibility value

                        end

                        Metric1=Metric1+(pP*iP*vP);

                    end

                    %ORIGINAL_Metrics_InitialLayouts(iLayout,2)=Metric1;   % Output from the 'Metric1_ExpectedProfit' program
                    CurrentPeriodMetrics(1,1)=Metric1; % First metric of this period stored in the correct format
                    %---------------------------------------------------------------------%


                    %------------------------ Metric#2 Calculation -----------------------%

                    % Metric#2 deals with 'historical' expected cost
                    % In our formulation, it was simply the objective of Step 3
                    % However, in this Flamand-GAP type of setup; historical impulses have nothing to do with the Step 3 objective function
                    % Therefore we need to go through the process of calculating them, here

                    Metric2=0;

            %       Metric2=Step3_Objective; % In our formulation, it was simply the objective of Step 3

                    for iAssignCase=[1:(size(Step3_AislesAssigned,1))]

                        AssignedDept=Step3_AislesAssigned(iAssignCase,1);
                        AssignedAisle=Step3_AislesAssigned(iAssignCase,2);

                        HistoricallyRelevantDepartment=ORIGINAL_AislesAssigned( (find(ORIGINAL_AislesAssigned(:,2)==AssignedAisle)),1 );  

                            if AssignedDept==HistoricallyRelevantDepartment

                                HabitualImpulse=0; % In this case, there's no 'habitually driven' impulse

                            else

                                HabitualImpulse=DEPARTMENT_SUPPORT_MATRIX(AssignedDept,HistoricallyRelevantDepartment); 

                            end 

                        tempHistoricalProfit=Step2_PSI(AssignedDept,AssignedAisle)*HabitualImpulse/1000; % /1000 is an adjustment to correct the earlier *1000 intrduced to ease AMPL work
                        Metric2=(Metric2+tempHistoricalProfit); %%%%%%%%%%%

                    end    
                    

                    CurrentPeriodMetrics(2,1)=Metric2; % Second metric of this period stored in the correct format
                    %---------------------------------------------------------------------%

                    
                    %------------------------ Metric#3 Calculation -----------------------%

                    % Considers the cross selling profitability of the new arrangement

                    Metric3=0;

                    for iSubCategory=(unique(FinalShelfAllocationforMetrics(:,2)))'  % for each subcategory


                        % Figure out the closest subcategories to iSubCategory

                        RelevantDept = unique(ProductGroupingDetails((find(ismember(ProductGroupingDetails(:,2), iSubCategory))),4));		% The department corresponding to iSubCategory
                        RelevantSubCategories = unique(ProductGroupingDetails((find(ismember(ProductGroupingDetails(:,4), RelevantDept))),2));			% The intra-department subcategories corresponding to iSubCategory

                        for SubCat1=iSubCategory

                         SUBCATandMinDIST=[];

                            for SubCat2=RelevantSubCategories'

                                if SubCat1==SubCat2

                                    SUBCATandMinDIST=SUBCATandMinDIST;

                                else

                                    Support=SUBCATEGORY_SUPPORT_MATRIX(SubCat1,SubCat2);
                                    AllDistances=[];

                                    for SubCat1seg=(FinalShelfAllocationforMetrics(find(FinalShelfAllocationforMetrics(:,2)==SubCat1),1))'
                                        for SubCat2seg=(FinalShelfAllocationforMetrics(find(FinalShelfAllocationforMetrics(:,2)==SubCat2),1))'

                                            AllDistances=[AllDistances Aisles_Inter_Segment_DISTANCES(SubCat1seg,SubCat2seg)];

                                        end
                                    end

                                    MinDistance=min(AllDistances);

                                    % Forming a matrix with each (relevant) subcategory index and their MinDistance to iSubCategory
                                    SUBCATandMinDIST=[SUBCATandMinDIST;SubCat2 MinDistance];

                                end

                            end
                        end

                        iSubCategory;
                        ClosestSUBCATs=SUBCATandMinDIST(find(SUBCATandMinDIST(:,2)==min(SUBCATandMinDIST(:,2))),1);		
                        % This may comprise SUBCATs that are directly opposite of iSubCategory, AND NOT IN THE SAME SHELF
                        % The following part takes care of this by isolating only those (if any) in the same shelf as iSubCategory

                        iSubCatSEGMENT=min(FinalShelfAllocationforMetrics(find(FinalShelfAllocationforMetrics(:,2)==iSubCategory),1)) ;
                        iSubCatSHELF=Aisles_and_Shelves(find(Aisles_and_Shelves(:,3)==iSubCatSEGMENT),2) ;
                        SEGMENTSiniSubCatSHELF=Aisles_and_Shelves(find(Aisles_and_Shelves(:,2)==iSubCatSHELF),3) ;

                        SUBCATSiniSubCatSHELF=[];
                        for i=SEGMENTSiniSubCatSHELF'

                            SUBCATSiniSubCatSHELF=[SUBCATSiniSubCatSHELF ; (FinalShelfAllocationforMetrics(find(FinalShelfAllocationforMetrics(:,1)==i),2))];

                        end
                        SUBCATSiniSubCatSHELF=unique(SUBCATSiniSubCatSHELF);

                        ClosestSUBCATs=intersect(SUBCATSiniSubCatSHELF,ClosestSUBCATs);		% Assumption: Only intra-shelf items are considered valid closest subcategories



                            if numel(ClosestSUBCATs)==0

                                Metric3=Metric3+0; % In case there's no neighboring subcategory within the same aisle, we ASSUME there no Cross Selling to account for

                            else

                                CrossSellingTEMP=0;

                                    for iNeighbor=ClosestSUBCATs'		
                                    % for the closest subcategories, do the following and get the AVERAGE (adding would give unnecessary emphasize in random shuffle) them all up	

                                        pP=Profitability_SUBCATEGORYwise(find(Profitability_SUBCATEGORYwise(:,1)==iNeighbor),2);     % The subcategory profitability
                                        iP=SUBCATEGORY_SUPPORT_MATRIX(iSubCategory,iNeighbor);   % The association between the two

                                        vP=0;

                                        % Adds visibility of neighbor and iSubcategory together (since this reflects practical consumer behavior)
                                        NeighborANDiSubCategory=[( find(FinalShelfAllocationforMetrics(:,2)==iSubCategory) );( find(FinalShelfAllocationforMetrics(:,2)==iNeighbor) )];

                                        for p=NeighborANDiSubCategory'

                                            vP=vP + (  (FinalShelfAllocationforMetrics(p,4))*(FinalShelfAllocationforMetrics(p,3))/(FinalShelfAllocationforMetrics(p,5))  ) ; % The visibility value

                                        end

                                        CrossSellingTEMP=CrossSellingTEMP+(pP*iP*vP); % Accumulating each neighbor's 'cross selling metric'

                                    end

                                Metric3=Metric3+(CrossSellingTEMP/(numel(ClosestSUBCATs)));  % Average of the 'cross selling with closest neighbor' values corresponding to current iSubCategory

                            end

                    end

                    %Analysis2_Metrics(iLayout,2)=Metric1;
                    CurrentPeriodMetrics(3,1)=Metric3; % Third metric of this period stored in the correct format                    
                    %---------------------------------------------------------------------%                    
                                                           
                    %%% --------------------------------------------------------------- %%%

                ORIGINAL_AislesAssigned=Step3_AislesAssigned; % Updating the 'original' aisle assignment with the current result, for the use of the next period
        
          end      
        
        %%--- Assigning this period's three metrics to 'Analysis4_MetricsOverTime'---%%
                
        Analysis4_MetricsOverTime{iTrial,(1+iPeriod)}=CurrentPeriodMetrics;
                
        %%---------------------------------------------------------------------------%%

        
        %% Storing the Aisle & Shelf Allocations for Visualization Purposes

        Analysis4_AislesAssigned{iTrial,(1+iPeriod)}=Step3_AislesAssigned;
        Analysis4_FinalShelfAllocation{iTrial,(1+iPeriod)}=Step3_FinalShelfAllocation;
        
        
        % Print the level of progress as a percentage, and the time elapsed so far
        fprintf('Current level of completion: %d percent.\n', round(((((iTrial-1)*NumberofPeriods)+iPeriod)/(NumberofTrials*NumberofPeriods)*100)))
        tTemp=(toc/3600);
        fprintf('Elapsed time: %d hours.\n',  tTemp)
        remTemp=((toc/3600)/((((iTrial-1)*NumberofPeriods)+iPeriod))*((NumberofTrials*NumberofPeriods)-((((iTrial-1)*NumberofPeriods)+iPeriod))));
        fprintf('Remaining time: %d hours.\n',  remTemp)
        
        end
end

end
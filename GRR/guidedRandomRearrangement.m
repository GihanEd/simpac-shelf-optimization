% guidedRandomRearrangement
% ---------------------------------------------------------
% Performs guided randomized rearrangement of department-subcategory 
% allocations across multiple time periods using profitability, 
% transaction, and support data.
%
% INPUTS:
%   NumberofTrials  — Integer. Number of randomized layout trials.
%   NumberofPeriods — Integer. Number of simulation periods (e.g., months).
%
% OUTPUTS (via global variables):
%   Analysis4_MetricsOverTime        — Metrics across trials and periods
%   Analysis4_AislesAssigned         — Aisle assignments per period
%   Analysis4_FinalShelfAllocation   — Final shelf allocation matrices
%
% USAGE:
%   guidedRandomRearrangement(1000, 20)
%
% REQUIRED DATA FILES (preloaded into workspace):
%   • AssociationData.mat
%   • InitialLayouts.mat
%   • Subcategory_FP_Growth.mat
%   • Profitability_Matrix.mat
%   • RawTransactionData.mat
%   • Department_FP_Growth.mat



function[]=guidedRandomRearrangement(NumberofTrials,NumberofPeriods)
tic

% Performs guided random layout optimization over multiple periods using randomized starting points.
%
% Inputs:
%   NumberofTrials - number of different layout trials
%   NumberofPeriods - number of periods (e.g., months/years)
%
% Outputs (in global scope):
%   Analysis4_MetricsOverTime
%   Analysis4_AislesAssigned
%   Analysis4_FinalShelfAllocation



% ---- INPUT GLOBALS ----
global ProductGroupingDetails Profitability_SUBCATEGORYwise DepartmentList  % the inputs from 'AssociationData'
global SUBCATEGORY_SUPPORT_MATRIX                                           % the inputs from the relevant FP Growth data file (subcategory) 
global Aisles_and_Shelves Aisles_Inter_Segment_DISTANCES SUBCATEGORY_SpaceReq AisleList    % the inputs from 'AssociationData'
global RawTransactionData                                   % the inputs from the relevant RawTransactionData file
global ORIGINAL_AislesAssigned ORIGINAL_FinalShelfAllocation ORIGINAL_AisleANDShelfAssignments   % the inputs from 'AssociationData'
global DEPARTMENT_SUPPORT_MATRIX                            % the inputs from the relevant FP Growth data file
global FULL1997SUMMARY_Profitability_SUBCATEGORYwise FULL1997SUMMARY_SUBCATEGORY_SUPPORT_MATRIX FULL1997SUMMARY_RawTransactionData FULL1997SUMMARY_DEPARTMENT_SUPPORT_MATRIX TWENTYPERIOD_ProfitabilityData TWENTYPERIOD_SUBCATEGORY_SUPPORT_MATRIX TWENTYPERIOD_RawTransactionData TWENTYPERIOD_DEPARTMENT_SUPPORT_MATRIX      % New variable as a result of the multi-period approach

% ---- INTERMEDIATE GLOBALS ----
global Step1_VirtualCategoryASSIGNMENT Step1_VirtualCategorySUPPORTS Step1_Profits AllX Alli_cat Obj % outputs of this program
global Step2_FeasibleIndicators Step2_PSI Step2_VCallocation4eachGROUPandAISLE AllYpe AMPLiVC % outputs of this program
global Step3_Objective Step3_AislesAssigned Step3_FinalShelfAllocation AllT AMPLi AMPLFEASIBLECASES AMPLINFEASIBLECASES    % intermediate outputs of this program
global ORIGINAL_Metrics_InitialLayouts         % Outputs of this program
global SEED_Random_AisleANDShelfAssignments    % Intermediate variable storing the 1000 or so original (seed) random allocations    
global Analysis2_AisleANDShelfAssignments Analysis2_Metrics

% ---- OUTPUT GLOBALS ----
global Analysis4_MetricsOverTime Analysis4_AislesAssigned Analysis4_FinalShelfAllocation     % Outputs of this program



% INSTRUCTIONS FOR USE
% -----------------------------------------------------------
% STEP 1: Load data in the following order (ensure compatibility)
%   a) AssociationData.mat
%   b) InitialLayouts.mat (e.g., 1000 layouts)
%   c) Subcategory_FP_Growth.mat
%   d) Profitability_Matrix.mat
%   e) RawTransactionData.mat
%   f) Department_FP_Growth.mat
%
% STEP 2: Run the main analysis script with parameters:
%   guidedRandomRearrangement(1000, 20)
%
% STEP 3: Save outputs:
%   Analysis4_MetricsOverTime.mat
%   Analysis4_AislesAssigned.mat
%   Analysis4_FinalShelfAllocation.mat
 

%------------------ Terminology Differences ------------------------------%
% TERMINOLOGY NOTE: 'Department' in MATLAB is referred to as 'Group' in AMPL;
% 'Subcategory' in MATLAB corresponds to 'Category' in AMPL.
%-------------------------------------------------------------------------%


% A cell array specifying all allocation details of the 1000 or so original (seed) random allocations
% Stores initial randomized layouts (used as baseline for period 1)
SEED_Random_AisleANDShelfAssignments=ORIGINAL_AisleANDShelfAssignments;

% A cell specifying the Analysis 4's metrics, for each period of each trial considered
% Column1: index of the trial
% Column2: [Metric1; Metric2; Metric3(neighbor)] for period 1
% Column3: [Metric1; Metric2; Metric3(neighbor)] for period 2
% .
% .
% .
% Column(NumberofPeriods+1): [Metric1; Metric2; Metric3(neighbor)] for period 'NumberofPeriods'
Analysis4_MetricsOverTime={};

% Stores the final department-to-aisle assignment for each period and trial.
% Structure:
%   - Column 1: Trial index
%   - Columns 2 onward: Cell arrays for each period
%     > Each cell contains:
%         • Column 1: Department (group) indices
%         • Column 2: Corresponding aisle assignment for each department
Analysis4_AislesAssigned = {};

% Stores the final shelf segment-to-subcategory allocation for each period and trial.
% Structure:
%   - Column 1: Trial index
%   - Columns 2 onward: Cell arrays for each period
%     > Each cell contains a matrix with:
%         • Column 1: Shelf segment IDs
%         • Column 2: Assigned pure subcategory
%         • Column 3: Space allocated for the subcategory on this shelf segment (e.g., meters or units)
%         • Column 4: 'k' value of the shelf segment (e.g., weight factor)
%         • Column 5: 'c' value (capacity) of the shelf segment
Analysis4_FinalShelfAllocation = {};



tic
for iTrial=[1:NumberofTrials]    
% Each trial contains multiple periods (see loop below)
% for each trial (the first two steps do not change based on the initial layout)
% They change only if the dataset (support values), the shelf layout (not allocation) and subcategory l,u values change
Analysis4_MetricsOverTime(iTrial,1)={iTrial};
Analysis4_AislesAssigned(iTrial,1)={iTrial};
Analysis4_FinalShelfAllocation(iTrial,1)={iTrial};

        for iPeriod=[1:NumberofPeriods]
        % for iLayout=[1:(size(ORIGINAL_AisleANDShelfAssignments,1))] 
        % Period 1: Original assignment (Flamand in this case)
        % Periods 2 to NumberofPeriods: result of a random shuffle
        % i.e. {1 initial layout} + {(NumberofPeriods-1) random schuffles}
        

        %Analysis2_AisleANDShelfAssignments(iLayout,1)={iLayout};
        %Analysis2_Metrics(iLayout,1)=iLayout;
        CurrentPeriodMetrics=zeros(3,1); % Stores the three metrics of this period, until we assign it to 'Analysis4_MetricsOverTime'

        
          if iPeriod==1     
					% Period 1: Use predefined (seed) layout assignments
					% Load full-year summary data (1997) as baseline
    
                    Profitability_SUBCATEGORYwise=FULL1997SUMMARY_Profitability_SUBCATEGORYwise;
                    SUBCATEGORY_SUPPORT_MATRIX=FULL1997SUMMARY_SUBCATEGORY_SUPPORT_MATRIX;
                    RawTransactionData=FULL1997SUMMARY_RawTransactionData;
                    DEPARTMENT_SUPPORT_MATRIX=FULL1997SUMMARY_DEPARTMENT_SUPPORT_MATRIX;
              
              
                    %% Step 3 - Assignment of Aisles

                    % Stores the FINAL DEPARTMENT-AISLE ASSIGNMENT
                    % Column1: group (department) indices ; Column2: the aisle to which the group (department) is newly assigned
                    Step3_AislesAssigned;

                    % Stores the FINAL SHELF SEGMENT-SUBCATEGORY ASSIGNMENT
                    % Format: [Column 1: Shelf segments, Column 2: Assigned PURE subcategory, Column 3: Space allocated for the subcategory on this shelf segment, Column 4: 'k value' of the shelf segment, Column 5: 'c value (capacity)' of the shelf segment]
                    Step3_FinalShelfAllocation;

                    % Specifies the ORIGINAL DEPARTMENT-AISLE ASSIGNMENT in the store
                    % Column1: group (department) indices ; Column2: the aisle to which the group (department) was ORIGINALLY assigned
                    % ORIGINAL_AislesAssigned=cell2mat(ORIGINAL_AisleANDShelfAssignments(iLayout,2));
                    ORIGINAL_AislesAssigned=[];
                    %%%%%% Needs adjusting in reshuffles (not here because this is the initial period) %%%%%% 


                    % A cell Specifying the ORIGINAL DEPARTMENT-AISLE ASSIGNMENT in the store
                    % Column1: index of the initial layout 
                    % Column2: the initial layout's DEPARTMENT-AISLE ASSIGNMENT - [Column1: group (department) indices ; Column2: the aisle to which the group (department) was ORIGINALLY assigned]
                    % Column3: the initial layout's SHELF SEGMENT-SUBCATEGORY ASSIGNMENT - [Column 1: Shelf segments, Column 2: Assigned PURE subcategory, Column 3: Space allocated for the subcategory on this shelf segment, Column 4: 'k value' of the shelf segment, Column 5: 'c value (capacity)' of the shelf segment]
                    ORIGINAL_AisleANDShelfAssignments={};

					% Assign all results to 'ORIGINAL_AisleANDShelfAssignments' - columns 2 & 3
                    %ORIGINAL_AisleANDShelfAssignments(iLayout,3)={FinalSegmentAllocationforLayout}; % Column 3 (segment allocations) - see top for description
                    %ORIGINAL_AisleANDShelfAssignments(iLayout,2)={AisleDeptAssignment};             % Column 2 (dept-aisle assignments) - see top for description
                    %Step3_AislesAssigned=AisleDeptAssignment; %%%%%%%%%%
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
					% METRIC 1: Visibility-driven Purchasing Potential
					% This metric evaluates the potential profit from visibility-enhanced purchasing behavior within a given period.				
					

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
                    CurrentPeriodMetrics(2,1)=Metric2; % Second metric of this period stored in the correct format
					% METRIC 2: Past-aisle Impulse Potential via Temporal Familiarity
					% This metric measures the likelihood of triggering impulse purchases based on aisle continuity across time.



                    %---------------------------------------------------------------------%

                    
                    %------------------------ Metric#3 Calculation -----------------------%

                    % Not applicable here, as Metric#4 deals with 'historical' expected
                    % cost. But since here we simply look at initial layouts (and no
                    % schuffling whatsoever).

                    Metric3=0;
					% METRIC 3: Present-shelf Impulse Effect
					% This metric captures the potential for impulse-driven profit increases through strategic shelf-level adjacencies.

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

                                Metric3=Metric3+0; % In case there's no neighboring subcategory within the same aisle, cross-selling potential is assumed to be zero

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
                    
                                         
                      
                %%--- Assigning this period's three metrics to 'Analysis4_MetricsOverTime'---%%

                Analysis4_MetricsOverTime{iTrial,(1+iPeriod)}=CurrentPeriodMetrics;

                %%---------------------------------------------------------------------------%%
                      
                
                %% Storing the Aisle & Shelf Allocations for Visualization Purposes
                        
                Analysis4_AislesAssigned{iTrial,(1+iPeriod)}=Step3_AislesAssigned;
                Analysis4_FinalShelfAllocation{iTrial,(1+iPeriod)}=Step3_FinalShelfAllocation;
                
                
                
                ORIGINAL_AislesAssigned=Step3_AislesAssigned; % Updating the 'original' aisle assignment with the current result, for the use of the next period         
                        
                        
          else  % In this case, for Periods 2 to NumberofPeriods - a new layout is randomly generated for each period


                % The values for each period (2 and beyond) is extracted from our saved database as follows:
                Profitability_SUBCATEGORYwise=TWENTYPERIOD_ProfitabilityData{(iPeriod-1)}{5};
                SUBCATEGORY_SUPPORT_MATRIX=TWENTYPERIOD_SUBCATEGORY_SUPPORT_MATRIX{(iPeriod-1)};
                RawTransactionData=TWENTYPERIOD_RawTransactionData{(iPeriod-1)};
                DEPARTMENT_SUPPORT_MATRIX=TWENTYPERIOD_DEPARTMENT_SUPPORT_MATRIX{(iPeriod-1)};
              
              
                % A cell Specifying the ORIGINAL DEPARTMENT-AISLE ASSIGNMENT in the store
                % Column1: index of the initial layout 
                % Column2: the initial layout's DEPARTMENT-AISLE ASSIGNMENT - [Column1: group (department) indices ; Column2: the aisle to which the group (department) was ORIGINALLY assigned]
                % Column3: the initial layout's SHELF SEGMENT-SUBCATEGORY ASSIGNMENT - [Column 1: Shelf segments, Column 2: Assigned PURE subcategory, Column 3: Space allocated for the subcategory on this shelf segment, Column 4: 'k value' of the shelf segment, Column 5: 'c value (capacity)' of the shelf segment]
                ORIGINAL_AisleANDShelfAssignments={};


                    %ORIGINAL_AisleANDShelfAssignments(iLayout,1)={iLayout};
                    %ORIGINAL_Metrics_InitialLayouts(iLayout,1)=iLayout;
                    AislesAssignedsoFar=[];
                    FinalSegmentAllocationforLayout=[];
                    AisleDeptAssignment=[];
    
    
                    % We create a sorted largest-to-smallest department list, to avoid bigger departments being considered later
                    % (when feasible aisles may no longer be available for assignment)
                    % This is a practical randomized approach to ensure that enough feasible aisles are available for the bigger departments when sampling
                    DeptandSize=[];
                    for Deptindex=DepartmentList'

                        RelevantSubCategories = unique(ProductGroupingDetails((find(ismember(ProductGroupingDetails(:,4), Deptindex))),2));
                        RelevantSubCategoriesCOUNT=size(RelevantSubCategories,1);
                        DeptandSize=[DeptandSize;Deptindex RelevantSubCategoriesCOUNT];

                    end

                    DeptandSize=sortrows(DeptandSize,2); 
                    SortedDepts=flip(DeptandSize(:,1));     % Sorting departments from heaviest (most subcategories in it) to lightest
                    %----------------------%


                    for CurrentDept=SortedDepts'    % We assign each group to a feasible aisle

                        AllocationofDEPARTMENTtoSEGMENTS=[]; % Output of this departmentwise run

                        Min_Dept_Space_Req=sum(SUBCATEGORY_SpaceReq(unique(ProductGroupingDetails((find(ismember(ProductGroupingDetails(:,4), CurrentDept))),2)),1));

                        CurrentAisle=datasample(setdiff(AisleList',AislesAssignedsoFar),1); % Randomly drawing a aisle index out of the aisles not assigned with departments, so far
                        Available_Aisle_Space=sum(Aisles_and_Shelves((find(Aisles_and_Shelves(:,1)==CurrentAisle)),4));

                        while (Min_Dept_Space_Req > Available_Aisle_Space)      % Considering ONLY feasible aisles
                            CurrentAisle=datasample(setdiff(AisleList',AislesAssignedsoFar),1);
                            Available_Aisle_Space=sum(Aisles_and_Shelves((find(Aisles_and_Shelves(:,1)==CurrentAisle)),4));
                        end
                        CurrentAisle;


                        RelevantSubCategories = unique(ProductGroupingDetails((find(ismember(ProductGroupingDetails(:,4), CurrentDept))),2));
                        RelevantSubCategoriesCOUNT=size(RelevantSubCategories,1);


                %         Available_Aisle_Space=64;RelevantSubCategoriesCOUNT=14; % N - Total aisle capacity; nE - number of subcategories in Group

                        if ((Available_Aisle_Space/RelevantSubCategoriesCOUNT)==4)  % If the aisle space divides evenly across subcategories, assign each the minimum required allocation.

                           v=4*(ones(1,RelevantSubCategoriesCOUNT));

                        elseif  (RelevantSubCategoriesCOUNT>10)   % This is a time-saving measure for departments with many subcategories

                            v1=4*(ones(1,10));

                            v2 = diff([0,sort(randperm((Available_Aisle_Space-(10*4))-1,(RelevantSubCategoriesCOUNT-10)-1)),(Available_Aisle_Space-(10*4))]);
                            while sum(v2(:) < 4) % Replace the 4 with the minimal 'l value' of the entire group
                                v2 = diff([0,sort(randperm((Available_Aisle_Space-(10*4))-1,(RelevantSubCategoriesCOUNT-10)-1)),(Available_Aisle_Space-(10*4))]);
                            end

                            v=[v1 v2];
                            v=v(randperm(length(v))); % randomly permute the allocation vector to diversify layouts

                        else    % We generate a set of random space allocation values that ADD UP TO THE AVAILABLE AISLE SPACE

                            v = diff([0,sort(randperm(Available_Aisle_Space-1,RelevantSubCategoriesCOUNT-1)),Available_Aisle_Space]);
                            while sum(v(:) < 4) % Replace the 4 with the minimal 'l value' of the entire group
                                v = diff([0,sort(randperm(Available_Aisle_Space-1,RelevantSubCategoriesCOUNT-1)),Available_Aisle_Space]);
                            end

                        end
                        v;

                        % Space allocations for the subcategories in the current department
                        % [Column 1: SubCat index , Column 2: Space reserved for that subcat] 
                        SpaceforSUBCATs=[];
                        SpaceforSUBCATs(:,1)=RelevantSubCategories;
                        SpaceforSUBCATs(:,2)=v';


                        % VCspacereservation holds the segments of the current aisle, their k values and the capacity of each segment...
                        % ... FOR THE CURRENT DEPARTMENT
                        Departmentspacereservation=[];
                        Departmentspacereservation(:,1)=Aisles_and_Shelves( find(Aisles_and_Shelves(:,1)==CurrentAisle),3 );

                        for seg=(Departmentspacereservation(:,1))'
                            Departmentspacereservation( find(Departmentspacereservation(:,1)==seg) ,2)=Aisles_and_Shelves( find(Aisles_and_Shelves(:,3)==seg) ,5);
                            Departmentspacereservation( find(Departmentspacereservation(:,1)==seg) ,3)=Aisles_and_Shelves( find(Aisles_and_Shelves(:,3)==seg) ,4);
                        end

                        Departmentspacereservation=sortrows(Departmentspacereservation,1); % Sorting by the first column            

                        % VCsubcatspaces holds the Department's subcategories, their l,u and phi values...
                        % ... FOR THE CURRENT DEPARTMENT
                        DeptSUBCATspaces=[];
                        DeptSUBCATspaces(:,1)=RelevantSubCategories;

                        for isubcat=(DeptSUBCATspaces(:,1))'
                            DeptSUBCATspaces( find(DeptSUBCATspaces(:,1)==isubcat) ,2 )=SUBCATEGORY_SpaceReq(isubcat,1) ;
                            DeptSUBCATspaces( find(DeptSUBCATspaces(:,1)==isubcat) ,3 )=SUBCATEGORY_SpaceReq(isubcat,2) ;
                            DeptSUBCATspaces( find(DeptSUBCATspaces(:,1)==isubcat) ,4 )=SUBCATEGORY_SpaceReq(isubcat,3) ;
                        end        

                        AllocationofDEPARTMENTtoSEGMENTS=[];
                        for CurrentSUBCAT=RelevantSubCategories' % for all subcats in the current dept

                                    SpaceforCurrentSubcat=SpaceforSUBCATs(find(SpaceforSUBCATs(:,1)==CurrentSUBCAT),2);    
                                    AllocationofSUBCATtoSEGMENTS=[]; % This is our output at each subcat run

                %                   SegmentSpaceAvailability = sortrows(Departmentspacereservation,1); 
                                    SegmentSpaceAvailability=Departmentspacereservation;

                                    for icase=[1:(size(AllocationofDEPARTMENTtoSEGMENTS,1))]

                                             seg=AllocationofDEPARTMENTtoSEGMENTS(icase,2);
                                             spacetaken=AllocationofDEPARTMENTtoSEGMENTS(icase,3);

                                             SegmentSpaceAvailability( find(SegmentSpaceAvailability(:,1)==seg),3 ) = (SegmentSpaceAvailability( find(SegmentSpaceAvailability(:,1)==seg),3 )) - spacetaken;
                                    end

                                    SegmentSpaceAvailability=SegmentSpaceAvailability(find(SegmentSpaceAvailability(:,3)),:); % Getting rid of 0 cases
                                    cumsum1=cumsum(SegmentSpaceAvailability(:,3));


                                    for ikval=1:size(cumsum1,1)

                                        if cumsum1(ikval)>SpaceforCurrentSubcat
                                        break
                                        end

                                    end
                                    ikval;


                                    if ikval==1

                                        AllocationofSUBCATtoSEGMENTS(:,2)=SegmentSpaceAvailability([1:(ikval)],1);
                                        AllocationofSUBCATtoSEGMENTS(:,1)=CurrentSUBCAT;
                                        AllocationofSUBCATtoSEGMENTS(ikval,3)=SpaceforCurrentSubcat;

                                    elseif cumsum1(ikval-1)==SpaceforCurrentSubcat

                                        AllocationofSUBCATtoSEGMENTS(:,2)=SegmentSpaceAvailability([1:(ikval-1)],1);
                                        AllocationofSUBCATtoSEGMENTS(:,1)=CurrentSUBCAT;
                                        AllocationofSUBCATtoSEGMENTS(:,3)=SegmentSpaceAvailability([1:(ikval-1)],3);

                                    elseif cumsum1(ikval-1)<SpaceforCurrentSubcat

                                        AllocationofSUBCATtoSEGMENTS(:,2)=SegmentSpaceAvailability([1:(ikval)],1);
                                        AllocationofSUBCATtoSEGMENTS(:,1)=CurrentSUBCAT;
                                        AllocationofSUBCATtoSEGMENTS([1:(ikval-1)],3)=SegmentSpaceAvailability([1:(ikval-1)],3);
                                        AllocationofSUBCATtoSEGMENTS(ikval,3)=SpaceforCurrentSubcat-sum(SegmentSpaceAvailability([1:(ikval-1)],3));

                                    end

                                    AllocationofDEPARTMENTtoSEGMENTS=[AllocationofDEPARTMENTtoSEGMENTS; AllocationofSUBCATtoSEGMENTS];


                          end

                        AislesAssignedsoFar=[AislesAssignedsoFar CurrentAisle]; % Record the current assignment for tracking
                        FinalSegmentAllocationforLayout=[FinalSegmentAllocationforLayout; AllocationofDEPARTMENTtoSEGMENTS];

                        AisleDeptAssignment=[AisleDeptAssignment; CurrentDept CurrentAisle];

                    end


                    % Processing the resultant matrix
                    % Format: [Column 1: Shelf segments, Column 2: Assigned PURE subcategory, Column 3: Space allocated for the subcategory on this shelf segment, Column 4: 'k value' of the shelf segment, Column 5: 'c value (capacity)' of the shelf segment]

                    FinalSegmentAllocationforLayout(:,[1 2]) = FinalSegmentAllocationforLayout(:,[2 1]); 	% Flipping the first two columns to fit the above format
                    FinalSegmentAllocationforLayout=sortrows(FinalSegmentAllocationforLayout,1);		% Sorting the matrix by the shelf segments in ascending order          

                    for isegmentrow=1:size(FinalSegmentAllocationforLayout,1)

                        isegment=FinalSegmentAllocationforLayout(isegmentrow,1);     % The shelf segment in this row

                        FinalSegmentAllocationforLayout(isegmentrow,4)=Aisles_and_Shelves(find(Aisles_and_Shelves(:,3)==isegment),5);  % Storing the 'k value' (traffic) of the shelf segment
                        FinalSegmentAllocationforLayout(isegmentrow,5)=Aisles_and_Shelves(find(Aisles_and_Shelves(:,3)==isegment),4);  % Storing the 'c value' (capacity) of the shelf segment

                    end


                    % Assign all results to 'ORIGINAL_AisleANDShelfAssignments' - columns 2 & 3
                    %ORIGINAL_AisleANDShelfAssignments(iLayout,3)={FinalSegmentAllocationforLayout}; % Column 3 (segment allocations) - see top for description
                    %ORIGINAL_AisleANDShelfAssignments(iLayout,2)={AisleDeptAssignment};             % Column 2 (dept-aisle assignments) - see top for description
                    Step3_AislesAssigned=AisleDeptAssignment;


                    %%% --------- Forming the Metrics relevant to the layouts --------- %%%

                    FinalShelfAllocationforMetrics=FinalSegmentAllocationforLayout;   % Input for the 'Metric1_ExpectedProfit' program

                    
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

                    Metric2=0;

            %       Metric2=Step3_Objective; % In our formulation, it was simply the objective of Step 3

                    for iAssignCase=[1:(size(Step3_AislesAssigned,1))]

                        AssignedDept=Step3_AislesAssigned(iAssignCase,1);
                        AssignedAisle=Step3_AislesAssigned(iAssignCase,2);

                        
                        % Habitual Impulse calculation
                        HistoricallyRelevantDepartment=ORIGINAL_AislesAssigned( (find(ORIGINAL_AislesAssigned(:,2)==AssignedAisle)),1 );  

                            if AssignedDept==HistoricallyRelevantDepartment

                                HabitualImpulse=0; % In this case, there's no 'habitually driven' impulse

                            else

                                HabitualImpulse=DEPARTMENT_SUPPORT_MATRIX(AssignedDept,HistoricallyRelevantDepartment); 

                            end 

                            
                        % ProfitofCurrentAssignment Impulse calculation (because using Step_PSI is invalid here)
                        ProfitofCurrentAssignment=0;

                        for iSubCategory=( unique(ProductGroupingDetails(find((ProductGroupingDetails(:,4))==AssignedDept),2)) )'  % for each subcategory

                            pP=Profitability_SUBCATEGORYwise(find(Profitability_SUBCATEGORYwise(:,1)==iSubCategory),2);     % The subcategory profitability
                            iP=SUBCATEGORY_SUPPORT_MATRIX(iSubCategory,iSubCategory);   % The purchase probability

                            vP=0;

                            for p=(find(FinalShelfAllocationforMetrics(:,2)==iSubCategory))'

                                vP=vP + (  (FinalShelfAllocationforMetrics(p,4))*(FinalShelfAllocationforMetrics(p,3))/(FinalShelfAllocationforMetrics(p,5))  ) ; % The visibility value

                            end

                           ProfitofCurrentAssignment=ProfitofCurrentAssignment+(pP*iP*vP);

                        end
                        
                        
                        
                        % Final Metric2 calculation    
                        tempHistoricalProfit=ProfitofCurrentAssignment*HabitualImpulse;
                        Metric2=(Metric2+tempHistoricalProfit); %%%%%%%%%%%

                    end    
                    

                    CurrentPeriodMetrics(2,1)=Metric2; % Second metric of this period stored in the correct format
                    %---------------------------------------------------------------------%

                    
                    %------------------------ Metric#3 Calculation -----------------------%

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

                                Metric3=Metric3+0; % In case there's no neighboring subcategory within the same aisle, we assume there is no cross-selling potential

                            else

                                CrossSellingTEMP=0;

                                    for iNeighbor=ClosestSUBCATs'		
                                    % for the closest subcategories, do the following and get the AVERAGE (adding would put unnecessary emphasis in random shuffle)	

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
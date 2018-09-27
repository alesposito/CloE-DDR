%% INIT
close all
clear all


% Input parameters - SIMULATION attributes

n_trees = 100;   % number of trees (initial cells) to be simulated
dt      = 1;     % time base in hours
t_max   = 24*30; % lenght of the simulation in hours
t_ncs   = 24;    % time at which DNA damage is 
bn      = 100;    % number of repeats

% Input parameter - EXPERIMENTALLY related

% Probability of mitosis will follow a normal distrubition with average
% time tm and width sm. 
tm = 24.5; % hrs
ts = 5;  % hrs (15)
div_dice = @(nc)randn(1,nc)*ts + tm;

p0_dead     = 0.0;  % rate of cell death before treatment
p_dead      = 0.55; % rate of cell death after treatment
t_dead_max  = 30;    % max delay before death
  

% Init storage variables
[memo_proliferation memo_exhausted ] = deal(zeros(ceil(t_max/dt+2),bn));

    
%%    
for bi=1:bn
    
    close all
    rng('shuffle'); 

    TreeStat = zeros(n_trees,round(t_max/dt));

    %% randomize initial birth and division times
    dt_birth = 10*(tm+3*ts);
    t0_birth = dt_birth*rand(1,n_trees);
    t0_mit   = div_dice(n_trees)+t0_birth;

    t = 0 ;
    while t<=dt_birth
        t = t + dt;

        idx = find((t>=t0_mit));
        if ~isempty(idx)
            t0_mit(idx) = t + div_dice(numel(idx));        
            t0_birth(idx) = t;
        end    
    end

    
    t0_mit   = t0_mit-max(t0_birth)+eps;
    t0_birth = t0_birth-max(t0_birth)+eps;    
    
    figure
    th = (-dt_birth:dt_birth);
    h0 = histogram(t0_birth,th);
    hold on
    hn = histogram(t0_mit,th);

    b_dead = rand(1,n_trees)<p0_dead;
    t0_mit(b_dead) = +inf;
    t0_dead = t_dead_max*rand(1,n_trees);
    t0_dead(~b_dead) = +inf;
    
    
    %% 
    t  = 0;
    ti = 1; 

    
    % Creates cell objects that can be analyzed with LinSense
    for id_cnt=1:n_trees

        CellList(id_cnt) = MockFucciCell;

        CellList(id_cnt).cellNumber = id_cnt;
        CellList(id_cnt).birthTime  = t0_birth(id_cnt); 
        CellList(id_cnt).tree_id    = id_cnt;        

        CellMit(id_cnt)   = t0_mit(id_cnt); % cell division time
        CellDead(id_cnt)  = t0_dead(id_cnt); % cell death time
        CellCheck(id_cnt) = 0; % this will become 1 or 2 when cell is either dead or divided

        TreeStat(id_cnt,ti) = 1; % each tree sart with a single cell
    end

    %% sim
    while t<=t_max

        t  = t + dt;
        ti = ti + 1;

        if t>=t_ncs
            pc_dead = p_dead;
        else
            pc_dead = p0_dead;
        end
        
        % Find active cells
        idx     = find(~CellCheck);
        n_cells = numel(idx);

        % Print to screen and define exit strategy
        if mod(ti,50)==0
            [t sum(TreeStat(:,ti-1))]
            if sum(TreeStat(:,ti-1))>2500
                TreeStat(:,ti:end)=[];
                ti = ti-1;
                break
            end
        end


        for ni=1:n_cells

           if t>= CellDead(idx(ni)); % Does cell die?
                CellList(idx(ni)).deathTime = ti;
                CellCheck(idx(ni)) = 1;                                       
           end

           if t>= CellMit(idx(ni)); % Does cell devide?

                %Update mother cell
                CellList(idx(ni)).divisionTime  = ti;
                CellList(idx(ni)).daughterCells = id_cnt+[1 2];

                %Create daughter cells
                CellList(id_cnt+1) = MockFucciCell;
                CellList(id_cnt+2) = MockFucciCell;
                CellList(id_cnt+1).cellNumber = id_cnt + 1;
                CellList(id_cnt+2).cellNumber = id_cnt + 2;                
                CellList(id_cnt+1).birthTime  = ti;
                CellList(id_cnt+2).birthTime  = ti;                    
                CellList(id_cnt+1).motherCell = idx(ni);
                CellList(id_cnt+2).motherCell = idx(ni);
                % Inherit lineage id
                CellList(id_cnt+1).tree_id = CellList(idx(ni)).tree_id;
                CellList(id_cnt+2).tree_id = CellList(idx(ni)).tree_id;            
                
               
                % Set devision time
                CellMit(id_cnt+[1 2])   = t + div_dice(2);
                CellDead(id_cnt+[1 2])  = +inf;
                
                % Set cell death
                if rand(1)<pc_dead; % first siste
                    CellMit(id_cnt+1) = +inf;
                    CellDead(id_cnt+1) = t + t_dead_max*rand(1);
                end
                if rand(1)<pc_dead; % first siste
                    CellMit(id_cnt+2) = +inf;
                    CellDead(id_cnt+2) = t + t_dead_max*rand(1);
                end
                
                
               
               CellCheck(id_cnt+[1 2]) = 0;          
               CellCheck(idx(ni)) = 2;
               id_cnt = id_cnt + 2;


            end % cell devide end



        end % scan cells

        % Update Tree statistics    
        TreeStat(:,ti) = hist([CellList(find(~CellCheck)).tree_id],(1:1:n_trees))';


    end




    %%
    figure
    subplot(1,2,1)
    imagesc(TreeStat)
    subplot(1,2,2)
    hist(TreeStat(:,end))

    exhausted_trees = nnz(TreeStat(:,end)==0)/n_trees
    proliferation = sum(TreeStat(:,end))/n_trees

    memo_proliferation(1:ti,bi) = sum(TreeStat)/n_trees;
    memo_exhausted(1:ti,bi) = sum(TreeStat==0)/n_trees;
    
    figure
    subplot(1,2,1)
    plot(sum(TreeStat)/n_trees)
    title('proliferation')
    subplot(1,2,2)
    plot(sum(TreeStat==0)/n_trees)
    title('exhausted trees')

    %%

    col_rnd = rand(n_trees,3).^.5

    sm = 100;
    time = dt*(0:ti-1)/24;
    showclone = @(ab,offset0,offset1,col)patch([time fliplr(time)],[0.5*medfilt1(ab,sm)/n_trees -0.5*fliplr(medfilt1(ab,sm))/n_trees]+[(offset0:(offset1-offset0)/(ti-1):offset1) fliplr(offset0:(offset1-offset0)/(ti-1):offset1)],col,'edgecolor','none','facealpha',.6)



    %offset_end = cumsum(TreeStat(:,end)/n_trees);

    idx_end = min(find(sum(TreeStat)==0));
    if isempty(idx_end)

        idx_end = ti;
    end


    offset_end = cumsum(TreeStat(:,idx_end)/n_trees);
    offset_end = mean(cat(2,[0; offset_end],[offset_end; max(offset_end)+eps]),2);


    offset_end = offset_end - max(offset_end)/2;
    offset_ini = (1:n_trees)/n_trees-0.5;

    TS0 = [zeros(1,ti);TreeStat];
    TS = TS0./repmat(sum(TreeStat,1)+eps,[n_trees+1 1]);
    TS0(1,:) = -sum(TS0,1)/2;

    hf1 = figure

        subplot(1,2,1)
        axis square
        box on
        set(gca,'xlim',[0 t_max/24])

        subplot(1,2,2)
        axis square
        box on
        set(gca,'xlim',[0 t_max/24],'ylim',[0 1])

    for ci=1:n_trees
        clone = TreeStat(ci,:);
        %showclone(clone,ci/n_trees,offset_end(ci),rand(1,3))
        %showclone(clone,offset_end(ci),offset_end(ci)+eps,rand(1,3))
        subplot(1,2,1)
        %showclone(clone,offset_ini(ci),offset_end(ci)+eps,col_rnd(ci,:))

        patch([time fliplr(time)], ...
              [medfilt1(sum(TS0(1:ci,:),1),sm) ....
              fliplr(medfilt1(sum(TS0(1:ci+1,:),1),sm))],col_rnd(ci,:),'edgecolor','none','facealpha',.6)

        subplot(1,2,2)

        patch([time fliplr(time)], ...
              [medfilt1(sum(TS(1:ci,:),1),sm) ....
              fliplr(medfilt1(sum(TS(1:ci+1,:),1),sm))],col_rnd(ci,:),'edgecolor','none','facealpha',.6)

        drawnow
    end


    subplot(1,2,1)
    axis square
    box on
    set(gca,'xlim',[0 t_max/24])

    axis off
    
    subplot(1,2,2)
    axis square
    box on
    set(gca,'xlim',[0 t_max/24],'ylim',[0 1])

    axis off
    
    saveas(hf1,[num2str(bi) '.png'])
    saveas(hf1,[num2str(bi) '.fig'])
    
     %export2corel(hf1,[num2str(bi)])

    clear CellList CellMit CellDead CellCheck
end


%%
time = dt*(0:size(memo_proliferation,1)-1)/24;
hf=figure
plot(time,memo_proliferation)
saveas(hf,['pr.png'])
saveas(hf,['pr.fig'])
saveas(hf,['pr.eps'])

hf=figure
plot(time,memo_exhausted)
saveas(hf,['ex.png'])
saveas(hf,['ex.fig'])
saveas(hf,['ex.eps'])
save 'backup.mat' memo_exhausted memo_proliferation


%%
hf = figure
hold all
plot(time,mean(memo_proliferation'),'b')
plot(time,mean(memo_proliferation')-std(memo_proliferation'),'b')
plot(time,mean(memo_proliferation')+std(memo_proliferation'),'b')
plot(time,mean(memo_exhausted'),'r')
plot(time,mean(memo_exhausted')-std(memo_exhausted'),'r')
plot(time,mean(memo_exhausted')+std(memo_exhausted'),'r')

title('Population growth relative to start of simulation')
ylabel('folds')
xlabel('days')

saveas(hf,['avg.fig'])
saveas(hf,['avg.png'])
saveas(hf,['avg.eps'])
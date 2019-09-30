%% THis code has been written by Akram Etemadi Amin. PhD student at Physics and Astronomy dept. 
%% University of New Mexico
%% All parameters used in this simulation has extracted from the experiments done in Neumann's lab. 


%% Parameters

%  all_dead_by_fret   the number of all dectins that have FRET
%  all_FRET_mono      combination of FRET efficiency from differnt
%  Acceptors surronding a donor
%  
%  all_Alive_Dectin     % the number of excited donors on time_step (tn)
%  all_dead             %location (coordination) and number of all donors that decay over two process









%%

if ispc
    clear all; close all;
    clc
    
    saveDir='.';
    dir=1;
    seed = randi(10000);
else
    dir=0;
end

if ~exist('seed', 'var')
    seed=1;
end
rng(seed*cputime);
%%

killed_by_FRET=0; % counting the donors that decay because of FREt
killed_by_exc=0;   %counting the donors that decay because of Emission
x_cont=[-0.4,0.4]*10^0; % container size
y_cont=[-0.4,0.4]*10^0;
x_range=max(x_cont)-min(x_cont); % container bunderies
y_range=max(x_cont)-min(y_cont);
area_cont= x_range .* y_range;

N=round(4612*area_cont); %Density and the number of particles
Diffusion=1.12; % um^2/s   Diffusion coefficient extract from experiments
fret_dist =10*10^-3; % the distance from donor that can be considered for FRET if there is an Acc
Forster_dis=5.24*10^-3;
frac=0.5; % Fraction of donors
dimer_frac =0; %The fraction of whole particles which are dimers
time_step=1e-10; % s
Pulse_dist=25e-9; % s
Pulsenum=1;  %the number of pulses over data aquiring
t_total=Pulse_dist*Pulsenum; % s   the pixel size are 0.188u*0.18u and the container size is 0.8um*0.8um
cont_Msize = 4;%contaner marker size
plot_Msize = 3;
Dead_Acc_frac=0;  % the percent to photobleach the acc that can be used when there is acc photobleaching
Alive_Dectin=[];
% Dimer_dist=0.003; %um
Pulse2plot=3;% the number of pulses for graphing 
time_cnt = 1;
% cnt_2=1;

makevid = 1; % for making a video should be 1
one_pulse_FRET=0;
one_pulse_Exc=0;
% dimmers_loc_offset = 5e-3;% the distance between dimers
pulse_cnt=0;
init_loc_type=2; % define the type of structure, patterned or random distribution (2 for random)
r=0.05;% the radious for patterns
time_life=[]; % donors life time over emission
tmp_FRET_mono=[]; % number of FRET happening via monomers
comb_acc_dist=[];
simulation_randomnum = randperm(1e6,1);
% tmp_loc_dead=[];

%%

if rem(Pulse_dist,time_step)~=0
    error('rem(Pulse_dist,time_step) should be zero')
end


%% making a structure

for pn=1:N
    
    Dectin(pn).type='Donor/Accepror'; % the particles difined as donor or acc
    
    Dectin(pn).x_loc=x_range*(rand(1))+min(x_cont);
    Dectin(pn).y_loc=y_range*(rand(1))+min(y_cont);
    Dectin(pn).phosphorate=[];
    Dectin(pn).state='dead';
    Dectin(pn).flou_decay=0;
    Dectin(pn).pair=[]; % if they get dimer they will have a pair for walking
    Dectin(pn).dimer='Monomer';
    
end




%%
Dectin_frac=round(N*frac);  % the fraction of Dectins which are donors
for i=1:Dectin_frac
    Dectin(i).type=['Donor'];
%     flag_dead_by_FRET(i)=0; % define a flag for donors can have FRET
    
    flag_don = 0;
    while flag_don==0
        tmp_loc=rand(2,1)*x_range+min(x_cont);% defining the pattern with donor and acc coordinations
        
        tmp_flag_acc = [];
        tmp_flag_don = [];
        cnt = 1;
        for l=min(x_cont)+0.1:0.2:max(x_cont)-0.1
            
            for m=min(y_cont)+0.1:0.2:max(y_cont)-0.1
                [tmp_flag_acc(cnt),tmp_flag_don(cnt)]=Myfunc(tmp_loc,[l m],r,init_loc_type); % using function to make the pattern l, m are center coordinate
                cnt=cnt+1;
            end
        end
        flag_don = any(tmp_flag_don);
    end
    
    Dectin(i).x_loc=tmp_loc(1);
    Dectin(i).y_loc=tmp_loc(2);
    all_loc_don(i,:) = tmp_loc;
    
    % display(['Don ' num2str(i) ' out of ' num2str(Dectin_frac)]);
end
%

for j=Dectin_frac+1:N
    Dectin(j).type=['Acceptor'] ;
    Acc_flag(j)=1;
    
    flag_acc = 0;
    while flag_acc==0
        tmp_loc=rand(2,1)*x_range+min(x_cont);
        
        tmp_flag_acc = [];
        tmp_flag_don = [];
        cnt = 1;
        for l=min(x_cont)+0.1:0.2:max(x_cont)-0.1
            
            for m=min(y_cont)+0.1:0.2:max(y_cont)-0.1
                [tmp_flag_acc(cnt),tmp_flag_don(cnt)]=Myfunc(tmp_loc,[l m],r,init_loc_type);
                cnt=cnt+1;
            end
        end
        flag_acc = all(tmp_flag_acc);
    end
    Dectin(j).x_loc=tmp_loc(1);
    Dectin(j).y_loc=tmp_loc(2);
    all_loc_acc(j-Dectin_frac,:) = tmp_loc;
    
    %     display(['Acc ' num2str(j-Dectin_frac) ' out of ' num2str(N-Dectin_frac)]);
    
end

%% Define dimer donors and acceptors
dimers=Dectin_frac*dimer_frac; % define the dimers
dimers_donors=randperm(Dectin_frac,dimers);% define donors that are dimers
dimers_acc=randperm((N-(Dectin_frac)),dimers)+Dectin_frac; % define the acc that are dimers
for i=1:length(dimers_donors)
    Dectin(dimers_donors(i)).dimer=['Dimer'];
    Dectin(dimers_acc(i)).dimer=['Dimer'];
    Dectin(dimers_donors(i)).pair= (dimers_acc(i));
    Dectin(dimers_acc(i)).pair=(dimers_donors(i));
end


%% movie

if makevid==true
    %     f = figure('visible','off','units','normalized','outerposition',[0 0 1 1]);
    f = figure('visible','off','units','normalized','outerposition',[0 0 .7 1]);
    
    if ~isdir('simulation_vidfiles')
        mkdir('simulation_vidfiles');
    end
    
    writerObj = VideoWriter(fullfile('simulation_vidfiles',['simulation1_' num2str(simulation_randomnum) '.avi'])); % Name it.
    close(writerObj);
    writerObj.FrameRate = 5; % How many frames per second.
    open(writerObj);
end
%% saving the innitials
for i = 1 : N
    DectinInit(i) = Dectin(i);
end
%%
time_cnt = 1;
All_dectine{time_cnt} = Dectin;
for tn=0:time_step:t_total
    
    tn/t_total;
    tmp_FRET_eff_mono = [];
    tmp_FRET_eff_dimer=[];
    dead_by_em = 0;
    dead_by_fret = 0;
    dead_donors=0;
    %     tmp_loc_dead=zeros(2,Dectin_frac);
    %     cnt_2=1;
    % clear tmp_loc_dead
    
    
    %% Diffusion part
    dx=randn(N,1).*sqrt(Diffusion.*time_step);  % making random steps
    dy=randn(N,1).*sqrt(Diffusion.*time_step);  % making random steps
    
    
    %%
    for pn=1:N
        %% random walk for donors and acceptors
        Dectin(pn).x_loc= Dectin(pn).x_loc + dx(pn);
        Dectin(pn).y_loc= Dectin(pn).y_loc + dy(pn);
        
        %% make sure dectines are not outside of container
        if Dectin(pn).x_loc < min(x_cont)
            Dectin(pn).x_loc = Dectin(pn).x_loc + x_range;
        elseif Dectin(pn).x_loc > max(x_cont)
            Dectin(pn).x_loc = Dectin(pn).x_loc - x_range;
        end
        
        if Dectin(pn).y_loc < min(y_cont)
            Dectin(pn).y_loc = Dectin(pn).y_loc + y_range;
        elseif Dectin(pn).y_loc > max(y_cont)
            Dectin(pn).y_loc = Dectin(pn).y_loc - y_range;
        end
    end
    %%
    for pn=1:N
        if strcmp(Dectin(pn).dimer,'Dimer') && strcmp(Dectin(pn).type,'Acceptor' )
            % this is a dimer and acceptor
            Dectin(pn).x_loc=Dectin(Dectin(pn).pair).x_loc + dimmers_loc_offset;
            Dectin(pn).y_loc=Dectin(Dectin(pn).pair).y_loc + dimmers_loc_offset;
        end
    end
    
    %% This section activate the pulse and excite the donors. while the pulse is happening
    pulse_flag = 0;
    if rem(tn,Pulse_dist)==0
        pulse_cnt;
        all_Acc_flag(int32(tn/Pulse_dist+1)) = length(find(Acc_flag)); % accumulating the acceptors that are bleached
        all_time_life{int32(tn/Pulse_dist+1)}=time_life;% donors life time
        %         all_FRET_mono{int32(tn/Pulse_dist+1)} = tmp_FRET_eff_mono; % fret efficiency
        %         all_dead_by_FRET{int32(tn/Pulse_dist+1)}=length(find(flag_dead_by_FRET));
        
        time_life=[];
        
        one_pulse_FRET/(one_pulse_FRET+one_pulse_Exc);
        %% bleaching the acc
        %         if tn>0
        %             acc_flag_1_ind = find(Acc_flag);
        %             Dead_Acc=randperm(numel(acc_flag_1_ind),round(numel(acc_flag_1_ind)*Dead_Acc_frac));
        %             Acc_flag(acc_flag_1_ind(Dead_Acc)) = 0;
        %             mean(Acc_flag(Dectin_frac+1:end));
        %
        %         end
        %%
        %         one_pulse_FRET=0;
        %         one_pulse_Exc=0;
        pulse_flag = 1;
        excited_idx=randperm(Dectin_frac,round(Dectin_frac*0.3));
        Alive_Dectin=length(excited_idx);
        decay_tmp=exprnd(2.4e-9, length(excited_idx),1);
        for i=1:length(excited_idx)
            Dectin(excited_idx(i)).state='alive';
            Dectin(excited_idx(i)).flou_decay=decay_tmp(i);
            
        end
        pulse_cnt=pulse_cnt+1;
        if mod(pulse_cnt,100)==0
            fprintf('pulse_cnt = %d\n', pulse_cnt);
        end
        
        
    else
        %% check for emision
        
        for j=1:Dectin_frac
            
            Dectin(j).flou_decay = Dectin(j).flou_decay - strcmp(Dectin(j).state,'alive')*time_step; %allocating a life time to every donor
            if Dectin(j).flou_decay<=0 & strcmp(Dectin(j).state,'alive') % finding the dectin which decay via emmision
                Dectin(j).state='dead';
                Alive_Dectin=Alive_Dectin-1;
                dead_by_em = dead_by_em + 1;
                one_pulse_Exc=one_pulse_Exc+1;
                dead_donors=dead_donors+1;
            end
            
        end
        
        %% FRET
        cnt_Fret=1;
        
        clear acc_x acc_y
        
        for jj=Dectin_frac+1:N
            acc_x(jj-Dectin_frac)=Dectin(jj).x_loc;
            acc_y(jj-Dectin_frac)=Dectin(jj).y_loc;
            
            Dectin(jj).state = 'dead';
        end
        tmp_FRET_eff_mono=zeros(1,Dectin_frac);
        comb_acc_dist=zeros(1, Dectin_frac);
        cnt_2=1;
        
        for j=1:Dectin_frac
            if strcmp(Dectin(j).state,'alive') % compare if the dectin.state is alive
                if strcmp(Dectin(j).dimer, 'Monomer') %compare if the dectin.state is monomer
                    
                    tmp_dis=pdist2([Dectin(j).x_loc; Dectin(j).y_loc]',[acc_x ;acc_y]'); %calculate the  distance betwee donor and acc
                    %                     tmp_dis = tmp_dis + (Acc_flag(Dectin_frac+1:end)==0)*2000;
                    if  any(tmp_dis<=fret_dist)% creteria for choosing the dectin in a certain redius
                        % Inverse_dist=1./tmp_dis;
                        [acc_winner]=find(tmp_dis<=fret_dist); %picking the acceptors in the certain radious
                        % Inverse_dist=1./tmp_dis(acc_winner);
                        %weighted_score=rand(1,length(acc_winner)).*Inverse_dist; %calculation the weigth distance
                        
                        weighted_score=rand(1,length(acc_winner));
                        %[~,I]=max(weighted_score);
                        % Dectin(acc_winner(I)+Dectin_frac).state='alive';% make alive the acc index
                        
                        for l=1:length(weighted_score)
                            Dectin(acc_winner(l)+Dectin_frac).state='alive';% make alive the acc index
                            comb_dis(l)=(tmp_dis(acc_winner(l))^6)/(Forster_dis^6);
                            comb_acc_dist(j)= comb_acc_dist(j)+comb_dis(l);
                            %tmp_FRET_mono(j)=tmp_FRET_mono(j)+tmp_FRET_eff_mono(l);
                        end
                        tmp_FRET_eff_mono(j)=1/(1+comb_acc_dist(j));
                        
                        clear  comb_dis
                        
                        Dectin(j).state='dead';
                        
                        Dectin(j).flou_decay=0;
                        dead_by_fret = dead_by_fret + 1;
                        dead_donors=dead_donors+1;
                        one_pulse_FRET= one_pulse_FRET+1;
                        Alive_Dectin=Alive_Dectin-1;
                        % tmp_FRET_eff_mono(cnt_Fret)=fret_dist.^6/(tmp_dis(acc_winner(I))^6+fret_dist^6);
                        time_cnt;
                        clear acc_winner weighted_score Inverse_dist
                        cnt_Fret=cnt_Fret+1;
                        
                        
                    end
                    
                else
                    Dectin(Dectin(j).pair).state='alive';% make alive the acc index
                    Dectin(j).state='dead';
                    Dectin(j).flou_decay=0;
                    dead_by_fret = dead_by_fret + 1;
                    one_pulse_FRET= one_pulse_FRET+1;
                    Alive_Dectin=Alive_Dectin-1;
                    tmp_FRET_eff_dimer(cnt_Fret)=Forster_dis^6/(Dimer_dist^6+Forster_dis^6);
                    cnt_Fret=cnt_Fret+1;
                    
                end
                
            end
            
            if  strcmp(Dectin(j).state,'dead') && strcmp(prev_state(j).state,'alive') && strcmp(Dectin(j).type,'Donor')
                time_life(j) = tn - floor(tn/Pulse_dist)*Pulse_dist;
                tmp_loc_dead(:,cnt_2)=[Dectin(j).x_loc; Dectin(j).y_loc];
                cnt_2=cnt_2+1;
            end
        end
        
        
        
    end
    
    prev_state = Dectin; % return dectin to their state
    
    %% plot tank
    if makevid==1
        s1 = subplot(3,5,[1 2 6 7 11 12]);
        
        %%
        dectine_plot_type = zeros(1,N);
        for pn=1:N
            if strcmp(Dectin(pn).type,'Donor' )  % compare the characters in Dectin.type to see if they are Donors
                if strcmp(Dectin(pn).state,'alive') && strcmp(Dectin(pn).dimer,'Dimer')
                    dectine_plot_type(pn)=1;
                elseif strcmp(Dectin(pn).state,'dead') && strcmp(Dectin(pn).dimer,'Dimer')
                    dectine_plot_type(pn)=2;
                elseif strcmp(Dectin(pn).state,'alive') && strcmp(Dectin(pn).dimer,'Monomer')
                    dectine_plot_type(pn)=3;
                else
                    dectine_plot_type(pn)=4;
                end
            else
                if strcmp(Dectin(pn).state,'alive') && strcmp(Dectin(pn).dimer,'Dimer')
                    dectine_plot_type(pn)=5;
                elseif strcmp(Dectin(pn).state,'dead') && strcmp(Dectin(pn).dimer,'Dimer')
                    dectine_plot_type(pn)=6;
                elseif strcmp(Dectin(pn).state,'alive') && strcmp(Dectin(pn).dimer,'Monomer')
                    dectine_plot_type(pn)=7;
                else
                    dectine_plot_type(pn)=8;
                end
            end
        end
        %
        tmp_struct = squeeze(struct2cell(Dectin));
        tmp_xloc = cell2mat(tmp_struct(2,:));
        tmp_yloc = cell2mat(tmp_struct(3,:));
        
        hold on
        p1 = plot(tmp_xloc(find(dectine_plot_type==1)),tmp_yloc(find(dectine_plot_type==1))...
            , 'dk','LineWidth',2,'MarkerFaceColor','k','MarkerSize',cont_Msize);
        p2 = plot(tmp_xloc(find(dectine_plot_type==2)),tmp_yloc(find(dectine_plot_type==2))...
            , 'dk','MarkerSize',cont_Msize);
        p3 = plot(tmp_xloc(find(dectine_plot_type==3)),tmp_yloc(find(dectine_plot_type==3))...
            , 'ob','MarkerFaceColor','b','MarkerSize',cont_Msize);
        p4 = plot(tmp_xloc(find(dectine_plot_type==4)),tmp_yloc(find(dectine_plot_type==4))...
            , 'ob','MarkerSize',cont_Msize);
        p5 = plot(tmp_xloc(find(dectine_plot_type==5)),tmp_yloc(find(dectine_plot_type==5))...
            , '*k','LineWidth',2,'MarkerFaceColor','k','MarkerSize',cont_Msize+5);
        p6 = plot(tmp_xloc(find(dectine_plot_type==6)),tmp_yloc(find(dectine_plot_type==6))...
            , 'ok','MarkerSize',cont_Msize);
        p7 = plot(tmp_xloc(find(dectine_plot_type==7)),tmp_yloc(find(dectine_plot_type==7))...
            , '*k','LineWidth',2,'MarkerFaceColor','r','MarkerSize',cont_Msize+5);
        p8 = plot(tmp_xloc(find(dectine_plot_type==8)),tmp_yloc(find(dectine_plot_type==8))...
            , 'or','MarkerSize',cont_Msize-1);
        hold off
        title(['tn=' num2str(tn)],'FontSize',20)
        
        %         lgnd = legend([p1,p2,p3,p4,p5,p6,p7],{'preformed donor-acceptor dimer','non-excited donor dimer',...
        %             'excited donors','non excited donors','non-excited acceptor dimer',...
        %             'excited acceptor','non-excited acceptor'},'Location','south');
        %
        %         lgnd.FontSize = 14;
        %         lgndpos = lgnd.Position;
        %         lgndpos(2) = .01;
        %         lgnd.Position = lgndpos;
        %
        xlim([min(x_cont),max(x_cont)])
        ylim([min(y_cont),max(y_cont)])
        xticks([])
        yticks([])
        bg_c = {'w',[.8 .8 .8]};
        set(gca,'Color',bg_c{pulse_flag+1})
        axis image
        s1pos = s1.Position;
        s1pos(3) = 3;
        s1pos(4) = 3;
        %         s1.Position = s1pos;
        
        clear tmp_dead tmp_struct tmp_type_acc tmp_xloc tmp_yloc
        clear tmp_dead_idx tmp_alive_idx dectine_plot_type lgndpos s1pos
        %%
        %             for pn=1:N
        %                 if strcmp(Dectin(pn).type,'Donor' )  % compare the characters in Dectin.type to see if they are Donors
        %                     if strcmp(Dectin(pn).state,'alive')
        %
        %                         plot(Dectin(pn).x_loc,Dectin(pn).y_loc, 'ob','MarkerFaceColor','b','MarkerSize',cont_Msize)
        %                     else
        %                         plot(Dectin(pn).x_loc,Dectin(pn).y_loc, 'ob','MarkerSize',cont_Msize)
        %                     end
        %                     hold on
        %                 else
        %                     if strcmp(Dectin(pn).state,'alive')
        %                         plot(Dectin(pn).x_loc,Dectin(pn).y_loc, 'or','MarkerFaceColor','r','MarkerSize',cont_Msize)
        %                     else
        %                         plot(Dectin(pn).x_loc,Dectin(pn).y_loc, 'or','MarkerSize',cont_Msize)
        %                     end
        %                     hold on
        %                 end
        %             end
        %             title(['tn=' num2str(tn)])
        %             hold off
        %
        %             xlim([min(x_cont),max(x_cont)])
        %             ylim([min(y_cont),max(y_cont)])
        %             xticks([])
        %             yticks([])
        %             %     hold on
        %             bg_c = {'w','k'};
        %             set(gca,'Color',bg_c{pulse_flag+1})
        %
        %%
        S2=subplot(3,5,[4 5]);
        hold on
        plot(rem(tn,Pulse2plot*Pulse_dist),Alive_Dectin,'bO','MarkerFaceColor','b','MarkerSize',plot_Msize)
        if rem(tn,Pulse_dist)==0
            pp = plot([tn tn],[0 400],'r','LineWidth',2);
            pp.Color(4) = .6;
        end
        xlim([0,Pulse2plot*Pulse_dist])
        ylim([0,Dectin_frac])
        ylabel('excited Donors','FontSize',14)
        hold off
        %         ax = get(gca);ax.XAxis.Exponent = 0;
        
        S3=subplot(3,5,[9 10]);
        hold on
        plot(rem(tn,Pulse2plot*Pulse_dist),log(dead_by_fret+1),'bO','MarkerFaceColor','b','MarkerSize',plot_Msize)
        xlim([0,Pulse2plot*Pulse_dist])
        ylim([0,10])
        ylabel('Number of FRET(log)','FontSize',14)
        hold off
        %         ax = get(gca);ax.XAxis.Exponent = 0;
        
        S4=subplot(3,5,[14 15]);
        hold on
        plot(rem(tn,Pulse2plot*Pulse_dist),log(dead_by_em+1),'bO','MarkerFaceColor','b','MarkerSize',plot_Msize)
        xlim([0,Pulse2plot*Pulse_dist])
        ylim([0,10])
        ylabel('Number of Photon(log)','FontSize',14)
        xlabel('time')
        hold off
        %         ax = get(gca);ax.XAxis.Exponent = 0;
        
        if rem(tn+time_step,Pulse2plot*Pulse_dist)==0
            delete(S2)
            delete(S3)
            delete(S4)
        end
        
        frame = getframe(f); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
        cla(s1);
    end
    %%
    all_dead_by_fret(time_cnt) = dead_by_fret; % num of FRET
    %time_cnt = time_cnt + 1;
    all_FRET_mono{time_cnt} =  tmp_FRET_eff_mono; % FRET Effic
    FRET_eff_dimer{time_cnt}= tmp_FRET_eff_dimer;
    all_Alive_Dectin{time_cnt}=Alive_Dectin;
    if exist('acc_winner')
        all_winner{time_cnt}=(acc_winner)
    end
    if exist('tmp_loc_dead')
        all_dead{time_cnt}=tmp_loc_dead;
        
        clear tmp_loc_dead
    end
    
    
    time_cnt = time_cnt + 1;
    
    %     total_eff{time_cnt}=FRET_eff_mono{time_cnt}+FRET_eff_dimer{time_cnt};
    %  All_dectine{time_cnt} = Dectin;
    
    %     clear tmp_FRET_eff_mono  tmp_FRET_eff_dimer
end

if makevid==true
    close(writerObj);
end

close all
disp('done')
SaveDir= 'D:\PhD\Fret_Simulation\new';

if ~isdir(fullfile(SaveDir,'simulation_matfiles'))
    mkdir(fullfile(SaveDir,'simulation_matfiles'));
end
fname = ['Pixel_' num2str(Pulsenum*(0.025/4.88)) '_initloctype_' num2str(init_loc_type) ...
    '_fretdist_' num2str(fret_dist*1000) '_x_range_' num2str(x_range*10) '_frac_' num2str(frac*10) '__' num2str(simulation_randomnum) '.mat'];

file_name = fullfile(SaveDir,'simulation_matfiles',fname);
save(file_name,'-v7.3')
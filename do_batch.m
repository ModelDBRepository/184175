function do_batch(pars, func, No_xpts)

% do a general batch 
%
% PARS is a cell array
%
% FUNC is a string defining the function to be used
% 
% No_xpts is the number of experiments to do
%

file_id = 1;    % could make this something else if we want to make a logfile of progress
                % could be useful in overcoming loss of connectivity to
                % cluster logged us out?

No_outputs = 0;  

fprintf(file_id, '\n\n\n starting batch \n\n\n');

%%%
sched = findResource('scheduler','type','generic');
set(sched, 'SubmitFcn', @submitFunc);
%%
                
Max_No_concurrent_tasks = 25;
xpt_No_start = 1;
FreedSpace = Max_No_concurrent_tasks;

check_interval = 90;    % seconds between checking for free space
Job_no = 1;


Running_tasks = [];
All_tasks = 1:No_xpts;
Pending_tasks = All_tasks;
Finished_tasks = [];

while ~isempty(Pending_tasks)
    % Break tasks into bactches
    if ~isempty(setdiff(All_tasks, Pending_tasks)) % not first batch

        % Cycle through running tasks until find tasks that are finished:

        fprintf(file_id, '\n\nWaiting for space to become available \n\n');

        finishedtasks = [];
        while isempty(finishedtasks)
            for z=1:length(Running_tasks)
                xpt_No = Running_tasks(z);
                if isequal(task_handle(xpt_No).State, 'finished')
                    finishedtasks = [finishedtasks xpt_No];
                end
            end
            pause(check_interval);	% Wait for a bit...
        end
        % update sets
        Running_tasks = setdiff(Running_tasks, finishedtasks);
        FreedSpace = length(finishedtasks);
        Finished_tasks = union(Finished_tasks, finishedtasks);
        
        % report
        log_string = [datestr(now) ' xpts {' num2str(finishedtasks) '} completed '];
        fprintf(file_id, '%s\n', log_string);
        
        % check for errors in finished tasks
        for j=1:length(finishedtasks)
            xpt_No = finishedtasks(j) ;
            err = get(task_handle(xpt_No), 'ErrorMessage');
            if ~isempty(err)
                fprintf(file_id, 'Error in Xpt %d\n\n', xpt_No);
                fprintf(file_id, '%s\n', err);
            end            
        end
    end
    
    % Create a new SubBatch
    xpt_No_end = xpt_No_start + FreedSpace - 1;
    if xpt_No_end > No_xpts
        xpt_No_end = No_xpts;
    end
    
    SubBatch = [xpt_No_start: xpt_No_end]; 
    
    xpt_No_start = xpt_No_end + 1;
    
    %%%%%%%%%%%%%%% Run Each SubBatch as a Seperate Job %%%%%%%%%%%%%%%%%

    jobs(Job_no) = createJob(sched);
    fprintf(file_id, 'creating tasks for job %d \n \n', Job_no)

    task_string = ['task_handle(xpt_no) = createTask(jobs(Job_no), @' func ', No_outputs, pars{xpt_no});'];
    for h = 1:length(SubBatch);
        xpt_no = SubBatch(h);
        eval(task_string);
        %task_handle(xpt_no) = createTask(jobs(Job_no), @STDE_shen_batch, No_outputs, pars{xpt_no});
    end
    submit(jobs(Job_no));
    
    log_string = [datestr(now) ' Submitting xpts{ ' num2str(SubBatch) '} as job ' num2str(Job_no)];
    fprintf(file_id, '%s\n', log_string);
    
    Running_tasks = union(Running_tasks, SubBatch);
    Pending_tasks = setdiff(Pending_tasks, Running_tasks);
    
    log_string = ['Summary - Pending tasks {' num2str(Pending_tasks) '}'];
    fprintf(file_id, '%s\n', log_string);
    log_string = ['Summary - Running tasks {' num2str(Running_tasks) '}'];
    fprintf(file_id, '%s\n', log_string);
    log_string = ['Summary - Completed tasks {' num2str(Finished_tasks) '}'];
    fprintf(file_id, '%s\n', log_string);
    
    Job_no = Job_no + 1;
end

%%%%%%%%%%%% wait for final running tasks to finish before proceeding  %%%%%%%%%%

while ~isempty(Running_tasks)

    % pick up any finished tasks
    finishedtasks = [];
    for z=1:length(Running_tasks)
        xpt_No = Running_tasks(z);
        if isequal(task_handle(xpt_No).State, 'finished')
            finishedtasks = [finishedtasks xpt_No];
        end

    end

    if ~isempty(finishedtasks)
        Finished_tasks = union(Finished_tasks, finishedtasks);
        Running_tasks = setdiff(Running_tasks, finishedtasks);

        log_string = [datestr(now) ' xpts ' num2str(finishedtasks) ' completed '];
        fprintf(file_id, '%s\n', log_string);
        log_string = ['Summary - Completed tasks {' num2str(Finished_tasks) '}'];
        fprintf(file_id, '%s\n', log_string);

        % check for errors
        for j=1:length(finishedtasks)
            xpt_No = finishedtasks(j) ;
            err = get(task_handle(xpt_No), 'ErrorMessage');
            if ~isempty(err)
                fprintf(file_id, 'Error in Xpt %d\n\n', xpt_No);
                fprintf(file_id, '%s\n', err);
            end
        end
    end

    pause(check_interval);	% Wait for a bit...
end

% final report
log_string = ['Summary - Pending tasks {' num2str(Pending_tasks) '}'];
fprintf(file_id, '%s\n', log_string);
log_string = ['Summary - Running tasks {' num2str(Running_tasks) '}'];
fprintf(file_id, '%s\n', log_string);
log_string = ['Summary - Completed tasks {' num2str(Finished_tasks) '}'];
fprintf(file_id, '%s\n', log_string);

pause(10); % just to be sure .....

all_jobs = get(sched,'Jobs');
destroy(all_jobs);


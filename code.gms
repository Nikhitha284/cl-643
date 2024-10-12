* Include the data file
$include new1.gms
Set
    k /1*5/;  
Parameter M;
M = 100;
Variables
    Cmax              ! Makespan (completion time of the last task)
    delta(t)          ! Duration of task t
    alpha(u,k)        ! Changeover time between positions k and k+1 on unit u
    alpha0(u)         ! Initial changeover time on unit u
    x(t,u,k) binary   ! 1 if task t is assigned to unit u at position k;

Equations
    obj               ! Objective function
    alloc1(u,k)       ! Each position contains at most one task
    alloc2(u,t)       ! Each task is assigned to one position
    continuity(u,k)   ! Continuous usage of event points
    demandeqn(s)         ! Demand satisfaction
    changeover(u,k,t,s) ! Estimating changeover time between tasks
    initial_changeover(u,t) ! Estimating initial changeover time
    processing_time(u)    ! Estimating total processing and changeover time
    task_duration(t)   ! Task duration fits within limits;

* Objective function
obj ..
    Cmax + 1 * sum(t, delta(t)) + 1 * sum(u, alpha0(u) + sum(k, alpha(u,k))) =e= Cmax;

* Constraints
alloc1(u,k) ..
    sum(t, x(t,u,k)) =l= 1;

alloc2(u,t) ..
    sum(k, x(t,u,k)) =l= 1;

continuity(u,k)$(ord(k) > 1) ..
    sum(t, x(t,u,k-1)) =g= sum(t, x(t,u,k));

* Use binary ProducedState(t, s)
demandeqn(s) ..
    sum(t, ProductionRate(t) * delta(t) * ProducedState(t, s)) =g= Demand(s);

* Use binary SuitableUnit(t, u)
changeover(u,k,t,s)$(ord(k) > 1) ..
    alpha(u,k) =g= sum(t, aij(t,s,u) * x(t,u,k-1) * SuitableUnit(t, u)) 
                 - M * (1 - x(s,u,k));

initial_changeover(u,t) ..
    alpha0(u) =g= sum(t, a0i(t,u) * x(t,u,k1) * SuitableUnit(t, u));

processing_time(u) ..
    sum(t, delta(t) * SuitableUnit(t, u)) + alpha0(u) + sum(k, alpha(u,k)) =l= Cmax;

task_duration(t) ..
    MinTime(t) * sum(k, x(t,SuitableUnit(t),k)) =l= delta(t) 
    =l= Tmax(t) * sum(k, x(t,SuitableUnit(t),k));

* Bounds
delta.lo(t) = 0;
alpha.lo(u,k) = 0;
alpha0.lo(u) = 0;

Model scheduling /all/;
Solve scheduling using mip minimizing Cmax;

if (modelStat <> 0) 
    display "Model is infeasible or unbounded.";

display Cmax, delta.l, alpha.l, alpha0.l, x.l;

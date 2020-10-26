
function disp_dr(dr,order,var_list)
% Display the decision rules
%
% INPUTS
%    dr [struct]:            decision rules
%    order [int]:            order of approximation
%    var_list [char array]:  list of endogenous variables for which the
%                            decision rules should be printed

% Copyright (C) 2001-2012 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global M_ options_

nx =size(dr.ghx,2);
nu =size(dr.ghu,2);
if options_.block
    k = M_.nspred;
    k1 = 1:M_.endo_nbr;
else
    k = find(dr.kstate(:,2) <= M_.maximum_lag+1);
    klag = dr.kstate(k,[1 2]);
    k1 = dr.order_var;
end;

if size(var_list,1) == 0
    var_list = M_.endo_names(1:M_.orig_endo_nbr, :);
end

nvar = size(var_list,1);

ivar=zeros(nvar,1);
for i=1:nvar
    i_tmp = strmatch(var_list(i,:),M_.endo_names(k1,:),'exact');
    if isempty(i_tmp)
        disp(var_list(i,:));
        error (['One of the variable specified does not exist']) ;
    else
        ivar(i) = i_tmp;
    end
end

disp('POLICY AND TRANSITION FUNCTIONS')
% variable names
str = '                        ';
for i=1:nvar
    str = [str sprintf('%16s',M_.endo_names(k1(ivar(i)),:))];
end
disp(str);
%
% constant
%
str = 'Constant            ';

decision = [];
for i=1:nvar
    x = dr.ys(k1(ivar(i)));
    if order > 1
        x = x + dr.ghs2(ivar(i))/2;
    end
    str = [str sprintf('%16.6f',x)];
    decision = [decision x];
end

disp(str)

if order > 1
    decisiontemp = [];
    str = '(correction)        ';
    for i=1:nvar
        x = dr.ghs2(ivar(i))/2;
        str = [str sprintf('%16.6f',x)];
        decisiontemp = [decisiontemp x];
    end
    disp(str)
    decision = [decision;decisiontemp];
end
%
% ghx
%
for k=1:nx
    decisiontemp = [];
    if isfield(dr,'state_var')
        str1 = subst_auxvar(dr.state_var(k),-1);
    else
        str1 = subst_auxvar(k1(klag(k,1)),klag(k,2)-M_.maximum_lag-2);
    end
    str = sprintf('%-20s',str1);
    for i=1:nvar
        x = dr.ghx(ivar(i),k);
        str = [str sprintf('%16.6f',x)];
        decisiontemp = [decisiontemp x];
    end
    disp(str)
    decision = [decision;decisiontemp];
end
%
% ghu
%
for k=1:nu
    decisiontemp = [];
    str = sprintf('%-20s',M_.exo_names(k,:));
    for i=1:nvar
        x = dr.ghu(ivar(i),k);
        str = [str sprintf('%16.6f',x)];
        decisiontemp = [decisiontemp x];
    end
    disp(str)
    decision = [decision;decisiontemp];
end

if order > 1
    % ghxx
    for k = 1:nx
        for j = 1:k
            decisiontemp = [];
            str1 = sprintf('%s,%s',subst_auxvar(k1(klag(k,1)),klag(k,2)-M_.maximum_lag-2), ...
                           subst_auxvar(k1(klag(j,1)),klag(j,2)-M_.maximum_lag-2));
            str = sprintf('%-20s',str1);
            for i=1:nvar
                if k == j
                    x = dr.ghxx(ivar(i),(k-1)*nx+j)/2;
                else
                    x = dr.ghxx(ivar(i),(k-1)*nx+j);
                end
                str = [str sprintf('%16.6f',x)];
                decisiontemp = [decisiontemp x];                
            end
            disp(str)
            decision = [decision;decisiontemp];
        end
    end
    %
    % ghuu
    %
    for k = 1:nu
        for j = 1:k
            decisiontemp = [];
            str = sprintf('%-20s',[M_.exo_names(k,:) ',' M_.exo_names(j,:)] );
            for i=1:nvar
                if k == j
                    x = dr.ghuu(ivar(i),(k-1)*nu+j)/2;
                else
                    x = dr.ghuu(ivar(i),(k-1)*nu+j);
                end
                
                str = [str sprintf('%16.6f',x)];
                decisiontemp = [decisiontemp x];
            end
            disp(str)
            decision = [decision;decisiontemp];
        end
    end
    %
    % ghxu
    %
    for k = 1:nx
        for j = 1:nu
            decisiontemp = [];
            str1 = sprintf('%s,%s',subst_auxvar(k1(klag(k,1)),klag(k,2)-M_.maximum_lag-2), ...
                           M_.exo_names(j,:));
            str = sprintf('%-20s',str1);
            for i=1:nvar
                x = dr.ghxu(ivar(i),(k-1)*nu+j);
                str = [str sprintf('%16.6f',x)];
                decisiontemp = [decisiontemp x];
            end
            disp(str)
            decision = [decision; decisiontemp];
        end
    end
end
save polfunction decision
end

% Given the index of an endogenous (possibly an auxiliary var), and a
% lead/lag, creates a string of the form "x(lag)".
% In the case of auxiliary vars for lags, replace by the original variable
% name, and compute the lead/lag accordingly.
function str = subst_auxvar(aux_index, aux_lead_lag)
global M_

if aux_index <= M_.orig_endo_nbr
    str = sprintf('%s(%d)', deblank(M_.endo_names(aux_index,:)), aux_lead_lag);
    return
end
for i = 1:length(M_.aux_vars)
    if M_.aux_vars(i).endo_index == aux_index
        switch M_.aux_vars(i).type
          case 1
            orig_name = deblank(M_.endo_names(M_.aux_vars(i).orig_index, :));
          case 3
            orig_name = deblank(M_.exo_names(M_.aux_vars(i).orig_index, :));
          case 4
            str = sprintf('EXPECTATION(%d)(...)', aux_lead_lag);
            return
          case 6
            str = sprintf('%s(%d)', ...
                          deblank(M_.endo_names(M_.aux_vars(i).endo_index, :)),aux_lead_lag);
            return
          otherwise
            error(sprintf('Invalid auxiliary type: %s', M_.endo_names(aux_index, :)))
        end
        str = sprintf('%s(%d)', orig_name, M_.aux_vars(i).orig_lead_lag+aux_lead_lag);
        return
    end
end
error(sprintf('Could not find aux var: %s', M_.endo_names(aux_index, :)))
end

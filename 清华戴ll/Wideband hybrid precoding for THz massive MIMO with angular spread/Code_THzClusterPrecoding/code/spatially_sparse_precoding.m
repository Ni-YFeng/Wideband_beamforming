function [F, F_BB, Theta]=spatially_sparse_precoding(Nr,N_s,F_opt,At)
F=[];
Theta = [];
F_res=F_opt;
for i=1:Nr
    pha=At'*F_res;
    [~,k]=max(diag(pha*pha'));
%     disp(max(diag(pha*pha')));
    Theta = [Theta k];
    F=[F At(:,k)];
    F_BB=inv(F'*F)*F'*F_opt;
    F_res=(F_opt-F*F_BB)/norm(F_opt-F*F_BB,'fro');
end
%F_BB = F_BB/norm(F*F_BB,'fro');
F_BB = sqrt(Nr)*F_BB/norm(F*F_BB,'fro');
error_sparse = norm(F_opt - F*F_BB,'fro')^2;
norm_opt = norm(F_opt,'fro')^2;



    
    

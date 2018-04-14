function plot_val(valap, Objectivefunc)
    ep = linspace(1, 10, 10);    
    figure(1)
    plot(ep, valap);
    title('Val_Ap Values');
    xlabel('Iterations');
    ylabel('ap');
    
    figure(2)
    plot(ep, Objectivefunc);
    title('Objectivefunc Values');
    xlabel('Iterations');
    ylabel('Objectivefunc');
end
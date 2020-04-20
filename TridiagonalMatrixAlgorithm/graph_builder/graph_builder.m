% graph builder
function y = Graph_builder(A)
    x = A(:,1:1);
   
    y = A(:,2:3);
    
    
    %disp(y);
    
    plot(x, y, 'LineWidth', 0.3);
   
   
    %axis([410 470 0.2 0.6]);
    
    
    str_arr = string(2);
    str_arr(1) = ' - Approximate solution';
    str_arr(2) = ' - Exact Solution';
   
    grid on;
    title('Lock')
    xlabel('threads')
    ylabel('time, sec')
    lgd = legend(str_arr, 'Location','northwest');
    lgd.FontSize = 8;
    
    print('Graphs.png','-dpng','-r300');  
end


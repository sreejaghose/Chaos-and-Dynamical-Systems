
clear all
chen();
function chen()

    % initial conditions
    x0 = -0.1;
    y0 = 0.5;
    z0 = -0.6;
    % numerically solve 
    [t y] = ode45(@chen_ode, [0 1e2], [x0; y0; z0]);
    
    % plot
    figure
    plot3(y(:,1),y(:,2),y(:,3));
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('Chen double scroll attractor')  
    
    function dudt = chen_ode(t,y)
        % parameters
        a = 40;
        b = 3;
        c = 28;
        
        dudt = [a*(y(2)-y(1));
                (c-a)*y(1) - y(1)*y(3) + c*y(2);
                y(1)*y(2) - b*y(3)];
    end
end



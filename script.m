clear; clc; format short g;

%% DATA

data = readtable('data.xlsx');
datapoints = 4 * height(data);
[C, S, R, L] = deal(zeros(datapoints,1)); % Curve, Sign, Radius, Length

dic1 = ["Αριστερή", "Αριστερη", "αριστερή", "αριστερη"];	% LEFT
dic2 = ["Δεξιά", "Δεξια", "δεξιά", "δεξια"];				% RIGHT

i=1;
while i <= height(data)
    % line
    C(4*i-3) = 1;
    if i == 1
        L(1) = 0; 
    else
        L(4*i-3) = (data{i,1}-data{i-1,2}) * 1000;
    end
    % cubpar 
    C(4*i-2) = 2; C(4*i-0) = 2;
    R(4*i-2) = data{i,3}; R(4*i-0) = data{i,3};
    L(4*i-2) = data{i,4}; L(4*i-0) = abs(data{i,6});
    % circle
    C(4*i-1) = 3;
    R(4*i-1) = data{i,3};
    L(4*i-1) = abs(data{i,5});
    if any(string(data{i,7}) == dic1)
        [S(4*i-0), S(4*i-1), S(4*i-2)] = deal(+1);
    elseif any(string(data{i,7}) == dic2)
        [S(4*i-0), S(4*i-1), S(4*i-2)] = deal(-1);
    else
        fprintf('\n'); error('ERROR: data{i,7} [data.xlsx]');
    end
    i = i + 1;
end

%% INI

gridpoints = datapoints + 1;

X = zeros(gridpoints,1); % longitude
Y = zeros(gridpoints,1); % latitude
E = zeros(gridpoints,1); % easting

X(1) = 0; Y(1) = 0; E(1) = 0; res = 20;
figure(); hold on; grid on; axis equal; axis tight;

%% SOL

for j = 1:datapoints
    switch C(j)
        case 1
            DT = L(j);
            DN = 0;
            DE = 0;
            dt = linspace(0,DT,res);
            dn = zeros(1,res);
            style = 'g:';
        case 2
            a = 6*L(j)*R(j);
            F = @(x) sqrt( 1 + (3/a * x.^2).^2 );
            G = @(x) integral(F,0,x) - L(j);
            DT = fzero(G,[0,L(j)]);
            DN = 1/a * DT^3;
            DE = rad2deg(3/a * DT^2);
            if mod(j,4) == 2			% ENTRY
            	dt = linspace(0,DT,res);
				dn = 1/a * dt.^3;
            else						% EXIT
                dt = linspace(DT,0,res);
                dn = 1/a * dt.^3;
            end
            style = 'c-.';
        case 3
            DE = rad2deg(L(j)/R(j));
            DT = R(j) * sind(DE);
            DN = R(j) * (1 - cosd(DE));
            dt = R(j) .* sind(linspace(0,DE,res));
            dn = R(j) .* (1 - cosd(linspace(0,DE,res)));
            center = rotz(E(i)) * [0; S(i)*R(i); 0];
            plot(X(i)+center(1), Y(i)+center(2), 'y+', 'MarkerSize', 2);
            style = 'y-';           
    end
    xy = [cosd(E(j)), -sind(E(j)); sind(E(j)), cosd(E(j))] * [dt; S(j)*dn];
    XY = [cosd(E(j)), -sind(E(j)); sind(E(j)), cosd(E(j))] * [DT; S(j)*DN];
    plot(X(j)+xy(1,:), Y(j)+xy(2,:), style);
    X(j+1) = X(j) + XY(1);
    Y(j+1) = Y(j) + XY(2);
    E(j+1) = E(j) + S(j)*DE;
end

plot(X, Y, 'ys', 'MarkerSize', 5, 'LineWidth', 1); % disp([X,Y,E]);
